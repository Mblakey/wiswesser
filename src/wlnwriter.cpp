#ifdef USING_OPENBABEL
/*********************************************************************
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef USING_OPENBABEL
#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/kekulize.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include "openbabel/stereo/stereo.h"
#include <openbabel/stereo/tetrahedral.h>

using namespace OpenBabel; 
#endif

#include "wlnparser.h"

// Element "magic numbers"
#define DUM   0
#define BOR   5
#define CAR   6
#define NIT   7
#define OXY   8
#define FLU   9
#define PHO   15
#define SUL   16
#define CHL   17
#define BRO   35
#define IOD   53


#define INT_TO_LOCANT(X) (X+64)
#define LOCANT_TO_INT(X) (X-64)


#if 0
struct LocantPos{
  unsigned int locant;  // lets branch locants be placed at the back of the array, past 255 indexing
  OBAtom *atom;         // atoms are mallocd by obabel so should be alive at all times
};

struct PathData {
  LocantPos *locant_path = 0;
  unsigned int path_size = 0;
  bool macro_ring = false;
};


/* allow a pass by reference for exiting the local sssr */
struct SubsetData{
  unsigned int path_size = 0;
  bool hetero = false;
  bool bridging = false; 
  bool multi = 0;
};


static void Fatal(const char *str){
  fprintf(stderr,"Fatal: %s\n",str);
  exit(1);
}


/**********************************************************************
                          Locant Path Functions
**********************************************************************/

template <typename T>
void burn_stack(std::stack<T> &stack){
  while(!stack.empty())
    stack.pop(); 
}


void copy_locant_path(LocantPos*new_path, LocantPos*locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++){
    new_path[i].atom = locant_path[i].atom;
    new_path[i].locant = locant_path[i].locant; 
  }
}

void zero_locant_path(LocantPos *locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++){
    locant_path[i].atom = 0;
    locant_path[i].locant = 0; 
  }
}


bool in_locant_path(OBAtom *atom,LocantPos*locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++){
    if(atom == locant_path[i].atom)
      return true; 
  }
  return false;
}


bool char_in_locant_path(unsigned int locant,LocantPos*locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++){
    if(locant == locant_path[i].locant)
      return true; 
  }
  return false;
}


unsigned int position_in_path(OBAtom *atom,LocantPos*locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++){
    if(atom == locant_path[i].atom)
      return i; 
  }

  Fatal("Error: atom not found in locant path");
  return 0; 
}


int fusion_locant(OBMol *mol,OBRing *ring, LocantPos*locant_path, unsigned int path_size){
  unsigned int lpos = path_size; 
  for(unsigned int i=0;i<ring->Size();i++){
    OBAtom *latom = mol->GetAtom(ring->_path[i]); 
    if(in_locant_path(latom,locant_path,path_size)){
      unsigned int pos = position_in_path(latom,locant_path,path_size); 

      if(pos < lpos)
        lpos = pos;
    }
  }
  return lpos; 
}


/* overall ring sum, combinates rule 30d and 30e for symmetrical structures */
unsigned int ring_sum(OBMol *mol, OBRing *ring, LocantPos*locant_path,unsigned int path_size){
  unsigned int rsum = 0;
  for(unsigned int i=0;i<ring->Size();i++){
    OBAtom *ratom = mol->GetAtom(ring->_path[i]);
    if(in_locant_path(ratom,locant_path, path_size))
      rsum += position_in_path(ratom, locant_path,path_size) + 1;
  }
  return rsum;
}


unsigned int fusion_sum(OBMol *mol, LocantPos*locant_path, unsigned int path_size, std::set<OBRing*> &local_SSSR){
  unsigned int ret = 0; 
  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end();riter++)
    ret += fusion_locant(mol,*riter,locant_path,path_size) + 1; // A=1, B=2 etc 
  return ret;
}





void write_ring_size(OBRing *ring, std::string &buffer){
  if(ring->Size() < 9)
    buffer += ring->Size() + '0'; 
  else{
    buffer += '-';
    buffer += std::to_string(ring->Size()); 
    buffer += '-';
  }

}

/* calculates the locant path A->X site for a broken locant */
unsigned char broken_parent_char(unsigned int locant){
  for(unsigned char ch = 'A'; ch < 'X';ch++){ // could expand this later on
    if(locant >= 128+(LOCANT_TO_INT(ch)*6) && locant < 128+(LOCANT_TO_INT(ch)*6) +6) 
      return ch; 
  }
  return 0; 
}


void static write_locant(unsigned int locant,std::string &buffer){
  if(locant >= 128){
    unsigned int loc_start = 0; 
    unsigned int offset = 0; 
    
    loc_start = broken_parent_char(locant); 
    if(!loc_start)
      Fatal("could not fetch off path parent for broken locant"); 

    offset = locant - (128 + (LOCANT_TO_INT(loc_start)*6)); // 0 = E-, 1 = E-&
    buffer += loc_start;

    switch (offset) {
      case 0:
        buffer +='-';
        break;
      case 1:
        buffer +="-&";
        break;
      case 2:
        buffer += "--";
        break;
      case 3:
        buffer += "--&";
        break;
      case 4:
        buffer += "-&-"; 
        break;
      case 5:
        buffer += "-&&";
        break;

      default:
        Fatal("broken locants exceeding tree limit of 6"); 
    }
  }
  else{
    if(locant < 'X')
      buffer += locant;
    else{
      unsigned int amps = 0;
      while(locant >= 'X'){
        amps++; 
        locant+= -23;
      }
      buffer += locant; 
      for(unsigned int i=0;i<amps;i++)
        buffer += '&';
    }
  }
}

static void print_locant_array(LocantPos* locant_path, unsigned int size,bool locants=false){
  fprintf(stderr,"[ ");
  for(unsigned int i=0; i<size;i++){
    if(!locant_path[i].atom)
      fprintf(stderr,"0 ");
    else{
      if(locants)
        fprintf(stderr,"%d ",locant_path[i].locant);
      else
        fprintf(stderr,"%d ",locant_path[i].atom->GetIdx());
    }
  }
    
  fprintf(stderr,"]\n");
}


void sort_locants(unsigned char *arr,unsigned int len){
	for (unsigned int j=1;j<len;j++){
		unsigned char key = arr[j];
		int i = j-1;
		while(i>=0 && arr[i] > key){
			arr[i+1] = arr[i];
			i--;
		}
		arr[i+1] = key;
	}
}


/* this does the same calculation as done in readwln, broken locants are given values
 * from 128 in the form of a binary tree, with a max tree depth of two, spawning 6
 * potential branches when needed, see readwln for tree structure and notes */
void update_broken_locants(OBMol* mol, LocantPos *ratom,
                          LocantPos *locant_path, unsigned int path_size,
                          LocantPos *off_branches, unsigned int off_branch_n)
{
  // if on E, then check for E-, if E- exists, create E-&
  // if either of these exist, make E--, E--& etc. 

  // since the array is layed out as follow: 
  // E-, E-&, E--, E--&, E-&-, E-&&
  for(unsigned int b=0;b<off_branch_n;b++){
    if( mol->GetBond(ratom->atom,off_branches[b].atom)){

      // this could be E- or E-&, where E- should always get set first. 
      unsigned int movement = 1;
      unsigned int offset = 1; // skip the E-& 
      unsigned int off_path_char = 128 + (LOCANT_TO_INT(ratom->locant) * 6);
      // check the atom, does it already have a E-, if so E-&, this has to be exact.
      if(off_branches[b].locant == off_path_char){
        off_path_char++; // move it to E-&;
        movement = 1; 
        offset = 0; // E-& E-&- and E-&& are all in a line
      }

      if(!off_branches[b].locant || off_branches[b].locant > off_path_char)
        off_branches[b].locant = off_path_char; 
      
      for(unsigned int bn=0;bn<off_branch_n;bn++){ // and update its children in the tree
       if(mol->GetBond(off_branches[b].atom, off_branches[bn].atom)){
          //if(!off_branches[bn].locant || off_branches[bn].locant > off_branches[b].locant)
            off_branches[bn].locant = off_branches[b].locant + offset + movement++;
        }
      }
    }
  }
 
}

/* writes the lowest locant in the ring given, an unspecified broken locant is given 'A'-1, since
 * this is impossible under normal operations, if this locant exists, it needs specifying, and then checking
 * whether its actually the lowest locant possible in that chain */
unsigned int lowest_ring_locant(OBMol*mol, OBRing *ring, LocantPos* locant_path, unsigned int plen){
  unsigned int lowest_locant = 0; 
  for(unsigned int i=0;i<plen;i++){
    if(locant_path[i].atom && ring->IsMember(locant_path[i].atom)
       && (!lowest_locant || locant_path[i].locant < lowest_locant)){
        lowest_locant = locant_path[i].locant;
    }
  }
  return lowest_locant; 
}


unsigned int lowest_ring_locantWB(OBMol*mol, OBRing *ring, 
                                  LocantPos* locant_path, unsigned int plen,
                                  LocantPos* branching_path, unsigned int blen){
  unsigned int lowest_locant = 0; 
  for(unsigned int i=0;i<plen;i++){
    if(locant_path[i].atom && ring->IsMember(locant_path[i].atom)
       && (!lowest_locant || locant_path[i].locant < lowest_locant)){
        lowest_locant = locant_path[i].locant;
    }
  }

  for(unsigned int i=0;i<blen;i++){
    if(branching_path[i].atom && ring->IsMember(branching_path[i].atom)){
      unsigned int parent_ch = broken_parent_char(branching_path->locant); 
      if(parent_ch && parent_ch < lowest_locant)
        lowest_locant = branching_path[i].locant; 
    }
  }
  
  return lowest_locant; 
}


bool IsConsecutiveLocants(LocantPos *a, LocantPos *b){
  if(a->locant == b->locant+1)
    return true;
  else if(a->locant == b->locant-1)
    return true;
  return false; 
}

/*
 * Performs a santity check on the ring set, if the starting locant is impossible
 * under path walk conditions, it must be a pseudo locant and be given a pseudo code - 
 * This acts as a fail safe for all compounds, but its useage should be rare. 
 * - This must be the full path size, detection made on locant value for bridging broken points
*/
void TestPathSequences( OBMol *mol,LocantPos*locant_path, unsigned int path_size, 
                        std::vector<OBRing*> &ring_order,std::map<OBAtom*,unsigned int> &atom_shares,
                        std::map<OBAtom*,bool> &bridge_atoms, std::string &buffer){
  
  std::map<OBBond*,bool> allowed_jumps; 
  std::map<unsigned int,unsigned int> allowed_connection; 

  // set up an allowed connections similar to read logic, but only needed on first char
  for(unsigned int i=0;i<path_size;i++){
    unsigned int a = locant_path[i].locant; 
    if(i==0 || i == path_size-1)
      allowed_connection[a] = 2;
    else
      allowed_connection[a] = 1;

    if(atom_shares[locant_path[i].atom] > 3)
      allowed_connection[a]++;

    if(bridge_atoms[locant_path[i].atom])
      allowed_connection[a]--; 
  }

  for(unsigned int r=0;r<ring_order.size();r++){
    OBRing *ring = ring_order[r]; 
    LocantPos *sequence = (LocantPos*)malloc(sizeof(LocantPos)*ring->Size()); 
    for(unsigned int i=0;i<ring->Size();i++){
      sequence[i].atom    = locant_path[position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)].atom; 
      sequence[i].locant  = locant_path[position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)].locant; 
    }

    unsigned int lowest = lowest_ring_locant(mol, ring, sequence, ring->Size()); 
    while(sequence[0].locant != lowest){
      LocantPos tmp = sequence[0]; 
      for(unsigned int i=0;i<ring->Size()-1;i++){
        sequence[i] = sequence[i+1]; 
      }
      sequence[ring->Size()-1] = tmp; 
    }

    while(allowed_connection[lowest] < 1)
      lowest++; 

    // always check that the ends first, as this takes highest priotrity due to fusion sum
    if(!IsConsecutiveLocants(&sequence[0], &sequence[ring->Size()-1]) && 
       !allowed_jumps[mol->GetBond(sequence[0].atom,sequence[ring->Size()-1].atom)])
    { 
      allowed_jumps[mol->GetBond(sequence[0].atom,sequence[ring->Size()-1].atom)] = true; 
    }

    for(unsigned int k=1;k<ring->Size();k++){
      if( mol->GetBond(sequence[k].atom,sequence[k-1].atom) && 
          !IsConsecutiveLocants(&sequence[k], &sequence[k-1]) && 
          !allowed_jumps[mol->GetBond(sequence[k].atom,sequence[k-1].atom)] 
          && (sequence[k].locant < 128 && sequence[k-1].locant < 128)){

        if(sequence[k-1].locant == lowest || sequence[k].locant==lowest)
          allowed_jumps[mol->GetBond(sequence[k].atom,sequence[k-1].atom)] = true; 
        else{
          if(sequence[k].locant < sequence[k-1].locant){
            buffer += '/'; 
            write_locant(sequence[k].locant, buffer); 
            write_locant(sequence[k-1].locant, buffer); 
          }
          else{
            buffer += '/'; 
            write_locant(sequence[k-1].locant, buffer); 
            write_locant(sequence[k].locant, buffer); 
          }
        }
      }
    }

    allowed_connection[lowest]--; 
    free(sequence);
  }
}


#if DEPRECATED
bool sequential_chain(  OBMol *mol,OBRing *ring, 
                        OBAtom **locant_path, unsigned int path_size,
                        std::map<unsigned char, bool> &in_chain)
{
  unsigned char *sequence = (unsigned char*)malloc(sizeof(unsigned char)*ring->Size()); 
  for(unsigned int i=0;i<ring->Size();i++){
    if(in_locant_path(mol->GetAtom(ring->_path[i]),locant_path, path_size)){
      unsigned int pos = position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size);
      sequence[i] = INT_TO_LOCANT(pos+1); 
    }
    else
      sequence[i] = 0; 
  }

  sort_locants(sequence,ring->Size());
  unsigned char prev = 0; 
  for(unsigned int k=0;k<ring->Size();k++){
    if(prev && prev != sequence[k]-1){
      // if we've made a jump, can we get to the jump via locants already wrapped
      for(unsigned char loc = prev+1;loc<sequence[k];loc++){
        if(!in_chain[loc]){
          free(sequence);
          return false; // if we cant do it, return false
        }
      }
    }
    prev = sequence[k];
  }
    
  free(sequence);
  return true;
}


unsigned int PseudoCheck( OBMol *mol, OBAtom **locant_path, unsigned int path_size,
                          std::vector<unsigned char>      &locant_order,
                          std::vector<OBRing*>            &ring_order,
                          std::map<OBAtom*,bool>          &bridge_atoms,
                          std::string &buffer)
{

  unsigned int pseudo_pairs = 0;

  // set up the read shadowing algorithm
  std::map<unsigned char,unsigned int> connections;
  std::map<unsigned char, unsigned char> highest_jump;

  for(unsigned int i=0;i<path_size;i++){
    unsigned char ch = INT_TO_LOCANT(i+1); 
    connections[ch] = 1;
    if(i==0 || i == path_size-1)
      connections[ch]++;
    
    if(bridge_atoms[locant_path[i]] && connections[ch])
      connections[ch]--;

    if(locant_path[i]->GetAtomicNum() == 6){
      unsigned int rbonds = 0; 
      for(unsigned int k=0;k<path_size;k++){
        if(mol->GetBond(locant_path[i],locant_path[k]))
          rbonds++;
      }

      if(rbonds==4)
        connections[INT_TO_LOCANT(i+1)] = 4;
    }
  }

  std::map<std::set<unsigned char>,bool> seen_nt;
  std::vector<std::set<unsigned char>> non_trivials;

  for(unsigned int i=0;i<path_size;i++){
    for(unsigned int j=i+2;j<path_size;j++){
      if(mol->GetBond(locant_path[i],locant_path[j]))
        non_trivials.push_back({INT_TO_LOCANT(i+1),INT_TO_LOCANT(j+1)});
    }
  }


  // shadow read graph traversal
  for(unsigned int i=0;i<ring_order.size();i++){
    unsigned int steps = ring_order[i]->Size()-1;
    unsigned char bind = locant_order[i];
    unsigned char locant = locant_order[i];
    std::vector<unsigned char> path;

    path.push_back(bind);
    for(unsigned int step = 0; step < steps;step++){
      if(highest_jump[locant])
        locant = highest_jump[locant];
      else if(locant < INT_TO_LOCANT(path_size))
        locant++;

      path.push_back(locant);
    }

    // add the loop back logic
    for(unsigned int i=0;i<path.size();i++){
      unsigned int tally = 1;
      if(path[i] == INT_TO_LOCANT(path_size)){
        for(unsigned int j=i+1;j<path.size();j++){
          if(path[j] == path[i]){
            path[j] += -tally; // looping back
            tally++;
          }
        }
      }
    }

    while(!connections[bind] && bind < INT_TO_LOCANT(path_size)){
      bind++; // increase bind_1
      bool found = false;
      for(unsigned int a=0;a<path.size();a++){
        if(path[a] == bind)
          found = true;
      }
      // if its already there, we dont change the path
      if(!found && !path.empty()){ 
        // if its not, we have to spawn it in by knocking one off
        path.pop_back();
        locant = path.back(); 
      }
    }
      
    highest_jump[bind] = locant;
    seen_nt[{bind,locant}] = true;


    if(connections[bind])
      connections[bind]--;
    if(connections[locant])
      connections[locant]--;
  }

  // now we check if this was possible without pseudo locants
  std::set<unsigned char>::iterator psd_iter;
  for(unsigned int i=0;i<non_trivials.size();i++){
    std::set<unsigned char> nt_bond = non_trivials[i];
    if(!seen_nt[nt_bond]){
      psd_iter = nt_bond.begin();
      buffer += '/';
      write_locant(*psd_iter,buffer); 
      write_locant(*(++psd_iter),buffer); 
      pseudo_pairs++;
    }
  }

  return pseudo_pairs;
}
#endif

bool ReachableFromEntry(OBAtom *entry, std::set<OBAtom*> &ring_atoms, std::set<OBAtom*> &seen){
  std::stack<OBAtom*> stack;  
  OBAtom *top = 0; 
  stack.push(entry);
  while(!stack.empty()){
    top = stack.top();
    seen.insert(top);
    stack.pop();

    FOR_NBORS_OF_ATOM(n ,top){
      OBAtom *nbor = &(*n); 
      for(std::set<OBAtom*>::iterator fiter = ring_atoms.begin(); fiter != ring_atoms.end();fiter++){
        if(*fiter == nbor)
          stack.push(nbor);
      }
    }
  }

  return (seen == ring_atoms);
}

// helper function whenever we need the ring bonds
void FillRingBonds(OBMol *mol,OBRing* obring, std::set<OBBond*> &ring_bonds){
   
  for(unsigned int i=0;i<obring->Size()-1;i++){
    OBAtom *satom = mol->GetAtom(obring->_path[i]);    
    OBAtom *eatom = mol->GetAtom(obring->_path[i+1]);
    OBBond *bond = mol->GetBond(satom,eatom); 
    if(bond)
      ring_bonds.insert(bond); 
  }

  OBAtom *satom = mol->GetAtom(obring->_path[0]);    
  OBAtom *eatom = mol->GetAtom(obring->_path[obring->Size()-1]);
  OBBond *bond = mol->GetBond(satom,eatom); 
  if(bond)
    ring_bonds.insert(bond); 
}



#if DEPRECATED
/* read locant path algorithm, we return the number of non consecutive blocks,
pseudo check will add determined pairs and check notation is viable for read
 - This is deprecated as rings are read on path read, but keep for posterity */
unsigned int ReadLocantPath(  OBMol *mol, OBAtom **locant_path, unsigned int path_size,
                      std::set<OBRing*>               &local_SSSR,
                      std::map<OBAtom*,bool>          &bridge_atoms,
                      std::vector<OBRing*>            &ring_order,
                      std::string &buffer,
                      bool verbose)
{  
  unsigned int arr_size = 0; 
  OBRing **ring_arr = (OBRing**)malloc(sizeof(OBRing*) * local_SSSR.size()); 
  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end();riter++)
    ring_arr[arr_size++]= *riter; 
  
  unsigned int assignment_score = 0;
  unsigned int rings_done = 0;
  std::vector<unsigned char> locant_order;
  std::map<unsigned char,bool> in_chain; 
  std::map<unsigned int,bool> pos_written; 
  while(rings_done < arr_size){
    
    unsigned int pos_to_write     = 0;
    unsigned char lowest_in_ring  = 255;
    unsigned char highest_in_ring  = 0;
    unsigned int lowest_rsum = UINT32_MAX; 
    
    bool updated = false;
    for(unsigned int i=0;i<arr_size;i++){
      OBRing *wring = ring_arr[i]; 
      if(!pos_written[i]){
        if(sequential_chain(mol,wring,locant_path,path_size,in_chain)){
          updated = true;
          unsigned char min_loc = 255; 
          unsigned char high_loc = 0; 
          unsigned int rsum = ring_sum(mol,wring,locant_path,path_size);
          for(unsigned int k=0;k<wring->Size();k++){
            if(in_locant_path(mol->GetAtom(wring->_path[k]),locant_path,path_size)){
              unsigned int pos = position_in_path(mol->GetAtom(wring->_path[k]),locant_path,path_size);
              unsigned char loc = INT_TO_LOCANT(pos+1);

              if(loc < min_loc)
                min_loc = loc; 
              if(loc > high_loc)
                high_loc = loc;
            }
          }

          if(min_loc < lowest_in_ring || 
            (min_loc == lowest_in_ring && rsum < lowest_rsum) || 
            (min_loc == lowest_in_ring && rsum == lowest_rsum && high_loc < highest_in_ring)){
            lowest_in_ring = min_loc;
            highest_in_ring = high_loc;
            lowest_rsum = rsum;
            pos_to_write = i; 
          }
        }
      }
    }

    // catch all whilst branching locant logic is not currently active
    if(!updated){
      free(ring_arr);
      return 255;
    }

    OBRing *to_write = ring_arr[pos_to_write];
    if(!to_write)
      Fatal("out of access locant path reading");
    

    for(unsigned int k=0;k<to_write->Size();k++){
      if(in_locant_path(mol->GetAtom(to_write->_path[k]),locant_path,path_size)){ 
        unsigned char loc = INT_TO_LOCANT(position_in_path(mol->GetAtom(to_write->_path[k]),locant_path,path_size)+1); 
        in_chain[loc] = true;
      }
    }

    if(OPT_DEBUG && verbose){
      fprintf(stderr,"  %d(%d): %c(%d) -",rings_done,pos_to_write,lowest_in_ring,lowest_in_ring);
      TestPathSequences(mol,to_write,locant_path,path_size,false);
    }

    if(lowest_in_ring != 'A'){
      buffer += ' ';
      write_locant(lowest_in_ring,buffer);
      assignment_score++;
    }
    else if(!rings_done)
      assignment_score++;
  
    write_ring_size(to_write, buffer); 

    locant_order.push_back(lowest_in_ring);
    ring_order.push_back(to_write);
    
    pos_written[pos_to_write] = true;
    rings_done++;
  }

  unsigned int pairs = PseudoCheck(mol,locant_path,path_size,locant_order,ring_order, bridge_atoms,buffer);
  for(unsigned int i=0;i<pairs;i++)
    assignment_score += 10; // this should mean we always try and avoid them

  free(ring_arr);
  return assignment_score;  
}
#endif



LocantPos *SingleWalk(OBMol *mol, unsigned int path_size,
                    std::set<OBRing*> &local_SSSR,
                    std::vector<OBRing*> &ring_order,
                    std::string &buffer )
{
  LocantPos *locant_path = (LocantPos*)malloc(sizeof(LocantPos) * path_size); 
  zero_locant_path(locant_path, path_size); 

  OBRing *mono = *(local_SSSR.begin());

  for(unsigned int i=0;i<mono->Size();i++){
    locant_path[i].atom = mol->GetAtom(mono->_path[i]);
    locant_path[i].locant = INT_TO_LOCANT(i+1); 
  }
  
  write_ring_size(mono, buffer); 
  ring_order.push_back(mono);
  return locant_path;
}



/* unless one of the atoms is multicyclic, crossing a ring junction is not allowed
 * a ring junction is determined if bond between 2 atoms is shared between rings
 * 
 * For multicyclics is a simple ask, you want the highest ring shares in the lowest position, 
 * therefore if all bonds are junctions, take the multicyclic point.
*/
bool IsRingJunction(OBMol*mol, OBAtom *curr, OBAtom *ahead, std::set<OBRing*>&local_SSSR){
  OBBond * bond = mol->GetBond(curr,ahead); 
  if(!bond){
    fprintf(stderr,"Error: bond does not exist\n"); 
    return false;
  }
    
  unsigned int shares = 0; 
  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end(); riter++){
    OBRing *obring = *riter; 
    if(obring->IsMember(bond)){
      shares++;
      if(shares >= 2)
        return true; 
    }
  }

  return false; 
}

bool IsRingComplete(OBRing *ring, LocantPos *locant_path, unsigned int path_len){
  unsigned int size = 0; 
  for(unsigned int i=0;i<path_len;i++){
    if(locant_path[i].atom && ring->IsMember(locant_path[i].atom))
      size++; 
  }
  
  return (size == ring->Size()); 
}

bool IsRingCompleteWB(OBRing *ring, 
                      LocantPos *locant_path, unsigned int path_len,
                      LocantPos *broken_path, unsigned int b_len){
  unsigned int size = 0; 
  for(unsigned int i=0;i<path_len;i++){
    if(locant_path[i].atom && ring->IsMember(locant_path[i].atom))
      size++; 
  }

  for(unsigned int i=0;i<b_len;i++){
    if(broken_path[i].atom && ring->IsMember(broken_path[i].atom))
      size++; 
  }

  return (size == ring->Size()); 
}

/*  a valid starting multicyclic point cannot be nested within others, it must only point to 
 *  1 other multicyclic, within the ring walk this is not true, and both routes must be taken
 *  */ 
unsigned int connected_multicycles(OBAtom *atom, std::map<OBAtom*,unsigned int>  &atom_shares){
  unsigned int m = 0;
  FOR_NBORS_OF_ATOM(a, atom){
    OBAtom *n = &(*a);
    if(atom_shares[n] >= 3)
      m++; 
  }
  return m;
}


/*
builds iteratively so ordering is correct, when a ring is filled, write the notation
there is some bit logic to speed all this up, - concepts first, optimisation later
Some rules to follow when walking the path,
as long as the sequential order is FORCED, we can take the max array size as 
max_path_size
*/
void write_complete_rings(  OBMol *mol, LocantPos *locant_path, unsigned int max_path_size, 
                            std::set<OBRing*> &local_SSSR, 
                            std::map<OBRing*,bool>  &handled_rings, 
                            std::vector<OBRing*> &ring_order, 
                            std::string &buffer)
{
  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end(); riter++){
    if(!handled_rings[*riter] && IsRingComplete(*riter, locant_path, max_path_size)){
      unsigned int lowest_locant = lowest_ring_locant(mol,*riter, locant_path, max_path_size);
      if(lowest_locant != 'A'){
        buffer += ' '; 
        write_locant(lowest_locant,buffer); 
      }
      
      write_ring_size(*riter, buffer); 
      handled_rings[*riter] = true;
      ring_order.push_back(*riter);
    }
  }
}

/* same as before but allow through a branching locant array to check solves */
void write_complete_ringsWB(  OBMol *mol, LocantPos *locant_path, unsigned int max_path_size, 
                            LocantPos *branching_locants, unsigned int branch_n,
                            std::set<OBRing*> &local_SSSR, std::map<OBRing*,bool> &handled_rings,
                            std::vector<OBRing*> &ring_order, 
                            std::string &buffer)
{

  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end(); riter++){
    if(!handled_rings[*riter] && IsRingCompleteWB(*riter, locant_path, max_path_size,branching_locants,branch_n)){
      unsigned int lowest_locant = lowest_ring_locantWB(mol,*riter, locant_path, max_path_size,branching_locants,branch_n);
      if(lowest_locant != 'A'){
        buffer += ' '; 
        write_locant(lowest_locant,buffer); 
      }
      
      write_ring_size(*riter, buffer); 
      handled_rings[*riter] = true;
      ring_order.push_back(*riter);
    }
  }
}

/*  standard ring walk, can deal with all standard polycyclics without an NP-Hard
    solution, fusion sum is the only filter rule needed here, for optimal branch, 
    create notation as the path is read to prove the concept for removal of NP flood fill. 

1) If pointing to something that is not a ring junction (if the choice is there)
   take the path with the highest atom share count

   * from a given starting point there will be two possible directions for a standard polycyclic
     call these direction A, and direction B, which can be scored by the fusion sum

2) Write the notation as rings loop back to the path, if a previous locant can be seen from the current
   position, the ring must of looped back, write the smallest locant in that subring. 
  
   this includes the final A if looping back to start 

3) In order to do this, keep a track of the last ring entered, for non-sharing atoms, this will update the ring
   * see UpdateCurrentRing 

4) When a stack point is hit, walk the locant path back to where the stack atom is, mark all visited nodes as 
   false in the walk back

3 and 4 are likely not needed for polycyclic, see ComplexWalk for implementation on multicyclics, bridges etc. 
*/
LocantPos *PathFinderIIIa(    OBMol *mol, unsigned int path_size,
                              std::set<OBAtom*>               &ring_atoms,
                              std::map<OBAtom*,unsigned int>  &atom_shares,
                              std::map<OBAtom*,bool>          &bridge_atoms,
                              std::set<OBRing*>               &local_SSSR,
                              std::vector<OBRing*>            &ring_order,
                              std::string                     &buffer)
{

  // create the path
  LocantPos *locant_path = (LocantPos*)malloc(sizeof(LocantPos) * path_size); 
  LocantPos *best_path = (LocantPos*)malloc(sizeof(LocantPos) * path_size); 
    
  zero_locant_path(locant_path, path_size);
  zero_locant_path(best_path, path_size); 

  std::vector<OBRing*> best_order; 
    
  OBAtom*                ratom  = 0; // ring
  OBAtom*                catom  = 0; // child
  OBAtom*                matom  = 0; // move atom
  unsigned int           lowest_sum       = UINT32_MAX;
  
  for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
    if(atom_shares[*aiter] == 2){ // these are the starting points 
     
      zero_locant_path(locant_path, path_size);

      std::set<OBBond*>      ring_junctions;
      std::map<OBAtom*,bool> visited; 
      unsigned int locant_pos = 0;
      
      ratom = *aiter; 
      for(;;){

        locant_path[locant_pos].atom = ratom; 
        locant_path[locant_pos].locant = INT_TO_LOCANT(locant_pos+1);        
        locant_pos++; 
        visited[ratom] = true;

        if(locant_pos >= path_size)
          break;
        
        matom = 0;  
        FOR_NBORS_OF_ATOM(a,ratom){ 
          catom = &(*a);   
          if(!visited[catom] && !IsRingJunction(mol, ratom, catom,local_SSSR)){
            if(!matom)
              matom = catom; 
            else if(atom_shares[catom] > atom_shares[matom])
              matom = catom; 
          }
        }
        
        if(!matom){
          fprintf(stderr,"Error: did not move in locant path walk!\n"); 
          return 0;
        }

        ratom = matom; 
      }
     
      unsigned int fsum = fusion_sum(mol,locant_path,path_size,local_SSSR);
      if(fsum < lowest_sum){ // rule 30d.
        lowest_sum = fsum;
        copy_locant_path(best_path,locant_path,path_size);
      }
    }
  }

  free(locant_path);
  

  std::map<OBRing*,bool> handled_rings; 
  for(unsigned int i=0;i<=path_size;i++){ // inner function are less than i, therefore <=
    write_complete_rings( mol,best_path, i, local_SSSR,handled_rings, 
                          ring_order,buffer); 
  }
  return best_path; 
}


void BackTrackWalk(  OBAtom *clear, LocantPos*locant_path, unsigned int path_size,unsigned int &locant_pos, std::map<OBAtom*,bool> &visited_atoms)
{

  // find the position in the path where the lowest one of these are, the other gets placed into locant path 
  unsigned int p=0;
  unsigned int q=0;
  for(p=0;p<path_size;p++,q++){
    if(locant_path[p].atom == clear)
      break;
  }

  // everything from p gets cleared, but not including
  for(++p;p<path_size;p++){
    visited_atoms[locant_path[p].atom] = 0;
    locant_path[p].atom   = 0;
    locant_path[p].locant = 0; 
  }
  locant_pos = q+1;
}


/*
Some rules to follow when walking the path:

1) If pointing to something that is not a ring junction (if the choice is there)
   take the path with the highest atom share count

   * from a given starting point there will be two possible directions for a standard polycyclic
     call these direction A, and direction B, which can be scored by the fusion sum

2) Write the notation as rings loop back to the path, if a previous locant can be seen from the current
   position, the ring must of looped back, write the smallest locant in that subring. 
  
   this includes the final A if looping back to start 

3) In order to do this, keep a track of the last ring entered, for non-sharing atoms, this will update the ring
   * see UpdateCurrentRing 

4) When a stack point is hit, walk the locant path back to where the stack atom is, mark all visited nodes as 
   false in the walk back

  - this significantly speeds up and removes an NP-hard algorithm 


  path size can now change to accomadate whether the locant path reduces due to branching locants
*/
LocantPos *PathFinderIIIb(  OBMol *mol,        unsigned int &path_size,
                            std::set<OBAtom*>               &ring_atoms,
                            std::map<OBAtom*,unsigned int>  &atom_shares,
                            std::map<OBAtom*,bool>          &bridge_atoms,
                            std::set<OBRing*>               &local_SSSR,
                            std::vector<OBRing*>            &ring_order,
                            std::string                     &buffer)
{

  unsigned int locant_pos   = 0; 
  unsigned int off_branch_n = 0; 
  
  LocantPos *locant_path = (LocantPos*)malloc(sizeof(LocantPos) * path_size); 
  LocantPos *best_path = (LocantPos*)malloc(sizeof(LocantPos) * path_size); 

  LocantPos *off_branches = (LocantPos*)malloc(sizeof(LocantPos) * 32); // hard limit 
  LocantPos *best_off_branches = (LocantPos*)malloc(sizeof(LocantPos) * 32); // hard limit 

  zero_locant_path(locant_path, path_size);
  zero_locant_path(best_path, path_size);

  zero_locant_path(off_branches, 32); 
  zero_locant_path(best_off_branches, 32); 

  OBAtom*                ratom  = 0; // ring
  OBAtom*                catom  = 0; // child
  OBAtom*                matom  = 0; // move atom
  unsigned int           lowest_sum         = UINT32_MAX;
  unsigned int           starting_path_size = path_size; // important if path size changes
  unsigned int           best_path_size     = 0; 
  unsigned int           best_off_branch_n  = 0; 

  for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
    // a multicyclic that connects to two other multicyclic points can never be the start, always take an edge case
    if( (atom_shares[*aiter] >= 3 && connected_multicycles(*aiter,atom_shares)<=1)  || bridge_atoms[*aiter]){ // these are the starting points 

      std::map<OBAtom*,bool>    visited; 
      std::stack<OBAtom*>       multistack; 
      std::stack<std::pair<OBAtom*,OBAtom*>> backtrack_stack;   // multicyclics have three potential routes, 
      
      locant_pos = 0;
      off_branch_n = 0; 
      path_size = starting_path_size; 
      zero_locant_path(locant_path, starting_path_size);
      zero_locant_path(off_branches, 32);
    
path_solve:        
      ratom = *aiter; 
      locant_pos = 0;
      do{

        if(!backtrack_stack.empty())
          backtrack_stack.pop(); 

        for(;;){
          
          locant_path[locant_pos].atom = ratom;
          locant_path[locant_pos].locant = INT_TO_LOCANT(locant_pos+1);
          locant_pos++; 
          
          // check if attached to broken, if yes, update their locants
          update_broken_locants(mol, &locant_path[locant_pos-1], locant_path, path_size, off_branches, off_branch_n); 
          visited[ratom] = true;

          if(locant_pos >= path_size)
            break;
          
          // here we allow multicyclics to cross ring junctions
          matom = 0; 
          FOR_NBORS_OF_ATOM(a,ratom){ 
            catom = &(*a);  
            if(!visited[catom]){
              // two things can happen, either we're at a ring junction or we're not
              if(IsRingJunction(mol, ratom, catom, local_SSSR)){
                // if its a ring junction, we can move if this is going to/from a multicyclic point,
                // if pointing at a multicyclic, or an edge atoms, try both
                if( (atom_shares[ratom]>=3 || atom_shares[catom]>=3) || (bridge_atoms[ratom] || bridge_atoms[catom]) ){
                  if(!matom)
                    matom = catom; 
                  else if(atom_shares[catom] > atom_shares[matom]){
                    // edge cases should be tried
                    if(atom_shares[matom] < 2)
                      backtrack_stack.push({ratom,matom}); 

                    matom = catom;
                  }
                  else 
                    backtrack_stack.push({ratom,catom}); 
                }
              }
              else if(atom_shares[catom] < 3){
                if(!matom)
                  matom = catom; 
                else if(atom_shares[catom] > atom_shares[matom]){
                  if(atom_shares[matom] < 2)
                    backtrack_stack.push({ratom,matom}); 
                  
                  matom = catom;
                }
                else if(atom_shares[catom] <= atom_shares[matom])
                  backtrack_stack.push({ratom,catom});
              }
              else if(atom_shares[ratom] < 2){
                if(!matom)
                  matom = catom;
              }
            }
          }
          
          if(!matom){
            // no locant path! add logic for off branches here! and pseudo locants
            //fprintf(stderr,"Error: did not move in locant path walk!\n"); 
            break;
          }
            
          ratom = matom;
          if(atom_shares[ratom ]>= 3 || bridge_atoms[ratom]) // ignores the first one
            multistack.push(ratom); 
        }
       
        if(locant_pos == path_size){
          unsigned int fsum = fusion_sum(mol,locant_path,path_size,local_SSSR);
          if(fsum < lowest_sum){ // rule 30d.
            lowest_sum = fsum;
            copy_locant_path(best_path,locant_path,starting_path_size);
            copy_locant_path(best_off_branches,off_branches,32);
            best_path_size = path_size;
            best_off_branch_n = off_branch_n; 
          }
        }
        
        if(!backtrack_stack.empty()){
          ratom = backtrack_stack.top().second;
          BackTrackWalk(backtrack_stack.top().first, locant_path, path_size,locant_pos,visited); 
        }
        else if(!best_path[0].atom && !multistack.empty()){ // once you find a branch path, take it!
          // this the where the broken locants happen, pop off a multistack atom 
          OBAtom *branch_locant = multistack.top();
          multistack.pop();
          for(unsigned int p=0;p<starting_path_size;p++)
            visited[locant_path[p].atom] = 0;
          
          burn_stack(backtrack_stack); 
          burn_stack(multistack);  
          visited[branch_locant] = true;
          
          off_branches[off_branch_n].atom   = branch_locant; 
          off_branches[off_branch_n].locant = 0; 
          off_branch_n++;
          path_size--; // decrement the path size, this is globally changed
          
          //fprintf(stderr,"removing atom: %d\n",branch_locant->GetIdx()); 
          zero_locant_path(locant_path, starting_path_size); 
          goto path_solve; 
        }
        else
          break;
        
      } while(!backtrack_stack.empty()) ; 
    }
  } 

  free(locant_path);
  free(off_branches); 

  std::map<OBRing*,bool> handled_rings; 
  for(unsigned int i=0;i<=path_size;i++){ // inner function are less than i, therefore <=
    write_complete_ringsWB( mol,best_path, i, 
                            best_off_branches,best_off_branch_n,
                            local_SSSR,handled_rings, 
                            ring_order,buffer); 
  }

  // add the branching locants at the back of the locant_array for indexing
  for(unsigned int i=0;i<best_off_branch_n;i++){
    best_path[best_path_size+i].atom = best_off_branches[i].atom; 
    best_path[best_path_size+i].locant = best_off_branches[i].locant; 
  }

  TestPathSequences(mol, best_path, starting_path_size,ring_order,atom_shares,bridge_atoms,buffer); 

  free(best_off_branches);
  path_size = best_path_size; 
  return best_path; 
}


#if DEPRECATED
  /* uses a flood fill style solution (likely NP-HARD), with some restrictions to 
  find a multicyclic path thats stable with disjoined pericyclic points */
OBAtom **PeriWalk(      OBMol *mol, unsigned int path_size,
                        std::set<OBAtom*>               &ring_atoms,
                        std::set<OBBond*>               &ring_bonds,
                        std::map<OBAtom*,unsigned int>  &atom_shares,
                        std::map<OBAtom*,bool>          &bridge_atoms, // rule 30f.
                        std::set<OBRing*>               &local_SSSR,
                        unsigned int recursion_tracker)
  {

    // parameters needed to seperate out the best locant path
    unsigned int           lowest_fsum       = UINT32_MAX;
    unsigned int           lowest_score      = UINT32_MAX;

    // multi atoms are the starting seeds, must check them all unfortuanately 
    std::vector<OBAtom*> seeds; 
    for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
      OBAtom *rseed = (*aiter);
      if(atom_shares[rseed] >= 1)
        seeds.push_back(rseed);
    }

    bool path_found = false;
    unsigned int found_path_size = path_size;

    OBAtom **locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * found_path_size); 
    OBAtom **best_path = (OBAtom**)malloc(sizeof(OBAtom*) * found_path_size); 
    for(unsigned int i=0;i<found_path_size;i++){
      locant_path[i] = 0;
      best_path[i] = 0; 
    }

    for(OBAtom *rseed : seeds){
      OBAtom*                catom  = 0;
      std::map<OBAtom*,bool> current; 
      std::vector<std::pair<OBAtom*,OBAtom*>> path; 
      path.push_back({rseed,0}); 

      unsigned int safety = 0;
      while(!path.empty()){ // need a loop safety
        safety++;
        OBAtom* ratom = path.back().first;
        OBAtom* next = path.back().second;

        current[ratom] = true;
        
        bool skipped = false;
        bool pushed = false;

        if(!next)
          skipped = true; 

        FOR_NBORS_OF_ATOM(a,ratom){ // this relies on this being a deterministic loop
          catom = &(*a);
  //        bool in_set = true;

          bool in_set = false; 
          for (std::set<OBAtom*>::iterator siter = ring_atoms.begin();siter != ring_atoms.end(); siter++) {
            if (catom == *siter)
              in_set = true;
          }

          if(in_set && atom_shares[catom]){
            if(catom == next)
              skipped = true;
            else if(!current[catom] && skipped && !pushed){
              path.push_back({catom,0});
              pushed = true; 
              break;
            }
          }
        }

        if(!pushed && !path.empty()){
          if(path.size() == found_path_size){
            path_found = true;
            for(unsigned int i=0;i<found_path_size;i++)
              locant_path[i] = path[i].first;
            
            std::vector<OBRing*> tmp; 
            std::string candidate_string; // unfortunate but necessary, very expensive procedure (we optimise later!)
            unsigned int score = ReadLocantPath(mol,locant_path,found_path_size,local_SSSR,bridge_atoms,tmp,candidate_string,false);
            unsigned int fsum = fusion_sum(mol,locant_path,found_path_size,local_SSSR);


#if WLNDEBUG
        fprintf(stderr, "%s - score: %d, fusion sum: %d\n", candidate_string.c_str(),score,fsum);       
#endif 

            if(score < lowest_score){ // rule 30(d and e).
              lowest_score = score;
              lowest_fsum = fsum;
              copy_locant_path(best_path,locant_path,found_path_size);
            }
            else if (score == lowest_score){
              if(fsum < lowest_fsum){
                lowest_fsum = fsum; 
                copy_locant_path(best_path,locant_path,found_path_size);
              }
            }

          }

          OBAtom *tmp = path.back().first; 
          path.pop_back();
          if(!path.empty()){
            path.back().second = tmp;
            current[tmp] = false; 
          }
        }

        // // super defensive temporary measure, this CANNOT be in the release
        // // guards agaisnt C60
        if(safety == 100000)
          break;
        
          
      }
    }

    for(unsigned int i=0;i<path_size;i++)
      locant_path[i] = 0;  
    free(locant_path);

    // call recursion here for path finding
    if(!path_found){
      free(best_path);
      best_path = 0; 

      // try and remove one ring and stop the recursion
      if(recursion_tracker == 0){
        
        unsigned int pos = 0;
        for(std::set<OBRing*>::iterator riter = local_SSSR.begin();riter != local_SSSR.end();riter++){
          OBRing *obring = *riter; 
          std::set<OBAtom*> local_atoms; 
          std::set<OBAtom*> difference;
          
          // its the difference ONLY if the atoms are ONLY contained in this ring
          for(unsigned int i=0;i<obring->Size();i++){
            OBAtom *latom = mol->GetAtom(obring->_path[i]);    
            if(atom_shares[latom] == 1)
              local_atoms.insert(latom);
          }
          
          if(!local_atoms.empty()){
            std::set_difference(ring_atoms.begin(), ring_atoms.end(), local_atoms.begin(), local_atoms.end(),
                                std::inserter(difference, difference.begin()));
            
            best_path = PeriWalk(mol, difference.size(), difference, ring_bonds,atom_shares, bridge_atoms, local_SSSR, 1); 
            if(best_path){
              // remove ring from the local SSSR, mark all the local_atoms set as non cyclic
              // and non-aromatic!
              for(std::set<OBAtom*>::iterator laiter=local_atoms.begin(); laiter != local_atoms.end();laiter++){
                (*laiter)->SetInRing(false);
                bridge_atoms[*laiter] = false;
              }
            
              // bridges can only be made on two ring intersections? 
              for(unsigned int i=0;i<obring->Size();i++){
                OBAtom *latom = mol->GetAtom(obring->_path[i]);    
                if(bridge_atoms[latom])
                  bridge_atoms[latom] = false;
              }
        
              std::set<OBRing*>::iterator it = std::next(local_SSSR.begin(), pos); 
              local_SSSR.erase(it);
              
              // remove any aromaticty contraints so bonds get written
              std::set<OBBond*> local_ring_bonds; 
              FillRingBonds(mol, obring, local_ring_bonds); 
              for(std::set<OBBond*>::iterator biter=local_ring_bonds.begin(); biter != local_ring_bonds.end();biter++){
                // these need to erased from the ring
                (*biter)->SetAromatic(false); 
              std::set<OBBond*>::iterator gpos = std::find(ring_bonds.begin(),ring_bonds.end(), *biter);   
              ring_bonds.erase(gpos); 
            }

            ring_atoms = difference;  // set the ring atoms to what we expect 
            return best_path;
          }

          pos++;
        }
      }      
    }  

  }
  
  return best_path;
}
#endif

/**********************************************************************
                         Debugging Functions
**********************************************************************/


void BabelDumptoDot(FILE *fp, OBMol *mol){

  fprintf(fp, "digraph BABELdigraph {\n");
  fprintf(fp, "  rankdir = LR;\n");
  FOR_ATOMS_OF_MOL(aiter,mol){
    OBAtom* atom = &(*aiter); 
    fprintf(fp, "  %d", atom->GetIdx());
    fprintf(fp, "[shape=circle,label=\"%s\"];\n", std::to_string(atom->GetIdx()).c_str());
  }

  FOR_BONDS_OF_MOL(biter,mol){
    OBBond* bond = &(*biter); 
    fprintf(fp, "  %d", bond->GetBeginAtom()->GetIdx() );
    fprintf(fp, " -> ");
    fprintf(fp, "%d\n", bond->GetEndAtom()->GetIdx() );

  }

  fprintf(fp, "}\n");
}


bool WriteBabelDotGraph(OBMol *mol){
  fprintf(stderr,"Dumping babel graph to babel-graph.dot:\n");
  FILE *fp = 0;
  fp = fopen("babel-graph.dot", "w");
  if (!fp)
  {
    fprintf(stderr, "Error: could not create dump .dot file\n");
    fclose(fp);
    return false;
  }
  else
    BabelDumptoDot(fp,mol);
  
  fclose(fp);
  fp = 0;
  fprintf(stderr,"  dumped\n");
  return true;
}


/**********************************************************************
                         BABEL Mol Functions
**********************************************************************/


// holds all the functions for WLN graph conversion, mol object is assumed ALIVE AT ALL TIMES
// uses old NM functions from previous methods: Copyright (C) NextMove Software 2019-present
struct BabelGraph{
  
  bool modern; 
  std::map<OBAtom*,bool> atoms_seen;
  std::map<OBAtom*,unsigned char> atom_chars;
  std::map<OBRing*,bool> rings_seen; 
  std::map<OBAtom*,int>  remaining_branches; // tracking for branch pop
  std::map<OBAtom*,unsigned int> string_position; // essential for writing post charges. 

  BabelGraph(){
    modern = 0; 
  };
  ~BabelGraph(){};
  

  // if modern, charges are completely independent apart from assumed K
  unsigned char WriteSingleChar(OBAtom* atom){

    if(!atom)
      Fatal("writing notation from dead atom ptr");
    
    unsigned int neighbours = atom->GetExplicitDegree(); 
    unsigned int orders = atom->GetExplicitValence(); 

    switch(atom->GetAtomicNum()){
      case 1:
        return 'H';

      case 5:
        if(orders > 3)
          return '*';
        else
          return 'B';

      case 6:
        if(neighbours <= 2)
          return '1';
        else if(neighbours == 3)
          return 'Y';
        else if(neighbours == 4)
          return 'X';
      
      case 7:
          if(orders <= 1 && atom->GetImplicitHCount() == 2){
            return 'Z'; 
          }
          else if(orders == 2 && atom->GetImplicitHCount() == 1)
              return 'M';
          else if(atom->GetFormalCharge() == +1 && orders == 4)
            return 'K';
          else if (orders >= 4)
            return '*';
          else 
            return 'N';
      
      case 8:
        if(neighbours == 1 && orders==1 && atom->GetFormalCharge() == 0)
          return 'Q';
        else if (neighbours == 0 && atom->GetFormalCharge() != -2)
          return 'Q';
        else if (atom->GetExplicitValence() > 2)
          return '*';
        else
          return 'O';
          
      case 9:
        if(neighbours > 1)
          return '*';
        else 
          return 'F';

      case 15:
        if(neighbours > 5)
          return '*'; 
        return 'P';

      case 16:
        return 'S';
  
      case 17:
        if(neighbours > 1)
          return '*';
        else
          return 'G';

      case 35:
        if(neighbours > 1)
          return '*';
        else
          return 'E';

      case 53:
        if(neighbours > 1)
          return '*';
        else
          return 'I';

      default:
        return '*';
        
    }

    return 0; 
  }
  
  void ModernCharge(OBAtom *atom, std::string &buffer){
    if(atom->GetFormalCharge() == 0)
      return; 
    
    if(abs(atom->GetFormalCharge())>1)
      buffer += std::to_string(abs(atom->GetFormalCharge()));  
    
    if(atom->GetFormalCharge() < 0)
      buffer += '-';
    else
     buffer += '+';
    return; 
  }



  bool CheckCarbonyl(OBAtom *atom){
    if(!atom)
      Fatal("checking for carbonyl on dead atom ptr");

    if(atom->GetAtomicNum() != 6)
      return false;

    FOR_NBORS_OF_ATOM(a,atom){
      OBAtom *nbor = &(*a);
      if(!atoms_seen[nbor] && !nbor->IsInRing() && nbor->GetAtomicNum() == 8){
        if(atom->GetBond(nbor)->GetBondOrder() == 2){
          atoms_seen[nbor] = true;
          return true;
        }
      }  
    }
    return false;
  }

  OBAtom *return_open_branch(std::stack<OBAtom*> &branch_stack){
    while(!branch_stack.empty()){
      if(remaining_branches[branch_stack.top()] > 0)
        return branch_stack.top();
      else
        branch_stack.pop();
    }
    return 0;
  }


  void AddPostCharges(OBMol *mol,std::string &buffer){
    
    bool working = true;
    while(working){
      working = false;
      FOR_ATOMS_OF_MOL(a,mol){
        OBAtom *atom = &(*a);
        if(atom->GetFormalCharge() != 0){


          if(atom->GetFormalCharge() > 0){
            buffer += ' ';
            buffer += '&';
            buffer += std::to_string(string_position[atom]);
            buffer += '/';
            buffer += '0';
            atom->SetFormalCharge(atom->GetFormalCharge()-1);
            working = true;
          }

          if(atom->GetFormalCharge() < 0){
            buffer += ' ';
            buffer += '&';
            buffer += '0';
            buffer += '/';
            buffer += std::to_string(string_position[atom]);
            atom->SetFormalCharge(atom->GetFormalCharge()+1);
            working = true;
          }
        }

      }
     
    }
  }

  /* parses the local ring system, return the size for creating the locant path with 
  non bonds to avoid */
  bool  ConstructLocalSSSR( OBMol *mol, OBAtom *ring_root,
                            std::set<OBAtom*>         &ring_atoms,
                            std::set<OBBond*>         &ring_bonds,
                            std::map<OBAtom*,bool>    &bridge_atoms,
                            std::map<OBAtom*,unsigned int> &atom_shares,
                            std::map<OBBond*,unsigned int> &bond_shares, 
                            std::set<OBRing*> &local_SSSR,
                            SubsetData &local_data)
  {

    if(!ring_root){
      fprintf(stderr,"Error: ring root is nullptr\n");
      return false;
    }

    OBAtom *ratom = 0; 
    OBAtom *prev = 0; 
    OBBond *bond = 0; 
    OBRing *obring = 0; 
    std::set<OBAtom*> tmp_bridging_atoms;
    std::set<OBAtom*> remove_atoms; 

    // get the seed ring and add path to ring_atoms
    FOR_RINGS_OF_MOL(r,mol){
      obring = &(*r);
      if(obring->IsMember(ring_root)){
        rings_seen[obring] = true;
        local_SSSR.insert(obring);
 
        prev = 0; 
        for(unsigned int i=0;i<obring->Size();i++){
          ratom = mol->GetAtom(obring->_path[i]);
          atom_shares[ratom]++;
          ring_atoms.insert(ratom); 
          if(ratom->GetAtomicNum() != 6)
            local_data.hetero = true;

          if(!prev)
            prev = ratom; 
          else{
            bond = mol->GetBond(prev,ratom);
            if(bond){
              ring_bonds.insert(bond); 
              bond_shares[bond]++; 
            }
            prev = ratom; 
          }
        }
        // get the last bond
        bond = mol->GetBond(mol->GetAtom(obring->_path.front()),mol->GetAtom(obring->_path.back()));
        if(bond){
          ring_bonds.insert(bond); 
          bond_shares[bond]++; 
        }
        break;
      }
    }
    
    bool running = true; 
    while(running){
      running = false;

      FOR_RINGS_OF_MOL(r,mol){
        obring = &(*r);

        if(!rings_seen[obring]){

          std::set<OBAtom*> ring_set; 
          std::set<OBAtom*> intersection; 
          bool all_ring = true;

          for(unsigned int i=0;i<obring->Size();i++){
            OBAtom *ratom = mol->GetAtom(obring->_path[i]);
            ring_set.insert(ratom);
            if(!ratom->IsInRing())
              all_ring = false;
          }

          std::set_intersection(ring_set.begin(), ring_set.end(), ring_atoms.begin(), ring_atoms.end(),
                                std::inserter(intersection, intersection.begin()));

          // intersection == 1 is a spiro ring, ignore,
          if(intersection.size() > 1 && all_ring){

            prev = 0;
            // if its enough to say that true bridges cannot have more than two bonds each?
            // yes but 2 bonds within the completed local SSSR,so this will needed filtering
            if(intersection.size() > 2){
              for(std::set<OBAtom*>::iterator iiter = intersection.begin(); iiter != intersection.end();iiter++){
                tmp_bridging_atoms.insert(*iiter);
              }
            }

            for(unsigned int i=0;i<obring->Size();i++){
              OBAtom *ratom = mol->GetAtom(obring->_path[i]);
              ring_atoms.insert(ratom);
              atom_shares[ratom]++;
              if(atom_shares[ratom] >= 3)
                local_data.multi = true; // set the multi classification
              if(ratom->GetAtomicNum() != 6)
                local_data.hetero = true;
                
              if(!prev)
                prev = ratom;
              else{
                bond = mol->GetBond(prev,ratom);
                if(bond){
                  ring_bonds.insert(bond); 
                  bond_shares[bond]++; 
                }
                prev = ratom; 
              }
            }
            // get the last bond
            bond = mol->GetBond(mol->GetAtom(obring->_path.front()),mol->GetAtom(obring->_path.back()));
            if(bond){
              ring_bonds.insert(bond); 
              bond_shares[bond]++; 
            }
            rings_seen[obring] = true; 
            local_SSSR.insert(obring);
            running = true;
          }
        }
      }
    }

    // filter out only the 2 bond bridge atoms
    unsigned int bridge_count = 0;
    if(!tmp_bridging_atoms.empty()){
      for(std::set<OBAtom*>::iterator brd_iter=tmp_bridging_atoms.begin(); brd_iter != tmp_bridging_atoms.end();brd_iter++){
        unsigned int inter_ring_bonds = 0;
        for(std::set<OBAtom*>::iterator aiter= ring_atoms.begin(); aiter != ring_atoms.end();aiter++){
          if(mol->GetBond(*brd_iter,*aiter))
            inter_ring_bonds++; 
        }
        if(inter_ring_bonds == 2){
          bridge_count++;
          bridge_atoms[*brd_iter] = true;
        }
      }
    }




    local_data.bridging = bridge_count;
    local_data.path_size = ring_atoms.size();  
    return true; 
  }

#if 0
  /* create the heteroatoms and locant path unsaturations where neccesary */
  bool ReadLocantAtomsBonds(  OBMol *mol, LocantPos* locant_path,unsigned int path_size,
                              std::vector<OBRing*> &ring_order,
                              std::set<OBBond*>   &ring_bonds,
                              std::string &buffer)
  {

    OBAtom *locant_atom = 0; 
    unsigned int locant_char = 0;
    unsigned int last_locant = 'A';
    std::map<OBBond*,bool> bonds_checked; 

    if(!std::isdigit(buffer.back()))
      last_locant = ' ';
    
    for(unsigned int i=0;i<path_size;i++){
      
      unsigned char het_char = 0;
      locant_atom = locant_path[i].atom; 
      locant_char = locant_path[i].locant;
      
      if(!locant_path[i].atom){
        print_locant_array(locant_path, path_size); 
        Fatal("dead locant path atom ptr in hetero read - atom");
      }

      if(!locant_path[i].locant)
        Fatal("dead locant path position in hetero read - locant");
      
      bool carbonyl = CheckCarbonyl(locant_atom);

      if( !carbonyl &&  
        locant_atom->GetAtomicNum() == 6 &&
        locant_atom->GetFormalCharge() == -1){
        // organometallics logic 
        if(locant_char != last_locant){
          buffer += ' ';
          write_locant(locant_char,buffer);
          last_locant = locant_char;
        } 

        buffer += '0';
        locant_atom->SetFormalCharge(0);
      }
      
      if(locant_atom->GetAtomicNum() == 6)
        atom_chars[locant_atom] = '1'; 

      if(carbonyl || locant_atom->GetAtomicNum() != 6){
        if(locant_char != last_locant){
          buffer += ' ';
          write_locant(locant_char,buffer);
          last_locant = locant_char;
        } 
        
        if(carbonyl){
          buffer += 'V';
          string_position[locant_atom] = buffer.size();
          last_locant++;
        }
        else{
          het_char = WriteSingleChar(locant_atom);
          atom_chars[locant_atom] = het_char; 
          if(het_char != '*'){
            if(het_char == 'K')
              locant_atom->SetFormalCharge(0);

            buffer+=het_char;
            string_position[locant_atom] = buffer.size();
          }
          else{
            WriteSpecial(locant_atom,buffer); 
          }

          last_locant++;
        }
      }

      if(locant_atom->GetAtomicNum() == 6){
        unsigned int rbonds = 0; 
        for(unsigned int k=0;k<path_size;k++){
          if(mol->GetBond(locant_atom,locant_path[k].atom))
            rbonds++;
        }
        if(rbonds == 4){
          if(locant_char != last_locant){
            buffer += ' ';
            write_locant(locant_char,buffer);
            last_locant = locant_char;
          } 
          buffer += 'X';
        }
      }

      // for now dont worry about positions, only that we get the double bonds
      // we can condense later
      OBAtom *first   = 0;
      OBAtom *second  = 0;
      // handles sequential locant unsaturations, when not aromatic
      first = locant_path[i].atom;
      if(i < path_size-1)
        second = locant_path[i+1].atom;
      else
        second = locant_path[0].atom;

      OBBond *locant_bond = first->GetBond(second);
      bonds_checked[locant_bond] = true;

      if(locant_bond && !locant_bond->IsAromatic() && locant_bond->GetBondOrder()>1){
        buffer += ' ';
        write_locant(locant_char,buffer);

        for(unsigned int b=1;b<locant_bond->GetBondOrder();b++)
          buffer += 'U';
      }
    }


    for(std::set<OBBond*>::iterator biter = ring_bonds.begin(); biter != ring_bonds.end();biter++){
      OBBond *fbond = *biter; 
      
      if(!bonds_checked[fbond] && fbond->GetBondOrder() > 1 && !fbond->IsAromatic()){
        unsigned int floc = locant_path[position_in_path(fbond->GetBeginAtom(),locant_path,path_size)].locant; 
        unsigned int bloc = locant_path[position_in_path(fbond->GetEndAtom(),locant_path,path_size)].locant; 
        
        buffer += ' ';
        write_locant(floc,buffer);
        for(unsigned int b=1;b<fbond->GetBondOrder();b++)
          buffer += 'U';
        buffer+='-';
        buffer+=' ';
        write_locant(bloc,buffer);
        break;
      }
    }
    
    return true;
  }

  void ReadMultiCyclicPoints( LocantPos *locant_path,unsigned int path_size, 
                              std::map<OBAtom*,unsigned int> &ring_shares,std::string &buffer)
  { 

    unsigned int count = 0;
    std::string append = ""; 
    for(unsigned int i=0;i<path_size;i++){
      if(ring_shares[locant_path[i].atom] > 2){
        count++; 
        write_locant(locant_path[i].locant,append);
      }
    }

    buffer += ' ';
    buffer+= std::to_string(count);
    buffer+= append;
  }

  /* constructs and parses a cyclic structure, locant path is returned with its path_size */
  void ParseCyclic(OBMol *mol, OBAtom *ring_root,OBAtom *spawned_from,bool inline_ring,PathData &pd,std::string &buffer){

    LocantPos*                      locant_path = 0; 
    std::set<OBRing*>               local_SSSR;
    std::set<OBAtom*>               ring_atoms;
    std::set<OBBond*>               ring_bonds;
    std::vector<OBRing*>            ring_order; 

    std::map<OBAtom*,bool>          bridge_atoms;
    std::map<OBAtom*,unsigned int>  atom_shares;
    std::map<OBBond*,unsigned int>  bond_shares;
    
    // can we get the local SSSR data out of this?
    SubsetData LocalSSRS_data;

    if(!ConstructLocalSSSR(mol,ring_root,ring_atoms,ring_bonds,bridge_atoms,atom_shares,bond_shares,local_SSSR, LocalSSRS_data)){
      Fatal("failed to contruct SSSR"); 
    }
    
    bool multi = LocalSSRS_data.multi;
    bool hetero = LocalSSRS_data.hetero;
    bool bridging = LocalSSRS_data.bridging;
    unsigned int path_size = LocalSSRS_data.path_size;
    bool macro_ring = false; 
    unsigned int branch_locants = 0; 
    std::string ring_segment; 

    if(local_SSSR.size() == 1)
      locant_path = SingleWalk(mol,path_size,local_SSSR,ring_order, ring_segment);
    else if(!multi && !bridging)
      locant_path = PathFinderIIIa(mol,path_size,ring_atoms,atom_shares,bridge_atoms,local_SSSR,ring_order,ring_segment);
    else
      locant_path = PathFinderIIIb(mol,path_size, ring_atoms, atom_shares, bridge_atoms, local_SSSR,ring_order,ring_segment); 
    if(!locant_path)
      return Fatal("no locant path could be determined");
    


    branch_locants = LocalSSRS_data.path_size - path_size; 
    

    if(inline_ring){
      buffer+= '-';
      bool spiro = false;
      unsigned char root_locant = 0;
      for(unsigned int i=0;i<path_size;i++){
        if(locant_path[i].atom == ring_root)
          root_locant = INT_TO_LOCANT(i+1);

        if(locant_path[i].atom == spawned_from){
          // must be spiro ring
          root_locant = INT_TO_LOCANT(i+1);
          spiro = true;
          break;
        }
      }

      if(spiro)
        buffer += '&';
      
      buffer += ' ';
      buffer += root_locant;
    }
    
    if(macro_ring)
      buffer += "T-"; 
 
    if(hetero)
      buffer += 'T';
    else
      buffer += 'L';
    

    buffer += ring_segment; 
    
    if(bridging){
      for(unsigned int i=0;i<LocalSSRS_data.path_size;i++){
        if(bridge_atoms[locant_path[i].atom]){
          buffer+= ' ';
          write_locant(locant_path[i].locant,buffer);
        }
      }
    }

    if(multi){
      ReadMultiCyclicPoints(locant_path,LocalSSRS_data.path_size,atom_shares,buffer);
      
      buffer += ' ';
      write_locant(INT_TO_LOCANT(path_size),buffer); // need to make the relative size
    }

    path_size = LocalSSRS_data.path_size; 
    ReadLocantAtomsBonds(mol,locant_path,path_size,ring_order,ring_bonds,buffer);

    // breaks incremented locant notation
    if(buffer.back() == '&')
      buffer+='-';
    
    unsigned int arom_state = 0; // 0 - null, 1 - all aromatic, 2- no aromatic, 3 mixed.
    for(unsigned int i=0;i<ring_order.size();i++){
      unsigned int arom = 0;
      if(ring_order[i]->IsAromatic())
        arom = 1;
      else
        arom = 2; 

      if(arom_state != 0){
        if(arom_state != arom){
          arom_state = 3; 
          break;
        }
      }
      
      arom_state = arom; 
    }

    if(arom_state == 2)
      buffer += 'T';
    else if(arom_state == 3){
      bool space_added = false;
      for(unsigned int i=0;i<ring_order.size();i++){
        if(ring_order[i]->IsAromatic()){
          if(!space_added && (buffer.back() >= 'A' && buffer.back() <= 'Z'))
            buffer+=' ';
      
          buffer+='&';
        }
        else
          buffer += 'T';

        space_added = true;
      }
    }

    buffer += 'J';
    
    pd.locant_path = locant_path;
    pd.path_size = path_size;
    pd.macro_ring = macro_ring; // could be zero, need check
  }
    


};


#if 0
bool 
parse_acyclic(char buffer, unsigned int i, OBMol *mol, OBAtom *atom) {
  
  PathData pd; 
  ParseCyclic(mol,atom,spawned_from,inline_ring,pd,buffer);
  if(!pd.locant_path){
    fprintf(stderr,"Error: failed on cyclic parse\n");
    return false;
  }
  
  // handle hydrogens here
  for(unsigned int i=0;i<pd.path_size;i++){
    atoms_seen[pd.locant_path[i].atom] = true;
    OBAtom *lc = pd.locant_path[i].atom;
    unsigned int locant = pd.locant_path[i].locant;
    
    switch(atom_chars[lc]){
      case '1':
        break; // ring carbons

      case 'M':
        for(unsigned int h=1;h<lc->GetImplicitHCount();h++){
          buffer += " ";
          write_locant(locant, buffer); 
          buffer += 'H';
        }
        break;


      case 'Z':
        for(unsigned int h=2;h<lc->GetImplicitHCount();h++){
          buffer += " ";
          write_locant(locant, buffer); 
          buffer += 'H';
        }
        break;

      case 'P':
        if(lc->GetExplicitValence() & 1)
          break;
        else{
          for(unsigned int h=0;h<lc->GetImplicitHCount();h++){
            buffer += " ";
            write_locant(locant, buffer); 
            buffer += 'H';
          }
        }
        break;

      case 'S':
        if( !(lc->GetExplicitValence() & 1))
          break;
        else{
          for(unsigned int h=0;h<lc->GetImplicitHCount();h++){
            buffer += " ";
            write_locant(locant, buffer); 
            buffer += 'H';
          }
        }
        break;

      default:
        for(unsigned int h=0;h<lc->GetImplicitHCount();h++){
          buffer += " ";
          write_locant(locant, buffer); 
          buffer += 'H';
        }
        break;
    }
  }
    
  for(unsigned int i=0;i<pd.path_size;i++){
    FOR_NBORS_OF_ATOM(iter,pd.locant_path[i].atom){
      OBAtom *latom = &(*iter);
      OBBond* lbond = pd.locant_path[i].atom->GetBond(latom);
      if(!atoms_seen[latom]){
        if(!ParseNonCyclic( mol,latom,pd.locant_path[i].atom,lbond,
                            pd.locant_path[i].locant,pd.locant_path,pd.path_size,buffer)){
          fprintf(stderr,"Error: failed on non-cyclic parse\n");
          return false;
        }

      }
    }

    // OM logic 
    if(pd.locant_path[i].atom->GetAtomicNum() == 6 && pd.locant_path[i].atom->GetFormalCharge() == -1){
      FOR_ATOMS_OF_MOL(om,mol){
        OBAtom *organometallic = &(*om);
        if( organometallic->GetAtomicNum() >= 20 && 
            organometallic->GetFormalCharge() > 1 &&
            organometallic->GetExplicitValence() == 0){
          
          unsigned int charge = organometallic->GetFormalCharge(); 
          if(!atoms_seen[organometallic]){
            buffer += ' ';
            buffer += '0';
            WriteSpecial(organometallic,buffer);
            atoms_seen[organometallic] = true;
            pd.locant_path[i].atom->SetFormalCharge(0);
            if(charge)
              charge--;

            // find and write the other rings based on the negative charges
            FOR_ATOMS_OF_MOL(negc,mol){
              OBAtom* next_pi =  &(*negc); 
              if(!atoms_seen[next_pi] && next_pi->GetAtomicNum() == 6
                  && next_pi->GetFormalCharge() == -1 && next_pi->IsInRing()){
                if(!ParseNonCyclic(mol,next_pi,pd.locant_path[i].atom,0,
                            '0',pd.locant_path,pd.path_size,buffer)){
                        

                  fprintf(stderr,"Error: failed on non-cyclic parse\n");
                  return false;
                }

                next_pi->SetFormalCharge(0);
                if(charge)
                  charge--;
                else
                  Fatal("Linking more pi bonded organometallics then charge allows\n");

                organometallic->SetFormalCharge(charge);
              }
            }
          }
        }
      }
    }

  }
  

  // lets try this
  buffer += '&'; // close the ring everytime 

  free(pd.locant_path);  
  return true;
}
#endif
#endif

#endif

unsigned char 
wln_character(OBAtom *atom) 
{
  unsigned int valence = atom->GetExplicitValence(); 
  unsigned int degree  = atom->GetExplicitDegree(); 
  switch (atom->GetAtomicNum()) {

    case BRO: return 'E';
    case IOD: return 'I';
    case FLU: return 'F';
    case CHL: return 'G';
    case BOR: return 'B';
    case PHO: return 'P';
    case SUL: return 'S';

    case OXY: 
      if (valence <= 1) 
        return 'Q';
      else return 'O';
    case NIT: 
      if (valence <= 1) 
        return 'Z';
      if (valence == 2) 
        return 'M';
      if (valence == 4)
        return 'K';
      else return 'N';

    case CAR:
      if (degree == 3) 
        return 'Y';
      if (degree == 4)
        return 'X';
  }
  return '*';  // parse as dash
}


static void 
wln_dash_character(char **buffer, OBAtom *atom)
{
  char *ptr = *buffer;
  *ptr++ = '-';
  switch (atom->GetAtomicNum()) {
    case 5:  *ptr++  = 'B'; break;
    case 7:  *ptr++  = 'N';  break;
    case 8:  *ptr++  = 'O'; break;
    case 9:  *ptr++  = 'F'; break;
    case 53: *ptr++  = 'I'; break;
    case 35: *ptr++  = 'E'; break;
    case 17: *ptr++  = 'G'; break; 
    case 15: *ptr++  = 'P'; break;
    case 89: *ptr++ = 'A';*ptr++ = 'C'; break;
    case 47: *ptr++ = 'A';*ptr++ = 'G'; break;
    case 13: *ptr++ = 'A';*ptr++ = 'L'; break;
    case 95: *ptr++ = 'A';*ptr++ = 'M'; break;
    case 18: *ptr++ = 'A';*ptr++ = 'R'; break;
    case 33: *ptr++ = 'A';*ptr++ = 'S'; break;
    case 85: *ptr++ = 'A';*ptr++ = 'T'; break;
    case 79: *ptr++ = 'A';*ptr++ = 'U'; break;
    case 56: *ptr++ = 'B';*ptr++ = 'A'; break;
    case 4: *ptr++ = 'B';*ptr++ = 'E'; break;
    case 107: *ptr++ = 'B';*ptr++ = 'H'; break;
    case 83: *ptr++ = 'B';*ptr++ = 'I'; break;
    case 97: *ptr++ = 'B';*ptr++ = 'K'; break;
    case 20: *ptr++ = 'C';*ptr++ = 'A'; break;
    case 48: *ptr++ = 'C';*ptr++ = 'D'; break;
    case 58: *ptr++ = 'C';*ptr++ = 'E'; break;
    case 98: *ptr++ = 'C';*ptr++ = 'F'; break;
    case 96: *ptr++ = 'C';*ptr++ = 'N'; break;
    case 112: *ptr++ = 'C';*ptr++ = 'N'; break;
    case 27: *ptr++ = 'C';*ptr++ = 'O'; break;
    case 24: *ptr++ = 'C';*ptr++ = 'R'; break;
    case 55: *ptr++ = 'C';*ptr++ = 'S'; break;
    case 29: *ptr++ = 'C';*ptr++ = 'U'; break;
    case 105: *ptr++ = 'D';*ptr++ = 'B'; break;
    case 110: *ptr++ = 'D';*ptr++ = 'S'; break;
    case 66: *ptr++ = 'D';*ptr++ = 'Y'; break;
    case 68: *ptr++ = 'E';*ptr++ = 'R'; break;
    case 99: *ptr++ = 'E';*ptr++ = 'S'; break;
    case 63: *ptr++ = 'E';*ptr++ = 'U'; break;
    case 26: *ptr++ = 'F';*ptr++ = 'E'; break;
    case 114: *ptr++ = 'F';*ptr++ = 'L'; break;
    case 100: *ptr++ = 'F';*ptr++ = 'M'; break;
    case 87: *ptr++ = 'F';*ptr++ = 'R'; break;
    case 31: *ptr++ = 'G';*ptr++ = 'A'; break;
    case 64: *ptr++ = 'G';*ptr++ = 'D'; break;
    case 32: *ptr++ = 'G';*ptr++ = 'E'; break;
    case 2: *ptr++ = 'H';*ptr++ = 'E'; break;
    case 72: *ptr++ = 'H';*ptr++ = 'F'; break;
    case 80: *ptr++ = 'H';*ptr++ = 'G'; break;
    case 67: *ptr++ = 'H';*ptr++ = 'O'; break;
    case 108: *ptr++ = 'H';*ptr++ = 'S'; break;
    case 49: *ptr++ = 'I';*ptr++ = 'N'; break;
    case 77: *ptr++ = 'I';*ptr++ = 'R'; break;
    case 36: *ptr++ = 'K';*ptr++ = 'R'; break;
    case 19: *ptr++ = 'K';*ptr++ = 'A'; break;
    case 57: *ptr++ = 'L';*ptr++ = 'A'; break;
    case 3: *ptr++ = 'L';*ptr++ = 'I'; break;
    case 103: *ptr++ = 'L';*ptr++ = 'R'; break;
    case 71: *ptr++ = 'L';*ptr++ = 'U'; break;
    case 116: *ptr++ = 'L';*ptr++ = 'V'; break;
    case 115: *ptr++ = 'M';*ptr++ = 'C'; break;
    case 101: *ptr++ = 'M';*ptr++ = 'D'; break;
    case 12: *ptr++ = 'M';*ptr++ = 'G'; break;
    case 25: *ptr++ = 'M';*ptr++ = 'N'; break;
    case 42: *ptr++ = 'M';*ptr++ = 'O'; break;
    case 109: *ptr++ = 'M';*ptr++ = 'T'; break;
    case 11: *ptr++ = 'N';*ptr++ = 'A'; break;
    case 41: *ptr++ = 'N';*ptr++ = 'B'; break;
    case 60: *ptr++ = 'N';*ptr++ = 'D'; break;
    case 10: *ptr++ = 'N';*ptr++ = 'E'; break;
    case 113: *ptr++ = 'N';*ptr++ = 'H'; break;
    case 28: *ptr++ = 'N';*ptr++ = 'I'; break;
    case 102: *ptr++ = 'N';*ptr++ = 'O'; break;
    case 93: *ptr++ = 'N';*ptr++ = 'P'; break;
    case 118: *ptr++ = 'O';*ptr++ = 'G'; break;
    case 76: *ptr++ = 'O';*ptr++ = 'S'; break;
    case 91: *ptr++ = 'P';*ptr++ = 'A'; break;
    case 82: *ptr++ = 'P';*ptr++ = 'B'; break;
    case 46: *ptr++ = 'P';*ptr++ = 'D'; break;
    case 61: *ptr++ = 'P';*ptr++ = 'M'; break;
    case 84: *ptr++ = 'P';*ptr++ = 'O'; break;
    case 59: *ptr++ = 'P';*ptr++ = 'R'; break;
    case 78: *ptr++ = 'P';*ptr++ = 'T'; break;
    case 94: *ptr++ = 'P';*ptr++ = 'U'; break;
    case 88: *ptr++ = 'R';*ptr++ = 'A'; break;
    case 37: *ptr++ = 'R';*ptr++ = 'B'; break;
    case 75: *ptr++ = 'R';*ptr++ = 'E'; break;
    case 104: *ptr++ = 'R';*ptr++ = 'F'; break;
    case 111: *ptr++ = 'R';*ptr++ = 'G'; break;
    case 45: *ptr++ = 'R';*ptr++ = 'H'; break;
    case 86: *ptr++ = 'R';*ptr++ = 'N'; break;
    case 44: *ptr++ = 'R';*ptr++ = 'U'; break;
    case 51: *ptr++ = 'S';*ptr++ = 'B'; break;
    case 21: *ptr++ = 'S';*ptr++ = 'C'; break;
    case 34: *ptr++ = 'S';*ptr++ = 'E'; break;
    case 106: *ptr++ = 'S';*ptr++ = 'G'; break;
    case 14: *ptr++ = 'S';*ptr++ = 'I'; break;
    case 62: *ptr++ = 'S';*ptr++ = 'M'; break;
    case 50: *ptr++ = 'S';*ptr++ = 'N'; break;
    case 38: *ptr++ = 'S';*ptr++ = 'R'; break;
    case 65: *ptr++ = 'T';*ptr++ = 'B'; break;
    case 43: *ptr++ = 'T';*ptr++ = 'C'; break;
    case 52: *ptr++ = 'T';*ptr++ = 'E'; break;
    case 90: *ptr++ = 'T';*ptr++ = 'H'; break;
    case 22: *ptr++ = 'T';*ptr++ = 'I'; break;
    case 81: *ptr++ = 'T';*ptr++ = 'L'; break;
    case 69: *ptr++ = 'T';*ptr++ = 'M'; break;
    case 117: *ptr++ = 'T';*ptr++ = 'S'; break;
    case 92: *ptr++ = 'U';*ptr++ = 'R'; break;
    case 23: *ptr++ = 'V';*ptr++ = 'A'; break;
    case 74: *ptr++ = 'W';*ptr++ = 'T'; break;
    case 54: *ptr++ = 'X';*ptr++ = 'E'; break;
    case 39: *ptr++ = 'Y';*ptr++ = 'T'; break;
    case 70: *ptr++ = 'Y';*ptr++ = 'B'; break;
    case 30: *ptr++ = 'Z';*ptr++ = 'N'; break;
    case 40: *ptr++ = 'Z';*ptr++ = 'R'; break;
  }

  *ptr++ = '-';
  *buffer = ptr; 
}

/* parses the acyclic substructures from a given start point */
bool parse_acyclic(char **buffer, OBMol *mol, OBAtom *atom) 
{
  char *ptr = *buffer; 
  unsigned char ch; 

  unsigned int carbon_chain = 0;
  unsigned int stack_ptr = 0; 

  struct node_edge {
    OBAtom *a;
    OBAtom *p;
    unsigned int order;
  } dfs_stack[256];  

#define add_traversal_node(s,t,_a,_p,o)\
  s[t].a = _a;\
  s[t].p = _p;\
  s[t].order = o;\
  
  unsigned char *seen = (unsigned char*)alloca(mol->NumAtoms()); 
  memset(seen, 0, mol->NumAtoms()); 
  
  ch = wln_character(atom);
  if (ch == '*')
    wln_dash_character(&ptr, atom); 
  seen[atom->GetId()] = ch; 

  FOR_BONDS_OF_ATOM(b, atom) {
    OBBond *bond = &(*b); 
    OBAtom *nxt = bond->GetNbrAtom(atom); 
    if (!seen[nxt->GetId()] && !nxt->IsInRing()) {
      add_traversal_node(dfs_stack, stack_ptr++, nxt, atom, bond->GetBondOrder()); 
    }
  }

  while (stack_ptr > 0) {
    OBAtom *p = dfs_stack[--stack_ptr].p;
    OBAtom *a = dfs_stack[stack_ptr].a;

    for (unsigned int i=1;i<dfs_stack[stack_ptr].order;i++)
      *ptr++ = 'U'; 
    
    ch = wln_character(a);
    if (ch == '*')
      wln_dash_character(&ptr, a); 
    seen[atom->GetId()] = ch; 

    FOR_BONDS_OF_ATOM(b, atom) {
      OBBond *bond = &(*b); 
      OBAtom *nxt = bond->GetNbrAtom(atom); 
      if (!seen[nxt->GetId()] && !nxt->IsInRing()) {
        add_traversal_node(dfs_stack, stack_ptr++, nxt, atom, bond->GetBondOrder()); 
      }
    }
  }

  return true; 
}

bool WriteWLN(char *buffer, OBMol* mol)
{   
  unsigned int i = 0; 
#ifdef USING_OPENBABEL
  bool cyclic = !mol->GetSSSR().empty();
#endif

  if (cyclic) FOR_RINGS_OF_MOL(r, mol) { 

  } 
  else FOR_ATOMS_OF_MOL(a, mol) {

  }
    
  return true; 
}


