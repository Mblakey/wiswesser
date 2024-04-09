/*********************************************************************
 
Author : Michael Blakey

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
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>

#include <set>
#include <vector>
#include <stack>
#include <map>
#include <string>

#include <utility> // std::pair
#include <iterator>
#include <algorithm>

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

#include "parser.h"

using namespace OpenBabel; 

#define WLNDEBUG 0
#define REASONABLE 1024
#define MACROTOOL 0
#define STEREO 0  


#define INT_TO_LOCANT(X) (X+64)
#define LOCANT_TO_INT(X) (X-64)


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


void static write_locant(unsigned char locant,std::string &buffer){
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


static void print_locant_array(LocantPos* locant_path, unsigned int size){
  fprintf(stderr,"[ ");
  for(unsigned int i=0; i<size;i++){
    if(!locant_path[i].atom)
      fprintf(stderr,"0 ");
    else
      fprintf(stderr,"%d ",locant_path[i].atom->GetIdx());
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





/**********************************************************************
                          Write Helper Functions
**********************************************************************/

void write_lowest_ring_locant(OBRing *ring, LocantPos* locant_path, unsigned int plen, std::string &buffer){
  for(unsigned int i=0;i<plen;i++){
    if(ring->IsMember(locant_path[i].atom)){
      if(i>0){
        buffer += ' '; 
        write_locant(INT_TO_LOCANT(i+1),buffer); 
      }
      break;
    }
  }
}


unsigned char highest_ring_locant(OBRing *ring, OBAtom **locant_path, unsigned int plen){
  unsigned char loc = 0; 
  for(unsigned int i=0;i<plen;i++){
    if(ring->IsMember(locant_path[i]))
      loc = i; 
  }
  return INT_TO_LOCANT(loc+1); 
}


void write_ring_size(OBRing *ring, std::string &buffer){
  if(ring->Size() < 9)
    buffer += ring->Size() + '0'; 
  else{
#if MODERN
    buffer += '-';
    buffer += std::to_string(ring->Size()); 
    buffer += '-';
#else
    buffer += '-';
    buffer += std::to_string(ring->Size()); 
    buffer += '-';
#endif
  }

}


/**********************************************************************
                          Locant Path Functions
**********************************************************************/


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


void print_ring_locants(OBMol *mol,OBRing *ring, LocantPos*locant_path, unsigned int path_size){
  unsigned char *sequence = (unsigned char*)malloc(sizeof(unsigned char)*ring->Size()); 
  
  for(unsigned int i=0;i<ring->Size();i++){
    sequence[i] = locant_path[position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)].locant; 
  }

  fprintf(stderr,"[ ");
  for(unsigned int k=0;k<ring->Size();k++)
    fprintf(stderr,"%c ",sequence[k]);
  fprintf(stderr,"]\n");

  free(sequence);
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
      print_ring_locants(mol,to_write,locant_path,path_size,false);
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
    if(ring->IsMember(locant_path[i].atom))
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
Some rules to follow when walking the path:
*/
void write_complete_rings(  LocantPos *locant_path, unsigned int locant_pos, 
                            std::set<OBRing*> &local_SSSR, std::map<OBRing*,bool> &handled_rings,
                            std::vector<OBRing*> &ring_order,std::string &buffer)
{
  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end(); riter++){
    if(!handled_rings[*riter] && IsRingComplete(*riter, locant_path, locant_pos)){
      write_lowest_ring_locant(*riter, locant_path, locant_pos, buffer);
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
LocantPos *PolyWalk(    OBMol *mol, unsigned int path_size,
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

  std::string best_notation; 
  std::vector<OBRing*> best_order; 
    
  OBAtom*                ratom  = 0; // ring
  OBAtom*                catom  = 0; // child
  OBAtom*                matom  = 0; // move atom
  unsigned int           lowest_sum       = UINT32_MAX;
  
  for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
    if(atom_shares[*aiter] == 2){ // these are the starting points 
     
      zero_locant_path(locant_path, path_size);

      std::string poly_buffer = "";
      std::vector<OBRing*> lring_order; 
      std::map<OBAtom*,bool> visited; 
      std::map<OBRing*,bool> handled_rings; 
      unsigned int locant_pos = 0;
      
      ratom = *aiter; 
      for(;;){

        locant_path[locant_pos].atom = ratom; 
        locant_path[locant_pos].locant = INT_TO_LOCANT(locant_pos+1);        
        locant_pos++; 

        visited[ratom] = true;

        write_complete_rings(locant_path, locant_pos, local_SSSR, handled_rings, lring_order,poly_buffer); 
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
        best_notation = poly_buffer; 
        best_order = lring_order; 
      }
    }
  }

  free(locant_path);
  for(unsigned int i=0;i<path_size;i++){
    if(!best_path[i].atom){
      free(best_path);
      return 0; 
    }
  }
  
  ring_order = best_order; 
  buffer = best_notation; 
  return best_path; 
}


void BackTrackWalk(  OBAtom *clear, LocantPos*locant_path, unsigned int path_size,unsigned int &locant_pos, std::map<OBAtom*,bool> &visited_atoms)
{

  // find the position in the path where the lowest one of these are, the other gets placed into locant path 
  unsigned int p=0;
  unsigned int q=0;
  for(p=0;p<path_size;p++,q++){
    if(locant_path[p].atom == clear){
      break;
    }
  }

  // everything from p gets cleared, but not including
  for(++p;p<path_size;p++){
    visited_atoms[locant_path[p].atom] = 0;
    locant_path[p].atom = 0; 
  }
  
  locant_pos = q+1;
}


bool duplicates_path(OBAtom **locant_path, unsigned int path_size){
  bool f = 0; 
  for(unsigned int i=0;i<path_size;i++){
    for(unsigned int j=i+1;j<path_size;j++){
      if(locant_path[i] && locant_path[j] && locant_path[i] == locant_path[j]){
        fprintf(stderr,"duplicate in path - pos %d == %d\n",i,j); 
        f = 1;
      }
    }
  }
  return f; 
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
LocantPos *PeriWalk2(   OBMol *mol,        unsigned int &path_size,
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

  OBAtom*                ratom  = 0; // ring
  OBAtom*                catom  = 0; // child
  OBAtom*                matom  = 0; // move atom
  unsigned int           lowest_sum = UINT32_MAX;
  unsigned int           starting_path_size = path_size; // important if path size changes
  std::string best_notation; 
  std::vector<OBRing*> best_order; 

  for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
    // a multicyclic that connects to two other multicyclic points can never be the start, always take an edge case
    if( (atom_shares[*aiter] >= 3 && connected_multicycles(*aiter,atom_shares)<=1)  || bridge_atoms[*aiter]){ // these are the starting points 

      std::string             peri_buffer; 
      std::vector<OBRing*>    lring_order; 
      std::map<OBAtom*,bool>  visited; 
      std::map<OBRing*,bool>  handled_rings;
      std::stack<OBAtom*>     multistack; 
      std::stack<std::pair<OBAtom*,OBAtom*>> backtrack_stack;   // multicyclics have three potential routes, 
      
      unsigned int locant_pos = 0;
      zero_locant_path(locant_path, path_size);
      
      
      // used for off branch locants when needed
      for(std::set<OBAtom*>::iterator multi_iter = ring_atoms.begin(); multi_iter != ring_atoms.end(); multi_iter++){
        if(atom_shares[*multi_iter]>=3 && *multi_iter != *aiter)
          multistack.push(*multi_iter); 
      }

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
          visited[ratom] = true;

          write_complete_rings(locant_path, locant_pos, local_SSSR, handled_rings, lring_order,peri_buffer); 
          
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
                if( (atom_shares[ratom]>=3 || atom_shares[catom]>=3)){
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
        }
       
        if(locant_pos == path_size){
          unsigned int fsum = fusion_sum(mol,locant_path,path_size,local_SSSR);
          if(fsum < lowest_sum){ // rule 30d.
            lowest_sum = fsum;
            copy_locant_path(best_path,locant_path,path_size);
            best_notation = peri_buffer; 
            best_order = lring_order;
          }
        }
        
        if(!backtrack_stack.empty()){
          ratom = backtrack_stack.top().second;
          BackTrackWalk(backtrack_stack.top().first, locant_path, path_size,locant_pos,visited); 
          
          // this is expensive but guarantees sequential ordeirng
          peri_buffer.clear();
          handled_rings.clear(); 
          lring_order.clear(); 
          for(unsigned int t=0;t<locant_pos;t++)
            write_complete_rings(locant_path, t, local_SSSR, handled_rings, lring_order,peri_buffer); 
        }
        else if(!best_path[0].atom && !multistack.empty()){
          // this the where the broken locants happen, pop off a multistack atom 
          OBAtom *branch_locant = multistack.top();
          
          multistack.pop();
          visited.clear(); 
          peri_buffer.clear();
          handled_rings.clear();
          lring_order.clear();
          visited[branch_locant] = true;
          
          // set the branch locant value here
          path_size--; // decrement the path size, this is globally changed
          locant_path[path_size].atom = branch_locant; 
          locant_path[path_size].locant = 'X';
          
          goto path_solve; 
        }

      } while(!backtrack_stack.empty()) ; 
    }
  } 

  free(locant_path);
  for(unsigned int i=0;i<path_size;i++){
    if(!best_path[i].atom){
      fprintf(stderr,"Error: locant path is missing a value at %d\n",i); 
      free(best_path);
      return 0; 
    }
  }
  
  ring_order = best_order; 
  buffer = best_notation; 
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

  void WriteSpecial(OBAtom *atom, std::string &buffer){
    if(!atom)
      Fatal("writing notation from dead atom ptr");
    // all special elemental cases
    //
    
#if MODERN
    buffer += "<";
#else
    buffer += "-";
#endif
    string_position[atom] = buffer.size()+1; // always first character 
    switch(atom->GetAtomicNum()){
#if MODERN
      case 5:
        buffer += "-";
        buffer += "B";
        break;

      case 7:
        buffer += "-";
        buffer += "N";
        break;
      
      case 8:
        buffer += "-";
        buffer += "O";
        break;

      case 9:
        buffer += "-";
        buffer += "F";
        break;

      case 53:
        buffer += "-";
        buffer += "I";
        break;

      case 35:
        buffer += "-";
        buffer += "E";
        break;

      case 17:
        buffer += "-";
        buffer += "G";
        break; 
    
      case 15:
        buffer += "-";
        buffer += "P";
        break;

#else
      case 5:
        buffer += "B";
        break;

      case 7:
        buffer += "N";
        break;
      
      case 8:
        buffer += "O";
        break;

      case 9:
        buffer += "F";
        break;

      case 53:
        buffer += "I";
        break;

      case 35:
        buffer += "E";
        break;

      case 17:
        buffer += "G";
        break; 
      
      case 15:
        buffer += "P";
        break;
#endif
      case 89:
        buffer += "AC";
        break;

      case 47:
        buffer += "AG";
        break;
    
      case 13:
        buffer += "AL";
        break;

      case 95:
        buffer += "AM";
        break;

      case 18:
        buffer += "AR";
        break;

      case 33:
        buffer += "AS";
        break;

      case 85:
        buffer += "AT";
        break;

      case 79:
        buffer += "AU";
        break;


      case 56:
        buffer += "BA";
        break;

      case 4:
        buffer += "BE";
        break;

      case 107:
        buffer += "BH";
        break;

      case 83:
        buffer += "BI";
        break;

      case 97:
        buffer += "BK";
        break;

      case 20:
        buffer += "CA";
        break;
      
      case 48:
        buffer += "CD";
        break;

      case 58:
        buffer += "CE";
        break;

      case 98:
        buffer += "CF";
        break;

      case 96:
        buffer += "CN";
        break;

      case 112:
        buffer += "CN";
        break;

      case 27:
        buffer += "CO";
        break;

      case 24:
        buffer += "CR";
        break;

      case 55:
        buffer += "CS";
        break;

      case 29:
        buffer += "CU";
        break;

      case 105:
        buffer += "DB";
        break;

      case 110:
        buffer += "DS";
        break;

      case 66:
        buffer += "DY";
        break;

      case 68:
        buffer += "ER";
        break;

      case 99:
        buffer += "ES";
        break;

      case 63:
        buffer += "EU";
        break;

      case 26:
        buffer += "FE";
        break;

      case 114:
        buffer += "FL";
        break;

      case 100:
        buffer += "FM";
        break;

      case 87:
        buffer += "FR";
        break;

      case 31:
        buffer += "GA";
        break;

      case 64:
        buffer += "GD";
        break;

      case 32:
        buffer += "GE";
        break;

      case 2:
        buffer += "HE";
        break;

      case 72:
        buffer += "HF";
        break;

      case 80:
        buffer += "HG";
        break;

      case 67:
        buffer += "HO";
        break;

      case 108:
        buffer += "HS";
        break;

      case 49:
        buffer += "IN";
        break;

      case 77:
        buffer += "IR";
        break;

      case 36:
        buffer += "KR";
        break;

      case 19:
        buffer += "KA";
        break;

      case 57:
        buffer += "LA";
        break;

      case 3:
        buffer += "LI";
        break;

      case 103:
        buffer += "LR";
        break;

      case 71:
        buffer += "LU";
        break;

      case 116:
        buffer += "LV";
        break;

      case 115:
        buffer += "MC";
        break;

      case 101:
        buffer += "MD";
        break;

      case 12:
        buffer += "MG";
        break;

      case 25:
        buffer += "MN";
        break;

      case 42:
        buffer += "MO";
        break;

      case 109:
        buffer += "MT";
        break;

      case 11:
        buffer += "NA";
        break;

      case 41:
        buffer += "NB";
        break;

      case 60:
        buffer += "ND";
        break;

      case 10:
        buffer += "NE";
        break;

      case 113:
        buffer += "NH";
        break;

      case 28:
        buffer += "NI";
        break;

      case 102:
        buffer += "NO";
        break;

      case 93:
        buffer += "NP";
        break;

      case 118:
        buffer += "OG";
        break;

      case 76:
        buffer += "OS";
        break;

      case 91:
        buffer += "PA";
        break;

      case 82:
        buffer += "PB";
        break;

      case 46:
        buffer += "PD";
        break;

      case 61:
        buffer += "PM";
        break;

      case 84:
        buffer += "PO";
        break;

      case 59:
        buffer += "PR";
        break;

      case 78:
        buffer += "PT";
        break;

      case 94:
        buffer += "PU";
        break;

      case 88:
        buffer += "RA";
        break;

      case 37:
        buffer += "RB";
        break;

      case 75:
        buffer += "RE";
        break;

      case 104:
        buffer += "RF";
        break;

      case 111:
        buffer += "RG";
        break;

      case 45:
        buffer += "RH";
        break;

      case 86:
        buffer += "RN";
        break;

      case 44:
        buffer += "RU";
        break;

      case 51:
        buffer += "SB";
        break;

      case 21:
        buffer += "SC";
        break;

      case 34:
        buffer += "SE";
        break;

      case 106:
        buffer += "SG";
        break;

      case 14:
        buffer += "SI";
        break;

      case 62:
        buffer += "SM";
        break;

      case 50:
        buffer += "SN";
        break;

      case 38:
        buffer += "SR";
        break;


      case 73:
        buffer += "TA";
        break;

      case 65:
        buffer += "TB";
        break;

      case 43:
        buffer += "TC";
        break;

      case 52:
        buffer += "TE";
        break;

      case 90:
        buffer += "TH";
        break;

      case 22:
        buffer += "TI";
        break;

      case 81:
        buffer += "TL";
        break;

      case 69:
        buffer += "TM";
        break;

      case 117:
        buffer += "TS";
        break;

      case 92:
        buffer += "UR";
        break;

      case 23:
        buffer += "VA";
        break;
      
      case 74:
        buffer += "WT";
        break;

      case 54:
        buffer += "XE";
        break;

      case 39:
        buffer += "YT";
        break;

      case 70:
        buffer += "YB";
        break;

      case 30:
        buffer += "ZN";
        break;

      case 40:
        buffer += "ZR";
        break;
    }

#if MODERN
    ModernCharge(atom, buffer); 
    buffer += ">";
#else
    buffer += "-";
#endif
    
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

  /* parse non-cyclic atoms DFS style - return last atom seen in chain */
  bool ParseNonCyclic(OBMol *mol, OBAtom* start_atom, OBAtom *spawned_from,OBBond *inc_bond,
                      unsigned char locant, LocantPos *locant_path, unsigned int path_size, 
                      std::string &buffer)
  {
    if(!start_atom)
      Fatal("writing notation from dead atom ptr");

    unsigned int border = 0; 
    unsigned char stereo = 0; 

    if(inc_bond){
      border = inc_bond->GetBondOrder();
      if(inc_bond->IsHash())
        stereo = 'D';
      else if (inc_bond->IsWedge())
        stereo =  'A'; 
    }

    if(locant && locant != '0' && border > 0){ // allows OM through
      buffer+=' ';
      write_locant(locant,buffer);
    }

    
#if MODERN
    if(stereo)
      buffer += stereo;
#endif

    for(unsigned int b=1;b<border;b++)
      buffer+='U';
    

    unsigned int carbon_chain = 0;
    unsigned char wln_character = 0; 

    OBAtom* atom = start_atom;
    OBAtom* prev = 0; 
    OBBond *bond = 0; 
   
    std::stack<OBAtom*> atom_stack; 
    std::stack<OBAtom*> branch_stack;
    std::map<OBAtom*,bool> branching_atom; 
    atom_stack.push(atom);

    bool require_macro_closure = false;

    while(!atom_stack.empty()){
      atom = atom_stack.top(); 
      atom_stack.pop();
      atoms_seen[atom] = true;
    
      if(prev){
        bond = mol->GetBond(prev,atom); 
        if(!bond && !branch_stack.empty()){
          
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }
          
          // a closure is any symbol that you would expect to have a symbol on either side
          // e.g 2N2 or even 2UUN<...> since N normally has two symbols either side it requires 
          // a closure, this can be stated with a degree check, terminators are not considered here
          if(!branching_atom[prev]  )
            buffer += '&'; // requires a closure 
          else if (!branch_stack.empty() && prev == branch_stack.top() && prev->GetExplicitDegree() == 1){
            buffer += '&'; 
            branch_stack.pop(); // returns immedietely
          }

          while(!branch_stack.empty()){
            prev = branch_stack.top();
            if(mol->GetBond(atom,prev)){
              bond = mol->GetBond(atom,prev); 
              break;
            }
            else{
              if (remaining_branches[prev] > 0)
                buffer += '&';
              
              branch_stack.pop();
            }
          }
        }

        remaining_branches[prev]--; // reduce the branches remaining  

#if MODERN && STEREO
        if(bond->IsHash())
          buffer += 'D';
        else if (bond->IsWedge())
          buffer += 'A'; 
#endif

        for(unsigned int i=1;i<bond->GetBondOrder();i++){
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }
          
          buffer += 'U';
          if(prev->GetAtomicNum() != 6)  // branches unaffected when carbon Y or X
            remaining_branches[prev]--;
        }
      }

      if(atom->IsInRing()){
        if(carbon_chain){
          buffer += std::to_string(carbon_chain);
          carbon_chain = 0;
        }

        if(locant == '0' && !border){
          buffer += '-';
          buffer += ' ';
          buffer += '0';
          if(!RecursiveParse(mol,atom,spawned_from,false,buffer))
            Fatal("failed to make pi bonded ring");
        }
        else{
          if(!RecursiveParse(mol,atom,spawned_from,true,buffer))
            Fatal("failed to make inline ring");
        }
        
        // this should count as a branch?, lets see - doesnt seem
        // if(prev) 
        //   remaining_branches[prev]--; // reduce the branches remaining  

        if(!branch_stack.empty())
          prev = return_open_branch(branch_stack);
        
        continue;
      }
        
       // remaining_branches are -1, we only look forward
      unsigned int correction = 0; 
      wln_character =  WriteSingleChar(atom);
      atom_chars[atom] = wln_character; 
    
      if(prev && bond)
        correction = bond->GetBondOrder() - 1;
      else if (border > 0)
        correction = border - 1;

      // last added char, not interested in the ring types of '-'
      switch(wln_character){
// oxygens
        case 'O':
        case 'V':
        case 'M':
        case 'W': // W is not actually seen as a prev,
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }
          
          prev = atom; 
#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else 
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          string_position[atom] = buffer.size(); 
          break;

// carbons 
        // alkyl chain 
        case '1':
          prev = atom; 
          if(CheckCarbonyl(atom)){
            if(carbon_chain){
              buffer += std::to_string(carbon_chain);
              carbon_chain = 0;
            }
            buffer += 'V';
            string_position[atom] = buffer.size();
          }
#if MODERN
          else if (atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += 'C';
            ModernCharge(atom, buffer); 
            buffer += '>'; 
            string_position[atom] = buffer.size()+1; // writen next so offset by 1
          }
#endif
          else{
            string_position[atom] = buffer.size()+1; // writen next so offset by 1
            carbon_chain++; 
          }
          
          break;

        case 'Y':
        case 'X':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom;
          if(CheckCarbonyl(atom))
            buffer += 'V';
          else{
#if MODERN
            if(atom->GetFormalCharge() != 0){
              buffer += '<'; 
              buffer += wln_character; 
              ModernCharge(atom, buffer); 
              buffer += '>'; 
            }
            else 
              buffer += wln_character; 
#else
            buffer += wln_character; 
#endif
            if(wln_character == 'X')
              remaining_branches[atom] += 3;
            else
              remaining_branches[atom] += 2; 

            branching_atom[atom] = true;
            branch_stack.push(atom);
          }
          string_position[atom] = buffer.size();
          break;

        case 'N':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else 
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          prev = atom; 
          string_position[atom] = buffer.size();
          
          for(unsigned int h=0;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 

          remaining_branches[atom] += 2 - correction; 
          branch_stack.push(atom);
          branching_atom[atom] = true;
          break;


        case 'B':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom; 

#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else 
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          string_position[atom] = buffer.size();
          
          for(unsigned int h=0;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 
         
          remaining_branches[atom] += 2 - correction; 
          branch_stack.push(atom);
          branching_atom[atom] = true;
          break;

        case 'K':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom;
          buffer += wln_character;
          string_position[atom] = buffer.size();

          for(unsigned int h=0;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 
          
          // K now given for all positive nitrogen - NO!
          remaining_branches[atom] += 3 - correction;
          branching_atom[atom] = true; 
          branch_stack.push(atom);
          atom->SetFormalCharge(0); // remove the charge, as this is expected 
          break;
        
        case 'P':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom;

#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else 
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          string_position[atom] = buffer.size();
          
          for(unsigned int h=0;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 
          
          remaining_branches[atom] += 4 - correction; 
          branching_atom[atom] = true;
          branch_stack.push(atom);
          break;

        case 'S':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom;

#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else 
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          string_position[atom] = buffer.size();
          for(unsigned int h=0;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 

          remaining_branches[atom] += 5 - correction; 
          branching_atom[atom] = true;
          branch_stack.push(atom);
          break;

        case '*':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom; 
          WriteSpecial(atom,buffer);
          for(unsigned int h=0;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 
          
          switch (atom->GetAtomicNum()) {
            case 8:
              remaining_branches[atom] += 2 - correction; // oxygen gets expanded to 3 only
              break;

            default:
              remaining_branches[atom] += 5 - correction; // octdhedral max geometry 
          }
            
          branching_atom[atom] = true;
          branch_stack.push(atom);
          break;
          
// terminators
        case 'Q':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          string_position[atom] = buffer.size();
          for(unsigned int h=1;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 

          if(!branch_stack.empty())
            prev = return_open_branch(branch_stack);

          if(atom->GetExplicitDegree() == 0)
            atom->SetFormalCharge(0); // notational implied, do not write ionic code

          break;

        case 'Z':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else  
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          string_position[atom] = buffer.size();
          for(unsigned int h=2;h<atom->GetImplicitHCount();h++)
            buffer += 'H'; 

          if(!branch_stack.empty())
            prev = return_open_branch(branch_stack);

          if(atom->GetExplicitDegree() == 0 && !atom->GetImplicitHCount())
            atom->SetFormalCharge(0); // notational implied, do not write ionic code
          break;

        case 'E':
        case 'F':
        case 'G':
        case 'I':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

#if MODERN
          if(atom->GetFormalCharge() != 0 && atom->GetHeteroDegree() > 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else  
            buffer += wln_character; 
#else
          buffer += wln_character; 
#endif
          string_position[atom] = buffer.size();

          if(!branch_stack.empty())
            prev = return_open_branch(branch_stack);

          if(atom->GetExplicitDegree() == 0)
            atom->SetFormalCharge(0); // notational implied, do not write ionic code
          
          if(atom->GetImplicitHCount())
            buffer += 'H'; 
          break;

        case 'H':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }
#if MODERN
          if(atom->GetFormalCharge() != 0){
            buffer += '<'; 
            buffer += wln_character; 
            ModernCharge(atom, buffer); 
            buffer += '>'; 
          }
          else  
            buffer += wln_character; 
#else
          buffer += wln_character;
#endif
          string_position[atom] = buffer.size();
          break;

        default:
          fprintf(stderr,"Error: unhandled char %c\n",wln_character); 
          return 0; 
      }

      // here we ask, is this bonded to a ring atom that is not 'spawned from'
      FOR_NBORS_OF_ATOM(a,atom){
        OBAtom *nbor = &(*a);
        if(nbor != spawned_from && nbor->IsInRing() && atoms_seen[nbor] == true){
  
          if(require_macro_closure){
            fprintf(stderr,"Error: macro-closure appearing more than once\n");
            return 0;
          }
          else{
            require_macro_closure = true;
            
            if(carbon_chain){
              buffer += std::to_string(carbon_chain);
              carbon_chain = 0;
            }

            OBBond *macro_bond = mol->GetBond(atom, nbor);
            if(branching_atom[atom])
              remaining_branches[atom]--; 

            if(macro_bond->GetBondOrder() > 1){
#if MODERN && STEREO      
              if(bond->IsHash())
                buffer += 'D';
              else if (bond->IsWedge())
                buffer += 'A'; 
#endif
              buffer += 'U'; 
              if(branching_atom[atom] && atom->GetAtomicNum() != 6)  // branches unaffected when carbon Y or X
                remaining_branches[prev]--;
            }

            buffer += '-';
            buffer += ' ';
            
            if(!locant_path){
              fprintf(stderr,"Error: no locant path to wrap back macro-closures\n");
              return 0;
            }
            else{
              for (unsigned int i=0;i<path_size;i++) {
                if(locant_path[i].atom == nbor){
                  write_locant(locant_path[i].locant, buffer);
                  break;
                }
              }
            }

            buffer += "-x-";
            // the wrapper needs to be handled on closure of all the locants in a ring... hmmm
            break;
          }
        }
      }

      FOR_NBORS_OF_ATOM(a,atom){
        if(!atoms_seen[&(*a)])
          atom_stack.push(&(*a));
      }
    }

    if(carbon_chain){
      buffer += std::to_string(carbon_chain);
      carbon_chain = 0;
    }

    if(require_macro_closure){
      buffer += 'J';
    }
    
    // burn the branch stack here, recursion takes care of the rest 
    // may be over cautious
        
    if(!branch_stack.empty()){
      OBAtom *top = return_open_branch(branch_stack);
      
      if(top && prev != top){
        switch (wln_character) {
            case 'E':
            case 'F':
            case 'G':
            case 'H':
            case 'I':
            case 'Q':
            case 'Z':
            case 'W':
              break;

            default:
              buffer += '&';
        }
      }

      while(!branch_stack.empty()){
        if(remaining_branches[branch_stack.top()] > 0){
          buffer += '&';
        }
        branch_stack.pop(); 
      }
    }

    return atom; 
  }

  void AddPostCharges(OBMol *mol,std::string &buffer){
    if(OPT_DEBUG)
      fprintf(stderr,"Post Charges\n");
    
    bool working = true;
    while(working){
      working = false;
      FOR_ATOMS_OF_MOL(a,mol){
        OBAtom *atom = &(*a);
        if(atom->GetFormalCharge() != 0){

          if(OPT_DEBUG)
            fprintf(stderr,"  adding charge %d to atomic num: %d\n",atom->GetFormalCharge(),atom->GetAtomicNum());

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
  bool  ConstructLocalSSRS( OBMol *mol, OBAtom *ring_root,
                            std::set<OBAtom*>         &ring_atoms,
                            std::set<OBBond*>         &ring_bonds,
                            std::map<OBAtom*,bool>    &bridge_atoms,
                            std::map<OBAtom*,unsigned int> &atom_shares,
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
            if(bond)
              ring_bonds.insert(bond); 
            
            prev = ratom; 
          }
        }
        // get the last bond
        bond = mol->GetBond(mol->GetAtom(obring->_path.front()),mol->GetAtom(obring->_path.back()));
        if(bond)
          ring_bonds.insert(bond); 

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
                if(bond)
                  ring_bonds.insert(bond); 
                
                prev = ratom; 
              }
            }
            // get the last bond
            bond = mol->GetBond(mol->GetAtom(obring->_path.front()),mol->GetAtom(obring->_path.back()));
            if(bond)
              ring_bonds.insert(bond); 
            
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



    if(OPT_DEBUG){
      fprintf(stderr,"  ring atoms: %lu\n",ring_atoms.size());
      fprintf(stderr,"  ring bonds: %lu\n",ring_bonds.size());
      fprintf(stderr,"  ring subcycles: %lu/%lu\n",local_SSSR.size(),mol->GetSSSR().size());
      if(bridge_count)
        fprintf(stderr,"  bridging atoms: %d\n",bridge_count);
      
    }


    local_data.bridging = bridge_count;
    local_data.path_size = ring_atoms.size();  
    return true; 
  }


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
      
      if(!locant_path[i].atom)
        Fatal("dead locant path atom ptr in hetero read");

      if(!locant_path[i].locant)
        Fatal("dead locant path position in hetero read");
      
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

#if MODERN
            if(locant_atom->GetFormalCharge() != 0){
              buffer += '<';
              buffer += het_char; 
              ModernCharge(locant_path[i], buffer);
              buffer += '>';
            }
            else 
              buffer += het_char; 
#else
            buffer+=het_char;
#endif
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

#if MODERN && STEREO      
        if(locant_bond->IsHash())
          buffer += 'D';
        else if (locant_bond->IsWedge())
          buffer += 'A'; 
#endif

        for(unsigned int b=1;b<locant_bond->GetBondOrder();b++)
          buffer += 'U';
      }
#if MODERN && STEREO
      else if(locant_bond && (locant_bond->IsWedge()||locant_bond->IsHash())){
        buffer += ' ';
        write_locant(locant,buffer);

        if(locant_bond->IsHash())
          buffer += 'D';
        else if (locant_bond->IsWedge())
          buffer += 'A'; 
      }
#endif
    }


    for(std::set<OBBond*>::iterator biter = ring_bonds.begin(); biter != ring_bonds.end();biter++){
      OBBond *fbond = *biter; 
      
      if(!bonds_checked[fbond] && fbond->GetBondOrder() > 1 && !fbond->IsAromatic()){
        unsigned int floc = locant_path[position_in_path(fbond->GetBeginAtom(),locant_path,path_size)].locant; 
        unsigned int bloc = locant_path[position_in_path(fbond->GetEndAtom(),locant_path,path_size)].locant; 
        
        buffer += ' ';
        write_locant(floc,buffer);
#if MODERN && STEREO
        if(fbond->IsHash())
          buffer += 'D';
        else if (fbond->IsWedge())
          buffer += 'A'; 
#endif
        for(unsigned int b=1;b<fbond->GetBondOrder();b++)
          buffer += 'U';
        buffer+='-';
        buffer+=' ';
        write_locant(bloc,buffer);
        break;
      }
#if MODERN && STEREO 
      else if(!bonds_checked[fbond] && fbond->IsWedgeOrHash()){

        unsigned char floc = INT_TO_LOCANT(position_in_path(fbond->GetBeginAtom(),locant_path,path_size)+1); 
        unsigned char bloc = INT_TO_LOCANT(position_in_path(fbond->GetEndAtom(),locant_path,path_size)+1); 
        buffer += ' ';
        write_locant(floc,buffer);
        
        if(fbond->IsHash())
          buffer += 'D';
        else if (fbond->IsWedge())
          buffer += 'A'; 
        
        buffer+='-';
        buffer+=' ';
        write_locant(bloc,buffer);
      }
#endif
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
        write_locant(INT_TO_LOCANT(i+1),append);
      }
    }

#if MODERN
    buffer += ' ';
    buffer+= std::to_string(count);
#else
    buffer += ' ';
    buffer+= std::to_string(count);
    buffer+= append;
#endif
  }


  /* constructs and parses a cyclic structure, locant path is returned with its path_size */
  void ParseCyclic(OBMol *mol, OBAtom *ring_root,OBAtom *spawned_from,bool inline_ring,PathData &pd,std::string &buffer){
    if(OPT_DEBUG)
      fprintf(stderr,"Reading Cyclic\n");

    LocantPos*                      locant_path = 0; 
    std::set<OBRing*>               local_SSSR;
    std::set<OBAtom*>               ring_atoms;
    std::set<OBBond*>               ring_bonds;
    std::vector<OBRing*>            ring_order; 

    std::map<OBAtom*,bool>          bridge_atoms;
    std::map<OBAtom*,unsigned int>  atom_shares;
    
    // can we get the local SSSR data out of this?
    SubsetData LocalSSRS_data;

    if(!ConstructLocalSSRS(mol,ring_root,ring_atoms,ring_bonds,bridge_atoms,atom_shares,local_SSSR, LocalSSRS_data)){
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
      locant_path = PolyWalk(mol,path_size,ring_atoms,atom_shares,bridge_atoms,local_SSSR,ring_order,ring_segment);
    else
      locant_path = PeriWalk2(mol,path_size, ring_atoms, atom_shares, bridge_atoms, local_SSSR,ring_order,ring_segment); 
    if(!locant_path)
      return Fatal("no locant path could be determined");

    branch_locants = LocalSSRS_data.path_size - path_size; 
    
    if(OPT_DEBUG){
      fprintf(stderr,"  multi-cyclic: %d\n",multi);
      fprintf(stderr,"  bridging:     %d\n",bridging);
      fprintf(stderr,"  off branches: %d\n",branch_locants);
    }

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
      for(unsigned int i=0;i<path_size;i++){
        if(bridge_atoms[locant_path[i].atom]){
          buffer+= ' ';
          unsigned char bloc = INT_TO_LOCANT(i+1);
          write_locant(bloc,buffer);
        }
      }
    }

    if(multi){
      ReadMultiCyclicPoints(locant_path,path_size,atom_shares,buffer);
      
#if MODERN
      write_locant(INT_TO_LOCANT(path_size),buffer); // need to make the relative size
#else
      buffer += ' ';
      write_locant(INT_TO_LOCANT(path_size),buffer); // need to make the relative size
#endif
    }

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
    

  bool RecursiveParse(OBMol *mol, OBAtom *atom, OBAtom *spawned_from, bool inline_ring,std::string &buffer){
    // assumes atom is a ring atom 
    
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


};



/**********************************************************************
                         API FUNCTION
**********************************************************************/

bool WriteWLN(std::string &buffer, OBMol* mol, bool modern)
{   
  
  OBMol *mol_copy = new OBMol(*mol); // performs manipulations on the mol object, copy for safety
  BabelGraph obabel;

#define PERCEPTION_DEBUG 1
#if PERCEPTION_DEBUG
    fprintf(stderr,"debugging ring perception:\n"); 
    
    unsigned int rc = 0;
    FOR_RINGS_OF_MOL(r,mol){
      OBRing *obring = &(*r);
      fprintf(stderr,"%d [ ",rc++);
      for(unsigned int i=0;i<obring->Size();i++){
        OBAtom *ratom = mol->GetAtom(obring->_path[i]);
        fprintf(stderr,"%d ",ratom->GetIdx());
      }
      fprintf(stderr,"]\n");
    }

#endif

#if MODERN
  StereoFrom0D(mol_copy); 

  OBStereoFacade stereo_facade = OBStereoFacade(mol_copy);
  FOR_ATOMS_OF_MOL(a,mol_copy){
    if(stereo_facade.HasTetrahedralStereo(a->GetId())){
    
      OBTetrahedralStereo *tet_centre = stereo_facade.GetTetrahedralStereo(a->GetId()); 
      OBTetrahedralStereo::Config cfg;
      cfg = tet_centre->GetConfig();
      
      unsigned int stereo = 0;
      for (unsigned long ref : cfg.refs){
        if(mol_copy->GetAtomById(ref)){
          OBBond * bond = mol_copy->GetBond(mol_copy->GetAtomById(cfg.center),mol_copy->GetAtomById(ref));

          if(stereo==1){
            bond->SetWedge(true); 
          }
          else if (stereo==0){
            bond->SetHash(true); 
          }
        }
        stereo++; 
      }
      
    }
  }
#endif

  unsigned int cyclic = 0;
  bool started = false; 
  FOR_RINGS_OF_MOL(r,mol_copy)
    cyclic++;

  if(OPT_DEBUG)
    WriteBabelDotGraph(mol_copy);

  if(!cyclic){
    FOR_ATOMS_OF_MOL(a,mol_copy){
      OBAtom *satom = &(*a); 
      if(!obabel.atoms_seen[satom] && (satom->GetExplicitDegree()==1 || satom->GetExplicitDegree() == 0) ){
        if(started)
          buffer += " &"; // ionic species
        if(!obabel.ParseNonCyclic(mol_copy,&(*a),0,0,0,0,0,buffer))
          Fatal("failed on recursive branch parse");

        started = true; 
      }
    }
  }
  else{
    FOR_RINGS_OF_MOL(r,mol_copy){
    // start recursion from first cycle atom
      if(!obabel.rings_seen[&(*r)]){
        if(started){
          buffer += " &"; // ionic species
        }
        
        if(!obabel.RecursiveParse(mol_copy,mol_copy->GetAtom( (&(*r))->_path[0]),0,false,buffer))
          Fatal("failed on recursive ring parse");

        started = true;
      }
    }
    
    // handles additional ionic atoms here
    FOR_ATOMS_OF_MOL(a,mol_copy){
      OBAtom *satom = &(*a); 
      if(!obabel.atoms_seen[satom] && (satom->GetExplicitDegree()==1 || satom->GetExplicitDegree() == 0) ){
        buffer += " &"; // ionic species
        if(!obabel.ParseNonCyclic(mol_copy,satom,0,0,0,0,0,buffer))
          Fatal("failed on recursive branch parse");
      }
    }
  }

#if !MODERN
  obabel.AddPostCharges(mol_copy,buffer); // add in charges where we can 
#endif 

  while(buffer.back() == '&')
    buffer.pop_back(); 

  delete mol_copy; 
  return true; 
}



