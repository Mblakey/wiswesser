/**********************************************************************
 
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include <set>
#include <deque>
#include <vector>
#include <stack>
#include <map>
#include <string>

#include <utility> // std::pair
#include <iterator>
#include <sstream>

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

#include "parser.h"

using namespace OpenBabel; 

#define REASONABLE 1024

// --- DEV OPTIONS  ---
static bool opt_debug = false;

static void Fatal(const char *str){
  fprintf(stderr,"Fatal: %s\n",str);
  exit(1);
}

unsigned char static int_to_locant(unsigned int i){
  return i + 64;
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


static void print_locant_array(OBAtom **locant_path, unsigned int size){
  fprintf(stderr,"[ ");
  for(unsigned int i=0; i<size;i++){
    if(!locant_path[i])
      fprintf(stderr,"0 ");
    else
      fprintf(stderr,"%d ",locant_path[i]->GetIdx());
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
                          Locant Path Functions
**********************************************************************/

void copy_locant_path(OBAtom ** new_path,OBAtom **locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++)
    new_path[i] = locant_path[i]; 
}

unsigned int position_in_path(OBAtom *atom,OBAtom**locant_path,unsigned int path_size,bool error=true){
  for(unsigned int i=0;i<path_size;i++){
    if(atom == locant_path[i])
      return i; 
  }

  if(error)
    fprintf(stderr,"Error: atom not found in locant path\n");
  return 0; 
}

unsigned int fusion_locant(OBMol *mol,OBRing *ring, OBAtom **locant_path, unsigned int path_size){
  unsigned int lpos = path_size; 
  for(unsigned int i=0;i<ring->Size();i++){
    OBAtom *latom = mol->GetAtom(ring->_path[i]); 
    unsigned int pos = position_in_path(latom,locant_path,path_size,false); 
    if(pos < lpos)
      lpos = pos; 
  }
  return lpos; 
}

/* overall ring sum, combinates rule 30d and 30e for symmetrical structures */
unsigned int ring_sum(OBMol *mol, OBRing *ring, OBAtom**locant_path,unsigned int path_size){
  unsigned int rsum = 0;
  for(unsigned int i=0;i<ring->Size();i++){
    OBAtom *ratom = mol->GetAtom(ring->_path[i]);
    rsum += position_in_path(ratom, locant_path,path_size,false) + 1; 
  }
  return rsum;
}

unsigned int fusion_sum(OBMol *mol, OBAtom **locant_path, unsigned int path_size, std::set<OBRing*> &local_SSSR){
  unsigned int ret = 0; 
  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end();riter++)
    ret += fusion_locant(mol,*riter,locant_path,path_size) + 1; // A=1, B=2 etc 
  return ret;
}


void print_ring_locants(OBMol *mol,OBRing *ring, OBAtom **locant_path, unsigned int path_size, bool sort=false){
  unsigned char *sequence = (unsigned char*)malloc(sizeof(unsigned char)*ring->Size()); 
  
  for(unsigned int i=0;i<ring->Size();i++)
    sequence[i] = int_to_locant(position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)+1); 
  
  if(sort)
    sort_locants(sequence,ring->Size());
  fprintf(stderr,"[ ");
  for(unsigned int k=0;k<ring->Size();k++)
    fprintf(stderr,"%c ",sequence[k]);
  fprintf(stderr,"]\n");

  free(sequence);
}

bool sequential_chain(  OBMol *mol,OBRing *ring, 
                        OBAtom **locant_path, unsigned int path_size,
                        std::map<unsigned char, bool> &in_chain)
{
  unsigned char *sequence = (unsigned char*)malloc(sizeof(unsigned char)*ring->Size()); 
  for(unsigned int i=0;i<ring->Size();i++)
    sequence[i] = int_to_locant(position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)+1); 

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
    unsigned char ch = int_to_locant(i+1); 
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
        connections[int_to_locant(i+1)] = 4;
    }
  }

  std::map<std::set<unsigned char>,bool> seen_nt;
  std::vector<std::set<unsigned char>> non_trivials;

  for(unsigned int i=0;i<path_size;i++){
    for(unsigned int j=i+2;j<path_size;j++){
      if(mol->GetBond(locant_path[i],locant_path[j]))
        non_trivials.push_back({int_to_locant(i+1),int_to_locant(j+1)});
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
      else if(locant < int_to_locant(path_size))
        locant++;

      path.push_back(locant);
    }

    // add the loop back logic
    for(unsigned int i=0;i<path.size();i++){
      unsigned int tally = 1;
      if(path[i] == int_to_locant(path_size)){
        for(unsigned int j=i+1;j<path.size();j++){
          if(path[j] == path[i]){
            path[j] += -tally; // looping back
            tally++;
          }
        }
      }
    }

    while(!connections[bind] && bind < int_to_locant(path_size)){
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


/* read locant path algorithm, we return the number of non consecutive blocks,
pseudo check will add determined pairs and check notation is viable for read */
unsigned int ReadLocantPath(  OBMol *mol, OBAtom **locant_path, unsigned int path_size,
                      std::set<OBRing*>               &local_SSSR,
                      std::map<OBAtom*,bool>          &bridge_atoms,
                      std::map<OBAtom*,OBAtom*>       &broken_atoms,
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
            unsigned char loc = int_to_locant(position_in_path(mol->GetAtom(wring->_path[k]),locant_path,path_size)+1); 
            if(loc < min_loc)
              min_loc = loc; 
            if(loc > high_loc)
              high_loc = loc;
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
      unsigned char loc = int_to_locant(position_in_path(mol->GetAtom(to_write->_path[k]),locant_path,path_size)+1); 
      in_chain[loc] = true;
    }

    if(opt_debug && verbose){
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
  
    if(to_write->Size() > 9){
      buffer+='-';
      buffer+= std::to_string(to_write->Size());
      buffer+='-';
    }
    else
      buffer+= std::to_string(to_write->Size());

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


OBAtom **MonoPath(OBMol *mol, unsigned int path_size,
                  std::set<OBRing*> &local_SSSR)
{
  OBAtom **locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
  for(unsigned int i=0;i<path_size;i++)
    locant_path[i] = 0; 

  OBRing *mono = *(local_SSSR.begin());

  for(unsigned int i=0;i<mono->Size();i++)
  locant_path[i] = mol->GetAtom(mono->_path[i]);
    
  return locant_path;
}

/*  standard ring walk, can deal with all standard polycyclics without an NP-Hard
    solution, fusion sum is the only filter rule needed here
*/
OBAtom **PLocantPath(   OBMol *mol, unsigned int path_size,
                        std::set<OBAtom*>               &ring_atoms,
                        std::set<OBBond*>               &ring_bonds,
                        std::map<OBAtom*,unsigned int>  &atom_shares,
                        std::map<OBAtom*,bool>          &bridge_atoms,
                        std::map<OBAtom*,OBAtom*>       &broken_atoms,
                        std::set<OBRing*>               &local_SSSR)
{

  // create the path
  OBAtom **locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
  OBAtom **best_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
  for(unsigned int i=0;i<path_size;i++){
    locant_path[i] = 0;
    best_path[i] = 0; 
  }

  // set up some non-trivial bonds
  std::vector<OBBond*> nt_bonds; 
  for(std::set<OBBond*>::iterator biter = ring_bonds.begin(); biter != ring_bonds.end(); biter++){
    OBBond *bond = (*biter);
    unsigned int share = 0; 
    for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end(); riter++){
      OBRing *obring = (*riter); 
      if(obring->IsMember(bond))
        share++; 
    }

    if(share > 1)
      nt_bonds.push_back(bond);   
  }

  unsigned int           lowest_sum       = UINT32_MAX;
  unsigned int           lowest_score       = UINT32_MAX;

  OBAtom*                ratom  = 0;
  OBAtom*                catom  = 0;
  OBBond*                bond   = 0; 
  std::map<OBBond*,bool> ignore_bond; 

  for(unsigned int i=0;i<nt_bonds.size();i++)
    ignore_bond[nt_bonds[i]] = true; 
  
  for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
    if(atom_shares[*aiter] == 2){
      std::stack<OBAtom*>    stack; 
      std::map<OBAtom*,bool> visited; 
      unsigned int locant_pos = 0;
      stack.push(*aiter);

      while(!stack.empty()){
        ratom = stack.top();
        stack.pop();
        locant_path[locant_pos++] = ratom; 
        visited[ratom] = true; 

        FOR_NBORS_OF_ATOM(a,ratom){ 
          catom = &(*a);   
          bond = mol->GetBond(ratom,catom); 
          if(atom_shares[catom] && !visited[catom] && !ignore_bond[bond]){
            stack.push(catom);
            break;
          }
        }
      }

      std::vector<OBRing*> tmp; 
      std::string candidate_string; // super annoying this has to go here, UPSET
      unsigned int score = ReadLocantPath(mol,locant_path,path_size,local_SSSR,bridge_atoms,broken_atoms,tmp,candidate_string,false);
      unsigned int fsum = fusion_sum(mol,locant_path,path_size,local_SSSR);

#define EXPERIMENT 0
#if EXPERIMENT
      fprintf(stderr,"%s - score: %d, fsum: %d\n",candidate_string.c_str(),score,fsum);
#endif
      
      if(score < lowest_score){
        lowest_sum = fsum;
        lowest_score = score; 
        copy_locant_path(best_path,locant_path,path_size);
      }
      else if (score == lowest_score){
        if(fsum < lowest_sum){ // rule 30d.
          lowest_sum = fsum;
          copy_locant_path(best_path,locant_path,path_size);
        }
      }
    }
  }

  free(locant_path);
  for(unsigned int i=0;i<path_size;i++){
    if(!best_path[i]){
      free(best_path);
      Fatal("no continous locant path was possible - currently unsupported - PAlgorithm\n");
    }
  }

  return best_path; 
}


/* uses a flood fill style solution (likely NP-HARD), with some restrictions to 
find a multicyclic path thats stable with disjoined pericyclic points */
OBAtom **NPLocantPath(      OBMol *mol, unsigned int path_size,
                            std::set<OBAtom*>               &ring_atoms,
                            std::map<OBAtom*,unsigned int>  &atom_shares,
                            std::map<OBAtom*,bool>          &bridge_atoms, // rule 30f.
                            std::map<OBAtom*,OBAtom*>       &broken_atoms,
                            std::set<OBRing*>               &local_SSSR)
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

  // decend logic idea from Tom Allam
  while(!path_found && found_path_size > 1){ 
    
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
          if(atom_shares[catom]){
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
            unsigned int score = ReadLocantPath(mol,locant_path,found_path_size,local_SSSR,bridge_atoms,broken_atoms,tmp,candidate_string,false);
            unsigned int fsum = fusion_sum(mol,locant_path,found_path_size,local_SSSR);
           
#if EXPERIMENT
      fprintf(stderr,"%s - score: %d, fsum: %d\n",candidate_string.c_str(),score,fsum);
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
        if(safety == 1000000)
          break;
        
          
      }
    }

    for(unsigned int i=0;i<path_size;i++)
      locant_path[i] = 0;
    
    if(!path_found){
      for(unsigned int i=0;i<path_size;i++)
        best_path[i] = 0;

      found_path_size--; // decrement the path size and see what we can do
    }
  }
  
  free(locant_path);
  if(!path_found){
    free(best_path);
    Fatal("no locant path could be generated, even with decrements - first stop point\n");
  }

  if(found_path_size < path_size){

    if(opt_debug)
      fprintf(stderr,"  found locant path with %d branches out\n",path_size-found_path_size);

    for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end();aiter++){
      if(position_in_path(*aiter,best_path,found_path_size,false) == 0 && best_path[0] != *aiter){
        OBAtom *branching = *aiter; 
        unsigned int lowest_pos = found_path_size; 
        FOR_NBORS_OF_ATOM(a,branching){
          unsigned int bpos = position_in_path(&(*a),best_path,found_path_size);
          if(bpos < lowest_pos)
            lowest_pos = bpos; 
        }

        if(opt_debug)
          fprintf(stderr,"  branching is bonded to: %c\n",int_to_locant(lowest_pos+1));
        
        if(broken_atoms[locant_path[lowest_pos]])
          Fatal("multiple broken locants per atom are currently unsupported");
        else
          broken_atoms[locant_path[lowest_pos]] = branching;
        
        if(found_path_size < path_size)
          best_path[found_path_size++] = branching; // stick at end of path
      }
    }
  }

  if(found_path_size < path_size){
    free(best_path);
    Fatal("no locant path could be generated, even with decrements - second stop point\n");
  }
  return best_path;
}


bool IsHeteroRing(OBAtom **locant_array,unsigned int size){
  for(unsigned int i=0;i<size;i++){
    if(locant_array[i]->GetAtomicNum() != 6)
      return true;
  }
  return false; 
}


unsigned int ClassifyRing(  std::set<OBAtom*>               &ring_atoms,
                            std::map<OBAtom*,unsigned int>  &atom_shares)
{
  unsigned int classification = 0;
  for(std::set<OBAtom*>::iterator iter = ring_atoms.begin();iter != ring_atoms.end();iter++){
    if(atom_shares[*iter] == 3){
      if(classification < 1)
        classification = 1; 
    }

    if(atom_shares[*iter] == 4)
      return 2;
  }
  return classification; 
}



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

  std::map<OBAtom*,bool> atoms_seen;
  std::map<OBRing*,bool> rings_seen; 
  std::map<OBAtom*,int>  remaining_branches; // tracking for branch pop
  std::map<OBAtom*,unsigned int> string_position; // essential for writing post charges. 

  // recursion tracking
  unsigned int cycle_count; 
  unsigned int last_cycle_seen;
  
  BabelGraph(){
    cycle_count = 0; 
    last_cycle_seen = 0; 
  };
  ~BabelGraph(){};

  unsigned char WriteSingleChar(OBAtom* atom){

    if(!atom)
      Fatal("writing notation from dead atom ptr");
    
    unsigned int neighbours = atom->GetExplicitDegree(); 
    unsigned int orders = atom->GetExplicitValence(); 

    switch(atom->GetAtomicNum()){
      case 1:
        return 'H';

      case 5:
        if(neighbours > 3)
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
        if(atom->GetFormalCharge() == +1)
          return 'K';
        
        if(orders == 0 || orders == 1)
            return 'Z'; 
        else if(orders == 2)
            return 'M';
        else if(orders == 3)
          return 'N';
        else 
          return 'K';
      
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

  void WriteSpecial(OBAtom *atom, std::string &buffer){

    if(!atom)
      Fatal("writing notation from dead atom ptr");
    // all special elemental cases
    switch(atom->GetAtomicNum()){
      case 5:
        buffer += "-B-";
        break;

      case 8:
        buffer += "-O-";
        break;

      case 9:
        buffer += "-F-";
        break;

      case 53:
        buffer += "-I-";
        break;

      case 35:
        buffer += "-E-";
        break;

      case 17:
        buffer += "-G-";
        break; 

      case 89:
        buffer += "-AC-";
        break;

      case 47:
        buffer += "-AG-";
        break;
    
      case 13:
        buffer += "-AL-";
        break;

      case 95:
        buffer += "-AM-";
        break;

      case 18:
        buffer += "-AR-";
        break;

      case 33:
        buffer += "-AS-";
        break;

      case 85:
        buffer += "-AT-";
        break;

      case 79:
        buffer += "-AU-";
        break;


      case 56:
        buffer += "-BA-";
        break;

      case 4:
        buffer += "-BE-";
        break;

      case 107:
        buffer += "-BH-";
        break;

      case 83:
        buffer += "-BI-";
        break;

      case 97:
        buffer += "-BK-";
        break;

      case 20:
        buffer += "-CA-";
        break;
      
      case 48:
        buffer += "-CD-";
        break;

      case 58:
        buffer += "-CE-";
        break;

      case 98:
        buffer += "-CF-";
        break;

      case 96:
        buffer += "-CN-";
        break;

      case 112:
        buffer += "-CN-";
        break;

      case 27:
        buffer += "-CO-";
        break;

      case 24:
        buffer += "-CR-";
        break;

      case 55:
        buffer += "-CS-";
        break;

      case 29:
        buffer += "-CU-";
        break;

      case 105:
        buffer += "-DB-";
        break;

      case 110:
        buffer += "-DS-";
        break;

      case 66:
        buffer += "-DY-";
        break;

      case 68:
        buffer += "-ER-";
        break;

      case 99:
        buffer += "-ES-";
        break;

      case 63:
        buffer += "-EU-";
        break;

      case 26:
        buffer += "-FE-";
        break;

      case 114:
        buffer += "-FL-";
        break;

      case 100:
        buffer += "-FM-";
        break;

      case 87:
        buffer += "-FR-";
        break;

      case 31:
        buffer += "-GA-";
        break;

      case 64:
        buffer += "-GD-";
        break;

      case 32:
        buffer += "-GE-";
        break;

      case 2:
        buffer += "-HE-";
        break;

      case 72:
        buffer += "-HF-";
        break;

      case 80:
        buffer += "-HG-";
        break;

      case 67:
        buffer += "-HO-";
        break;

      case 108:
        buffer += "-HS-";
        break;

      case 49:
        buffer += "-IN-";
        break;

      case 77:
        buffer += "-IR-";
        break;

      case 36:
        buffer += "-KR-";
        break;

      case 19:
        buffer += "-KA-";
        break;

      case 57:
        buffer += "-LA-";
        break;

      case 3:
        buffer += "-LI-";
        break;

      case 103:
        buffer += "-LR-";
        break;

      case 71:
        buffer += "-LU-";
        break;

      case 116:
        buffer += "-LV-";
        break;

      case 115:
        buffer += "-MC-";
        break;

      case 101:
        buffer += "-MD-";
        break;

      case 12:
        buffer += "-MG-";
        break;

      case 25:
        buffer += "-MN-";
        break;

      case 42:
        buffer += "-MO-";
        break;

      case 109:
        buffer += "-MT-";
        break;

      case 11:
        buffer += "-NA-";
        break;

      case 41:
        buffer += "-NB-";
        break;

      case 60:
        buffer += "-ND-";
        break;

      case 10:
        buffer += "-NE-";
        break;

      case 113:
        buffer += "-NH-";
        break;

      case 28:
        buffer += "-NI-";
        break;

      case 102:
        buffer += "-NO-";
        break;

      case 93:
        buffer += "-NP-";
        break;

      case 118:
        buffer += "-OG-";
        break;

      case 76:
        buffer += "-OS-";
        break;

      case 91:
        buffer += "-PA-";
        break;

      case 82:
        buffer += "-PB-";
        break;

      case 46:
        buffer += "-PD-";
        break;

      case 61:
        buffer += "-PM-";
        break;

      case 84:
        buffer += "-PO-";
        break;

      case 59:
        buffer += "-PR-";
        break;

      case 78:
        buffer += "-PT-";
        break;

      case 94:
        buffer += "-PU-";
        break;

      case 88:
        buffer += "-RA-";
        break;

      case 37:
        buffer += "-RB-";
        break;

      case 75:
        buffer += "-RE-";
        break;

      case 104:
        buffer += "-RF-";
        break;

      case 111:
        buffer += "-RG-";
        break;

      case 45:
        buffer += "-RH-";
        break;

      case 86:
        buffer += "-RN-";
        break;

      case 44:
        buffer += "-RU-";
        break;

      case 51:
        buffer += "-SB-";
        break;

      case 21:
        buffer += "-SC-";
        break;

      case 34:
        buffer += "-SE-";
        break;

      case 106:
        buffer += "-SG-";
        break;

      case 14:
        buffer += "-SI-";
        break;

      case 62:
        buffer += "-SM-";
        break;

      case 50:
        buffer += "-SN-";
        break;

      case 38:
        buffer += "-SR-";
        break;


      case 73:
        buffer += "-TA-";
        break;

      case 65:
        buffer += "-TB-";
        break;

      case 43:
        buffer += "-TC-";
        break;

      case 52:
        buffer += "-TE-";
        break;

      case 90:
        buffer += "-TH-";
        break;

      case 22:
        buffer += "-TI-";
        break;

      case 81:
        buffer += "-TL-";
        break;

      case 69:
        buffer += "-TM-";
        break;

      case 117:
        buffer += "-TS-";
        break;

      case 92:
        buffer += "-UR-";
        break;

      case 23:
        buffer += "-VA-";
        break;

      case 54:
        buffer += "-XE-";
        break;

      case 39:
        buffer += "-YT-";
        break;

      case 70:
        buffer += "-YB-";
        break;

      case 30:
        buffer += "-ZN-";
        break;

      case 40:
        buffer += "-ZR-";
        break;
    }
  }

  unsigned int CountDioxo(OBAtom *atom){
    if(!atom)
      Fatal("count dioxo on dead atom ptr");

    unsigned int Ws = 0; 
    unsigned int carbonyls = 0;
    unsigned int oxo_ions = 0; 
    std::vector<OBAtom*> seen; 
    FOR_NBORS_OF_ATOM(a,atom){
      OBAtom *nbor = &(*a);
      if(!atoms_seen[nbor] && !nbor->IsInRing() && nbor->GetAtomicNum() == 8){
 
        if(atom->GetBond(nbor)->GetBondOrder() == 2){
          carbonyls++;
          seen.push_back(nbor);
        }
        else if(nbor->GetFormalCharge() == -1){
          oxo_ions++;
          seen.push_back(nbor);
        }
          
        if(carbonyls == 2 || (oxo_ions == 1 && carbonyls == 1)){
          Ws++;
          atoms_seen[seen[0]] = true;
          atoms_seen[seen[1]] = true;

          for(OBAtom *pcharge : seen){
            if(pcharge->GetFormalCharge() == -1)
              pcharge->SetFormalCharge(0); // remove the charge, as this is expected 
          }

          carbonyls = 0;
          oxo_ions = 0;
          seen.clear();
        }
      }  
    }
    return Ws;
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
  bool ParseNonCyclic(OBAtom* start_atom, OBAtom *spawned_from, unsigned int b_order,
                      OBMol *mol, std::string &buffer, 
                      unsigned int cycle_num, unsigned char locant){
    if(!start_atom)
      Fatal("writing notation from dead atom ptr");
    //##################################
    //      INDIRECT RECURSION TRACKING

    if(last_cycle_seen > cycle_num){
      for(unsigned int i=0;i<(last_cycle_seen-cycle_num);i++){
        buffer+='&';
        if(cycle_count)
          cycle_count--; // once a ring is closed can you ever get back? - GOOD
      }
    }
    last_cycle_seen = cycle_num;

    if(locant && locant != '0' && b_order > 0){ // allows OM through
      buffer+=' ';
      write_locant(locant,buffer);
    }

    for(unsigned int b=1;b<b_order;b++)
      buffer+='U';
    //##################################

    unsigned int carbon_chain = 0;

    OBAtom* atom = start_atom;
    OBAtom* prev = 0; 
    OBBond *bond = 0; 

    std::stack<OBAtom*> atom_stack; 
    std::stack<OBAtom*> branch_stack;
    std::map<OBAtom*,bool> branching_atom; 
    atom_stack.push(atom); 

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

          if(!branching_atom[prev])
            buffer += '&'; // requires a closure 
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

        if(!bond)
          Fatal("failure to read branched bond segment");

        remaining_branches[prev]--; // reduce the branches remaining  
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

        cycle_count++;
        if(locant == '0' && b_order == 0){
          buffer += '-';
          buffer += ' ';
          buffer += '0';
          if(!RecursiveParse(atom,spawned_from,mol,false,buffer,cycle_count))
            Fatal("failed to make pi bonded ring");
        }
        else{
          if(!RecursiveParse(atom,spawned_from,mol,true,buffer,cycle_count))
            Fatal("failed to make inline ring");
        }

        if(!atom_stack.empty()){

          if(last_cycle_seen > cycle_num){
            for(unsigned int i=0;i<(last_cycle_seen-cycle_num);i++){
              buffer+='&';
              if(cycle_count)
                cycle_count--; // once a ring is closed can you ever get back? - GOOD
            }
          }

          last_cycle_seen = cycle_count;
          if(!branch_stack.empty())
            prev = return_open_branch(branch_stack);
        }
        continue;
      }

       // remaining_branches are -1, we only look forward
      unsigned int correction = 0; 
      unsigned char wln_character =  WriteSingleChar(atom);
      unsigned int Wgroups = CountDioxo(atom);

    
      if(prev && bond)
        correction = bond->GetBondOrder() - 1;
      else if (b_order > 0)
        correction = b_order - 1;

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
          buffer += wln_character; 
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
          if(!Wgroups && CheckCarbonyl(atom))
            buffer += 'V';
          else{
            buffer += wln_character;
            if(wln_character == 'X')
              remaining_branches[atom] = 3;
            else
              remaining_branches[atom] = 2; 

            branching_atom[atom] = true;
            branch_stack.push(atom);
          }
          string_position[atom] = buffer.size();
          break;

        case 'N':
        case 'B':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom; 
          buffer += wln_character;

          if(atom->GetTotalDegree() > 1){
            remaining_branches[atom] = 2 - correction; 
            branch_stack.push(atom);
          }
          string_position[atom] = buffer.size();
          break;

        case 'K':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom;
          buffer += wln_character;

          // K now given for all positive nitrogen
          string_position[atom] = buffer.size();
          if(atom->GetExplicitValence() < 4){
            for(unsigned int i=atom->GetExplicitValence();i<4;i++){
              buffer += 'H';
              correction++;
            }
          }

          if(atom->GetTotalDegree() > 1){
            remaining_branches[atom] = 3 - correction;
            branching_atom[atom] = true; 
            branch_stack.push(atom);
          }
          atom->SetFormalCharge(0); // remove the charge, as this is expected 
          break;
        
        case 'P':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom;
          buffer += wln_character;
          string_position[atom] = buffer.size();
          if(atom->GetExplicitValence() < 2)
            buffer += 'H';
          
          if(atom->GetTotalDegree() > 1){
            remaining_branches[atom] = 4 - correction; 
            branching_atom[atom] = true;
            branch_stack.push(atom);
          }
          break;

        case 'S':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom;
          buffer += wln_character;
          string_position[atom] = buffer.size();
          if(atom->GetExplicitValence() < 2)
            buffer += 'H';

          if(atom->GetTotalDegree() > 1){
            remaining_branches[atom] = 5 - correction; 
            branching_atom[atom] = true;
            branch_stack.push(atom);
          }          
          break;

        case '*':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          prev = atom; 
          string_position[atom] = buffer.size()+2;          
          WriteSpecial(atom,buffer);
          if(atom->GetTotalDegree() > 1){
            remaining_branches[atom] = 5 - correction; 
            branching_atom[atom] = true;
            branch_stack.push(atom);
          }

          for(unsigned int i=0;i<atom->GetImplicitHCount();i++)
            buffer += 'H';
          
          break;
          
// terminators
        case 'Q':
        case 'Z':
        case 'E':
        case 'F':
        case 'G':
        case 'I':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          buffer += wln_character; 
          string_position[atom] = buffer.size();
          if(atom->GetExplicitValence() == 0 && atom->GetFormalCharge() == 0)
            buffer += 'H';

          if(!branch_stack.empty())
            prev = return_open_branch(branch_stack);

          if(atom->GetExplicitDegree() == 0)
            atom->SetFormalCharge(0); // notational implied, do not write ionic code

          break;

        case 'H':
          if(carbon_chain){
            buffer += std::to_string(carbon_chain);
            carbon_chain = 0;
          }

          buffer += wln_character;
          if(atom->GetExplicitValence() == 0 && atom->GetFormalCharge() == 0)
            buffer += 'H';

          string_position[atom] = buffer.size();
          break;

        
        default:
          fprintf(stderr,"Error: unhandled char %c\n",wln_character); 
          return 0; 
      }

      if(Wgroups){
        for(unsigned int i=0;i<Wgroups;i++){
          buffer+='W';
          remaining_branches[atom] += -3;
        }
        OBAtom *ret = return_open_branch(branch_stack);
        if(ret)
          prev = ret; 
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

    return atom; 
  }

  void AddPostCharges(OBMol *mol,std::string &buffer){
    if(opt_debug)
      fprintf(stderr,"Post Charges\n");
    
    bool working = true;
    while(working){
      working = false;
      FOR_ATOMS_OF_MOL(a,mol){
        OBAtom *atom = &(*a);
        if(atom->GetFormalCharge() != 0){

          if(opt_debug)
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
  unsigned int ConstructLocalSSSR(  OBMol *mol, OBAtom *ring_root,
                                    std::set<OBAtom*>         &ring_atoms,
                                    std::set<OBBond*>         &ring_bonds,
                                    std::map<OBAtom*,bool>    &bridge_atoms,
                                    std::map<OBAtom*,unsigned int> &atom_shares,
                                    std::set<OBRing*> &local_SSSR)
  {

    if(!ring_root){
      fprintf(stderr,"Error: ring root is nullptr\n");
      return 0;
    }

    OBAtom *ratom = 0; 
    OBAtom *prev = 0; 
    OBBond *bond = 0; 
    OBRing *obring = 0; 
    std::set<OBAtom*> tmp_bridging_atoms;

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

          for(unsigned int i=0;i<obring->Size();i++){
            OBAtom *ratom = mol->GetAtom(obring->_path[i]);
            ring_set.insert(ratom);
          }

          std::set_intersection(ring_set.begin(), ring_set.end(), ring_atoms.begin(), ring_atoms.end(),
                                std::inserter(intersection, intersection.begin()));

          // intersection == 1 is a spiro ring, ignore,
          if(intersection.size() > 1){
            rings_seen[obring] = true; 
            local_SSSR.insert(obring);
            prev = 0;

            // if its enough to say that true bridges cannot have more than two bonds each?
            // yes but 2 bonds within the completed local SSSR,so this will needed filtering
            if(intersection.size() > 2){
              for(std::set<OBAtom*>::iterator iiter = intersection.begin(); iiter != intersection.end();iiter++)
                tmp_bridging_atoms.insert(*iiter); 
            }

            for(unsigned int i=0;i<obring->Size();i++){
              OBAtom *ratom = mol->GetAtom(obring->_path[i]);
              ring_atoms.insert(ratom);
              atom_shares[ratom]++;

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

    if(opt_debug){
      fprintf(stderr,"  ring atoms: %lu\n",ring_atoms.size());
      fprintf(stderr,"  ring bonds: %lu\n",ring_bonds.size());
      fprintf(stderr,"  ring subcycles: %lu/%lu\n",local_SSSR.size(),mol->GetSSSR().size());
      if(bridge_count)
        fprintf(stderr,"  bridging atoms: %d\n",bridge_count);
      
    }
      
    return ring_atoms.size(); 
  }


  /* create the heteroatoms and locant path unsaturations where neccesary */
  bool ReadLocantAtomsBonds(  OBMol *mol, OBAtom** locant_path,unsigned int path_size,
                              std::vector<OBRing*> &ring_order,
                              std::set<OBBond*>   &ring_bonds,
                              std::string &buffer)
  {

    
    unsigned char locant = 0;
    unsigned char last_locant = 'A'; 
    std::map<OBBond*,bool> bonds_checked; 

    if(!std::isdigit(buffer.back()))
      last_locant = ' ';
    

    for(unsigned int i=0;i<path_size;i++){
      if(!locant_path[i])
        Fatal("dead locant path atom ptr");

      unsigned char het_char = 0;
      locant = int_to_locant(i+1);
      unsigned int Wgroups = CountDioxo(locant_path[i]);
      bool carbonyl = CheckCarbonyl(locant_path[i]);

      if( !carbonyl && !Wgroups && 
        locant_path[i]->GetAtomicNum() == 6 &&
        locant_path[i]->GetFormalCharge() == -1){
        // organometallics logic 
        if(locant != last_locant){
          buffer += ' ';
          write_locant(locant,buffer);
          last_locant = locant;
        } 

        buffer += '0';
        locant_path[i]->SetFormalCharge(0);
      }

      if(carbonyl || Wgroups || locant_path[i]->GetAtomicNum() != 6){
        if(locant != last_locant){
          buffer += ' ';
          write_locant(locant,buffer);
          last_locant = locant;
        } 
        if(Wgroups){
          het_char = WriteSingleChar(locant_path[i]);
          
          if(het_char != '*'){
            if(het_char == 'K')
              locant_path[i]->SetFormalCharge(0);

            buffer+=het_char; 
            string_position[locant_path[i]] = buffer.size();
          }else{
            WriteSpecial(locant_path[i],buffer); 
            string_position[locant_path[i]] = buffer.size()+2;
          }
          for(unsigned int w=0;w<Wgroups;w++)
            buffer+='W';

          last_locant++;
        }
        else if(carbonyl){
          buffer += 'V';
          string_position[locant_path[i]] = buffer.size();
          last_locant++;
        }
        else{
          het_char = WriteSingleChar(locant_path[i]);
          if(het_char != '*'){
            if(het_char == 'K')
              locant_path[i]->SetFormalCharge(0);

            buffer+=het_char;
            string_position[locant_path[i]] = buffer.size(); 
          }
          else{
            WriteSpecial(locant_path[i],buffer); 
            string_position[locant_path[i]] = buffer.size()+2;
          }

          last_locant++;
        }
      }

      if(locant_path[i]->GetAtomicNum() == 6){
        unsigned int rbonds = 0; 
        for(unsigned int k=0;k<path_size;k++){
          if(mol->GetBond(locant_path[i],locant_path[k]))
            rbonds++;
        }
        if(rbonds == 4){
          if(locant != last_locant){
            buffer += ' ';
            write_locant(locant,buffer);
            last_locant = locant;
          } 
          buffer += 'X';
        }
      }

      // for now dont worry about positions, only that we get the double bonds
      // we can condense later
      bool bonds = false;
      OBAtom *first   = 0;
      OBAtom *second  = 0;
      // handles sequential locant unsaturations, when not aromatic
      first = locant_path[i];
      if(i < path_size-1)
        second = locant_path[i+1];
      else
        second = locant_path[0];

      OBBond *locant_bond = first->GetBond(second);
      
      if(locant_bond){
        bonds_checked[locant_bond] = true;
        for(unsigned int k=0;k<ring_order.size();k++){
          if(!ring_order[k]->IsAromatic() && ring_order[k]->IsMember(locant_bond))
            bonds = true; 
        }

        if(bonds && locant_bond && locant_bond->GetBondOrder() > 1){
          buffer += ' ';
          write_locant(locant,buffer);
          for(unsigned int b=1;b<locant_bond->GetBondOrder();b++)
            buffer += 'U';
        }
      }
    }

    for(std::set<OBBond*>::iterator biter = ring_bonds.begin(); biter != ring_bonds.end();biter++){
      OBBond *fbond = *biter; 
      if(!bonds_checked[fbond] && fbond->GetBondOrder() > 1){

        for(unsigned int k=0;k<ring_order.size();k++){
          if(!ring_order[k]->IsAromatic() && ring_order[k]->IsMember(fbond)){
            
            unsigned char floc = int_to_locant(position_in_path(fbond->GetBeginAtom(),locant_path,path_size)+1); 
            unsigned char bloc = int_to_locant(position_in_path(fbond->GetEndAtom(),locant_path,path_size)+1); 
            
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
      }
      bonds_checked[fbond] = true;
    }
    
    return true;
  }


  void ReadMultiCyclicPoints( OBAtom**locant_path,unsigned int path_size, 
                              std::map<OBAtom*,unsigned int> &ring_shares,std::string &buffer)
  { 
    unsigned int count = 0;
    std::string append = ""; 
    for(unsigned int i=0;i<path_size;i++){
      if(ring_shares[locant_path[i]] > 2){
        count++; 
        write_locant(int_to_locant(i+1),append);
      }
    }

    buffer += ' ';
    buffer+= std::to_string(count);
    buffer+= append; 
  }


  /* constructs and parses a cyclic structure, locant path is returned with its path_size */
  std::pair<OBAtom **,unsigned int> ParseCyclic(OBAtom *ring_root,OBAtom *spawned_from,OBMol *mol, bool inline_ring,std::string &buffer){
    if(opt_debug)
      fprintf(stderr,"Reading Cyclic\n");

    OBAtom **                       locant_path = 0; 
    std::set<OBRing*>               local_SSSR;
    std::set<OBAtom*>               ring_atoms;
    std::set<OBBond*>               ring_bonds;
    std::vector<OBRing*>            ring_order; 

    std::map<OBAtom*,bool>          bridge_atoms;
    std::map<OBAtom*,OBAtom*>       broken_atoms;
    std::map<OBAtom*,unsigned int>  atom_shares;
  
    unsigned int path_size   =  ConstructLocalSSSR(mol,ring_root,ring_atoms,ring_bonds,bridge_atoms,atom_shares,local_SSSR); 
    if(!path_size)
      Fatal("failed to write ring");


    bool multi = ClassifyRing(ring_atoms,atom_shares); 

    if(opt_debug)
      fprintf(stderr,"  multi classification: %d\n",multi);

    if(local_SSSR.size() == 1)
      locant_path = MonoPath(mol,path_size,local_SSSR);
    else if(!multi && bridge_atoms.empty())
      locant_path = PLocantPath(mol,path_size,ring_atoms,ring_bonds,atom_shares,bridge_atoms,broken_atoms,local_SSSR);
    else 
      locant_path = NPLocantPath(mol,path_size,ring_atoms,atom_shares,bridge_atoms,broken_atoms,local_SSSR);

    if(inline_ring){
      buffer+= '-';
      bool spiro = false;
      unsigned char root_locant = 0;
      for(unsigned int i=0;i<path_size;i++){
        if(locant_path[i] == ring_root)
          root_locant = int_to_locant(i+1);

        if(locant_path[i] == spawned_from){
          // must be spiro ring
          root_locant = int_to_locant(i+1);
          spiro = true;
          break;
        }
      }

      if(spiro)
        buffer += '&';
      
      buffer += ' ';
      buffer += root_locant;
    }
    
    if(IsHeteroRing(locant_path,path_size))
      buffer += 'T';
    else
      buffer += 'L';

    ReadLocantPath( mol,locant_path,path_size,
                    local_SSSR,bridge_atoms,broken_atoms,ring_order,
                    buffer,true); 
    
    if(!bridge_atoms.empty()){
      for(unsigned int i=0;i<path_size;i++){
        if(bridge_atoms[locant_path[i]]){
          buffer+= ' ';
          unsigned char bloc = int_to_locant(i+1);
          write_locant(bloc,buffer);
        }
      }
    }

    if(multi){
      ReadMultiCyclicPoints(locant_path,path_size,atom_shares,buffer);
      buffer += ' ';
      write_locant(int_to_locant(path_size),buffer); // need to make the relative size
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
    return {locant_path,path_size};
  }
    

  bool RecursiveParse(OBAtom *atom, OBAtom *spawned_from, OBMol *mol, bool inline_ring,std::string &buffer, unsigned int cycle_num){
    // assumes atom is a ring atom 
    last_cycle_seen = cycle_num;
    std::pair<OBAtom**,unsigned int> path_pair; 

    path_pair = ParseCyclic(atom,spawned_from,mol,inline_ring,buffer);

    if(!path_pair.first){
      fprintf(stderr,"Error: failed on cyclic parse\n");
      return false;
    }

    for(unsigned int i=0;i<path_pair.second;i++)
      atoms_seen[path_pair.first[i]] = true;
      
    for(unsigned int i=0;i<path_pair.second;i++){
      FOR_NBORS_OF_ATOM(iter,path_pair.first[i]){
        OBAtom *latom = &(*iter);
        OBBond* lbond = path_pair.first[i]->GetBond(latom);
        if(!atoms_seen[latom]){
          if(!ParseNonCyclic( latom,path_pair.first[i],lbond->GetBondOrder(),
                              mol,buffer,
                              cycle_num,int_to_locant(i+1))){
            fprintf(stderr,"Error: failed on non-cyclic parse\n");
            return false;
          }

        }
      }

      // OM logic 
      if(path_pair.first[i]->GetAtomicNum() == 6 && path_pair.first[i]->GetFormalCharge() == -1){
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
              path_pair.first[i]->SetFormalCharge(0);
              if(charge)
                charge--;

              // find and write the other rings based on the negative charges
              FOR_ATOMS_OF_MOL(negc,mol){
                OBAtom* next_pi =  &(*negc); 
                if(!atoms_seen[next_pi] && next_pi->GetAtomicNum() == 6
                    && next_pi->GetFormalCharge() == -1 && next_pi->IsInRing()){
                  if(!ParseNonCyclic(next_pi,path_pair.first[i],0,
                              mol,buffer,
                              cycle_num,'0')){
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
    
    free(path_pair.first);
    return true;
  }


};



/**********************************************************************
                         API FUNCTION
**********************************************************************/

bool WriteWLN(std::string &buffer, OBMol* mol)
{   
  
  OBMol *mol_copy = new OBMol(*mol); // performs manipulations on the mol object, copy for safety

  BabelGraph obabel; 
  unsigned int cyclic = 0;
  bool started = false; 
  FOR_RINGS_OF_MOL(r,mol_copy)
    cyclic++;

  if(opt_debug)
    WriteBabelDotGraph(mol_copy);

  if(!cyclic){
    
    FOR_ATOMS_OF_MOL(a,mol_copy){
      OBAtom *satom = &(*a); 
      if(!obabel.atoms_seen[satom] && (satom->GetExplicitDegree()==1 || satom->GetExplicitDegree() == 0) ){
        if(started)
          buffer += " &"; // ionic species
        if(!obabel.ParseNonCyclic(&(*a),0,0,mol_copy,buffer,0,0))
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
          obabel.cycle_count = 0;
          obabel.last_cycle_seen = 0;
        }
        
        if(!obabel.RecursiveParse(mol_copy->GetAtom( (&(*r))->_path[0]),0,mol_copy,false,buffer,0))
          Fatal("failed on recursive ring parse");

        started = true;
      }
    }
    
    // handles additional ionic atoms here
    obabel.cycle_count = 0;
    obabel.last_cycle_seen = 0;
    FOR_ATOMS_OF_MOL(a,mol_copy){
      OBAtom *satom = &(*a); 
      if(!obabel.atoms_seen[satom] && (satom->GetExplicitDegree()==1 || satom->GetExplicitDegree() == 0) ){
        buffer += " &"; // ionic species
        if(!obabel.ParseNonCyclic(satom,0,0,mol_copy,buffer,0,0))
          Fatal("failed on recursive branch parse");
      }
    }
  }

  obabel.AddPostCharges(mol_copy,buffer); // add in charges where we can 

  delete mol_copy; 
  return true; 
}



