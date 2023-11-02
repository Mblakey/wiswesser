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

using namespace OpenBabel; 

#define REASONABLE 1024

const char *cli_inp;
const char *format; 

// --- options ---
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

void shift_locants_left(OBAtom **locant_path,unsigned int path_size){
  OBAtom *tmp = locant_path[0];
  for(unsigned int i=0;i<path_size-1;i++)
    locant_path[i] = locant_path[i+1]; 
  locant_path[path_size-1] = tmp; 
}

void copy_locant_path(OBAtom ** new_path,OBAtom **locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++)
    new_path[i] = locant_path[i]; 
}

unsigned int position_in_path(OBAtom *atom,OBAtom**locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++){
    if(atom == locant_path[i])
      return i; 
  }
  fprintf(stderr,"Error: atom not found in locant path\n");
  return 0; 
}


void print_ring_locants(OBMol *mol,OBRing *ring, OBAtom **locant_path, unsigned int path_size, bool sort=false){
  unsigned char *sequence = (unsigned char*)malloc(sizeof(unsigned char)*ring->Size()); 
  for(unsigned int i=0;i<ring->Size();i++)
    sequence[i] = int_to_locant(position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)+1); 
  
  if(sort)
    sort_locants(sequence,ring->Size());
  
  for(unsigned int k=0;k<ring->Size();k++)
    fprintf(stderr,"%c ",sequence[k]); 
  
  free(sequence);
}

bool BondedEndPoints(   OBMol *mol,OBRing *ring, 
                        OBAtom **locant_path, unsigned int path_size,
                        std::set<OBAtom*>              &bridge_atoms,
                        std::map<OBAtom*,unsigned int> &atom_shares,
                        std::string &buffer){

  unsigned char *sequence = (unsigned char*)malloc(sizeof(unsigned char)*ring->Size()); 
  OBAtom **ring_array = (OBAtom**)malloc(sizeof(OBAtom*)*ring->Size()); 
  for(unsigned int i=0;i<ring->Size();i++){
    sequence[i] = int_to_locant(position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)+1); 
    ring_array[i] = mol->GetAtom(ring->_path[i]); 
  }

  for (unsigned int j=1;j<ring->Size();j++){
		unsigned char key = sequence[j];
    OBAtom *ptr = ring_array[j];
		int i = j-1;
		while(i>=0 && sequence[i] > key){
			sequence[i+1] = sequence[i];
      ring_array[i+1] = ring_array[i];
			i--;
		}
		sequence[i+1] = key;
    ring_array[i+1] = ptr;
	}

  // need a new approach
  bool sequential_read = true;
  unsigned int last = ring->Size()-1;


  // for(unsigned int i=0;i<ring->Size();i++)
  //   fprintf(stderr,"%c ",sequence[i]);
  // fprintf(stderr,"\n");

  unsigned char bind = 0;
  if(!mol->GetBond(ring_array[0],ring_array[last]) && !(bridge_atoms.find(ring_array[last]) != bridge_atoms.end()) ){
    unsigned char loc = sequence[0];
    for(unsigned int i=1;i<last;i++){
      if(sequence[i] != loc+1)
        sequential_read = false;
      
      if(mol->GetBond(ring_array[i],ring_array[last])){
        bind = sequence[i];
        break;
      }
        
      
  
      loc = sequence[i];
    }
    if(!sequential_read){
      buffer += '/';
      write_locant(bind,buffer);
      write_locant(sequence[last],buffer);
    }
  }
 
  free(ring_array);
  free(sequence);
  return true;
}


/*  standard ring walk, can deal with all standard polycyclics without an NP-Hard
    solution 
*/
OBAtom **PLocantPath(   OBMol *mol, unsigned int path_size,
                        std::set<OBAtom*>              &ring_atoms,
                        std::set<OBBond*>              &ring_bonds,
                        std::map<OBAtom*,unsigned int> &atom_shares,
                        std::set<OBRing*>              &local_SSSR)
{

  // create the path
  unsigned int locant_pos = 0; 
  OBAtom **locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
  for(unsigned int i=0;i<path_size;i++)
    locant_path[i] = 0;

  // choose a seed with the highest degree of ring shares
  OBAtom *rseed = 0; 
  for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
    if(!rseed)
      rseed = (*aiter);
    else if(atom_shares[(*aiter)] > atom_shares[rseed])
      rseed = (*aiter);
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

  OBAtom*                ratom  = 0;
  OBAtom*                catom  = 0;
  OBBond*                bond   = 0; 
  std::map<OBAtom*,bool> atoms_in_lp; 
  std::map<OBBond*,bool> ignore_bond; 
  std::stack<OBAtom*>    stack; 

  for(unsigned int i=0;i<nt_bonds.size();i++){
    if(opt_debug)
      fprintf(stderr,"  fused ring junction: %d --> %-d\n",nt_bonds[i]->GetBeginAtomIdx(),nt_bonds[i]->GetEndAtomIdx());
    ignore_bond[nt_bonds[i]] = true; 
  } 
  
  stack.push(rseed);
  while(!stack.empty()){
    ratom = stack.top();
    stack.pop();
    locant_path[locant_pos++] = ratom; 
    atoms_in_lp[ratom] = true; 

    FOR_NBORS_OF_ATOM(a,ratom){ 
      catom = &(*a);   
      bond = mol->GetBond(ratom,catom); 
      if(atom_shares[catom] && !atoms_in_lp[catom] && !ignore_bond[bond]){
        stack.push(catom);
        break;
      }
    }
  }

  for(unsigned int i=0;i<path_size;i++){
    if(!locant_path[i]){
      free(locant_path);
      Fatal("no continous locant path was possible - currently unsupported\n");
    }
  }

  return locant_path; 
}

/* uses a flood fill style solution (likely NP-HARD), with some restrictions to 
find a multicyclic path thats stable with disjoined pericyclic points */
OBAtom **NPLocantPath(      OBMol *mol, unsigned int path_size,
                            std::set<OBAtom*>               &ring_atoms,
                            std::set<OBAtom*>               &bridge_atoms,
                            std::map<OBAtom*,unsigned int>  &atom_shares,
                            std::set<OBRing*>               &local_SSSR)
{

  // create the path
  OBAtom **locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
  OBAtom **best_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
  for(unsigned int i=0;i<path_size;i++){
    locant_path[i] = 0;
    best_path[i] = 0; 
  }

  // if there are bridges, set up a bridge map
  std::map<OBAtom*,bool> is_bridging; 
  for(std::set<OBAtom*>::iterator biter = bridge_atoms.begin(); biter != bridge_atoms.end(); biter++)
    is_bridging[*biter] = true; 

  // parameters needed to seperate out the best locant path
  unsigned int           lowest_sum       = UINT32_MAX;
  unsigned char          lowest_non_multi = int_to_locant(path_size);
  unsigned char          lowest_multi     = int_to_locant(path_size); 

  // multi atoms are the starting seeds, must check them all unfortuanately 
  for(std::set<OBAtom*>::iterator aiter = ring_atoms.begin(); aiter != ring_atoms.end(); aiter++){
    OBAtom *rseed = (*aiter);
    if(atom_shares[rseed] == 3 || is_bridging[rseed]){

      OBAtom*                catom  = 0;
      std::map<OBAtom*,bool> current; 
      std::deque<std::pair<OBAtom*,OBAtom*>> path; 
      path.push_back({rseed,0}); 
      while(!path.empty()){

        OBAtom* ratom = path.back().first;
        OBAtom* next = path.back().second;  

        if(!current[ratom])
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
            else if(!current[catom] && skipped){
              path.push_back({catom,0});
              pushed = true; 
              break;
            }
          }
        }

        if(!pushed){
          if(path.size() == path_size){

            for(unsigned int i=0;i<path_size;i++)
              locant_path[i] = path[i].first;
            
            // calculate the fusion sum here, expensive but necessary
            unsigned int fusion_sum = 0;
            for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end();riter++){
              OBRing *obring = (*riter); 

              unsigned int min_loc = path_size;
              for(unsigned int i=0;i<obring->Size();i++){
                unsigned int pos = position_in_path( mol->GetAtom(obring->_path[i]) ,locant_path,path_size);
                if(pos < min_loc)
                  min_loc = pos;
              }

              fusion_sum+=min_loc;
            }

            // trying to tease out 30e without a notation write and compare
            unsigned char earliest_non_multi = 0;
            unsigned char highest_multi = 0;
            for(int i=0;i< (int)path_size;i++){
              OBAtom *src = locant_path[i];
              if(atom_shares[src] == 3)
                highest_multi = int_to_locant(i+1);
              for(int j=0;j<i;j++){
                OBAtom *trg = locant_path[j]; 
                if(mol->GetBond(src,trg)){
                  if( atom_shares[locant_path[i]]==2 &&
                      atom_shares[locant_path[j]] == 2 &&
                      !earliest_non_multi)
                    earliest_non_multi = int_to_locant(i+1);
                }
                  
                               
              }
            }

            if(fusion_sum < lowest_sum){ // rule 30d.
              lowest_sum = fusion_sum;
              lowest_non_multi = earliest_non_multi; 
              lowest_multi = highest_multi;
              copy_locant_path(best_path,locant_path,path_size);
              if(opt_debug)
                fprintf(stderr,"  set on fs:  fusion sum for path: %-2d, lowest non-multi: %c, highest_multi: %c\n",lowest_sum,lowest_non_multi,lowest_multi);
            }
            else if(fusion_sum == lowest_sum && earliest_non_multi && earliest_non_multi < lowest_non_multi){ // rule 30e.
              lowest_non_multi = earliest_non_multi; 
              lowest_multi = highest_multi;
              copy_locant_path(best_path,locant_path,path_size);
              if(opt_debug)
                fprintf(stderr,"  set on enm: fusion sum for path: %-2d, lowest non-multi: %c, highest_multi: %c\n",lowest_sum,lowest_non_multi,lowest_multi);
            }
            else if(fusion_sum == lowest_sum && earliest_non_multi == lowest_non_multi && highest_multi < lowest_multi){ // not officially a rule but a good filter
              lowest_multi = highest_multi;
              copy_locant_path(best_path,locant_path,path_size);
              if(opt_debug)
                fprintf(stderr,"  set on hm:  fusion sum for path: %-2d, lowest non-multi: %c, highest_multi: %c\n",lowest_sum,lowest_non_multi,lowest_multi);
            }
            else {
              // if(opt_debug)
              //   fprintf(stderr,"  other:      fusion sum for path: %-2d, lowest non-multi: %c, highest_multi: %c\n",fusion_sum,earliest_non_multi,highest_multi);
            }


          }
          
          OBAtom *tmp = path.back().first; 
          path.pop_back();

          if(!path.empty()){
            path.back().second = tmp;
            current[tmp] = false; 
          }
          
        }
            
      }
    }
  }

  // do a count check here or|else return null for unsuccessful locant path, - future 
  free(locant_path);
  locant_path=0;

  for(unsigned int i=0;i<path_size;i++){
    if(!best_path[i]){
      free(best_path);
      Fatal("no continous locant path was possible - currently unsupported\n");
    }
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



struct RingWrapper{
  unsigned char loc_a;
  unsigned char loc_b; 
  OBRing *ring; 
};

void write_wrapper(RingWrapper *wrapper, std::string &buffer){
  
  if(wrapper->loc_a != 'A'){
    buffer += ' ';
    buffer += wrapper->loc_a;
  }

  if(wrapper->ring->Size() > 9){
    buffer+='-';
    buffer+= std::to_string(wrapper->ring->Size());
    buffer+='-';
  }
  else
    buffer+= std::to_string(wrapper->ring->Size());
}

void sort_wrapper(RingWrapper **wrapper_stack,unsigned int stack_size){
	for (unsigned int j=1;j<stack_size;j++){
		unsigned char key = wrapper_stack[j]->loc_b;
    RingWrapper *ptr = wrapper_stack[j];
		int i = j-1;
		while(i>=0 && wrapper_stack[i]->loc_b > key){
			wrapper_stack[i+1] = wrapper_stack[i];
			i--;
		}
		wrapper_stack[i+1] = ptr;
	}
}

bool valid_pair(unsigned char loc_a,unsigned char loc_b, RingWrapper **ring_stack,unsigned int stack_size){
  for(unsigned int i=0;i<stack_size;i++){
    if(ring_stack[i]->loc_a == loc_a && ring_stack[i]->loc_b == loc_b)
      return true;
  }
  return false;
}

unsigned char add_to_path(OBMol *mol,OBRing *ring, OBAtom **locant_path, unsigned int path_size, std::map<unsigned char,bool> &in_path){
  unsigned char highest_added = 'A';
  for(unsigned int i=0;i<ring->Size();i++){
    unsigned char loc = int_to_locant(position_in_path(mol->GetAtom(ring->_path[i]),locant_path,path_size)+1); 
    in_path[loc] = true;

    if(loc > highest_added)
      highest_added = loc; 
  }
  return highest_added;
}



bool ReadLocantPath(  OBMol *mol, OBAtom **locant_path, unsigned int path_size,
                      std::set<OBAtom*>               &bridge_atoms,
                      std::map<OBAtom*,unsigned int>  &atom_shares, 
                      std::set<OBRing*>               &local_SSSR,
                      std::vector<OBRing*>            &ring_order,
                      std::string &buffer)
{  
  
  unsigned int stack_size = 0;
  RingWrapper **ring_stack = (RingWrapper**)malloc(sizeof(RingWrapper) * path_size); // path size will be a hard limit
  for(unsigned int i=0;i<path_size;i++)
    ring_stack[i] = 0;

  std::map<OBRing*,bool> rings_checked;
  std::map<unsigned char,unsigned int> char_in_stack;

  for(int i=0;i<(int)path_size;i++){
    OBAtom *src = locant_path[i]; 
    if(atom_shares[src] > 1){
      for(int j=0;j<i;j++){
        OBAtom *trg = locant_path[j];

        if(atom_shares[trg] > 1 && mol->GetBond(src,trg) && ( (j == i-1 && atom_shares[trg] == 3) || j < i-1)){ 
          RingWrapper *wrapped = (RingWrapper*)malloc(sizeof(RingWrapper));
          wrapped->loc_a = int_to_locant(j+1);
          wrapped->loc_b = int_to_locant(i+1);
          wrapped->ring = 0; 
          ring_stack[stack_size++] = wrapped; 
        }

      }
    }
  }

  for(unsigned int i=0;i<stack_size;i++){
    // for multicyclic which ring contains the position after loc_b, should only be one
    // for non multicyclic, which ring contains the position after loc_a, again only one
    OBAtom *src = locant_path[ ring_stack[i]->loc_a-'A'];
    OBAtom *trg = locant_path[ ring_stack[i]->loc_b-'A'];
    OBBond *lbond = mol->GetBond(src,trg); 

    OBAtom *find = 0;
    if(atom_shares[src]==3)
      find = locant_path[ ring_stack[i]->loc_b-'A' +1];
    else
      find = locant_path[ ring_stack[i]->loc_a-'A' +1];

    OBRing *obring = 0;
    for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end();riter++){
      if(!rings_checked[(*riter)] && (*riter)->IsMember(lbond) && (*riter)->IsMember(find)){
        obring = (*riter);
        break;
      }
    } 

    if(!obring){
      free(ring_stack[i]);
      ring_stack[i] = 0;
    }
    else{

      // what if we also change loc_b here to the highest in ring

      if(opt_debug)
        fprintf(stderr,"  assigning (%c -> %c) ring with %c\n",ring_stack[i]->loc_a,ring_stack[i]->loc_b,int_to_locant(1+position_in_path(find,locant_path,path_size)));

      char_in_stack[ring_stack[i]->loc_a]++;

      // assign the new loc_a to lowest locant seen in the assigned ring
      unsigned int lowest_loc = path_size; 
      unsigned int highest_loc = 0;
      for(unsigned int j=0;j<obring->Size();j++){
        unsigned int npos = position_in_path(mol->GetAtom(obring->_path[j]),locant_path,path_size);
        if(npos < lowest_loc)
          lowest_loc = npos;

        if(npos > highest_loc)
          highest_loc = npos;  
      }

      ring_stack[i]->loc_b = int_to_locant(highest_loc+1); // BETA

      unsigned char pa = int_to_locant(lowest_loc+1);
      if(pa != ring_stack[i]->loc_a){

        if(opt_debug)
          fprintf(stderr,"    shifting to %c\n",pa);

        char_in_stack[pa]++;
        char_in_stack[ring_stack[i]->loc_a]--;

        ring_stack[i]->loc_a = pa; 
      }

      ring_stack[i]->ring = obring;
      rings_checked[obring] = true;
    }
  }

  // clean up stage 
  unsigned int idx = 0;
  for(unsigned int i=0;i<stack_size;i++){
    RingWrapper *r = ring_stack[i];
    ring_stack[i] = 0;
    if(r)
      ring_stack[idx++] = r;
  }
  stack_size = idx;

  // check if we've left any behind, always seems to happen on final looping closure when on a flat angle to 'A'
  unsigned int missed = 0;
  for(std::set<OBRing*>::iterator riter = local_SSSR.begin(); riter != local_SSSR.end();riter++){
    OBRing *cring = *(riter);
    if(!rings_checked[cring]){
      unsigned int lowest_loc = path_size; 
      unsigned int highest_loc = 0;
      for(unsigned int j=0;j<cring->Size();j++){
       unsigned int npos = position_in_path(mol->GetAtom(cring->_path[j]),locant_path,path_size);
        if(npos < lowest_loc)
          lowest_loc = npos; 
        if(npos > highest_loc)
          highest_loc = npos;  
      }

      RingWrapper *wrapped = (RingWrapper*)malloc(sizeof(RingWrapper));
      wrapped->loc_a = int_to_locant(lowest_loc+1);
      wrapped->loc_b = int_to_locant(highest_loc+1);
      char_in_stack[wrapped->loc_a]++;

      wrapped->ring = cring; 
      ring_stack[stack_size++] = wrapped; 
      missed++;
    }
  }

  if(missed > 1)
    Fatal("more than one ring seems to have been missed in read path");


  if(opt_debug){
    fprintf(stderr,"  pre-read:\n");
    for(unsigned int i=0;i<stack_size;i++){
      RingWrapper *wrapper = ring_stack[i];
      fprintf(stderr,"    %c --> %c [%d] contains: [ ",wrapper->loc_a,wrapper->loc_b,char_in_stack[wrapper->loc_a]);
      print_ring_locants(mol,wrapper->ring,locant_path,path_size);
      fprintf(stderr,"]\n");
    } 
  }

  // path walk writing algorithm
  RingWrapper **local_stack = (RingWrapper**)malloc(sizeof(RingWrapper) * stack_size); // should be a hard limit, write all at once
  for(unsigned int i=0;i<stack_size;i++)
    local_stack[i] = 0;
  
  for(int i=0;i< (int)path_size;i++){
    unsigned char src_char = int_to_locant(i+1); 

    for(int j=0;j<i;j++){
      unsigned char trg_char = int_to_locant(j+1); 

      if(char_in_stack[trg_char] && valid_pair(trg_char,src_char,ring_stack,stack_size)){
        
        // we can sort here to guarentee that pseudo bridges are placed correctly
        unsigned int seen = 0;
        for(unsigned int k=0;k<stack_size;k++){
          RingWrapper *wrapper = ring_stack[k];
          if(wrapper->loc_a == trg_char && wrapper->loc_b <= src_char)
            local_stack[seen++] = wrapper; 
        }

        if(char_in_stack[trg_char] == seen){
          sort_wrapper(local_stack,seen);
          for(unsigned int p=0;p<seen;p++){
            write_wrapper(local_stack[p],buffer);
            // we can check for pseudo bridges here - are there any non consecutive bonds in the sorted path?
            BondedEndPoints(mol,local_stack[p]->ring,locant_path,path_size,bridge_atoms,atom_shares,buffer);
            ring_order.push_back(local_stack[p]->ring);
          }
            
        }
        // clear the local stack
        for(unsigned int i=0;i<stack_size;i++)
          local_stack[i] = 0;
      }
    }

  }

  for(unsigned int i=0;i<stack_size;i++)
    free(ring_stack[i]);

  free(ring_stack);
  free(local_stack);
  return true;  
}


/**********************************************************************
                          Reduction Functions
**********************************************************************/

std::string CondenseCarbonylChains(std::string &buffer){
  
  bool special = false;
  unsigned int counter = 0; 
  std::string condensed = {}; 
  for(unsigned int i=0;i<buffer.size();i++){
    
    if(buffer[i] == '-'){
      if(special)
        special = false;
      else
        special = true;
    }
    else if (buffer[i] == ' ')
      special = false;
  
    if(buffer[i] == '1' && !special)
      counter++; 
    else{
      if(counter){
        condensed += std::to_string(counter);
        counter = 0; 
      }
      condensed += buffer[i];
    }
  }
  if(counter)
    condensed += std::to_string(counter);
  return condensed;
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
  std::map<OBRing*, bool> rings_seen; 

  std::map<OBAtom*,unsigned char> given_char; // post lookup
  std::map<OBAtom*,unsigned int>  evaluated_branches; // tracking for branch pop

  std::vector<OBAtom**> stored_paths; // mem pool

  // recursion tracking
  unsigned int cycle_count; 
  unsigned int last_cycle_seen;
  
  BabelGraph(){
    cycle_count = 0; 
    last_cycle_seen = 0; 
  };
  ~BabelGraph(){
    for(unsigned int i=0;i<stored_paths.size();i++){
      free(stored_paths[i]);
      stored_paths[i] = 0; 
    }
  };

  bool WriteBabelAtom(OBAtom* atom, std::string &buffer){

    if(!atom){
      fprintf(stderr,"Error: nullptr OpenBabel Atom*\n");
      return false; 
    }

    unsigned int neighbours = atom->GetExplicitDegree(); 
    unsigned int orders = atom->GetExplicitValence(); 

    switch(atom->GetAtomicNum()){
      case 1:
        buffer += 'H';
        break; 

      case 5:
        buffer += 'B';
        break;

      case 6:
        if(neighbours <= 2)
          buffer += '1';
        else if(neighbours == 3)
          buffer += 'Y';
        else if(neighbours == 4)
          buffer += 'X';
        break;
      
      case 7:
        switch(orders){
          case 0:
          case 1:
            buffer += 'Z'; 
            break;
          case 2:
            buffer += 'M';
            break;
          case 3:
            buffer += 'N';
            break;
          case 4:
            buffer += 'K';
            break;
          default: 
            fprintf(stderr,"Error: Unrecognised nitrogen bonding enviroment, orders - %d\n",orders);
            return 0;

        }   
        break;
      
      case 8:
        if(atom->GetExplicitValence() < 2 && atom->GetFormalCharge() != -1)
          buffer += 'Q';
        else
          buffer += 'O';
        break;
      
      case 9:
        if(atom->GetExplicitValence() > 1)
          buffer += "-F-";
        else
          buffer += 'F';
        break;

      case 15:
        buffer += 'P';
        if(neighbours == 1 && orders == 1)
          buffer+='H';
        break;

      case 16:
        buffer += 'S';
        if(neighbours == 1 && orders == 1)
          buffer+='H';
        break;

      case 17:
        if(atom->GetExplicitValence() > 1)
          buffer += "-G-";
        else
          buffer += 'G';
        break;

      case 35:
        if(atom->GetExplicitValence() > 1)
          buffer += "-E-";
        else
          buffer += 'E';
        break;

      case 53:
        if(atom->GetExplicitValence() > 1)
          buffer += "-I-";
        else
          buffer += 'I';
        break;



  // all special elemental cases

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
      

      default:
        fprintf(stderr,"Error: unhandled element for WLNSymbol formation\n");
        return false;
    }

    given_char[atom] = buffer.back(); // assign the last char, - will mean hypervalent

    return true; 
  }

  unsigned int CountDioxo(OBAtom *atom){
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
          
        if(carbonyls == 2 || oxo_ions == 2 || (oxo_ions == 1 && carbonyls == 1)){
          Ws++;
          atoms_seen[seen[0]] = true;
          atoms_seen[seen[1]] = true;
          carbonyls = 0;
          oxo_ions = 0;
          seen.clear();
        }
      }  
    }
    return Ws;
  }

  bool CheckCarbonyl(OBAtom *atom){
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

  /* parse non-cyclic atoms DFS style - return last atom seen in chain */
  bool ParseNonCyclic(OBAtom* start_atom, unsigned int b_order,
                      OBMol *mol, std::string &buffer, 
                      unsigned int cycle_num, unsigned char locant){
    
    //##################################
    //      INDIRECT RECURSION TRACKING

    if(last_cycle_seen > cycle_num){
      for(unsigned int i=0;i<(last_cycle_seen-cycle_num);i++){
        buffer+='&';
        
        if(cycle_count)
          cycle_count+= -1; // once a ring is closed can you ever get back? - GOOD
      }
    }
    last_cycle_seen = cycle_num;

    if(locant){
      buffer+=' ';
      buffer+= locant;
    }

    for(unsigned int b=1;b<b_order;b++)
      buffer+='U';
    //##################################


    unsigned int Wgroups = 0;
  
    OBAtom* atom = start_atom;
    OBAtom* prev = 0; 
    OBBond *bond = 0; 

    std::stack<OBAtom*> atom_stack; 
    std::stack<OBAtom*> branch_stack; 
    atom_stack.push(atom); 

    while(!atom_stack.empty()){
      atom = atom_stack.top(); 
      atom_stack.pop();
      atoms_seen[atom] = true;

      if(prev){
        bond = mol->GetBond(prev,atom); 
        if(!bond){

          if(!branch_stack.empty() && prev != branch_stack.top()){
            prev = branch_stack.top();
            buffer += '&';
          }

          // need to have a method where the branches are tracked 
          while(!branch_stack.empty() && !mol->GetBond(prev,atom)){
            branch_stack.pop();

            switch(given_char[prev]){ // sorts out last closures
              
              case 'Y':
              case 'B':
              case 'N':
                if(evaluated_branches[prev] < 2 )
                  buffer += '&';
                break;

              case 'X':
              case 'K':
                if(evaluated_branches[prev] < 3 )
                  buffer += '&';
                break;

              case 'P':
                if(evaluated_branches[prev] < 4 )
                  buffer += '&';
                break;

              case 'S':
              case '-':
                if(evaluated_branches[prev] < 5 )
                  buffer += '&';
                break;

              default:
                fprintf(stderr,"Error: unrecognised branching assignment - %c\n",given_char[prev]);
                return 0;
            }
      
            prev = branch_stack.top();
          }

          bond = mol->GetBond(prev,atom); 
        }

        if(!bond)
          Fatal("failure to read branched bond segment");

        for(unsigned int i=1;i<bond->GetBondOrder();i++)
          buffer += 'U';

        evaluated_branches[prev]++;
      }

      if(atom->IsInRing()){
        buffer+= '-';
        buffer+= ' '; 

        cycle_count++;
        if(!RecursiveParse(atom,mol,true,buffer,cycle_count))
          Fatal("failed to make inline ring");

        if(!atom_stack.empty()){
          buffer+='&';
          if(cycle_count)
            cycle_count--;

          prev = branch_stack.top();
        }
        continue;
      }

      if(!WriteBabelAtom(atom,buffer))
        Fatal("failed to write atom");
    
      // last added char, not interested in the ring types of '-'
      switch(buffer.back()){
        
// oxygens
        case 'O':
        case 'V':
        case 'M':
        case 'W': // W is not actually seen as a prev,
          prev = atom; 
          break;

// carbons 
        // alkyl chain 
        case '1':
          prev = atom; 
          if(CheckCarbonyl(atom))
            buffer.back() = 'V';
          break;

        case 'Y':
        case 'X':
          prev = atom;
          Wgroups = CountDioxo(atom);
          for(unsigned int i=0;i<Wgroups;i++){
            buffer+='W';
            evaluated_branches[atom]+=3;
          }
            
          if(!Wgroups && CheckCarbonyl(atom))
            buffer.back() = 'V';
          else
            branch_stack.push(atom);
          break;

        case 'N':
        case 'K':
        case 'B':
        case 'S':
        case 'P':
        case '-':
          Wgroups = CountDioxo(atom);
          for(unsigned int i=0;i<Wgroups;i++){
            buffer+='W';
            evaluated_branches[atom]+=3;
          }
          
          prev = atom;
          branch_stack.push(atom);
          break;
          
      
// halogens
        case 'Q':
        case 'Z':
        case 'E':
        case 'F':
        case 'G':
        case 'I':
          if(atom->GetExplicitValence() == 0 && atom->GetFormalCharge() == 0)
            buffer += 'H';

          if(!branch_stack.empty())
            prev = branch_stack.top();
          break;

        case 'H':
          prev = atom;
          break;


        default:
          fprintf(stderr,"Error: unhandled char %c\n",buffer.back()); 
          return 0; 
      }

      FOR_NBORS_OF_ATOM(a,atom){
        if(!atoms_seen[&(*a)])
          atom_stack.push(&(*a));
      }

    }

    return atom; 
  }

  /* parses the local ring system, return the size for creating the locant path with 
  non bonds to avoid */
  unsigned int ConstructLocalSSSR(  OBMol *mol, OBAtom *ring_root,
                                    std::set<OBAtom*>  &ring_atoms,
                                    std::set<OBBond*>  &ring_bonds,
                                    std::set<OBAtom*>  &bridging_atoms,
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

          // intersection 1 is a spiro ring, ignore, 
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
    if(!tmp_bridging_atoms.empty()){
      for(std::set<OBAtom*>::iterator brd_iter=tmp_bridging_atoms.begin(); brd_iter != tmp_bridging_atoms.end();brd_iter++){

        unsigned int inter_ring_bonds = 0;
        for(std::set<OBAtom*>::iterator aiter= ring_atoms.begin(); aiter != ring_atoms.end();aiter++){
          if(mol->GetBond(*brd_iter,*aiter))
            inter_ring_bonds++; 
        }
        if(inter_ring_bonds == 2)
          bridging_atoms.insert(*brd_iter); 
      }
    }

    if(opt_debug){
      fprintf(stderr,"  ring atoms: %lu\n",ring_atoms.size());
      fprintf(stderr,"  ring bonds: %lu\n",ring_bonds.size());
      if(!bridging_atoms.empty())
        fprintf(stderr,"  bridging atoms: %lu\n",bridging_atoms.size());
      
    }
      
    return ring_atoms.size(); 
  }


  /* create the heteroatoms and locant path unsaturations where neccesary */
  bool ReadLocantAtomsBonds(  OBAtom** locant_path,unsigned int path_size,
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

      locant = int_to_locant(i+1);
      unsigned int Wgroups = CountDioxo(locant_path[i]);
      bool carbonyl = CheckCarbonyl(locant_path[i]);

      if(carbonyl || Wgroups || locant_path[i]->GetAtomicNum() != 6){
        if(locant != last_locant){
          buffer += ' ';
          write_locant(locant,buffer);
          last_locant = locant;
        } 
        if(Wgroups){
          if(!WriteBabelAtom(locant_path[i],buffer))
            return false; 

          for(unsigned int w=0;w<Wgroups;w++)
            buffer+='W';

          last_locant++;
        }
        else if(carbonyl){
          buffer += 'V';
          last_locant++;
        }
        else{
          WriteBabelAtom(locant_path[i],buffer);
          last_locant++;
        }
      }

    
#define WIP 1
#if WIP
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
#endif
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


  /* constructs and parses a cyclic structre, locant path is returned with its path_size */
  std::pair<OBAtom **,unsigned int> ParseCyclic(OBAtom *ring_root,OBMol *mol, bool inline_ring,std::string &buffer){
    if(opt_debug)
      fprintf(stderr,"Reading Cyclic\n");

    OBAtom **                       locant_path = 0; 
    std::set<OBRing*>               local_SSSR;
    std::set<OBAtom*>               ring_atoms;
    std::set<OBBond*>               ring_bonds;
    std::set<OBAtom*>               bridge_atoms;
    
    std::vector<OBRing*>            ring_order; 
    
    std::map<OBAtom*,unsigned int>  atom_shares;
  
    unsigned int path_size   =  ConstructLocalSSSR(mol,ring_root,ring_atoms,ring_bonds,bridge_atoms,atom_shares,local_SSSR); 
    if(!path_size)
      Fatal("failed to write ring");


    bool multi = ClassifyRing(ring_atoms,atom_shares); 
    
    if(multi || !bridge_atoms.empty())
      locant_path = NPLocantPath(mol,path_size,ring_atoms,bridge_atoms,atom_shares,local_SSSR);
    else
      locant_path = PLocantPath(mol,path_size,ring_atoms,ring_bonds,atom_shares,local_SSSR);
      
    for(unsigned int i=0;i<path_size;i++){
      if(inline_ring && locant_path[i] == ring_root)
        buffer += int_to_locant(i+1);
    }

    if(IsHeteroRing(locant_path,path_size))
      buffer += 'T';
    else
      buffer += 'L';

    ReadLocantPath( mol,locant_path,path_size,
                    bridge_atoms,atom_shares,local_SSSR,ring_order,
                    buffer); 
    
    if(!bridge_atoms.empty()){
      for(std::set<OBAtom*>::iterator biter = bridge_atoms.begin(); biter != bridge_atoms.end(); biter++){
        buffer+= ' ';
        unsigned char bloc = int_to_locant(position_in_path(*biter,locant_path,path_size)+1);
        write_locant(bloc,buffer);
      }
    }

    unsigned int size_check = 0;
    if(multi){
      ReadMultiCyclicPoints(locant_path,path_size,atom_shares,buffer);
      buffer += ' ';
      write_locant(int_to_locant(path_size),buffer); // need to make the relative size
      size_check = buffer.size();
    }

    ReadLocantAtomsBonds(locant_path,path_size,ring_order,ring_bonds,buffer);

    if(buffer.size() == size_check)
      buffer += ' '; // has to be a choice whether to add the space here
    
    bool space_added = false;
    bool dash_added = false;
    for(unsigned int i=0;i<ring_order.size();i++){
      // check aromaticity
      if(ring_order[i]->IsAromatic()){
        if(!dash_added && buffer.back() == '&')
          buffer+='-';
        else if(!space_added && (buffer.back() >= 'A' && buffer.back() <= 'Z'))
          buffer+=' ';
        
        buffer+='&';
      }
      else
        buffer += 'T';

      space_added = true;
      dash_added = true;
    }

    buffer += 'J';
    return {locant_path,path_size};
  }
    

  bool RecursiveParse(OBAtom *atom, OBMol *mol, bool inline_ring,std::string &buffer, unsigned int cycle_num){
    // assumes atom is a ring atom 
    last_cycle_seen = cycle_num;
    std::pair<OBAtom**,unsigned int> path_pair;  
    path_pair = ParseCyclic(atom,mol,inline_ring,buffer);

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
          if(!ParseNonCyclic( latom,lbond->GetBondOrder(),
                              mol,buffer,
                              cycle_num,int_to_locant(i+1))){
            fprintf(stderr,"Error: failed on non-cyclic parse\n");
            return false;
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
 
  BabelGraph obabel; 
  unsigned int cyclic = 0;
  bool started = false; 
  FOR_RINGS_OF_MOL(r,mol)
    cyclic++;

  if(opt_debug)
    WriteBabelDotGraph(mol);

  if(!cyclic){
    
    FOR_ATOMS_OF_MOL(a,mol){
      if(!obabel.atoms_seen[&(*a)]){
        if(started)
          buffer += " &"; // ionic species
        
        if(!obabel.ParseNonCyclic(&(*a),0,mol,buffer,0,0))
          Fatal("failed on recursive branch parse");

        started = true; 
      }
    }
  }
  else{
    FOR_RINGS_OF_MOL(r,mol){
    // start recursion from first cycle atom
      if(!obabel.rings_seen[&(*r)]){
        if(started)
          buffer += " &"; // ionic species

        if(!obabel.RecursiveParse(mol->GetAtom( (&(*r))->_path[0]),mol,false,buffer,0))
          Fatal("failed on recursive ring parse");

        started = true;
      }
    }
    
    // handles additional ionic atoms here
    obabel.cycle_count = 0;
    obabel.last_cycle_seen = 0;
    FOR_ATOMS_OF_MOL(a,mol){
      if(!obabel.atoms_seen[&(*a)]){
        buffer += " &"; // ionic species
        if(!obabel.ParseNonCyclic(&(*a),0,mol,buffer,0,0))
          Fatal("failed on recursive branch parse");
      }
    }
    
  }

  buffer = CondenseCarbonylChains(buffer); 
  return true; 
}


static void DisplayUsage()
{
  fprintf(stderr, "writewln <options> -i<format> -s <input (escaped)>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -d                    print debug messages to stderr\n");
  fprintf(stderr, "  -h                    show the help for executable usage\n");
  fprintf(stderr, "  -i                    choose input format (-ismi, -iinchi, -ican)\n");
  fprintf(stderr, "  -w                    dump wln trees & babel graphs to dot files in [build]\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser writes to wiswesser\n"
                  " line notation (wln) from smiles/inchi, the parser is built on OpenBabels\n"
                  " toolkit and will return the minimal WLN string\n");
  DisplayUsage();
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i;

  cli_inp = (const char *)0;
  format = (const char *)0;

  if (argc < 2)
    DisplayUsage();

  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){

        case 'd':
          opt_debug = true;
          break;

        case 'h':
          DisplayHelp();


        case 'i':
          if (!strcmp(ptr, "-ismi"))
          {
            format = "smi";
            break;
          }
          else if (!strcmp(ptr, "-iinchi"))
          {
            format = "inchi";
            break;
          }
          else if (!strcmp(ptr, "-ican"))
          {
            format = "can";
            break;
          }
          else{
            fprintf(stderr,"Error: unrecognised format, choose between ['smi','inchi','can']\n");
            DisplayUsage();
          }

        case 's':
          if(i+1 >= argc){
            fprintf(stderr,"Error: must add string after -s\n");
            DisplayUsage();
          }
          else{
            cli_inp = argv[i+1];
            i++;
          }
          break;

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
  }

  if(!format){
    fprintf(stderr,"Error: no input format selected\n");
    DisplayUsage();
  }

  if(!cli_inp){
    fprintf(stderr,"Error: no input string entered\n");
    DisplayUsage();
  }

  return;
}

int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);
  
  std::string res;
  OBMol mol;
  OBConversion conv;

  conv.SetInFormat(format);
  res = conv.ReadString(&mol,cli_inp);

  std::string buffer;
  buffer.reserve(1000);
  if(!WriteWLN(buffer,&mol))
    return 1;
  
  std::cout << buffer << std::endl;

  return 0;
}


