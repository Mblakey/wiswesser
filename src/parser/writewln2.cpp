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
static bool opt_wln2dot = false;
static bool opt_debug = false;


unsigned char static int_to_locant(unsigned int i){
  return i + 64;
}

unsigned int static locant_to_int(unsigned char loc){
  return loc - 64;
}

unsigned char static create_relative_position(unsigned char parent){
  // A = 129
  unsigned int relative = 128 + locant_to_int(parent);
  if(relative > 252){
    fprintf(stderr,"Error: relative position is exceeding 252 allowed space - is this is suitable molecule for WLN notation?\n");
    return '\0';
  }
  else
    return relative;
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


/**********************************************************************
                          Ring Construction Functions
**********************************************************************/

/* constructs the local rings system, return highest ring share value */
unsigned int ConstructLocalSSSR(  OBAtom *ring_root, OBMol *mol, 
                                  std::set<OBAtom*> &ring_atoms,
                                  std::map<OBAtom*,unsigned int> &ring_shares,
                                  std::set<OBRing*> &local_SSSR){

  if(!ring_root){
    fprintf(stderr,"Error: ring root is nullptr\n");
    return 0; 
  }

  OBRing *seed_ring = 0;
  OBRing *obring = 0; 

  unsigned int fuses = 0;
  unsigned int multicyclic = 0;
  unsigned int branching = 0;   

  // get the seed ring and add path to ring_atoms
  FOR_RINGS_OF_MOL(r,mol){
    obring = &(*r);
    if(obring->IsMember(ring_root)){
      seed_ring = obring; // so we do not consider again
      local_SSSR.insert(obring);
      for(unsigned int i=0;i<obring->Size();i++){
        OBAtom *ratom = mol->GetAtom(obring->_path[i]);
        ring_atoms.insert(ratom);
        ring_shares[ratom]++;
      }
      break;
    }
  }
  

  FOR_RINGS_OF_MOL(r,mol){
    obring = &(*r);
    if(obring != seed_ring){
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
        local_SSSR.insert(obring);
        for(unsigned int i=0;i<obring->Size();i++){
          OBAtom *ratom = mol->GetAtom(obring->_path[i]);
          ring_atoms.insert(ratom);
          ring_shares[ratom]++;
        }
      }
    }
  }

  for(std::set<OBAtom*>::iterator iter=ring_atoms.begin(); iter != ring_atoms.end(); iter++){
    OBAtom *atom = *iter; 
    if(ring_shares[atom] == 2)
      fuses++; 
    else if(ring_shares[atom] == 3)
      multicyclic++;
    else if(ring_shares[atom] == 4)
      branching++;
  }


  if(opt_debug){
    fprintf(stderr,"  SSSR for system:    ");
    for(std::set<OBRing*>::iterator set_iter = local_SSSR.begin();set_iter != local_SSSR.end();set_iter++)
      fprintf(stderr,"%ld(%c) ",(*set_iter)->Size(), (*set_iter)->IsAromatic()?'a':'s');
    fprintf(stderr,"\n");

    fprintf(stderr,"  ring size:          %d\n",(unsigned int)ring_atoms.size());
    fprintf(stderr,"  fuse points:        %d\n",fuses);
    fprintf(stderr,"  multicyclic points: %d\n",multicyclic);
    fprintf(stderr,"  branching points:   %d\n",branching);
  }

  if(branching){
    fprintf(stderr,"NON-SUPPORTED: branching cyclics\n");
    return 0;
  }
  else if(multicyclic)
    return 3; 
  else 
    return 2; 
}


// get all potential seeds for locant path start
void GetSeedAtoms(  std::set<OBAtom*> &ring_atoms,
                    std::map<OBAtom*,unsigned int> &ring_shares,
                    std::vector<OBAtom*> &seed_atoms,
                    unsigned int target_shares)
{ 
  for(std::set<OBAtom*>::iterator iter = ring_atoms.begin(); iter != ring_atoms.end();iter++){
    if(ring_shares[(*iter)] == target_shares)
      seed_atoms.push_back((*iter));
  }
}


/* construct locant paths without hamiltonians - return the new locant pos */
unsigned int ShiftandAddLocantPath( OBMol *mol, OBAtom **locant_path,
                            unsigned int locant_pos,unsigned int path_size,
                            unsigned int hp_pos, OBRing *obring,
                            std::map<OBAtom*,bool> &atoms_seen,
                            std::vector<std::pair<OBAtom*,OBAtom*>> &nt_pairs,
                            std::vector<unsigned int> &nt_sizes)
                            
{
  
  bool seen = false;
  OBAtom *ratom = 0; 
  OBAtom *insert_start  =  locant_path[hp_pos];
  OBAtom *insert_end    =  locant_path[hp_pos+1]; 
  std::deque<int> path;

  for(unsigned int i=0;i<obring->Size();i++){
    path.push_back(obring->_path[i]);
    if(insert_end->GetIdx() == obring->_path[i])
      seen = true;
  }

  if(!seen){
    insert_start = locant_path[locant_pos-1];
    insert_end = locant_path[0];
  }
    

  while(path[0] != insert_start->GetIdx()){
    unsigned int tmp = path[0];
    path.pop_front();
    path.push_back(tmp);
  }

    
  // standard clockwise and anti clockwise additions to the path
  if(seen){
  
    // must be anti-clockwise - shift and reverse
    if(insert_start && path[1] == insert_start->GetIdx()){
      unsigned int tmp = path[0];
      path.pop_front();
      path.push_back(tmp);
      std::reverse(path.begin(),path.end());
    }

    if(opt_debug)
      fprintf(stderr,"  non-trivial bonds:  %-2d <--> %-2d from size: %ld\n",locant_path[hp_pos]->GetIdx(),locant_path[hp_pos+1]->GetIdx(),obring->Size());

    // add nt pair + size
    nt_pairs.push_back({locant_path[hp_pos],locant_path[hp_pos+1]});
    nt_sizes.push_back(obring->Size()); 

    // spit the locant path between hp_pos and hp_pos + 1, add elements
    unsigned int j=0;
    for(unsigned int i=0;i<path.size();i++){
      ratom = mol->GetAtom(path[i]);
      if(!atoms_seen[ratom]){
        // shift
        for(int k=path_size-1;k>hp_pos+j;k--) // potential off by 1 here. 
          locant_path[k]= locant_path[k-1];
          
        locant_path[hp_pos+1+j] = ratom;
        atoms_seen[ratom] = true;
        j++;
        locant_pos++;
      }
    }
  }

  // must be a ring wrap on the locant path, can come in clockwise or anticlockwise
  else{

    // this is the reverse for the ring wrap
    if(path[1] == insert_end->GetIdx()){ // end due to shift swap
      unsigned int tmp = path[0];
      path.pop_front();
      path.push_back(tmp);
      std::reverse(path.begin(),path.end());
    }

    // just add to the back, no shift required
    for(unsigned int i=0;i<path.size();i++){
      ratom = mol->GetAtom(path[i]);
      if(!atoms_seen[ratom]){
        locant_path[locant_pos++] = ratom;
        atoms_seen[ratom] = true;
      }
    }

    if(opt_debug)
      fprintf(stderr,"  non-trivial ring wrap:  %-2d <--> %-2d from size: %ld\n",locant_path[0]->GetIdx(),locant_path[locant_pos-1]->GetIdx(),obring->Size());


    // ending wrap condition 
    nt_pairs.push_back({locant_path[0],locant_path[locant_pos-1]});
    nt_sizes.push_back(obring->Size()); 

  } 

  return locant_pos;
}

/* works on priority, and creates locant path via array shifting, the locant path
currently working for MONO,POLY, unsupported: PERI,BRANCH */
OBAtom ** CreateLocantPath(   OBMol *mol, std::set<OBRing*> &local_SSSR, 
                              std::map<OBAtom*,unsigned int> &ring_shares,
                              std::vector<OBRing*> &ring_order,
                              std::vector<std::pair<OBAtom*,OBAtom*>> &nt_pairs,
                              std::vector<unsigned int> &nt_sizes,
                              unsigned int path_size,
                              OBAtom *seed_atom)
{

  OBAtom **locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
  for(unsigned int i=0;i<path_size;i++)
    locant_path[i] = 0;


  OBAtom *ratom = 0;
  OBRing *obring = 0;

  // get the ring with the seed atom present
  for(std::set<OBRing*>::iterator iter = local_SSSR.begin(); iter != local_SSSR.end(); iter++){
    for(unsigned int i=0;i<(*iter)->_path.size();i++){
      if(mol->GetAtom((*iter)->_path[i]) == seed_atom){
        obring = (*iter);
        break;
      }
    }
  }

  if(!obring){
    fprintf(stderr,"Error: seed atom could not be found in local SSSR\n");
    return 0; 
  }
  else
    ring_order.push_back(obring);

  
  unsigned int locant_pos = 0;
  std::map<OBRing*,bool> rings_seen; 
  std::map<OBAtom*,bool> atoms_seen; 

  // add into the array directly and shift so seed is guareented in position 0
  for(unsigned int i=0;i<obring->_path.size();i++){
    ratom = mol->GetAtom(obring->_path[i]);
    locant_path[locant_pos++] = ratom;
    atoms_seen[ratom] = true;
  }

  while(locant_path[0] != seed_atom){
    locant_path[locant_pos] = locant_path[0];
    for(unsigned int i=0;i<path_size-1;i++)
      locant_path[i] = locant_path[i+1];
  }

  if(opt_debug)
    fprintf(stderr,"  non-trivial bonds:  %-2d <--> %-2d from size: %ld\n",locant_path[0]->GetIdx(),locant_path[locant_pos-1]->GetIdx(),obring->Size());
  

  nt_pairs.push_back({locant_path[0],locant_path[locant_pos-1]});
  nt_sizes.push_back(obring->Size());


  // get next ring in locant order
  for(unsigned int rings_handled = 0; rings_handled < local_SSSR.size()-1;rings_handled++){
    rings_seen[obring] = true;
    unsigned int hp_pos = 0; 
    for(unsigned int i=0;i<locant_pos;i++){
      bool found = false;
      ratom = locant_path[i];
      if(ring_shares[ratom] > 1){
        for(std::set<OBRing*>::iterator iter = local_SSSR.begin(); iter != local_SSSR.end(); iter++){
          if(!rings_seen[(*iter)] && (*iter)->IsInRing(ratom->GetIdx())){
            hp_pos = i;
            obring = (*iter);
            ring_order.push_back(obring);
            found = true;
            break;
          }
        }
        if(found)
          break;
      }
    }

    locant_pos =  ShiftandAddLocantPath(  mol,locant_path,
                                          locant_pos,path_size,hp_pos,obring,
                                          atoms_seen,nt_pairs,nt_sizes);
    if(!locant_pos)
      return 0;
  }
  

  
  return locant_path;
}




bool IsHeteroRing(OBAtom **locant_array,unsigned int size){
  for(unsigned int i=0;i<size;i++){
    if(locant_array[i]->GetAtomicNum() != 6)
      return true;
  }
  return false; 
}


void UpdateReducedPath( OBAtom **reduced_path, OBAtom** locant_path, unsigned int size,
                        std::map<OBAtom*,unsigned int> &ring_shares){
  for(unsigned int i=0;i<size;i++){
    if(ring_shares[locant_path[i]] > 1)
      reduced_path[i] = locant_path[i];
    else
      reduced_path[i] = 0; 
  }
}


std::string ReadLocantPath(OBAtom **locant_path,unsigned int path_size,
                            std::map<OBAtom*,unsigned int> ring_shares, // copy unavoidable 
                            std::vector<std::pair<OBAtom*,OBAtom*>> &nt_pairs,
                            std::vector<unsigned int> &nt_sizes,
                            unsigned int expected_rings)
{
  
  std::string ring_str; 
  if(IsHeteroRing(locant_path,path_size))
    ring_str += 'T';
  else
    ring_str += 'L';


  // can we take an interrupted walk between the points, if so, write ring size 
  // and remove 

  // create a reduced array 
  OBAtom **reduced_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size);
  UpdateReducedPath(reduced_path,locant_path,path_size,ring_shares);

  if(opt_debug){
    fprintf(stderr,"  locant path:  ");
    print_locant_array(locant_path,path_size); 
    fprintf(stderr,"  reduced path: ");
    print_locant_array(reduced_path,path_size);
  }
  
  
  unsigned int safety = 0;
  while(!nt_pairs.empty() && safety < expected_rings){
  
    for(unsigned int i=0;i<nt_pairs.size();i++){
      OBAtom *first =  nt_pairs[i].first; 
      OBAtom *second = nt_pairs[i].second;

      // find the position of first in the array
      unsigned int pos = 0; 
      for(pos=0;pos < path_size;pos++){
        if(locant_path[pos] == first)
          break;
      }

      // can we go to the second without interuption
      bool popped = false;
      for(unsigned int j=pos+1;j < path_size;j++){
        if(reduced_path[j] && reduced_path[j] != second){
          // interuption - pair cannot be handled in this iteration
          // break out of search, search next pair
          break;
        }
        else if(reduced_path[j] && reduced_path[j] == second){
          // write the ring, and pop nt_pair and nt_ring at position 
          if(pos){
            ring_str+= ' ';
            ring_str+= int_to_locant(pos+1);
          }

          if(nt_sizes[i] > 9){
            ring_str += '-';
            ring_str += std::to_string(path_size); 
            ring_str += '-';
          } 
          else
            ring_str += std::to_string(nt_sizes[i]);

          nt_pairs.erase(nt_pairs.begin() + i);
          nt_sizes.erase(nt_sizes.begin() + i);

          // update the reduced locant path based on ring_shares
          ring_shares[first]--;
          ring_shares[second]--;
          UpdateReducedPath(reduced_path,locant_path,path_size,ring_shares);

          // requires reset to while loop 
          popped = true; 
          break;
        }
      }

      if(popped) // resets to while
        break;
    }
    
    safety++;
  }

  if( nt_pairs[0].first == locant_path[0] 
      && nt_pairs[0].second == locant_path[path_size-1])
  {
    // last implied ring wrap
    ring_str += std::to_string(nt_sizes[0]);
    nt_pairs.clear();
    nt_sizes.clear();
  }
  else{
    fprintf(stderr,"Error: safety caught on reduced locant loop\n");
    return {};
  }

  free(reduced_path);
  reduced_path=0;
  return ring_str;
}




/**********************************************************************
                          Canonicalisation Function
**********************************************************************/

unsigned int highest_unbroken_numerical_chain(std::string &str){
  unsigned int highest_chain = 0; 
  unsigned int current_chain = 0; 
  for(unsigned int i=0;i<str.size();i++){
    if(str[i] <= '9' && str[i] >= '0')
      current_chain++;
    else{
      if(current_chain > highest_chain)
        highest_chain = current_chain;
      
      current_chain = 0; 
    }
  }

  if(current_chain > highest_chain)
    highest_chain = current_chain;

  return highest_chain; 
}

unsigned char first_locant_seen(std::string &str){
  // ignore the L|T
  for(unsigned int i=1;i<str.size();i++){
    if(str[i] != ' ' && !(str[i] <= '9' && str[i] >= '0'))
      return str[i];
  }
  return 0;
}


/* returns the index of the minimal ring ring
- unbroken numerical chain count and lowest locant sum? */
unsigned int MinimalWLNRingNotation(std::vector<std::string> &ring_strings){

  unsigned int highest_chain = 0; 
  unsigned char lowest_loc  =  0;
  unsigned int return_index = 0;  
  for(unsigned int i=0;i<ring_strings.size();i++){
    unsigned int chain =  highest_unbroken_numerical_chain(ring_strings[i]); 
    unsigned char loc =  first_locant_seen(ring_strings[i]); 

    if(chain > highest_chain){
      highest_chain = chain;
      lowest_loc = loc; 
      return_index = i; 
    }
    else if (chain == highest_chain && lowest_loc > loc){
      lowest_loc = loc; 
      return_index = i;
    }
  }

  return return_index;
}


std::string CondenseCarbonylChains(std::string &buffer){
  
  unsigned int counter = 0; 
  std::string condensed = {}; 
  for(unsigned int i=0;i<buffer.size();i++){
    if(buffer[i] == '1')
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
        else if(neighbours > 2){
          if(orders == 3)
            buffer += 'Y';
          else
            buffer += 'X';
        }
        else
          buffer += 'C';
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
        break;

      case 16:
        buffer += 'S';
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
    //      INDIRECT RECURSION TRACK

    if(opt_debug)
      fprintf(stderr,"non-cyclic level: %d, last seen - %d\n", cycle_num, last_cycle_seen);

    if(last_cycle_seen > cycle_num){
      for(unsigned int i=0;i<(last_cycle_seen-cycle_num);i++){
        buffer+='&';
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

          if(prev != branch_stack.top()){
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
                if(evaluated_branches[prev] != 2 )
                  buffer += '&';
                break;

              case 'X':
              case 'K':
                if(evaluated_branches[prev] != 3 )
                  buffer += '&';
                break;

              case 'P':
              case 'S':
              case '-':
                if(evaluated_branches[prev] != 5 )
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

        if(!bond){
          fprintf(stderr,"Error: failure to read branched bond segment\n");
          return 0;
        }

        for(unsigned int i=1;i<bond->GetBondOrder();i++)
          buffer += 'U';

        evaluated_branches[prev]++;
      }

      if(atom->IsInRing()){
        buffer+= '-';
        buffer+= ' '; 

        cycle_count++;
        if(!RecursiveParse(atom,mol,true,buffer,cycle_count)){
          fprintf(stderr,"Error: failed to make inline ring\n");
          return false;
        }
        if(!atom_stack.empty()){
          buffer+='&';
          cycle_count--;
        }
        continue;
      }

      if(!WriteBabelAtom(atom,buffer))
        return 0;

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
          break;


        case 'Y':
        case 'X':
          prev = atom;
          Wgroups = CountDioxo(atom);
          for(unsigned int i=0;i<Wgroups;i++)
            buffer+='W';
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
          for(unsigned int i=0;i<Wgroups;i++)
            buffer+='W';

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

  /* create the hetero atoms where neccesary */
  bool ReadLocantBondsxAtoms( OBAtom** locant_path,unsigned int path_size,
                              std::vector<OBRing*> &ring_order,
                              std::string &buffer){

    unsigned char last_locant = 0; 
    unsigned char locant = 0;
    

    OBAtom *first = 0;
    OBAtom *second = 0; 
    for(unsigned int i=0;i<path_size;i++){
      unsigned int Wgroups = 0; 
      locant = int_to_locant(i+1); 
      if(locant_path[i]->GetAtomicNum() != 6){
        // handles 'A' starting and consecutive locants
        if(i > 0 && locant-1 != last_locant){
          buffer += ' ';
          buffer += locant;
        }

        Wgroups = CountDioxo(locant_path[i]);
        if(Wgroups){
          if(!WriteBabelAtom(locant_path[i],buffer))
            return false; 

          for(unsigned int w=0;w<Wgroups;w++)
            buffer+='W';
        }
        else if(CheckCarbonyl(locant_path[i]))
          buffer += 'V';
        else if(!WriteBabelAtom(locant_path[i],buffer))
          return false; 
        
        last_locant = locant; 
      }

      // handles sequential locant unsaturations, when not aromatic
      first = locant_path[i];
      if(i < path_size-1)
        second = locant_path[i+1];
      else
        second = locant_path[0];

      OBBond *locant_bond = first->GetBond(second);
      if(locant_bond->GetBondOrder() > 1){
        bool arom_skip = false; 
        for(unsigned int j=0;j<ring_order.size();j++){
          OBRing *ring = ring_order[j];
          if( (ring->IsMember(first) || ring->IsMember(second)) 
              && ring->IsAromatic()){
            arom_skip = true;
            break;
          }
        }
        if(!arom_skip){
          if(i > 0 && locant != last_locant){
            buffer += ' ';
            buffer += locant;
          }

          for(unsigned int b=1;b<locant_bond->GetBondOrder();b++)
            buffer += 'U';
        }
      }

      

      
    }
    
    return true;
  }


  void ReadAromaticity(std::vector<OBRing*> &ring_order,std::string &buffer){

    std::string append = ""; 
    for(unsigned int i=0;i<ring_order.size();i++){
      // mark all the rings as seen
      rings_seen[ring_order[i]] = true;

      // check aromaticity
      if(ring_order[i]->IsAromatic())
        append += '&';
      else
        append += 'T';
    }

    // check if all aromatic or not
    if(append.find_first_not_of(append[0]) == std::string::npos){

      if(append[0] != '&')
        buffer+='T';
    }
    else
      buffer+= append; 
  }
  
  /* constructs and parses a cyclic structre, locant path is returned with its path_size */
  std::pair<OBAtom **,unsigned int> ParseCyclic(OBAtom *ring_root,OBMol *mol, bool inline_ring,std::string &buffer){
    if(opt_debug)
      fprintf(stderr,"Reading Cyclic\n");
   
    OBAtom **                               locant_path = 0;
    std::set<OBRing*>                       local_SSSR; 
    std::set<OBAtom*>                       ring_atoms; 
    std::map<OBAtom*,unsigned int>          ring_shares; 
    std::vector<std::pair<OBAtom*,OBAtom*>> nt_pairs;
    std::vector<unsigned int>               nt_sizes;  
    
    // needed for minimising the locant path
    std::vector<std::string>                cyclic_strings;
    std::vector<OBAtom**>                   locant_paths;  
    std::vector<std::vector<OBRing*>>       cycle_orders;
    unsigned int                            minimal_index = 0; 

    unsigned int expected_rings = 0;
    unsigned int ring_type = ConstructLocalSSSR(ring_root,mol,ring_atoms,ring_shares,local_SSSR); 
    unsigned int path_size = ring_atoms.size(); 
    
    // monocyclic
    if(local_SSSR.size() == 1){
      locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size); 
      for(unsigned int i=0;i<path_size;i++)
        locant_path[i] = 0;
      stored_paths.push_back(locant_path);

      OBRing* monoring = *(local_SSSR.begin());
      rings_seen[monoring] = true;
      for(unsigned int i=0;i<monoring->Size();i++){
        locant_path[i] = mol->GetAtom(monoring->_path[i]);
        if(inline_ring && locant_path[i] == ring_root)
          buffer += int_to_locant(i+1);
      }
      cycle_orders.push_back({monoring});
        
      if(IsHeteroRing(locant_path,path_size))
        buffer += 'T';
      else
        buffer += 'L';

      if(path_size > 9){
        buffer += '-';
        buffer += std::to_string(path_size); 
        buffer += '-';
      }
      else
        buffer += std::to_string(path_size); 
      
      ReadLocantBondsxAtoms(locant_path,path_size,cycle_orders[0],buffer);
      if(!monoring->IsAromatic())
        buffer += 'T';
      buffer += 'J';
      return {locant_path,path_size};
    }
    
    else if(path_size && ring_type){
      // generate seeds
      std::vector<OBAtom*> seed_atoms; 
      expected_rings = local_SSSR.size(); 
      GetSeedAtoms(ring_atoms,ring_shares,seed_atoms,ring_type);
      if(seed_atoms.empty()){
        fprintf(stderr,"Error: no seeds found to build locant path\n");
        return {0,0};
      }

      for(unsigned int i=0;i<seed_atoms.size();i++){
        
        // create a new path per string
        std::vector<OBRing*> ring_order;
        locant_path = CreateLocantPath(   mol,
                                          local_SSSR,ring_shares,ring_order,
                                          nt_pairs,nt_sizes,
                                          path_size,
                                          seed_atoms[i]);
        if(!locant_path)
          return {0,0}; 
        else
          stored_paths.push_back(locant_path);

        std::string cyclic_str = ReadLocantPath(  locant_path,path_size,ring_shares,
                                                  nt_pairs,nt_sizes,
                                                  expected_rings);
        if(cyclic_str.empty())
          return {0,0};

        cyclic_strings.push_back(cyclic_str);
        locant_paths.push_back(locant_path);  
        cycle_orders.push_back(ring_order);

        if(opt_debug)
          std::cout << "  produced: " << cyclic_str << "\n\n";
      }


      // get the minal WLN cyclic system from the notations generated
      minimal_index = MinimalWLNRingNotation(cyclic_strings); 
      buffer+= cyclic_strings[minimal_index]; // add to the buffer

 
      // add any hetero atoms at locant positions
      ReadLocantBondsxAtoms(locant_paths[minimal_index],path_size,cycle_orders[minimal_index],buffer);
      
      // aromatics
      ReadAromaticity(cycle_orders[minimal_index],buffer);


      // close ring
      buffer += 'J';
      return {locant_paths[minimal_index],path_size};
    }
    
    return {0,0}; 
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
          return false;

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
          return false; 

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
          return false;
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

        case 'w':
          opt_wln2dot = true;
          break;

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


