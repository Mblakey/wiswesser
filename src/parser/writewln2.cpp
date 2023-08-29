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
                          Locant Path Functions
**********************************************************************/

void shift_locant_left(OBAtom **locant_path,unsigned int path_size){
  OBAtom *tmp = locant_path[0];
  for(unsigned int i=0;i<path_size-1;i++)
    locant_path[i] = locant_path[i+1]; 
  locant_path[path_size-1] = tmp; 
}

void SweepReducedPath(OBAtom** reduced_path,unsigned int path_size){
  unsigned int idx = 0; 
  for(unsigned int i=0;i<path_size;i++){
    if(reduced_path[i])
      reduced_path[idx++] = reduced_path[i];
  }
  for(unsigned int i=idx;i<path_size;i++)
    reduced_path[i] = 0; 
}

unsigned int position_in_path(OBAtom *atom,OBAtom**locant_path,unsigned int path_size){
  for(unsigned int i=0;i<path_size;i++){
    if(atom == locant_path[i])
      return i; 
  }
  fprintf(stderr,"Error: atom not found in locant path\n");
  return 0; 
}



/* standard ring walk, since nt_bonds are pre-calculated*/
OBAtom **CreateLocantPath3( OBMol *mol, unsigned int path_size,
                            std::set<OBAtom*>              &ring_atoms,
                            std::set<OBBond*>              &ring_bonds,
                            std::map<OBAtom*,unsigned int> &atom_shares,
                            std::map<OBBond*,unsigned int> &bond_shares,
                            std::set<OBRing*>    &local_SSSR, 
                            std::vector<OBBond*> &nt_bonds)
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

  OBAtom*                ratom  = 0;
  OBAtom*                catom  = 0;
  OBBond*                bond   = 0; 
  std::map<OBAtom*,bool> atoms_in_lp; 
  std::map<OBBond*,bool> ignore_bond; 
  std::stack<OBAtom*> stack; 

  for(unsigned int i=0;i<nt_bonds.size();i++)
    ignore_bond[nt_bonds[i]] = true; 

  stack.push(rseed);
  while(!stack.empty()){
    ratom = stack.top();
    stack.pop();
    locant_path[locant_pos++] = ratom; 
    atoms_in_lp[ratom] = true; 

    FOR_NBORS_OF_ATOM(a,ratom){
      catom = &(*a);
      bond = mol->GetBond(ratom,catom); 
      if(!ignore_bond[bond] && !atoms_in_lp[catom]){
        stack.push(catom); 
        break;
      }
    }
  }

  if(opt_debug){
    fprintf(stderr,"  ");
    print_locant_array(locant_path,path_size);
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


bool IsMultiCyclic(OBAtom **locant_array,unsigned int size, std::map<OBAtom*,unsigned int> &atom_shares){
  for(unsigned int i=0;i<size;i++){
    if(atom_shares[locant_array[i]] == 3)
      return true;
  }
  return false; 
}

/* 0 = not, 1 = true, 2 = partial */
unsigned int IsAromatic(std::set<OBRing*> &local_SSSR){
  unsigned int arom = 3; // just an init value  
  for(std::set<OBRing*>::iterator iter = local_SSSR.begin(); iter != local_SSSR.end(); iter++){
    if((*iter)->IsAromatic()){
      if(arom == 3)
        arom = 1; 
      else if(arom == 0){
        arom = 2;
        return arom;
      }
    }
    else{
      if(arom == 3)
        arom = 0; 
      else if(arom == 1){
        arom = 2;
        return arom;
      }
    }
  }
  return arom; 
}

/* take the created locant path and writes WLN ring assignments,
probably O(n^4) but should be limited to below 100 rings */
std::string ReadLocantPath2(  OBMol *mol, OBAtom **locant_path, unsigned int path_size,
                              std::vector<OBBond*> nt_bonds, unsigned int expected_rings)
{  
  std::string ring_str; 
  if(IsHeteroRing(locant_path,path_size))
    ring_str += 'T';
  else
    ring_str += 'L';

  // before reading the path, the last bond from the end must be calculated
  OBAtom *ratom = locant_path[path_size-1];
  for(unsigned int i=0;i<path_size;i++){
    OBAtom *catom = locant_path[i];
    OBBond *bond = mol->GetBond(ratom,catom); 
    if(bond){
      nt_bonds.push_back(bond);
      if(opt_debug)
        fprintf(stderr,"  ending path bond: %d --> %d\n",bond->GetBeginAtomIdx(),bond->GetEndAtomIdx());
      break;
    }
  }

  // working copy
  OBAtom ** reduced_path = (OBAtom**)malloc(sizeof(OBAtom*) * path_size);
  for(unsigned int i=0;i<path_size;i++)
    reduced_path[i] = locant_path[i]; 

  // sort on the lowest atom value in nt_bond, resolve ties with the other aotm
  for(unsigned int i=1;i<nt_bonds.size();i++){
    // insertion sort with tie breaks using end atom
    OBBond *key = nt_bonds[i]; 
    unsigned int first_beg  =  position_in_path(nt_bonds[i]->GetBeginAtom(),locant_path,path_size);
    int j = i - 1; 
    while(j>=0){
      unsigned int second_beg  =  position_in_path(nt_bonds[j]->GetBeginAtom(),locant_path,path_size);
      if(first_beg < second_beg)
        nt_bonds[j+1] = nt_bonds[j];
      else if (first_beg == second_beg){
        unsigned int first_end  =  position_in_path(nt_bonds[i]->GetEndAtom(),locant_path,path_size);
        unsigned int second_end  =  position_in_path(nt_bonds[j]->GetEndAtom(),locant_path,path_size);
        if(first_end < second_end)
          nt_bonds[j+1] = nt_bonds[j];
        else
          break;
      }
      else
        break;
      j--;
    }
    nt_bonds[j+1] = key;
  }

#define PATH_SHIFT 1
#ifdef PATH_SHIFT
  // move the first nt_bond to the back as WLN uses L66TJ as an init structure. 
  OBBond*tmp = nt_bonds[0];
  for(unsigned int i=0;i<nt_bonds.size()-1;i++)
    nt_bonds[i] = nt_bonds[i+1];
  nt_bonds[nt_bonds.size()-1] = tmp;
#endif

  // can we take an interrupted walk between the points, if so, write ring size 
  // and decremenent the active state
  std::map<OBAtom*,unsigned int> active_atoms; 
  for(unsigned int i=0;i<nt_bonds.size();i++){
    if(opt_debug)
      fprintf(stderr,"  bond %d: %d --> %d\n",i,nt_bonds[i]->GetBeginAtomIdx(),nt_bonds[i]->GetEndAtomIdx());
    
    active_atoms[nt_bonds[i]->GetBeginAtom()]++;
    active_atoms[nt_bonds[i]->GetEndAtom()]++;
  }

  unsigned int loops = 0; 
  while(loops < expected_rings){
    for(unsigned int i=0;i<nt_bonds.size();i++){
      OBBond *bond = nt_bonds[i]; 
      OBAtom *first = 0;
      OBAtom *second = 0; 
      unsigned int rsize = 0; 
      unsigned char first_locant; 
      for(unsigned int j=0;j<path_size;j++){
        // find the starting point in nt_bond
        if(!first || !second){
          if(reduced_path[j] == bond->GetBeginAtom()){
            first = bond->GetBeginAtom();
            second = bond->GetEndAtom();
            first_locant = int_to_locant( position_in_path(first,locant_path,path_size) +1);
            rsize++;
          }
          else if (reduced_path[j] == bond->GetEndAtom()){
            first   = bond->GetEndAtom();
            second  = bond->GetBeginAtom();
            first_locant = int_to_locant( position_in_path(first,locant_path,path_size) +1);
            rsize++;
          }
        }
        else if (first && second){

          if(reduced_path[j] != second && active_atoms[reduced_path[j]])
            break; // path interupted
          else if(reduced_path[j] == second){
            rsize++;

            if(first_locant != 'A'){
              ring_str+= ' ';
              ring_str+= first_locant;
            }

            if(rsize > 9){
              ring_str += '-';
              ring_str += std::to_string(rsize); 
              ring_str += '-';
            } 
            else
              ring_str += std::to_string(rsize);

            nt_bonds.erase(nt_bonds.begin() + i);
            // go back rsize+1 spaces and reduce the array

            for(unsigned int k = j-rsize+1;k<=j;k++){
              if(active_atoms[reduced_path[k]])
                active_atoms[reduced_path[k]]--;
              else
                reduced_path[k] = 0; 
            }
            SweepReducedPath(reduced_path,path_size);
          }
          else
            rsize++;
        }
      }

    }
    loops++;
  }

  if(opt_debug)
    fprintf(stderr,"  produced %s\n",ring_str.c_str());
  
  free(reduced_path);
  reduced_path = 0; 
  return ring_str;  
}


/**********************************************************************
                          Reduction Functions
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


/* compare ring notation to find the minimal ring direction */
bool MinimalWLNRingNotation(std::string &first_str, std::string &second_str){

  unsigned int  chain_f  =  highest_unbroken_numerical_chain(first_str); 
  unsigned char loc_f   =   first_locant_seen(first_str); 

  unsigned int  chain_s  =  highest_unbroken_numerical_chain(second_str); 
  unsigned char loc_s   =   first_locant_seen(second_str); 

  if(chain_s > chain_f)
    return true; 
  else if ((chain_s == chain_f) && loc_s < loc_f)
    return true;
  else
    return false;
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
    //      INDIRECT RECURSION TRACKING

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

  /* parses the local ring system, return the size for creating the locant path with 
  non bonds to avoid */
  unsigned int ConstructLocalSSSR(  OBMol *mol, OBAtom *ring_root,
                                    std::set<OBAtom*> &ring_atoms,
                                    std::set<OBBond*> &ring_bonds,
                                    std::map<OBAtom*,unsigned int> &atom_shares,
                                    std::map<OBBond*,unsigned int> &bond_shares,
                                    std::set<OBRing*> &local_SSSR,
                                    std::vector<OBBond*> &nt_bonds)
  {

    if(!ring_root){
      fprintf(stderr,"Error: ring root is nullptr\n");
      return 0; 
    }

    OBAtom *ratom = 0; 
    OBAtom *prev = 0; 
    OBBond *bond = 0; 
    OBRing *obring = 0; 
  
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
            if(bond){
              bond_shares[bond]++;
              ring_bonds.insert(bond); 
            }
            prev = ratom; 
          }
        }
        // get the last bond
        bond = mol->GetBond(mol->GetAtom(obring->_path.front()),mol->GetAtom(obring->_path.back()));
        if(bond){
          bond_shares[bond]++;
          ring_bonds.insert(bond); 
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

          for(unsigned int i=0;i<obring->Size();i++){
            OBAtom *ratom = mol->GetAtom(obring->_path[i]);
            ring_set.insert(ratom);
          }

          std::set_intersection(ring_set.begin(), ring_set.end(), ring_atoms.begin(), ring_atoms.end(),
                                std::inserter(intersection, intersection.begin()));

          // intersection 1 is a spiro ring, ignore, 
          if(intersection.size() > 1){
            rings_seen[obring]; 
            local_SSSR.insert(obring);

            prev = 0;
            for(unsigned int i=0;i<obring->Size();i++){
              OBAtom *ratom = mol->GetAtom(obring->_path[i]);
              ring_atoms.insert(ratom);
              atom_shares[ratom]++;

              if(!prev)
                prev = ratom;
              else{
                bond = mol->GetBond(prev,ratom);
                if(bond){
                  bond_shares[bond]++;
                  ring_bonds.insert(bond); 
                }
                prev = ratom; 
              }
            }
            // get the last bond
            bond = mol->GetBond(mol->GetAtom(obring->_path.front()),mol->GetAtom(obring->_path.back()));
            if(bond){
              bond_shares[bond]++;
              ring_bonds.insert(bond); 
            }
              
            rings_seen[obring] = true;
            running = true;
          }
        }
      }
    }

    for(std::set<OBBond*>::iterator biter = ring_bonds.begin(); biter != ring_bonds.end(); biter++){
      bond = (*biter);
      if(bond_shares[bond] > 1){
        nt_bonds.push_back(bond);
        if(opt_debug)
          fprintf(stderr,"  locant path bond: %d --> %d\n",bond->GetBeginAtomIdx(),bond->GetEndAtomIdx());
      }
    }

    if(opt_debug){
      fprintf(stderr,"  ring atoms: %lu\n",ring_atoms.size());
      fprintf(stderr,"  ring bonds: %lu\n",ring_bonds.size());
    }
      
    return ring_atoms.size(); 
  }


  /* create the hetero atoms where neccesary */
  bool ReadLocantAtomsBonds( OBAtom** locant_path,unsigned int path_size,std::string &buffer)
  {
    unsigned char locant = 0;
    unsigned char last_locant = 0; 
    
    OBAtom *first   = 0;
    OBAtom *second  = 0; 
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

#ifdef WIP
      // handles sequential locant unsaturations, when not aromatic
      first = locant_path[i];
      if(i < path_size-1)
        second = locant_path[i+1];
      else
        second = locant_path[0];

      OBBond *locant_bond = first->GetBond(second);
      if(locant_bond && locant_bond->GetBondOrder() > 1){
        if(i > 0 && locant != last_locant){
          buffer += ' ';
          buffer += locant;
        }
        for(unsigned int b=1;b<locant_bond->GetBondOrder();b++)
          buffer += 'U';
      }
#endif
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
        append+= int_to_locant(i+1);
      }
    }

    buffer += ' ';
    buffer+= std::to_string(count);
    buffer+= append; 
  }


  void ReadAromaticity(std::vector<OBRing*> &ring_order,std::string &buffer){

    std::string append = ""; 
    for(unsigned int i=0;i<ring_order.size();i++){
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
   
    std::set<OBRing*>               local_SSSR; 
    std::set<OBAtom*>               ring_atoms;
    std::set<OBBond*>               ring_bonds;
    std::vector<OBBond*>            nt_bonds;
    std::map<OBAtom*,unsigned int>  atom_shares;
    std::map<OBBond*,unsigned int>  bond_shares;
    
    unsigned int path_size   =  ConstructLocalSSSR(mol,ring_root,ring_atoms,ring_bonds,atom_shares,bond_shares,local_SSSR,nt_bonds); 
    OBAtom **    locant_path =  CreateLocantPath3(mol,path_size,ring_atoms,ring_bonds,atom_shares,bond_shares,local_SSSR,nt_bonds);
    
    // monocyclic
    if(local_SSSR.size() == 1){
      OBRing* monoring = *(local_SSSR.begin());
      for(unsigned int i=0;i<monoring->Size();i++){
        locant_path[i] = mol->GetAtom(monoring->_path[i]);
        if(inline_ring && locant_path[i] == ring_root)
          buffer += int_to_locant(i+1);
      }
        
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
      
      ReadLocantAtomsBonds(locant_path,path_size,buffer);
      if(!monoring->IsAromatic())
        buffer += 'T';
      buffer += 'J';
      return {locant_path,path_size};
    }
    // theres a slightly different procedure for minimising Multicyclic compounds
    else if(!IsMultiCyclic(locant_path,path_size,atom_shares)){
      
      OBAtom *start = locant_path[0];
      std::string cyclic_str = ReadLocantPath2(mol,locant_path,path_size,nt_bonds,local_SSSR.size()); 
      
      do{
        shift_locant_left(locant_path,path_size);
        
      }while(atom_shares[locant_path[0]] != 2);
      print_locant_array(locant_path,path_size);

      cyclic_str = ReadLocantPath2(mol,locant_path,path_size,nt_bonds,local_SSSR.size()); 
      
      buffer += cyclic_str;

      ReadLocantAtomsBonds(locant_path,path_size,buffer);
      
      unsigned int aromatic = IsAromatic(local_SSSR);
      if(!aromatic)
        buffer += "TJ";
      else if (aromatic == 1)
        buffer += 'J';
      else{
        fprintf(stderr,"Warning: mixed aromaticity not implemented yet\n");
        buffer += 'J';
      }

      return {locant_path,path_size};
    }
    else{
      fprintf(stderr,"Error: could not form locant path\n");
      return {0,0};
    }
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


