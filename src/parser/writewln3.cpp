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

#include <set>
#include <deque>
#include <vector>
#include <stack>
#include <map>

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

// --- macros ---
#define REASONABLE 1024


// --- inputs ---
const char *cli_inp;
const char *dotfile;

// --- options ---
static bool opt_wln2dot = false;
static bool opt_allow = false;
static bool opt_debug = false;
static bool opt_convert = false;

// --- globals ---
const char *wln;
struct WLNSymbol;
struct WLNEdge; 
struct WLNRing;

unsigned int edge_count   = 0;
unsigned int symbol_count = 0;
unsigned int ring_count   = 0;

WLNSymbol *SYMBOLS[REASONABLE];
WLNEdge   *EDGES  [REASONABLE];
WLNRing   *RINGS  [REASONABLE];

std::map<WLNSymbol *, unsigned int> index_lookup;
std::map<unsigned int, WLNSymbol *> symbol_lookup;
std::map<unsigned int,OpenBabel::OBAtom*> babel_atom_lookup;

unsigned int glob_index = 1; // babel starts from 1, keep consistent  

// ionic parsing
std::map<unsigned int,WLNSymbol*> string_positions; 
std::map<WLNSymbol*,int> charge_additions;


enum WLNTYPE
{
  STANDARD = 0,
  CHAIN = 1, 
  LOCANT = 2,   
  RING = 3,     
  ELEMENT = 4
};


// rule 2 - hierarchy - rules have diverged due to end terminator char, also use for locant setting from 14
std::map<unsigned char, unsigned int> char_hierarchy =
    {
      {' ', 1}, {'-', 2}, {'/', 3}, {'0', 4}, {'1', 5}, {'2', 6}, {'3', 7}, {'4', 8}, {'5', 9}, {'6', 10}, {'7', 11}, {'8', 12}, {'9', 13}, {'A', 14}, {'B', 15}, {'C', 16}, {'D', 17}, {'E', 18}, {'F', 19}, {'G', 20}, {'H', 21}, {'I', 22}, {'J', 23}, {'K', 24}, {'L', 25}, {'M', 26}, {'N', 27}, {'O', 28}, {'P', 29}, {'Q', 30}, {'R', 31}, {'S', 32}, {'T', 33}, {'U', 34}, {'V', 35}, {'W', 36}, {'X', 37}, {'Y', 38}, {'Z', 40}, {'&', 41}};


unsigned char static int_to_locant(unsigned int i){
  return i + 64;
}

unsigned int static locant_to_int(unsigned char loc){
  return loc - 64;
}


std::string get_notation(unsigned int s, unsigned int e)
{
  std::string res; 
  for (unsigned int i = s; i <= e; i++)
  {
    res.push_back(wln[i]);
  }
  return res; 
}

static void Fatal(unsigned int pos)
{
  fprintf(stderr, "Fatal: %s\n", wln);
  fprintf(stderr, "       ");

  for (unsigned int i = 0; i < pos; i++)
    fprintf(stderr, " ");

  fprintf(stderr, "^\n");

  exit(1);
}



struct WLNEdge{
  WLNSymbol *parent;
  WLNSymbol *child;
  WLNEdge *nxt;

  bool aromatic;
  unsigned int order;

  WLNEdge(){
    parent = 0;
    child = 0;
    aromatic = 0;
    order = 0;
    nxt = 0;
  }
};


struct WLNSymbol
{
  unsigned char ch;
  std::string special; // string for element, or ring, if value = '*'

  unsigned int type;
  unsigned int allowed_edges;
  unsigned int num_edges;

  WLNSymbol *previous;
  WLNEdge   *bonds; // array of bonds

  // if default needed
  WLNSymbol()
  {
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
    previous = 0;
  }
  ~WLNSymbol(){};

  void set_edge_and_type(unsigned int e, unsigned int t=STANDARD){
    allowed_edges = e;
    type = t;
  }

  void add_special(unsigned int s, unsigned int e)
  {
    for (unsigned int i = s; i <= e; i++)
      special.push_back(wln[i]);
  }

};


WLNSymbol *AllocateWLNSymbol(unsigned char ch)
{

  
  symbol_count++;
  if(symbol_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wlln symbols - is this reasonable?\n");
    return 0;
  }
  
  WLNSymbol *wln = new WLNSymbol;

  SYMBOLS[symbol_count] = wln;
  
  wln->ch = ch;
  index_lookup[wln] = glob_index;
  symbol_lookup[glob_index] = wln;
  glob_index++;
  return wln;
}


WLNEdge *AllocateWLNEdge(WLNSymbol *child, WLNSymbol *parent){

  if(!child || !parent){
    fprintf(stderr,"Error: attempting bond of non-existent symbols - %s|%s is dead\n",child ? "":"child",parent ? "":"parent");
    return 0;
  }

  edge_count++;
  if(edge_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
    return 0;
  }
  
  if ((child->num_edges + 1) > child->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges+1, child->allowed_edges);
    return 0;
  }
  
  if ((parent->num_edges + 1) > parent->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges+1, parent->allowed_edges);
    return 0;
  }

  WLNEdge *edge = new WLNEdge;
  EDGES[edge_count] = edge;

  // use a linked list to store the bond, can also check if it already exists
  if(parent->bonds){
    WLNEdge *curr = parent->bonds;

    while(curr->nxt){
      if(curr->child == child){
        fprintf(stderr,"Error: trying to bond already bonded symbols\n");
        return 0;
      }
      curr = curr->nxt;
    }
      
    curr->nxt = edge;
  }
  else
    parent->bonds = edge; 

  // set the previous for look back
  child->previous = parent; 

  child->num_edges++;
  parent->num_edges++;

  edge->parent = parent; 
  edge->child = child;
  edge->order = 1;
  return edge;
}

WLNEdge *search_edge(WLNSymbol *child, WLNSymbol*parent){
  WLNEdge *edge = 0;
  for (edge=parent->bonds;edge;edge = edge->nxt){
    if(edge->child == child)
      return edge;
  }
  fprintf(stderr,"Error: could not find edge in search\n");
  return 0;
}

WLNEdge *unsaturate_edge(WLNEdge *edge,unsigned int n){
  if(!edge){
    fprintf(stderr,"Error: unsaturating non-existent edge\n");
    return 0;
  }

  edge->order += n; 
  edge->parent->num_edges += n;
  edge->child->num_edges+= n;

  if(edge->parent->num_edges > edge->parent->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->parent->ch,edge->parent->num_edges, edge->parent->allowed_edges);
    return 0;
  }

  if(edge->child->num_edges > edge->child->allowed_edges){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->child->ch,edge->child->num_edges, edge->child->allowed_edges);
    return 0;
  }

  return edge;
}


bool remove_edge(WLNSymbol *head,WLNEdge *edge){
  if(!head || !edge){
    fprintf(stderr,"Error: removing bond of non-existent symbols\n");
    return false;
  }

  if(head->bonds == edge){
    head->bonds = 0;
    return true;
  }

  bool found = false;
  WLNEdge *search = head->bonds;

  WLNEdge *prev = 0;
  while(search){
    if(search == edge){ 
      found = true;
      break;
    }
    prev = search; 
    search = search->nxt;
  }


  if(!found){
    fprintf(stderr,"Error: trying to remove bond from wln character[%c] - bond not found\n",head->ch);
    return false;
  }
  else{
    WLNEdge *tmp = edge->nxt;
    prev->nxt = tmp;
    // dont null the edges as we use the mempool to release them
  }

  return true;
}


WLNEdge* add_methyl(WLNSymbol *head){

  WLNSymbol *carbon = AllocateWLNSymbol('C');
  carbon->set_edge_and_type(4); // used for hydrogens

  for(unsigned int i=0;i<3;i++){
    WLNSymbol *hydrogen = AllocateWLNSymbol('H');
    WLNEdge   *edge = AllocateWLNEdge(hydrogen,carbon);
  }

  WLNEdge *bond = AllocateWLNEdge(carbon,head);
  return bond; 
}


bool add_diazo(WLNSymbol *head){

  WLNEdge *edge = 0;
  WLNSymbol *oxygen = 0;

  oxygen = AllocateWLNSymbol('O');
  oxygen->set_edge_and_type(2,head->type);
  charge_additions[oxygen] = -1;

  edge = AllocateWLNEdge(oxygen,head);
  
  
  oxygen = AllocateWLNSymbol('O');
  oxygen->set_edge_and_type(2,head->type);
  edge = AllocateWLNEdge(oxygen,head);
  
  edge = unsaturate_edge(edge,1);
  if(!edge)
    return false;
  
  return true;
}


/* resolve carbon methyl assumptions */
bool resolve_methyls(WLNSymbol *target){

  switch(target->ch){

    case 'Y':
    case 'X':
    case 'K':
      while(target->num_edges < target->allowed_edges){
        if(!add_methyl(target))
          return false;
      }
      target->num_edges = target->allowed_edges;
      break;

    default:
      fprintf(stderr,"Error: resolving methyls performed on invalid symbol: %c\n",target->ch);
      return false;
  }

  return true;
}


WLNSymbol* define_hypervalent_element(unsigned char sym,WLNSymbol *override=0){

  WLNSymbol *new_symbol = 0;
  if(!override)
    new_symbol = AllocateWLNSymbol(sym);
  else{
    new_symbol = override;
    new_symbol->ch = sym;
  }

  switch(sym){
    
    case 'P':
    case 'S':
    case 'G':
    case 'E':
    case 'I':
    case 'F':
      new_symbol->set_edge_and_type(6);            // allows FCl6
      return new_symbol;

    default:
      fprintf(stderr,"Error: character %c does not need - notation for valence expansion, please remove -\n",sym);
      return 0;
  }
  
  return new_symbol;
}

/* allocate new or override exisiting node*/
WLNSymbol* define_element(std::string special,WLNSymbol *remake = 0){
    
  WLNSymbol *created_wln = 0;

  if(!remake)
    created_wln = AllocateWLNSymbol('*');
  else{
    created_wln = remake;
    created_wln->ch = '*';
  }
   

  created_wln->allowed_edges = 18; // allow anything for now;

  switch (special[0]){

    case 'A':
      switch(special[1]){
        case 'C':
        case 'G':
        case 'L':
        case 'M':
        case 'R':
        case 'S':
        case 'T':
        case 'U':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'B':
      switch(special[1]){
        case 'A':
        case 'E':
        case 'H':
        case 'I':
        case 'K':
        case 'R':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      

    case 'C':
      switch(special[1]){
        case 'A':
        case 'D':
        case 'E':
        case 'F':
        case 'M':
        case 'N':
        case 'O':
        case 'R':
        case 'S':
        case 'U':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      
    case 'D':
      switch(special[1]){
        case 'B':
        case 'S':
        case 'Y':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'E':
      switch(special[1]){
        case 'R':
        case 'S':
        case 'U':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'F':
      switch(special[1]){
        case 'E':
        case 'L':
        case 'M':
        case 'R':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'G':
      switch(special[1]){
        case 'A':
        case 'D':
        case 'E':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'H':
      switch(special[1]){
        case 'E':
        case 'F':
        case 'G':
        case 'O':
        case 'S':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'I':
      switch(special[1]){
        case 'N':
        case 'R':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'K':
      if (special[1] == 'R')
        created_wln->special = "Kr";
      else if (special[1] == 'A'){
        created_wln->ch = 'K';
        created_wln->set_edge_and_type(4);
        return created_wln;
      }
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;

    case 'L':
      switch(special[1]){
        case 'A':
        case 'I':
        case 'R':
        case 'U':
        case 'V':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'M':
      switch(special[1]){
        case 'C':
        case 'D':
        case 'G':
        case 'N':
        case 'O':
        case 'T':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'N':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'D':
        case 'E':
        case 'H':
        case 'I':
        case 'O':
        case 'P':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }


    case 'O':
      switch(special[1]){
        case 'O':
        case 'G':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'P':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'D':
        case 'M':
        case 'O':
        case 'R':
        case 'T':
        case 'U':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'R':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'n':
        case 'U':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
     

    case 'S':
      switch(special[1]){
        case 'B':
        case 'C':
        case 'E':
        case 'G':
        case 'I':
        case 'M':
        case 'N':
        case 'R':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }


    case 'T':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'C':
        case 'E':
        case 'H':
        case 'I':
        case 'L':
        case 'M':
        case 'S':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }

    case 'U':
      if(special[1] == 'R')
        created_wln->special = "UR";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;

    case 'V':
      if (special[1] == 'A')
        created_wln->special = "VA";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;
    
    case 'W':
      if(special[1] == 'T')
        created_wln->special = "WT";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;
    

    case 'X':
      if (special[1] == 'E')
        created_wln->special = "XE";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;

    case 'Y':
      switch(special[1]){
        case 'B':
        case 'T':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'Z':
      switch(special[1]){
        case 'N':
        case 'R':
          created_wln->special = special;
          return created_wln;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    default:
      fprintf(stderr, "Error: invalid character in special definition switch\n");
      return (WLNSymbol *)0;
  }

  return 0;
}


/* checks are already made, this should just return*/
unsigned int special_element_atm(std::string &special){
  
  switch (special[0]){

    case 'A':
      if (special[1] == 'C')
        return 89;
      else if (special[1] == 'G')
        return 47;
      else if (special[1] == 'L')
        return 13;
      else if (special[1] == 'M')
        return 95;
      else if (special[1] == 'R')
        return 18;
      else if (special[1] == 'S')
        return 33;
      else if (special[1] == 'T')
        return 85;
      else if (special[1] == 'U')
        return 79;
      break;

    case 'B':
      if (special[1] == 'A')
        return 56;
      else if (special[1] == 'E')
        return 4;
      else if (special[1] == 'H')
        return 107;
      else if (special[1] == 'I')
        return 83;
      else if (special[1] == 'K')
        return 97;
      else if (special[1] == 'R')
        return 35;
      break;

    case 'C':
      if (special[1] == 'A')
        return 20;
      else if (special[1] == 'D')
        return 48;
      else if (special[1] == 'E')
        return 58;
      else if (special[1] == 'F')
        return 98;
      else if (special[1] == 'M')
        return 96;
      else if (special[1] == 'N')
        return 112;
      else if (special[1] == 'O')
        return 27;
      else if (special[1] == 'R')
        return 24;
      else if (special[1] == 'S')
        return 55;
      else if (special[1] == 'U')
        return 29;
      break;

    case 'D':
      if (special[1] == 'B')
        return 105;
      else if (special[1] == 'S')
        return 110;
      else if (special[1] == 'Y')
        return 66;
      break;

    case 'E':
      if (special[1] == 'R')
        return 68;
      else if (special[1] == 'S')
        return 99; 
      else if (special[1] == 'U')
        return 63;
      break;

    case 'F':
      if (special[1] == 'E')
        return 26;
      else if (special[1] == 'L')
        return 114;
      else if (special[1] == 'M')
        return 100;
      else if (special[1] == 'R')
        return 97;
      break;

    case 'G':
      if (special[1] == 'A')
        return 31;
      else if (special[1] == 'D')
        return 64;
      else if (special[1] == 'E')
        return 32;
      break;

    case 'H':
      if (special[1] == 'E')
        return 2;
      else if (special[1] == 'F')
        return 72;
      else if (special[1] == 'G')
        return 80;
      else if (special[1] == 'O')
        return 67;
      else if (special[1] == 'S')
        return 108;

      break;

    case 'I':
      if (special[1] == 'N')
        return 49;
      else if (special[1] == 'R')
        return 77;
      break;

    case 'K':
      if (special[1] == 'R')
        return 39;
      break;

    case 'L':
      if (special[1] == 'A')
        return 57;
      else if (special[1] == 'I')
        return 3;
      else if (special[1] == 'R')
        return 103;
      else if (special[1] == 'U')
        return 71;
      else if (special[1] == 'V')
        return 116;
      break;

    case 'M':
      if (special[1] == 'C')
        return 115;
      else if (special[1] == 'D')
        return 101;
      else if (special[1] == 'G')
        return 12;
      else if (special[1] == 'N')
        return 25;
      else if (special[1] == 'O')
        return 42;
      else if (special[1] == 'T')
        return 109;
      break;

    case 'N':
      if (special[1] == 'A')
       return 11;
      else if (special[1] == 'B')
        return 41;
      else if (special[1] == 'D')
        return 60;
      else if (special[1] == 'E')
        return 10;
      else if (special[1] == 'H')
        return 113;
      else if (special[1] == 'I')
        return 28;
      else if (special[1] == 'O')
        return 102;
      else if (special[1] == 'P')
        return 93;
      break;

    case 'O':
      if (special[1] == 'G')
        return 118;
      else if (special[1] == 'S')
        return 76;
      break;

    case 'P':
      if (special[1] == 'A')
        return 91;       
      else if (special[1] == 'B')
        return 82;
      else if (special[1] == 'D')
        return 46;
      else if (special[1] == 'M')
        return 61;
      else if (special[1] == 'O')
        return 84;
      else if (special[1] == 'R')
        return 59;
      else if (special[1] == 'T')
        return 97;
      else if (special[1] == 'U')
        return 94;
      
      break;

    case 'R':
      if (special[1] == 'A')
        return 88;
      else if (special[1] == 'B')
        return 37;
      else if (special[1] == 'E')
        return 75;
      else if (special[1] == 'F')
        return 104;
      else if (special[1] == 'G')
        return 111;
      else if (special[1] == 'H')
        return 45;
      else if (special[1] == 'N')
        return 86;
      else if (special[1] == 'U')
        return 44;
      break;

    case 'S':
      if (special[1] == 'B')
        return 51;
      else if (special[1] == 'C')
        return 21;
      else if (special[1] == 'E')
        return 34;
      else if (special[1] == 'G')
        return 106;
      else if (special[1] == 'I')
        return 14;
      else if (special[1] == 'M')
        return 62;
      else if (special[1] == 'N')
        return 50;
      else if (special[1] == 'R')
        return 38;
      
      break;

    case 'T':
      if (special[1] == 'A')
        return 73;
      else if (special[1] == 'B')
        return 65;
      else if (special[1] == 'C')
        return 43;
      else if (special[1] == 'E')
        return 52;
      else if (special[1] == 'H')
        return 90;
      else if (special[1] == 'I')
        return 22;
      else if (special[1] == 'L')
        return 81;
      else if (special[1] == 'M')
        return 69;
      else if (special[1] == 'S')
        return 117;

      break;

    case 'U':
      if(special[1] == 'R')
        return 92;
      break;

    case 'V':
      if(special[1] == 'A')
        return 23;
      break;

    case 'X':
      if (special[1] == 'E')
        return 54;
      break;

    case 'Y':
      if(special[1] == 'T')
        return 39;
      else if (special[1] == 'B')
        return 70;
      break;

    case 'Z':
      if (special[1] == 'N')
        return 30;
      else if (special[1] == 'R')
        return 40;
  
      break;

    default:
      fprintf(stderr, "Error: invalid character in special definition switch\n");
      return 0;
  }

  return 0;
}


/* struct to hold pointers for the wln ring */
struct WLNRing
{
  std::vector<unsigned int> rings;
  std::map<unsigned char, WLNSymbol *> locants;
  std::map<WLNSymbol*,unsigned char> locants_ch;

  // some bridge notation is index dependent
  struct indexed_pair{
    unsigned char bind_1 = '\0';
    unsigned char bind_2 = '\0';
    unsigned int index   = 0;

    void set(unsigned char a, unsigned char b, unsigned int p){
      bind_1 = a;
      bind_2 = b;
      index = p;
    }

  };

  WLNRing(){}
  ~WLNRing(){};

  // both lookups needed for QOL in ring building
  WLNSymbol* assign_locant(unsigned char loc,unsigned char type){
    WLNSymbol *locant = 0; 
    locant = AllocateWLNSymbol(type);
    locants[loc] = locant; 
    locants_ch[locant] = loc;
    
    locant->type = RING;
    return locant; 
  }

  // both lookups needed for QOL in ring building
  WLNSymbol* assign_locant(unsigned char loc,WLNSymbol *locant){
    
    if(!locant)
      return 0;
    
    locants[loc] = locant; 
    locants_ch[locant] = loc;
    locant->type = RING;
    return locant; 
  }  

  void debug_locants(){
    std::map<unsigned char, WLNSymbol *>::iterator locant_itr;
    fprintf(stderr,"alive locants: ");
    for (locant_itr = locants.begin(); locant_itr != locants.end(); locant_itr++){
      if(locant_itr->second){
        fprintf(stderr," %c",locant_itr->first);
      }
    }
    fprintf(stderr,"\n");
  }


  /* creates poly rings, aromaticity is defined in reverse due to the nature of notation build */
  bool CreatePolyCyclic(std::vector<std::pair<unsigned int,unsigned char>> &ring_assignments, 
                  std::vector<bool> &aromaticity)
    {
     

    unsigned int local_size = 0; 
    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i]; 

      if(local_size)
        local_size += component.first - 2;
      else
        local_size = component.first;
    }

    // create all the nodes in a large straight chain
    
    WLNSymbol *curr= 0; 
    WLNSymbol *prev = 0; 
    for (unsigned int i=1;i<=local_size;i++){
      unsigned char loc = int_to_locant(i);
      if(!locants[loc]){
        curr = assign_locant(loc,'C');
        curr->set_edge_and_type(4,RING);
      }
      else
        curr = locants[loc];

      if(prev){
        WLNEdge *edge = AllocateWLNEdge(curr,prev);
        if(!edge)
          return false;
      }
      prev = curr;
    }


    // calculate bindings and then traversals round the loops
    unsigned int comp_size = 0;
    unsigned char bind_1 = '\0';
    unsigned char bind_2 = '\0';
    unsigned int fuses = 0; 
    bool aromatic = false; 


    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i];
      comp_size = component.first;
      bind_1 = component.second;
      aromatic = aromaticity[i];

      std::deque<unsigned char> ring_path;
      // first pair can be calculated directly without a path travel
      if(!fuses){
        bind_2 = bind_1 + comp_size - 1; // includes start atom
        for (unsigned int i=0; i<comp_size;i++)
          ring_path.push_back(bind_1+i);
      }
      else{
        //there needs to be a graph travel here taking the longest locant

        // 1. starting on bind_1, travel n-1 places through the maximal locant path, to calculate fuse
        
        // annoyingly n2 ... 
        WLNSymbol *path = locants[bind_1];
        unsigned char highest_loc = '\0';
        for (unsigned int i=0;i<comp_size - 1; i++){
          ring_path.push_back(locants_ch[path]);

          WLNEdge *lc = 0;
          for(lc = path->bonds;lc;lc = lc->nxt){
            WLNSymbol *child = lc->child;
            unsigned char child_loc = locants_ch[child];
            if(child_loc > highest_loc)
              highest_loc = child_loc;
          }    
          path = locants[highest_loc];
        }

        ring_path.push_back(locants_ch[path]); // add the last symbol
        bind_2 = highest_loc;
      }

      if(opt_debug){
        fprintf(stderr,"  %d  fusing: %c <-- %c   [",fuses,bind_2,bind_1);
        for (unsigned char ch : ring_path){
          fprintf(stderr," %c(%d)",ch,ch);
        }
        fprintf(stderr," ]\n");
      }

      WLNEdge *edge = AllocateWLNEdge(locants[bind_2],locants[bind_1]);
      if(!edge)
        return false;

      // if(aromatic){
      //   if(!AssignAromatics(ring_path))
      //     return false;
      // }
        
      fuses++;
    }

    return true; 
  }


  /* interesting here that the multicyclic points are not explicitly used */
  bool CreateMultiCyclic(std::vector<std::pair<unsigned int,unsigned char>> &ring_assignments, 
                  std::vector<bool> &aromaticity,
                  std::vector<unsigned char> &multicyclic_locants,
                  std::vector<indexed_pair> &pseudo_locants,
                  std::set<unsigned char> &broken_locants,
                  unsigned char size_designator)
  {
    
    // create a chain size of ring designator
    unsigned int local_size = locant_to_int(size_designator);
    std::map<unsigned char, unsigned int> rings_shared;
  
    // create all the nodes in a large straight chain
    WLNSymbol *curr = 0; 
    WLNSymbol *prev = 0; 
    for (unsigned int i=1;i<=local_size;i++){
      unsigned char loc = int_to_locant(i);
      if(!locants[loc]){
        curr = assign_locant(loc,'C');
        curr->set_edge_and_type(4,RING);
      }
      else
        curr = locants[loc];

      if(prev){
        WLNEdge *edge = AllocateWLNEdge(curr,prev);
        if(!edge)
          return false;
      }
      prev = curr;
    }


    // have these as indexed lookups in the component pass

    // kick back if pairs are occupying the same index
    std::map<unsigned int,std::vector<indexed_pair>> pseudo_lookup;
    for (indexed_pair psd_pair : pseudo_locants){
      pseudo_lookup[psd_pair.index].push_back(psd_pair);
    }

    // bind the broken locant to its parents character 'J- will bond to 'J'

    // broken locants have an intergration marker, they can only be considered in 
    // the path after all prior rings have been evaulated to that point

    // parent -> all dead ends, e.g 'B' --> {B-, B-&, B--, B--&}
    std::map<unsigned char,std::vector<unsigned char>> broken_lookup;
    std::map<unsigned char,bool>                       resolved;

    if(!broken_locants.empty()){
      // create the atoms, 
      for (unsigned char loc_broken : broken_locants){
        unsigned char calculate_origin = loc_broken;
        unsigned int pos = 0;
        while( (calculate_origin - 23) > 128){
          calculate_origin += -23;
          pos++;
        }

        // position here decodes where to link them

        unsigned char parent = '\0';
        parent = int_to_locant(128 + calculate_origin); // relative positioning
        if(pos == 2 || pos == 3)
          parent = locant_to_int(parent) + 128;
        else if(pos > 3){
          fprintf(stderr,"Error: non-locant links past a two-level tree are unsuitable for this parser\n");
          return false;
        }

        if(opt_debug){
          fprintf(stderr,"  ghost linking %d to parent %d\n",loc_broken,parent);
        }

        if(!locants[loc_broken]){
          WLNSymbol *broken = assign_locant(loc_broken,'C');
          broken->set_edge_and_type(4,RING);

          broken_lookup[parent].push_back(loc_broken);
          resolved[loc_broken] = false;
        }
        else{
          fprintf(stderr,"Error: branching locants are overlapping created elements already in the locant path\n");
          return false;
        }
      }
    }

    // calculate bindings and then traversals round the loops
    unsigned int comp_size = 0;
    unsigned char bind_1 = '\0';
    unsigned char bind_2 = '\0';
    unsigned int fuses = 0; 
    bool aromatic = false;

    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i];
      comp_size = component.first;
      bind_1 = component.second;

      aromatic = aromaticity[i];
      WLNSymbol *path = locants[bind_1];

      std::deque<unsigned char> ring_path; // need push and pop function at both ends
      unsigned int predefined = 1;

      if(path->num_edges > 2)
        predefined++; 


      // GIVEN BRIDGE LOCANTS ONLY
      
      if(!pseudo_lookup[i].empty()){
        indexed_pair psd_pair = pseudo_lookup[i].front(); // should only be 1
        
        bind_1 = psd_pair.bind_1;
        bind_2 = psd_pair.bind_2;

        WLNEdge *edge = AllocateWLNEdge(locants[bind_2],locants[bind_1]);
        if(!edge)
          return false;
        
        // a ring path then needs to be calculated for aromaticity assignment
        // should be only one path

        WLNSymbol *path = locants[bind_1];
        unsigned char highest_loc = '\0'; // highest of each child iteration 

        ring_path.push_back(locants_ch[path]);
        for (unsigned int i=0;i<comp_size - predefined; i++){
          WLNEdge *lc = 0;
          for(lc = path->bonds;lc;lc = lc->nxt){
            WLNSymbol *child = lc->child;
            unsigned char child_loc = locants_ch[child];

            // cannot be the bound symbol
            if(child_loc == bind_2)
              continue;

            if(child_loc >= highest_loc)
              highest_loc = child_loc;
          }

          path = locants[highest_loc];
          ring_path.push_back(locants_ch[path]);
        }

        if(ring_path.back() != bind_2)
          ring_path.push_back(bind_2);

        if(opt_debug){
          fprintf(stderr,"  %d  fusing: %c <-- %c   [",fuses,bind_2,bind_1);
          for (unsigned char ch : ring_path){
            fprintf(stderr," %c(%d)",ch,ch);
          }
          fprintf(stderr," ]\n");
        }

        fuses++;
       
        // skip the proper algorithmic block
        continue;      
      }
      


      // MULTI ALGORITHM

      ring_path.push_back(locants_ch[path]);

      unsigned char highest_loc = '\0'; 
      for (unsigned int i=0;i<comp_size - predefined; i++){  // 2 points are already defined
        
        highest_loc = '\0'; // highest of each child iteration 

        WLNEdge *lc = 0;
        for(lc = path->bonds;lc;lc = lc->nxt){
          WLNSymbol *child = lc->child;
          unsigned char child_loc = locants_ch[child];

          if(child_loc >= highest_loc)
            highest_loc = child_loc;
        }

        if(!highest_loc){
          
          if(locant_to_int(locants_ch[path]) == local_size)
            highest_loc = locants_ch[path];
          else{
            fprintf(stderr,"Error: locant path formation is broken in ring definition - '%c'\n",locants_ch[path]);
            return false;
          }
        }

        path = locants[highest_loc];
        ring_path.push_back(locants_ch[path]);
      }

      bind_2 = highest_loc; 

      if(locants[bind_1]->num_edges > 2){
        bool shift = true;
        if(!broken_lookup[bind_1].empty()){
          // check if maxed out, if so - shift as standard
          for (unsigned char extra :broken_lookup[bind_1]){
            if(locants[extra]->num_edges < 3){
              shift = false;
              break;
            }
          }
        }
        if(shift){
          bind_1 += 1; // handles all normal  multicyclic denototions
          while(locants[bind_1]->num_edges > 2){
            bind_1 += 1;
          }
          ring_path.push_front(bind_1);
        }
      }

      // check are we going to make this a multi point with a look up?
      if(locants[bind_1]->num_edges >= 2 && !broken_lookup[bind_1].empty()){

        // while loop allows multiple decends if needed
        while(!broken_lookup[bind_1].empty()){
          
          WLNSymbol *broken = 0;
          for (unsigned int k=0;k<broken_lookup[bind_1].size();k++){
            broken = locants[broken_lookup[bind_1].at(k)];
            if(!resolved[broken_lookup[bind_1].at(k)]){
              resolved[broken_lookup[bind_1].at(k)] = true;
              break;
            }
              
          }


          WLNEdge *edge = AllocateWLNEdge(broken,locants[bind_1]);
          if(!edge)
            return false;

          // since we've added a symbol to path we need to push and pop
          ring_path.push_front(locants_ch[broken]);
          bind_1 = locants_ch[broken];
        }

        while(ring_path.size() != comp_size)
          ring_path.pop_back();
      
        // the new path should be set
        bind_1 = ring_path.front();
        bind_2 = ring_path.back();
      }


      // annoying catch needed for bridge notation that is 'implied' 
      if(i == ring_assignments.size() - 1 && bind_2 != int_to_locant(local_size)){
        unsigned char back = ring_path.back();
        while(back != int_to_locant(local_size)){
          back++;
          ring_path.pop_front();
        }
        bind_2 = back;
        bind_1 = ring_path.front();
      }
      

      if(opt_debug){
        fprintf(stderr,"  %d  fusing: %c <-- %c   [",fuses,bind_2,bind_1);
        for (unsigned char ch : ring_path){
          fprintf(stderr," %c(%d)",ch,ch);
        }
        fprintf(stderr," ]\n");
      }

      WLNEdge *edge = AllocateWLNEdge(locants[bind_2],locants[bind_1]);
      if(!edge)
        return false;
        
      // if(aromatic){
      //   if(!AssignAromatics(ring_path))
      //     return false;
      // }

      fuses++;
    }


    return true; 
  }

  unsigned int FindAromatics(const char *wln_ptr, unsigned int len){
    int i = 0;
    for(i=len-1;i>-1;i--){
      unsigned char ch = wln_ptr[i];
      if(ch == 'J')
        continue;
      else if(ch == 'T')
        continue;
      else if(ch == '&')
        continue;
      else if(ch == '-')
        return i; // this so we can ignore it in the notation 
      else 
        break;
    }
    return i+1;
  };


  unsigned char create_relative_position(unsigned char parent){
    // A = 129
    unsigned int relative = 128 + locant_to_int(parent);
    if(relative > 252){
      fprintf(stderr,"Error: relative position is exceeding 252 allowed space - is this is suitable molecule for WLN notation?\n");
      return '\0';
    }
    else
      return relative;
  }

  
  void FormWLNRing(std::string &block, unsigned int start){

  
    enum RingType{ POLY=1, PERI=2, BRIDGED=3, PSDBRIDGED = 4}; 
    const char* ring_strings[] = {"MONO","POLY","PERI","BRIDGED","PSDBRIDGED"};
    unsigned int ring_type = POLY;   // start in mono and climb up

    bool warned             = false;  // limit warning messages to console
    bool heterocyclic       = false;  // L|T designator can throw warnings

    // -- paths -- // 
    // int allows way more description in states

    unsigned int state_multi          = 0; // 0 - closed, 1 - open multi notation, 2 - expect size denotation
    unsigned int state_pseudo         = 0; 
    unsigned int state_bridge         = 0;
    unsigned int state_aromatics      = 0;

    bool implied_assignment_used       = 0; // allows a shorthand if wanted, but not mixing
   
    unsigned int expected_locants     = 0;

    unsigned int  evaluating_break      = 0;
    unsigned char ring_size_specifier   = '\0';
    unsigned char positional_locant     = '\0';

    std::string special;  

    std::vector<bool> aromaticity; 
    std::vector<std::pair<unsigned char, unsigned char>>  bond_increases; 

    std::vector<unsigned char> pseudo_locants;
    std::vector<unsigned int> pseudo_positions; 
    std::vector<unsigned char> bridge_locants;
    std::vector<unsigned char> multicyclic_locants;
    std::set<unsigned char> broken_locants;
    
    // broken locants start at A = 129 for extended ascii 
    // first is the standard 'X-' second is 'X-&', third is 'X--', fourth is 'X--&' and so on

    std::vector<std::pair<unsigned int, unsigned char>>  ring_components;
    std::vector<indexed_pair>                            indexed_bindings;  
    
  
    unsigned int i = 0; 
    const char *block_str = block.c_str();    
    unsigned int len = block.size();
    unsigned int aromatic_position = FindAromatics(block_str,len); // should be a copy so no effect on pointer

    unsigned char ch = *block_str++;


    while(ch){

      if(i >= aromatic_position)
        state_aromatics = 1;

      switch(ch){

        // specials

        case ' ':
          if(evaluating_break){
            broken_locants.insert(positional_locant);

            if(state_multi == 1 && expected_locants){
              multicyclic_locants.back() = positional_locant;
              state_multi = 2;
              expected_locants--;
            }
            else if(state_pseudo == 1 && expected_locants){
              pseudo_locants.back() = positional_locant;
              expected_locants--;
            }

            evaluating_break = 0;
          }
          if(expected_locants){
            fprintf(stderr,"Error: %d locants expected before space character\n",expected_locants);
            Fatal(i+start);
          }
          else if(state_multi == 1){
            state_multi = 2;
          }
          state_pseudo = 0;
          positional_locant = '\0'; // hard resets on spaces
          break;

        case '&':
          if (state_aromatics){
            aromaticity.push_back(1);
            break;
          }
          else if (state_multi == 3){
            ring_size_specifier += 23;
          }
          else if(positional_locant){
            // can only be an extension of a positional multiplier for a branch
            positional_locant += 23;
          }
          else{
            if(ch > 252){
              fprintf(stderr,"Error: creating molecule with atoms > 252, is this a reasonable for WLN?\n");
              Fatal(i+start);
            }

            ch = block[i-1] + 23; 

            if(positional_locant){
              positional_locant = ch;
              break;
            }
          }
          break;

        case '/':
          if(state_aromatics){
            fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
            Fatal(i+start);
          }

          if(!pseudo_positions.empty() && pseudo_positions.back() == ring_components.size() -1){
            for (unsigned int p=0;p<pseudo_positions.size();p++){
              pseudo_positions[p] += - 1; // back shift;
            }
          }
          pseudo_positions.push_back(ring_components.size() -1); // which ring has the pseudo?
          expected_locants = 2; 
          state_pseudo = 1;
          break; 
        

        // turn this into a look ahead type bias in order to significantly tidy this up
        case '-':{

          if(state_aromatics)
            break;

          // gives us a local working copy
          char local_arr [strlen(block_str)+1]; 
          memset(local_arr,'\0',strlen(block_str)+1);
          memcpy(local_arr,block_str,strlen(block_str));
          const char *local = local_arr;

          unsigned char local_ch = *(local)++; // moved it over
          unsigned int gap = 0; 
          bool found_next = false;
          
          while(local_ch != '\0'){
            if(local_ch == ' ')
              break;
            if(local_ch == '-'){
              // this calculates the gap to the next '-'
              found_next = true;
              break;
            }
            special.push_back(local_ch);
            gap++;
            local_ch = *local++;
          }

          // this will change on metallocenes defintions
          if( (state_multi || state_pseudo) && expected_locants){
           gap = 0;
          }

          // could ignore gap zeros as they will come back round again, therefore need a size check

          if(found_next){
            
            // pointer is moved by 1 at the bottom, so its positions -1 
            switch(gap){
              case 0:
                // we resolve only the first one
                evaluating_break = 1; // on a space, number or character, we push to broken_locants

                if(positional_locant){
                  if(positional_locant < 128){
                    positional_locant = create_relative_position(positional_locant); // i believe breaking modifier will then get removed
                    if(!positional_locant)
                      Fatal(i+start);
                  }
                  else{
                    // this means its already been moved, so we move the locant 23+23 across
                    if(positional_locant + 46 > 252){
                      fprintf(stderr,"Error: branching locants are exceeding the 252 space restriction on WLN notation, is this a reasonable molecule?\n");
                      Fatal(start+i);
                    }
                    positional_locant += 46;
                    // no need to move the global pointer here
                  }
                }
                else{
                  fprintf(stderr,"Error: trying to branch out character without starting point\n");
                  Fatal(start+i);  
                }

                break;
              case 1:
                if(!implied_assignment_used){
                  implied_assignment_used = true;
                  positional_locant = 'A';
                }
                // this can only be hypervalent element
                if(positional_locant){
                  WLNSymbol* new_locant = assign_locant(positional_locant,define_hypervalent_element(special[0]));  // elemental definition 
                  if(!new_locant)
                    Fatal(i+start);

                  string_positions[start+i + 1] = new_locant; // attaches directly

                  if(opt_debug)
                    fprintf(stderr,"  assigning hypervalent %c to position %c\n",special[0],positional_locant);

                  positional_locant++; // allow inline definition
                }
                else{
                  fprintf(stderr,"Error: trying to assign element without starting point\n");
                  Fatal(start+i);  
                }
                block_str+=2; 
                i+=2; 
                break;
              case 2:
                if(!implied_assignment_used){
                  implied_assignment_used = true;
                  positional_locant = 'A';
                }

                if(std::isdigit(special[0])){
                  for(unsigned char dig_check : special){
                    if(!std::isdigit(dig_check)){
                      fprintf(stderr,"Error: mixing numerical and alphabetical special defintions is not allowed\n");
                      Fatal(start+i);
                    }
                  }
                  if(positional_locant)
                    ring_components.push_back({std::stoi(special),positional_locant}); //big ring
                  else
                    ring_components.push_back({std::stoi(special),'A'});
                }
                else{
                  if(positional_locant){
                    WLNSymbol* new_locant = assign_locant(positional_locant,define_element(special));  // elemental definition
                    if(!new_locant)
                      Fatal(i+start);

                    string_positions[start+i + 1] = new_locant; // attaches directly to the starting letter

                    if(opt_debug)
                      fprintf(stderr,"  assigning element %s to position %c\n",special.c_str(),positional_locant);

                    positional_locant++; // allow inline definition
                  }
                  else{
                    fprintf(stderr,"Error: trying to assign element without starting point\n");
                    Fatal(start+i);  
                  }
                }
                block_str+=3; 
                i+=3;              
                break;
              default:
                fprintf(stderr,"Error: %d numerals incased in '-' brackets is unreasonable for WLN to create\n",gap);
                Fatal(start+i);
            }

          }
          else{
            // if there wasnt any other symbol, it must be a notation extender
            evaluating_break = 1; // on a space, number or character, we push to broken_locants

            if(positional_locant){
              if(positional_locant < 128){
                positional_locant = create_relative_position(positional_locant); // i believe breaking modifier will then get removed
                if(!positional_locant)
                  Fatal(i+start);
              }
              else{
                // this means its already been moved, so we move the locant 23+23 across
                if(positional_locant + 46 > 252){
                  fprintf(stderr,"Error: branching locants are exceeding the 252 space restriction on WLN notation, is this a reasonable molecule?\n");
                  Fatal(start+i);
                }
                positional_locant += 46;
                // no need to move the global pointer here
              }
            }
            else{
              fprintf(stderr,"Error: trying to branch out character without starting point\n");
              Fatal(start+i);  
            }

          }
          special.clear();
          break;
        }

        // numerals - easy access

        case '0':
          fprintf(stderr,"Error: Metallocene and Catenane compounds are valid within WLN notation, however\n"
                         "converting between common formats (smi & InChI) leads to undefined and undesirable\n"
                         "behaviour, see reconnected InChi's for a modern way of representing these compounds\n"
                         "as a line notation. For now, these will be unsupported alongside WLN 'uncertainties'\n");
          Fatal(i+start);

        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if(state_aromatics){
            fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
            Fatal(i+start);
          }

          if(evaluating_break){
            broken_locants.insert(positional_locant);

            if(state_multi == 1 && expected_locants){
              multicyclic_locants.back() = positional_locant;
              expected_locants--;
            }
            else if(state_pseudo == 1 && expected_locants){
              pseudo_locants.back() = positional_locant;
              expected_locants--;
            }

            evaluating_break = 0;
          }
          if (i > 1 && block[i-1] == ' '){
            state_multi   = 1; // enter multi state
            expected_locants = ch - '0';
          }  
          else{

            if(positional_locant)
              ring_components.push_back({ch-'0',positional_locant});
            else
              ring_components.push_back({ch-'0','A'});

            positional_locant = '\0';
          }
          break;

        case 'A':
        case 'B':
        case 'C':
        case 'D':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'I':
        case 'K':
        case 'M':
        case 'N':
        case 'O':
        case 'P':
        case 'Q':
        case 'R':
        case 'S':
        case 'U':
        case 'V':
        case 'W':
        case 'X':
        case 'Y':
        case 'Z':

          if(state_aromatics){
            fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
            Fatal(i+start);
          }

          if(evaluating_break){
            broken_locants.insert(positional_locant);

            if(state_multi == 1 && expected_locants)
              multicyclic_locants.back() = positional_locant;
            else if(state_pseudo == 1 && expected_locants)
              pseudo_locants.back() = positional_locant;
            
            evaluating_break = 0;
          }

          if(expected_locants){

            positional_locant = ch; // use for look back
            expected_locants--;

            if(state_multi)
              multicyclic_locants.push_back(ch);
            else if (state_pseudo)
              pseudo_locants.push_back(ch);
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }
            break;
          }
          else if(state_multi == 2){
            ring_size_specifier = ch;
            state_multi = 3;
          }
          else if (positional_locant){
           
            if (opt_debug)
              fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

            WLNSymbol *new_locant = 0; 

            switch(ch){
              case 'S':
              case 'P':
                if(!heterocyclic)
                  warned = true;

                if(locants[positional_locant])
                  positional_locant++; 

                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edge_and_type(5,RING);
                break;

              case 'Y':
                if(locants[positional_locant])
                  positional_locant++; 
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edge_and_type(3,RING);
                break;
              case 'N':
                if(!heterocyclic)
                  warned = true;

                if(locants[positional_locant])
                  positional_locant++; 

                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edge_and_type(3,RING);
                break;

              case 'V':
                if(locants[positional_locant])
                  positional_locant++; 
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edge_and_type(2,RING);
                break;

              case 'M':
              case 'O':
                if(!heterocyclic)
                  warned = true;

                if(locants[positional_locant])
                  positional_locant++; 

                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edge_and_type(2,RING);
                break;

              case 'X':
                if(locants[positional_locant])
                  positional_locant++; 

                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edge_and_type(4,RING);
                break;
              case 'K':
                if(!heterocyclic)
                  warned = true;

                if(locants[positional_locant])
                  positional_locant++; 

                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edge_and_type(4,RING);
                break;

              case 'U':
                // no need to put this in implied, it has to be specified
                if(i < len - 3 && block[i+1] == '-' && block[i+2] == ' '){
                  bond_increases.push_back({positional_locant,block[i+3]});
                  block_str += 3;
                  i += 3;
                }
                else
                  bond_increases.push_back({positional_locant,positional_locant+1});
                break;

              case 'W':
                switch(locants[positional_locant]->ch){
                  case 'K':
                    locants[positional_locant]->allowed_edges++;
                    break;
                  default:
                    break;
                }
                if(!add_diazo(locants[positional_locant]))
                  Fatal(i+start);
                break;

              default:
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
            }
            string_positions[start+i] = new_locant;
          }
          else{

            if(i>0 && block[i-1] == ' '){
              positional_locant = ch;
            }
            else{
              implied_assignment_used = true;
              positional_locant = 'A';

              if (opt_debug)
                fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

              WLNSymbol *new_locant = 0; 

              switch(ch){
                case 'S':
                case 'P':
                  if(!heterocyclic)
                    warned = true;

                  if(locants[positional_locant])
                    positional_locant++;     

                  new_locant = assign_locant(positional_locant,ch);
                  new_locant->set_edge_and_type(5,RING);
                  break;

                case 'Y':
                  if(locants[positional_locant])
                    positional_locant++; 

                  new_locant = assign_locant(positional_locant,ch);
                  new_locant->set_edge_and_type(3,RING);
                  break;
                case 'N':
                  if(!heterocyclic)
                    warned = true;
                  
                  if(locants[positional_locant])
                    positional_locant++; 

                  new_locant = assign_locant(positional_locant,ch);
                  new_locant->set_edge_and_type(3,RING);
                  break;

                case 'V':
                  if(locants[positional_locant])
                      positional_locant++; 
                  new_locant = assign_locant(positional_locant,ch);
                  new_locant->set_edge_and_type(2,RING);
                  break;

                case 'M':
                case 'O':
                  if(!heterocyclic)
                    warned = true;
                  
                  if(locants[positional_locant])
                    positional_locant++; 

                  new_locant = assign_locant(positional_locant,ch);
                  new_locant->set_edge_and_type(2,RING);
                  break;

                case 'X':
                  if(locants[positional_locant])
                    positional_locant++; 

                  new_locant = assign_locant(positional_locant,ch);
                  new_locant->set_edge_and_type(4,RING);
                  break;

                case 'K':
                  if(!heterocyclic)
                    warned = true;

                  if(locants[positional_locant])
                    positional_locant++; 

                  new_locant = assign_locant(positional_locant,ch);
                  new_locant->set_edge_and_type(4,RING);
                  break;

                case 'U':

                  if(i < len )
                  
                  bond_increases.push_back({positional_locant,positional_locant+1});
                break;

                default:
                  fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                  Fatal(start+i);
              }

              string_positions[start+i] = new_locant;
            }
          }

          break;


        case 'L':
          if(state_aromatics){
            fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
            Fatal(i+start);
          }

          if(evaluating_break){
            broken_locants.insert(positional_locant);

            if(state_multi == 1 && expected_locants)
              multicyclic_locants.back() = positional_locant;
            else if(state_pseudo == 1 && expected_locants)
              pseudo_locants.back() = positional_locant;

            evaluating_break = 0;
          }

          if(i==0){
            heterocyclic = false; 
            break;
          }
          if(expected_locants){

            positional_locant = ch; // use for look back
            expected_locants--;

            if(state_multi)
              multicyclic_locants.push_back(ch);
            else if (state_pseudo)
              pseudo_locants.push_back(ch);
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }

            break;
          }
          else if(state_multi == 2){
            ring_size_specifier = ch;
            state_multi = 3;
          }
          else{
            if(i>0 && block[i-1] == ' '){
              positional_locant = ch;
            }
            else{
              fprintf(stderr,"Error: symbol '%c' is in an unhandled state, please raise issue if this notation is 100%% correct\n",ch);
              Fatal(i+start);
            }
          }
        
          break;


        case 'T':
          if(state_aromatics){
            aromaticity.push_back(0);
            break;
          }
       
          if(evaluating_break){
            broken_locants.insert(positional_locant);

            if(state_multi == 1 && expected_locants)
              multicyclic_locants.back() = positional_locant;
            else if(state_pseudo == 1 && expected_locants)
              pseudo_locants.back() = positional_locant;

            evaluating_break = 0;
          }

          if(i==0){
            heterocyclic = true; 
            break;
          }

          if(expected_locants){

            positional_locant = ch; // use for look back
            expected_locants--;

            if(state_multi)
              multicyclic_locants.push_back(ch);
            else if (state_pseudo)
              pseudo_locants.push_back(ch);
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }
            break;
          }
          else if(state_multi == 2){
            ring_size_specifier = ch;
            state_multi = 3;
          }
          else{
            if(i>0 && block[i-1] == ' '){
              positional_locant = ch;
            }
            else{
              fprintf(stderr,"Error: symbol '%c' is in an unhandled state, please raise issue if this notation is 100%% correct\n",ch);
              Fatal(i+start);
            }
          }
            

          break;

        
        // CLOSE

        case 'J':
          if(state_aromatics)
            state_aromatics = 0;
          

          if(evaluating_break){
            broken_locants.insert(positional_locant);

            if(state_multi == 1 && expected_locants)
              multicyclic_locants.back() = positional_locant;
            else if(state_pseudo == 1 && expected_locants)
              pseudo_locants.back() = positional_locant;

            evaluating_break = 0;
          }

          if (i == block.size()-1){

            if(ring_components.empty()){
              fprintf(stderr,"Error: error in reading ring components, check numerals in ring notation\n");
              Fatal(start+i);
            }

            if (pseudo_locants.size() > 0)
              ring_type = PSDBRIDGED;

            if (multicyclic_locants.size() > 0 && ring_type < PSDBRIDGED)
              ring_type = PERI;

            if (aromaticity.size() == 1 && aromaticity[0] == false){
              while(aromaticity.size() < ring_components.size())
                aromaticity.push_back(false);
            }
            else if (aromaticity.empty()){
              while(aromaticity.size() < ring_components.size())
                aromaticity.push_back(true);
            }

            // perform the aromatic denotion check
            if (ring_components.size() != aromaticity.size()){
              fprintf(stderr,"Error: mismatch between number of rings and aromatic assignments - %ld vs expected %ld\n",aromaticity.size(),ring_components.size());
              Fatal(i+start);
            }

            // create the bindings needed for pseudo bridges
            for (unsigned int i=0; i< pseudo_positions.size();i++){
              indexed_pair pseudo; 
              pseudo.set(pseudo_locants[i+i],pseudo_locants[i+i+1],pseudo_positions[i]);
              indexed_bindings.push_back(pseudo);
            }

            break;
          }
          if(expected_locants){

            positional_locant = ch; // use for look back
            expected_locants--;

            if(state_multi)
              multicyclic_locants.push_back(ch);
            else if (state_pseudo)
              pseudo_locants.push_back(ch);
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }
            break;
          }
          else if(state_multi == 2){
            ring_size_specifier = ch;
            state_multi = 3;
          }
          else{
            if(i>0 && block[i-1] == ' '){
              positional_locant = ch;
            }
            else{
              fprintf(stderr,"Error: symbol '%c' is in an unhandled state, please raise issue if this notation is 100%% correct\n",ch);
              Fatal(i+start);
            }
          }
          
          break;

        default:
          // these can only be expanded chars
          fprintf(stderr,"WARNING: SWITCH UNCLOSED\n");
      }
      
      i++;
      ch = *(block_str++);
    }

    if(warned)
      fprintf(stderr,"Warning: heterocyclic ring notation required for inter atom assignment, change starting 'L' to 'T'\n");
    

    // debug here
    if (opt_debug){
      
      fprintf(stderr,"  ring type: %s\n",ring_strings[ring_type]);

      fprintf(stderr,"  ring components: ");
      for (std::pair<unsigned int, unsigned char> comp : ring_components){
        
        if(comp.second > 'Z')
          fprintf(stderr,"%d(%d) ",comp.first,comp.second);
        else
          fprintf(stderr,"%d(%c) ",comp.first,comp.second);
      } 
        
      fprintf(stderr,"\n");

      fprintf(stderr,"  aromaticity: ");
      for (bool aromatic : aromaticity)
        fprintf(stderr,"%d ",aromatic);
      fprintf(stderr,"\n");

      fprintf(stderr,"  multicyclic points: ");
      for (unsigned char loc : multicyclic_locants){
        if(loc > 'Z')
          fprintf(stderr,"%d ",loc);
        else
          fprintf(stderr,"%c ",loc);
      }
      fprintf(stderr,"\n");

      fprintf(stderr,"  broken path points: ");
      for (unsigned char loc : broken_locants){
        fprintf(stderr,"%d ",loc);
      }
      fprintf(stderr,"\n");

      fprintf(stderr,"  bridge points: ");
      for (unsigned char loc : bridge_locants){
        fprintf(stderr,"%c ",loc == ' ' ? '_':loc);
      }
      fprintf(stderr,"\n");

      fprintf(stderr,"  pseudo bridge points: ");
      for (unsigned int i=0; i< pseudo_positions.size();i++){
        fprintf(stderr,"(%d)[%c <-- %c] ",pseudo_positions[i],pseudo_locants[i+i],pseudo_locants[i+i+1]);
      }
      fprintf(stderr,"\n");
  
      fprintf(stderr,"  size denotion: %d\n",ring_size_specifier ? locant_to_int(ring_size_specifier) : 0);
      fprintf(stderr,"  heterocyclic: %s\n", heterocyclic ? "yes":"no");
    }
    
    bool state = true;
    switch(ring_type){
      case POLY:
        state = CreatePolyCyclic(ring_components,aromaticity);
        break;
      case PERI:
      case PSDBRIDGED:
        state = CreateMultiCyclic(ring_components,aromaticity,
                                  multicyclic_locants,indexed_bindings,
                                  broken_locants,
                                  ring_size_specifier);
        break;
      case BRIDGED:
        break;
    }

    if (!state)
      Fatal(start+i);

    for (std::pair<unsigned char, unsigned char> bond_pair : bond_increases){
      WLNEdge *increase_edge = search_edge(locants[bond_pair.second],locants[bond_pair.first]);
      increase_edge = unsaturate_edge(increase_edge,1);
      if(!increase_edge)
        Fatal(start+i);
    }
  }

};

WLNRing *AllocateWLNRing()
{
  ring_count++;
  if(ring_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln rings - is this reasonable?\n");
    return 0;
  }

  WLNRing *wln_ring = new WLNRing;
  RINGS[ring_count] = wln_ring;
  return wln_ring;
}


struct WLNGraph
{
  
  WLNSymbol *root;

  WLNGraph() : root{(WLNSymbol *)0} {};
  ~WLNGraph()
  {
    for (unsigned int i = 0; i < 1024;i++){
      if(SYMBOLS[i])
        delete SYMBOLS[i];
      if(EDGES[i])
        delete EDGES[i];
      if(RINGS[i])
        delete RINGS[i];
    }
  };


  bool expand_carbon_chain(WLNSymbol *head,unsigned int size){

    if (size > REASONABLE)
      fprintf(stderr,"Warning: making carbon chain over 1024 long, reasonable molecule?\n");

    head->ch = 'C';
    head->set_edge_and_type(4);

    if(size == 1)
      return true;
    
    // leave the original node where it is, and expand out
    WLNEdge *tmp = 0;
    WLNSymbol *bonded = 0;
    unsigned int bonded_order;

    // if the chain has any children
    if(head->bonds){
      // hold it
      tmp = head->bonds;

      bonded = tmp->child;
      bonded_order = tmp->order;
      if(!remove_edge(head,tmp))
        return false;
    }

    WLNEdge *edge = 0;
    WLNSymbol *prev = head;
    for(unsigned int i=0;i<size-1;i++){
      WLNSymbol* carbon = AllocateWLNSymbol('C');
      carbon->set_edge_and_type(4); // allows hydrogen resolve
      edge = AllocateWLNEdge(carbon,prev);
      prev = carbon;
    } 

    if(bonded){
      edge = AllocateWLNEdge(bonded,prev);
      if(!edge)
        return false;

      if(bonded_order > 2)
        edge = unsaturate_edge(edge,bonded_order-1);
    }

    return true;
  }



  /* must be performed before sending to obabel graph*/
  bool ExpandWLNGraph(){

 
    unsigned int stop = symbol_count;
    for (unsigned int i=1;i<=stop;i++){
      WLNSymbol *sym = SYMBOLS[i];

      if(sym->type == LOCANT){
        // floating locants get a methyl
        if(!sym->bonds){
          if(!add_methyl(sym))
            return false;
        }
        else 
          continue;
      }
        
      switch(sym->ch){

        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if (!sym->special.empty()){
            if(!expand_carbon_chain(sym,std::stoi(sym->special))){
              fprintf(stderr,"Error: error in expanding out numeral to carbon chain\n");
              return false;
            }
          }
          else{
            if(!expand_carbon_chain(sym,sym->ch - '0')){
              fprintf(stderr,"Error: error in expanding out numeral to carbon chain\n");
              return false;
            }
          }
          break;
        
        case 'Y':
        case 'X':
        case 'K':
          resolve_methyls(sym);
          break;

        case 'V':{
          sym->ch = 'C';
          sym->set_edge_and_type(4);
          WLNSymbol *oxygen = AllocateWLNSymbol('O');
          oxygen->set_edge_and_type(2);
          WLNEdge *e = AllocateWLNEdge(oxygen,sym);
          e = unsaturate_edge(e,1);
          if(!e)
            return false;
          break;
        }
        
        case 'W':
          sym->ch = 'C';
          sym->set_edge_and_type(4);
          if(!add_diazo(sym))
            return false;
          break;


        default:
          break; // ignore
      }
    }
    
    return true; 
  }


  // this one pops based on bond numbers
  WLNSymbol *return_open_branch(std::stack<WLNSymbol *> &branch_stack){
     // only used for characters that can 'act' like a '&' symbol
    
    WLNSymbol *top = 0;
    if (branch_stack.empty())
      return top;

    while(!branch_stack.empty()){
      top = branch_stack.top();
      if(top->num_edges == top->allowed_edges)
        branch_stack.pop();
      else
        return top;
    }

    return top;
  }



  /*  bonds a locant to the ring, returns the edge */
  WLNEdge* assign_locant_to_ring(WLNSymbol *curr,unsigned int bond_modifier,std::stack<WLNRing *> &ring_stack)
  {
    if(!curr)
      return 0;

    // locant should be allocated already here
    WLNRing *s_ring = 0;
    WLNEdge *edge   = 0;

    if (ring_stack.empty()){
      fprintf(stderr, "Error: no rings to assign locants to\n");
      return 0;
    }
    else
      s_ring = ring_stack.top();

    if (s_ring->locants[curr->ch]){
      edge = AllocateWLNEdge(curr,s_ring->locants[curr->ch]);
      if(bond_modifier)
        edge = unsaturate_edge(edge,bond_modifier);
    } 
    else {
      fprintf(stderr, "Error: assigning locant outside of ring\n");
      return 0;
    }

    return edge; 
  }

  /* backwards search for tentative ionic rule procedures */
  unsigned int search_ionic(const char *wln_ptr, unsigned int len,
                           std::vector<std::pair<unsigned int, int>> &charges)
  {
    unsigned int first_instance = 0;
  
    for (unsigned int i=0;i<len;i++){

      // these are required in blocks of 5
      if(wln_ptr[i] == ' ' && wln_ptr[i+1] == '&')
      {
        
        std::string position_1;
        std::string position_2;

        unsigned int local_search = i+2;

        if(std::isdigit(wln_ptr[i+2])){

          while(std::isdigit(wln_ptr[local_search])){
            position_1.push_back(wln_ptr[local_search]);
            local_search++;

            if(local_search > len)
              return first_instance;
          }
        }
        else 
          continue;

        // local search should now be pointing at the '\'
        if(wln_ptr[local_search] == '/')
          local_search++;
        else
          continue;

        if(std::isdigit(wln_ptr[local_search])){
          while(std::isdigit(wln_ptr[local_search])){
            position_2.push_back(wln_ptr[local_search]);
            local_search++;

            if(local_search > len)
              return first_instance;
          }
        }
        else 
          continue;

        
        if(std::stoi(position_1) != 0)
          charges.push_back({std::stoi(position_1),1});
        
        if(std::stoi(position_2) != 0)
          charges.push_back({std::stoi(position_2),-1});

        if(!first_instance)
          first_instance = i;
      }
    }

    return first_instance;
  }

  /* uses the global position map */
  bool AssignCharges(std::vector<std::pair<unsigned int, int>> &charges){
    if(charges.empty())
      return true;

    for (std::pair<unsigned int, int> pos_charge : charges){
      WLNSymbol *assignment = string_positions[pos_charge.first - 1]; // reindex as wln 1 is string 0
      if(!assignment){
        fprintf(stderr,"Error: trying to assign ionic charge to unavaliable element, check that character %d is avaliable for assignment\n",pos_charge.first);
        return false;
      }
      else if(assignment->type == LOCANT){
        fprintf(stderr,"Error: trying to assign charge to a locant character\n");
        return false;
      }
      else{
        charge_additions[assignment] += pos_charge.second;

        if(opt_debug){
          fprintf(stderr, "  character at position [%d] has the following charge addition - %d\n",pos_charge.first,pos_charge.second);
        }
      }
    }
    return true;
  }


  // performs standard sequence multiplers from a bind point
  WLNSymbol *sequence_multiplier(WLNSymbol *prev,bool mirror,unsigned int n, std::string &sequence,unsigned int bond_ticks=0){
    
    // return the last symbol created
    return 0;
  }


  
  /* returns the head of the graph, parse all normal notation */
  bool ParseWLNString(std::string wln_string) 
  {

    if (opt_debug)
      fprintf(stderr, "Parsing WLN notation: %s\n",wln_string.c_str());

    std::stack<WLNRing *>   ring_stack;   // access through symbol
    std::stack<WLNSymbol *> branch_stack; // between locants, clean branch stack
    std::stack<WLNSymbol *> linker_stack; // used for branching ring systems
    
    std::vector<std::pair<unsigned int, int>> ionic_charges;
    
    WLNSymbol *curr   = 0;
    WLNSymbol *prev   = 0;
    WLNEdge   *edge   = 0;
    WLNRing   *ring   = 0;

    bool pending_locant           = false;
    bool pending_J_closure        = false;
    bool pending_inline_ring      = false;
    bool pending_spiro            = false;
    bool pending_diazo            = false;
    bool pending_linker           = false;
    unsigned int pending_unsaturate = 0;    // 'U' style bonding
   
    std::string special;  // special elemental definitions
    
    // allows consumption of notation after block parses
    unsigned int block_start = 0;
    unsigned int block_end = 0;
 
    
    unsigned int len = wln_string.length();
    const char * wln_ptr = wln_string.c_str();

    unsigned int zero_position = search_ionic(wln_ptr,len,ionic_charges);
    
    unsigned int i=0;
    unsigned char ch = *wln_ptr;
    
    while(ch)
    {  
      
      // dont read any ionic notation
      if(zero_position && zero_position == i)
        break;

      switch (ch)
      {

      case '0': // cannot be lone, must be an addition to another num
        if (pending_J_closure)
          break;
        
        if(prev && std::isdigit(prev->ch)){
          prev->special.push_back(ch);
        }
        else
          Fatal(i);
        break;

      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        if (pending_J_closure)
          break;
        else if(pending_locant){  // handle all multiplier contractions
          
          // block to hunt the multiplier maximum --> lookahead
          // look for spaces followed by ints. 
 
          std::string int_sequence;
          int_sequence.push_back(ch);

          while(i < len - 2){
            if(wln_string[i+1] == ' ' && std::isdigit(wln_string[i+2])){
              int_sequence.push_back(wln_string[i+2]);
              i+=2;
              wln_ptr += 2;
            }
            else
              break;
          }

          // pointer is moved to the last number, multiplier value is calculated
          unsigned int multiplier_value = std::stoi(int_sequence);

          fprintf(stderr,"Error: multipliers are not currently supported\n");
          Fatal(i);
          
          pending_linker = true;
          pending_locant = false;
        }
        else{

          if(pending_diazo){
              // do a hydroxy transform here
              
            curr = prev; // might be overkill 
            curr->set_edge_and_type(4);

            if(!add_diazo(curr))
              Fatal(i-1);
              
            curr->ch = ch;
            pending_diazo = false;
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(4);
            
            if(prev)
              edge = AllocateWLNEdge(curr,prev);
          }
          
          // moves naturally, so end on the last number
    
          curr->special.push_back(ch);

          while(*(wln_ptr+1)){
            if(!std::isdigit(*(wln_ptr+1)))
              break;

            curr->special.push_back(*wln_ptr);
            wln_ptr++;
            i++;
          }

          if(pending_linker){
          }

          pending_unsaturate = 0;
          prev = curr;
          break;
        }
        break;
     

      case 'Y':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);
          
          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            curr = prev; 
            curr->set_edge_and_type(3);

            if(!add_diazo(curr))
              Fatal(i-1);
              
            curr->ch = ch;
            pending_diazo = false;
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(3);
            if(prev)
              edge = AllocateWLNEdge(curr,prev);
          }

          branch_stack.push(curr);

          string_positions[i] = curr;
          pending_unsaturate = 0;
          prev = curr;
        }
        break;

      case 'X':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            curr = prev; 
            curr->set_edge_and_type(4);

            if(!add_diazo(curr))
              Fatal(i-1);
              
            curr->ch = ch;
            pending_diazo = false;
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(4);
            if(prev){
              edge = AllocateWLNEdge(curr,prev);
              if(!edge)
                Fatal(i);

              if(pending_unsaturate){
                edge = unsaturate_edge(edge,pending_unsaturate);
                pending_unsaturate = 0;
              }
            }
          }

          branch_stack.push(curr);
          string_positions[i] = curr;
          prev = curr;
        }
        break;

        // oxygens

      case 'O':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            fprintf(stderr,"Error: diazo assignment to an oxygen is a disallowed bond type\n");
            Fatal(i);
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2);

          branch_stack.push(curr);

          if(prev){
            edge = AllocateWLNEdge(curr,prev);
            if(!edge)
              Fatal(i);

            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
          }

          string_positions[i] = curr;
          prev = curr;
        }
        break;

      case 'Q':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {

          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if(pending_diazo){
            fprintf(stderr,"Error: diazo assignment to an oxygen is a disallowed bond type\n");
            Fatal(i);
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(1);

          if(prev){
            edge = AllocateWLNEdge(curr,prev);
            
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }

            if(!edge)
              Fatal(i);
          }

          string_positions[i] = curr;
          pending_unsaturate = 0;
          prev = return_open_branch(branch_stack);
          if(!prev)
            prev = curr;
        }
        break;

      case 'V':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);
      
          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            fprintf(stderr,"Error: diazo assignment to an carbonyl is a disallowed bond type\n");
            Fatal(i);
          }
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,STANDARD);
          if(prev){
            edge = AllocateWLNEdge(curr,prev);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            
            if(!edge)
              Fatal(i);
          }

          string_positions[i] = curr;
          prev = curr;
        }
        break;

      case 'W':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            fprintf(stderr,"Error: double diazo assignment is a disallowed bond type\n");
            Fatal(i);
          }

          if(prev){
            prev->allowed_edges++;
            if(!add_diazo(prev))
              Fatal(i);
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(2);
            if(prev){
              edge = AllocateWLNEdge(curr,prev);
              if(pending_unsaturate){
                edge = unsaturate_edge(edge,pending_unsaturate);
                pending_unsaturate = 0;
              }
              if(!edge)
                Fatal(i);
            }
            pending_diazo = true;
          }

          string_positions[i] = curr;
          pending_unsaturate = 0;
          prev = curr;
        }
        break;

        // nitrogens

      case 'N':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          
          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
        
            curr = prev;
            curr->set_edge_and_type(4); // 'x-NW' is allowed 

            if(!add_diazo(curr))
              Fatal(i-1);
            
            curr->ch = ch;
            pending_diazo = false;
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(3);

            if(prev){
              edge = AllocateWLNEdge(curr,prev);
              if(pending_unsaturate){
                edge = unsaturate_edge(edge,pending_unsaturate);
                pending_unsaturate = 0;
              }
              if(!edge)
                Fatal(i);
              }
            }
         
          branch_stack.push(curr);

          string_positions[i] = curr;
          pending_unsaturate = 0;
          prev = curr;
        }
        break;

      case 'M':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            fprintf(stderr,"Error: diazo assignment to NH is a disallowed bond type\n");
            Fatal(i);
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2);

          if(prev){
            edge = AllocateWLNEdge(curr,prev);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }

          string_positions[i] = curr;
          pending_unsaturate = 0;
          prev = curr;
        }
        break;

      case 'K':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            curr = prev; 
            curr->set_edge_and_type(5);

            if(!add_diazo(curr))
              Fatal(i-1);
              
            curr->ch = ch;
            pending_diazo = false;
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(4);
            if(prev){
              edge = AllocateWLNEdge(curr,prev);
              if(pending_unsaturate){
                edge = unsaturate_edge(edge,pending_unsaturate);
                pending_unsaturate = 0;
              }
              if(!edge)
                Fatal(i);
            }
          }

          branch_stack.push(curr);
          string_positions[i] = curr;
          prev = curr;
        }
        break;

      case 'Z':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            fprintf(stderr,"Error: diazo assignment to NH2 is a disallowed bond type\n");
            Fatal(i);
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(1);

          if(prev){
            edge = AllocateWLNEdge(curr,prev);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            
            if(!edge)
              Fatal(i);
          }

          string_positions[i] = curr;
          pending_unsaturate = 0;
          prev = return_open_branch(branch_stack);
          if(!prev)
            prev = curr;
        }
        break;

        // halogens - need to add rules for semi allowed hyper valence in ionions

      case 'E':
      case 'G':
      case 'F':
      case 'I':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            fprintf(stderr,"Error: diazo assignment to a non expanded valence halogen is a disallowed bond type\n");
            Fatal(i);
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(1);

          if(prev){
            edge = AllocateWLNEdge(curr,prev);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            
            if(!edge)
              Fatal(i);
          }

          string_positions[i] = curr;
          pending_unsaturate = 0;
          prev = return_open_branch(branch_stack);
          if(!prev)
            prev = curr;
        }
        break;

        // inorganics

      case 'B':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring && prev)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            curr = prev; 
            curr->set_edge_and_type(3);

            if(!add_diazo(curr))
              Fatal(i-1);
              
            curr->ch = ch;
            pending_diazo = false;
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(3);
            if(prev){
              edge = AllocateWLNEdge(curr,prev);
              if(pending_unsaturate){
                edge = unsaturate_edge(edge,pending_unsaturate);
                pending_unsaturate = 0;
              }
            
              if(!edge)
                Fatal(i);
            }         
          }

          branch_stack.push(curr);
          string_positions[i] = curr;
          prev = curr;
        }
        break;

      case 'P':
      case 'S':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if(pending_diazo){
            curr = prev; 
            curr->set_edge_and_type(6); // might be overkill 

            if(!add_diazo(curr))
              Fatal(i-1);
            
            curr->ch = ch;
            pending_diazo = false;
          }
          else{
            curr = AllocateWLNSymbol(ch);
            curr->set_edge_and_type(6);

            if(prev){
              edge = AllocateWLNEdge(curr,prev);
              if(pending_unsaturate){
                edge = unsaturate_edge(edge,pending_unsaturate);
                pending_unsaturate = 0;
              }
          
              if(!edge)
                Fatal(i);
            }
          }

          branch_stack.push(curr);
          string_positions[i] = curr;
          prev = curr;
        }
        break;

        // locants only?

      case 'A':
      case 'C':
      case 'D':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);
          
            

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else{
          fprintf(stderr,"Error: locant only symbol used in atomic definition\n");
          Fatal(i);
        }
        break;
          
          
      // hydrogens explicit

      case 'H':
        if (pending_J_closure)
          break;
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);
 
          prev = curr;
          pending_locant = false;
        }
        else{
          // explicit hydrogens
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(1);

          // will add with more examples
          switch(prev->ch){
            case 'Z':
              charge_additions[prev]++;
              break;
            default:
              break;
          }

          if(prev){
            edge = AllocateWLNEdge(curr,prev);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }

          string_positions[i] = curr;
          prev = return_open_branch(branch_stack);
          if(!prev)
            prev = curr;
        }
        break;

        // ring notation

      case 'J':
        if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(curr,prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }

        else if (pending_J_closure && ( (i<len-1 && wln[i+1] == ' ') || i == len -1))
        {
          block_end = i;
          
          ring = AllocateWLNRing();
          std::string r_notation = get_notation(block_start,block_end);

          ring->FormWLNRing(r_notation,block_start);
          ring_stack.push(ring);

          block_start = 0;
          block_end = 0;

          if (pending_spiro)
          {
            prev->type = LOCANT; // spiros are normal rings with dual linker notation
            prev->previous->type = LOCANT;
            pending_spiro = false;
          }

          // does the incoming locant check
          if (prev)
          {
            if (ring->locants[prev->ch]){
              if(prev){
                edge = AllocateWLNEdge(ring->locants[prev->ch],prev);
                if(pending_unsaturate){
                  edge = unsaturate_edge(edge,pending_unsaturate);
                  pending_unsaturate = 0;
                }
                if(!edge)
                  Fatal(i);
              }
            }   
            else
            {
              fprintf(stderr, "Error: attaching inline ring with out of bounds locant assignment\n");
              Fatal(i);
            }
          }

          pending_J_closure = false;
        }
        
        break;

      case 'L':
      case 'T':
        if (pending_J_closure)
          break;
        
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(ring->locants[prev->ch],prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if (i == 0)
            pending_inline_ring = true;
          
        
          if (!pending_inline_ring)
          {
            fprintf(stderr, "Error: ring notation started without '-' denotion\n");
            Fatal(i);
          }
          else
            pending_inline_ring = false;

          block_start = i;
          pending_J_closure = true;
        }
        break;

      case 'R':
        if (pending_J_closure)
          break;
        
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(ring->locants[prev->ch],prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          
          prev = curr;
          pending_locant = false;
        }
        else
        {
          ring = AllocateWLNRing();

          std::string r_notation = "L6J";
          ring->FormWLNRing(r_notation,i);
          ring_stack.push(ring);

          if( !(wln_ptr + 1))
            ring->locants['A']->num_edges++; // will place a minus charge on the centre carbon

          curr = ring->locants['A'];
          if(prev){
            edge = AllocateWLNEdge(ring->locants[prev->ch],prev);
            
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }

          string_positions[i] = curr;
          prev = curr;;
        }
        break;

        // bonding

      case 'U':
        if (pending_J_closure)
          break;
        
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_edge_and_type(2,LOCANT); // locants always have two edges

          if (pending_inline_ring)
            edge = AllocateWLNEdge(ring->locants[prev->ch],prev);
          else
            edge = assign_locant_to_ring(curr, pending_unsaturate,ring_stack);

          if(!edge)
            Fatal(i);

          prev = curr;
          pending_locant = false;
        }
        else if (pending_diazo){
          fprintf(stderr,"Error: diazo assignment followed by a bond increase is a disallowed bond type\n");
          Fatal(i);
        }
        else
          pending_unsaturate++;
        
        break;

        // specials

      case ' ':
        if (pending_J_closure)
          break;
        
        else if (pending_diazo){
          fprintf(stderr,"Error: diazo assignment followed by a space seperator is a disallowed bond type\n");
          Fatal(i);
        }
        
        if (pending_inline_ring)
        {
          // send the linker into its own stack
          if (!branch_stack.empty() && (branch_stack.top()->num_edges < branch_stack.top()->allowed_edges))
            linker_stack.push(branch_stack.top());
        }
 
        // clear the branch stack betweek locants and ions
        while (!branch_stack.empty())
          branch_stack.pop();
        
        pending_locant = true;
        break;



      case '&':
        if (pending_diazo){
          fprintf(stderr,"Error: diazo assignment followed by a branch terminator is a disallowed bond type\n");
          Fatal(i);
        }

        if (pending_J_closure)
          break;
        
        if (pending_inline_ring)
        {
          // spiro notation open
          pending_spiro = true;
        }
        else if (pending_locant)
        {
          // ionic species or spiro, reset the linkings
          prev = 0;
          pending_locant = false;
        }
        else if(curr->type == LOCANT)
          curr->ch += 23;
        else if(i < len - 1 && wln_string[i+1] == ' '){

          ring = 0;
          if(!ring_stack.empty())
            ring_stack.pop();
          else{
            fprintf(stderr,"Error: popping too many rings, check '&' count\n");
            Fatal(i);
          }
          
          if(!ring_stack.empty())
            ring = ring_stack.top();
          else{
            fprintf(stderr,"Error: popping too many rings, check '&' count\n");
            Fatal(i);
          }
          break;
        }
        else
        {
          WLNSymbol *top = 0;
          if(!branch_stack.empty())
            top = branch_stack.top();
          
          if(!top){
            fprintf(stderr,"Error: '&' punctuation outside of branching characters is disallowed notation\n");
            Fatal(i);
          }

          // this means a <Y|X|..>'&' so handle methyl
          if(prev && prev == top){
            
            switch(prev->ch){
              // methyl contractions
              case 'X':
              case 'Y':
              case 'K':
                if(prev->num_edges < prev->allowed_edges){
                  if(!add_methyl(prev))
                    Fatal(i);

                  prev = return_open_branch(branch_stack);
                }
                else{ 
                  // we pop it as all contractions / bonds are made
                  prev = return_open_branch(branch_stack);
                }
                break;

              // no contractions possible, we pop the stack
              default:
                prev = return_open_branch(branch_stack);
                break;
            }

          }
          else{
            // means a closure is done, we return to the first avaliable symbol on the branch stack
            prev = return_open_branch(branch_stack);
          }

        
        }
        break;


      case '-':
        if (pending_J_closure)
          break;

        else if (pending_inline_ring)
        {
          fprintf(stderr, "Error: only one pending ring can be active, check closures\n");
          Fatal(i);
        }
        else{

          // look ahead and consume the special

          char local_arr [strlen(wln_ptr)+1]; 
          memset(local_arr,'\0',strlen(wln_ptr)+1);
          memcpy(local_arr,wln_ptr,strlen(wln_ptr));
          const char *local = local_arr;
          
          unsigned char local_ch = *(++local); // moved it over
          unsigned int gap = 0; 
          bool found_next = false;

          while(local_ch != '\0'){
            if(local_ch == ' ')
              break;
            if(local_ch == '-'){
              // this calculates the gap to the next '-'
              found_next = true;
              break;
            }
            special.push_back(local_ch);
            gap++;
            local_ch = *(++local);
          }
          
          if(!found_next)
            pending_inline_ring = true;
          else{
            if(gap == 1){
              curr = define_hypervalent_element(special[0]);
              if(!curr)
                Fatal(i);
            }
            else if(gap == 2){
              curr = define_element(special);
              if(!curr)
                Fatal(i); 
            }
            else{
              fprintf(stderr,"Error: special '-' must be either 1 or 2 symbols - %d seen\n",gap);
              Fatal(i);
            }

            if(prev){
              edge = AllocateWLNEdge(ring->locants[prev->ch],prev);
              if(!edge)
                Fatal(i);
            }

            i+= gap+1;
            wln_ptr+= gap+1;

            pending_unsaturate = 0;
            prev = curr;
          }
        }
        break;


      // do a lookahead here as well for multiplier symbols

      case '/':
        if (pending_J_closure){
          break;
        }
        else if (pending_diazo){
          fprintf(stderr,"Error: diazo assignment followed by a multiplier is a disallowed bond type\n");
          Fatal(i);
        }
        else{

          fprintf(stderr,"Error: multipliers are not currently supported\n");
          Fatal(i);
          
          if(pending_linker){
            // we are expected to look for a linker to evaulate

            // some trick recusion is going to go here, but should generalise the linking

            char local_arr [strlen(wln_ptr)+1]; 
            memset(local_arr,'\0',strlen(wln_ptr)+1);
            memcpy(local_arr,wln_ptr,strlen(wln_ptr));
            const char *local = local_arr;
            
            unsigned char local_ch = *(++local); // moved it over
            unsigned int gap = 0; 
            bool found_next = false;

            while(local_ch != '\0'){
              if(local_ch == ' ')
                break;
              if(local_ch == '-'){
                // this calculates the gap to the next '-'
                found_next = true;
                break;
              }
              special.push_back(local_ch);
              gap++;
              local_ch = *(++local);
            }

          }
        }
        
        break;

      default:
        fprintf(stderr, "Error: unallowed character! - alphabet: [A-Z][0-1][&-/' ']\n");
        Fatal(i);
      }

      i++;
      ch = *(++wln_ptr);
    }
    

    // if a multiplier is here, we handle it
    if(pending_linker){

      fprintf(stderr,"will catch here\n");

    }



    if (pending_J_closure)
    {
      fprintf(stderr, "Error: expected 'J' to close ring\n");
      Fatal(len);
    }

    if (pending_locant)
    {
      fprintf(stderr, "Error: expected locant to attach to ring\n");
      Fatal(len);
    }

    if (pending_inline_ring)
    {
      fprintf(stderr, "Error: expected inline ring to be defined\n");
      Fatal(len);
    }

    if (pending_spiro)
    {
      fprintf(stderr, "Error: expected sprio ring to be defined\n");
      Fatal(len);
    }

    if(!AssignCharges(ionic_charges))
      Fatal(len);

      // use this for recursion on multipliers
    return true;
  }

  /* dump wln tree to a dotvis file */
  void WLNDumpToDot(FILE *fp)
  {  
    fprintf(fp, "digraph WLNdigraph {\n");
    fprintf(fp, "  rankdir = LR;\n");
    for (unsigned int i=0; i<=symbol_count;i++)
    {
      WLNSymbol *node = SYMBOLS[i];
      if(!node)
        continue;
  
      fprintf(fp, "  %d", index_lookup[node]);
      if (node->ch == '*')
        fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
      else if(node->type == LOCANT)
        fprintf(fp, "[shape=circle,label=\"%c\",color=blue];\n", node->ch);
      else if (node->type == RING)
        fprintf(fp, "[shape=circle,label=\"%c\",color=green];\n", node->ch);
      else{
        if(std::isdigit(node->ch)){
          if (!node->special.empty())
            fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
          else
            fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->ch);
        } 
        else
          fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->ch);
      }
    
        
      WLNEdge *edge = 0;
      for (edge = node->bonds;edge;edge = edge->nxt){

        WLNSymbol *child = edge->child;
        unsigned int bond_order = edge->order;

        // aromatic
        if(bond_order == 4){
          fprintf(fp, "  %d", index_lookup[node]);
          fprintf(fp, " -> ");
          fprintf(fp, "%d [arrowhead=none,color=red]\n", index_lookup[child]);
        }
        else if (bond_order > 1){
          for (unsigned int k=0;k<bond_order;k++){
            fprintf(fp, "  %d", index_lookup[node]);
            fprintf(fp, " -> ");
            fprintf(fp, "%d [arrowhead=none]\n", index_lookup[child]);
          }
        }
        else{
          fprintf(fp, "  %d", index_lookup[node]);
          fprintf(fp, " -> ");
          fprintf(fp, "%d [arrowhead=none]\n", index_lookup[child]);
        }
      }
    }

    fprintf(fp, "}\n");
  }
};



// holds all the functions for WLN graph conversion, mol object is assumed ALIVE AT ALL TIMES
// uses old NM functions from previous methods: Copyright (C) NextMove Software 2019-present
struct BabelGraph{


  BabelGraph(){};
  ~BabelGraph(){};


  OpenBabel::OBAtom* NMOBMolNewAtom(OpenBabel::OBMol* mol, unsigned int elem,unsigned int charge=0,unsigned int hcount=0)
  {

    OpenBabel::OBAtom* result = mol->NewAtom();
    
    result->SetAtomicNum(elem);
    result->SetImplicitHCount(hcount);

    if(charge)
      result->SetFormalCharge(charge);

    return result;
  }

  void NMOBAtomSetAromatic(OpenBabel::OBAtom* atm, bool arom)
  {
    OpenBabel::OBMol* mol = (OpenBabel::OBMol*)atm->GetParent();
    if (mol && !mol->HasAromaticPerceived())
        mol->SetAromaticPerceived();

    atm->SetAromatic(arom);
  }


  bool NMOBMolNewBond(OpenBabel::OBMol* mol,
                      OpenBabel::OBAtom* s,
                      OpenBabel::OBAtom* e,
                      unsigned int order, bool arom)
  {
    
    if(!s || !e){
      fprintf(stderr,"Error: could not find atoms in bond, bond creation impossible\n");
      return false;
    }

    if(opt_debug)
      fprintf(stderr,"  bonding: atoms %3d --> %3d [%d]\n",s->GetIdx(),e->GetIdx(),order);
    

    if (!mol->AddBond(s->GetIdx(), e->GetIdx(), order)){
      fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n",s->GetIdx(),e->GetIdx());
      return false;
    }
        
    OpenBabel::OBBond* bptr = mol->GetBond(mol->NumBonds() - 1);
    if(!bptr){
      fprintf(stderr,"Error: could not re-return bond for checking\n");
      return false;
    }

    if (arom){
      bptr->SetAromatic();
      NMOBAtomSetAromatic(s,true);
      NMOBAtomSetAromatic(e,true);
    }
    return true;
  }


  bool NMOBSanitizeMol(OpenBabel::OBMol* mol)
  {
    
    mol->SetAromaticPerceived(true);

    if(!OBKekulize(mol)){
      fprintf(stderr,"Error: failed on kekulize mol\n");
      return false;
    }
      
    
    // WLN has no inherent stereochemistry, this can be a flag but should be off by default
    mol->SetChiralityPerceived(true);
    
    return true;
  }


  bool ConvertFromWLN(OpenBabel::OBMol* mol,WLNGraph &wln_graph){

    // aromaticity is handled by reducing the edges

    if(opt_debug)
      fprintf(stderr,"Converting wln to obabel mol object: \n");

    // set up atoms
    for (unsigned int i=1; i<=symbol_count;i++){
      WLNSymbol *sym = SYMBOLS[i];
      if(sym->type != LOCANT){

        OpenBabel::OBAtom *atom = 0;

        unsigned int atomic_num = 0;
        unsigned int charge = 0; 
        unsigned int hcount = 0;

        switch(sym->ch){

          case 'H':
            atomic_num = 1;
            hcount = 0;
            break; 

          case 'B':
            atomic_num = 5;
            break;

          case 'C':
            atomic_num = 6;
            while(sym->num_edges < sym->allowed_edges){
              hcount++;
              sym->num_edges++;
            }
            break;

          case 'X':
            atomic_num = 6;
            hcount = 0; // this must have 4 + methyl
            break;

          case 'Y':
            atomic_num = 6;
            hcount = 1;
            break;

          case 'N':
            atomic_num = 7;
            while(sym->num_edges < sym->allowed_edges){
              hcount++;
              sym->num_edges++;
            }
            break;

          case 'M':
            atomic_num = 7;
            hcount = 1;
            break;

          case 'Z':
            atomic_num = 7; 
            hcount = 2;
            break;

          case 'K':
            atomic_num = 7;
            charge = 1; 
            hcount = 0;
            break;

          case 'O':
            atomic_num = 8;
            if(!sym->num_edges)
              charge = -1;
            break;

          case 'Q':
            atomic_num = 8;
            hcount = 1;
            break;

          case 'F':
            atomic_num = 9;
            if(!sym->num_edges)
              charge = -1;
            break;
          
          case 'P':
            atomic_num = 15;
            while(sym->num_edges < 3){
              hcount++;
              sym->num_edges++;
            }
            break;
          
          case 'S':
            atomic_num = 16;
            while(sym->num_edges < 3){
              hcount++;
              sym->num_edges++;
            }
            break;

          case 'G':
            atomic_num = 17;
            if(!sym->num_edges)
              charge = -1;
            break;

          case 'E':
            atomic_num = 35;
            if(!sym->num_edges)
              charge = -1;
            break;

          case 'I':
            atomic_num = 53;
            if(!sym->num_edges)
              charge = -1;
            break;
        
          case '*':
            atomic_num = special_element_atm(sym->special);
            break;

          default:
            fprintf(stderr,"Error: unrecognised WLNSymbol* char in obabel mol build - %c\n",sym->ch);
            return false;
        }

        // ionic notation - overides any given formal charge
        if(charge_additions[sym]){
          charge = charge_additions[sym];
        }

        atom = NMOBMolNewAtom(mol,atomic_num,charge,hcount);
        if(!atom){
          fprintf(stderr,"Error: formation of obabel atom object\n");
          return false;
        }

        if(sym->type == RING)
          atom->SetInRing();

        babel_atom_lookup[index_lookup[sym]] = atom;
        if(opt_debug)
          fprintf(stderr,"  created: atom[%d] - atomic num(%d), charge(%d)\n",atom->GetIdx(),atomic_num,charge);
      }

    }

    // set bonds
    for(unsigned int i=1;i<=symbol_count;i++){
      WLNSymbol *parent = SYMBOLS[i];
      
      // ignore parents; 
      if(parent->type == LOCANT)
        continue;

      unsigned int parent_id = index_lookup[parent];
      OpenBabel::OBAtom *par_atom = babel_atom_lookup[parent_id];

      WLNEdge *e = 0;
      for (e = parent->bonds;e;e = e->nxt){
        
        WLNSymbol *child = e->child;
      
        // skip across locants
        if(child->type == LOCANT){
          e = child->bonds;
          child = e->child;
        }
          

        unsigned int bond_order = e->order;  
   

        unsigned int child_id = index_lookup[child];
        OpenBabel::OBAtom *chi_atom = babel_atom_lookup[child_id];
        if(bond_order == 4){
          if(!NMOBMolNewBond(mol,par_atom,chi_atom,1,true))
            return false;
        }
        else{
          if(!NMOBMolNewBond(mol,par_atom,chi_atom,bond_order,false))
            return false;
        }
      }
    }

    return true;
  }


};



bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol)
{   
  if(!ptr){
    fprintf(stderr,"Error: could not read wln string pointer\n");
    return false;
  }
  else 
    wln = ptr; 

  WLNGraph wln_graph;
  BabelGraph obabel; 

  if(!wln_graph.ParseWLNString(std::string(ptr))){
    fprintf(stderr,"Error: string pass was successful but return nullptr for wln graph\n");
    return false;
  }
  

  // create an optional wln dotfile
  if (opt_wln2dot)
  {
    fprintf(stderr,"Dumping wln graph to wln-graph.dot:\n");
    FILE *fp = 0;
    fp = fopen("wln-graph.dot", "w");
    if (!fp)
    {
      fprintf(stderr, "Error: could not create dump .dot file\n");
      fclose(fp);
      return 1;
    }
    else
    {
      wln_graph.WLNDumpToDot(fp);
      fclose(fp);
    }
    fprintf(stderr,"  dumped\n");
  }

  if(!wln_graph.ExpandWLNGraph())
    return false;

  if(!obabel.ConvertFromWLN(mol,wln_graph))
    return false;
  

  if(!obabel.NMOBSanitizeMol(mol))
    return false; 

  return true;
}


static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser reads and evaluates wiswesser\n"
                  " line notation (wln), the parser is native\n"
                  " and will can return either a reformatted string*\n"
                  " *if rules do not parse exactly, and the connection\n"
                  " table which can be used in other libraries\n");
  exit(1);
}

static void DisplayUsage()
{
  fprintf(stderr, "wln-writer <options> < input (escaped) >\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -a | --allow-changes          allow changes to notation to allow parsing\n");
  fprintf(stderr, "  -c | --convert                convert the wln graph into SCT table\n");
  fprintf(stderr, "  -d | --debug                  print debug messages to stderr\n");
  fprintf(stderr, "  -h | --help                   print debug messages to stderr\n");
  fprintf(stderr, "  -w | --wln2dot                dump wln trees to dot file in [build]\n");
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i, j;

  cli_inp = (const char *)0;
  dotfile = (const char *)0;

  if (argc < 2)
    DisplayUsage();

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];

    if (ptr[0] == '-' && ptr[1])
      switch (ptr[1])
      {

      case 'a':
        opt_allow = true;
        break;

      case 'c':
        opt_convert = true;
        break;

      case 'd':
        opt_debug = true;
        break;

      case 'h':
        DisplayHelp();

      case 'w':
        opt_wln2dot = true;
        break;

      case '-':
        if (!strcmp(ptr, "--allow-changes"))
        {
          opt_allow = true;
          break;
        }
        else if (!strcmp(ptr, "--convert"))
        {
          opt_convert = true;
          break;
        }
        else if (!strcmp(ptr, "--debug"))
        {
          opt_debug = true;
          break;
        }
        else if (!strcmp(ptr, "--help"))
        {
          DisplayHelp();
        }
        else if (!strcmp(ptr, "--wln2dot"))
        {
          opt_wln2dot = true;
          break;
        }

      default:
        fprintf(stderr, "Error: unrecognised input %s\n", ptr);
        DisplayUsage();
      }

    else
      switch (j++)
      {
      case 0:
        cli_inp = ptr;
        break;
      default:
        break;
      }
  }

  return;
}

int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);
  
  std::string res;
  OpenBabel::OBMol* mol = new OpenBabel::OBMol;
  if(!ReadWLN(cli_inp,mol))
    return 1;


  OpenBabel::OBConversion conv;
  conv.SetOutFormat("smi");
  res = conv.WriteString(mol);

  std::cout << res;

  delete mol; 
  return 0;
}