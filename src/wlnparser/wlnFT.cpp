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

#include "huffman.h"
#include <cctype>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <set>
#include <stack>
#include <deque>
#include <queue>
#include <string>
#include <vector>
#include <map>

#include <utility> // std::pair
#include <iterator>


#define STRUCT_COUNT 1024
#define MAX_EDGES 8

// --- DEV OPTIONS  ---
#define OPT_DEBUG 1
#define OPT_CORRECT 0

const char *cli_inp;

const char *wln_string;
struct WLNSymbol;
struct WLNEdge; 
struct WLNRing;
struct WLNGraph;
struct ObjectStack;

// ##############################################################


int isNumber(const std::string& str)
{
  char* ptr;
  unsigned int val = strtol(str.c_str(), &ptr, 10);
  if(*ptr != '\0')
    return -1;
  else
    return val;
}

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
    res.push_back(wln_string[i]);
  }
  return res; 
}

bool Fatal(unsigned int pos, const char *message)
{ 

#if ERRORS == 1
  fprintf(stderr,"%s\n",message);
  fprintf(stderr, "Fatal: %s\n", wln_string);
  fprintf(stderr, "       ");
  for (unsigned int i = 0; i < pos; i++)
    fprintf(stderr, " ");

  fprintf(stderr, "^\n");
#endif

  return false;
}


/**********************************************************************
                          STRUCT DEFINTIONS
**********************************************************************/
 

struct WLNEdge{
  WLNSymbol *parent;
  WLNSymbol *child;
  WLNEdge *reverse; // points to its backwards edge 
  unsigned char order;
  bool aromatic;

  WLNEdge(){
    parent   = 0;
    child    = 0;
    order    = 0;
    aromatic = 0;
  }
  ~WLNEdge(){};
};


struct WLNSymbol
{
  unsigned short int id;
  unsigned int str_position; 
  short int charge;

  unsigned char ch;
  std::string special; // string for element, or ring, if value = '*'
  
  bool aromatic; 
  WLNRing *inRing; // allows quick lookback 
  unsigned char allowed_edges; // save the bits here as well
  unsigned char num_edges;
  
  unsigned char barr_n; // doesnt need to be large
  unsigned char parr_n; // doesnt need to be large
  WLNEdge bond_array[MAX_EDGES]; // move to undirected 
  WLNEdge prev_array[MAX_EDGES]; // points backwards 

  // if default needed
  WLNSymbol()
  {
    id = 0;
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
    inRing = 0;
    aromatic = 0;
    charge = 0;
    str_position = 0; 
    barr_n = 0; 
    parr_n = 0; 
  }
  ~WLNSymbol(){};

  void add_special(unsigned int s, unsigned int e)
  {
    for (unsigned int i = s; i <= e; i++)
      special.push_back(wln_string[i]);
  }

};

struct WLNRing
{
  unsigned int rsize;
  unsigned int aromatic_atoms;
  unsigned int *adj_matrix; 
  
  std::vector<unsigned char> assignment_locants;
  std::vector<unsigned int> assignment_digits;
  std::map<unsigned char, WLNSymbol *>      locants; 
  std::map<WLNSymbol*,unsigned char>        locants_ch;

  bool spiro;
  unsigned int ranking; // saves the ranking for inline canonical choice
  unsigned int multi_points;
  unsigned int pseudo_points;
  unsigned int bridge_points;
  unsigned int loc_count;
  std::string str_notation; // used for write back

  WLNRing(){
    rsize = 0;
    aromatic_atoms = 0;
    adj_matrix = 0;
    spiro = 0;
    multi_points = 0;
    pseudo_points = 0;
    bridge_points = 0; 
    loc_count = 0; 
    ranking = 0; 
  }
  ~WLNRing(){
    if(adj_matrix)
      free(adj_matrix);
    adj_matrix = 0;
  };

  bool FillAdjMatrix(){

    aromatic_atoms = 0;
    adj_matrix = (unsigned int*)malloc(sizeof(unsigned int) * (rsize*rsize)); 
    if(!adj_matrix)
      return false;
    
    for (unsigned int i = 0; i< rsize;i++){
      for (unsigned int j=0; j < rsize;j++){
        adj_matrix[i * rsize + j] = 0; 
      }
    }

    for (unsigned int i = 0; i< rsize;i++){
      unsigned int r = i;
      unsigned char loc_a = int_to_locant(i+1);
      WLNSymbol *rsym = locants[loc_a]; 
      if(rsym->ch == 'S') // for now lets see
        continue;

      if(rsym->aromatic && rsym->num_edges < rsym->allowed_edges){
        for(unsigned int ei = 0;ei < rsym->barr_n;ei++){
          WLNEdge *redge = &rsym->bond_array[ei];
          WLNSymbol *csym = redge->child; 

          if(csym->ch == 'S' || redge->order > 1)
            continue;
        
          if(csym->aromatic && redge->aromatic && csym->num_edges < csym->allowed_edges){
            unsigned char loc_b = locants_ch[csym];
            unsigned int c = locant_to_int(loc_b-1);
            adj_matrix[r * rsize + c] = 1; 
            adj_matrix[c * rsize + r] = 1; 
            aromatic_atoms++;
          }
        } 
      }
    }

    return true;
  }
};


// handles all memory and 'global' vars
struct WLNGraph
{
  WLNSymbol *root;

  unsigned int symbol_count;
  unsigned int ring_count;

  WLNSymbol *SYMBOLS[STRUCT_COUNT];
  WLNRing   *RINGS  [STRUCT_COUNT];
  
  unsigned int last_cycle_seen=0;  
  
  WLNGraph(){
    symbol_count = 0;
    ring_count   = 0;

    // pointer safety
    root = 0;
    for (unsigned int i = 0; i < STRUCT_COUNT;i++){
      SYMBOLS[i] = 0;
      RINGS[i] = 0;
    }
  };

  ~WLNGraph(){
    for (unsigned int i = 0; i < STRUCT_COUNT;i++){
      delete SYMBOLS[i];
      delete RINGS[i];
    }
  }
};


// needs to be able to hold both a WLNSymbol and WLNRing for branch returns
struct ObjectStack{  
  std::vector<std::pair<WLNRing*,WLNSymbol*>> stack; // vector so i can iterate down and instant access
  WLNRing   *ring;
  WLNSymbol *branch;
  unsigned int size;

  ObjectStack(){
    ring = 0;
    branch = 0;
    size = 0;
  }

  void reserve(unsigned int n){
    stack.reserve(n);
  }

  bool peek(){
    if(!size){
      fprintf(stderr,"Error: peeking empty stack\n");
      return false;
    }
    else{
      fprintf(stderr,"top: ring: %p   branch: %p\n",stack.back().first,stack.back().second);
      return true;
    }
     
  }

  bool pop(){
    if(!size){
#if ERRORS == 1
      fprintf(stderr,"Error: popping empty stack\n");
#endif
      return false;
    }

    stack.pop_back();
    size--;

    ring = 0;
    branch = 0;

    if(stack.empty())
      return false;
    
    for (int i=size-1;i>-1;i--){
      if(!ring && stack[i].first)
        ring = stack[i].first;
      
      if(!branch && stack[i].second)
        branch = stack[i].second;        
    }

    return true; 
  }

  void push(std::pair<WLNRing*,WLNSymbol*> pair,bool verbose = false){
    stack.push_back(pair);
    if(pair.first)
      ring = pair.first;

    if(pair.second)
      branch = pair.second;

    if(verbose){
      fprintf(stderr,"pushed: ring: %p    branch: %p\n",pair.first,pair.second);
    }

    size++;
  }

  void debug_stack(){
    for(unsigned int i=0;i<size;i++){
      fprintf(stderr,"%p,",stack[i].first);
      if(stack[i].second)
        fprintf(stderr,"%c)\n",stack[i].second->ch);
      else
        fprintf(stderr,"%p)\n",stack[i].second);
    }
  }

  bool empty(){
    if (stack.empty())
      return true;
    else 
      return false;
  }

  void clear_all(){
    ring = 0;
    branch = 0;
    while(!stack.empty())
      stack.pop_back();
    size = 0;
  }

  std::pair<WLNRing*,WLNSymbol*> top(){
    if(stack.empty())
      return {0,0};
    else
      return stack.back();
  }

  // cleans branches
  bool branch_avaliable(){
    if(branch && branch->num_edges < branch->allowed_edges)
      return true;
    else
      return false;
  }

  WLNRing *pop_to_ring(){
    std::pair<WLNRing*,WLNSymbol*> t;
    t = top();
    while(!t.first && !stack.empty()){
      pop();
      t = top();
    }
    return t.first;
  }

};



/**********************************************************************
                         WLNSYMBOL Functions
**********************************************************************/


WLNSymbol *AllocateWLNSymbol(unsigned char ch, WLNGraph &graph)
{

  if(graph.symbol_count >= STRUCT_COUNT){
#if ERRORS == 1
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
#endif
    return 0;
  }

  WLNSymbol *wln = new WLNSymbol;
  wln->id = graph.symbol_count++;
  graph.SYMBOLS[wln->id] = wln;
  wln->ch = ch;
 
  return wln;
}

bool IsTerminator(WLNSymbol *symbol){
  switch(symbol->ch){
    case 'E':
    case 'F':
    case 'G':
    case 'H':
    case 'I':
    case 'Q':
    case 'Z':
      return true;
  }
  return false;
}

WLNSymbol* define_hypervalent_element(unsigned char sym, WLNGraph &graph){

  WLNSymbol *new_symbol = 0;

  if(!sym)
    return new_symbol;
  
  switch(sym){
    case 'O':
      new_symbol = AllocateWLNSymbol(sym,graph);
      if(new_symbol)
        new_symbol->allowed_edges = 3;
      break;

    case 'P':
    case 'S':
      new_symbol = AllocateWLNSymbol(sym,graph);
      if(new_symbol)
        new_symbol->allowed_edges = 8;
      break;

    case 'G':
    case 'E':
    case 'I':
    case 'F':
    case 'B':
      new_symbol = AllocateWLNSymbol(sym,graph);
      if(new_symbol)
        new_symbol->allowed_edges = 6;         // allows FCl6
      break;
  }
  
  return new_symbol;
}

/* allocate new or override exisiting node*/
WLNSymbol* define_element(std::string special, WLNGraph &graph){
    
  WLNSymbol *created_wln = 0;
  
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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'B':
      switch(special[1]){
        case 'A':
        case 'E':
        case 'H':
        case 'I':
        case 'K':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;
      

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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;
      
    case 'D':
      switch(special[1]){
        case 'B':
        case 'S':
        case 'Y':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'E':
      switch(special[1]){
        case 'R':
        case 'S':
        case 'U':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'F':
      switch(special[1]){
        case 'E':
        case 'L':
        case 'M':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'G':
      switch(special[1]){
        case 'A':
        case 'D':
        case 'E':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'H':
      switch(special[1]){
        case 'E':
        case 'F':
        case 'G':
        case 'O':
        case 'S':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'I':
      switch(special[1]){
        case 'N':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'K':
      switch(special[1]){
        case 'R':
        case 'A':
          created_wln = AllocateWLNSymbol('*',graph);
          break;

        default:
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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'M':
      switch(special[1]){
        case 'C':
        case 'D':
        case 'G':
        case 'N':
        case 'O':
        case 'T':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;


    case 'O':
      switch(special[1]){
        case 'O':
        case 'G':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'R':
      switch(special[1]){
        case 'A':
        case 'B':
        case 'E':
        case 'F':
        case 'G':
        case 'H':
        case 'N':
        case 'U':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;
     

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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;


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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'U':
      if(special[1] == 'R')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;
      
      break;

    case 'V':
      if (special[1] == 'A')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;
      
      break;
    
    case 'W':
      if(special[1] == 'T')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;
      
      break;
    

    case 'X':
      if (special[1] == 'E')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;

      break;

    case 'Y':
      switch(special[1]){
        case 'B':
        case 'T':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'Z':
      switch(special[1]){
        case 'N':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
           
        default:
          return (WLNSymbol *)0;
      }
      break;

    default:
      return (WLNSymbol *)0;
  }

  created_wln->special = special;
  created_wln->allowed_edges = 8; // allow anything for now;
  return created_wln;
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
        return 87;
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
        return 36;
      else if(special[1] == 'A')
        return 19;
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
        return 78;
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
      return 0;
  }

  return 0;
}

unsigned int count_children(WLNSymbol *sym){
  unsigned int count = 0;
  count += sym->barr_n; 
  count += sym->parr_n;
  return count; 
}


// this one pops based on bond numbers
WLNSymbol *return_object_symbol(ObjectStack &branch_stack){
  WLNSymbol *top = 0;
  while(!branch_stack.empty()){
    top = branch_stack.top().second;
    if(!top)
      return top; // only iterate to the next
    else if(top->ch == 'Y' && count_children(top)==3)
      branch_stack.pop();
    else if(top->num_edges == top->allowed_edges)
      branch_stack.pop();
    else
      return top;
  }
  return top;
}


/**********************************************************************
                          WLNEdge Functions
**********************************************************************/

/* adds an edge between two symbols, arrays should be held by each symbol */
bool AddEdge(WLNSymbol *child, WLNSymbol *parent)
{
  if(!child || !parent || child == parent){
#if ERRORS == 1 
    fprintf(stderr,"Error: binding invalid nodes\n"); 
#endif
    return 0;
  }
  
  if(parent->barr_n >= MAX_EDGES){
#if ERRORS == 1 
    fprintf(stderr,"Error: creating more %d bonds on a singular symbol - is this reasonable?\n",MAX_EDGES);
#endif
    return 0;
  }
  
  if ( ((child->num_edges + 1) > child->allowed_edges)){
#if ERRORS == 1
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges+1, child->allowed_edges);
#endif
    return 0;
  }
  
  if ( ((parent->num_edges + 1) > parent->allowed_edges)){
#if ERRORS == 1
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges+1, parent->allowed_edges);
#endif
    return 0;
  }

  child->num_edges++;
  parent->num_edges++;
  
  WLNEdge *forward = &parent->bond_array[parent->barr_n];
  forward->child = child;
  forward->parent = parent;
  forward->order = 1; 
  parent->barr_n++;
  
  WLNEdge *backward = &child->prev_array[child->parr_n]; 
  // store the reverse in the prev array
  backward->child = parent; 
  backward->parent = child; 
  backward->order = 1;
  child->parr_n++; 

  forward->reverse = backward;
  backward->reverse = forward; 
  
  return true;
}


WLNEdge *search_edge(WLNSymbol *child, WLNSymbol*parent){
  if(!child || !parent)
    return 0;
  
  for (unsigned int ei=0;ei<child->barr_n;ei++){
    if(child->bond_array[ei].child == parent)
      return &child->bond_array[ei];
  }

  for (unsigned int ei=0;ei<parent->barr_n;ei++){
    if(parent->bond_array[ei].child == child)
      return &parent->bond_array[ei];
  }

  return 0;
}

bool unsaturate_edge(WLNEdge *edge,unsigned int n, unsigned int pos=0){
  if(!edge)
    return 0;
  
  edge->order += n;
  edge->reverse->order = edge->order; 
  edge->parent->num_edges += n;
  edge->child->num_edges+= n;

  if( (edge->child->num_edges > edge->child->allowed_edges)){
#if ERRORS == 1
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->child->ch,edge->child->num_edges, edge->child->allowed_edges);
#endif
    return false;
  }

  if( (edge->parent->num_edges > edge->parent->allowed_edges)){
#if ERRORS == 1
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->parent->ch,edge->parent->num_edges, edge->parent->allowed_edges);
#endif
    return false;
  }

  return true;
}

bool saturate_edge(WLNEdge *edge,unsigned int n){
  if(!edge)
    return false;
  
  if(edge->order < 2)
    return true;
  
  edge->order -= n; 
  edge->reverse->order = edge->order; 
  edge->parent->num_edges -= n;
  edge->child->num_edges -= n;

  return true;
}


bool add_methyl(WLNSymbol *head, WLNGraph &graph){
  WLNSymbol *carbon = AllocateWLNSymbol('1',graph);
  if(carbon)
    carbon->allowed_edges = 4;
  else 
    return false;

  return AddEdge(carbon, head); 
}

bool has_dioxo(WLNSymbol *node){
  if(node->parr_n && node->prev_array[0].child->ch == 'W')
    return true;
  
  for(unsigned int ei=0;ei<node->barr_n;ei++)
    if(node->bond_array[ei].child->ch == 'W')
      return true;
  
  return false;
}

bool add_dioxo(WLNSymbol *head,WLNGraph &graph){

  WLNEdge *edge = 0;
  WLNSymbol *oxygen = 0;
  WLNSymbol *binded_symbol = 0; 

  // 2 options here, either the symbol is binded 'W' -> sym
  // or the symbol is binded sym -> 'W'

  if(head->barr_n){
    binded_symbol = head->bond_array[0].child;
    edge = &head->bond_array[0]; 
  }
  else if(head->parr_n){
    binded_symbol = head->prev_array[0].child;
    edge = &head->prev_array[0]; 
  }

  if(!binded_symbol || edge->order != 3){
    return false; 
  }
  
  head->ch = 'O'; // change the W symbol into the first oxygen
  head->allowed_edges = 2;

  // should still be double bonded. 
  oxygen = AllocateWLNSymbol('O',graph);
  oxygen->allowed_edges = 2;
  if(!saturate_edge(edge,1))
    return false;

  if(!AddEdge(oxygen, binded_symbol))
    return false;

  WLNEdge *sedge = &binded_symbol->bond_array[binded_symbol->barr_n-1]; 
  
  if(binded_symbol->num_edges < binded_symbol->allowed_edges)
    if(!unsaturate_edge(sedge,1))
      return false;
  
  if(binded_symbol->ch == 'N'){
    binded_symbol->charge++;
  }

  if(!edge || !sedge)
    return false;
  
  return true;
}


/* resolve carbon methyl assumptions */
bool resolve_methyls(WLNSymbol *target, WLNGraph &graph){

  switch(target->ch){

    case 'X':
    case 'K':
      while(target->num_edges < target->allowed_edges){
        if(!add_methyl(target,graph))
          return false;
      }
      target->num_edges = target->allowed_edges;
      break;

    case 'Y':
      while(count_children(target) < 3){
        if(!add_methyl(target,graph))
          return false;
      }
      target->num_edges = target->allowed_edges;
      break;

    default:
      return false;
  }

  return true;
}



/**********************************************************************
                          WLNRing Functions
**********************************************************************/

WLNRing *AllocateWLNRing(WLNGraph &graph)
{
  if(graph.ring_count >= STRUCT_COUNT)
    return 0;
  
  WLNRing *wln_ring = new WLNRing;
  graph.RINGS[graph.ring_count++] = wln_ring;
  return wln_ring;
}


// both lookups needed for QOL in ring building
WLNSymbol* assign_locant(unsigned char loc,WLNSymbol *locant, WLNRing *ring){
  if(!locant)
    return 0;
  ring->locants[loc] = locant; 
  ring->locants_ch[locant] = loc;
  locant->inRing = ring;
  return locant; 
}  


bool set_up_broken( WLNRing *ring, WLNGraph &graph,
                    std::set<unsigned char> &broken_locants,
                    std::map<unsigned char,std::deque<unsigned char>> &broken_lookup,
                    std::map<unsigned char, bool> &spawned_broken, 
                    std::map<unsigned char,unsigned int> &allowed_connections)
{
  // broken locants
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
      else if(pos > 3)
        return false;
      
      if(OPT_DEBUG)
        fprintf(stderr,"  ghost linking %d to parent %c\n",loc_broken,parent);
      
      if(!ring->locants[loc_broken]){
        // bond them in straight away
        allowed_connections[loc_broken] = 3; 
        if(allowed_connections[parent])
          allowed_connections[parent]--;

        WLNSymbol *broken = AllocateWLNSymbol('C',graph);
        broken->inRing = ring;
        broken->allowed_edges = 4;
        broken = assign_locant(loc_broken,broken,ring);
        broken_lookup[parent].push_back(loc_broken);
        if(!AddEdge(ring->locants[loc_broken], ring->locants[parent]))
          return false;
      }
      else
        return false;
      
    }
  }

  return true; 
}

bool set_up_pseudo( WLNRing *ring, WLNGraph &graph,
                            std::vector<unsigned char> &pseudo_locants,
                            std::map<unsigned char,unsigned char> &pseudo_lookback)
{ 

  if(!pseudo_locants.empty()){

    if(pseudo_locants.size() % 2 != 0)
      return false;
    

    for(unsigned int i=0;i<pseudo_locants.size()-1;i+=2){
      unsigned int bind_1 = pseudo_locants[i];
      unsigned int bind_2 = pseudo_locants[i+1];
      pseudo_lookback[bind_2] = bind_1;
    }
  }
  
  return true;
}

/* interesting here that the multicyclic points are not explicitly used */
unsigned int BuildCyclic( std::vector<std::pair<unsigned int,unsigned char>>  &ring_assignments, 
                          std::vector<bool>                                   &aromaticity,
                          std::vector<unsigned char>                          &multicyclic_locants,
                          std::vector<unsigned char>                          &pseudo_locants,
                          std::set<unsigned char>                             &broken_locants,
                          std::map<unsigned char,unsigned int>                &bridge_locants,
                          unsigned char size_designator,
                          WLNRing *ring,
                          WLNGraph &graph) 
{
  unsigned int local_size = 0;
  if(!size_designator){
 
    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i]; 
      if(local_size)
        local_size += component.first - 2;
      else
        local_size = component.first;
    }

    for (unsigned int i=0;i<252;i++){
      if(bridge_locants[i])
        local_size+= -1; 
    }

    local_size+= - broken_locants.size();
    if(OPT_DEBUG)
      fprintf(stderr,"  calculated size: %c(%d)\n",int_to_locant(local_size),local_size);
  }
  else
    local_size = locant_to_int(size_designator);

#if BUILD_RING
  // create all the nodes in a large straight chain and assign how many bonds
  // each atom is allowed to take

  WLNSymbol *curr = 0; 
  WLNSymbol *prev = 0; 
  std::map<unsigned char,unsigned int> allowed_connections; 

  for (unsigned int i=1;i<=local_size;i++){
    unsigned char loc = int_to_locant(i);

    if(i == 1 || i == local_size)
      allowed_connections[loc] = 2; 
    else
      allowed_connections[loc] = 1; 

    if(!ring->locants[loc]){
      curr = AllocateWLNSymbol('C',graph);
      if(!curr)
        return 0;

      curr->allowed_edges = 4;
      curr->inRing = ring;
      curr = assign_locant(loc,curr,ring);
    }
    else{
      curr = ring->locants[loc];
      if(curr->ch == 'X')
        allowed_connections[loc]++;
      else if(curr->ch == '*')
        allowed_connections[loc] = 6; // allow octahedral geometry 

      if(!ring->locants_ch[curr])
        ring->locants_ch[curr] = loc;
    }

    if(bridge_locants[loc] && allowed_connections[loc])
      allowed_connections[loc]--; 
      
    if(prev){
      if(!AddEdge(curr, prev))
        return false;
    }
    prev = curr;
  }

  std::map<unsigned char,unsigned char>             pseudo_lookup;  
  std::map<unsigned char,std::deque<unsigned char>> broken_lookup;    // want to pop front
  std::map<unsigned char, bool>                     spawned_broken; 
  
  std::map<unsigned char, bool>                     shortcuts; // take when possible?

  // broken locant map spawn + pseudo locants map spawn 
  if(!set_up_broken(ring,graph,broken_locants,broken_lookup,spawned_broken,allowed_connections) || 
     !set_up_pseudo(ring,graph,pseudo_locants,pseudo_lookup))
    return false;

  // calculate bindings and then traversals round the loops
  bool aromatic             = false;
  unsigned char bind_1      = '\0';
  unsigned char bind_2      = '\0';
  unsigned int fuses        = 0; 
  unsigned int comp_size    = 0;
  unsigned int pseudo_pairs = pseudo_locants.size()/2;
  
  for (unsigned int i=0;i<ring_assignments.size();i++){
    std::pair<unsigned int, unsigned char> component = ring_assignments[i];
    comp_size = component.first;
    bind_1    = component.second;
    aromatic = aromaticity[i];
    WLNSymbol *path = ring->locants[bind_1];

    if(!path){
#if ERRORS == 1
      fprintf(stderr,"Error: out of bounds locant access in cyclic builder\n");
#endif
      return 0;
    }

    // If we need to catch fuse, we do
    if(i==ring_assignments.size()-1 && pseudo_pairs){
      bool caught = false;
      for(unsigned int s = 1; s<=local_size;s++){
        if(pseudo_lookup[int_to_locant(s)]){
          unsigned char pbind_2 = int_to_locant(s); 
          unsigned char pbind_1 = pseudo_lookup[pbind_2];

          if(OPT_DEBUG)
            fprintf(stderr,"  %d  catch fusing: %c <-- %c\n",fuses,pbind_2,pbind_1);
          
          if(!search_edge(ring->locants[pbind_2],ring->locants[pbind_1])){
            if(!AddEdge(ring->locants[pbind_2], ring->locants[pbind_1]))
              return false;
            fuses++;
            caught = true;
          }
        }
      }
      if(caught)
        break;
    }

    // --- MULTI ALGORITHM --- 
    unsigned int path_size = 0; 
    unsigned char *ring_path = (unsigned char*)malloc(sizeof(unsigned char)*comp_size);
    for (unsigned int a=0;a<comp_size;a++)
      ring_path[a] = 0;
    
    ring_path[path_size++] = ring->locants_ch[path];
    unsigned char highest_loc = '\0'; 

    while(path_size < comp_size){ 
      // this should now overshoot

      highest_loc = '\0'; // highest of each child iteration 
      for(unsigned int ei=0;ei<path->barr_n;ei++){
        WLNSymbol *child = path->bond_array[ei].child;
        unsigned char child_loc = ring->locants_ch[child];

        if(child_loc > 128 && !spawned_broken[child_loc])  // skip the broken child if not yet included in a ring
          continue;
        else if(shortcuts[child_loc]){
          highest_loc = child_loc;
          break;
        }
        else if(child_loc >= highest_loc)
          highest_loc = child_loc;  
      }

      if(!highest_loc){
        if(locant_to_int(ring->locants_ch[path]) == local_size)
          highest_loc = ring->locants_ch[path];
        else{
#if ERRORS == 1
          fprintf(stderr,"Error: locant path formation is broken in ring definition - '%c(%d)'\n",ring->locants_ch[path],ring->locants_ch[path]);
#endif    
          if(ring_path)
            free(ring_path);
          ring_path = 0;
          return false;
        }
      }

      path = ring->locants[highest_loc];
      ring_path[path_size++] = highest_loc;   

      // let catch fuse take care of last pseudo pairs
      if(pseudo_lookup[highest_loc] != '\0' && path_size < comp_size && i!=ring_assignments.size()-1){
        // lets get the bonds right and then worry about the path 
        
        bind_1 = pseudo_lookup[highest_loc];
        bind_2 = highest_loc; 
        path_size = comp_size;
        for(unsigned int a = 0;a<comp_size;a++)
          ring_path[a] = 0;
        
        // just reset
        pseudo_lookup[highest_loc] = 0;
        if(bind_1 > 128)
          spawned_broken[bind_1] = true;
        
        shortcuts[bind_1] = true;
        if(pseudo_pairs)
          pseudo_pairs--;
      }

      bind_2 = highest_loc;
    }

    // sanitize the ring path, if we've got duplicates (i.e hit the end of the path),
    // we decrement and loop back round
    for(unsigned int i=0;i<path_size;i++){
      if(ring_path[i] == int_to_locant(local_size)){
        unsigned int tally = 1;
        for(unsigned int j=i+1;j<path_size;j++){
          if(ring_path[j] == ring_path[i]){
            ring_path[j] += -tally; // looping back
            tally++;
          } 
        }
      }
    }

    // shifting now performed here should be more stable
    for(;;){
      if(!broken_lookup[bind_1].empty()){
        
        while(spawned_broken[broken_lookup[bind_1].front()])
          broken_lookup[bind_1].pop_front();
        
        if(broken_lookup[bind_1].empty())
          continue;

        unsigned char bloc = broken_lookup[bind_1].front();
        broken_lookup[bind_1].pop_front();
        bind_1 = bloc;
        
        for(unsigned int a=path_size-1;a>0;a--)
          ring_path[a] = ring_path[a-1];
        
        ring_path[0] = bind_1;
        spawned_broken[bind_1] = true;

        // pseudo check, will remove later
        if(ring_path[path_size-1])
          bind_2 = ring_path[path_size-1]; 
      }
      else if(allowed_connections[bind_1]){
        // very rare case where we get the right path a different way to normal
        
        while((!allowed_connections[bind_2] || bind_2 == bind_1) && bind_2 < int_to_locant(local_size))
          ring_path[path_size-1] = ++bind_2;
        
        if(OPT_DEBUG){
          fprintf(stderr,"  %d  fusing (%d): %c <-- %c   [",fuses,comp_size,bind_2,bind_1);
          for (unsigned int a=0;a<path_size;a++)
            fprintf(stderr," %c(%d)",ring_path[a],ring_path[a]);
          fprintf(stderr," ]\n");
        }

        if(!AddEdge(ring->locants[bind_2], ring->locants[bind_1])){
          if(ring_path)
            free(ring_path);
          ring_path = 0;
          return false;
        }
         
        allowed_connections[bind_1]--;
        if(allowed_connections[bind_2])
          allowed_connections[bind_2]--;

        break;
      }
      else{

        bind_1++; // increase bind_1
        if(bind_1 > int_to_locant(local_size)+1)
           break;

        bool found = false;
        for(unsigned int a=0;a<path_size;a++){
          if(ring_path[a] == bind_1)
            found = true;
        }

        // if its already there, we dont change the path
        if(!found){ 
          // if its not, we have to spawn it in by knocking one off
          for(unsigned int a=path_size-1;a>0;a--)
            ring_path[a] = ring_path[a-1];

          ring_path[0] = bind_1; 
          bind_2 = ring_path[path_size-1]; 
        }
      }
    }


    if(aromatic){
      for(unsigned int a=0;a<path_size;a++){
        WLNSymbol *arom = ring->locants[ring_path[a]]; 
        if(arom){
          arom->aromatic = true;
          ring->aromatic_atoms = 1;
        }
      }
    
      // add the edges based on statement before n^2 but should work
      for(unsigned int a=0;a<path_size;a++){
        WLNSymbol *src = ring->locants[ring_path[a]];
        for(unsigned int b=a+1;b<path_size;b++){
          WLNSymbol *trg = ring->locants[ring_path[b]];
          if(src && trg && src->aromatic && trg->aromatic){
            WLNEdge *edge = search_edge(src,trg); 
            if(edge)
              edge->aromatic = true;
          }
        }
      }

    }

    if(ring_path)
      free(ring_path);
    ring_path = 0; 
    fuses++;
  }
#endif
  return local_size; 
}


unsigned char create_relative_position(unsigned char parent){
  // A = 129
  unsigned int relative = 128 + locant_to_int(parent);
  if(relative > 252){
#if ERRORS == 1
    fprintf(stderr,"Error: relative position is exceeding 252 allowed space - is this is suitable molecule for WLN notation?\n");
#endif
    return '\0';
  }
  else
    return relative;
}


// try to handle if errors are occuring due to path changes
bool post_unsaturate(std::vector<std::pair<unsigned char, unsigned char>> &bonds, 
                        unsigned int final_size,
                        WLNRing *ring){

  // post unsaturate bonds
  for (std::pair<unsigned char, unsigned char> bond_pair : bonds){
    
    unsigned char loc_1 = bond_pair.first;
    unsigned char loc_2 = bond_pair.second;

    if(loc_2 > int_to_locant(final_size)){
      loc_1 = 'A';
      loc_2--;
    }

    WLNEdge *edge = search_edge(ring->locants[loc_2],ring->locants[loc_1]);
    if(!edge)
      return false;
    else if(!unsaturate_edge(edge,1))
      return false;

    edge->aromatic = 0;
    // edge->child->aromatic = 0;
    // edge->parent->aromatic = 0;
  }

  return true;
}

// try to handle if errors are occuring due to path changes
bool post_saturate( std::vector<std::pair<unsigned char, unsigned char>> &bonds, 
                    unsigned int final_size,
                    WLNRing *ring){

  // post saturate bonds
  for (std::pair<unsigned char, unsigned char> bond_pair : bonds){
    
    unsigned char loc_1 = bond_pair.first;
    unsigned char loc_2 = bond_pair.second;

    if(loc_2 > int_to_locant(final_size)){
      loc_1 = 'A';
      loc_2--;
    }

    WLNEdge *edge = search_edge(ring->locants[loc_2],ring->locants[loc_1]);
    if(!edge)
      return false;
    else
      edge->aromatic = false; // removes from kekulize consideration
  }

  return true;
}

/* parse the WLN ring block, use ignore for already predefined spiro atoms */
bool FormWLNRingLazy(WLNRing *ring,std::string &block, unsigned int start, WLNGraph &graph,unsigned char spiro_atom='\0'){
  ring->str_notation = block;

  bool warned             = false;  // limit warning messages to console
  bool heterocyclic       = false;  // L|T designator can throw warnings

  unsigned int state_multi          = 0; // 0 - closed, 1 - open multi notation, 2 - expect size denotation
  unsigned int state_pseudo         = 0; 
  unsigned int state_aromatics      = 0;

  unsigned int  expected_locants      = 0;
  unsigned char ring_size_specifier   = '\0';
  
  bool locant_attached                = false;  // more sensible way of modifying the locants 
  unsigned char positional_locant     = 'A';    // have A as a default, lets tidy around this

  std::string str_buffer;  

  std::vector<bool> aromaticity; 
  std::vector<std::pair<unsigned char, unsigned char>>  unsaturations;
  std::vector<std::pair<unsigned char, unsigned char>>  saturations;

  std::vector<unsigned char>                multicyclic_locants;
  std::vector<unsigned char>                pseudo_locants;
  std::map<unsigned char,unsigned int>      bridge_locants;
  std::set<unsigned char>                   broken_locants;
  
  std::vector<std::pair<unsigned int, unsigned char>>  ring_components;

  // broken locants start at A = 129 for extended ascii 
  // first is the standard 'X-' second is 'X-&', third is 'X--', fourth is 'X--&' and so on

  unsigned int i = 0;    
  unsigned int len = block.size();
  const char *block_str = block.c_str();
  unsigned char ch = *block_str;

  while(ch){

    switch(ch){
      case ' ':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(state_multi == 3)
          state_multi = 0;

        if(expected_locants)
          return Fatal(i+start,"Error: not enough locants before space character");
        
        else if(state_multi == 1)
          state_multi = 2;
        else if (state_pseudo)
          state_pseudo = 0;
        else if(positional_locant && locant_attached){
          if(ring_components.empty())
            return Fatal(start+i,"Error: assigning bridge locants without a ring");
          else
            bridge_locants[positional_locant] = true;
        }

        positional_locant = '\0'; // hard resets on spaces
        locant_attached = false;
        break;


      case '&':
        if (state_aromatics){
          aromaticity.push_back(1);
          break;
        }
        else if (state_multi == 3){
          ring_size_specifier += 23;
        }
        else if (state_pseudo){
          pseudo_locants.back() += 23;
        }
        else if(positional_locant && locant_attached){
          positional_locant += 23;
        }
        else{
          // if unhandled, then it must be an aromatic start
          state_aromatics = 1;
          aromaticity.push_back(1);
        }
        break;

      case '/':
        if(state_aromatics)
          return Fatal(i+start,"Error: invalid character in the aromaticity assignment block");
             
        expected_locants = 2; 
        state_pseudo = 1;
        break; 
      

      // turn this into a look ahead type bias in order to significantly tidy this up
      case '-':{
        
        str_buffer.clear();
        bool found_next = false;
        unsigned int k = i+1; // use the block string and iterate

        if(!expected_locants){
          while(k < block.size()){
            if(block[k] == ' ' || block[k] == '&')
              break;
            if(block[k] == '-'){
              // this calculates the gap to the next '-'
              if(k != i+1)
                found_next = true; // if there was a double '--' this will have gap zero
              break;
            }
            str_buffer.push_back(block[k]);
            k++;
          }
        }

        if(!found_next){
          if(i > 0 && block[i-1] == '&')
            state_aromatics = 1;
          else{

            if(positional_locant && locant_attached){
        
              if(positional_locant < 128){
                positional_locant = create_relative_position(positional_locant); // i believe breaking modifier will then get removed
                if(!positional_locant)
                  return Fatal(i+start, "Error: failed to make expanded locant position");

              }
              else{
                // this means its already been moved, so we move the locant 23+23 across
                if(positional_locant + 46 > 252)
                  return Fatal(start+i,"Error: branching locants are exceeding the 252 space restriction on WLN notation, is this a reasonable molecule?");
                
                positional_locant += 46;
              }

              if(state_pseudo && !pseudo_locants.empty())
                  pseudo_locants.back() = positional_locant;
              if(state_multi == 1 && !multicyclic_locants.empty())
                multicyclic_locants.back() = positional_locant;
            }
          }
        }
        else{
          
          if(str_buffer.size() == 1){
            
            if(positional_locant != spiro_atom){
              WLNSymbol* new_locant = assign_locant(positional_locant,define_hypervalent_element(str_buffer[0],graph),ring);  // elemental definition 
              if(!new_locant)
                return Fatal(i+start, "Error: could not create hypervalent element");
              
              new_locant->str_position = start+i+1+1;
              if(OPT_DEBUG)
                fprintf(stderr,"  assigning hypervalent %c to position %c\n",str_buffer[0],positional_locant);
            }
            else 
              positional_locant++;

            block_str+=2; 
            i+=2; 
            positional_locant++;
            locant_attached = false;
          }
          else if (str_buffer.size() == 2){
    
            if(std::isdigit(str_buffer[0])){
              for(unsigned char dig_check : str_buffer){
                if(!std::isdigit(dig_check))
                  return Fatal(start+i,"Error: mixing numerical and alphabetical special defintions is not allowed");
              }

              int big_ring = isNumber(str_buffer);
              if(big_ring < 0)
                return Fatal(start+i,"Error: non numeric value entered as ring size");
              

              ring_components.push_back({big_ring,positional_locant}); //big ring
              positional_locant = 'A';
              locant_attached = false;
            }
            else if(positional_locant != spiro_atom){

              // must be a special element 
              WLNSymbol* new_locant = assign_locant(positional_locant,define_element(str_buffer,graph),ring);  // elemental definition
              if(!new_locant)
                return Fatal(i+start, "Error: could not create periodic code element");

              new_locant->str_position = start+i + 1+1; // attaches directly to the starting letter
    
              if(OPT_DEBUG)
                fprintf(stderr,"  assigning element %s to position %c\n",str_buffer.c_str(),positional_locant);
              
              positional_locant++;
            }
            else
              positional_locant++;
            
            locant_attached = false;
            
            block_str+=3; 
            i+=3;
          }
          else
            return Fatal(start+i,"Error: ended in an unexpected state due to '-' characters");
        }
      
        break;
      }

      case '0':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(!ring_components.empty()){
          if(!positional_locant)
            positional_locant = 'A';

          if(OPT_DEBUG)
            fprintf(stderr,"  placing pi bond charge on locant - %c\n",positional_locant);

          
//          WLNSymbol *zero_carbon = AllocateWLNSymbol('C',graph);
//          zero_carbon->allowed_edges = 3; 
//          zero_carbon = assign_locant(positional_locant++,zero_carbon,ring);
//          zero_carbon->str_position = (start+i+1); 
//          zero_carbon->charge--; 
        }
        locant_attached = false;
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
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(state_aromatics)
          return Fatal(i+start,"Error: invalid character in the aromaticity assignment block");
        

        if (i > 1 && block[i-1] == ' '){
          state_multi   = 1; // enter multi state
          expected_locants = ch - '0';
        }
        else{
          ring_components.push_back({ch-'0',positional_locant});
          positional_locant = 'A';
          locant_attached = false;
        }
        break;

      // bring this out as chelating notation is seperate
      case 'D':
        if(i == 0){
          heterocyclic = true;
          if(OPT_DEBUG)
            fprintf(stderr,"  opening chelating notation\n");
        }

        if(state_aromatics)
          return Fatal(i+start,"Error: invalid character in the aromaticity assignment block");
        

        if(expected_locants){
          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            return Fatal(start+i,"Error: unhandled locant rule");
          

          positional_locant = ch; // use for look back
          locant_attached = true;
          expected_locants--;
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if(i>0 && block[i-1] == ' '){
          positional_locant = ch;
          locant_attached = true;
        }
        break;

      case 'A':
      case 'B':
      case 'C':
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
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(state_aromatics)
          return Fatal(i+start,"Error: invalid character in the aromaticity assignment block");
        

        if(expected_locants){
          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            Fatal(start+i,"Error: unhandled locant rule");

          positional_locant = ch; // use for look back
          locant_attached = true;
          expected_locants--;
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if (spiro_atom && positional_locant == spiro_atom){
          positional_locant++;
          locant_attached = false;
        }
        else if (positional_locant){
          
          if (OPT_DEBUG)
            fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

          WLNSymbol *new_locant = 0; 
          switch(ch){
            
            case 'S':
            case 'P':
              if(!heterocyclic)
                warned = true;

#if ALLOCATE_HETERO
              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant++,new_locant,ring);
              new_locant->str_position = (start+i+1); 

              if(ch == 'P')
                new_locant->allowed_edges = 5;
              else
                new_locant->allowed_edges = 6;

              new_locant->inRing = ring;
#endif
              break;

            case 'Y':
            case 'X':
            case 'K':
              if(!heterocyclic && ch=='K')
                warned = true;

#if ALLOCATE_HETERO
              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant++,new_locant,ring);
              new_locant->allowed_edges = 4;
              new_locant->inRing = ring;
              new_locant->str_position = (start+i+1); 
#endif
              break;

            case 'Z': // treat as NH2
            case 'N':
            case 'B':
              if(!heterocyclic)
                warned = true;

#if ALLOCATE_HETERO
              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant++,new_locant,ring);
              new_locant->allowed_edges = 3;
              new_locant->inRing = ring;
              new_locant->str_position = (start+i+1); 
#endif
              break;

            case 'M':
            case 'O':
            case 'V':
              if(!heterocyclic && (ch == 'M' || ch == 'O'))
                warned = true;

#if ALLOCATE_HETERO
              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant++,new_locant,ring);
              new_locant->allowed_edges = 2;
              new_locant->inRing = ring;
              new_locant->str_position = (start+i+1); 
#endif
              break;

        
            case 'U':{
              // no need to put this in implied, it has to be specified
#if ALLOCATE_HETERO
              if(i < len - 3 && block[i+1] == '-' && block[i+2] == ' '){
                
                // can double bond to a amped locant
                unsigned int k = 1;
                unsigned char dloc = block[i+3]; 
                while(block[k+i+3] == '&'){
                  dloc+=23;
                  k++;
                } 
                unsaturations.push_back({positional_locant,dloc});
                block_str += 2+k;
                i += 2+k;
              }
              else
                unsaturations.push_back({positional_locant,positional_locant+1});
              
              if(i < len && block[i+1] != 'U')
                positional_locant++;
#endif
              break;
            }

            // externally bonded to the symbol as a locant symbol
            case 'W':{
              if(!heterocyclic)
                warned = true;

#if ALLOCATE_HETERO
              if(positional_locant > 'A')
                positional_locant--;

              // the carbon might not be created yet
              if(!ring->locants[positional_locant]){
                new_locant = AllocateWLNSymbol('C',graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->allowed_edges = 2;
                new_locant->inRing = ring;
                new_locant->str_position = (start+i+1); 
              }
              else
                new_locant = ring->locants[positional_locant];
              
              if(new_locant->ch == 'N')
                new_locant->allowed_edges++;

              WLNSymbol *dioxo = AllocateWLNSymbol('W',graph);
              dioxo->allowed_edges = 3;
              dioxo->inRing = ring;
              if(!AddEdge(dioxo, new_locant))
                return Fatal(start+i,"Error: failed to create bond");

              WLNEdge *e = &new_locant->bond_array[new_locant->barr_n-1];  
              if(!unsaturate_edge(e,2))
                return Fatal(start+i,"Error: failed to unsaturate edge");
#endif              
              break;
            }

            // has the effect of unaromatising a bond, remove from edge consideration
            case 'H':
              //saturations.push_back({positional_locant,positional_locant+1});
              break;


            default:
              return Fatal(start+i,"Error: invalid character in atom assignment within ring notation");
          }

          locant_attached = false; // locant is no longer primary
        }
        else if(i>0 && block[i-1] == ' '){
          positional_locant = ch;
          locant_attached = true;
        }
        break;

      case 'L':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(state_aromatics)
          return Fatal(i+start,"Error: invalid character in the aromaticity assignment block");
        

        if(i==0){
          heterocyclic = false; 
          break;
        }
        if(expected_locants){

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            return Fatal(start+i,"Error: unhandled locant rule");
          

          positional_locant = ch; // use for look back
          locant_attached = true;
          expected_locants--;
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else{
          if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            locant_attached = true;
          }
          else
            return Fatal(i+start,"Error: symbol is in an unhandled state, please raise issue if this notation is 100%% correct");
          
        }
      
        break;


      case 'T':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(state_aromatics){
          aromaticity.push_back(0);
          break;
        }
      
        if(i==0){
          heterocyclic = true; 
          break;
        }

        if(expected_locants){
          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            return Fatal(start+i,"Error: unhandled locant rule");
          

          positional_locant = ch; // use for look back
          locant_attached = true;
          expected_locants--;
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if(positional_locant && locant_attached){
          if(ring_components.empty())
            return Fatal(start+i,"Error: assigning bridge locants without a ring");
          else
            bridge_locants[positional_locant] = true;

          state_aromatics = 1;
          aromaticity.push_back(0);
        }
        else{
          if(i>0 && block[i-1] == ' ' && block[i+1] != 'J'){
            positional_locant = ch;
            locant_attached = true;
          }
          else{
            // this must be an aromatic state right?
            state_aromatics = 1;
            aromaticity.push_back(0);
          }
        }
        break;

      
      // CLOSE
      case 'J':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(state_aromatics)
          state_aromatics = 0;
        
        if (i == block.size()-1){
          
          if(ring_components.empty())
            return Fatal(start+i,"Error: error in reading ring components, check numerals in ring notation");
          

          if (aromaticity.size() == 1 && aromaticity[0] == false){
            while(aromaticity.size() < ring_components.size())
              aromaticity.push_back(false);
          }
          else if (aromaticity.empty()){
            while(aromaticity.size() < ring_components.size())
              aromaticity.push_back(true);
          }

          // perform the aromatic denotion check
          if (ring_components.size() != aromaticity.size())
            return Fatal(i+start,"Error: mismatch between number of rings and aromatic assignments");
          
          break;
        }
        if(expected_locants){

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            return Fatal(start+i,"Error: unhandled locant rule");
        
          positional_locant = ch; // use for look back
          locant_attached = true;
          expected_locants--;
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if(positional_locant && locant_attached){
          if(ring_components.empty())
            return Fatal(start+i,"Error: assigning bridge locants without a ring");
          
          else
            bridge_locants[positional_locant] = true;
        }
        else{
          if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            locant_attached = true;
          }
          else
            return Fatal(i+start,"Error: symbol is in an unhandled state, please raise issue if this notation is 100%% correct");
        }
        break;

      default:
        break;
    }
    
    i++;
    ch = *(++block_str);
  }

  if(OPT_DEBUG && warned)
    fprintf(stderr,"Warning: heterocyclic ring notation required for inter atom assignment, change starting 'L' to 'T'\n");
  

  // debug here
  if (OPT_DEBUG){
    
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
    for (unsigned int i=0;i<252;i++){
      if(bridge_locants[i])
        fprintf(stderr,"%c ",i);
    }
    fprintf(stderr,"\n");

    if(!pseudo_locants.empty()){
      fprintf(stderr,"  pseudo locants: ");
      for (unsigned int i=0; i< pseudo_locants.size()-1;i+=2)
        fprintf(stderr,"[%c <-- %c] ",pseudo_locants[i],pseudo_locants[i+1]);
      fprintf(stderr,"\n");
    }

    fprintf(stderr,"  multi size: %c(%d)\n",ring_size_specifier ,ring_size_specifier ? locant_to_int(ring_size_specifier) : 0);
    fprintf(stderr,"  heterocyclic: %s\n", heterocyclic ? "yes":"no");
  }
  
  unsigned int final_size = 0; 
  final_size = BuildCyclic(ring_components,aromaticity,
                                multicyclic_locants,pseudo_locants,
                                broken_locants,
                                bridge_locants,
                                ring_size_specifier,
                                ring,
                                graph);

  ring->rsize = final_size;
  ring->multi_points = multicyclic_locants.size(); 
  ring->pseudo_points = pseudo_locants.size(); 
  
    for (unsigned int i=0;i<252;i++){
      if(bridge_locants[i])
        ring->bridge_points++; 
    }

  for (std::pair<unsigned int, unsigned char> comp : ring_components){    
      ring->assignment_locants.push_back(comp.second); 
      ring->assignment_digits.push_back(comp.first); 
    } 

  if (!final_size)
    return Fatal(start+i, "Error: failed to build WLN cycle unit");
 
#if BUILD_CYCLE
  if(!post_unsaturate(unsaturations,final_size,ring) || !post_saturate(saturations,final_size,ring))
    return Fatal(start+i, "Error: failed on post ring bond (un)/saturation");
#endif 

  return true;
}

bool multiply_carbon(WLNSymbol *sym){
  if(!sym->parr_n || !sym->barr_n)
    return false;
  
  WLNSymbol *back     = sym->prev_array[0].child; 
  WLNEdge   *fedge    = &sym->bond_array[0];
  
  WLNSymbol *forward  = fedge->child;
  WLNEdge *bedge = 0;
  for(unsigned int ei=0;ei<back->barr_n;ei++){      
    if(back->bond_array[ei].child == sym){
      bedge = &back->bond_array[ei];
      break;
    }
  }

  if(!forward||!bedge)
    return false;
  

  unsigned int back_edges = back->allowed_edges - back->num_edges; 
  unsigned int forward_edges = forward->allowed_edges - forward->num_edges;

  // it seems to be the convention
  // that if the forward facing symbol can take a triple bond, we do it
  
  // this cannot be done for alkyl numbers, so.. turn off the chain edges
  if(back->ch == '#')
    back_edges = 1;
  if(forward->ch == '#')
    forward_edges = 1;

  // experimental rule, if a triple bond will completely saturate an, atom, 
  // we should always take it. 

  if(forward->num_edges== 1 && forward->num_edges+2 == forward->allowed_edges){
    if(!unsaturate_edge(fedge,2))
      return false;
  }
  else if(back->num_edges == 1 && back->num_edges+2 == back->allowed_edges){
    if(!unsaturate_edge(bedge,2))
      return false;
  }
  else{
    // first cases, forward edges force the triple
    if(forward_edges >= 2){
      if(!unsaturate_edge(fedge,2))
        return false;
    }
    else if (forward_edges == 1 && back_edges >=1){
      if(!unsaturate_edge(bedge,1) || !unsaturate_edge(fedge,1))
        return false;
    }
  }

  return true; 
}

/* adds in assumed double bonding that some WLN forms take
- it leads to ambiguety on how some structures could be represented */
bool ResolveHangingBonds(WLNGraph &graph){
  for(unsigned int i=0;i<graph.symbol_count;i++){
    WLNSymbol *sym = graph.SYMBOLS[i];
    WLNEdge *edge = 0; 

    if( ( sym->ch == 'O'  ||
          sym->ch == 'N'  ||
          sym->ch ==  'P' || 
          sym->ch ==  'S') &&
          sym->num_edges == 1 && sym->charge == 0)
    {
      edge = &sym->bond_array[0]; 
      if(edge && edge->order == 1){
        while((sym->num_edges < sym->allowed_edges) && 
        (edge->child->num_edges < edge->child->allowed_edges)){
          if(!unsaturate_edge(edge,1))
            return false;
        }
      }
    } 
    else{
      for(unsigned int ei=0;ei<sym->barr_n;ei++){
        edge = &sym->bond_array[ei]; 
        if( (edge->child->ch == 'O' ||
            edge->child->ch ==  'P'  || 
            edge->child->ch ==  'N'  || 
            edge->child->ch ==  'S') &&
            edge->child->num_edges == 1 && edge->child->charge == 0)
        {
          while((sym->num_edges < sym->allowed_edges) && 
            (edge->child->num_edges < edge->child->allowed_edges)){
            if(!unsaturate_edge(edge,1))
              return false;
          }
        }
      }
    }
  }
  return true;
}



/* must be performed before sending to obabel graph*/
bool ExpandWLNSymbols(WLNGraph &graph, unsigned int len){

  unsigned int stop = graph.symbol_count;
  // dioxo must be handled before 
  for (unsigned int i=0;i<stop;i++){
    WLNSymbol *sym = graph.SYMBOLS[i];
    if(sym->ch == 'W' && !add_dioxo(sym,graph))
      return Fatal(len, "Error: failed on past handling of W dioxo symbol");

    // unsaturated carbons with C
    if(sym->ch == 'c'){
      sym->ch = 'C';
      if(!multiply_carbon(sym))
        return Fatal(len,"Error: failed on post handling of multiplier carbon");
    }
  }

  stop = graph.symbol_count;
  for (unsigned int i=0;i<stop;i++){
    WLNSymbol *sym = graph.SYMBOLS[i];

    switch(sym->ch){
      case 'Y':
      case 'X':
      case 'K':
        if(!resolve_methyls(sym,graph))
          return Fatal(len,"Error: failed on post handling of undefined methyl groups");
        break;

      case 'V':{
        sym->ch = 'C';
        sym->allowed_edges = 4;
        
        WLNSymbol *oxygen = AllocateWLNSymbol('O',graph);
        oxygen->allowed_edges = 2;
        
        if(!AddEdge(oxygen, sym))
          return Fatal(len,"Error: failed on post expansion on 'V' symbol");
      
        WLNEdge *e = &sym->bond_array[sym->barr_n-1]; 
        if(!unsaturate_edge(e,1))
          return Fatal(len,"Error: failed on post expansion on 'V' symbol");
        
        break;
      }
      
      default:
        break; // ignore
    }
  }

  return ResolveHangingBonds(graph); 
}


/* uses dfs style traversal to see what is reachble from a given start
 * will use both types of edges */
void Reachable(WLNSymbol *node, std::set<WLNSymbol*> &reachable){
  std::stack<WLNSymbol*> stack;
  std::map<WLNSymbol*,bool> seen;
  stack.push(node); 
  while(!stack.empty()){
    WLNSymbol *top = stack.top(); 
    stack.pop(); 
    seen[top] = true;
    reachable.insert(top); 

    for(unsigned int ei=0;ei<top->barr_n;ei++){
      if(!seen[top->bond_array[ei].child])
        stack.push(top->bond_array[ei].child); 
    }

    for(unsigned int ei=0;ei<top->parr_n;ei++){
      if(!seen[top->prev_array[ei].child])
        stack.push(top->prev_array[ei].child); 
    }
  }
}

/**********************************************************************
                        WLN Ring Kekulize
**********************************************************************/

bool IsBipartite(WLNRing *ring){

  WLNSymbol *top = ring->locants['A'];
  if(!top){
    fprintf(stderr,"Error: graph is empty\n");
    return false; 
  }

  std::deque<WLNSymbol*> queue; 
  std::map<WLNSymbol*,unsigned int> color; // 0 un assigned, 1 first color, 2 second color

  color[top] = 1;

  queue.push_back(top); 
  while(!queue.empty()){

    top = queue.back();
    queue.pop_back();

    for(unsigned int ei=0;ei<top->barr_n;ei++){
      WLNSymbol *child = top->bond_array[ei].child; 
      if(!ring->locants_ch[child])
        continue;

      if(!color[child]){
        // find all non-coloured and assign alternative colour
        if(color[top] == 1)
          color[child] = 2;
        else
          color[child] = 1;
        
        queue.push_front(child);
      }
      else if(color[child] == color[top])
        return false;
      else if(child == top){
        // self loops impossible
        return false; 
      }
    }
  }

  return true;
}

bool AdjMatrixBFS(WLNRing *ring, unsigned int src, unsigned int sink, int *path){
  
  bool *visited = (bool*)malloc(sizeof(bool) * ring->rsize); // square so doesnt matter rows vs cols
  for(unsigned int i=0;i<ring->rsize;i++)
    visited[i] = false;
  
  std::deque<unsigned int> queue; 

  unsigned int u = src; 
  path[src] = -1; // stop condition 

  queue.push_back(u);

  while(!queue.empty()){
    u = queue.back();
    queue.pop_front();
    visited[u] = true; 

    for(unsigned int v=0;v<ring->rsize;v++){

      if((u * ring->rsize + v) >= ring->rsize*ring->rsize){
        if(visited)
          free(visited);
        visited = 0;
        return false;
      }
        
      if(!visited[v] && ring->adj_matrix[u * ring->rsize + v] > 0){
        path[v] = u;
        if(v == sink){
          if(visited)
            free(visited);
          visited = 0;
          return true;
        }

        queue.push_front(v);
      }
    }
  }

  if(visited)
    free(visited);
  visited = 0;
  return false;
}

bool BPMatching(WLNRing *ring, unsigned int u, bool *seen, int *MatchR){
  for(unsigned int v=0;v < ring->rsize;v++){
    if((u * ring->rsize + v) >= ring->rsize*ring->rsize)
      return false;

    if(ring->adj_matrix[u * ring->rsize + v] > 0 && !seen[v]){
      seen[v] = true;

      if(MatchR[v] < 0 || BPMatching(ring,MatchR[v],seen,MatchR)){
        MatchR[v] = u;
        return true;
      }
    }
  }
  return false; 
}

bool WLNRingBPMaxMatching(WLNRing *ring, int *MatchR){
  bool  *seen = (bool*)malloc(sizeof(bool) * ring->rsize);
  for (unsigned int i=0;i<ring->rsize;i++)
    seen[i] = false;
  
  for(unsigned int u=0; u< ring->rsize;u++)
    BPMatching(ring,u,seen,MatchR);
  
  if(seen)
    free(seen);
  seen = 0;
  return true;
}


/**********************************************************************
                         High Level Parser Functions
**********************************************************************/


/* returns the head of the graph, parse all normal notation */
bool ParseWLNLazy(const char *wln_ptr, WLNGraph &graph) 
{
  // keep the memory alive
  if (OPT_DEBUG)
    fprintf(stderr, "Parsing WLN notation: %s\n",wln_ptr);

  ObjectStack branch_stack;   // access to both rings and symbols
  //branch_stack.reserve(512);  // reasonable size given
  
  std::map<WLNRing*,WLNSymbol*> lazy_ring_lookup; 

  WLNSymbol *curr       = 0;
  WLNSymbol *prev       = 0;
  WLNEdge   *edge       = 0;
  WLNRing   *ring       = 0;
  WLNRing   *wrap_ring  = 0;

  bool cleared = true; // resets on ionics
  bool pending_locant           = false;
  bool pending_J_closure        = false;
  bool pending_inline_ring      = false;
  bool pending_spiro            = false;
  bool pending_ring_in_ring     = false; // rings in rings
  bool pending_rir_closure      = false;
  bool pending_negative_charge  = false; // lets get rid of a lot of waste
  bool pending_carbon_chain     = false;

  bool no_shift = false; // stop shifting if already done
  std::string str_buffer; 
  std::string digits_buffer; 

  unsigned char on_locant = '\0';         // locant tracking
  unsigned int pending_unsaturate = 0;    // 'U' style bonding
  bool j_skips = false;                   // handle skipping of 'J' if in cyclic notation legitimately 
  
  // allows consumption of notation after block parses
  unsigned int block_start = 0;
  unsigned int block_end = 0;

  unsigned int len = strlen(wln_ptr);

  unsigned int i=0;
  unsigned char ch = *wln_ptr;
  
  while(ch)
  {  

    // this will need to resolved at the end as well
    if(pending_carbon_chain && (ch < '0' || ch > '9') && ch != '/' ){
      if(digits_buffer.empty() || digits_buffer[0] == '0')
        return Fatal(i,"Error: zero mark opened for chain length, invalid numeral"); 

      int carbon_len = isNumber(digits_buffer);
      digits_buffer.clear(); 
      if(carbon_len < 0)
        return Fatal(i, "Error: non-numeric value entered for carbon length");
      else if (carbon_len > 100)
        return Fatal(i,"Error: creating a carbon chain > 100 long, is this reasonable for WLN?");

      // create the head 
      WLNSymbol *curr = AllocateWLNSymbol('#', graph);
      curr->str_position = i;
      curr->special = std::to_string(carbon_len); 
      curr->allowed_edges = 4; 

      if(prev){
        if(prev == branch_stack.branch){
          while(!branch_stack.top().second && !branch_stack.empty())
            branch_stack.pop();
        }
        
        if(!AddEdge(curr, prev))
          return Fatal(i, "Error: failed to bond to previous symbol");

        edge = &prev->bond_array[prev->barr_n-1]; 
        if(pending_unsaturate){
          if(!unsaturate_edge(edge,pending_unsaturate))
            return Fatal(i, "Error: failed to unsaturate bond"); 
          pending_unsaturate = 0;
        }
      }

      pending_carbon_chain = false;
      prev = curr; 
      cleared = false;
    }

    switch (ch)
    {

    case '0': // cannot be lone, must be an addition to another num
      if(pending_J_closure)
        break;

      else if (pending_locant){
        
        if(pending_inline_ring && prev && !prev->inRing)
          prev->charge++;

        prev = 0;
        on_locant = '0';
        pending_locant = false;
      }
      else if (pending_carbon_chain){
        digits_buffer+= ch;
        break;
      }
      else if(pending_negative_charge){
        digits_buffer += ch;
      }
      else if (cleared){
        // null effect for positive charge
        digits_buffer += ch;
        pending_carbon_chain = true; 
      }
      else
        return Fatal(i,"Error: a lone zero mark is not allowed without positive numerals either side");
      
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
      if (pending_J_closure){
        // small addition to allow J handling in points
        if(i > 0 && wln_string[i-1] == ' ')
          j_skips = true;
        
        break;
      }
        
      else if(pending_locant){  // handle all multiplier contractions
        return Fatal(i,"Error: multipliers are not currently supported");
        pending_locant = false;
        on_locant = ch;
      }
      else if(pending_ring_in_ring && pending_inline_ring){
        // onlocant holds the char needed to wrap the ring back, 
        
        if(on_locant != '0'){
          curr = wrap_ring->locants[on_locant];
          if(!curr)
            return Fatal(i,"Error: cannot access looping ring structure");
          
          if(prev){

            if(prev == branch_stack.branch){
              while(!branch_stack.top().second && !branch_stack.empty())
                branch_stack.pop();
            }
        
            if(!AddEdge(curr, prev))
              return Fatal(i, "Error: failed to bond to previous symbol");
            
            edge = &prev->bond_array[prev->barr_n-1]; 
            if(pending_unsaturate){
              if(!unsaturate_edge(edge,pending_unsaturate))
                return Fatal(i, "Error: failed to unsaturate bond"); 
              pending_unsaturate = 0;
            }

          }
          else
            return Fatal(i,"Error: no previous symbol for inline ring defintion");

          on_locant = '\0';
        }
        
        // last notation is not neccessary
        while(wln_ptr){
          if(*wln_ptr == 'J')
            break;
          wln_ptr++;
          i++;
        }

        pending_ring_in_ring = false;
        pending_inline_ring = false;
        cleared = false;
      }
      else if (pending_carbon_chain){
        digits_buffer += ch;  // simple resolve
        break;
      }
      else if (pending_negative_charge){
        digits_buffer += ch; 
      }
      else{
        on_locant = '\0';
        pending_carbon_chain = true;
        digits_buffer += ch; 
      }
      break;
    

    case 'Y':
      if (pending_J_closure)
        break;
      else if (pending_locant)
        return Fatal(i,"Error: 'Y' cannot be a locant assignment, please expand [A-W] with &\n");
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 4; // change methyl addition

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }
          
          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }
        
        branch_stack.push({0,curr});
        pending_unsaturate = 0;
        prev = curr;
      }
      cleared = false;
      break;

    case 'X':
      if (pending_J_closure)
        break;
      else if (pending_locant){
        return Fatal(i, "Error: Wiswesser Uncertainities lead to runaway outcomings");
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 4;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }
        
        branch_stack.push({0,curr});
        prev = curr;
      }
      cleared = false;
      break;

      // oxygens

    case 'O':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 2;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }

        prev = curr;
      }
      cleared = false;
      break;

    case 'Q':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 1;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }

        pending_unsaturate = 0;
        prev = return_object_symbol(branch_stack);

        if(!prev)
          prev = curr;
      }
      cleared = false;
      break;

    case 'V':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 2;
        
        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          } 
        }

        prev = curr;
      }
      cleared = false;
      break;

    case 'W':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 3;
        curr->str_position = i+1;

        if(prev){

          if(prev->ch == 'N')
            prev->allowed_edges++;

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(!unsaturate_edge(edge, 2))
            return Fatal(i, "Error: failed to attach W symbol"); 

          if(pending_unsaturate)
            return Fatal(i,"Error: a bond unsaturation followed by dioxo is undefined notation");
        }
        else
          pending_unsaturate = 2;

        // if there is no prev, this is then a character we go from, 
        // else do not update.
        if(!prev)
          prev = curr;
        else
          prev = return_object_symbol(branch_stack);

      }
      cleared = false;
      break;

      // nitrogens

    case 'N':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 3;

        if(prev){
          if(prev->ch == 'W')
            curr->allowed_edges++;

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }
          
          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }

        branch_stack.push({0,curr});
        pending_unsaturate = 0;
        prev = curr;
      }
      cleared = false;
      break;

    case 'M':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 2;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){ 
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }

        pending_unsaturate = 0;
        prev = curr;
      }
      cleared = false;
      break;

    case 'K':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 4;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
          
        }
        
        branch_stack.push({0,curr});
        prev = curr;
      }
      cleared = false;
      break;

    case 'Z':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];          
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
      
        pending_locant = false;
        on_locant = ch;
      }
      else
      { 
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 1;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
          
        }

        pending_unsaturate = 0;
        prev = return_object_symbol(branch_stack);
        if(!prev)
          prev = curr;
      }
      cleared = false;
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
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 1;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }

        pending_unsaturate = 0;
        prev = return_object_symbol(branch_stack);
        if(!prev)
          prev = curr;
      }
      cleared = false;
      break;

      // inorganics

    case 'B':
      if (pending_J_closure)  
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }

          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      { 
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 3;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }         
        
        branch_stack.push({0,curr});
        prev = curr;
      }
      cleared = false;
      break;

    case 'P':
    case 'S':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        
        if(ch == 'P')
          curr->allowed_edges = 5;
        else
          curr->allowed_edges = 6;
          
        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }
        branch_stack.push({0,curr});
        prev = curr;
      }
      cleared = false;
      break;

    // multiply bonded carbon, therefore must be at least a double bond
    case 'C':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        // set lower case for multiplier carbon
        curr = AllocateWLNSymbol('c',graph);
        curr->str_position = i+1;
        curr->allowed_edges = 4;

        if(prev && i < len - 1){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }

        prev = curr;
      }
      cleared = false;
      break;

    case 'A':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
        return Fatal(i,"Error: locant only symbol used in atomic definition");
      
      cleared = false;
      break;
        
    // this can start a chelating ring compound, so has the same block as 'L\T'
    case 'D':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {

        if(i < len - 2 && wln_string[i+1] == '-' && (wln_string[i+2] == 'T' || wln_string[i+2] == 'L')){
          pending_ring_in_ring = true;

          i++;
          wln_ptr++;
          pending_inline_ring = true;
          break;
        }
          
        if (i == 0)
          pending_inline_ring = true;

        if (!pending_inline_ring)
          return Fatal(i,"Error: chelating ring notation started without '-' denotion");
        
        
        pending_inline_ring = false;
        block_start = i;
        pending_J_closure = true;
      }
      cleared = false;
      break;
        
        
    // hydrogens explicit

    case 'H':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else{
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 1;
  
        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }

          switch(prev->ch){
            case 'Z':
                //graph.charge_additions[prev]++;
                prev->allowed_edges++;
                break;
            default:
              break;
          }
        }

        if(prev && (prev->ch == 'V' || prev->ch == 'M')){
          curr = prev;
        }
        else
          prev = return_object_symbol(branch_stack);
        
        if(!prev)
          prev = curr; // failsafe to starting hydrogen
      }
      cleared = false;
      break;

      // ring notation

    case 'J':
      if(pending_rir_closure){
        wrap_ring = 0;
        pending_rir_closure = false;
        break;
      }
      if(pending_J_closure && j_skips)
        break;
      if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else if (pending_J_closure 
              && ( (i<len-1 && (wln_string[i+1] == ' ' || wln_string[i+1] == '&') && wln_string[i-1] != ' ') 
              || i == len -1)
              )     
      {
        block_end = i;
        
        ring = AllocateWLNRing(graph);
        std::string r_notation = get_notation(block_start,block_end);

        if(pending_spiro){
          
          if(!prev)
            Fatal(i,"Error: sprio notation opened without a previous atom");
          else
            ring->locants[on_locant] = prev;
          
          if(!FormWLNRingLazy(ring,r_notation,block_start,graph,on_locant)){
            return Fatal(len,"Error: failed on ring evaluation"); 
          }

          curr = AllocateWLNSymbol('&', graph); 
          curr->inRing = ring; 
          curr->allowed_edges = 2; // default spiro locant
          AddEdge(curr, prev);
        }
        else{
          if(!FormWLNRingLazy(ring,r_notation,block_start,graph))
            return Fatal(len,"Error: failed on ring evaluation"); 
        }

        WLNSymbol *ring_symbol = AllocateWLNSymbol('*',graph); 
        ring_symbol->inRing = ring;
        ring_symbol->special = ring->str_notation;
        ring_symbol->allowed_edges = ring->rsize*2; // allow each locant to have 2 bonds for now 
        lazy_ring_lookup[ring] = ring_symbol; 

        if(pending_ring_in_ring && !wrap_ring)
          wrap_ring = ring; // instant back access

        branch_stack.push({ring,0});
        block_start = 0;
        block_end = 0;

        // does the incoming locant check
        if(pending_spiro){  
          AddEdge(curr, lazy_ring_lookup[ring]);
          pending_spiro = false;
        }
        else if (prev && on_locant && on_locant != '0')
        {
          
          if(!ring->locants[on_locant]){
            ring->locants[on_locant] = AllocateWLNSymbol(on_locant, graph);
            curr = ring->locants[on_locant]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }

          if(!AddEdge(ring->locants[on_locant], prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
          ring->loc_count++; //in-line locant assumed
        }   
          
        on_locant = '\0';
        pending_J_closure = false;
      }
    //  else
     //   return Fatal(i, "Error: J character in unrecognised state"); 
      
      cleared = false;
      break;

    case 'L':
    case 'T':
      if (pending_J_closure || pending_rir_closure)
        break;
   //   else if (pending_ring_in_ring)
   //     break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        if(i < len - 2 && wln_string[i+1] == '-' && (wln_string[i+2] == 'T' || wln_string[i+2] == 'L')){
          pending_ring_in_ring = true;

          // if pending inline ring doesnt get a L or T, we know it
          // the ring wrap symbol. 

          i++;
          wln_ptr++;
          pending_inline_ring = true;
          break;
        }
          
        if (cleared)
          pending_inline_ring = true;
          
        if (!pending_inline_ring)
          return Fatal(i,"Error: ring notation started without '-' denotion");
        
        
        pending_inline_ring = false;
        block_start = i;
        pending_J_closure = true;
      }
      cleared = false;
      break;

    case 'R':
      if (pending_J_closure)
        break;
      
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        ring = AllocateWLNRing(graph);

        std::string r_notation = "L6J";
        FormWLNRingLazy(ring,r_notation,i,graph);
        branch_stack.push({ring,0});

        curr = ring->locants['A'];
        if(prev){
          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          if(pending_unsaturate){
            if(!unsaturate_edge(edge,pending_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            pending_unsaturate = 0;
          }
        }

        prev = curr;
      }
      cleared = false;
      break;

      // bonding

    case 'U':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            ring->locants[ch] = AllocateWLNSymbol(ch, graph);
            curr = ring->locants[ch]; 
            curr->inRing = ring; 
            curr->allowed_edges = 3; // default locant
            AddEdge(curr, lazy_ring_lookup[ring]); 
          }
          ring->loc_count++;  
          prev = curr;
        }
        
        on_locant = ch;
      }
      else if(cleared)
        return Fatal(i, "Error: floating double bond after ionic clear");
      
      else{
        on_locant = '\0';
        pending_unsaturate++;
      }
      break;

      // specials

    case ' ':
      if (pending_J_closure){
        j_skips = false;
        break;
      }
      else if(pending_negative_charge){
        int negative_index = isNumber(digits_buffer); 
        digits_buffer.clear(); 

        if(negative_index < 0)
          return Fatal(i, "Error: assigning non-numerical value to charge index");
        else if(negative_index != 0){
          // find the symbol and increment its charge + 1
          bool found = false;
          for(unsigned int cs = 0;cs<graph.symbol_count;cs++){
            if(graph.SYMBOLS[cs]->str_position == (unsigned int)negative_index){
              graph.SYMBOLS[cs]->charge--; 
#if OPT_DEBUG
              fprintf(stderr,"assigning %c charge %d\n",graph.SYMBOLS[cs]->ch,graph.SYMBOLS[cs]->charge); 
#endif
              found = true;
              break;
            }
          }
          
          if(!found)
            return Fatal(i, "Error: negative charge index out of range, check letter index");
        }
        pending_negative_charge = false;
      }

      if(!branch_stack.empty() && !pending_inline_ring)
        branch_stack.pop_to_ring();
      
      if( (i < len - 1 && wln_string[i+1] == '&') || branch_stack.ring){
        pending_locant = true;

        // single letter methyl branches
        if(on_locant && !pending_inline_ring){
          if(!branch_stack.ring || !add_methyl(branch_stack.ring->locants[on_locant],graph))
            return Fatal(i, "Error: could not attach implied methyl to ring");
          
          branch_stack.ring->loc_count++;  
          on_locant = '\0';
        }
        
      }
      else
        return Fatal(i, "Error: space used outside ring and ionic notation");
      
      // only burn the stacks now on ionic clearance
      break;


    case '&':
      if (pending_J_closure)
        break;
      
      if (pending_inline_ring && !pending_ring_in_ring)
      {
        // spiro notation open
        pending_spiro = true;
      }
      else if (pending_locant)
      { 
        // ionic species or spiro, reset the linkings
        prev = 0;
        curr = 0;
        ring = 0;
        pending_locant = false;
        cleared = true;
        branch_stack.clear_all(); // burn stack
      }
      else if(on_locant){
        if(curr && curr == ring->locants[on_locant]){
          on_locant += 23;
          curr = ring->locants[on_locant];  
          if(!curr)
            return Fatal(i, "Error: could not fetch expanded locant position - out of range");
          
          prev = curr;
        }
      }
#if DEPRECATED
      else if (i < len-1 && wln_string[i+1] == ' '){
        // this must be a ring pop, no matter what
        // technically not true!
        
        if(branch_stack.empty() || !branch_stack.ring)
          return Fatal(i, "Error: '&' followed by a space indicates a ring pop, are there any rings?"); 
        
        else if(!branch_stack.empty()){
          branch_stack.pop_to_ring();
          branch_stack.pop(); // pop whats open
          ring = branch_stack.ring; // assign the previous ring
          
          prev = return_object_symbol(branch_stack);
          if(!prev)
            prev = branch_stack.branch;
        }
      }
#endif
      else if(!branch_stack.empty())
      {

        if(branch_stack.top().first){
          branch_stack.pop();
          prev = return_object_symbol(branch_stack);
          if(!prev)
            prev = branch_stack.branch;
          ring = branch_stack.ring;
        }

        else if(branch_stack.top().second){
          WLNSymbol *top = 0;
          top = branch_stack.top().second;

            // this means a <Y|X|..>'&' so handle methyl
          if(prev && prev == top){
            switch(prev->ch){
              // methyl contractions

              case 'Y':
                if(count_children(prev) < 3){
                  if(!add_methyl(prev,graph))
                    return Fatal(i, "Error: failed to add methyl group on methyl contraction");

                  prev = return_object_symbol(branch_stack);
                }
                else{ 
                  // we pop,
                  branch_stack.pop();
                  prev = branch_stack.branch; 
                }
                break;

              case 'X':
              case 'K':
                if(prev->num_edges < prev->allowed_edges){
                  if(!add_methyl(prev,graph))
                    return Fatal(i, "Error: failed to add methyl group on methyl contraction");

                  prev = return_object_symbol(branch_stack);
                }
                else{ 
                  // we pop, 
                  branch_stack.pop();
                  prev = branch_stack.branch; 
                }
                break;

              // no contractions possible, we pop the stack
              default:
                // default pop
                branch_stack.pop();
                prev = return_object_symbol(branch_stack); // if prev is nulled, then a ring is active
                if(!prev)
                  prev = branch_stack.branch;

                break;
            }
          }
          else{
            // means a closure is done, we return to the first avaliable symbol on the branch stack
            prev = return_object_symbol(branch_stack);
            if(branch_stack.top().first)
              branch_stack.pop();
            
              
            if(!prev)
              prev = branch_stack.branch; // catches branching ring closure
          }
        }
      }
      else
        return Fatal(i, "Error: popping too many rings|symbols, check '&' count");
        
      break;


    case '-':{
      if (pending_J_closure)
        break;

      else if (pending_inline_ring)
      { 
        
        if(pending_ring_in_ring){
          // onlocant holds the char needed to wrap the ring back, 
          
          if(!wrap_ring)
            return Fatal(i, "Error: wrap ring is not active");

          curr = wrap_ring->locants[on_locant];
          if(!curr)
            return Fatal(i, "Error: cannot access looping ring structure");
        
          if(prev){  
            if(prev == branch_stack.branch){
              while(!branch_stack.top().second && !branch_stack.empty())
                branch_stack.pop();
            }

            if(!AddEdge(curr, prev))
              return Fatal(i, "Error: failed to bond to previous symbol");

            edge = &prev->bond_array[prev->barr_n-1]; 
            if(pending_unsaturate){
              if(!unsaturate_edge(edge,pending_unsaturate))
                return Fatal(i, "Error: failed to unsaturate bond"); 
              pending_unsaturate = 0;
            }
          }
          else
            return Fatal(i,"Error: no previous symbol for inline ring defintion");
      

          // if we remove this, we can allow branching inline defintiions 
          // last notation is not neccessary, so we eat two positions
          unsigned int hit = 0;
          while(wln_ptr){
            if(*wln_ptr == 'J')
              return Fatal(i, "Error: macro-notation requires closure with the ring size in two dashes e.g -6-");
            
            if(*wln_ptr == '-'){
              hit++;
              if(hit == 2)
                break;
            }
            wln_ptr++;
            i++;
          }
          
          curr = prev; // set back to prev
          on_locant = '\0';
          pending_ring_in_ring = false;
          pending_inline_ring = false;
          pending_rir_closure = true;
        }
        else
          return Fatal(i, "Error: only one pending ring can be active, check closures");
        
      }
      else{

        str_buffer.clear();
        bool closed_dash = false;
        unsigned int first_dash = i;

        // move off the first dash
        i++;
        ch = *(++wln_ptr);
        while(ch != '\0'){
          if(ch == '-'){
            closed_dash = true;
            break;
          }
          else if (ch == ' '){
            no_shift = true;
            break; // positions an open ring
          }
          else if (ch == '&'){
            pending_spiro = true;
            break;
          }
          else
            str_buffer.push_back(ch);
          
          i++;
          ch = *(++wln_ptr);
        }

        if(!closed_dash){
          pending_inline_ring = true;
          return_object_symbol(branch_stack);
          if(branch_stack.branch && !prev){
            // prev must be at top of the branch stack
            while(branch_stack.top().second != branch_stack.branch)
              branch_stack.pop();

            prev = return_object_symbol(branch_stack);
          }    
        }
        else{

          if(str_buffer.size() == 1){
            curr = define_hypervalent_element(str_buffer[0],graph);
            if(!curr)
              return Fatal(i, "Error: failed to define hypervalent element");
            curr->str_position = i+1;
          }
          else if(str_buffer.size() == 2){
            curr = define_element(str_buffer,graph);
            if(!curr)
              return Fatal(i, "Error: failed to define periodic element"); 
            
            curr->str_position = i+1;

            if(on_locant == '0'){
              curr->charge++;
              //on_locant = '\0';
            }
          }
          else
            return Fatal(i, "Error: special '-' must be either 1 or 2 symbols");
          
        
          if(prev){
            if(str_buffer.empty() && ring){
              if(!AddEdge(ring->locants[prev->ch], prev))
                return Fatal(i, "Error: failed to bond to previous symbol");
              edge = &prev->bond_array[prev->barr_n-1]; 
            }
            else{
              if(prev == branch_stack.branch){
                while(!branch_stack.top().second && !branch_stack.empty())
                  branch_stack.pop();
              }
              if(!AddEdge(curr, prev))
                return Fatal(i, "Error: failed to bond to previous symbol");
              edge = &prev->bond_array[prev->barr_n-1]; 
            }
              
            if(!edge)
              return Fatal(i, "Error: failed to bond to previous symbol");
            else if(pending_unsaturate){
              if(!unsaturate_edge(edge,pending_unsaturate))
                return Fatal(i, "Error: failed to unsaturate bond"); 
              pending_unsaturate = 0;
            }
            
          }
          on_locant = '\0';
          branch_stack.push({0,curr});
          
          curr->str_position = first_dash+1+1;
          pending_unsaturate = 0;
          prev = curr;
        }
      }
      cleared = false;
      break;
    }
    
    case '/':
      if (pending_J_closure){
        j_skips = true;
        break;
      }
      else if(pending_carbon_chain){ // state that this must be a charge 

        if(!cleared)
          return Fatal(i, "Error: opening post charge assignment without proper syntax [ &x/x ]");
        
        int positive_index = isNumber(digits_buffer); 
        digits_buffer.clear(); 
        
        if(positive_index < 0)
          return Fatal(i, "Error: assigning non-numerical value to charge index");
        else if (positive_index != 0){
          // find the symbol and increment its charge + 1
          bool found = false;
          for(unsigned int cs = 0;cs<graph.symbol_count;cs++){
            if(graph.SYMBOLS[cs]->str_position == (unsigned int)positive_index){
              graph.SYMBOLS[cs]->charge++;
#if OPT_DEBUG
              fprintf(stderr,"assigning %c charge %d\n",graph.SYMBOLS[cs]->ch,graph.SYMBOLS[cs]->charge); 
#endif
              found = true;
              break;
            }
          }
          
          if(!found)
            return Fatal(i, "Error: positive charge index out of range, check letter index");
        }

        pending_carbon_chain = false;
        pending_negative_charge = true;
      }
      else
        return Fatal(i,"Error: multipliers are not currently supported");
      cleared = false;
      break;

    default:
      return Fatal(i,"Error: unallowed character! - alphabet: [A-Z][0-1][&-/' ']");
    }

    if(!no_shift){
      i++;
      ch = *(++wln_ptr);
    }
    no_shift = false;
  }


  if(pending_carbon_chain){        
    if(digits_buffer.empty() || digits_buffer[0] == '0')
      return Fatal(i,"Error: zero mark opened for chain length, invalid numeral"); 
    int carbon_len = isNumber(digits_buffer);

    digits_buffer.clear(); 
    if(carbon_len < 0)
      return Fatal(i, "Error: non-numeric value entered for carbon length");
    else if (carbon_len > 100)
      return Fatal(i,"Error: creating a carbon chain > 100 long, is this reasonable for WLN?");
    

    // create the head 
    WLNSymbol *curr = AllocateWLNSymbol('#', graph);
    curr->str_position = i;
    curr->special = std::to_string(carbon_len); 
    curr->allowed_edges = 4; 
    
    if(prev){
      if(prev == branch_stack.branch){
        while(!branch_stack.top().second && !branch_stack.empty())
          branch_stack.pop();
      }

      if(!AddEdge(curr, prev))
        return Fatal(i, "Error: failed to bond to previous symbol");
      edge = &prev->bond_array[prev->barr_n-1]; 
      
      if(pending_unsaturate){
        if(!unsaturate_edge(edge,pending_unsaturate))
          return Fatal(i, "Error: failed to unsaturate bond"); 
        pending_unsaturate = 0;
      }
    }

    pending_carbon_chain = false;
    prev = curr; 
  }

  // single letter methyl branches
  if(on_locant && on_locant != '0' && !pending_inline_ring && !branch_stack.empty()){
    if(!add_methyl(branch_stack.ring->locants[on_locant],graph))
      return Fatal(i, "Error: could not attach implied methyl to ring");
    
    on_locant = 0;
  }
  
  if(pending_negative_charge){
    int negative_index = isNumber(digits_buffer); 
    digits_buffer.clear();

    if(negative_index < 0)
      return Fatal(i, "Error: assigning non-numerical value to charge index");
    else if (negative_index != 0){
      // find the symbol and increment its charge + 1
      bool found = false;
      for(unsigned int cs = 0;cs<graph.symbol_count;cs++){
        if(graph.SYMBOLS[cs]->str_position == (unsigned int)negative_index){
          graph.SYMBOLS[cs]->charge--; 
#if OPT_DEBUG
          fprintf(stderr,"assigning %c charge %d\n",graph.SYMBOLS[cs]->ch,graph.SYMBOLS[cs]->charge); 
#endif
          found = true;
          break;
        }
      }
      
      if(!found)
        return Fatal(i, "Error: negative charge index out of range, check letter index");
    }
  }

  if (pending_J_closure)
    return Fatal(len, "Error: ring open at end of notation, inproper closure");
  

  if (pending_locant)
    return Fatal(len, "Error: locant open at end of notation, inproper closure");
  

  if (pending_inline_ring)
    return Fatal(len, "Error: inline ring expected at end of notation, inproper closure");
  

  if (pending_spiro)
    return Fatal(len, "Error: spiro ring expected at end of notation, inproper closure");

    // use this for recursion on multipliers
  return true;
}


/* dump wln tree to a dotvis file */
void WLNDumpToDot(FILE *fp, WLNGraph &graph)
{  
  fprintf(fp, "digraph WLNdigraph {\n");
  for (unsigned int i=0; i< graph.symbol_count;i++)
  {
    WLNSymbol *node = graph.SYMBOLS[i];
    if(!node)
      continue;

    fprintf(fp, "  %d", node->id);
    if (node->ch == '*' || node->ch == '#')
      fprintf(fp, "[shape=circle,label=\"* %s\"];\n", node->special.c_str());
    else if (node->inRing)
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
  
      
    for (unsigned int ei=0;ei<node->barr_n;ei++){
      WLNEdge *edge = &node->bond_array[ei];
      WLNSymbol *child = edge->child;
      unsigned int bond_order = edge->order;

      if (bond_order > 1){
        for (unsigned int k=0;k<bond_order;k++){
          fprintf(fp, "  %d", node->id);
          fprintf(fp, " -> ");

          if(edge->aromatic)
            fprintf(fp, "%d [color=red]\n", child->id);
          else
            fprintf(fp, "%d\n", child->id);
        }
      }
      else{
        fprintf(fp, "  %d", node->id);
        fprintf(fp, " -> ");
        if(edge->aromatic)
          fprintf(fp, "%d [color=red]\n", child->id);
        else
          fprintf(fp, "%d\n", child->id);
      }
    }
  }

  fprintf(fp, "}\n");
}

bool WriteGraph(WLNGraph &graph,const char*filename){
  fprintf(stderr,"Dumping wln tree to %s:\n",filename);
  FILE *fp = 0;
  fp = fopen(filename, "w");
  if (!fp)
  {
    fprintf(stderr, "Error: could not create dump .dot file\n");
    fclose(fp);
    return false;
  }
  else
    WLNDumpToDot(fp,graph);
  
  fclose(fp);
  fp = 0;
  fprintf(stderr,"  dumped\n");
  return true;
}


static void DisplayUsage()
{
  fprintf(stderr, "wlntree <options> <input (escaped)>\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser feature trees ---\n\n");
  fprintf(stderr, " This parser reads and evaluates wiswesser\n"
                  " line notation (wln) and builds feature tree\n"
                  " objects for graph edit search, the parser is native C\n"
                  " with a plug in function to OpenBabel\n"
        );
  DisplayUsage();
}


static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i;
  unsigned int j = 0;

  cli_inp = (const char *)0;

  if (argc < 2)
    DisplayUsage();

  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];

    if (ptr[0] == '-' && ptr[1]){

      if(ptr[1] >= 'A' && ptr[1] <= 'Z' && !j){
        cli_inp = ptr;
        j++; 
      }
    
      else{
        switch (ptr[1])
        {
        case 'h':
          DisplayHelp();

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
        }
      }
    }
    else{
      switch (j)
      {
      case 0:
        cli_inp = ptr;
        break;
      
      default:
        fprintf(stderr,"Error: wln string already set - %s\n",cli_inp);
        DisplayUsage();
      }
      j++;
    }
  }

  if(!cli_inp){
    fprintf(stderr,"Error: no input string entered\n");
    DisplayUsage();
  }
  else 
    wln_string = ptr;

  return;
}

int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);

  WLNGraph graph; 
  if(!ParseWLNLazy(cli_inp, graph)){
    fprintf(stderr,"Error: failed building feature tree\n"); 
    return 1;
  }
  WriteGraph(graph, "wln-tree.dot"); 
  return 0;
}
