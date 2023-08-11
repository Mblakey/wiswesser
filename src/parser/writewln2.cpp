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

#define REASONABLE 1024

const char *cli_inp;
const char *format; 

// --- options ---
static bool opt_wln2dot = false;
static bool opt_debug = false;


struct WLNSymbol;
struct WLNEdge; 
struct WLNRing;
struct WLNGraph;
struct ObjectStack;


enum WLNTYPE
{
  STANDARD = 0,
  RING = 1,     
  SPECIAL = 2  // for now this is only going to be the pi bond
};

unsigned char static int_to_locant(unsigned int i){
  return i + 64;
}

unsigned int static locant_to_int(unsigned char loc){
  return loc - 64;
}



/**********************************************************************
                          STRUCT DEFINTIONS
**********************************************************************/
 

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
  ~WLNEdge(){};
};


struct WLNSymbol
{
  unsigned char ch;
  std::string special; // string for element, or ring, if value = '*'
  
  unsigned int type;
  unsigned int allowed_edges;
  unsigned int num_edges;

  unsigned int num_children; // specifically for forward edges

  WLNSymbol *previous;
  WLNEdge   *bonds; // array of bonds

  // if default needed
  WLNSymbol()
  {
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
    type = 0;
    previous = 0;
    bonds = 0;
    num_children = 0; 
  }
  ~WLNSymbol(){};

  void set_edge_and_type(unsigned int e, unsigned int t=STANDARD){
    allowed_edges = e;
    type = t;
  }

};

struct WLNRing
{
  std::vector<unsigned int> rings;
  std::map<unsigned char, WLNSymbol *> locants; 
  std::map<WLNSymbol*,unsigned char> locants_ch;
  std::vector<std::pair<unsigned char,int>> post_charges; 
  
  WLNRing(){}
  ~WLNRing(){};
};


// handles all memory and 'global' vars
struct WLNGraph
{
  
  WLNSymbol *root;

  unsigned int edge_count;
  unsigned int symbol_count;
  unsigned int ring_count;

  WLNSymbol *SYMBOLS[REASONABLE];
  WLNEdge   *EDGES  [REASONABLE];
  WLNRing   *RINGS  [REASONABLE];

  std::map<WLNSymbol *, unsigned int> index_lookup;
  std::map<unsigned int, WLNSymbol *> symbol_lookup;

  unsigned int glob_index; // babel starts from 1, keep consistent  

    // ionic parsing
  std::map<unsigned int,WLNSymbol*> string_positions; 
  std::map<WLNSymbol*,int> charge_additions;

  WLNGraph(){
    edge_count   = 0;
    symbol_count = 0;
    ring_count   = 0;
    glob_index   = 1; // babel atoms are +1 indexed

    // pointer safety
    root = 0;
    for (unsigned int i = 0; i < REASONABLE;i++){
      SYMBOLS[i] = 0;
      EDGES[i] = 0;
      RINGS[i] = 0;
    }
  };

  ~WLNGraph(){
    for (unsigned int i = 0; i < REASONABLE;i++){
      delete SYMBOLS[i];
      delete EDGES[i];
      delete RINGS[i];
    }
  }
};

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
      fprintf(stderr,"Error: peeking empty ring stack\n");
      return false;
    }
    else{
      fprintf(stderr,"top: ring: %p   branch: %p\n",stack.back().first,stack.back().second);
      return true;
    }
     
  }

  bool pop(){
    stack.pop_back();
    size--;

    ring = 0;
    branch = 0;

    if(stack.empty()){
      fprintf(stderr,"Error: popping empty ring stack\n");
      return false;
    }
     
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


  bool empty(){
    if (stack.empty())
      return true;
    else 
      return false;
  }

  void clear_all(){
    ring = 0;
    branch = 0;
    stack.clear();
  }

  std::pair<WLNRing*,WLNSymbol*> & top(){
    return stack.back();
  }

  // cleans branches
  bool branch_avaliable(){
    if(branch && branch->num_edges < branch->allowed_edges)
      return true;
    else
      return false;
  }

};



/**********************************************************************
                         WLNSYMBOL Functions
**********************************************************************/


WLNSymbol *AllocateWLNSymbol(unsigned char ch, WLNGraph &graph)
{

  graph.symbol_count++;
  if(graph.symbol_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
    return 0;
  }

  if(!ch){
    fprintf(stderr,"Error: null char used to symbol creation\n");
    return 0;
  }

  WLNSymbol *wln = new WLNSymbol;
  graph.SYMBOLS[graph.symbol_count] = wln;
  
  wln->ch = ch;
  graph.index_lookup[wln] = graph.glob_index;
  graph.symbol_lookup[graph.glob_index] = wln;
  graph.glob_index++;
  return wln;
}

WLNSymbol* define_hypervalent_element(unsigned char sym, WLNGraph &graph){

  if(!sym){
    fprintf(stderr,"Error: null char used for hypervalent element allocation\n");
    return 0;
  }

  WLNSymbol *new_symbol = 0;
  
  switch(sym){
    
    case 'P':
    case 'S':
    case 'G':
    case 'E':
    case 'I':
    case 'F':
      new_symbol = AllocateWLNSymbol(sym,graph);
      if(new_symbol)
        new_symbol->set_edge_and_type(6);            // allows FCl6
      break;

    default:
      fprintf(stderr,"Error: character %c does not need - notation for valence expansion, please remove -\n",sym);
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    case 'U':
      if(special[1] == 'R')
        created_wln = AllocateWLNSymbol('*',graph);
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;

    case 'V':
      if (special[1] == 'A')
        created_wln = AllocateWLNSymbol('*',graph);
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;
    
    case 'W':
      if(special[1] == 'T')
        created_wln = AllocateWLNSymbol('*',graph);
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
        return (WLNSymbol *)0;
      }
      break;
    

    case 'X':
      if (special[1] == 'E')
        created_wln = AllocateWLNSymbol('*',graph);
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
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
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
          fprintf(stderr, "Error: invalid element symbol in special definition - %s\n",special.c_str());
          return (WLNSymbol *)0;
      }
      break;

    default:
      fprintf(stderr, "Error: invalid character in special definition switch\n");
      return (WLNSymbol *)0;
  }

  created_wln->special = special;
  created_wln->allowed_edges = 8; // allow anything for now;
  return created_wln;
}



/**********************************************************************
                          WLNEdge Functions
**********************************************************************/


WLNEdge *AllocateWLNEdge(WLNSymbol *child, WLNSymbol *parent,WLNGraph &graph){

  if(!child || !parent){
    fprintf(stderr,"Error: attempting bond of non-existent symbols - %s|%s is dead\n",child ? "":"child",parent ? "":"parent");
    return 0;
  }

  graph.edge_count++;
  if(graph.edge_count > REASONABLE){
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
  graph.EDGES[graph.edge_count] = edge;

  // use a linked list to store the bond, can also check if it already exists

  WLNEdge *curr = parent->bonds;
  if(curr){
    
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

  parent->num_children++;
  return edge;
}


void debug_edge(WLNEdge *edge){
  if(!edge)
    fprintf(stderr,"Error: debugging nullptr edge\n");  
  else
    fprintf(stderr,"%c -- %d --> %c\n",edge->parent->ch, edge->order ,edge->child->ch);
}


WLNEdge *search_edge(WLNSymbol *child, WLNSymbol*parent, bool verbose=true){
  if(!child || !parent){
    fprintf(stderr,"Error: searching edge on nullptrs\n");
    return 0;
  }
  
  WLNEdge *edge = 0;
  for (edge=parent->bonds;edge;edge = edge->nxt){
    if(edge->child == child)
      return edge;
  }
  if(verbose)
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
  
  head->num_edges--;
  edge->child->num_edges--;

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



/**********************************************************************
                          WLNRing Functions
**********************************************************************/

WLNRing *AllocateWLNRing(WLNGraph &graph)
{
  graph.ring_count++;
  if(graph.ring_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln rings - is this reasonable?\n");
    return 0;
  }

  WLNRing *wln_ring = new WLNRing;
  graph.RINGS[graph.ring_count] = wln_ring;
  return wln_ring;
}


// both lookups needed for QOL in ring building
WLNSymbol* assign_locant(unsigned char loc,WLNSymbol *locant, WLNRing *ring){
    
  if(!locant)
    return 0;
  
  ring->locants[loc] = locant; 
  ring->locants_ch[locant] = loc;
  locant->type = RING;
  return locant; 
}  




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



/**********************************************************************
                         High Level Parser Functions
**********************************************************************/


/* dump wln tree to a dotvis file */
void WLNDumpToDot(FILE *fp, WLNGraph &graph)
{  
  fprintf(fp, "digraph WLNdigraph {\n");
  fprintf(fp, "  rankdir = LR;\n");
  for (unsigned int i=0; i<=graph.symbol_count;i++)
  {
    WLNSymbol *node = graph.SYMBOLS[i];
    if(!node)
      continue;

    fprintf(fp, "  %d", graph.index_lookup[node]);
    if (node->ch == '*')
      fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
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
      if (bond_order > 1){
        for (unsigned int k=0;k<bond_order;k++){
          fprintf(fp, "  %d", graph.index_lookup[node]);
          fprintf(fp, " -> ");
          fprintf(fp, "%d\n", graph.index_lookup[child]);
        }
      }
      else{
        fprintf(fp, "  %d", graph.index_lookup[node]);
        fprintf(fp, " -> ");
        fprintf(fp, "%d\n", graph.index_lookup[child]);
      }
    }
  }

  fprintf(fp, "}\n");
}

bool WriteGraph(WLNGraph &graph){
  fprintf(stderr,"Dumping wln graph to wln-graph.dot:\n");
  FILE *fp = 0;
  fp = fopen("wln-graph.dot", "w");
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




/**********************************************************************
                         BABEL Mol Functions
**********************************************************************/


// holds all the functions for WLN graph conversion, mol object is assumed ALIVE AT ALL TIMES
// uses old NM functions from previous methods: Copyright (C) NextMove Software 2019-present
struct BabelGraph{

  std::map<OpenBabel::OBAtom*, WLNSymbol*>  atom_symbol_map; 
  std::map<WLNSymbol*,OpenBabel::OBAtom*>   symbol_atom_map; 
  
  bool CYCLIC = 0;


  BabelGraph(){};
  ~BabelGraph(){};


  WLNSymbol* CreateWLNNode(OpenBabel::OBAtom* atom,  WLNGraph &graph){

    if(!atom){
      fprintf(stderr,"Error: nullptr OpenBabel Atom*\n");
      return 0; 
    }

    WLNSymbol *node = 0;
    switch(atom->GetAtomicNum()){
      case 1:
        node = AllocateWLNSymbol('H',graph);
        break; 

      case 5:
        node = AllocateWLNSymbol('B',graph);
        break;

      case 6:
        node = AllocateWLNSymbol('C',graph);
        break;

      case 7:
        node = AllocateWLNSymbol('N',graph);
        break;
      
      case 8:
        if(atom->GetExplicitValence() == 1 && atom->GetFormalCharge() != -1)
          node = AllocateWLNSymbol('Q',graph);
        else
          node = AllocateWLNSymbol('O',graph);
        break;
      
      case 9:
        node = AllocateWLNSymbol('F',graph);
        break;

      case 15:
        node = AllocateWLNSymbol('P',graph);
        break;

      case 16:
        node = AllocateWLNSymbol('S',graph);
        break;

      case 17:
        node = AllocateWLNSymbol('G',graph);
        break;

      case 35:
        node = AllocateWLNSymbol('E',graph);
        break;

      case 53:
        node = AllocateWLNSymbol('I',graph);
        break;

      default:
        fprintf(stderr,"Error: unhandled element for WLNSymbol formation\n");
        return 0;
    }
    
    if(!graph.root)
      graph.root = node; 

    return node; 
  }


  /* add the starting atom to build a tree for a locant position */
  WLNSymbol* BuildWLNTree(OpenBabel::OBAtom* start_atom, OpenBabel::OBMol *mol,WLNGraph &graph){

    // has to be done as DFS in order to keep bond direction on the tree
    WLNSymbol *root = 0; 
    WLNSymbol *node   = 0; 
    WLNSymbol *child  = 0; 
    WLNEdge   *edge   = 0; 

    OpenBabel::OBAtom* atom = start_atom;
    std::map<OpenBabel::OBAtom*,bool> visited; 
    std::stack<OpenBabel::OBAtom*> atom_stack; 
    atom_stack.push(atom); 

    while(!atom_stack.empty()){
      atom = atom_stack.top(); 
      atom_stack.pop();
      visited[atom] = true;

      // negative oxygen should be given a W character, therefore pointing in
      if(atom->GetFormalCharge() == -1 && atom->GetAtomicNum() == 8){
        FOR_NBORS_OF_ATOM(iterator, atom){
          OpenBabel::OBAtom *neighbour = &(*iterator);
          if(!visited[neighbour])
            atom_stack.push(neighbour); 
        }
        continue;
      }

      // create the first atom if needed
      if(!atom_symbol_map[atom]){
        node = CreateWLNNode(atom,graph); 

        if(!root)
          root = node;

        if(!node){
          fprintf(stderr,"Error: could not create node in BuildWLNTree\n");
          return 0;
        }
        node->set_edge_and_type(atom->GetTotalValence(),STANDARD); // allow smiles to dictate
        atom_symbol_map[atom] = node; 
        symbol_atom_map[node] = atom;
      }
      else
        node = atom_symbol_map[atom]; 

      // this will look back, so order is important to maintain
      FOR_NBORS_OF_ATOM(iterator, atom){
        OpenBabel::OBAtom *neighbour = &(*iterator);
        if(!atom_symbol_map[neighbour]){
          child = CreateWLNNode(neighbour,graph); 
          if(!child){
            fprintf(stderr,"Error: could not create node in BuildWLNTree\n");
            return 0;
          }
          child->set_edge_and_type(neighbour->GetTotalValence(),STANDARD); // allow smiles to dictate
          atom_symbol_map[neighbour] = child; 
          symbol_atom_map[child] = neighbour; 

          // bond here, and don't consider the symbol if the atom is already made 
          OpenBabel::OBBond *bond = atom->GetBond(neighbour); 
          if(!bond){
            fprintf(stderr,"Error: accessing non-existent bond in BuildWLNTree\n");
            return 0;
          }
          unsigned int order = bond->GetBondOrder(); 
          edge = AllocateWLNEdge(child,node,graph); 
          if(order > 1)
            edge = unsaturate_edge(edge,order-1);
        }
          
        if(!visited[neighbour])
          atom_stack.push(neighbour); 
      } 

    }
    
    return root;
  }


  // reads the babel graph into the appropriate wln graph
  std::vector<WLNSymbol*> ReadBabelGraph(OpenBabel::OBMol *mol,WLNGraph &graph){
    // no cycles, build from root node and output
    std::vector<WLNSymbol*> roots; 
    if(mol->GetSSSR().empty()){
      OpenBabel::OBAtomAtomIter iter; 
      FOR_ATOMS_OF_MOL(iter,mol){
        OpenBabel::OBAtom *atom = &(*iter);
        if(!atom_symbol_map[atom]){
          WLNSymbol *root_node = BuildWLNTree (atom,mol,graph); 
          if(!root_node){
            fprintf(stderr,"Error: failure in building tree from source atom\n");
            return {}; 
          }
          else
            roots.push_back(root_node); 
        }
      }
    }
    else
      CYCLIC = true; // set the flag, temp for now
    
    return roots; 
  }

  bool ParseWLNGraph(std::vector<WLNSymbol*> &roots, WLNGraph &graph, std::string &buffer){
    if(!CYCLIC){
      for (unsigned int i=0;i<roots.size();i++){
        if(!WriteWLNFromNode(roots[i],graph,buffer))
          return false;

        if(i<roots.size()-1)
          buffer += " &"; // add for ionic
      }
    }
    return true; 
  }



  // will also add to handled
  bool CheckOXO(WLNSymbol *sym, std::map<WLNSymbol*,bool> &handled){

    WLNEdge *edge = 0; 
    WLNEdge *e_oxo_1 = 0; 
    WLNEdge *e_oxo_2 = 0; 

    for(edge=sym->bonds;edge;edge=edge->nxt){
      if((edge->child->ch == 'O') && (edge->order == 2 || symbol_atom_map[edge->child]->GetFormalCharge() == -1)){
        if(!e_oxo_1)
          e_oxo_1 = edge; 
        else
          e_oxo_2 = edge; 
      }
    }

    if(!e_oxo_1 || !e_oxo_2)
      return false;
    else{
      handled[e_oxo_1->child] = true;
      handled[e_oxo_2->child] = true;
      return true;
    } 
  }

  unsigned int CheckCarbonChain(WLNSymbol *sym, std::map<WLNSymbol*,bool> &handled){
    unsigned int carbons = 0; 
    WLNSymbol *carbon_sym = sym;

    while(carbon_sym->num_children == 1){
      carbons++;
      handled[carbon_sym] = true;

      if(carbon_sym->bonds && carbon_sym->bonds->order == 1 && carbon_sym->bonds->child->ch == 'C')
        carbon_sym = carbon_sym->bonds->child;
      else 
        break;
    }

    if(!carbon_sym->num_children && !handled[carbon_sym]){
      carbons++; 
      handled[carbon_sym] = true;
    }

    return carbons; 
  }



  bool WriteWLNFromNode(WLNSymbol *root,WLNGraph &graph,std::string &buffer){

    // dfs style notational build from a given root, use to build the notation 

    WLNSymbol *top = 0; 
    WLNSymbol *prev = 0; 
    WLNEdge   *edge = 0; // for iterating
    
    std::stack<std::pair<WLNSymbol*,unsigned int>> stack; 
    std::stack <WLNSymbol*> branch_stack; // this is for non-cyclics, simpler than writing
    std::map<WLNSymbol*,bool> visited;
    std::map<WLNSymbol*,bool> handled; // if a character corresponds to multiple symbols 

    
    unsigned int order = 0; 
    unsigned int carbon_chain = 0;   


    stack.push({root,0});

    while(!stack.empty()){
      top = stack.top().first;
      order = stack.top().second; 
      

// branching returns
      if( top->previous && top->previous != prev && 
          !branch_stack.empty() && !handled[top]){

        while(!branch_stack.empty() && prev != branch_stack.top()){
          buffer += '&';
          branch_stack.pop();
        }
      }
      
      stack.pop();
      visited[top] = true; 
      prev = top; // use for checks on branching stack

      if(!handled[top]){

// bond unsaturations
        if(order == 2)
          buffer+='U';
        if(order == 3)
          buffer+="UU";        

        switch (top->ch){

// oxygens
          case 'O':
            buffer += 'O';
            break;

          case 'Q':
            buffer += 'Q';
            if(!branch_stack.empty())
              prev = branch_stack.top(); 
            break;

// carbons 
          case 'C':
            // special V case
            if(top->num_children == 2){
              
              if(CheckOXO(top, handled)){
                buffer += 'W'; 
                break;  // immediate break out
              } 
              else{
                // V case
                for(edge=top->bonds;edge;edge=edge->nxt){
                  if(edge->order == 2 && edge->child->ch == 'O'){
                    buffer+='V';
                    handled[edge->child] = true; 
                    goto out_switch;
                  }
                }
              }


            }
            
            if(top->num_children > 1){
              if(top->num_edges == 4)
                buffer += 'X';
              else if(top->num_edges == 3)
                buffer += 'Y'; 

              if(CheckOXO(top, handled))
                buffer += 'W'; 
              
              branch_stack.push(top);
            }
            else{
              carbon_chain = CheckCarbonChain(top,handled);
              char digits[4] = {0}; 
              sprintf (digits,"%d",carbon_chain);
              buffer += digits;
            }
            break;

// nitrogen
          case 'N':
            if(top->num_edges == 1)
              buffer += 'Z';
            else if(top->num_edges == 2)
              buffer += 'M';
            else if (top->num_edges == 3){
              buffer += 'N';
              if(CheckOXO(top, handled))
                buffer += 'W'; 

              branch_stack.push(top);
            }
            else{
              if(CheckOXO(top, handled)){
                buffer += 'N';
                buffer += 'W';
              }
              else{
                buffer += 'K';
                branch_stack.push(top); 
              }
            }
            break;


// halogens
          case 'E':
          case 'F':
          case 'G':
          case 'I':
            buffer += top->ch;
            if(!top->num_edges && symbol_atom_map[top]->GetFormalCharge() == 0)
              buffer += 'H';

            if(!branch_stack.empty())
              prev = branch_stack.top(); 
            break;

// branching heteroatoms 
          case 'S':
            branch_stack.push(top);
            buffer += 'S';
            if(CheckOXO(top, handled))
              buffer += 'W'; 
            break;

          default:
            fprintf(stderr,"Error: unhandled WLN char %c\n",top->ch); 
            return false; 
        }
      }

      out_switch:
      for(edge=top->bonds;edge;edge=edge->nxt){
        if(!visited[edge->child])
          stack.push({edge->child,edge->order}); 
      }

    }
    return true;
  }

};



/**********************************************************************
                         API FUNCTION
**********************************************************************/


bool WriteWLN(std::string &buffer, OpenBabel::OBMol* mol)
{   
 
  WLNGraph wln_graph;
  BabelGraph obabel; 

  bool state = true;
  std::vector<WLNSymbol*> roots; 

  roots = obabel.ReadBabelGraph(mol,wln_graph);
  if(roots.empty() && !obabel.CYCLIC)
    return false; 
  
  // create an optional wln dotfile
  if (opt_wln2dot)
    WriteGraph(wln_graph);

  if(state)
    state = obabel.ParseWLNGraph(roots,wln_graph,buffer);

  return state;
}


static void DisplayUsage()
{
  fprintf(stderr, "writewln <options> -i<format> -s <input (escaped)>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -d                    print debug messages to stderr\n");
  fprintf(stderr, "  -h                    show the help for executable usage\n");
  fprintf(stderr, "  -i                    choose input format (-ismi, -iinchi, -ican)\n");
  fprintf(stderr, "  -w                    dump wln trees to dot file in [build]\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser writes to wiswesser\n"
                  " line notation (wln) from smiles/inchi, the parser is native\n"
                  " and will can return either a reformatted string*\n"
                  " *if rules do not parse exactly, and the connection\n"
                  " table which can be used in other libraries\n");
  DisplayUsage();
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i, j;

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
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;

  conv.SetInFormat(format);
  res = conv.ReadString(&mol,cli_inp);

  std::string buffer;
  buffer.reserve(1000);
  if(!WriteWLN(buffer,&mol))
    return 1;
  
  std::cout << buffer << std::endl;

  return 0;
}


