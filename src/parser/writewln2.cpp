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

  unsigned int num_children;  // specifically for forward edges
  unsigned int on_child;      // which branch are we on for stack tracking

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
    on_child = 0; 
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

  std::map<OBAtom*, WLNSymbol*>  atom_symbol_map; 
  std::map<WLNSymbol*,OBAtom*>   symbol_atom_map; 
  
  BabelGraph(){};
  ~BabelGraph(){};


  WLNSymbol* CreateWLNNode(OBAtom* atom, WLNGraph &graph){

    if(!atom){
      fprintf(stderr,"Error: nullptr OpenBabel Atom*\n");
      return 0; 
    }

    unsigned int neighbours = 0; 
    unsigned int orders = 0; 
    OBAtom *neighbour = 0; 
    OBBond *bond = 0; 

    WLNSymbol *node = 0;
    switch(atom->GetAtomicNum()){
      case 1:
        node = AllocateWLNSymbol('H',graph);
        node->set_edge_and_type(1);
        break; 

      case 5:
        node = AllocateWLNSymbol('B',graph);
        node->set_edge_and_type(3);
        break;

      case 6:
        FOR_NBORS_OF_ATOM(iterator, atom){
          neighbour = &(*iterator);
          bond = atom->GetBond(neighbour);
          orders += bond->GetBondOrder(); 
          neighbours++;
        }
        if(neighbours <= 2){
          node = AllocateWLNSymbol('1',graph);
          node->set_edge_and_type(4);
        }
        else if(neighbours > 2){
          if(orders == 3){
            node = AllocateWLNSymbol('Y',graph);
            node->set_edge_and_type(3);
          }
          else{
            node = AllocateWLNSymbol('X',graph);
            node->set_edge_and_type(4);
          }
        }
        else{
          node = AllocateWLNSymbol('C',graph);
          node->set_edge_and_type(4);
        }
        break;
      
      case 7:
        node = AllocateWLNSymbol('N',graph);
        node->set_edge_and_type(atom->GetExplicitValence());
        break;
      
      case 8:
        if(atom->GetExplicitValence() < 2 && atom->GetFormalCharge() != -1){
          node = AllocateWLNSymbol('Q',graph);
          node->set_edge_and_type(1);
        }
        else{
          node = AllocateWLNSymbol('O',graph);
          node->set_edge_and_type(2);
        }
        break;
      
      case 9:
        node = AllocateWLNSymbol('F',graph);
        node->set_edge_and_type(atom->GetExplicitValence());
        break;

      case 15:
        node = AllocateWLNSymbol('P',graph);
        node->set_edge_and_type(6);
        break;

      case 16:
        node = AllocateWLNSymbol('S',graph);
        node->set_edge_and_type(6);
        break;

      case 17:
        node = AllocateWLNSymbol('G',graph);
        node->set_edge_and_type(atom->GetExplicitValence());
        break;

      case 35:
        node = AllocateWLNSymbol('E',graph);
        node->set_edge_and_type(atom->GetExplicitValence());
        break;

      case 53:
        node = AllocateWLNSymbol('I',graph);
        node->set_edge_and_type(atom->GetExplicitValence());
        break;



// all special elemental cases

      case 89:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AC";
        break;

      case 47:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AG";
        break;
    
      case 13:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AL";
        break;

      case 95:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AM";
        break;

      case 18:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AR";
        break;

      case 33:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AS";
        break;

      case 85:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AT";
        break;

      case 79:
        node = AllocateWLNSymbol('*',graph);
        node->special += "AU";
        break;


      case 56:
        node = AllocateWLNSymbol('*',graph);
        node->special += "BA";
        break;

      case 4:
        node = AllocateWLNSymbol('*',graph);
        node->special += "BE";
        break;

      case 107:
        node = AllocateWLNSymbol('*',graph);
        node->special += "BH";
        break;

      case 83:
        node = AllocateWLNSymbol('*',graph);
        node->special += "BI";
        break;

      case 97:
        node = AllocateWLNSymbol('*',graph);
        node->special += "BK";
        break;

      case 20:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CA";
        break;
      
      case 48:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CD";
        break;

      case 58:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CE";
        break;

      case 98:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CF";
        break;

      case 96:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CN";
        break;

      case 112:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CN";
        break;

      case 27:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CO";
        break;

      case 24:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CR";
        break;

      case 55:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CS";
        break;

      case 29:
        node = AllocateWLNSymbol('*',graph);
        node->special += "CU";
        break;

      case 105:
        node = AllocateWLNSymbol('*',graph);
        node->special += "DB";
        break;

      case 110:
        node = AllocateWLNSymbol('*',graph);
        node->special += "DS";
        break;

      case 66:
        node = AllocateWLNSymbol('*',graph);
        node->special += "DY";
        break;

      case 68:
        node = AllocateWLNSymbol('*',graph);
        node->special += "ER";
        break;

      case 99:
        node = AllocateWLNSymbol('*',graph);
        node->special += "ES";
        break;

      case 63:
        node = AllocateWLNSymbol('*',graph);
        node->special += "EU";
        break;

      case 26:
        node = AllocateWLNSymbol('*',graph);
        node->special += "FE";
        break;

      case 114:
        node = AllocateWLNSymbol('*',graph);
        node->special += "FL";
        break;

      case 100:
        node = AllocateWLNSymbol('*',graph);
        node->special += "FM";
        break;

      case 87:
        node = AllocateWLNSymbol('*',graph);
        node->special += "FR";
        break;

      case 31:
        node = AllocateWLNSymbol('*',graph);
        node->special += "GA";
        break;

      case 64:
        node = AllocateWLNSymbol('*',graph);
        node->special += "GD";
        break;

      case 32:
        node = AllocateWLNSymbol('*',graph);
        node->special += "GE";
        break;

      case 2:
        node = AllocateWLNSymbol('*',graph);
        node->special += "HE";
        break;

      case 72:
        node = AllocateWLNSymbol('*',graph);
        node->special += "HF";
        break;

      case 80:
        node = AllocateWLNSymbol('*',graph);
        node->special += "HG";
        break;

      case 67:
        node = AllocateWLNSymbol('*',graph);
        node->special += "HO";
        break;

      case 108:
        node = AllocateWLNSymbol('*',graph);
        node->special += "HS";
        break;

      case 49:
        node = AllocateWLNSymbol('*',graph);
        node->special += "IN";
        break;

      case 77:
        node = AllocateWLNSymbol('*',graph);
        node->special += "IR";
        break;

      case 36:
        node = AllocateWLNSymbol('*',graph);
        node->special += "KR";
        break;

      case 19:
        node = AllocateWLNSymbol('*',graph);
        node->special += "KA";
        break;

      case 57:
        node = AllocateWLNSymbol('*',graph);
        node->special += "LA";
        break;

      case 3:
        node = AllocateWLNSymbol('*',graph);
        node->special += "LI";
        break;

      case 103:
        node = AllocateWLNSymbol('*',graph);
        node->special += "LR";
        break;

      case 71:
        node = AllocateWLNSymbol('*',graph);
        node->special += "LU";
        break;

      case 116:
        node = AllocateWLNSymbol('*',graph);
        node->special += "LV";
        break;

      case 115:
        node = AllocateWLNSymbol('*',graph);
        node->special += "MC";
        break;

      case 101:
        node = AllocateWLNSymbol('*',graph);
        node->special += "MD";
        break;

      case 12:
        node = AllocateWLNSymbol('*',graph);
        node->special += "MG";
        break;

      case 25:
        node = AllocateWLNSymbol('*',graph);
        node->special += "MN";
        break;

      case 42:
        node = AllocateWLNSymbol('*',graph);
        node->special += "MO";
        break;

      case 109:
        node = AllocateWLNSymbol('*',graph);
        node->special += "MT";
        break;

      case 11:
        node = AllocateWLNSymbol('*',graph);
        node->special += "NA";
        break;

      case 41:
        node = AllocateWLNSymbol('*',graph);
        node->special += "NB";
        break;

      case 60:
        node = AllocateWLNSymbol('*',graph);
        node->special += "ND";
        break;

      case 10:
        node = AllocateWLNSymbol('*',graph);
        node->special += "NE";
        break;

      case 113:
        node = AllocateWLNSymbol('*',graph);
        node->special += "NH";
        break;

      case 28:
        node = AllocateWLNSymbol('*',graph);
        node->special += "NI";
        break;

      case 102:
        node = AllocateWLNSymbol('*',graph);
        node->special += "NO";
        break;

      case 93:
        node = AllocateWLNSymbol('*',graph);
        node->special += "NP";
        break;


      case 118:
        node = AllocateWLNSymbol('*',graph);
        node->special += "OG";
        break;

      case 76:
        node = AllocateWLNSymbol('*',graph);
        node->special += "OS";
        break;


      case 91:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PA";
        break;

      case 82:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PB";
        break;

      case 46:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PD";
        break;

      case 61:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PM";
        break;

      case 84:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PO";
        break;

      case 59:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PR";
        break;

      case 78:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PT";
        break;

      case 94:
        node = AllocateWLNSymbol('*',graph);
        node->special += "PU";
        break;

      case 88:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RA";
        break;

      case 37:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RB";
        break;

      case 75:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RE";
        break;

      case 104:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RF";
        break;

      case 111:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RG";
        break;

      case 45:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RH";
        break;

      case 86:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RN";
        break;

      case 44:
        node = AllocateWLNSymbol('*',graph);
        node->special += "RU";
        break;

      case 51:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SB";
        break;

      case 21:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SC";
        break;

      case 34:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SE";
        break;

      case 106:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SG";
        break;

      case 14:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SI";
        break;

      case 62:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SM";
        break;

      case 50:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SN";
        break;

      case 38:
        node = AllocateWLNSymbol('*',graph);
        node->special += "SR";
        break;


      case 73:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TA";
        break;

      case 65:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TB";
        break;

      case 43:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TC";
        break;

      case 52:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TE";
        break;

      case 90:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TH";
        break;

      case 22:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TI";
        break;

      case 81:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TL";
        break;

      case 69:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TM";
        break;

      case 117:
        node = AllocateWLNSymbol('*',graph);
        node->special += "TS";
        break;

      case 92:
        node = AllocateWLNSymbol('*',graph);
        node->special += "UR";
        break;

      case 23:
        node = AllocateWLNSymbol('*',graph);
        node->special += "VA";
        break;

      case 54:
        node = AllocateWLNSymbol('*',graph);
        node->special += "XE";
        break;

      case 39:
        node = AllocateWLNSymbol('*',graph);
        node->special += "YT";
        break;

      case 70:
        node = AllocateWLNSymbol('*',graph);
        node->special += "YB";
        break;

      case 30:
        node = AllocateWLNSymbol('*',graph);
        node->special += "ZN";
        break;

      case 40:
        node = AllocateWLNSymbol('*',graph);
        node->special += "ZR";
        break;
      

      default:
        fprintf(stderr,"Error: unhandled element for WLNSymbol formation\n");
        return 0;
    }
    
    if(!graph.root)
      graph.root = node; 

    if(!node->allowed_edges)
      node->set_edge_and_type(8);

    return node; 
  }


  /* add the starting atom to build a tree for a locant position */
  WLNSymbol* BuildWLNTree(OBAtom* start_atom, OBMol *mol,WLNGraph &graph){

    // has to be done as DFS in order to keep bond direction on the tree
    WLNSymbol *root = 0; 
    WLNSymbol *node   = 0; 
    WLNSymbol *child  = 0; 
    WLNEdge   *edge   = 0; 

    OBAtom* atom = start_atom;
    std::map<OBAtom*,bool> visited; 
    std::stack<OBAtom*> atom_stack; 
    atom_stack.push(atom); 

    while(!atom_stack.empty()){
      atom = atom_stack.top(); 
      atom_stack.pop();
      visited[atom] = true;

      // negative oxygen should be given a W character, therefore pointing in
      if(atom->GetFormalCharge() == -1 && atom->GetAtomicNum() == 8){
        FOR_NBORS_OF_ATOM(iterator, atom){
          OBAtom *neighbour = &(*iterator);
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
        
        atom_symbol_map[atom] = node; 
        symbol_atom_map[node] = atom;
      }
      else
        node = atom_symbol_map[atom]; 

      // this will look back, so order is important to maintain
      FOR_NBORS_OF_ATOM(iterator, atom){
        OBAtom *neighbour = &(*iterator);
        if(!atom_symbol_map[neighbour]){
          child = CreateWLNNode(neighbour,graph); 
          if(!child){
            fprintf(stderr,"Error: could not create node in BuildWLNTree\n");
            return 0;
          }

          atom_symbol_map[neighbour] = child; 
          symbol_atom_map[child] = neighbour; 

          // bond here, and don't consider the symbol if the atom is already made 
          OBBond *bond = atom->GetBond(neighbour); 
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


  /* reads the babel graph into the appropriate wln graph
      either return roots for ionic species to build tree,
      or if cyclic, return ring objects to parse  */


  std::vector<WLNSymbol*> ReadBabelNonCyclic(OBMol *mol,WLNGraph &graph){
    // no cycles, build from root node and output
    
    std::vector<WLNSymbol*> roots; 
    OBAtomAtomIter iter; 
    FOR_ATOMS_OF_MOL(iter,mol){
      OBAtom *atom = &(*iter);
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
  
    return roots; 
  }


  /* handles non cyclic and ionic species from roots vector */
  bool ParseNonCyclicWLN(std::vector<WLNSymbol*> &roots, WLNGraph &graph, std::string &buffer){    
    for (unsigned int i=0;i<roots.size();i++){
      if(!WriteWLNFromNode(roots[i],graph,buffer))
        return false;

      if(i<roots.size()-1)
        buffer += " &"; // add for ionic
    }
    return true; 
  }

  // will also add to handled
  bool CheckCarbonyl(WLNSymbol *sym, std::map<WLNSymbol*,bool> &visited){
    WLNEdge *edge = 0; 
    WLNEdge *oxygen = 0; 
  
    for(edge=sym->bonds;edge;edge=edge->nxt){
      if((edge->child->ch == 'O') && (edge->order == 2 || symbol_atom_map[edge->child]->GetFormalCharge() == -1)){
        oxygen = edge; 
        break;
      }
    }
    if(!oxygen)
      return false;
    else{
      visited[oxygen->child] = true;
      return true;
    } 
  }


  // will also add to handled
  bool CheckDIOXO(WLNSymbol *sym, std::map<WLNSymbol*,bool> &visited){

    WLNEdge *edge = 0; 
    // needs to a be a priority, so taking a double bond =O over a =O + -O-
    std::deque<WLNSymbol*> oxygens;
    for(edge=sym->bonds;edge;edge=edge->nxt){
      // highest priority double bond =O
      if(edge->child->ch == 'O' && edge->order == 2)
        oxygens.push_front(edge->child); 
      // lower priority  =O + -O-
      else if(edge->child->ch == 'O' && symbol_atom_map[edge->child]->GetFormalCharge() == -1)
        oxygens.push_back(edge->child); 
    }

    if(oxygens.size() < 2)
      return false;
    else{
      visited[oxygens[0]] = true;
      visited[oxygens[1]] = true;
      return true;
    } 
  }

  // writes to the buffer
  WLNSymbol* WriteCarbonChain(WLNSymbol *sym, std::string &buffer){
    
    unsigned int carbons = 1; 
    WLNSymbol *carbon_sym = sym;

    while(carbon_sym->bonds && carbon_sym->bonds->child->ch == '1' && carbon_sym->bonds->order == 1){
      carbons++;
      carbon_sym = carbon_sym->bonds->child;
    }

    buffer += std::to_string(carbons);
    return carbon_sym;
  }



  bool WriteWLNFromNode(WLNSymbol *root,WLNGraph &graph,std::string &buffer){

    // dfs style notational build from a given root, use to build the notation 
    WLNSymbol *top = 0; 
    WLNSymbol *prev = 0; 
    WLNEdge   *edge = 0; // for iterating
    
    std::stack<std::pair<WLNSymbol*,unsigned int>> stack; 
    std::stack <WLNSymbol*> branch_stack; 
    std::map<WLNSymbol*,bool> visited;
    
    bool following_terminator = false;
    unsigned int order = 0;  

    stack.push({root,0});

    while(!stack.empty()){
      top = stack.top().first;
      order = stack.top().second; 
      
// branching returns
      if( (top->previous && prev) && top->previous != prev && 
          !branch_stack.empty()){

        prev = top->previous;
        
        if(opt_debug)
          fprintf(stderr,"%c is on branch: %d\n",prev->ch,prev->on_child);

        // distinction between a closure and a pop
        if(!following_terminator)
          buffer += '&';

        WLNSymbol *branch_top = 0; 
        while(!branch_stack.empty() && prev != branch_stack.top()){
          branch_top = branch_stack.top();
          if(opt_debug)
            fprintf(stderr,"stack_top: %c - %d\n",branch_top->ch,branch_top->on_child);
          
          if( (branch_top->num_children != branch_top->on_child) 
              || branch_top->num_edges < branch_top->allowed_edges)
            buffer += '&';

          branch_stack.pop();
        }

        prev->on_child++;
      }
      else if(prev)
        prev->on_child++;

      following_terminator = false;
      
      stack.pop();
      visited[top] = true;
      prev = top;      

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
          if(!top->num_edges)
            buffer += 'H';

          if(!branch_stack.empty()){
            prev = branch_stack.top();
            following_terminator = true;
          }
          break;

// carbons 
        // alkyl chain 
        case '1':
          top = WriteCarbonChain(top,buffer);
          prev = top;
          break;

        case 'Y':
        case 'X':
          if(CheckDIOXO(top, visited)){
            buffer += top->ch;
            buffer += 'W'; 
          }
          else if(CheckCarbonyl(top,visited))
            buffer += 'V'; 
          else{
            buffer += top->ch;
            branch_stack.push(top);
          }
          break;


// nitrogen
        case 'N':
          if(top->num_edges < 2){
            buffer += 'Z';
            if(!top->num_edges)
              buffer += 'H';

            if(!branch_stack.empty()){
              prev = branch_stack.top();
              following_terminator = true;
            }
          }
          else if(top->num_children < 2 && top->num_edges < 3)
            buffer += 'M';
          else if (top->num_children < 3 && top->num_edges < 4){
            buffer += 'N';
            if(CheckDIOXO(top, visited))
              buffer += 'W'; 

            branch_stack.push(top);
          }
          else{
            if(CheckDIOXO(top, visited)){
              buffer += 'N';
              buffer += 'W';
            }
            else{
              buffer += 'K';
              branch_stack.push(top); // implied methyl, must add to branch
            }
          }
          break;


// halogens
        case 'E':
        case 'F':
        case 'G':
        case 'I':
          if(top->num_edges > 1){
            buffer += '-';
            buffer += top->ch;
            buffer += '-';
            if(CheckDIOXO(top, visited))
              buffer += 'W'; 

            branch_stack.push(top);
          }
          else{
            buffer += top->ch;
            if(!top->num_edges && symbol_atom_map[top]->GetFormalCharge() == 0)
              buffer += 'H';

            if(!branch_stack.empty()){
              prev = branch_stack.top();
              following_terminator = true;
            }
          }
          break;

// branching heteroatoms 
        case 'B':
        case 'S':
        case 'P':
          buffer += top->ch;
          if(CheckDIOXO(top, visited))
            buffer += 'W'; 

          if(top->num_children > 0)
            branch_stack.push(top);
          
          break;


// specials 
        case '*':
          buffer += '-';
          buffer += top->special;
          buffer += '-';
          if(!top->num_edges && symbol_atom_map[top]->GetFormalCharge() == 0)
            buffer += 'H';
          else if(top->num_children > 0)
            branch_stack.push(top);
          
          break;

        default:
          fprintf(stderr,"Error: unhandled WLN char %c\n",top->ch); 
          return false; 
      }
    

      for(edge=top->bonds;edge;edge=edge->nxt){
        if(!visited[edge->child])
          stack.push({edge->child,edge->order}); 
      }

    }
    return true;
  }

  /* works on priority, and creates locant path via array shifting */
  OBAtom **CreateLocantArray(OBMol *mol, std::set<OBRing*> &local_SSSR, std::map<OBAtom*,unsigned int> &ring_shares, unsigned int size){

    OBAtom **locant_path = (OBAtom**)malloc(sizeof(OBAtom*) * size); 
    for(unsigned int i=0;i<size;i++)
      locant_path[i] = 0; 

   
    OBAtom *ratom = 0;
    unsigned int rings_handled = 0; 
    unsigned int locant_pos = 0;

    std::map<OBRing*,bool> rings_seen; 
    std::map<OBAtom*,bool> atoms_seen; 
    OBRing *obring = *(local_SSSR.begin());

     // add into the array directly for precondition, shift until highest priority is found  
    unsigned int init_priority = 0;
    for(unsigned int i=0;i<obring->_path.size();i++){
      ratom = mol->GetAtom(obring->_path[i]);
      locant_path[locant_pos++] = ratom;
      atoms_seen[ratom] = true;
      if(init_priority < ring_shares[ratom])
        init_priority = ring_shares[ratom];  
    }

    // shift so A is guaranteed a ring share
    while(ring_shares[locant_path[0]] != init_priority){
      locant_path[locant_pos] = locant_path[0];
      for(unsigned int i=0;i<size-1;i++)
        locant_path[i] = locant_path[i+1];
    }

    if(opt_debug){
      fprintf(stderr,"  locant path: ");
      print_locant_array(locant_path,size); 
    }

    while(rings_handled < local_SSSR.size()-1){

      rings_seen[obring] = true;
      rings_handled++; 

      // find the first point seen with an external ring condition
      OBAtom *next_seed = 0;
      unsigned int hp_pos = 0; 

      // this needs to be done on the array, with a external ring check
      for(unsigned int i=0;i<locant_pos;i++){
        bool found = false;
        ratom = locant_path[i];
        if(ring_shares[ratom] > 1){
          // find out if its pointing at a ring we havent yet considered
          for(std::set<OBRing*>::iterator iter = local_SSSR.begin(); iter != local_SSSR.end(); iter++){
            if(!rings_seen[(*iter)] && (*iter)->IsInRing(ratom->GetIdx())){
              next_seed = ratom;
              hp_pos = i;
              obring = (*iter);
              found = true;
              break;
            }
          }
          if(found)
            break;
        }
      }

      if(opt_debug)
        fprintf(stderr,"  shift %d from position %d\n",obring->_path.size(),hp_pos);
      
      // --- shift and add procedure ---
      // rings atoms must flow clockwise in _path to form the locant path correctly
      // is this an obabel default? --> algorithm relies on it
      unsigned int j=0;
      for(unsigned int i=0;i<obring->_path.size();i++){
        ratom = mol->GetAtom(obring->_path[i]);
        if(!atoms_seen[ratom]){
          // shift
          for(int k=size-1;k>hp_pos+j;k--)
            locant_path[k]= locant_path[k-1];
            
          locant_path[hp_pos+1+j] = ratom;
          atoms_seen[ratom] = true;
          j++;
          locant_pos++;
        }
      }

      if(opt_debug){
        fprintf(stderr,"  locant path: ");
        print_locant_array(locant_path,size); 
      }


    }

    return locant_path;
  }


  /* Works from a starting ring atom, build up local SSSR from there*/
  WLNRing* ReadBabelCyclic(OBAtom *ring_root, std::string &buffer,OBMol *mol, WLNGraph &graph){

    if(opt_debug)
    fprintf(stderr,"Building ring\n");

    if(!ring_root){
      fprintf(stderr,"Error: ring root is nullptr\n");
      return 0; 
    }

    WLNRing *wln_ring = 0; 
    wln_ring = AllocateWLNRing(graph);

    std::map<OBAtom*,bool> visited; 
    std::stack<OBAtom*> atom_stack; 
    
    bool hetero = false; 
    unsigned int size = 0; 
    unsigned int in_rings = 0; 
    OBAtom *atom = 0; 
    OBAtom *neighbour = 0; 


    unsigned int fuses = 0;
    unsigned int multicyclic = 0;
    unsigned int branching = 0;   
    std::set<OBRing*> local_SSSR; 
    std::map<OBAtom*,unsigned int> ring_shares; 

    atom_stack.push(ring_root); 
    while(!atom_stack.empty()){
      in_rings = 0;
      atom = atom_stack.top();
      atom_stack.pop();
      visited[atom] = true; 
      size++;

      if(atom->GetAtomicNum() != 6)
        hetero = true;

      FOR_RINGS_OF_MOL(r,mol){
        OBRing *obring = &(*r);
        if(obring->IsMember(atom)){
          in_rings++;
          local_SSSR.insert(obring); // use to get the SSSR size into WLNRing 
        }
      }
      ring_shares[atom] = in_rings; 

      FOR_NBORS_OF_ATOM(aiter,atom){
        neighbour = &(*aiter); 
        if(neighbour->IsInRing() && !visited[neighbour]){
          atom_stack.push(neighbour); 
          visited[neighbour] = true;
        }
      }

      if(in_rings > 3)
        branching++;
      else if(in_rings == 3)
        multicyclic++;
      else if (in_rings == 2)
        fuses++;
    }

    if(opt_debug){
      fprintf(stderr,"  SSSR for system:    ");
      for(std::set<OBRing*>::iterator set_iter = local_SSSR.begin();set_iter != local_SSSR.end();set_iter++)
        fprintf(stderr,"%ld(%c) ",(*set_iter)->Size(), (*set_iter)->IsAromatic()?'a':'s');
      fprintf(stderr,"\n");

      fprintf(stderr,"  ring size:          %d\n",size);
      fprintf(stderr,"  fuse points:        %d\n",fuses);
      fprintf(stderr,"  multicyclic points: %d\n",multicyclic);
      fprintf(stderr,"  branching points:   %d\n",branching);
    }

    OBAtom **locant_path = 0;
    if(!branching)
      locant_path = CreateLocantArray(mol,local_SSSR,ring_shares,size);


    if(locant_path){
      free(locant_path);
      locant_path = 0; 
    }
    return wln_ring;   
  }
  

};



/**********************************************************************
                         API FUNCTION
**********************************************************************/


bool WriteWLN(std::string &buffer, OBMol* mol)
{   
 
  WLNGraph wln_graph;
  BabelGraph obabel; 

  unsigned int cyclic = 0;
  bool state = true;
  std::vector<WLNSymbol*> roots; 
  WLNRing* start_ring = 0; 

  FOR_RINGS_OF_MOL(r,mol)
    cyclic++;

  if(!cyclic){
    roots = obabel.ReadBabelNonCyclic(mol,wln_graph);
    if(roots.empty())
      return false; 

    // create an optional wln dotfile
    if (opt_wln2dot)
      WriteGraph(wln_graph);

    if(state)
      state = obabel.ParseNonCyclicWLN(roots,wln_graph,buffer);
  }
  else{

    // get the start ring, and then use that as the jump point
    start_ring = obabel.ReadBabelCyclic(mol->GetAtom(mol->GetSSSR()[0]->_path[0]),buffer,mol ,wln_graph);
    
  }

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


