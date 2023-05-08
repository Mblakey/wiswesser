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

// should have the same structs, we just go from the babel graph
// to the wln graph, and then write back the notation


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
const char *inp_format; // only needed for outside babel build
const char *dotfile;

// --- options ---
static bool opt_wln2dot = false;
static bool opt_debug = false;

// --- globals ---
const char *wln;
struct WLNSymbol;
struct WLNEdge; 
struct WLNRing;
struct ObjectStack;


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
  RING = 1,     
  SPECIAL = 2  // for now this is only going to be the pi bond
};

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

void Fatal(unsigned int pos)
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
    bonds = 0;
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
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
    exit(0);
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
  return edge;
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
};


/* graph now interprets all babel mol bonds and atoms, created wln graph */
struct BabelGraph{

    BabelGraph(){};
    ~BabelGraph(){};


    bool BuildWLNGraph(OpenBabel::OBMol *mol,WLNGraph &wln_graph){

      // we can do the reverse to create the babel graph

      OpenBabel::OBAtom* atom = 0;
      FOR_ATOMS_OF_MOL(atom,mol){
        
        if(opt_debug)
          fprintf(stderr,"  created: atom[%d] - atomic num(%d), charge(%d)\n",atom->GetIdx(),atom->GetAtomicNum(),atom->GetFormalCharge());

        WLNSymbol *node = 0;
        switch(atom->GetAtomicNum()){
          case 1:
            node = AllocateWLNSymbol('H');
            break; 

          case 5:
            node = AllocateWLNSymbol('B');
            break;

          case 6:
            node = AllocateWLNSymbol('C');
            break;

          case 7:
            node = AllocateWLNSymbol('N');
            break;
          
          case 8:
            node = AllocateWLNSymbol('O');
            break;
          
          case 9:
            node = AllocateWLNSymbol('F');
            break;

          case 15:
            node = AllocateWLNSymbol('P');
            break;

          case 16:
            node = AllocateWLNSymbol('S');
            break;

          case 17:
            node = AllocateWLNSymbol('G');
            break;

          case 35:
            node = AllocateWLNSymbol('E');
            break;

          case 53:
            node = AllocateWLNSymbol('I');
            break;

          default:
            fprintf(stderr,"Error: unhandled element for WLNSymbol formation\n");
            return false;

        }

        if(node)
          node->set_edge_and_type(atom->GetTotalValence(),STANDARD); // allow smiles to dictate

        if(!wln_graph.root)
          wln_graph.root = node; 
        
      }


      // set up the bonds
      OpenBabel::OBBond* bond = 0;
      FOR_BONDS_OF_MOL(bond,mol){

        unsigned int b_idx = bond->GetBeginAtomIdx();
        unsigned int e_idx = bond->GetEndAtomIdx();
        unsigned int order = bond->GetBondOrder();

        if(opt_debug)
          fprintf(stderr,"  bonding: atoms %3d --> %3d [%d]\n",b_idx,e_idx,order);

        WLNEdge *edge = 0; 
        edge = AllocateWLNEdge(symbol_lookup[b_idx],symbol_lookup[e_idx]);
        if(order > 1){
          for (unsigned int i=1;i<order;i++)
            edge = unsaturate_edge(edge,1);
        }

      }

      return true;
    }



};

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

bool WriteGraph(){
  fprintf(stderr,"Dumping wln graph to wln-write-graph.dot:\n");
  FILE *fp = 0;
  fp = fopen("wln-write-graph.dot", "w");
  if (!fp)
  {
    fprintf(stderr, "Error: could not create dump .dot file\n");
    fclose(fp);
    return false;
  }
  else
  {
    WLNDumpToDot(fp);
    fclose(fp);
  }
  fprintf(stderr,"  dumped\n");
  return true;
}



bool WriteWLN(const char *inp,OpenBabel::OBMol *mol){

  return true;
}


static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser reads and evaluates smiles/inchi\n"
                  " and converts to wisswesser line notation (wln)\n"
        ); 
  exit(1);
}

static void DisplayUsage()
{
  fprintf(stderr, "writewln <options> < input (escaped) >\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -i<format>                    choose the input format for reading\n");
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

      case 'd':
        opt_debug = true;
        break;

      case 'h':
        DisplayHelp();
        break;

      case 'w':
        opt_wln2dot = true;
        break;

      case '-':

        if (!strcmp(ptr, "--debug"))
          opt_debug = true;
        else if (!strcmp(ptr, "--help"))
          DisplayHelp();
        else if (!strcmp(ptr, "--wln2dot"))
          opt_wln2dot = true;
  
      break;

      case 'i':
        if (!strcmp(ptr, "-ismi"))
          inp_format = "smi";
        else if (!strcmp(ptr, "-iinchi"))
          inp_format = "inchi";
        break;
        
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

  if(!inp_format){
    fprintf(stderr,"Error: please select an input format for string\n");
    DisplayUsage();
  }

  return;
}



int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);
  
  OpenBabel::OBMol* mol = new OpenBabel::OBMol;

  OpenBabel::OBConversion conv;
  conv.SetInFormat(inp_format);
  conv.ReadString(mol,cli_inp);
  
  BabelGraph obabel; 
  WLNGraph wln_graph;

  obabel.BuildWLNGraph(mol,wln_graph);

  WriteGraph();

  delete mol; 
  return 0;
}