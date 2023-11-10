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
#include <openbabel/graphsym.h>

using namespace OpenBabel; 

#define REASONABLE 1024

const char *cli_inp;
const char *format; 

// --- options ---
static bool opt_debug = false;
static bool opt_correct = false; 

const char *wln_string;
struct WLNSymbol;
struct WLNEdge; 
struct WLNRing;
struct WLNGraph;
struct ObjectStack;


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

void Fatal(unsigned int pos)
{
  fprintf(stderr, "Fatal: %s\n", wln_string);
  fprintf(stderr, "       ");
  for (unsigned int i = 0; i < pos; i++)
    fprintf(stderr, " ");

  fprintf(stderr, "^\n");

  exit(1);
}


/**********************************************************************
                          STRUCT DEFINTIONS
**********************************************************************/
 

struct WLNEdge{
  WLNSymbol *parent;
  WLNSymbol *child;
  WLNEdge *nxt;
  unsigned int order;
  bool aromatic;

  WLNEdge(){
    parent   = 0;
    child    = 0;
    order    = 0;
    nxt      = 0;
    aromatic = 0;
  }
  ~WLNEdge(){};
};


struct WLNSymbol
{
  unsigned int id;
  unsigned char ch;
  std::string special; // string for element, or ring, if value = '*'
  
  bool aromatic; 
  bool inRing;
  unsigned int allowed_edges;
  unsigned int num_edges;

  WLNSymbol *previous;
  WLNEdge   *bonds; // array of bonds

  // if default needed
  WLNSymbol()
  {
    id = 0;
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
    inRing = 0;
    previous = 0;
    bonds = 0;
    aromatic = 0;
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

  std::map<unsigned char, WLNSymbol *>      locants; 
  std::map<WLNSymbol*,unsigned char>        locants_ch;
  std::vector<std::pair<unsigned char,int>> post_charges; 

  WLNRing(){
    rsize = 0;
    aromatic_atoms = 0;
    adj_matrix = 0;
  }
  ~WLNRing(){
    if(adj_matrix)
      free(adj_matrix);
  };

  bool FillAdjMatrix(){

    aromatic_atoms = 0;
    adj_matrix = (unsigned int*)malloc(sizeof(unsigned int) * (rsize*rsize)); 
    
    if(!adj_matrix){
      fprintf(stderr,"Error: fatal memory allocation on AdjMatrix\n");
      return false;
    }
    
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
        WLNEdge *redge = 0;
        for(redge=rsym->bonds;redge;redge=redge->nxt){
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

  void print_matrix(){
    for (unsigned int i = 0; i< rsize;i++){
      fprintf(stderr,"[ ");
      for (unsigned int j=0; j < rsize;j++)
        fprintf(stderr,"%d ",adj_matrix[i * rsize + j]);
      fprintf(stderr,"]\n");
    }
  }

};

struct WLNBlossom{
  int n; 
  int m; 

  std::vector<int> mate; 
  std::vector<int> p,d,bl;
  std::vector<std::vector<int>> b,g; 

  WLNBlossom(int _n){
    n = _n;
    m = n+n / 2;
    mate.assign(n,-1); 

    b.resize(m);
    p.resize(m);
    d.resize(m);
    bl.resize(m);
    g.assign(m,std::vector<int>(m,-1));
  }

  void add_edge(int u, int v){
    g[u][v] = u;
    g[v][u] = v;
  }

  void match(int u, int v){
    g[u][v] = -1; 
    g[v][u] = -1;
    mate[u] = v;
    mate[v] = u;
  }

  std::vector<int> trace(int x){
    std::vector<int> vx; 
    for(;;){
      while(bl[x] != x)
        x = bl[x];
      
      if(!vx.empty() && vx.back() == x)
        break;
      
      vx.push_back(x);
      x = p[x];
    }

    return vx;
  }

  void contract(int c, int x, int y, std::vector<int> &vx, std::vector<int> &vy){
    b[c].clear();
    int r = vx.back();
    while(!vx.empty() && !vy.empty() && vx.back() == vy.back()){
      r = vx.back();
      vx.pop_back();
      vy.pop_back();
    }

    b[c].push_back(r);
    b[c].insert(b[c].end(),vx.rbegin(),vx.rend());
    b[c].insert(b[c].end(),vy.rbegin(),vy.rend());

    for(int i = 0; i <= c; i++) 
      g[c][i] = g[i][c] = -1;
      
    for(int z : b[c]) {
      bl[z] = c;
      for(int i = 0; i < c; i++) {
        if(g[z][i] != -1) {
          g[c][i] = z;
          g[i][c] = g[i][z];
        }
      }
    }
  }

  std::vector<int> lift(std::vector<int> &vx) {
    std::vector<int> A;
    while(vx.size() >= 2) {
      int z = vx.back(); 
      vx.pop_back();
      if(z < n) {
        A.push_back(z);
        continue;
      }
      int w = vx.back();
      int i = (A.size() % 2 == 0 ? std::find(b[z].begin(), b[z].end(), g[z][w]) - b[z].begin() : 0);
      int j = (A.size() % 2 == 1 ? std::find(b[z].begin(), b[z].end(), g[z][A.back()]) - b[z].begin() : 0);
      int k = b[z].size();
      int dif = (A.size() % 2 == 0 ? i % 2 == 1 : j % 2 == 0) ? 1 : k - 1;
      while(i != j) {
        vx.push_back(b[z][i]);
        i = (i + dif) % k;
      }
      vx.push_back(b[z][i]);
    }
    return A;
  }

  int solve() {
    for(int ans = 0; ; ans++) {
      fill(d.begin(), d.end(), 0);
      std::queue<int> queue;

      for(int i = 0; i < m; i++) 
        bl[i] = i;

      for(int i = 0; i < n; i++) {
        if(mate[i] == -1) {
          queue.push(i);
          p[i] = i;
          d[i] = 1;
        }
      }

      int c = n;
      bool aug = false;
      while(!queue.empty() && !aug) {
        int x = queue.front(); 
        queue.pop();
          
        if(bl[x] != x) 
          continue;

        for(int y = 0; y < c; y++) {
          
          if(bl[y] == y && g[x][y] != -1) {
           
            if(d[y] == 0){
              p[y] = x;
              d[y] = 2;
              p[mate[y]] = y;
              d[mate[y]] = 1;
              queue.push(mate[y]);
            }
            else if(d[y] == 1){
              std::vector<int> vx = trace(x);
              std::vector<int> vy = trace(y);
              if(vx.back() == vy.back()) {
                contract(c, x, y, vx, vy);
                queue.push(c);
                p[c] = p[b[c][0]];
                d[c] = 1;
                c++;
              }
              else {
                aug = true;
                vx.insert(vx.begin(), y);
                vy.insert(vy.begin(), x);
                std::vector<int> A = lift(vx);
                std::vector<int> B = lift(vy);
                A.insert(A.end(), B.rbegin(), B.rend());
                for(int i = 0; i < (int) A.size(); i += 2) {
                  match(A[i], A[i + 1]);
                  if(i + 2 < (int) A.size()) 
                    add_edge(A[i + 1], A[i + 2]);
                }
              }
              break;
            }
          }
        }
      }

      if(!aug) 
        return ans;
    }
  }

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


  // ionic parsing
  std::map<unsigned int,WLNSymbol*> string_positions; 
  std::map<WLNSymbol*,int> charge_additions;

  WLNGraph(){
    edge_count   = 0;
    symbol_count = 0;
    ring_count   = 0;

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

  if(graph.symbol_count >= REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
    return 0;
  }

  if(!ch){
    fprintf(stderr,"Error: null char used to symbol creation\n");
    return 0;
  }

  WLNSymbol *wln = new WLNSymbol;
  wln->id = graph.symbol_count++;
  graph.SYMBOLS[wln->id] = wln;
  wln->ch = ch;
 
  return wln;
}

WLNSymbol* define_hypervalent_element(unsigned char sym, WLNGraph &graph){

  if(!sym){
    fprintf(stderr,"Error: null char used for hypervalent element allocation\n");
    return 0;
  }

  WLNSymbol *new_symbol = 0;
  
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
      fprintf(stderr, "Error: invalid character in special definition switch\n");
      return 0;
  }

  return 0;
}


unsigned int count_children(WLNSymbol *sym){
  WLNEdge *edge = 0; 
  unsigned int count = 0;
  for(edge=sym->bonds;edge;edge=edge->nxt)
    count++;

  if(sym->previous)
    count++;

  if(sym->num_edges == sym->allowed_edges)
    return sym->num_edges;

  return count; 
}

// this one pops based on bond numbers
WLNSymbol *return_object_symbol(ObjectStack &branch_stack){
  WLNSymbol *top = 0;
  while(!branch_stack.empty()){
    top = branch_stack.top().second;
    if(!top)
      return top; // only iterate to the next
    
    if(top->ch == 'Y' && count_children(top)==3)
      branch_stack.pop();
    else if(top->num_edges == top->allowed_edges)
      branch_stack.pop();
    else
      return top;
  }
  return top;
}


bool RaiseBranchingSymbol(WLNSymbol *sym){
  if(!opt_correct)
    return false;

  switch(sym->ch){

    case 'M':
      fprintf(stderr,"Warning: M branches are exceeding 2, raising to N\n");
      sym->allowed_edges++;
      sym->ch = 'N';
      break;

    case 'N':
      if(!sym->inRing){
        fprintf(stderr,"Warning: N branches are exceeding 3, raising to K\n");
        sym->allowed_edges++;
        sym->ch = 'K';
      }
      break;

    default:
      fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", sym->ch,sym->num_edges+1, sym->allowed_edges);
      return 0;
  }
  return true;
}

/**********************************************************************
                          WLNEdge Functions
**********************************************************************/


WLNEdge *AllocateWLNEdge(WLNSymbol *child, WLNSymbol *parent,WLNGraph &graph){

  if(!child || !parent){
    fprintf(stderr,"Error: attempting bond of non-existent symbols - %s|%s is dead\n",child ? "":"child",parent ? "":"parent");
    return 0;
  }
  else if (child == parent){
    fprintf(stderr,"Error: making bond to self is impossible\n");
    return 0;
  }


  graph.edge_count++;
  if(graph.edge_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln symbols - is this reasonable?\n");
    return 0;
  }
  
  if ( ((child->num_edges + 1) > child->allowed_edges) && !RaiseBranchingSymbol(child) ){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges+1, child->allowed_edges);
    return 0;
  }
  
  if ( ((parent->num_edges + 1) > parent->allowed_edges) && !RaiseBranchingSymbol(parent)){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges+1, parent->allowed_edges);
    return 0;
  }

  WLNEdge *edge = new WLNEdge;
  graph.EDGES[graph.edge_count] = edge;
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


WLNEdge *search_edge(WLNSymbol *child, WLNSymbol*parent){
  if(!child || !parent){
    fprintf(stderr,"Error: searching edge on nullptrs\n");
    return 0;
  }
  
  WLNEdge *edge = 0;
  for (edge=parent->bonds;edge;edge = edge->nxt){
    if(edge->child == child)
      return edge;
  }

  for (edge=child->bonds;edge;edge = edge->nxt){
    if(edge->child == parent)
      return edge;
  }

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

  if( (edge->child->num_edges > edge->child->allowed_edges) && !RaiseBranchingSymbol(edge->child)){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->child->ch,edge->child->num_edges, edge->child->allowed_edges);
    return 0;
  }

  if( (edge->parent->num_edges > edge->parent->allowed_edges) && !RaiseBranchingSymbol(edge->parent) ){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", edge->parent->ch,edge->parent->num_edges, edge->parent->allowed_edges);
    return 0;
  }

  return edge;
}

WLNEdge *saturate_edge(WLNEdge *edge,unsigned int n){
  if(!edge){
    fprintf(stderr,"Error: saturating non-existent edge\n");
    return 0;
  }

  if(edge->order < 2)
    return edge;
  
  edge->order -= n; 
  edge->parent->num_edges -= n;
  edge->child->num_edges -= n;

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


WLNEdge* add_methyl(WLNSymbol *head, WLNGraph &graph){

  WLNSymbol *carbon = AllocateWLNSymbol('C',graph);
  WLNSymbol *hydrogen = 0;
  WLNEdge   *edge = 0;
  
  if(carbon)
    carbon->allowed_edges = 4;
  else 
    return 0;
  
  for(unsigned int i=0;i<3;i++){
    hydrogen = AllocateWLNSymbol('H',graph);
    if(hydrogen)
      hydrogen->allowed_edges = 1;
    else
      return 0;
    edge = AllocateWLNEdge(hydrogen,carbon,graph);
    if(!edge)
      return 0;
  }

  WLNEdge *bond = AllocateWLNEdge(carbon,head,graph);
  return bond; 
}


WLNSymbol* create_carbon_chain(WLNSymbol *head,unsigned int size, WLNGraph &graph){

  if (size > REASONABLE){
    fprintf(stderr,"Error: making carbon chain over 1024 long, reasonable molecule?\n");
    return 0;
  }
    
  head->ch = '1';
  head->allowed_edges = 4;

  if(size == 1)
    return head;
  
  WLNEdge *edge = 0;
  WLNSymbol *prev = head;
  for(unsigned int i=0;i<size-1;i++){
    WLNSymbol* carbon = AllocateWLNSymbol('1',graph);
    carbon->allowed_edges = 4; // allows hydrogen resolve
    if(!carbon)
      return 0;
    edge = AllocateWLNEdge(carbon,prev,graph);
    if(!edge)
      return 0;
    prev = carbon;
  } 

  return prev;
}


/* W is going to need to be post resolved, triple bond W will be unsaturated 
with an maximal group dioxo added */
bool add_dioxo(WLNSymbol *head,WLNGraph &graph){

  WLNEdge *edge = 0;
  WLNSymbol *oxygen = 0;
  WLNSymbol *binded_symbol = 0; 

  // 2 options here, either the symbol is binded 'W' -> sym
  // or the symbol is binded sym -> 'W'

  if(head->bonds){
    binded_symbol = head->bonds->child;
    edge = head->bonds; 
  }
  else{
    binded_symbol = head->previous;
    WLNEdge *e = 0; 
    for(e=binded_symbol->bonds;e;e=e->nxt){
      if(e->child == head)
        edge = e;
    }
  }

  if(!binded_symbol || edge->order != 3){
    fprintf(stderr,"Error: dioxo seems to be unbound\n");
    return false; 
  }

  
  head->ch = 'O'; // change the W symbol into the first oxygen
  head->allowed_edges = 2;

  // should still be double bonded. 
  oxygen = AllocateWLNSymbol('O',graph);
  oxygen->allowed_edges = 2;
  edge = saturate_edge(edge,1);


  WLNEdge *sedge = AllocateWLNEdge(oxygen,binded_symbol,graph);
  
  if(binded_symbol->num_edges < binded_symbol->allowed_edges)
    sedge = unsaturate_edge(sedge,1); 
  
  if(binded_symbol->ch == 'N')
    graph.charge_additions[binded_symbol]++;

  if(!edge || !sedge){
    fprintf(stderr,"Error: failure on W post handling\n");
    return false;
  }
    

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
      fprintf(stderr,"Error: resolving methyls performed on invalid symbol: %c\n",target->ch);
      return false;
  }

  return true;
}



/**********************************************************************
                          WLNRing Functions
**********************************************************************/

WLNRing *AllocateWLNRing(WLNGraph &graph)
{
  if(graph.ring_count > REASONABLE){
    fprintf(stderr,"Error: creating more than 1024 wln rings - is this reasonable?\n");
    return 0;
  }

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
  locant->inRing = true;
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
      else if(pos > 3){
        fprintf(stderr,"Error: non-locant links past a two-level tree are unsuitable for this parser\n");
        return false;
      }

      if(opt_debug)
        fprintf(stderr,"  ghost linking %d to parent %c\n",loc_broken,parent);
      
      if(!ring->locants[loc_broken]){
        // bond them in straight away
        allowed_connections[loc_broken] = 3; 
        if(allowed_connections[parent])
          allowed_connections[parent]--;

        WLNSymbol *broken = AllocateWLNSymbol('C',graph);
        broken->inRing = true;
        broken->allowed_edges = 4;
        broken = assign_locant(loc_broken,broken,ring);
        broken_lookup[parent].push_back(loc_broken);
        WLNEdge *edge = AllocateWLNEdge(ring->locants[loc_broken],ring->locants[parent],graph);
        if(!edge)
          return false; 
      }
      else{
        fprintf(stderr,"Error: branching locants are overlapping created elements already in the locant path\n");
        return false;
      }
    }
  }

  return true; 
}

bool set_up_pseudo( WLNRing *ring, WLNGraph &graph,
                            std::vector<unsigned char> &pseudo_locants,
                            std::map<unsigned char,unsigned char> &pseudo_lookback)
{ 

  if(!pseudo_locants.empty()){

    if(pseudo_locants.size() % 2 != 0){
      fprintf(stderr,"Error: uneven pairs read for pseudo locants - ignoring designation\n");
      return false;
    }

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
    if(opt_debug)
      fprintf(stderr,"  calculated size: %c(%d)\n",int_to_locant(local_size),local_size);
  }
  else
    local_size = locant_to_int(size_designator);


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
      curr->inRing = true;
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
      WLNEdge *edge = AllocateWLNEdge(curr,prev,graph);
      if(!edge)
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
      fprintf(stderr,"Error: out of bounds locant access in cyclic builder\n");
      return 0;
    }

    // If we need to catch fuse, we do
    if(i==ring_assignments.size()-1 && pseudo_pairs){
      bool caught = false;
      for(unsigned int s = 1; s<=local_size;s++){
        if(pseudo_lookup[int_to_locant(s)]){
          unsigned char pbind_2 = int_to_locant(s); 
          unsigned char pbind_1 = pseudo_lookup[pbind_2];

          if(opt_debug)
            fprintf(stderr,"  %d  catch fusing: %c <-- %c\n",fuses,pbind_2,pbind_1);
          
          if(!search_edge(ring->locants[pbind_2],ring->locants[pbind_1])){
            WLNEdge *edge = AllocateWLNEdge(ring->locants[pbind_2],ring->locants[pbind_1],graph);
            if(!edge)
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
      WLNEdge *lc = 0;
      for(lc = path->bonds;lc;lc = lc->nxt){
        WLNSymbol *child = lc->child;
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
          fprintf(stderr,"Error: locant path formation is broken in ring definition - '%c(%d)'\n",ring->locants_ch[path],ring->locants_ch[path]);
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
        while(!allowed_connections[bind_2] || bind_2 == bind_1)
          ring_path[path_size-1] = ++bind_2;
        

        if(opt_debug){
          fprintf(stderr,"  %d  fusing (%d): %c <-- %c   [",fuses,comp_size,bind_2,bind_1);
          for (unsigned int a=0;a<path_size;a++)
            fprintf(stderr," %c(%d)",ring_path[a],ring_path[a]);
          fprintf(stderr," ]\n");
        }

        WLNEdge *edge = AllocateWLNEdge(ring->locants[bind_2],ring->locants[bind_1],graph);
        if(!edge)
          return false;

        allowed_connections[bind_1]--;
        if(allowed_connections[bind_2])
          allowed_connections[bind_2]--;

        break;
      }
      else{

        bind_1++; // increase bind_1 
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

    free(ring_path);
    ring_path = 0; 
    fuses++;
  }

  return local_size; 
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
    else
      edge = unsaturate_edge(edge,1);

    if(edge){
      edge->aromatic = 0;
      // edge->child->aromatic = 0;
      // edge->parent->aromatic = 0;
    }
    else
      return false;
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
void FormWLNRing(WLNRing *ring,std::string &block, unsigned int start, WLNGraph &graph,unsigned char spiro_atom='\0'){

  bool warned             = false;  // limit warning messages to console
  bool heterocyclic       = false;  // L|T designator can throw warnings

  unsigned int state_multi          = 0; // 0 - closed, 1 - open multi notation, 2 - expect size denotation
  unsigned int state_pseudo         = 0; 
  unsigned int state_aromatics      = 0;
  unsigned int state_chelate        = 0;
  bool implied_assignment_used      = false; // allows a shorthand if wanted, but not mixing
  
  unsigned int  expected_locants       = 0;
  unsigned int  evaluating_break      = 0;
  unsigned char ring_size_specifier   = '\0';
  unsigned char positional_locant     = '\0';
  unsigned int last_locant_position   = 0;

  std::string special;  

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
  unsigned char ch = *block_str++;

  while(ch){

    switch(ch){

      // specials

      case ' ':

        if(state_multi == 3)
          state_multi = 0;

        if(evaluating_break){
          broken_locants.insert(positional_locant);
          if(state_multi >= 1){
            multicyclic_locants.back() = positional_locant;
            state_multi = 2;
          }
          else if(state_pseudo == 1)
            pseudo_locants.back() = positional_locant;
          else
            bridge_locants[positional_locant] = true;
          

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
        else if (state_pseudo){
          pseudo_locants.back() += 23;
        }
        else if(positional_locant){
          
          // can only be an extension of a positional multiplier for a branch
          if (last_locant_position && last_locant_position == i-1)
            positional_locant += 23;
          else{
            state_aromatics = 1;
            aromaticity.push_back(1);
          }
        }
        else{
          // if unhandled, then it must be an aromatic start
          state_aromatics = 1;
          aromaticity.push_back(1);
        }
        break;

      case '/':
        if(state_aromatics){
          fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
          Fatal(i+start);
        }      
        expected_locants = 2; 
        state_pseudo = 1;
        break; 
      

      // turn this into a look ahead type bias in order to significantly tidy this up
      case '-':{

        // gives us a local working copy
        char *local_arr = new char [strlen(block_str)+1]; 
        memset(local_arr,'\0',sizeof(char) * (strlen(block_str)+1) );
        strcpy(local_arr,block_str);
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

        if(local_arr){
          delete [] local_arr;
          local = 0;
          local_arr = 0;
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
                  last_locant_position = i;
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
                  last_locant_position = i;
                  // no need to move the global pointer here
                }
              }
              else{
                fprintf(stderr,"Error: trying to branch out character without starting point\n");
                Fatal(start+i);  
              }

              break;
            case 1:{
              if(!implied_assignment_used && !positional_locant){
                implied_assignment_used = true;
                positional_locant = 'A';
              }
              // this can only be hypervalent element
  
              if(positional_locant){
                if(spiro_atom){
                  if(positional_locant == spiro_atom){
                    positional_locant++;
                    block_str+=2; 
                    i+=2;
                    break;
                  }
                  else if(ring->locants[positional_locant]){
                    positional_locant++;
                    if(positional_locant == spiro_atom){
                      positional_locant++;
                      block_str+=2; 
                      i+=2;
                      break;
                    }
                  }
                }
                else if(ring->locants[positional_locant])
                  positional_locant++;


                WLNSymbol* new_locant = assign_locant(positional_locant,define_hypervalent_element(special[0],graph),ring);  // elemental definition 
                if(!new_locant)
                  Fatal(i+start);

                graph.string_positions[start+i + 1] = new_locant; // attaches directly

                if(opt_debug)
                  fprintf(stderr,"  assigning hypervalent %c to position %c\n",special[0],positional_locant);
              }
              else{
                fprintf(stderr,"Error: trying to assign element without starting point\n");
                Fatal(start+i);  
              }
              block_str+=2; 
              i+=2; 
              break;
            }
            case 2:{
              if(!implied_assignment_used && !positional_locant){
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

                positional_locant = '\0';
              }
              else{

                if(positional_locant){

                  if(spiro_atom){
                    if(positional_locant == spiro_atom){
                      positional_locant++;
                      block_str+=3; 
                      i+=3;
                      break;
                    }
                    else if(ring->locants[positional_locant]){
                      positional_locant++;
                      if(positional_locant == spiro_atom){
                        positional_locant++;
                        block_str+=3; 
                        i+=3;
                        break;
                      }
                    }
                  }
                  else if(ring->locants[positional_locant])
                    positional_locant++;
                  
                  WLNSymbol* new_locant = assign_locant(positional_locant,define_element(special,graph),ring);  // elemental definition
                  if(!new_locant)
                    Fatal(i+start);

                  graph.string_positions[start+i + 1] = new_locant; // attaches directly to the starting letter

                  if(opt_debug)
                    fprintf(stderr,"  assigning element %s to position %c\n",special.c_str(),positional_locant);
                }
                else{
                  fprintf(stderr,"Error: trying to assign element without starting point\n");
                  Fatal(start+i);  
                }
              }

              block_str+=3; 
              i+=3;              
              break;
            }
            default:
              fprintf(stderr,"Error: %d numerals incased in '-' brackets is unreasonable for WLN to create\n",gap);
              Fatal(start+i);
          }


        }
        else if(i > 0 && block[i-1] == '&')
          state_aromatics = 1;
        else{

          // if there wasnt any other symbol, it must be a notation extender
          evaluating_break = 1; // on a space, number or character, we push to broken_locants

          if(positional_locant){
            if(positional_locant < 128){
              positional_locant = create_relative_position(positional_locant); // i believe breaking modifier will then get removed
              last_locant_position = i;
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
              last_locant_position = i;
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
        // place the minus charges on the last ring seen
        if(ring_components.size() == 1){
          ring->post_charges.push_back({'B',-1});
        }else{
          unsigned int track = 0;
          for (unsigned int rn = 0; rn<ring_components.size()-1;rn++)
            track += ring_components[rn].first;
          ring->post_charges.push_back({int_to_locant(track + 1),-1});  // post is needed as this is before pointer assignment
        }
        
        
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
        if(state_aromatics){
          fprintf(stderr,"Error: character '%c' cannot be in the aromaticity assignment block\n",ch);
          Fatal(i+start);
        }

        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi == 1)
            multicyclic_locants.back() = positional_locant;
          else if(state_pseudo == 1)
            pseudo_locants.back() = positional_locant;
  
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
        if(i == 0 && ch == 'D'){
          state_chelate = 1;
          heterocyclic = true;
          break;
        }

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
          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }

          positional_locant = ch; // use for look back
          expected_locants--;

          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if (positional_locant){

          if(spiro_atom && positional_locant == spiro_atom){
            positional_locant++;
            break;
          }
            
          WLNSymbol *new_locant = 0; 

          switch(ch){
            
            case 'D':
              if(!state_chelate){
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
              }
              // means open chelating bond
              break;

            case 'S':
            case 'P':
              if(!heterocyclic)
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++; 
              
              if(spiro_atom && positional_locant == spiro_atom){
                positional_locant++;
                break;
              }  

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              
              if(ch == 'P')
                new_locant->allowed_edges = 5;
              else
                new_locant->allowed_edges = 6;

              new_locant->inRing = true;
              break;

            case 'Y':
            case 'X':
            case 'K':

              if(!heterocyclic && ch=='K')
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++;

              if(spiro_atom && positional_locant == spiro_atom){
                positional_locant++;
                break;
              }  
              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->allowed_edges = 4;
              new_locant->inRing = true;
              break;

            case 'Z': // treat as NH2
            case 'N':
            case 'B':
              if(!heterocyclic)
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++;

              if(spiro_atom && positional_locant == spiro_atom){
                positional_locant++;
                break;
              }  

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->allowed_edges = 3;
              new_locant->inRing = true;
              break;

            case 'M':
            case 'O':
            case 'V':

              if(!heterocyclic && (ch == 'M' || ch == 'O'))
                warned = true;

              if(ring->locants[positional_locant])
                positional_locant++; 

              new_locant = AllocateWLNSymbol(ch,graph);
              new_locant = assign_locant(positional_locant,new_locant,ring);
              new_locant->allowed_edges = 2;
              new_locant->inRing = true;
              break;

        
            case 'U':{
              // no need to put this in implied, it has to be specified
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

                fprintf(stderr,"triggering- k=%d\n",k);
              }
              else
                unsaturations.push_back({positional_locant,positional_locant+1});
              break;
            }

            // externally bonded to the symbol as a locant symbol
            case 'W':{
              if(!heterocyclic)
                warned = true;

              // the carbon might not be created yet
              if(!ring->locants[positional_locant]){
                new_locant = AllocateWLNSymbol('C',graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->allowed_edges = 2;
                new_locant->inRing = true;
              }
              else
                new_locant = ring->locants[positional_locant];
              
              if(new_locant->ch == 'N')
                new_locant->allowed_edges++;

              WLNSymbol *dioxo = AllocateWLNSymbol('W',graph);
              dioxo->allowed_edges = 3;
              dioxo->inRing = true;
              WLNEdge *e = AllocateWLNEdge(dioxo,new_locant,graph);
              e = unsaturate_edge(e,2);
              if(!e)
                Fatal(start+i);
              
              break;
            }

            // has the effect of unaromatising a bond, remove from edge consideration
            case 'H':
              saturations.push_back({positional_locant,positional_locant+1});
              break;


            default:
              fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
              Fatal(start+i);
          }

          if (opt_debug)
            fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

          graph.string_positions[start+i] = new_locant;
        }
        else{
          if( (i>0 && i<len-1 && block[i-1] == ' ') && (block[i+1] == ' ' || block[i+1] == 'T' || block[i+1] == 'J') ){
            if(ring_components.empty()){
              fprintf(stderr,"Error: assigning bridge locants without a ring\n");
              Fatal(start+i);
            }
            else
              bridge_locants[ch] = true;
          }
          else if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            last_locant_position = i;
          }
          else{
            implied_assignment_used = true;
            positional_locant = 'A';

            if(spiro_atom && positional_locant == spiro_atom){
              positional_locant++;
              break;
            }


            WLNSymbol *new_locant = 0; 

            switch(ch){
              
              case 'D':
              if(!state_chelate){
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
              }
              break;

              case 'S':
              case 'P':
                if(!heterocyclic)
                  warned = true;

                if(ring->locants[positional_locant])
                  positional_locant++; 

                if(spiro_atom && positional_locant == spiro_atom){
                  positional_locant++;
                  break;
                }  

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                
                 if(ch == 'P')
                  new_locant->allowed_edges = 5;
                else
                  new_locant->allowed_edges = 6;

                new_locant->inRing = true;
                break;

              case 'Y':
              case 'X':
              case 'K':
                if(!heterocyclic && ch=='K')
                  warned = true;

                if(ring->locants[positional_locant])
                  positional_locant++;

                if(spiro_atom && positional_locant == spiro_atom){
                  positional_locant++;
                  break;
                }  

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->allowed_edges = 4;
                new_locant->inRing = true;
                break;


              case 'Z':
              case 'N':
              case 'B':
                if(!heterocyclic)
                  warned = true;
                
                if(ring->locants[positional_locant])
                  positional_locant++; 

                if(spiro_atom && positional_locant == spiro_atom){
                  positional_locant++;
                  break;
                }  

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->allowed_edges = 3;
                new_locant->inRing = true;
                break;

              case 'V':
              case 'M':
              case 'O':
                if(!heterocyclic && (ch == 'M' || ch == 'O'))
                  warned = true;

                if(ring->locants[positional_locant])
                  positional_locant++; 

                new_locant = AllocateWLNSymbol(ch,graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->allowed_edges = 2;
                new_locant->inRing = true;
                break;


              case 'U':
                unsaturations.push_back({positional_locant,positional_locant+1});
                break;

              // externally bonded to the symbol as a locant symbol
            case 'W':{
              if(!heterocyclic)
                warned = true;

              // the carbon might not be created yet
              if(!ring->locants[positional_locant]){
                new_locant = AllocateWLNSymbol('C',graph);
                new_locant = assign_locant(positional_locant,new_locant,ring);
                new_locant->allowed_edges = 2;
                new_locant->inRing = true;
              }
              else
                new_locant = ring->locants[positional_locant];

              if(new_locant->ch == 'N')
                new_locant->allowed_edges++;

              WLNSymbol *dioxo = AllocateWLNSymbol('W',graph);
              dioxo->allowed_edges = 3;
              dioxo->inRing = true;
              WLNEdge *e = AllocateWLNEdge(dioxo,new_locant,graph);
              e = unsaturate_edge(e,2);
              if(!e)
                Fatal(start+i);
              
              break;
              }
                // has the effect of unsaturating a bond
              case 'H':
                saturations.push_back({positional_locant,positional_locant+1});
                break;

              default:
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
            }

            if (opt_debug)
              fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

            graph.string_positions[start+i] = new_locant;
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

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }

          positional_locant = ch; // use for look back
          expected_locants--;

          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if( (i>0 && i<len-1 && block[i-1] == ' ') && (block[i+1] == ' ' || block[i+1] == 'T' || block[i+1] == 'J') ){
            if(ring_components.empty()){
              fprintf(stderr,"Error: assigning bridge locants without a ring\n");
              Fatal(start+i);
            }
            else
              bridge_locants[ch] = true;
        }
        else{
          if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            last_locant_position = i;
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

          if(state_multi >= 1 && expected_locants)
            multicyclic_locants.back() = positional_locant;
          else if(state_pseudo == 1 && expected_locants)
            pseudo_locants.back() = positional_locant;
          else
            bridge_locants[positional_locant] = true;

          evaluating_break = 0;
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
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }

          positional_locant = ch; // use for look back
          expected_locants--;
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if( (i>0 && i<len-1 && block[i-1] == ' ') && (block[i+1] == ' ' || block[i+1] == 'T' || block[i+1] == 'J') ){
            if(ring_components.empty()){
              fprintf(stderr,"Error: assigning bridge locants without a ring\n");
              Fatal(start+i);
            }
            else
              bridge_locants[ch] = true;
        }
        else{
          if(i>0 && block[i-1] == ' ' && block[i+1] != 'J'){
            positional_locant = ch;
            last_locant_position = i;
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
        if(state_aromatics)
          state_aromatics = 0;
        
        if(evaluating_break){
          broken_locants.insert(positional_locant);

          if(state_multi >= 1 && expected_locants)
            multicyclic_locants.back() = positional_locant;
          else if(state_pseudo == 1 && expected_locants)
            pseudo_locants.back() = positional_locant;
          else
            bridge_locants[positional_locant] = true;

          evaluating_break = 0;
        }

        if (i == block.size()-1){
          
          if(ring_components.empty()){
            fprintf(stderr,"Error: error in reading ring components, check numerals in ring notation\n");
            Fatal(start+i);
          }

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

          break;
        }
        if(expected_locants){

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else{
            fprintf(stderr,"Error: unhandled locant rule\n");
            Fatal(start+i);
          }

          positional_locant = ch; // use for look back
          expected_locants--;

          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if( (i>0 && i<len-1 && block[i-1] == ' ') && (block[i+1] == ' ' || block[i+1] == 'T' || block[i+1] == 'J') ){
            if(ring_components.empty()){
              fprintf(stderr,"Error: assigning bridge locants without a ring\n");
              Fatal(start+i);
            }
            else
              bridge_locants[ch] = true;
        }
        else{
          if(i>0 && block[i-1] == ' '){
            positional_locant = ch;
            last_locant_position = i;
          }
          else{
            fprintf(stderr,"Error: symbol '%c' is in an unhandled state, please raise issue if this notation is 100%% correct\n",ch);
            Fatal(i+start);
          }
        }
        
        break;

      default:
        break;
    }
    
    i++;
    ch = *(block_str++);
  }

  if(warned)
    fprintf(stderr,"Warning: heterocyclic ring notation required for inter atom assignment, change starting 'L' to 'T'\n");
  

  // debug here
  if (opt_debug){
    
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
  if (!final_size)
    Fatal(start+i);
  
  for (std::pair<unsigned char,int> &post : ring->post_charges)
    graph.charge_additions[ring->locants[post.first]] += post.second;
  
  if(!post_unsaturate(unsaturations,final_size,ring) || !post_saturate(saturations,final_size,ring)){
    fprintf(stderr,"Error: failed on post ring bond (un)/saturation\n");
    Fatal(start+i);
  }
    
}

bool multiply_carbon(WLNSymbol *sym){
  
  WLNEdge   *edge     = 0;
  WLNSymbol *back     = sym->previous; 
  WLNEdge   *fedge    = sym->bonds;
  
  if(!back||!fedge){
    fprintf(stderr,"Error: multiplier carbon must have surrounding symbols, use H to resolve?\n");
    return false;
  }
  
  WLNSymbol *forward  = fedge->child;
  WLNEdge *bedge = 0;
  for(edge=back->bonds;edge;edge=edge->nxt){      
    if(edge->child == sym){
      bedge = edge;
      break;
    }
  }

  if(!forward||!bedge){
    fprintf(stderr,"Error: multiplier carbon must have surrounding symbols, use H to resolve?\n");
    return false;
  }

  unsigned int back_edges = back->allowed_edges - back->num_edges; 
  unsigned int forward_edges = forward->allowed_edges - forward->num_edges;

  // it seems to be the convention
  // that if the forward facing symbol can take a triple bond, we do it
  
  // this cannot be done for alkyl numbers, so.. turn off the chain edges
  if(std::isdigit(back->ch))
    back_edges = 1;
  if(std::isdigit(forward->ch))
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
          sym->num_edges == 1 && graph.charge_additions[sym] == 0)
    {
      edge = sym->bonds; 
      if(edge && edge->order == 1){
        
        while((sym->num_edges < sym->allowed_edges) && 
        (edge->child->num_edges < edge->child->allowed_edges)){
          if(!unsaturate_edge(edge,1))
            return false;
        }
      }
    } 
    else{
      for(edge=sym->bonds;edge;edge = edge->nxt){
        if( (edge->child->ch == 'O' ||
            edge->child->ch ==  'P'  || 
            edge->child->ch ==  'N'  || 
            edge->child->ch ==  'S') &&
            edge->child->num_edges == 1 && graph.charge_additions[edge->child] == 0)
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
bool ExpandWLNSymbols(WLNGraph &graph){

  unsigned int stop = graph.symbol_count;
  // dioxo must be handled before 
  for (unsigned int i=0;i<stop;i++){
    WLNSymbol *sym = graph.SYMBOLS[i];
    if(sym->ch == 'W' && !add_dioxo(sym,graph))
      return false;

    // unsaturated carbons with C
    if(sym->ch == 'c'){
      sym->ch = 'C';
      if(!multiply_carbon(sym))
        return false;
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
          return false;
        break;

      case 'V':{
        sym->ch = 'C';
        sym->allowed_edges = 4;
        
        WLNSymbol *oxygen = AllocateWLNSymbol('O',graph);
        oxygen->allowed_edges = 2;
        
        WLNEdge *e = AllocateWLNEdge(oxygen,sym,graph);
        e = unsaturate_edge(e,1);
        if(!e)
          return false;
        
        break;
      }
      
      default:
        break; // ignore
    }
  }

  return ResolveHangingBonds(graph); 
}


/* backwards search for tentative ionic rule procedures */
unsigned int SearchIonic(const char *wln_ptr, unsigned int len,
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
bool AssignCharges(std::vector<std::pair<unsigned int, int>> &charges,WLNGraph &graph){
  if(charges.empty())
    return true;

  for (std::pair<unsigned int, int> pos_charge : charges){
    WLNSymbol *assignment = graph.string_positions[pos_charge.first - 1]; // reindex as wln 1 is string 0
    if(!assignment){
      fprintf(stderr,"Error: trying to assign ionic charge to unavaliable element, check that character %d is avaliable for assignment\n",pos_charge.first);
      return false;
    }
    else{
      graph.charge_additions[assignment] += pos_charge.second;

      if(opt_debug){
        fprintf(stderr, "  character at position [%d] has the following charge addition - %d\n",pos_charge.first,pos_charge.second);
      }
    }
  }
  return true;
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

    WLNEdge *e = 0 ;
    for(e=top->bonds;e;e=e->nxt){
      if(!ring->locants_ch[e->child])
        continue;

      if(!color[e->child]){
        // find all non-coloured and assign alternative colour
        if(color[top] == 1)
          color[e->child] = 2;
        else
          color[e->child] = 1;
        
        queue.push_front(e->child);
      }
      else if(color[e->child] == color[top])
        return false;
      else if(e->child == top){
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
      if(!visited[v] && ring->adj_matrix[u * ring->rsize + v] > 0){
        path[v] = u;
        if(v == sink)
          return true;

        queue.push_front(v);
      }
    }
  }

  free(visited);
  return false;
}

bool BPMatching(WLNRing *ring, unsigned int u, bool *seen, int *MatchR){
  for(unsigned int v=0;v < ring->rsize;v++){

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
      
  free(seen);
  return true;
}

/* provides methods for `kekulising` wln ring structures, using blossums to maximise pairs */
bool WLNKekulize(WLNGraph &graph){
  for(unsigned int i=0;i<graph.ring_count;i++){
    WLNRing *wring = graph.RINGS[i]; 
    if(wring->aromatic_atoms){

      int   *MatchR = (int*)malloc(sizeof(int) * wring->rsize);
      if(!wring->FillAdjMatrix() || !MatchR){
        fprintf(stderr,"Error: failed to make aromatic matrix\n");
        return false;
      }

      for (unsigned int i=0;i<wring->rsize;i++)
        MatchR[i] = -1;
  

      if(IsBipartite(wring) && !WLNRingBPMaxMatching(wring,MatchR)){
        
        
        free(MatchR);
        MatchR = 0;
        return false;
      }
      else{
        WLNBlossom B(wring->rsize);

        for(unsigned int u = 0; u < wring->rsize; u++){
          for(unsigned int v = 0; v < wring->rsize; v++){
            if(wring->adj_matrix[u * wring->rsize + v] > 0)
              B.add_edge(u,v);
          }
        }

        B.solve();

        for(int i = 0; i < (int)wring->rsize; i++) {
          if(i < B.mate[i]) {
            MatchR[i] = B.mate[i];
            //MatchR[B.mate[i]] = i;
          }
        }
      }

      for(unsigned int i = 0; i<wring->rsize;i++){
        if(MatchR[i] > 0){
          WLNSymbol *f = wring->locants[int_to_locant(i+1)];
          WLNSymbol *s = wring->locants[int_to_locant(MatchR[i]+1)];
          if(f && s){
            WLNEdge *edge = search_edge(f,s);
            if(edge && edge->order == 1)
              edge = unsaturate_edge(edge,1);
            
            if(!edge){
              fprintf(stderr,"Error: failed to unsaturate bond in kekulize\n");
              return false;
            }

            MatchR[MatchR[i]] = 0; // remove from matching
          }
        }
      }
      
      free(MatchR);
      MatchR = 0;
    }
  }

  return true;
}

/**********************************************************************
                         High Level Parser Functions
**********************************************************************/


/* returns the head of the graph, parse all normal notation */
bool ParseWLNString(const char *wln_ptr, WLNGraph &graph) 
{
  
  // keep the memory alive

  if (opt_debug)
    fprintf(stderr, "Parsing WLN notation: %s\n",wln_ptr);

  ObjectStack branch_stack;   // access to both rings and symbols
  branch_stack.reserve(100);  // reasonable size given

  std::vector<std::pair<unsigned int, int>> ionic_charges;
  
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


  unsigned char on_locant = '\0';         // locant tracking
  unsigned int pending_unsaturate = 0;    // 'U' style bonding
  bool j_skips = false;                   // handle skipping of 'J' if in cyclic notation legitimately 

  std::string special;  // special elemental definitions
  
  // allows consumption of notation after block parses
  unsigned int block_start = 0;
  unsigned int block_end = 0;

  unsigned int len = strlen(wln_ptr);
  unsigned int zero_position = SearchIonic(wln_ptr,len,ionic_charges);

  unsigned int i=0;
  unsigned char ch = *wln_ptr;
  
  while(ch)
  {  
    // fprintf(stderr,"char: %c onlocant: %c\n",ch,on_locant);
    
    // if(on_locant && ring)
    //   prev = ring->locants[on_locant];
    
    // dont read any ionic notation
    if(zero_position && zero_position == i)
      break;

    switch (ch)
    {

    case '0': // cannot be lone, must be an addition to another num
      if(pending_J_closure)
        break;

      else if (pending_locant){
        
        if(prev && !prev->inRing)
          graph.charge_additions[prev]++;

        prev = 0;
        on_locant = '0';
        pending_locant = false;
      }
      else{
        fprintf(stderr,"Error: a lone zero mark is not allowed without positive numerals either side\n");
        Fatal(i);
      }
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

        fprintf(stderr,"Error: multipliers are not currently supported\n");
        Fatal(i);

        pending_locant = false;
        on_locant = ch;
      }
      else if(pending_ring_in_ring && pending_inline_ring){
          // onlocant holds the char needed to wrap the ring back, 
        
        if(on_locant != '0'){
          curr = wrap_ring->locants[on_locant];
          if(!curr){
            fprintf(stderr,"Error: cannot access looping ring structure\n");
            Fatal(i);
          }

          if(prev){

            if(prev == branch_stack.branch){
              while(!branch_stack.top().second && !branch_stack.empty())
                branch_stack.pop();
            }

            edge = AllocateWLNEdge(curr,prev,graph);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }
          else
            Fatal(i);

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
      } 
      else{
        on_locant = '\0';
        curr = AllocateWLNSymbol('1',graph);
        curr->allowed_edges = 4;
        graph.string_positions[i] = curr;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(!edge)
            Fatal(i);

          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
        }

        std::string int_sequence;
        int_sequence.push_back(ch);
        while(*(wln_ptr+1)){
          if(!std::isdigit(*(wln_ptr+1)))
            break;

          int_sequence.push_back(*(wln_ptr+1));
          wln_ptr++;
          i++;
        }

        curr = create_carbon_chain(curr,std::stoi(int_sequence),graph);
        if(!curr){
          fprintf(stderr,"Error: error in creating carbon chain, raise algorithm issue\n");
          Fatal(i);
        }

        prev = curr;
      }
      cleared = false;
      break;
    

    case 'Y':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        fprintf(stderr,"Error: '%c' cannot be a locant assignment, please expand [A-W] with &\n",ch);
        Fatal(i);
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 4; // change methyl addition

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(!edge)
            Fatal(i);

          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
        pending_unsaturate = 0;
        prev = curr;
      }
      cleared = false;
      break;

    case 'X':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        fprintf(stderr, "Wiswesser Uncertainities will produce multiple smiles per X entry\n"
                        "since the number of these is at least the size of the ring system\n"
                        "its likely to blow memory allocations, as such they are not supported\n");
        Fatal(i);
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 4;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(!edge)
            Fatal(i);

          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 2;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 1;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 2;
        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
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
        graph.string_positions[i] = curr;

        if(prev){

          if(prev->ch == 'N')
            prev->allowed_edges++;

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          edge = unsaturate_edge(edge,2); // at minimum dioxo must take 3 bonds
          if(pending_unsaturate){
            fprintf(stderr,"Error: a bond unsaturation followed by dioxo is undefined notation\n");
            Fatal(i);
          }
          if(!edge)
            Fatal(i);
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
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

        if(prev){
          if(prev->ch == 'W')
            curr->allowed_edges++;

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }
          
          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 2;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 4;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
          
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
      
        pending_locant = false;
        on_locant = ch;
      }
      else
      { 
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 1;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
          
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 1;

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          
          if(!edge)
            Fatal(i);
          
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
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

        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
        
          if(!edge)
            Fatal(i);
        }         
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        
        if(ch == 'P')
          curr->allowed_edges = 5;
        else
          curr->allowed_edges = 6;
          
        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
      
          if(!edge)
            Fatal(i);
          
        }
        
        branch_stack.push({0,curr});
        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
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
        curr->allowed_edges = 4;

        if(prev && i < len - 1){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }

          prev = curr;
        }

        pending_locant = false;
        on_locant = ch;
      }
      else{
        fprintf(stderr,"Error: locant only symbol used in atomic definition\n");
        Fatal(i);
      }
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
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
        {
          fprintf(stderr, "Error: chelating ring notation started without '-' denotion\n");
          Fatal(i);
        }
        
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else{
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->allowed_edges = 1;
  
        if(prev){

          if(prev == branch_stack.branch){
            while(!branch_stack.top().second && !branch_stack.empty())
              branch_stack.pop();
          }

          edge = AllocateWLNEdge(curr,prev,graph);
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);

          switch(prev->ch){
            case 'Z':
                //graph.charge_additions[prev]++;
                prev->allowed_edges++;
                break;
            default:
              break;
          }
        }

        graph.string_positions[i] = curr;

        if(prev && prev->num_edges < prev->allowed_edges)
          curr = prev;
        else
          prev = return_object_symbol(branch_stack);
        
        if(!prev)
          prev = curr; // failsafe to starting hydrogen
      }
      cleared = false;
      break;

      // ring notation

    case 'J':
      if(pending_J_closure && j_skips)
        break;
      if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
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

          ring->locants[on_locant] = prev;

          // check for an aromaticity bond move?
          if(prev->allowed_edges - prev->num_edges < 2){

            // spiro would not be possible here, check if a double bond can be shifted
            WLNEdge *e = 0;
            WLNSymbol *shift = 0;
            for (e = prev->bonds;e;e = e->nxt){
              if (e->order == 2){
                e = saturate_edge(e,1);
                if(!e)
                  Fatal(i);

                shift = e->child;
                break;
              }
            }
            unsigned char next_loc = branch_stack.ring->locants_ch[shift]+1;
            if(!next_loc)
              next_loc = 'A'; // must of done the full loop

            e = search_edge(branch_stack.ring->locants[next_loc],shift);
            e = unsaturate_edge(e,1);
            if(!e)
              Fatal(i);
          }
          
          FormWLNRing(ring,r_notation,block_start,graph,on_locant);
        }
        else
          FormWLNRing(ring,r_notation,block_start,graph);
        

        if(pending_ring_in_ring && !wrap_ring)
          wrap_ring = ring; // instant back access


        branch_stack.push({ring,0});

        block_start = 0;
        block_end = 0;

        // does the incoming locant check

        if(pending_spiro)
          pending_spiro = false;
        else if (prev && on_locant && on_locant != '0')
        {
          if (ring->locants[on_locant]){
            edge = AllocateWLNEdge(ring->locants[on_locant],prev,graph);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }   
          else
          {
            fprintf(stderr, "Error: attaching inline ring with out of bounds locant assignment\n");
            Fatal(i);
          }
        }

        on_locant = '\0';
        pending_J_closure = false;
      }
      cleared = false;
      break;

    case 'L':
    case 'T':
      if (pending_J_closure)
        break;
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr){
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
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
        {
          fprintf(stderr, "Error: ring notation started without '-' denotion\n");
          Fatal(i);
        }
        
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
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
        FormWLNRing(ring,r_notation,i,graph);
        branch_stack.push({ring,0});

        curr = ring->locants['A'];
        if(prev){
          edge = AllocateWLNEdge(curr,prev,graph);
          
          if(pending_unsaturate){
            edge = unsaturate_edge(edge,pending_unsaturate);
            pending_unsaturate = 0;
          }
          if(!edge)
            Fatal(i);
        }

        graph.string_positions[i] = curr;
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
            fprintf(stderr,"Error: accessing locants out of range\n");
            Fatal(i);
          }
          
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else if(cleared){
        fprintf(stderr,"Error: floating double bond after ionic clear\n");
        Fatal(i);
      }
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

      if(!branch_stack.empty() && !pending_inline_ring)
        branch_stack.pop_to_ring();
      
      if( (i < len - 1 && wln_string[i+1] == '&') || branch_stack.ring){
        pending_locant = true;

        // single letter methyl branches
        if(on_locant && !pending_inline_ring){
          if(!add_methyl(branch_stack.ring->locants[on_locant],graph)){
            fprintf(stderr,"Error: could not attach implied methyl to ring\n");
            Fatal(i);
          }
          on_locant = '\0';
        }
        
      }
      else if (!opt_correct){
        fprintf(stderr,"Error: space used outside ring and ionic notation\n");
        Fatal(i);
      }
      // only burn the stacks now on ionic clearance
      break;


    case '&':
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
        curr = 0;
        pending_locant = false;
        cleared = true;
        branch_stack.clear_all(); // burn stack
      }
      else if(on_locant){
        if(curr && curr == ring->locants[on_locant]){
          on_locant += 23;
          curr = ring->locants[on_locant];  
          if(!curr){
            fprintf(stderr,"Error: could not fetch expanded locant position - out of range\n");
            Fatal(i);
          }
          prev = curr;
        }
      }
      else if (i < len-1 && wln_string[i+1] == ' '){
        // this must be a ring pop, no matter what
        
        if(branch_stack.empty() || !branch_stack.ring){
          fprintf(stderr,"Error: '&' followed by a space indicates a ring pop, are there any rings?\n");
          Fatal(i); 
        }
        else if(!branch_stack.empty()){
          branch_stack.pop_to_ring();
          branch_stack.pop(); // pop whats open
          ring = branch_stack.ring; // assign the previous ring
          
          prev = return_object_symbol(branch_stack);
          if(!prev)
            prev = branch_stack.branch;
        }
      }
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
                    Fatal(i);

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
                    Fatal(i);

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
      else{
          fprintf(stderr,"Error: popping too many rings|symbols, check '&' count\n");
          Fatal(i);
        }
      break;


    case '-':{
      if (pending_J_closure)
        break;

      else if (pending_inline_ring)
      { 
        if(pending_ring_in_ring){
          // onlocant holds the char needed to wrap the ring back, 
          curr = wrap_ring->locants[on_locant];
          if(!curr){
            fprintf(stderr,"Error: cannot access looping ring structure\n");
            Fatal(i);
          }

          if(prev){
            
            if(prev == branch_stack.branch){
              while(!branch_stack.top().second && !branch_stack.empty())
                branch_stack.pop();
            }

            edge = AllocateWLNEdge(curr,prev,graph);
            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);
          }
          else
            Fatal(i);
    
          // last notation is not neccessary
          while(wln_ptr){
            if(*wln_ptr == 'J')
              break;
            wln_ptr++;
            i++;
          }

          on_locant = '\0';
          pending_ring_in_ring = false;
          pending_inline_ring = false;
        }
        else{
          fprintf(stderr, "Error: only one pending ring can be active, check closures\n");
          Fatal(i);
        }
      }
#ifdef HIGH_LEVEL_ASSIGNMENTS
      else if(on_locant){
        if(prev->ch < 128){
          prev->ch += 128;
        }
        else
          prev->ch += 46; 
      }
#endif
      else{

        // ahh variable length array.. gotcha
        char *local_arr  = new char [strlen(wln_ptr)+1]; 
        memset(local_arr,'\0',strlen(wln_ptr)+1 * sizeof(char));
        strcpy(local_arr,wln_ptr);
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

        if(local_arr){
          delete [] local_arr;
          local = 0;
          local_arr = 0;
        }
        
        if(!found_next){
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

          if(gap == 1){
            curr = define_hypervalent_element(special[0],graph);
            if(!curr)
              Fatal(i);
            special.clear();
          }
          else if(gap == 2){
            curr = define_element(special,graph);
            if(!curr)
              Fatal(i); 
            
            if(on_locant == '0'){
              graph.charge_additions[curr]++;
              //on_locant = '\0';
            }
            special.clear();
          }
          else{
            fprintf(stderr,"Error: special '-' must be either 1 or 2 symbols - %d seen\n",gap);
            Fatal(i);
          }

          if(prev){
            if(!gap && ring)
              edge = AllocateWLNEdge(ring->locants[prev->ch],prev,graph);
            else{
              if(prev == branch_stack.branch){
                while(!branch_stack.top().second && !branch_stack.empty())
                  branch_stack.pop();
              }
              edge = AllocateWLNEdge(curr,prev,graph);
            }
              

            if(pending_unsaturate){
              edge = unsaturate_edge(edge,pending_unsaturate);
              pending_unsaturate = 0;
            }
            if(!edge)
              Fatal(i);            
            
          }
          on_locant = '\0';
          branch_stack.push({0,curr});

          i+= gap+1;
          wln_ptr+= gap+1;

          graph.string_positions[i-gap] = curr;
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
      else{
        fprintf(stderr,"Error: multipliers are not currently supported\n");
        Fatal(i);
      }
      
      
      break;

    default:
      fprintf(stderr, "Error: unallowed character! - alphabet: [A-Z][0-1][&-/' ']\n");
      Fatal(i);
    }

    i++;
    ch = *(++wln_ptr);
  }

  // single letter methyl branches
  if(on_locant && on_locant != '0' && !pending_inline_ring && !branch_stack.empty()){
    if(!add_methyl(branch_stack.ring->locants[on_locant],graph)){
      fprintf(stderr,"Error: could not attach implied methyl to ring\n");
      Fatal(i);
    }
    on_locant = 0;
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


  if(!AssignCharges(ionic_charges,graph))
    Fatal(len);

    // use this for recursion on multipliers
  return true;
}


/* dump wln tree to a dotvis file */
void WLNDumpToDot(FILE *fp, WLNGraph &graph)
{  
  fprintf(fp, "digraph WLNdigraph {\n");
  fprintf(fp, "  rankdir = LR;\n");
  for (unsigned int i=0; i< graph.symbol_count;i++)
  {
    WLNSymbol *node = graph.SYMBOLS[i];
    if(!node)
      continue;

    fprintf(fp, "  %d", node->id);
    if (node->ch == '*')
      fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
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
  
      
    WLNEdge *edge = 0;
    for (edge = node->bonds;edge;edge = edge->nxt){

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
  fprintf(stderr,"Dumping wln graph to %s:\n",filename);
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



/**********************************************************************
                         BABEL Mol Functions
**********************************************************************/


// holds all the functions for WLN graph conversion, mol object is assumed ALIVE AT ALL TIMES
// uses old NM functions from previous methods: Copyright (C) NextMove Software 2019-present
struct BabelGraph{

  std::map<unsigned int,OBAtom*> babel_atom_lookup;

  BabelGraph(){};
  ~BabelGraph(){};


  OBAtom* NMOBMolNewAtom(OBMol* mol, unsigned int elem,int charge,unsigned int hcount)
  {
    OBAtom* result = mol->NewAtom();
    result->SetAtomicNum(elem);
    result->SetImplicitHCount(hcount);
    result->SetFormalCharge(charge);
    return result;
  }


  OBBond* NMOBMolNewBond( OBMol* mol,
                          OBAtom* s,
                          OBAtom* e,
                          unsigned int order)
  {
    
    OBBond* bptr = 0; 
    if(!s || !e){
      fprintf(stderr,"Error: could not find atoms in bond, bond creation impossible\n");
      return bptr;
    }

    if (!mol->AddBond(s->GetIdx(), e->GetIdx(), order)){
      fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n",s->GetIdx(),e->GetIdx());
      return bptr;
    }
        
    bptr = mol->GetBond(mol->NumBonds() - 1);
    return bptr;
  }


  bool NMOBSanitizeMol(OBMol* mol)
  {
    // WLN has no inherent stereochemistry, this can be a flag but should be off by default
    mol->SetChiralityPerceived(true);
    mol->SetAromaticPerceived(false);

#ifdef REPLACED
    if(!OBKekulize(mol)){
      fprintf(stderr,"Error: failed to kekulize mol\n");
      if(!opt_debug) // if we cant kekulise lets see why
        return false; 
    }
#endif

    mol->DeleteHydrogens();
    return true;
  }


  bool ConvertFromWLN(OBMol* mol,WLNGraph &graph){

    if(opt_debug)
      fprintf(stderr,"Converting wln to obabel mol object: \n");

    // set up atoms
    for (unsigned int i=0; i<graph.symbol_count;i++){
      WLNSymbol *sym = graph.SYMBOLS[i];
      OBAtom *atom = 0;

      int charge = 0; 
      unsigned int atomic_num = 0;
      unsigned int hcount = 0;

      switch(sym->ch){

        case 'H':
          atomic_num = 1;
          hcount = 0;
          break; 

        case 'B':
          atomic_num = 5;
          break;

        case '1': 
        case 'C':
          atomic_num = 6;
          while(sym->num_edges < sym->allowed_edges){
            hcount++;
            sym->num_edges++;
          }
          break;
        
        case 'X':
          atomic_num = 6;
          break;

        case 'Y':{
          atomic_num = 6;
          WLNEdge *e = 0;
          unsigned int orders = 0; 
          if(!sym->inRing){
            for(e = sym->bonds;e;e=e->nxt)
              orders += e->order;
            
            if(sym->previous){
              e = search_edge(sym,sym->previous);
              orders += e->order;
            }
              

            if(orders < 4)
              hcount = 1;
          }
          break;
        }

        case 'N':
          atomic_num = 7;
          if(sym->inRing)
            sym->allowed_edges = 3;

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
          if(sym->num_edges == 1)
            charge = -1;
          if(!sym->num_edges)
            charge = -2; 
          break;

        case 'Q':
          if(sym->num_edges == 0)
            charge = -1;
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
          while(sym->num_edges % 2 == 0){ // 3 and 5 valence
            hcount++;
            sym->num_edges++;
          }
          break;
        
        case 'S':
          atomic_num = 16;
          while(sym->num_edges % 2 != 0){ // 2 and 6 valence
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

      // ionic notation - overrides any given formal charge
      if(graph.charge_additions[sym]){
        charge = graph.charge_additions[sym];
        if(charge < 0 && hcount)
          hcount --; // let the charges relax the hydrogens 
      }
        
      atom = NMOBMolNewAtom(mol,atomic_num,charge,hcount);
      if(!atom){
        fprintf(stderr,"Error: formation of obabel atom object\n");
        return false;
      }

      babel_atom_lookup[sym->id] = atom;
    }

    // create edges 
    std::map<WLNEdge*, OBBond*> bond_map; 

    for(unsigned int i=0;i<graph.symbol_count;i++){
      WLNSymbol *parent = graph.SYMBOLS[i];
      WLNEdge *e = 0;
      if(parent->bonds){
        for (e = parent->bonds;e;e = e->nxt){
          WLNSymbol *child = e->child;
          unsigned int bond_order = e->order;  
          OBBond *bptr = NMOBMolNewBond(mol,babel_atom_lookup[parent->id],babel_atom_lookup[child->id],bond_order);
          if(!bptr)  
            return false;
          bond_map[e] = bptr; 
        }
      }
    }

    return true;
  }

};


/**********************************************************************
                         API FUNCTION
**********************************************************************/

bool ReadWLN(const char *ptr, OBMol* mol)
{   
  if(!ptr){
    fprintf(stderr,"Error: could not read wln string pointer\n");
    return false;
  }
  else 
    wln_string = ptr; 

  WLNGraph wln_graph;
  BabelGraph obabel; 

  if(!ParseWLNString(ptr,wln_graph))
    return false;

  if (opt_debug)
    WriteGraph(wln_graph,"wln-graph.dot");
  
    // needs to be this order to allow K to take the methyl groups
  if(!WLNKekulize(wln_graph))
    return false;

  if(!ExpandWLNSymbols(wln_graph))
    return false;

  if(!obabel.ConvertFromWLN(mol,wln_graph))
    return false;

  if(!obabel.NMOBSanitizeMol(mol))
    return false;

  return true;
}

static void DisplayUsage()
{
  fprintf(stderr, "readwln <options> -o<format> -s <input (escaped)>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, " -c                   allow run-time spelling correction where possible\n");
  fprintf(stderr, " -d                   print debug messages to stderr\n");
  fprintf(stderr, " -h                   show the help for executable usage\n");
  fprintf(stderr, " -o                   choose output format (-osmi, -oinchi, -ocan)\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser reads and evaluates wiswesser\n"
                  " line notation (wln), the parser is native\n"
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
      switch (ptr[1])
      {

      case 'c':
        opt_correct = true;
        break;

      case 'd':
        opt_debug = true;
        break;

      case 'h':
        DisplayHelp();

      case 'o':
        if (!strcmp(ptr, "-osmi"))
        {
          format = "smi";
          break;
        }
        else if (!strcmp(ptr, "-oinchi"))
        {
          format = "inchi";
          break;
        }
        else if (!strcmp(ptr, "-ocan"))
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
    fprintf(stderr,"Error: no output format selected\n");
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
  if(!ReadWLN(cli_inp,&mol))
    return 1;
  
  OBConversion conv;
  conv.AddOption("h",OBConversion::OUTOPTIONS);
  conv.SetOutFormat(format);

  res = conv.WriteString(&mol);
  std::cout << res;
  return 0;
}