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

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <set>
#include <deque>
#include <string>
#include <strings.h>
#include <vector>
#include <map>

#include <utility> // std::pair
#include <iterator>

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

#include "huffman.h"
#include "openbabel/parsmart.h"
#include "openbabel/stereo/stereo.h"
#include "parser.h"

using namespace OpenBabel; 

#define STRUCT_COUNT 1024
#define MAX_EDGES 8
#define AMPERSAND_EXPAND 23
#define BROKEN_TREE_LIMIT 6

// --- DEV OPTIONS  ---
#define OPT_CORRECT 0


const char *wln_input;
struct WLNSymbol;
struct WLNEdge; 
struct WLNRing;
struct WLNGraph;
struct ObjectStack;

#define INT_TO_LOCANT(X) (X+64)
#define LOCANT_TO_INT(X) (X-64)

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


bool Fatal(unsigned int pos, const char *message)
{ 

#if ERRORS == 1
  fprintf(stderr,"%s\n",message);
  fprintf(stderr, "Fatal: %s\n", wln_input);
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
  unsigned int order;
  bool aromatic;
  unsigned int stereo; // 0- none, 1 - descend, 2- accend

  WLNEdge(){
    parent   = 0;
    child    = 0;
    order    = 0;
    aromatic = 0;
    stereo = 0; 
  }
  ~WLNEdge(){};
};


struct WLNSymbol
{
  unsigned short int id;
  unsigned int str_position; 
  short int charge;
  unsigned int explicit_H; // if explicit, take the value, else we calculate the fill valence?

  unsigned char ch;
  std::string special; // string for element, or ring, if value = '*'
  
  bool aromatic;
  bool spiro;
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
    spiro = 0; 
    charge = 0;
    explicit_H = 0; 
    str_position = 0; 
    barr_n = 0; 
    parr_n = 0; 
  }
  ~WLNSymbol(){};

  void add_special(const char* ptr, unsigned int s, unsigned int e)
  {
    for (unsigned int i = s; i <= e; i++)
      special.push_back(ptr[i]);
  }

};

struct WLNRing
{
  unsigned int rsize;
  unsigned int aromatic_atoms;
  unsigned int *adj_matrix; 
  
  // quick ring perception
  std::vector<unsigned char> assignment_locants;
  std::vector<unsigned int>   assignment_digits;
  
  std::map<unsigned int, WLNSymbol *>     locants; // int allows effectively infinite broken positions 
  std::map<WLNSymbol*,unsigned int>       locants_ch;
  std::map<WLNSymbol*,unsigned int>       position_offset; // gives hetero position in the string

  bool spiro;
  WLNEdge *macro_return;

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
    macro_return = 0; 
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
      unsigned char loc_a = INT_TO_LOCANT(i+1);
      WLNSymbol *rsym = locants[loc_a]; 
      if(rsym->ch == 'S' || (rsym->ch == 'N' && rsym->charge < 0)) // for now lets see
        continue;
      
      if(rsym->aromatic && rsym->num_edges < rsym->allowed_edges){
        
        for(unsigned int ei = 0;ei < rsym->barr_n;ei++){
          WLNEdge *redge = &rsym->bond_array[ei];
          WLNSymbol *csym = redge->child; 

          if(csym->ch == 'S' || redge->order > 1 || (csym->ch == 'N' && csym->charge < 0))
            continue;
        
          if(csym->aromatic && redge->aromatic && csym->num_edges < csym->allowed_edges){
            unsigned char loc_b = locants_ch[csym];
            unsigned int c = LOCANT_TO_INT(loc_b-1);
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


// see notes in build cyclic for usuage. 
struct LocantPos{
  bool active;
  int  allowed_connections; 
  WLNSymbol *locant; // look up character based on the ring_map
  WLNSymbol *broken_a; 
  WLNSymbol *broken_b; 

  LocantPos(){
    active = false; 
    locant = 0; 
    broken_a = 0; 
    broken_b = 0; 
    allowed_connections = 0; 
  }
}; 


struct LocantPair{
  unsigned int first; 
  unsigned int second;
  unsigned int stereo; 
  LocantPair(){
    first = 0; 
    second = 0;
    stereo = 0;  
  }

  LocantPair(unsigned char f, unsigned char s){
    first = f; 
    second = s; 
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
    m = +n+ (n / 2);
    mate.assign(n,-1); 

    b.resize(m);
    p.resize(m);
    d.resize(m);
    bl.resize(m);
    g.assign(m,std::vector<int>(m,-1));
  }

  void add_edge(int u, int v){
    if(u >= g.size())
      return;
    if(v >= g[u].size())
      return;
    
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

      if(z >= b.size())
        return {};

      int i = (A.size() % 2 == 0 ? std::find(b[z].begin(), b[z].end(), g[z][w]) - b[z].begin() : 0);
      int j = (A.size() % 2 == 1 ? std::find(b[z].begin(), b[z].end(), g[z][A.back()]) - b[z].begin() : 0);
      int k = b[z].size();
      int dif = (A.size() % 2 == 0 ? i % 2 == 1 : j % 2 == 0) ? 1 : k - 1;
      
      unsigned int safety = 10000;
      while(i != j) {
        vx.push_back(b[z][i]);
        i = (i + dif) % k;
        safety--;
        if(!safety){
          return {};
        }
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

  unsigned int symbol_count;
  unsigned int ring_count;

  WLNSymbol *SYMBOLS[STRUCT_COUNT];
  WLNRing   *RINGS  [STRUCT_COUNT];
  
  WLNGraph(){
    symbol_count    = 0;
    ring_count      = 0;
    root            = 0;

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


// needs to be able to hold both a WLNSymbol and WLNRing for branch returns,
// first is the WLNRing, second is the branching WLN symbol 
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
    case 'I':
    case 'Q':
    case 'Z':
      return true;
  }
  return false;
}

bool IsBranching(WLNSymbol *symbol){
  switch (symbol->ch) {
    case 'S':
    case 'P':
    case 'Y':
    case 'X':
    case 'K':
    case 'N':
    case 'B':
    case '*':
      return true;

    case 'G':
    case 'F':
    case 'I':
      if(symbol->allowed_edges > 1)
        return true;
    case 'O':
      if(symbol->allowed_edges > 2)
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

    case 'N':
      new_symbol = AllocateWLNSymbol(sym,graph);
      if(new_symbol)
        new_symbol->allowed_edges = 6; // tautomeric possibility 
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
        new_symbol->allowed_edges = 8;         // allows FCl6, Ive seen [CL](=O)(=O)(=O)(O-)
      break;
  }
  
  return new_symbol;
}

/* allocate new or override exisiting node*/
WLNSymbol* define_element(unsigned char special_1, unsigned char special_2, WLNGraph &graph){
    
  WLNSymbol *created_wln = 0;
  
  switch (special_1){

    case 'A':
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
        case 'N':
        case 'R':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'K':
      switch(special_2){
        case 'R':
        case 'A':
          created_wln = AllocateWLNSymbol('*',graph);
          break;

        default:
          return (WLNSymbol *)0;
      }
      break;
      

    case 'L':
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
        case 'O':
        case 'G':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'P':
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      switch(special_2){
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
      if(special_2 == 'R')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;
      
      break;

    case 'V':
      if (special_2 == 'A')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;
      
      break;
    
    case 'W':
      if(special_2 == 'T')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;
      
      break;
    

    case 'X':
      if (special_2 == 'E')
        created_wln = AllocateWLNSymbol('*',graph);
      else
        return (WLNSymbol *)0;

      break;

    case 'Y':
      switch(special_2){
        case 'B':
        case 'T':
          created_wln = AllocateWLNSymbol('*',graph);
          break;
          
        default:
          return (WLNSymbol *)0;
      }
      break;

    case 'Z':
      switch(special_2){
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

  created_wln->special+= special_1;
  created_wln->special+= special_2;
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
  
    case 'W':
      if(special[1] == 'T')
        return 74;
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
    else if(top->num_edges >= top->allowed_edges)
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
    fprintf(stderr,"Error: binding invalid nodes\n"); 
    return 0;
  }
  
  // dont make the same bond twice
  for(unsigned int i=0;i<parent->barr_n;i++)
    if(parent->bond_array[i].child == child)
      return true; 

  
  if(parent->barr_n >= MAX_EDGES){
    fprintf(stderr,"Error: creating more %d bonds on a singular symbol - is this reasonable?\n",MAX_EDGES);
    return 0;
  }
  
  if ( ((child->num_edges + 1) > child->allowed_edges)){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges+1, child->allowed_edges);
    return 0;
  }
  
  if ( ((parent->num_edges + 1) > parent->allowed_edges)){
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges+1, parent->allowed_edges);
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
#if ERRORS == 1
  if(edge->order == 3){
    fprintf(stderr,"Error: attempting a quadruple bond - not allowed\n"); 
    return false;
  }
#endif 

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
  WLNSymbol *carbon = AllocateWLNSymbol('#',graph);
  carbon->special = '1'; 
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
  
  if(binded_symbol->ch == 'N' && binded_symbol->allowed_edges == 4){
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
      
      while( (target->num_edges + target->explicit_H) < target->allowed_edges){
        if(!add_methyl(target,graph))
          return false;
      }
      break;

    case 'Y':
      while(count_children(target) < 3 && target->num_edges < target->allowed_edges){
        if(!add_methyl(target,graph))
          return false;
      }
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



bool assign_locant(unsigned char loc,WLNSymbol *locant, WLNRing *ring){
  if(!locant)
    return 0;
  ring->locants[loc] = locant; 
  ring->locants_ch[locant] = loc;
  locant->inRing = ring;
  return 1; 
}  


/*
 The broken locant conventions are as follows, the first broken position is given a value of 
 128, its extremely unlikly that a reasonable chemical will have 128 non broken cyclic atoms:

 limiting the tree depth to 2, for now, which gives 6 values per parent symbol
 - Each normal locant can have two broken locants, one at - and one at -&, e.g B- and B-&
 
 Access is given by a binary tree sequence in an array starting at 128 and indexed by locant
 gets you the start node: examples: 
  
      B
     / \
    B- B-&
   /    \
  B--   B--&
 
 128 + (O * 6) = 'A' node, 128 = A-, 129 = A--, 130 = A--&, 131 = A-&, 132 = A-&-, 133 = A-&&
 128 + (1 * 6) = 'B' node, 134 = B-, 135 = B--, 136 = B--&, 137 = B-&, 138 = B-&-, 139 = B-&&
 128 + (2 * 6) = 'C' node, 140... 
*/
WLNSymbol* assign_broken_locant(unsigned int loc,WLNSymbol *locant, WLNRing *ring){
  if(!locant)
    return 0;
  ring->locants[loc] = locant; 
  ring->locants_ch[locant] = loc;
  locant->inRing = ring;
  return locant; 
}  


#if BROKEN
bool set_up_broken( WLNRing *ring, WLNGraph &graph,
                    std::set<unsigned char> &broken_locants,
                    LocantPos *locant_path)
{
  if(broken_locants.empty())
    return true; 

  for (unsigned char loc_broken : broken_locants){
    unsigned char calculate_origin = loc_broken;
    unsigned int pos = 0;
    while( (calculate_origin - 23) > 128){
      calculate_origin += -23;
      pos++;
    }


    // position here decodes where to link them
    unsigned char parent = '\0';
    parent = INT_TO_LOCANT(128 + calculate_origin); // relative positioning
    
                                                    
    if(pos == 2 || pos == 3)
      parent = LOCANT_TO_INT(parent) + 128;
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

  return true; 
}
#endif

/*
 Some rules for solving ring systems in WLN. 

 1) From a given locant, the path walks sequentially forward to the next position 
    C6 = C,D,E,F,G,H , bond C->H

    * if the locant already has saturated bonds, move the next avaliable locant
    A6, when A already has all its bonds, A, J, K, L, M (EOC), B->M. 
 
 2) The lowest locant in that ring should be the one specified.
    If a ring has A,B,C,D,E,F, bonding is A6

 4) The walk is taken by the highest set of locants, If A->F, then A6 goes A, F, G, H, I, J

 5) If a locant contains a branch off the locant path, it is always involved in its parent ring
    and taken last. 
    M has M-:  M6 = M, N, O, P, Q, M-
    sequentally, the broken locants sits between its parent and the next symbol e.g M < M- < N

 6) Each normal position can have maximum two broken children, H-, H-&, of which these can have two etc, 
    given with appending another dash H- can have H-- and H--&. 

 7) Broken locants are only "alive" when their parent locant has been targeteted.
    B has B-, L666 A->F, A->J, A->M, A is saturated therefore B->M, B has B-, backtrack: B- -> L
 
 8) The effect of bridging is that the bonds allowed is lowered by 1.
    Bridging A brings its allowed connections to 2, since it is at the start of the chain
    bridging B brings its allowed connections to 1, since its in the middle of the chain. 
*/
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

    local_size+= - broken_locants.size(); // shouldnt be possible, lets see.
    if(OPT_DEBUG)
      fprintf(stderr,"  calculated size: %c(%d)\n",INT_TO_LOCANT(local_size),local_size);
  }
  else
    local_size = LOCANT_TO_INT(size_designator);
  
  if(OPT_DEBUG)
    fputc('\n',stderr); 

  // create all the nodes in a large straight chain and assign how many bonds
  // each atom is allowed to take
  
  // just use the contructors with new
  LocantPos *locant_path = new LocantPos[local_size];  

  WLNSymbol *curr = 0; 
  WLNSymbol *prev = 0; 

  for (unsigned int i=0;i<local_size;i++){
    unsigned char loc = INT_TO_LOCANT(i+1);
    
    // default ring chain connections
    if(i == 0 || i == local_size-1)
      locant_path[i].allowed_connections = 2; 
    else
      locant_path[i].allowed_connections = 1; 

    if(!ring->locants[loc]){
      curr = AllocateWLNSymbol('C',graph);
      if(!curr)
        return 0;

      curr->allowed_edges = 4;
      assign_locant(loc,curr,ring);
      locant_path[i].locant = curr; 
    }
    else{
      curr = ring->locants[loc];
      locant_path[i].locant = curr; 
      if(curr->ch == 'X')
        locant_path[i].allowed_connections++;
      else if(curr->ch == '*')
        locant_path[i].allowed_connections = 6;

    }

    if(bridge_locants[loc])
      locant_path[i].allowed_connections--; // decrement any bridging positions
      
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

#if DEPRECATED
  // broken locant map spawn + pseudo locants map spawn 
  if(!set_up_broken(ring,graph,broken_locants,broken_lookup,spawned_broken,allowed_connections) || 
     !set_up_pseudo(ring,graph,pseudo_locants,pseudo_lookup))
    return false;
#endif

  // calculate bindings and then traversals round the loops
  unsigned int pseudo_pairs = pseudo_locants.size()/2;
  unsigned char max_locant = INT_TO_LOCANT(local_size);
  
  for (unsigned int i=0;i<ring_assignments.size();i++){
    
    std::pair<unsigned int, unsigned char> component = ring_assignments[i];
    
    bool aromatic           = aromaticity[i];
    unsigned int comp_size  = component.first;
    unsigned int start_char  = component.second;
    ring->aromatic_atoms = ring->aromatic_atoms ? 1:aromatic; 

    if(start_char > max_locant){
      fprintf(stderr,"Error: out of bounds locant access in cyclic builder\n");
      return 0;
    }
    
#if PSEUDO_PAIRS
    // If we need to catch fuse, we do
    if(i==ring_assignments.size()-1 && pseudo_pairs){
      bool caught = false;
      for(unsigned int s = 1; s<=local_size;s++){
        if(pseudo_lookup[INT_TO_LOCANT(s)]){
          unsigned char pbind_2 = INT_TO_LOCANT(s); 
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
#endif

    // --- MULTI ALGORITHM, see notes on function for rules and use --- 
    
    unsigned int  path_size   = 0;  
    unsigned char end_char    = 0; 
    unsigned int  over_shoot  = 0; // simplification on the end of chain logic 

    LocantPos *start_locant   = &locant_path[ LOCANT_TO_INT(start_char-1) ]; 
    LocantPos *curr_locant    = &locant_path[ LOCANT_TO_INT(start_char-1) ]; 
    start_locant->locant->aromatic = start_locant->locant->aromatic==1 ? 1:aromatic;

    // -1 as one locant is already given in start
    while(path_size < comp_size-1){
      WLNEdge*      edge_taken  = 0; 
      unsigned char highest_loc = 0; // highest of each child iteration 

      // walk the ring
      for(unsigned int ei=0;ei<curr_locant->locant->barr_n;ei++){
        WLNSymbol *child = curr_locant->locant->bond_array[ei].child;
        unsigned char child_loc = ring->locants_ch[child];
        if(child_loc > highest_loc){
          highest_loc = child_loc;
          edge_taken = &curr_locant->locant->bond_array[ei]; 
        }
      }

      if(!highest_loc){
        // if at the end of the path, its okay, break the loop and do post shifting 
        if(curr_locant->locant != ring->locants[max_locant]){
          fprintf(stderr,"Error: highest locant not found in path walk\n");
          return false;
        }

        over_shoot++;
        path_size++; 
      }
      else {
        curr_locant = &locant_path[LOCANT_TO_INT(highest_loc-1)];
        curr_locant->locant->aromatic = curr_locant->locant->aromatic ? 1:aromatic;
        edge_taken->aromatic = edge_taken->aromatic==1 ? 1:aromatic;
        end_char = highest_loc; 
        path_size++; 
      }
    }
      
    for(;;){

      if(start_locant->allowed_connections > 0){
#if RARE
        // very rare case where we get the right path a different way to normal
        while((!allowed_connections[bind_2] || bind_2 == bind_1) && bind_2 < INT_TO_LOCANT(local_size))
          ring_path[path_size-1] = ++bind_2;
#endif
        if(OPT_DEBUG)
          fprintf(stderr,"  fusing (%d): %c --> %c\n",comp_size,start_char,end_char);

        if(!AddEdge(curr_locant->locant,start_locant->locant)){
          fprintf(stderr,"Error: failed to bond locant path edge\n");    
          return false;
        }
        
        WLNEdge *new_edge = &start_locant->locant->bond_array[start_locant->locant->barr_n-1];
        new_edge->aromatic = new_edge->aromatic ? 1: aromatic;
        start_locant->allowed_connections--; 
        break;
      }
      else{
        // increase the start char and move the path locant
        start_char++;
        start_locant = &locant_path[LOCANT_TO_INT(start_char-1)]; 
        
        // the current then moves back by 1

        if(over_shoot)
          over_shoot--;
        else
          end_char--; 
        
        curr_locant = &locant_path[LOCANT_TO_INT(end_char-1)]; 
      }
    }


#if PSEUDO
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
#endif

    // sanitize the ring path, if we've got duplicates (i.e hit the end of the path),
    // we decrement and loop back round

    // shifting now performed here should be more stable
#if DEPRECATED
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
        
        while((!allowed_connections[bind_2] || bind_2 == bind_1) && bind_2 < INT_TO_LOCANT(local_size))
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
        if(bind_1 > INT_TO_LOCANT(local_size)+1)
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
#endif
  }
  
  delete [] locant_path;
  return local_size; 
}


unsigned char create_relative_position(unsigned char parent){
  // A = 129
  unsigned int relative = 128 + LOCANT_TO_INT(parent);
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
bool post_saturate( std::vector<LocantPair> &bonds, 
                    unsigned int final_size,
                    WLNRing *ring){

  // post saturate bonds
  for (LocantPair bond_pair : bonds){
    
    unsigned char loc_1 = bond_pair.first;
    unsigned char loc_2 = bond_pair.second;

    if(loc_2 > INT_TO_LOCANT(final_size)){
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
bool FormWLNRing(WLNRing *ring, const char *wln_block,unsigned int i, unsigned int len,WLNGraph &graph,unsigned char spiro_atom='\0'){

  bool warned             = false;  // limit warning messages to console
  bool heterocyclic       = false;  // L|T designator can throw warnings

  unsigned int state_multi          = 0; // 0 - closed, 1 - open multi notation, 2 - expect size denotation
  unsigned int state_pseudo         = 0; 
  unsigned int state_aromatics      = 0;
  
  int pending_charge      = 0; 
  unsigned char inline_unsaturate    = 0; 
  unsigned int  expected_locants      = 0;
  unsigned char ring_size_specifier   = '\0';
  
  bool locant_attached                = false;  // more sensible way of modifying the locants 
  unsigned char positional_locant     = 'A';    // have A as a default, lets tidy around this

  std::vector<bool> aromaticity; 
  std::vector<LocantPair>  saturations;
  
  // make the symbols as we go, or overwrite them
  WLNSymbol *locant_a = 0; 
  WLNSymbol *locant_b = 0;
  WLNEdge   *edge     = 0; 

#if MODERN
  std::vector<LocantPair>  stereo;
#endif

  std::vector<unsigned char>                multicyclic_locants;
  std::vector<unsigned char>                pseudo_locants;
  std::map<unsigned char,unsigned int>      bridge_locants;
  std::set<unsigned char>                   broken_locants;
  
  std::vector<std::pair<unsigned int, unsigned char>>  ring_components;

  // broken locants start at A = 129 for extended ascii 
  // first is the standard 'X-' second is 'X-&', third is 'X--', fourth is 'X--&' and so on

  unsigned char ch = wln_block[i];
  while(i < len){

character_start_ring:
    if(state_multi == 3 && ch != '-' && ch != '&'){
      state_multi = 0; 
      positional_locant = 'A';
    }

    if(inline_unsaturate && positional_locant && (ch != '-' || ch != '&')){
      
      if(!ring->locants[inline_unsaturate]){
        locant_a = AllocateWLNSymbol('C',graph);
        locant_a->allowed_edges = 4; 
        assign_locant(inline_unsaturate, locant_a,ring); 
      }
      else 
        locant_a = ring->locants[inline_unsaturate]; 
      
      if(!ring->locants[positional_locant]){
        locant_b = AllocateWLNSymbol('C',graph);
        locant_b->allowed_edges = 4; 
        assign_locant(positional_locant, locant_b,ring); 
      }
      else 
        locant_b = ring->locants[positional_locant]; 
      
      AddEdge(locant_b, locant_a); // graph flows smaller to larger
      edge = &locant_a->bond_array[locant_a->barr_n-1]; 
      
      if(!unsaturate_edge(edge, 1))
        return false;
      if(wln_block[i+1] == 'U'){
        if(!unsaturate_edge(edge, 1))
          return false;
        i++; 
      }

      inline_unsaturate = 0; 
      positional_locant = 0; 
    }

    switch(ch){
      case ' ':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(expected_locants)
          return Fatal(i,"Error: not enough locants before space character");
        
        else if(state_multi == 1)
          state_multi = 2;
        else if (state_pseudo)
          state_pseudo = 0;
        else if(positional_locant && locant_attached){
          if(ring_components.empty())
            return Fatal(i,"Error: assigning bridge locants without a ring");
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
          ring_size_specifier += AMPERSAND_EXPAND;
        }
        else if (state_pseudo){
          pseudo_locants.back() += AMPERSAND_EXPAND;
        }
        else if(positional_locant && locant_attached){
          positional_locant += AMPERSAND_EXPAND;
        }
        else{
          // if unhandled, then it must be an aromatic start
          state_aromatics = 1;
          aromaticity.push_back(1);
        }
        break;

#if MODERN
      case '<':{
          str_buffer.clear();
          bool found_close = false; 
          unsigned int k = i+1; // use the block string and iterate
          unsigned int stop = 0; 

          std::string dash_buffer; // since this could now be multiple things
    
          while(k < block.size()){
            if(block[k] == '>'){
              // this calculates the gap to the next '-'
              found_close = true; // if there was a double '--' this will have gap zero
              break;
            }

            dash_buffer.push_back(block[k]);
            k++;
            stop++; 
          }

          if(!found_close)
            return Fatal(i+start, "Error: unbalanced < > characters"); 
          else if (dash_buffer.empty())
            return Fatal(i+start, "Error: empty < > characters"); 

          // split the dash buffer up into large ring or not large ring
          //
          unsigned int mode = 0;  // 1- ring_size, 2 = element of some sort 
          int charge = 0;
          bool hypervalent = false;
          std::string ch_store; 
      
          for(unsigned int c=0;c<dash_buffer.size();c++){
            unsigned char ch = dash_buffer[c];    
            if(ch == '>')
              break;
            else{
              if(ch <= 'Z' && ch >= 'A'){
                if(!mode)
                  mode = 2;
                if(mode == 2)
                  str_buffer.push_back(ch);
              }
              else if(ch >= '0' && ch <= '9'){
                ch_store += ch; 
                if(!mode)
                  mode = 1; 
                if(mode == 1) 
                  str_buffer.push_back(ch); // for large rings
              }
              else if(ch == '-'){
                if(ch_store.empty() && str_buffer.empty())
                  hypervalent = true;
                else if (ch_store.empty())
                  charge = -1; 
                else
                 charge = -isNumber(ch_store); 
              }
              else if(ch == '+'){
                if(ch_store.empty())
                  charge = 1;
                else
                  charge = isNumber(ch_store); 
              }
              else 
                return Fatal(i, "Error: invalid character in special '< >'");
              
            }
          }
          

          if(mode == 1){
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
          else if (mode == 2){
            if(str_buffer.size() == 1){
              
              if(hypervalent){
                if(positional_locant != spiro_atom){
                  WLNSymbol* new_locant = assign_locant(positional_locant,define_hypervalent_element(str_buffer[0],graph),ring);  // elemental definition 
                  if(!new_locant)
                    return Fatal(i+start, "Error: could not create hypervalent element");
                  
                  new_locant->str_position = start+i+1+1;
                  new_locant->charge = charge; 
                  ring->position_offset[new_locant] = i+1;
                  if(OPT_DEBUG)
                    fprintf(stderr,"  assigning hypervalent %c to position %c\n",str_buffer[0],positional_locant);
                }
                else 
                  positional_locant++;

                positional_locant++;
                locant_attached = false;
              }
              else {
                ch = str_buffer[0];
                pending_charge = charge; 
                i+=stop;
                goto character_start_ring;
              }

            }
            else if (str_buffer.size() == 2){
      
              if(positional_locant != spiro_atom){

                // must be a special element 
                WLNSymbol* new_locant = assign_locant(positional_locant,define_element(str_buffer,graph),ring);  // elemental definition
                if(!new_locant)
                  return Fatal(i+start, "Error: could not create periodic code element");
                
                new_locant->charge = charge; 
                new_locant->str_position = start+i+1+1;
                ring->position_offset[new_locant] = i+1;
                if(OPT_DEBUG)
                  fprintf(stderr,"  assigning element %s to position %c\n",str_buffer.c_str(),positional_locant);
                
                positional_locant++;
              }
              else
                positional_locant++;
              
              locant_attached = false;
            }
          }

          i+=stop;
          break;
        }
        break; 
#endif

      case '/':
        if(state_aromatics)
          return Fatal(i,"Error: invalid character in the aromaticity assignment block");
             
        expected_locants = 2; 
        state_pseudo = 1;
        break; 
      

      // turn this into a look ahead type bias in order to significantly tidy this up
      
/*
      Some logic here; 
      1)  If size is followed by a '-', it immediately closes size and starts the aromatic state
          e.g B&-TJ

      2)  If dash is started on a locant, and either is not closed or closed with no gap before 
          the 2 character limit or a space is seen, it must be a branching locant

      3)  If a branching locant is deemed, the tree logic is that - is left, and & is right, with a limit of 
          2. A buffer to consume and move all the positions is needed, and can track the tree size
*/
      case '-':
        if(state_multi){
          if(state_multi == 1 && !multicyclic_locants.empty())
            multicyclic_locants.back() = positional_locant;
          else if(state_multi == 3){
            state_multi = 0;
            state_aromatics = true;
          }
        }
        else if(!expected_locants){
          
          // conditions for - within a ring system:
          // i+2 == '-' for hypervalent
          // i+3 == '-' for element and large ring with 2 digits
          // 
          // branching locants - come to these another time

          // two things to do in the lookahead here, determine whether branching or non-branching, 
          // and then access the branching tree
          
          if( i + 3 < len && wln_block[i+3] == '-'){
             
            if(wln_block[i+1] >= '0' && wln_block[i+1] <= '9' && wln_block[i+2] >= '0' && wln_block[i+2] <= '9'){
              unsigned int big_ring = ( (wln_block[i+1]-'0') * 10) + wln_block[i+2]-'0';
              if(!big_ring)
                return Fatal(i,"Error: non numeric value entered as ring size");

              ring_components.push_back({big_ring,positional_locant}); //big ring
              positional_locant = 'A';
              locant_attached = false;
              i+= 3; 
            }
            else if (wln_block[i+1] >= 'A' && wln_block[i+1] <= 'Z' && wln_block[i+2] >= 'A' && wln_block[i+2] <= 'Z'){
              
              if(positional_locant != spiro_atom){

                // must be a special element
                WLNSymbol *new_locant =  define_element(wln_block[i+1],wln_block[i+2],graph); 
                if(!new_locant)
                  return Fatal(i, "Error: could not create periodic code element");

                assign_locant(positional_locant,new_locant,ring);  // elemental definition

                new_locant->str_position = i + 2; // attaches directly to the starting letter
                ring->position_offset[new_locant] = i+1;
      
                if(OPT_DEBUG)
                  fprintf(stderr,"  assigning element %c%c to position %c\n",wln_block[i+1],wln_block[i+2],positional_locant);
                
                positional_locant++;
              }
              else
                positional_locant++;
              
              locant_attached = false;
              i+= 3; 
            }
          }
          else if (i+2 < len && wln_block[i+2] == '-'){
            if(positional_locant != spiro_atom){
              WLNSymbol *new_locant = define_hypervalent_element(wln_block[i+1],graph); 
              if(!new_locant)
                return Fatal(i, "Error: could not create hypervalent element");
              
              assign_locant(positional_locant,new_locant,ring);
              
              new_locant->str_position = i+1+1;
              ring->position_offset[new_locant] = i+1;
              if(OPT_DEBUG)
                fprintf(stderr,"  assigning hypervalent %c to position %c\n",wln_block[i+1],positional_locant);
            }
            else 
              positional_locant++;

            positional_locant++;
            locant_attached = false;
            i+= 2; 
          } 
        }
        break;
      

      case '0':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(!ring_components.empty()){
          if(!positional_locant)
            positional_locant = 'A';

          if(OPT_DEBUG)
            fprintf(stderr,"  placing pi bond charge on locant - %c\n",positional_locant);

          
          WLNSymbol *zero_carbon = AllocateWLNSymbol('C',graph);
          zero_carbon->allowed_edges = 3; 
          assign_locant(positional_locant++,zero_carbon,ring);
          zero_carbon->str_position = i+1;
          ring->position_offset[zero_carbon] = i; 
          zero_carbon->charge--; 
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
          return Fatal(i,"Error: invalid character in the aromaticity assignment block");
        
        if (i > 1 && wln_block[i-1] == ' '){
#if MODERN
          state_multi = 2;
#else
          state_multi   = 1; // enter multi state
          expected_locants = ch - '0';
#endif
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
          return Fatal(i,"Error: invalid character in the aromaticity assignment block");
        

        if(expected_locants){
          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            return Fatal(i,"Error: unhandled locant rule");
          

          positional_locant = ch; // use for look back
          locant_attached = true;
          expected_locants--;
          break;
        }
        else if(state_multi == 2){
          ring_size_specifier = ch;
          state_multi = 3;
        }
        else if(i>0 && wln_block[i-1] == ' '){
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
          return Fatal(i,"Error: invalid character in the aromaticity assignment block");
        
        if(expected_locants){
          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            Fatal(i,"Error: unhandled locant rule");

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

          switch(ch){
            
            case 'S':
            case 'P':
              if(!heterocyclic)
                warned = true;
              
              if(!ring->locants[positional_locant]){
                locant_a = AllocateWLNSymbol(ch,graph);
                assign_locant(positional_locant,locant_a,ring);
              }
              else 
                locant_a = ring->locants[positional_locant]; 

              locant_a->str_position = (i+1); 
              locant_a->charge = pending_charge; 
              ring->position_offset[locant_a] = i; 
              
              if(ch == 'P')
                locant_a->allowed_edges = 5;
              else
                locant_a->allowed_edges = 6;
              
              pending_charge = 0; 
              positional_locant++; 
              break;

            case 'Y':
            case 'X':
            case 'K':
              if(!heterocyclic && ch=='K')
                warned = true;
            
              if(!ring->locants[positional_locant]){
                locant_a = AllocateWLNSymbol(ch,graph);
                assign_locant(positional_locant,locant_a,ring);
              }
              else 
                locant_a = ring->locants[positional_locant]; 

              locant_a->str_position = (i+1); 
              locant_a->charge = pending_charge; 
              ring->position_offset[locant_a] = i; 

              locant_a->allowed_edges = 4;
              
              pending_charge = 0; 
              positional_locant++; 
              break;

            case 'Z': // treat as NH2
              if(!heterocyclic)
                warned = true;

              if(!ring->locants[positional_locant]){
                locant_a = AllocateWLNSymbol(ch,graph);
                assign_locant(positional_locant,locant_a,ring);
              }
              else 
                locant_a = ring->locants[positional_locant]; 

              locant_a->str_position = (i+1); 
              locant_a->charge = pending_charge; 
              ring->position_offset[locant_a] = i; 

              locant_a->allowed_edges = 3;
              locant_a->explicit_H = 2; 
              
              pending_charge = 0; 
              positional_locant++; 
              break;

            case 'N':
            case 'B':
              if(!heterocyclic)
                warned = true;


              if(!ring->locants[positional_locant]){
                locant_a = AllocateWLNSymbol(ch,graph);
                assign_locant(positional_locant,locant_a,ring);
              }
              else 
                locant_a = ring->locants[positional_locant]; 

              locant_a->str_position = (i+1); 
              locant_a->charge = pending_charge; 
              ring->position_offset[locant_a] = i; 

              locant_a->allowed_edges = 3;
              
              pending_charge = 0; 
              positional_locant++; 
              break;

            case 'M':
              if(!heterocyclic)
                warned = true;

              if(!ring->locants[positional_locant]){
                locant_a = AllocateWLNSymbol(ch,graph);
                assign_locant(positional_locant,locant_a,ring);
              }
              else 
                locant_a = ring->locants[positional_locant]; 

              locant_a->str_position = (i+1); 
              locant_a->charge = pending_charge; 
              ring->position_offset[locant_a] = i; 

              locant_a->allowed_edges = 2;
              locant_a->explicit_H = 1;
              
              pending_charge = 0; 
              positional_locant++; 
              break;

            case 'O':
            case 'V':
              if(!heterocyclic && ch == 'O')
                warned = true;

              if(!ring->locants[positional_locant]){
                locant_a = AllocateWLNSymbol(ch,graph);
                assign_locant(positional_locant,locant_a,ring);
              }
              else 
                locant_a = ring->locants[positional_locant]; 

              locant_a->str_position = (i+1); 
              locant_a->charge = pending_charge; 
              ring->position_offset[locant_a] = i; 

              locant_a->allowed_edges = 2;
              
              pending_charge = 0; 
              positional_locant++; 
              break;

        
            case 'U':
              if(i+1 < len && wln_block[i+1] != '-'){

                if(!ring->locants[positional_locant]){
                  locant_a = AllocateWLNSymbol('C',graph);
                  locant_a->allowed_edges = 4; 
                  assign_locant(positional_locant, locant_a,ring); 
                }
                else 
                  locant_a = ring->locants[positional_locant]; 
                
                if(!ring->locants[positional_locant+1]){
                  locant_b = AllocateWLNSymbol('C',graph);
                  locant_b->allowed_edges = 4; 
                  assign_locant(positional_locant+1, locant_b,ring); 
                }
                else 
                  locant_b = ring->locants[positional_locant+1]; 
                
                AddEdge(locant_b, locant_a); // graph flows smaller to larger
                edge = &locant_a->bond_array[locant_a->barr_n-1]; 
                
                if(!unsaturate_edge(edge, 1))
                  return false;
                if(wln_block[i+1] == 'U'){
                  if(!unsaturate_edge(edge, 1))
                    return false;
                  i++; 
                }
                positional_locant++;
              }
              else if (i+1 < len && wln_block[i+1] == '-'){
                inline_unsaturate = positional_locant; 
                positional_locant = 0; 
                i++; // skip the dash
              }
              break;
            
#if MODERN
            case 'A':
            case 'D':
              // no need to put this in implied, it has to be specified
              if(i < len - 3 && block[i+1] == '-' && block[i+2] == ' '){
                
                // can double bond to a amped locant
                unsigned int k = 1;
                unsigned char dloc = block[i+3]; 
                while(block[k+i+3] == '&'){
                  dloc+=AMPERSAND_EXPAND;
                  k++;
                } 

                LocantPair loc_pair(positional_locant,dloc);
                if(ch == 'A')
                  loc_pair.stereo = 2;
                else
                  loc_pair.stereo = 1;

                stereo.push_back({positional_locant,dloc});
                block_str += 2+k;
                i += 2+k;
              }
              else{  
                LocantPair loc_pair(positional_locant,positional_locant+1); 
                if(ch == 'A')
                  loc_pair.stereo = 2;
                else
                  loc_pair.stereo = 1;
                stereo.push_back(loc_pair);
              }
              positional_locant++;
              break;
            
#endif
            // externally bonded to the symbol as a locant symbol
            case 'W':
              if(!heterocyclic)
                warned = true;

              if(positional_locant > 'A')
                positional_locant--;

              // the carbon might not be created yet
              if(!ring->locants[positional_locant]){
                locant_a = AllocateWLNSymbol('C',graph);
                assign_locant(positional_locant,locant_a,ring);
                locant_a->allowed_edges = 2;
                locant_a->str_position = (i+1); 
                ring->position_offset[locant_a] = i; 
              }
              else
                locant_a = ring->locants[positional_locant];
              
              if(locant_a->ch == 'N' && locant_a->allowed_edges == 3)
                locant_a->allowed_edges++;

              locant_b = AllocateWLNSymbol('W',graph);
              locant_b->allowed_edges = 3;
              locant_b->inRing = ring;
              ring->position_offset[locant_b] = i; 
              if(!AddEdge(locant_b, locant_a))
                return Fatal(i,"Error: failed to create bond");

              edge = &locant_a->bond_array[locant_a->barr_n-1];  
              if(!unsaturate_edge(edge,2))
                return Fatal(i,"Error: failed to unsaturate edge");
              
              positional_locant++; 
              break;

            // has the effect of unaromatising a bond, remove from edge consideration
            case 'H':{
              LocantPair loc_pair(positional_locant,positional_locant+1);  
              saturations.push_back(loc_pair);
              break;
            }

            default:
              return Fatal(i,"Error: invalid character in atom assignment within ring notation");
          }

          locant_attached = false; // locant is no longer primary
        }
        else if(i>0 && wln_block[i-1] == ' '){
          if(ring_size_specifier && ch > ring_size_specifier)
            return Fatal(i, "Error: specifying locants outside of allowed range"); 

          positional_locant = ch;
          locant_attached = true;
        }
        break;

      case 'L':
        if(positional_locant >= 128)
          broken_locants.insert(positional_locant);

        if(state_aromatics)
          return Fatal(i,"Error: invalid character in the aromaticity assignment block");
        

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
            return Fatal(i,"Error: unhandled locant rule");
          

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
          if(i>0 && wln_block[i-1] == ' '){
            if(ring_size_specifier && ch > ring_size_specifier)
              return Fatal(i, "Error: specifying locants outside of allowed range"); 
            
            positional_locant = ch;
            locant_attached = true;
          }
          else
            return Fatal(i,"Error: symbol is in an unhandled state, please raise issue if this notation is 100%% correct");
          
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
            return Fatal(i,"Error: unhandled locant rule");
          

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
            return Fatal(i,"Error: assigning bridge locants without a ring");
          else
            bridge_locants[positional_locant] = true;

          state_aromatics = 1;
          aromaticity.push_back(0);
        }
        else{
          if(i>0 && wln_block[i-1] == ' ' && wln_block[i+1] != 'J'){
            if(ring_size_specifier && ch > ring_size_specifier)
              return Fatal(i, "Error: specifying locants outside of allowed range"); 
            
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
        
        if (i == len-1){
          
          if(ring_components.empty())
            return Fatal(i,"Error: error in reading ring components, check numerals in ring notation");
          

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
            return Fatal(i,"Error: mismatch between number of rings and aromatic assignments");
          
          break;
        }
        if(expected_locants){

          if(state_multi)
            multicyclic_locants.push_back(ch);
          else if (state_pseudo)
            pseudo_locants.push_back(ch);
          else
            return Fatal(i,"Error: unhandled locant rule");
        
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
            return Fatal(i,"Error: assigning bridge locants without a ring");
          
          else
            bridge_locants[positional_locant] = true;
        }
        else{
          if(i>0 && wln_block[i-1] == ' '){
            positional_locant = ch;
            locant_attached = true;
          }
          else
            return Fatal(i,"Error: symbol is in an unhandled state, please raise issue if this notation is 100%% correct");
        }
        break;

      default:
        break;
    }
    
    ch = wln_block[++i];
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

    fprintf(stderr,"  multi size: %c(%d)\n",ring_size_specifier ,ring_size_specifier ? LOCANT_TO_INT(ring_size_specifier) : 0);
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
    return Fatal(i, "Error: failed to build WLN cycle unit");
  
  if(! post_saturate(saturations,final_size,ring))
    return Fatal(i, "Error: failed on post ring bond (un)/saturation");
  
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

    switch (sym->ch){
      case 'X':
      case 'Y':
      case '#':
      case 'C':
      case '1': // same for ring carbons
        break;

      default:
        for(unsigned int ei=0;ei<sym->barr_n;ei++){
          edge = &sym->bond_array[ei]; 
          if( (edge->child->ch == 'O' ||
              edge->child->ch ==  'P'  || 
              edge->child->ch ==  'N'  || 
              edge->child->ch ==  'S') &&
              edge->child->num_edges == 1 && edge->child->charge == 0)
          {
            while((sym->num_edges < sym->allowed_edges && edge->order < 3) && 
              (edge->child->num_edges < edge->child->allowed_edges)){
              if(!unsaturate_edge(edge,1))
                return false;
            }
          }
        }

        for(unsigned int ei=0;ei<sym->parr_n;ei++){
          edge = &sym->prev_array[ei]; 
          if( (edge->child->ch == 'O' ||
              edge->child->ch ==  'P'  || 
              edge->child->ch ==  'N'  || 
              edge->child->ch ==  'S') &&
              edge->child->num_edges == 1 && edge->child->charge == 0)
          {
            while((sym->num_edges < sym->allowed_edges && edge->order < 3) && 
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

/* provides methods for `kekulising` wln ring structures, using blossums to maximise pairs */
bool WLNKekulize(WLNGraph &graph){
  for(unsigned int i=0;i<graph.ring_count;i++){
    WLNRing *wring = graph.RINGS[i]; 
    if(wring->aromatic_atoms){

      int   *MatchR = (int*)malloc(sizeof(int) * wring->rsize);
      if(!wring->FillAdjMatrix() || !MatchR)
        return false;
      

      for (unsigned int i=0;i<wring->rsize;i++)
        MatchR[i] = -1;
  

      if(IsBipartite(wring) && !WLNRingBPMaxMatching(wring,MatchR)){
        if(MatchR)
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
          unsigned char floc = INT_TO_LOCANT(i+1);
          unsigned char sloc = INT_TO_LOCANT(MatchR[i]+1);
          
          WLNSymbol *f = wring->locants[floc];
          WLNSymbol *s = wring->locants[sloc];
          if(f && s){
            WLNEdge *edge = search_edge(f,s);
            if(edge && edge->order == 1 && !unsaturate_edge(edge,1))
              return false;
          }
        }
        MatchR[MatchR[i]] = 0; // remove from matching
      }
      
      if(MatchR)
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
  if (OPT_DEBUG)
    fprintf(stderr, "Parsing WLN notation: %s\n",wln_ptr);

  ObjectStack branch_stack;   // access to both rings and symbols
  //branch_stack.reserve(512);  // reasonable size given
  
  WLNSymbol *curr       = 0;
  WLNSymbol *prev       = 0;
  WLNSymbol *last       = 0; // last handled, ignores terminators
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
  bool pending_numbers          = false;
  bool pending_locant_skips     = false; // special case handling 
  
  int pending_charge = 0; 
  unsigned int pending_stereo = false; // 0 = none, 1 = dashed, 2 = wedged 
  unsigned int inline_unsaturate = 0;    // 'U' style bonding
    
  // hold three sequential digits, error if 4, uses the ascii values '1', '2' etc 
  unsigned char d1 = 0; 
  unsigned char d2 = 0; 
  unsigned char d3 = 0; 
  
  unsigned char on_locant = '\0';         // locant tracking
  unsigned int locant_skips = false;                   // handle skipping of 'J' if in cyclic notation legitimately 
  
  // allows consumption of notation after block parses
  unsigned int block_start = 0;
  unsigned int block_end = 0;

  unsigned int i=0;
  unsigned int len = strlen(wln_ptr);
  unsigned char ch = wln_ptr[i];
  
  while(ch)
  {  
character_start:
    // this will need to resolved at the end as well
    if(pending_numbers && (ch < '0' || ch > '9') && ch != '/' ){
      
      // bases are given by the number reads
      unsigned int carbon_len = 0; 
      if(d3){
        carbon_len += (100*(d1-'0')); 
        carbon_len += (10*(d2-'0')); 
        carbon_len += (1*(d3-'0')); 
      }
      else if(d2){
        carbon_len += (10*(d1-'0')); 
        carbon_len += (1*(d2-'0')); 
      }      
      else
        carbon_len += (1*(d1-'0')); 
      
      d1=0;
      d2=0;
      d3=0;

      // create the head 
      WLNSymbol *curr = AllocateWLNSymbol('#', graph);
      curr->str_position = i;
      curr->special = std::to_string(carbon_len); 
      
      if(carbon_len > 1)
        curr->allowed_edges = 6; // can hold a triple bond on either side
      else
        curr->allowed_edges = 4;

      if(prev){
        
        if(!AddEdge(curr, prev))
          return Fatal(i, "Error: failed to bond to previous symbol");

        edge = &prev->bond_array[prev->barr_n-1]; 
        edge->stereo = pending_stereo; 
        pending_stereo = 0; 

        if(inline_unsaturate){
          if(!unsaturate_edge(edge,inline_unsaturate))
            return Fatal(i, "Error: failed to unsaturate bond"); 
          inline_unsaturate = 0;
        }
      }


      pending_numbers = false;
      prev = curr;
      last = curr; 
      last = prev; 
      cleared = false;
    }
    
    // slight logic change to J variation
    if(pending_locant_skips && (ch < '0' || ch > '9') ){
#if MODERN
      locant_skips = 1;
#else
      if(d3){
        locant_skips += (100*(d1-'0')); 
        locant_skips += (10*(d2-'0')); 
        locant_skips += (1*(d3-'0')); 
      }
      else if(d2){
        locant_skips += (10*(d1-'0')); 
        locant_skips += (1*(d2-'0')); 
      }      
      else
        locant_skips += (1*(d1-'0')); 
      
      d1=0;
      d2=0;
      d3=0;
#endif
      pending_locant_skips = false; 
    }

    switch (ch)
    {

    case '0': // cannot be lone, must be an addition to another num
      if(pending_J_closure){
        if(pending_locant_skips){
          if(d1 && d2 && d3)
            return Fatal(i, "Error: specifying a number greater than 3 digits - WLN isn't meant for this!");
          else if (d1 && d2)
            d3 = ch; 
          else if(d1)
            d2 = ch; 
          else 
            d1 = ch; 
        }
        break;
      }

      else if (pending_locant){
        
        if(pending_inline_ring && prev && !prev->inRing)
          prev->charge++;

        prev = 0;
        on_locant = '0';
        pending_locant = false;
      }
      else if (pending_numbers || pending_negative_charge || cleared){
        if(d1 && d2 && d3)
          return Fatal(i, "Error: specifying a number greater than 3 digits - WLN isn't meant for this!");
        else if (d1 && d2)
          d3 = ch; 
        else if(d1)
          d2 = ch; 
        else 
          d1 = ch; 

        if(cleared)
          pending_numbers = true; 
        break;
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
        if(i > 0 && wln_ptr[i-1] == ' '){
          pending_locant_skips = true; 
        }
        
        if(pending_locant_skips){
          if(d1 && d2 && d3)
            return Fatal(i, "Error: specifying a number greater than 3 digits - WLN isn't meant for this!");
          else if (d1 && d2)
            d3 = ch; 
          else if(d1)
            d2 = ch; 
          else 
            d1 = ch; 
        }
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
    
            if(!AddEdge(curr, prev))
              return Fatal(i, "Error: failed to bond to previous symbol");
            
            edge = &prev->bond_array[prev->barr_n-1];
            wrap_ring->macro_return = edge; 
            
            edge->stereo = pending_stereo; 
            pending_stereo = 0; 

            if(inline_unsaturate){
              if(!unsaturate_edge(edge,inline_unsaturate))
                return Fatal(i, "Error: failed to unsaturate bond"); 
              inline_unsaturate = 0;
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
      else {
        if(d1 && d2 && d3)
          return Fatal(i, "Error: specifying a number greater than 3 digits - WLN isn't meant for this!");
        else if (d1 && d2)
          d3 = ch; 
        else if(d1)
          d2 = ch; 
        else 
          d1 = ch; 
        
        
        on_locant = '\0';
        pending_numbers = true; 
        break;
      }
      break;
    

    case 'Y':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
        return Fatal(i,"Error: 'Y' cannot be a locant assignment, please expand [A-W] with &\n");
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 4; // change methyl addition
 
        curr->charge = pending_charge; 
        pending_charge = 0;
        
        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 

          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }
        
        branch_stack.push({0,curr});
        inline_unsaturate = 0;
        prev = curr;
        last = prev; 
      }
      cleared = false;
      break;

    case 'X':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant){
        return Fatal(i, "Error: Wiswesser Uncertainities lead to runaway outcomings");
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        curr->allowed_edges = 4;

        curr->charge = pending_charge; 
        pending_charge = 0;
        
        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }
        
        branch_stack.push({0,curr});
        prev = curr;
        last = prev; 
      }
      cleared = false;
      break;

      // oxygens

    case 'O':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          ring->loc_count++;  
          prev = curr;
          last = prev; 
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

        curr->charge = pending_charge; 
        pending_charge = 0;
        
        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }

        prev = curr;
        last = prev; 
      }
      cleared = false;
      break;

    case 'Q':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = prev; 
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
        curr->explicit_H = 1; 

        curr->charge = pending_charge; 
        pending_charge = 0;
        
        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }

        inline_unsaturate = 0;
        last = curr; 
        prev = return_object_symbol(branch_stack);

        if(!prev)
          prev = curr;
      }
      cleared = false;
      break;

    case 'V':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          } 
        }

        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;

    case 'W':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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
          
          if(prev->ch == 'N' && prev->allowed_edges == 3)
            prev->allowed_edges++;

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(!unsaturate_edge(edge, 2))
            return Fatal(i, "Error: failed to attach W symbol"); 

          if(inline_unsaturate)
            return Fatal(i,"Error: a bond unsaturation followed by dioxo is undefined notation");
        }
        else
          inline_unsaturate = 2;

        // if there is no prev, this is then a character we go from, 
        // else do not update.
        
        last = curr; 
        if(!prev)
          prev = curr;
        else
          prev = return_object_symbol(branch_stack);

      }
      cleared = false;
      break;

      // nitrogens

    case 'N':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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

        curr->charge = pending_charge; 
        pending_charge = 0;
        
        if(prev){
          if(prev->ch == 'W' && curr->allowed_edges == 3)
            curr->allowed_edges++;

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }

        branch_stack.push({0,curr});
        inline_unsaturate = 0;
        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;

    case 'M':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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
        curr->explicit_H = 1; 
        
        curr->charge = pending_charge; 
        pending_charge = 0;

        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){ 
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }

        inline_unsaturate = 0;
        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;

    case 'K':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
        
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
          
        }
        
        branch_stack.push({0,curr});
        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;

    case 'Z':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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
        curr->explicit_H = 2; // this should be a better way

        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
          
        }

        inline_unsaturate = 0;
        last = curr; 
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
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 

          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }

        inline_unsaturate = 0;
        last = curr; 
        prev = return_object_symbol(branch_stack);
        if(!prev)
          prev = curr;
      }
      cleared = false;
      break;

      // inorganics

    case 'B':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
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

        curr->charge = pending_charge; 
        pending_charge = 0;

        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }         
        
        branch_stack.push({0,curr});
        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;

    case 'P':
    case 'S':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        curr = AllocateWLNSymbol(ch,graph);
        curr->str_position = i+1;
        
        curr->charge = pending_charge; 
        pending_charge = 0;
        
        if(ch == 'P')
          curr->allowed_edges = 5;
        else
          curr->allowed_edges = 6;
          
        if(prev){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }
        branch_stack.push({0,curr});
        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;

    // multiply bonded carbon, therefore must be at least a double bond
    case 'C':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
#if MODERN
        // in modern notation, mutliplier carbons are removed
        on_locant = '\0';
        WLNSymbol *curr = AllocateWLNSymbol('#', graph);
        curr->str_position = i;
        curr->special = '1'; 
        
        curr->allowed_edges = 4;
        curr->charge = pending_charge;
        pending_charge = 0; 


        if(prev){
          
          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 

          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }
        
        prev = curr;
        last = curr; 
#else
        on_locant = '\0';
        // set lower case for multiplier carbon
        curr = AllocateWLNSymbol('c',graph);
        curr->str_position = i+1;
        curr->allowed_edges = 4;

        if(prev && i < len - 1){

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }

        prev = curr;
        last = curr; 
#endif
      }
      cleared = false;
      break;

    case 'A':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {

        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");

          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }

        pending_locant = false;
        on_locant = ch;
      }
      else{
#if MODERN
        pending_stereo = 2; 
#else
        return Fatal(i,"Error: locant only symbol used in atomic definition");
#endif
      }
      cleared = false;
      break;
        
    // this can start a chelating ring compound, so has the same block as 'L\T'
    case 'D':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {

#if MODERN
        pending_stereo = 1;     
#else
        if(i < len - 2 && wln_ptr[i+1] == '-' && (wln_ptr[i+2] == 'T' || wln_ptr[i+2] == 'L')){
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
#endif 
      }
      cleared = false;
      break;
        
    // hydrogens explicit
    case 'H':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }
        pending_locant = false;
        on_locant = ch;
      }
      else{
        // treat like the W or U symbol, modify the last symbol if possible, else create H explicit
        on_locant = '\0';
        if(!prev){
          curr = AllocateWLNSymbol(ch,graph);
          
          curr->charge = pending_charge;
          pending_charge = 0; 

          curr->str_position = i+1;
          curr->allowed_edges = 1;
          prev = curr; 
          last = curr; 
        }
        else if (prev && prev->ch == 'c'){
          curr = AllocateWLNSymbol(ch,graph);
          curr->str_position = i+1;
          curr->allowed_edges = 1;
          if(!AddEdge(curr, prev)) 
            return Fatal(i, "Error: failed to bond to previous symbol");
          prev = curr; 
          last = curr; 
        }
        else if (last)
          last->explicit_H++; 
        
        break;
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
      if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          last = curr; 
          prev = curr;
        }
        pending_locant = false;
        on_locant = ch;
      }
      else if (pending_J_closure)
      {
        
        if(locant_skips)
          locant_skips--;
        else if(i>0 && wln_ptr[i-1] != ' '){
          block_end = i;
          
          ring = AllocateWLNRing(graph);
          if(pending_spiro){
            
            if(!prev)
              Fatal(i,"Error: sprio notation opened without a previous atom");
            else{
              ring->locants[on_locant] = prev;
              ring->locants_ch[prev] = on_locant; 
              prev->spiro = true; 
            }
            // check for an aromaticity bond move?
            if(prev && (prev->allowed_edges - prev->num_edges) < 2){

              // spiro would not be possible here, check if a double bond can be shifted
              WLNSymbol *shift = 0;
              WLNEdge *e = 0; 
              for (unsigned int ei=0;ei<prev->barr_n;ei++){
                e = &prev->bond_array[ei];
                if (e->order == 2){
                  if(!saturate_edge(e,1))
                    return Fatal(i, "Error: could not shift aromaticity for spiro ring addition");

                  shift = e->child;
                  break;
                }
              }

              if(!branch_stack.ring)
                return Fatal(i, "Error: ring stack is empty, nothing to fetch");

              unsigned char next_loc = branch_stack.ring->locants_ch[shift]+1;
              if(!next_loc)
                next_loc = 'A'; // must of done the full loop

              e = search_edge(branch_stack.ring->locants[next_loc],shift);
              if(!e && !unsaturate_edge(e, 1))
                return Fatal(i, "Error: failed to re-aromatise previous ring");
            }
            
            // +1 gets the J included
            if(!FormWLNRing(ring,wln_ptr,block_start,i+1,graph,on_locant))
              return false;
          }
          else{
            if(!FormWLNRing(ring,wln_ptr,block_start,i+1,graph))
              return false;
          }
          

          if(pending_ring_in_ring && !wrap_ring){
            wrap_ring = ring; // instant back access
          }

          branch_stack.push({ring,0});
          block_start = 0;
          block_end = 0;

          // does the incoming locant check
          if(pending_spiro)
            pending_spiro = false;
          else if (prev && on_locant && on_locant != '0')
          {
            if (ring->locants[on_locant]){
                      
              if(!AddEdge(ring->locants[on_locant], prev))
                return Fatal(i, "Error: failed to bond to previous symbol");

              edge = &prev->bond_array[prev->barr_n-1]; 
              edge->stereo = pending_stereo; 
              pending_stereo = 0; 

              if(inline_unsaturate){
                if(!unsaturate_edge(edge,inline_unsaturate))
                  return Fatal(i, "Error: failed to unsaturate bond"); 
                inline_unsaturate = 0;
              }
              ring->loc_count++; //in-line locant assumed
            }   
            else
              return Fatal(i,"Error: attaching inline ring with out of bounds locant assignment");
            
          }
          

          on_locant = '\0';
          pending_J_closure = false;
        }

      }
      
      cleared = false;
      break;

    case 'L':
    case 'T':
      if (pending_J_closure || pending_rir_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){
          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        if(i < len - 2 && wln_ptr[i+1] == '-' && (wln_ptr[i+2] == 'T' || wln_ptr[i+2] == 'L')){
          pending_ring_in_ring = true;

          // if pending inline ring doesnt get a L or T, we know it
          // the ring wrap symbol. 

          i++;
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
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }
        pending_locant = false;
        on_locant = ch;
      }
      else
      {
        on_locant = '\0';
        ring = AllocateWLNRing(graph);

        const char* benzyl = "L6J";
        FormWLNRing(ring,benzyl,0,3,graph);
        branch_stack.push({ring,0});

        curr = ring->locants['A'];
        if(prev){
          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 

          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }

        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;

      // bonding

    case 'U':
      if (pending_J_closure){
        if(locant_skips)
          locant_skips--; 
        break;
      }
      else if (pending_locant)
      {
        if(!pending_inline_ring){

          ring = branch_stack.ring;
          curr = ring->locants[ch];
          if(!curr)
            return Fatal(i,"Error: accessing locants out of range");
          
          ring->loc_count++;  
          prev = curr;
          last = curr; 
        }
        pending_locant = false;
        on_locant = ch;
      }
      else if(cleared)
        return Fatal(i, "Error: floating double bond after ionic clear");
      
      else{
        on_locant = '\0';
        inline_unsaturate++;
      }
      break;

      // specials

    case ' ':
      if (pending_J_closure){
        locant_skips = 0; 
        break;
      }
      else if(pending_negative_charge){
        
        int negative_index = 0; 
        if(d3){
          negative_index += (100*(d1-'0')); 
          negative_index += (10*(d2-'0')); 
          negative_index += (1*(d3-'0')); 
        }
        else if(d2){
          negative_index += (10*(d1-'0')); 
          negative_index += (1*(d2-'0')); 
        }      
        else
          negative_index += (1*(d1-'0')); 
      
        d1=0;
        d2=0;
        d3=0;

        if(negative_index < 0)
          return Fatal(i, "Error: assigning non-numerical value to charge index");
        else if(negative_index != 0){
          // find the symbol and increment its charge + 1
          bool found = false;
          for(unsigned int cs = 0;cs<graph.symbol_count;cs++){
            if(graph.SYMBOLS[cs]->str_position == (unsigned int)negative_index){
              graph.SYMBOLS[cs]->charge--; 
#if OPT_DEBUG
              fprintf(stderr,"  assigning %c charge %d\n",graph.SYMBOLS[cs]->ch,graph.SYMBOLS[cs]->charge); 
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
      
      if( (i < len - 1 && wln_ptr[i+1] == '&') || branch_stack.ring){
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

    // spiro is handled on dash
    case '&':
      if (pending_J_closure)
        break;
      
      if (pending_locant)
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
          on_locant += AMPERSAND_EXPAND;
          curr = ring->locants[on_locant];  
          if(!curr)
            return Fatal(i, "Error: could not fetch expanded locant position - out of range");
          
          prev = curr;
          last = curr; 
        }
      }
      else if(!branch_stack.empty())
      {
        // for implied closures:
        // If we are on the last branch that a symbol could have, a '&' also pops it off the stack and returns the next item
        // ring or not ring

        // options that happen here; 
        // 1 ) close a branch and return to the top of the branch stack
        // 2 ) we are already at the top of the branch stack, in which case we pop it off
        // 3 ) there is no symbol in the branch stack, in which case we pop a ring
        
        // check whats going on. 
        WLNSymbol *branch_top = return_object_symbol(branch_stack); 
        
        if(branch_top){
          // this must be a pop, unless its a methyl contraction
          if(prev == branch_top){
            
            switch(prev->ch){
              case 'X':
              case 'K':
                if( (prev->num_edges + prev->explicit_H) < prev->allowed_edges){
                  if(!add_methyl(prev,graph))
                    return Fatal(i, "Error: failed to add methyl group on methyl contraction");
                
                  prev = return_object_symbol(branch_stack); // if its the last one pop it 
                }
                else{ 
                  branch_stack.pop();
                  prev = return_object_symbol(branch_stack);
                }
                break;
               
              case 'Y':
                if(count_children(prev) < 3){
                  if(!add_methyl(prev,graph))
                    return Fatal(i, "Error: failed to add methyl group on methyl contraction");

                  prev = return_object_symbol(branch_stack);
                }
                else{ 
                  // we pop,
                  branch_stack.pop();
                  prev = return_object_symbol(branch_stack);
                }
                break;

              default:
                branch_stack.pop(); 
                prev = return_object_symbol(branch_stack);
            }
          }
          else // else its a simple return
            prev = branch_top; 
        }
        else if (!branch_stack.empty() && branch_stack.top().first){
          // there either must be a ring here to pop, or the branch stack is empty and we error
          branch_stack.pop(); // pop off the ring
          prev = return_object_symbol(branch_stack); // return to any open branches
          ring = branch_stack.ring; // assign the ring
        }
        else{
          return Fatal(i, "Error: popping too many rings|symbols, check '&' count");
        }
      }
      else
        return Fatal(i, "Error: popping too many rings|symbols, check '&' count");
      break;

#if MODERN 
    case '<':{
       if (pending_J_closure)
          break;
        
        bool found_close = false; 
        str_buffer.clear();
        unsigned int first_angle = i;
        int charge = 0;
        bool hypervalent = false;
        std::string ch_store; 

        // move off the first dash
        i++;
        ch = *(++wln_ptr);
        while(ch != '\0'){
          if(ch == '>'){
            found_close = true; 
            break;
          }
          
          else{
            if(ch <= 'Z' && ch >= 'A'){
              str_buffer.push_back(ch);
            }
            else if(ch >= '0' && ch <= '9')
              ch_store += ch; 
            else if(ch == '-'){
              if(ch_store.empty() && str_buffer.empty())
                hypervalent = true;
              else if (ch_store.empty())
                charge = -1; 
              else
                charge = -isNumber(ch_store); 
            }
            else if(ch == '+'){
              if(ch_store.empty())
                charge = 1;
              else
                charge = isNumber(ch_store); 
            }
            else 
              return Fatal(i, "Error: invalid character in special '< >'");
            
          }
          i++;
          ch = *(++wln_ptr);
        }
          
        if(!found_close)
          return Fatal(i, "Error: unbalanced < > characters"); 
        else if (str_buffer.empty())
          return Fatal(i, "Error: empty < > characters"); 
        
        if(str_buffer.size() == 1){
          if(hypervalent)
            curr = define_hypervalent_element(str_buffer[0],graph);
          else{
            ch = str_buffer[0]; 
            pending_charge = charge; 
            goto character_start; 
          }
          if(!curr)
            return Fatal(i, "Error: failed to define hypervalent element");
          curr->str_position = i+1;
        }
        else if(str_buffer.size() == 2){
          curr = define_element(str_buffer,graph);
          if(!curr)
            return Fatal(i, "Error: failed to define periodic element"); 
          
          curr->str_position = i+1;
        }
        else
          return Fatal(i, "Error: special '< >' must be either 1 or 2 symbols");
        
        curr->charge = charge; 
      
        if(prev){
          if(str_buffer.empty() && ring){
            if(!AddEdge(ring->locants[prev->ch], prev))
              return Fatal(i, "Error: failed to bond to previous symbol");
            edge = &prev->bond_array[prev->barr_n-1]; 
          }
          else{
            if(!AddEdge(curr, prev))
              return Fatal(i, "Error: failed to bond to previous symbol");
            edge = &prev->bond_array[prev->barr_n-1]; 
          }
            
          if(!edge)
            return Fatal(i, "Error: failed to bond to previous symbol");
          
          edge->stereo = pending_stereo; 
          pending_stereo = 0; 

          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        } 
        on_locant = '\0';
        branch_stack.push({0,curr});
        
        curr->str_position = first_angle+1+1;
        inline_unsaturate = 0;
        prev = curr;
        last = curr; 
      }
      cleared = false;
      break;


    case '>':
    case '+':
      if(pending_J_closure)
        break; 
      else
       return Fatal(i, "Error: modern read outside of known rules");
#endif 

    
    case '-':
      if (pending_J_closure)
        break;
      
      else if (pending_inline_ring && pending_ring_in_ring) // macro_closure
      { 
        
        // onlocant holds the char needed to wrap the ring back, 
        if(!wrap_ring)
          return Fatal(i, "Error: wrap ring is not active");

        curr = wrap_ring->locants[on_locant];
        if(!curr)
          return Fatal(i, "Error: cannot access looping ring structure");
      
        if(prev){  

          if(!AddEdge(curr, prev))
            return Fatal(i, "Error: failed to bond to previous symbol");

          edge = &prev->bond_array[prev->barr_n-1]; 
          edge->stereo = pending_stereo; 
          pending_stereo = 0;

          wrap_ring->macro_return = edge; 
          if(inline_unsaturate){
            if(!unsaturate_edge(edge,inline_unsaturate))
              return Fatal(i, "Error: failed to unsaturate bond"); 
            inline_unsaturate = 0;
          }
        }
        else
          return Fatal(i,"Error: no previous symbol for inline ring defintion");
        
        if(i+3 < len && wln_ptr[i+3] == '-')
          i+=3;
        else if(i+2 < len && wln_ptr[i+2] == '-')
          i+=2;
        else
          return Fatal(i, "Error: macro-notation requires closure with the ring size in two dashes e.g -6-");

        curr = prev; // set back to prev
        on_locant = '\0';
        pending_ring_in_ring = false;
        pending_inline_ring = false;
        pending_rir_closure = true;
      }
      else{
        
        // conditions for '-' 
        //
        //  1) closed at position +3 away -NA- - elemental definition
        //  2) closed at position +2 away -K- - hypervalent element
        //  3) branch at pos +1 and space at +2 = spiro
        //  else) not closed - inline ring
        //
        //  1/2*) if +1 = '&'|'-' branching locant defintion e.g E--P, E-&P
        //  must eliminate the spaces every time
        //
        //  Not specified here, macro opens e.g L- T- are handled on L/T 

        // 1)
        if( i + 3 < len && wln_ptr[i+3] == '-' && wln_ptr[i+1] != ' '){

          curr = define_element(wln_ptr[i+1],wln_ptr[i+2],graph);
          if(!curr)
            return Fatal(i, "Error: failed to define periodic element"); 
          
          if(on_locant == '0'){
            curr->charge++;
            //on_locant = '\0';
          }

          if(prev){
            if(!AddEdge(curr, prev))
              return Fatal(i, "Error: failed to bond to previous symbol");
            edge = &prev->bond_array[prev->barr_n-1];


            if(!edge)
              return Fatal(i, "Error: failed to bond to previous symbol");
            
            edge->stereo = pending_stereo; 
            pending_stereo = 0; 

            if(inline_unsaturate){
              if(!unsaturate_edge(edge,inline_unsaturate))
                return Fatal(i, "Error: failed to unsaturate bond"); 
              inline_unsaturate = 0;
            }

          }

          on_locant = '\0';
          branch_stack.push({0,curr});
          curr->str_position = i+2;
          inline_unsaturate = 0;
          prev = curr;
          last = curr; 
          i+=3; 
        }
        
        // 2 and 3 
        else if(i+2 < len && wln_ptr[i+2] =='-' && wln_ptr[i+1] != ' '){
          
          curr = define_hypervalent_element(wln_ptr[i+1],graph);
          if(!curr)
            return Fatal(i, "Error: failed to define hypervalent element");
          
          if(prev){
            if(!AddEdge(curr, prev))
              return Fatal(i, "Error: failed to bond to previous symbol");
            edge = &prev->bond_array[prev->barr_n-1];
            if(!edge)
              return Fatal(i, "Error: failed to bond to previous symbol");
            
            edge->stereo = pending_stereo; 
            pending_stereo = 0; 

            if(inline_unsaturate){
              if(!unsaturate_edge(edge,inline_unsaturate))
                return Fatal(i, "Error: failed to unsaturate bond"); 
              inline_unsaturate = 0;
            }
          }

          on_locant = '\0';
          branch_stack.push({0,curr});
          curr->str_position = i+2;
          inline_unsaturate = 0;
          prev = curr;
          last = curr; 

          i+=2; 
        }

        // 3)
        else if( i + 2 < len && wln_ptr[i+1] == '&' && wln_ptr[i+2] == ' '){
          pending_spiro = true;
          pending_inline_ring = true;
          i+=1; 
        }
        else{

          if(pending_inline_ring)
            return Fatal(i,"Error: previous in-line ring definition not finished\n");

          pending_inline_ring = true;
          return_object_symbol(branch_stack);
          if(branch_stack.branch && !prev){
            // prev must be at top of the branch stack
            while(branch_stack.top().second != branch_stack.branch)
              branch_stack.pop();
            
            prev = return_object_symbol(branch_stack);
          }    
        }
      }
      cleared = false;
      break;
    
    
    case '/':
      if (pending_J_closure){
        locant_skips = 2;
        break;
      }
      else if(pending_numbers){ // state that this must be a charge 

        if(!cleared)
          return Fatal(i, "Error: opening post charge assignment without proper syntax [ &x/x ]");
        
        int positive_index = 0; 
        if(d3){
          positive_index += (100*(d1-'0')); 
          positive_index += (10*(d2-'0')); 
          positive_index += (1*(d3-'0')); 
        }
        else if(d2){
          positive_index += (10*(d1-'0')); 
          positive_index += (1*(d2-'0')); 
        }      
        else
          positive_index += (1*(d1-'0')); 
      
        d1=0;
        d2=0;
        d3=0;
       

#if OPT_DEBUG
        fprintf(stderr,"  attempting +1 charge on index %d\n",positive_index); 
#endif

        if (positive_index != 0){
          // find the symbol and increment its charge + 1
          bool found = false;
          for(unsigned int cs = 0;cs<graph.symbol_count;cs++){
            if(graph.SYMBOLS[cs]->str_position == (unsigned int)positive_index){
              graph.SYMBOLS[cs]->charge++;
              found = true;
              break;
            }
          }
          
          if(!found)
            return Fatal(i, "Error: positive charge index out of range, check letter index");
        }

        pending_numbers = false;
        pending_negative_charge = true;
      }
      else
        return Fatal(i,"Error: multipliers are not currently supported");
      cleared = false;
      break;

    default:
      return Fatal(i,"Error: unallowed character! - alphabet: [A-Z][0-1][&-/' ']");
    }

    ch = wln_ptr[++i];
  }


  if(pending_numbers){        
    int carbon_len = 0;
    if(d3){
      carbon_len += (100*(d1-'0')); 
      carbon_len += (10*(d2-'0')); 
      carbon_len += (1*(d3-'0')); 
    }
    else if(d2){
      carbon_len += (10*(d1-'0')); 
      carbon_len += (1*(d2-'0')); 
    }      
    else
      carbon_len += (1*(d1-'0')); 
  
    d1=0;
    d2=0;
    d3=0;

    if(carbon_len < 0)
      return Fatal(i, "Error: non-numeric value entered for carbon length");
    else if (carbon_len > 100)
      return Fatal(i,"Error: creating a carbon chain > 100 long, is this reasonable for WLN?");
    

    // create the head 
    WLNSymbol *curr = AllocateWLNSymbol('#', graph);
    curr->str_position = i;
    curr->special = std::to_string(carbon_len); 
    curr->allowed_edges = 6; // can hold a triple bond on either side
    
    if(prev){

      if(!AddEdge(curr, prev))
        return Fatal(i, "Error: failed to bond to previous symbol");
      edge = &prev->bond_array[prev->barr_n-1]; 
      edge->stereo = pending_stereo; 
      pending_stereo = 0; 
      
      if(inline_unsaturate){
        if(!unsaturate_edge(edge,inline_unsaturate))
          return Fatal(i, "Error: failed to unsaturate bond"); 
        inline_unsaturate = 0;
      }
    }

    pending_numbers = false;
    prev = curr; 
    last = curr; 

  }

  // single letter methyl branches
  if(on_locant && on_locant != '0' && !pending_inline_ring && !branch_stack.empty()){
    if(!add_methyl(branch_stack.ring->locants[on_locant],graph))
      return Fatal(i, "Error: could not attach implied methyl to ring");
    
    on_locant = 0;
  }
  
  if(pending_negative_charge){
    int negative_index = 0; 
    if(d3){
      negative_index += (100*(d1-'0')); 
      negative_index += (10*(d2-'0')); 
      negative_index += (1*(d3-'0')); 
    }
    else if(d2){
      negative_index += (10*(d1-'0')); 
      negative_index += (1*(d2-'0')); 
    }      
    else
      negative_index += (1*(d1-'0')); 
  
    d1=0;
    d2=0;
    d3=0;

    if(negative_index < 0)
      return Fatal(i, "Error: assigning non-numerical value to charge index");
    else if (negative_index != 0){
      // find the symbol and increment its charge + 1
      bool found = false;
      for(unsigned int cs = 0;cs<graph.symbol_count;cs++){
        if(graph.SYMBOLS[cs]->str_position == (unsigned int)negative_index){
          graph.SYMBOLS[cs]->charge--; 
#if OPT_DEBUG
          fprintf(stderr,"  assigning %c charge %d\n",graph.SYMBOLS[cs]->ch,graph.SYMBOLS[cs]->charge); 
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
  fprintf(fp, "  rankdir = LR;\n");
  for (unsigned int i=0; i< graph.symbol_count;i++)
  {
    WLNSymbol *node = graph.SYMBOLS[i];
    if(!node)
      continue;

    fprintf(fp, "  %d", node->id);
    if (node->ch == '*' || node->ch == '#')
      fprintf(fp, "[shape=circle,label=\"*:%s\"];\n", node->special.c_str());
    else if (node->spiro)
      fprintf(fp, "[shape=circle,label=\"%c\",color=blue];\n", node->ch);
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

  BabelGraph(){};
  ~BabelGraph(){};


  OBAtom* NMOBMolNewAtom(OBMol* mol, unsigned int elem,int charge,unsigned int hcount)
  {
    OBAtom* result = mol->NewAtom();
    if(!result)
      return 0;

    result->SetAtomicNum(elem);
    result->SetFormalCharge(charge);
    result->SetImplicitHCount(hcount);
    return result;
  }


  OBBond* NMOBMolNewBond( OBMol* mol,
                          OBAtom* s,
                          OBAtom* e,
                          unsigned int order)
  {
    
    OBBond* bptr = 0; 
    if(!s || !e){
#if ERRORS == 1
      fprintf(stderr,"Error: could not find atoms in bond, bond creation impossible\n");
#endif      
      return bptr;
    }

    if (!mol->AddBond(s->GetIdx(), e->GetIdx(), order)){

#if ERRORS == 1
      fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n",s->GetIdx(),e->GetIdx());
#endif
      return bptr;
    }
        
    bptr = mol->GetBond(mol->NumBonds() - 1);
    return bptr;
  }
  
  /* creates a carbon chain, from an inited node. */
  OBAtom *OBMolCarbonChain(OBMol *mol, OBAtom *head, unsigned int size){
    OBAtom *prev = head; 
    OBAtom *carbon = 0; 
    for(unsigned int i=0;i<size;i++){
      carbon = NMOBMolNewAtom(mol, 6, 0, 2); 
      NMOBMolNewBond(mol, prev, carbon, 1); 
      prev = carbon; 
    }

    carbon->SetImplicitHCount(0); 
    return carbon; // head of the chain
  }


  void NMOBSanitizeMol(OBMol* mol)
  {
    // WLN has no inherent stereochemistry, this can be a flag but should be off by default
#if MODERN
    //StereoFrom2D(mol); 
    mol->SetAromaticPerceived(false);
    mol->SetChiralityPerceived(true);
#else
    mol->SetChiralityPerceived(true);
    mol->SetAromaticPerceived(false);
#endif 

    mol->DeleteHydrogens();
  }
  
  OBAtom *WLNSymbolToAtom(OBMol* mol, WLNSymbol*sym){
    int charge = 0; 
    unsigned int atomic_num = 0;
    unsigned int hcount = sym->explicit_H;

    switch(sym->ch){
      case 'H':
        atomic_num = 1;
        hcount = 0;
        charge = sym->charge; 
        break; 

      case 'B':
        atomic_num = 5;
        charge = sym->charge; 
        break;

      case '1': 
      case 'C':
        atomic_num = 6;
        hcount = 4- sym->num_edges;
        charge = sym->charge; 
        if(charge < 0)
          hcount+= charge; 
        break; 

        if(sym->charge < 0 && hcount > 0)
          hcount--; 
        break;
      
      case 'X':
        atomic_num = 6;
        charge = sym->charge; 
        break;

      case 'Y':
        atomic_num = 6;
        charge = sym->charge; 
        hcount = 4 - sym->num_edges + sym->explicit_H;
        if(charge > 0){
          for(unsigned int c=0;c<charge;c++)
            if(hcount > 0)
              hcount--; 
        }
        break;

      case 'N':
        atomic_num = 7;
        if(sym->inRing)
          sym->allowed_edges = 3;
        charge = sym->charge; 

        if(!hcount && sym->aromatic && sym->num_edges < sym->allowed_edges && sym->charge==0)
          hcount = 1; // allows allow 1 H in aromatic species 
        break;

      case 'M':
        atomic_num = 7;
        charge = sym->charge; 
        break;

      case 'Z':
        atomic_num = 7; 
        charge = sym->charge;
        break;

      case 'K':
        atomic_num = 7;
        charge = 1; 
        break;

      case 'O':
        atomic_num = 8;
        if(!sym->charge){
          if(sym->num_edges==1 && !hcount)
            charge = -1;
          if(!sym->num_edges && !hcount)
            charge = -2;
        }
        else {
          charge = sym->charge; 
        }
        break;

      case 'Q':
        atomic_num = 8;
        if(!sym->num_edges && hcount == 1)
          charge = -1;
        break;

      case 'F':
        atomic_num = 9;
        if(!sym->num_edges && !hcount)
          charge = -1;
        break;
      
      case 'P':
        atomic_num = 15; 
        charge = sym->charge; 
        if(!hcount && sym->aromatic && sym->num_edges < 3)
          hcount = 3-sym->num_edges; // allows allow 1 H in aromatic species 
        break;
      
      case 'S':
        atomic_num = 16;
        charge = sym->charge; 
        if(!hcount && sym->aromatic && sym->num_edges < 2)
          hcount = 2-sym->num_edges; // allows allow 1 H in aromatic species 
        break;

      case 'G':
        atomic_num = 17;
        if(!sym->num_edges && !hcount)
          charge = -1;
        break;

      case 'E':
        atomic_num = 35;
        if(!sym->num_edges && !hcount)
          charge = -1;
        break;

      case 'I':
        atomic_num = 53;
        if(!sym->num_edges && !hcount)
          charge = -1;
        break;
    
      case '*':
        atomic_num = special_element_atm(sym->special);
        charge = sym->charge; 
        break;
      
      case '#': // should make a dummy atom
        hcount = 0; 
        break;

      default:
        return 0;
    }
    
    OBAtom *atom = NMOBMolNewAtom(mol,atomic_num,charge,hcount);
    return atom; 
  }


  bool ConvertFromWLN(OBMol* mol,WLNGraph &graph, unsigned int len){
    if(OPT_DEBUG)
      fprintf(stderr,"Converting wln to obabel mol object: \n");
    // doing this depth first might make more sense now
    
    std::map<WLNSymbol*,std::pair<OBAtom*,OBAtom*>> chain_pairs; // ooffh, STL is probably super slow 
    std::map<WLNSymbol*, OBAtom*> atom_map;

    for (unsigned int i=0; i<graph.symbol_count;i++){
      WLNSymbol *sym = graph.SYMBOLS[i];
      
      // handle any '#' carbon chains
      if(sym->ch == '#' && isNumber(sym->special) > 1 ){
        OBAtom *chain_head = NMOBMolNewAtom(mol, 6, 0, 0); // might need to be upped if terminal
        OBAtom *chain_end = OBMolCarbonChain(mol, chain_head, isNumber(sym->special)-1);  
        
        int order = 3; 
        for(unsigned int h=0;h<sym->parr_n;h++){
          order -= sym->prev_array[h].order;
        }
        if(order >= 0)
          chain_head->SetImplicitHCount(order); 
       
        // for bonds going forward, we assume head bond, for backward end bond, 
        order = 3; // assume a single bond is made to the chain 
        for(unsigned int h=0;h<sym->barr_n;h++)
          order -= sym->bond_array[h].order;
        if(order >= 0)
          chain_end->SetImplicitHCount(order); 
        
        chain_pairs[sym] = {chain_head,chain_end}; 
      }
      else if (sym->ch == '#'){
        sym->ch = '1';
        OBAtom *atom = WLNSymbolToAtom(mol, sym);
        atom_map[sym] = atom;  
      }
      else{
        OBAtom *atom = WLNSymbolToAtom(mol, sym);
        atom_map[sym] = atom;  
      }
    }

    // create edges 
    for(unsigned int i=0;i<graph.symbol_count;i++){
      WLNSymbol *parent = graph.SYMBOLS[i];
      if(parent->barr_n){
        for (unsigned int ei=0;ei<parent->barr_n;ei++){
          WLNEdge *e = &parent->bond_array[ei];
          WLNSymbol *child = e->child;
          unsigned int bond_order = e->order;  
          
          OBAtom *catom = 0;
          OBAtom *patom = 0; 
          if(parent->ch == '#')
            patom = chain_pairs[parent].second;
          
          if(child->ch == '#')
            catom = chain_pairs[child].first; 

          if(!catom)
            catom = atom_map[child];
          if(!patom)
            patom = atom_map[parent]; 

          OBBond *bptr = NMOBMolNewBond(mol,patom,catom,bond_order);
          
          if(e->stereo==1){
            bptr->SetHash();
          }
          else if (e->stereo == 2){
            bptr->SetWedge(); 
          }
          
          if(bptr->IsWedge())
            fprintf(stderr,"stereo - set\n"); 

          if(!bptr)  
            return false;
        }
      }
    }



    return true;
  }

};



/**********************************************************************
                          Canonical Algorithms
**********************************************************************/


unsigned int methyl_contract(WLNSymbol *sym){
  // methyl contraction 
  if(sym->special != "1")
    return 0; 

  if(sym->barr_n == 0 && sym->parr_n == 1 && sym->prev_array[0].order == 1){
    switch(sym->prev_array[0].child->ch){
      case 'X':
      case 'Y':
      case 'K':
        return 1; 
        break;

      default:
        return 0; 
    }
  }
  else if (sym->barr_n == 1 && sym->parr_n == 0 && sym->bond_array[0].order==1){
    switch(sym->bond_array[0].child->ch){
      case 'X':
      case 'Y':
      case 'K':
        return 1; 
      default:
        return 0; 
    }
  }

  return 0; 
}


// unfortuantely quite expensive, dioxo will always be forward facing
// 0 - none, 1 = (=O)(=O), 2 = ([O-])(=O)
// if the symbol is not saturated and type 2, return 0, as this cannot be implied
unsigned int check_dioxo_type(WLNSymbol *node, std::map<WLNSymbol*,bool> &seen_symbols, std::string &buffer){
  WLNSymbol* double_oxygen_1 = 0; 
  WLNSymbol* double_oxygen_2 = 0; 
  WLNSymbol *oxo_ion = 0; 

  for(unsigned int ei=0;ei<node->barr_n;ei++){
    if( node->bond_array[ei].child->ch == 'O'){
      if(node->bond_array[ei].child->charge == -1 && node->bond_array[ei].child->num_edges == 1
          && !oxo_ion &&  !seen_symbols[node->bond_array[ei].child]){
       
        oxo_ion = node->bond_array[ei].child;
      }
      else if (node->bond_array[ei].order == 2){
        if(!double_oxygen_1 &&  !seen_symbols[node->bond_array[ei].child])
          double_oxygen_1 = node->bond_array[ei].child;
        else if (double_oxygen_1 && !seen_symbols[node->bond_array[ei].child])
          double_oxygen_2 = node->bond_array[ei].child;
      }
    }
  }

  
  // fprintf(stderr,"Atom %c: %p %p %p\n",node->ch, double_oxygen_1,double_oxygen_2,oxo_ion); 

  // need to flag that this charge is no longer needed, as implicitly given
  // buffer type is a resetable flag
  if(double_oxygen_1 && double_oxygen_2){
    seen_symbols[double_oxygen_1] = true;
    seen_symbols[double_oxygen_2] = true;
    buffer += 'W';
    double_oxygen_1->str_position = buffer.size(); 
    double_oxygen_2->str_position = buffer.size(); 
    return 1;
  }
  else if (double_oxygen_1 && oxo_ion && node->num_edges == node->allowed_edges){
    seen_symbols[double_oxygen_1] = true;
    seen_symbols[oxo_ion] = true;
    buffer += 'W';
    double_oxygen_1->str_position = buffer.size(); 
    oxo_ion->str_position = buffer.size(); 
    return 2; 
  }
  else
    return 0;
}

/* used to perform a radix style sort to order bond stack pushes */
struct ChainScore{
  WLNEdge *e;
  std::string chunk; 
  bool terminates;
  unsigned int ring_ranking; 
  bool has_branch; 
};

struct SortedEdges{
  WLNEdge *edges[MAX_EDGES*2];
  unsigned char e_n;
  unsigned char e_max;
};


void debug_score(ChainScore *score){
  fprintf(stderr,"%c --> %c: %s",score->e->parent->ch, score->e->child->ch,score->chunk.c_str()); 
  fprintf(stderr,"\tterm:%d, branch: %d, ring:%d\n",score->terminates,score->has_branch,score->ring_ranking); 
}


void SortByTerminal(ChainScore **arr,unsigned int len){
  for (unsigned int j=1;j<len;j++){
    unsigned int key = 0; 
    
    ChainScore *s = arr[j];
    key = s->terminates;
    
		int i = j-1;
    while(i>=0){
      unsigned int val = 0;      
      val = arr[i]->terminates;
      if(val <= key)
        break;

      arr[i+1] = arr[i];
      i--;
    }
		arr[i+1] = s;
	}
}

void SortByBranch(ChainScore **arr,unsigned int len){
  for (unsigned int j=1;j<len;j++){
    unsigned int key = 0; 
    
    ChainScore *s = arr[j];
    key = s->has_branch;
    
		int i = j-1;
    while(i>=0){
      unsigned int val = 0;      
      val = arr[i]->has_branch;
      if(val <= key)
        break;

      arr[i+1] = arr[i];
      i--;
    }
		arr[i+1] = s;
	}
}

void SortByRing(ChainScore **arr,unsigned int len){
  for (unsigned int j=1;j<len;j++){
    unsigned int key = 0; 
    
    ChainScore *s = arr[j];
    key = s->ring_ranking;
    
		int i = j-1;
    while(i>=0){
      unsigned int val = 0;      
      val = arr[i]->ring_ranking;
      if(val >= key) // this does opposite ordering 
        break;

      arr[i+1] = arr[i];
      i--;
    }
		arr[i+1] = s;
	}
}

// sort based on the letter ordering, only if symbol lens match
void SortByRule2(ChainScore **arr,unsigned int len){
  for (unsigned int j=1;j<len;j++){
    ChainScore *s = arr[j];
		int i = j-1;
    while(i>=0){
      unsigned int k=0;
      bool _break = false;
      
      // if one string is longer, global rule 2 will sort the strings
      while(k<s->chunk.size() && k<arr[i]->chunk.size()){
        if(s->chunk[k] != arr[i]->chunk[k]){
          if(arr[i]->chunk[k] < s->chunk[k]){
            _break = true;
          }
          break;
        }

        k++;
      }
      if(_break)
        break;
      arr[i+1] = arr[i];
      i--;
    }
		arr[i+1] = s;
	}
}


/* run the chain until either a ring atom/branch point/EOC is seen,
 * prototype does have a map copy which is irrating */
ChainScore *RunChain(WLNEdge *edge,std::map<WLNSymbol*,bool> seen){
  
  ChainScore *score = new ChainScore; // must use new for string 
  score->e = edge;
  score->chunk = ""; 
  score->terminates = 0;
  score->ring_ranking = 0; // these should be placed last when possible 
  score->has_branch = 0; // also takes a lower priority?
  
  WLNSymbol *node = edge->child; 
  // no stack needed as looking ahead in a linear fashion
  for(;;){
    seen[node] = true;
    if(node->inRing){
      score->ring_ranking = node->inRing->ranking;
      return score; // immediate
    }

    switch(node->ch){
    case '#':
      score->chunk += node->special;  
      break;

      // these must branch
      case 'Y':
      case 'X':
      case 'K':
        score->chunk += node->ch;
        score->has_branch = true;
        return score;

      // these could branch, only split if they have
      case '*':
      case 'P':
      case 'S':
      case 'B':
      case 'N': // further canonical rule to add to minimise this
        score->chunk += node->ch;
        score->has_branch = true;
        return score; 
      
      // terminators end the branch lookahead immediately
      case 'E':
      case 'F':
      case 'G':
      case 'I':
      case 'Q':
      case 'Z':
        score->chunk += node->ch;
        score->terminates = true;
        return score; 

      default:
        score->chunk += node->ch; 
    }
   
    if(node->barr_n && !seen[node->bond_array[0].child]){
      for(unsigned int i=1;i<node->bond_array[0].order;i++) // let symbol orderer sort here
        score->chunk += 'U';

      node = node->bond_array[0].child; // no need to iterate, it should only have 1.
    }
    else if (node->parr_n && !seen[node->prev_array[0].child]){
      for(unsigned int i=1;i<node->prev_array[0].order;i++) // let symbol orderer sort here
        score->chunk += 'U';

      node = node->prev_array[0].child; // no need to iterate, it should only have 1.
    }
    else 
      return score;

  }

  return score; 
}


// forward declaration 
bool CanonicalWLNChain(WLNSymbol *node, WLNGraph &graph, WLNSymbol *ignore, std::string &buffer); 
bool CanonicalWLNRing(WLNSymbol *node, WLNGraph &graph, WLNSymbol *ignore, std::string &buffer);


/* sorts the bonds of any given chain, function will look through both previous and forward bonds in order to continue the parse, 
as such, it requires a seen map to avoid looping back to areas its previously been  */
SortedEdges* ArrangeBonds( WLNSymbol *sym, std::map<WLNSymbol*,bool> &seen, WLNSymbol *ignore)
{
  SortedEdges *se = (SortedEdges*)malloc(sizeof(SortedEdges));
  for(unsigned int i=0;i<MAX_EDGES*2;i++)
    se->edges[i] = 0; 

  unsigned int l = 0;
  ChainScore *scores[64] = {0};
  
  for(unsigned int ei=0;ei<sym->barr_n;ei++){
    WLNEdge *e =  &sym->bond_array[ei]; 
    if(!seen[e->child] && e->child != ignore)
      scores[l++] = RunChain(e,seen); // score each chain run
  }

  for(unsigned int ei=0;ei<sym->parr_n;ei++){
    WLNEdge *e =  &sym->prev_array[ei];
    if(!seen[e->child] && e->child != ignore)
      scores[l++] = RunChain(e,seen); // score each chain run
  }

  SortByRule2(scores, l); 
  SortByBranch(scores, l); 
  SortByTerminal(scores, l); // sort by terminals, prefer them
  SortByRing(scores, l); 
  
  unsigned char a = 0; 
  for(int i=l-1;i>=0;i--){ // sort the chains (radix style) to get high priorities first
    se->edges[a++] = scores[i]->e; 
    delete scores[i];
  }
    

  // this all valence full by definition
  if(sym->ch != 'X' && sym->ch != 'K' && sym->ch != 'Y'){
    if(IsBranching(sym) && sym->num_edges < sym->allowed_edges)
      scores[a++] = 0; // adds a forcable pop 
  }

  se->e_n = 0; 
  se->e_max = a;
  return se; 
}

SortedEdges* ArrangeRingBonds( WLNSymbol *locant,WLNRing *ring, std::map<WLNSymbol*,bool> seen,WLNSymbol *ignore)
{

  SortedEdges *se = (SortedEdges*)malloc(sizeof(SortedEdges));

  unsigned int l = 0;
  ChainScore *scores[64] = {0};
  
  for (unsigned int ei=0;ei<locant->barr_n;ei++){
    WLNEdge *fe = &locant->bond_array[ei];


    if( fe->child->inRing != ring 
        && fe->child != ignore 
        && !seen[fe->child])
    {
      if(fe->child->ch == 'H' && !locant->explicit_H)
        ;
      else
        scores[l++] = RunChain(fe, seen);
    }
  }

  for (unsigned int ei=0;ei<locant->parr_n;ei++){
    WLNEdge *be = &locant->prev_array[ei];
    if(be->child->inRing != ring && be->child != ignore && !seen[be->child]){

      if(be->child->ch == 'H' && !locant->explicit_H)
        ;
      else
        scores[l++] = RunChain(be, seen);
    }
  }

  SortByRule2(scores, l); 
  SortByBranch(scores, l); 
  SortByRing(scores, l); 

  // sort by ring will cluster the spiro rings together, remove 1
  if(locant->spiro){
    for(unsigned int i=0;i<l-1;i++){
      if(scores[i]->ring_ranking == scores[i+1]->ring_ranking){
        delete scores[i];
        scores[i] = 0;
        break;
      }
    }
  }

  unsigned char a = 0; 
  for(int i=l-1;i>=0;i--){ // sort the chains (radix style) to get high priorities first
    if(scores[i]){
      se->edges[a++] = scores[i]->e; 
      delete scores[i];
    }
  }

  se->e_n = 0; 
  se->e_max = a; 
  return se; 
}

void static write_locant(unsigned char locant,std::string &buffer){
  if(locant < 'X')
    buffer += locant;
  else{
    unsigned int amps = 0;
    while(locant >= 'X'){
      amps++; 
      locant+= -AMPERSAND_EXPAND;
    }
    buffer += locant; 
    for(unsigned int i=0;i<amps;i++)
      buffer += '&';
  }
}


// pass through a methyl contraction 

void WriteCharacter(WLNSymbol *sym, std::string &buffer, WLNGraph &graph){
  unsigned int modifier = 0; 
  switch (sym->ch) {
    
    case '#': 
      sym->str_position = buffer.size()+1; // place on first number 
      buffer += sym->special;
      break;

    case '*':
      buffer += '-';
      sym->str_position = buffer.size()+1; // place on first letter
      buffer += sym->special;
      buffer += '-';
      break;

    case 'E':
    case 'F':
    case 'G':
    case 'H':
    case 'I':
      if(sym->allowed_edges>1){
        buffer += '-';
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
        buffer += '-';
      }
      else{
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
      }
      break;

    case 'O':
      if(sym->allowed_edges>2){
        buffer += '-';
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
        buffer += '-';
      }
      else{
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
      }
      break;

    case 'B':
      if(sym->allowed_edges>3){
        buffer += '-';
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
        buffer += '-';
      }
      else{
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
      }
      break;
    
    case 'N':
      if(sym->allowed_edges>3){
        buffer += '-';
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
        buffer += '-';
      }
      else{
        buffer += sym->ch;
        sym->str_position = buffer.size(); 
      }
      break;

    case 'M':
      buffer += sym->ch;
      sym->str_position = buffer.size(); 
      modifier = 1; 
      break;

    case 'Q':
      buffer += sym->ch;
      sym->str_position = buffer.size(); 
      modifier = 1; 
      break;

    case 'Z':
      buffer += sym->ch;
      sym->str_position = buffer.size(); 
      modifier = 2; 
      break;
    
    case 'c':
      buffer += 'C';
      sym->str_position = buffer.size(); 
      break;

    default:
      buffer += sym->ch;
      sym->str_position = buffer.size(); 
  };

  for(unsigned int i=modifier;i<sym->explicit_H;i++)
    buffer += 'H'; 

}



bool CanonicalWLNChain(WLNSymbol *node, WLNGraph &graph, WLNSymbol *ignore, std::string &buffer){
  std::map<WLNSymbol*,SortedEdges*> sorted_edges;
  std::set<SortedEdges*> se_memory_hold; 
  std::map<WLNSymbol*,bool> seen_symbols; 
  std::stack<WLNSymbol*> branching_symbols;

  WLNEdge *edge = 0;
  SortedEdges *se = 0; 

  WriteCharacter(node,buffer,graph);
  unsigned int dioxo_write = check_dioxo_type(node,seen_symbols,buffer); 

  seen_symbols[node] = true;
  sorted_edges[node] = ArrangeBonds(node, seen_symbols, ignore);

  if(IsBranching(node)){
    if(dioxo_write == 1)
      sorted_edges[node]->edges[sorted_edges[node]->e_max++] = 0; // adds pop since W will go from 4 to 3  
    
    if(dioxo_write==2 && node->allowed_edges==4){
      dioxo_write = 0;
    }
    else
      branching_symbols.push(node);
  }

  for(;;){
    if(!node)
      break;
    

    se = sorted_edges[node];
    se_memory_hold.insert(se); 
    if(se->e_max > 0 && se->e_n < se->e_max){
      edge = se->edges[se->e_n++];

      if(!edge){
        if(!branching_symbols.empty()){
          buffer += '&';
          branching_symbols.pop(); 
        }

        while(  !branching_symbols.empty() && 
                sorted_edges[branching_symbols.top()]->e_n == sorted_edges[branching_symbols.top()]->e_max)
        {
          branching_symbols.pop(); 
        }

        if(!branching_symbols.empty())
          node = branching_symbols.top();
        else 
          node = 0; 
      }
      else{ 
        for(unsigned int i=1;i<edge->order;i++)
          buffer += 'U';

        // benzene special case; 
        if(edge->child->inRing && edge->child->inRing->str_notation == "L6J"){
          // incoming == 'A' locant need to be rearranged          
          
          // locants are rotated from the incoming position, C = A, therefore D = B etc
          unsigned char incoming_char = edge->child->inRing->locants_ch[edge->child]; 
          
          std::vector<std::pair<WLNSymbol*,unsigned char>> new_positions; 
          if(incoming_char != 'A'){
          // wrap the locants around 
            WLNRing * benzene =edge->child->inRing;
            for (unsigned char ch = 'A'; ch <= 'F'; ch++){
              if(ch < incoming_char){
                unsigned char new_loc = (ch-incoming_char) + 'F'+1;   
                new_positions.push_back({benzene->locants[ch],new_loc}); 
              }
              else{
                unsigned char new_loc = (ch-incoming_char) + 'A';   
                new_positions.push_back({benzene->locants[ch],new_loc}); 
              }
            }
            
            // seems a bit like bullying, but does work
            benzene->locants.clear();
            benzene->locants_ch.clear(); 
            for (std::pair<WLNSymbol*,unsigned char> p : new_positions){
              benzene->locants[p.second] = p.first;
              benzene->locants_ch[p.first] = p.second; 
            }
          }

          CanonicalWLNRing(edge->child, graph, edge->parent, buffer);
          while(  !branching_symbols.empty() && 
                  sorted_edges[branching_symbols.top()]->e_n == sorted_edges[branching_symbols.top()]->e_max)
          {
            //buffer += '&';
            branching_symbols.pop(); 
          }

          if(!branching_symbols.empty())
            node = branching_symbols.top();
          else
            node = 0; 
        }
        else if(edge->child->inRing){
          
          buffer += '-';
          buffer += ' ';
          write_locant(edge->child->inRing->locants_ch[edge->child], buffer); 
          
          if(edge == edge->child->inRing->macro_return || edge->reverse == edge->child->inRing->macro_return)
            buffer += "-x-J";
          else 
            CanonicalWLNRing(edge->child, graph, edge->parent, buffer);

          
          while(  !branching_symbols.empty() && 
                  sorted_edges[branching_symbols.top()]->e_n == sorted_edges[branching_symbols.top()]->e_max)
          {
            //buffer += '&';
            branching_symbols.pop(); 
          }

          if(!branching_symbols.empty())
            node = branching_symbols.top();
          else
            node = 0; 
        }
        else{
          node = edge->child;
          if(methyl_contract(node)){
            buffer += '&'; 
            seen_symbols[node] = true;
            sorted_edges[node] = ArrangeBonds(node, seen_symbols, ignore);
          }
          else{
            WriteCharacter(node, buffer,graph); 
            unsigned int dioxo_write = check_dioxo_type(node,seen_symbols,buffer); 

            seen_symbols[node] = true;
            sorted_edges[node] = ArrangeBonds(node, seen_symbols, ignore);
            
            if(IsBranching(node)){
              if(dioxo_write == 1)
                sorted_edges[node]->edges[sorted_edges[node]->e_max++] = 0; // adds pop since W will go from 4 to 3  
              
              if(dioxo_write==2 && node->allowed_edges==4){
                dioxo_write = 0;
              }
              else
                branching_symbols.push(node);
            }
          }
        }
      }
    }
    else if (!branching_symbols.empty()){
      // return to something on the branch stack here

      while(  !branching_symbols.empty() && 
              sorted_edges[branching_symbols.top()]->e_n == sorted_edges[branching_symbols.top()]->e_max)
      {
        branching_symbols.pop(); 
      }

      if(!branching_symbols.empty()){
        if( !(IsTerminator(node) || node->ch == 'W' || (!buffer.empty() && buffer.back() == 'W') ) ){ 
            
          if(!methyl_contract(node))
            buffer+= '&';
        }

        node = branching_symbols.top();
      }
      else
        node = 0; 
    }
    else
      node = 0; 
  }
  
  for(std::set<SortedEdges*>::iterator miter=se_memory_hold.begin();miter!=se_memory_hold.end();miter++){
    SortedEdges *free_edge = (*miter);
    free(free_edge);
    free_edge = 0; 
  }
  
  return true; 
}

bool CanonicalWLNRing(WLNSymbol *node, WLNGraph &graph,WLNSymbol *ignore, std::string &buffer){

  // expect the node to be within a ring, fetch ring and write the cycle
  WLNRing *ring = node->inRing; 
  
  if(ring->macro_return){
    buffer += ring->str_notation[0];
    buffer += '-'; 
  }
 
  // reaassign the positions for charges within the ring --
  //
  for(unsigned char ch = 'A';ch < 'A'+ring->rsize;ch++){
    ring->locants[ch]->str_position = buffer.size() + ring->position_offset[ring->locants[ch]] + 1; 
  } 

  if(!strcmp(ring->str_notation.c_str(),"L6J"))
    buffer += 'R'; // incoming locant will always be in position A, requires a bit of a re-jig in order to work   
  else
    buffer += ring->str_notation; 

  // add all the locants to the seen map, prevents back tracking within the ring structure
  std::map<WLNSymbol*,bool> seen_locants; 
  for(unsigned char ch = 'A';ch < 'A'+ring->rsize;ch++){ 
    // add the hydrogens here
    WLNSymbol *lc = ring->locants[ch]; 
    seen_locants[lc] = true; 
    switch(lc->ch){
      case '1':
        break; // ring carbons

      case 'M':
        for(unsigned int h=1;h<lc->explicit_H;h++){
          buffer += " ";
          write_locant(ch, buffer); 
          buffer += 'H';
        }
        break;


      case 'Z':
        for(unsigned int h=2;h<lc->explicit_H;h++){
          buffer += " ";
          write_locant(ch, buffer); 
          buffer += 'H';
        }
        break;

      case 'P':
        if(lc->explicit_H & 1)
          break;
        else{
          for(unsigned int h=0;h<lc->explicit_H;h++){
            buffer += " ";
            write_locant(ch, buffer); 
            buffer += 'H';
          }
        }
        break;

      case 'S':
        if( !(lc->explicit_H & 1))
          break;   
        else{   
          for(unsigned int h=0;h<lc->explicit_H;h++){
            buffer += " ";
            write_locant(ch, buffer); 
            buffer += 'H';
          }
        }
        break;

      default:
        for(unsigned int h=0;h<lc->explicit_H;h++){
          buffer += " ";
          write_locant(ch, buffer); 
          buffer += 'H';
        }
        break;
      }
  }
  
  // spiro can then be sorted with a ignore clause on the locant atom, as the first ring instance can handle
  // the R groups. 

    // From the arranged locants, arrange the bonds
  for(unsigned char ch = 'A';ch < 'A'+ring->rsize;ch++) // locant sorting can be done after spiro
  {
    
    WLNSymbol *locant = ring->locants[ch]; 
    if(locant->spiro && locant == ignore){
      continue;
    }

    SortedEdges *ring_se = ArrangeRingBonds(locant,ring,seen_locants ,ignore);
    for(unsigned int i=0;i<ring_se->e_max;i++){
      WLNEdge *e = ring_se->edges[i];
      if(e == ring->macro_return || e->reverse == ring->macro_return){
        continue; 
      }

      if(!e->child->inRing){
        buffer += ' '; 
        write_locant(ch, buffer); 
        
        for(unsigned int i=1;i<e->order;i++){
          buffer += 'U';
        }

        CanonicalWLNChain(e->child, graph,locant,buffer);
      }
      else if(locant->spiro && e->child->inRing->locants_ch[locant]){
        buffer += ' ';
        write_locant(ch, buffer); 
        buffer += "-&"; 
        buffer += ' '; 
        buffer += e->child->inRing->locants_ch[locant];

        CanonicalWLNRing(e->child, graph,locant,buffer); // ignore where we've come from
      }
      else if(e->child->inRing != node->inRing && e->child->inRing->str_notation == "L6J"){
        buffer += ' '; 
        write_locant(ch, buffer); 
        
        for(unsigned int i=1;i<e->order;i++)
          buffer += 'U';

        // incoming == 'A' locant need to be rearranged          
        
        // locants are rotated from the incoming position, C = A, therefore D = B etc
        unsigned char incoming_char = e->child->inRing->locants_ch[e->child]; 
        
        std::vector<std::pair<WLNSymbol*,unsigned char>> new_positions; 
        if(incoming_char != 'A'){
        // wrap the locants around 
          WLNRing * benzene =e->child->inRing;
          for (unsigned char ch = 'A'; ch <= 'F'; ch++){
            if(ch < incoming_char){
              unsigned char new_loc = (ch-incoming_char) + 'F'+1;   
              new_positions.push_back({benzene->locants[ch],new_loc}); 
            }
            else{
              unsigned char new_loc = (ch-incoming_char) + 'A';   
              new_positions.push_back({benzene->locants[ch],new_loc}); 
            }
          }
          
          // seems a bit like bullying, but does work
          benzene->locants.clear();
          benzene->locants_ch.clear(); 
          for (std::pair<WLNSymbol*,unsigned char> p : new_positions){
            benzene->locants[p.second] = p.first;
            benzene->locants_ch[p.first] = p.second; 
          }
        }

        CanonicalWLNRing(e->child, graph, e->parent, buffer);
      }
      else if (e->child->inRing != node->inRing){
        buffer += ' '; 
        write_locant(ch, buffer); 
        
        for(unsigned int i=1;i<e->order;i++)
          buffer += 'U';
      
        buffer += '-';
        buffer += ' ';
        write_locant(e->child->inRing->locants_ch[e->child], buffer); 

        if(e != e->child->inRing->macro_return && e->reverse != e->child->inRing->macro_return){
          CanonicalWLNRing(e->child, graph,locant,buffer); // ignore where we've come from
        }
        else{
          buffer += "-x-J"; 
        }
      }
    }

    free(ring_se); 
  }
  
  buffer += '&'; // better logic following writer
  return true;
}

void WritePostCharges(WLNGraph &wln_graph, std::string &buffer){
  // handle post charges, no need to check ring here, solid function
  for(unsigned int i=0;i<wln_graph.symbol_count;i++){
    WLNSymbol *pos = wln_graph.SYMBOLS[i]; 
    if(pos->charge > 0 && pos->ch != 'K'){  
      for(unsigned int i=0;i<abs(pos->charge);i++){
        buffer += " &";
        buffer += std::to_string(pos->str_position);
        buffer += "/0"; 
      }
    }
    else if (pos->charge < 0 && !(pos->charge == -1 && pos->inRing && pos->ch == 'C') && buffer[pos->str_position-1] != 'W'){
      for(unsigned int i=0;i<abs(pos->charge);i++){
        buffer += " &0/";
        buffer += std::to_string(pos->str_position); 
      }
    }
  }
}

bool ChainOnlyCanonicalise(WLNGraph &wln_graph, std::set<WLNSymbol*> &whole_set,std::string &store){
  bool ion_write = false;
  

  for (unsigned int i=0;i<wln_graph.symbol_count;i++){
    WLNSymbol *node = wln_graph.SYMBOLS[i];
    if( (!node->barr_n || !node->parr_n) && !node->inRing && !(whole_set.count(node)) ){

      // create a local set, iterate through the set
      if(ion_write){ // ion condition
        while(store.back() == '&') // trail cleaning
          store.pop_back(); 
        store += " &";
      }
      
      std::set<WLNSymbol*> symbol_set; 
      Reachable(node, symbol_set); 
      
      std::string last_chain;
      WLNSymbol *best_start_point = 0; 

      for(std::set<WLNSymbol*>::iterator siter = symbol_set.begin();siter != symbol_set.end(); siter++){
        if( (!(*siter)->barr_n || !(*siter)->parr_n) && !(*siter)->inRing){
          std::string new_chain; 
          CanonicalWLNChain(*siter, wln_graph, 0,new_chain); 
          while(new_chain.back() == '&') // trail cleaning
            new_chain.pop_back(); 
          
          // i think charges must be preserved, so writing here is a must

          if(new_chain.size() < last_chain.size() || last_chain.empty()){
            last_chain = new_chain;
            best_start_point = (*siter); 
          }
          else if(new_chain.size() == last_chain.size()){ // take the highest ascii character, actually follows rule 2 order
            for(unsigned int j=0;j<new_chain.size();j++){
              if(new_chain[j] > last_chain[j]){
                last_chain = new_chain;
                best_start_point = (*siter); 
                break;
              }
              else if(new_chain[j] < last_chain[j])
                break;
            }
          }
        }
      }

      // write out the final string, this sets all the charge positions to where they should be
      CanonicalWLNChain(best_start_point, wln_graph, 0,store); 
      whole_set.insert(symbol_set.begin(),symbol_set.end()); 
      ion_write = true;
    }
  }
  
 while(store.back() == '&') // trail cleaning
   store.pop_back(); 

  // tidy up any remaining closures that do not need to be added
  return true; 
}


std::string FullCanonicalise(WLNGraph &graph){

  std::string last_chain;
  std::string store; 
  bool first_write = false;

  unsigned int r = 0; 
  unsigned int b = 0; 
  WLNRing *sorted_rings[128] = {0};
  WLNRing *benzyl[128] = {0}; 
  for (unsigned int i=0;i<graph.ring_count;i++){
    if(graph.RINGS[i]->str_notation != "L6J"){
      sorted_rings[r++] = graph.RINGS[i];
      graph.RINGS[i]->ranking = r; 
    }
    else {
      if(graph.RINGS[i]->loc_count <=1){
        for(unsigned char loc='A';loc<='F';loc++){
          if(graph.RINGS[i]->locants[loc]->num_edges < 4){
            WLNSymbol *eHydrogen = AllocateWLNSymbol('H', graph);
            eHydrogen->allowed_edges = 1; 
            AddEdge(eHydrogen, graph.RINGS[i]->locants[loc]);
            break;
          }
        }
      }
      benzyl[b++] = graph.RINGS[i];
    }
  }
  
  // sort by how many locants, inverse order will likely give the shorter string
  for (unsigned int j=1;j<r;j++){
    WLNRing *s = sorted_rings[j];
    unsigned int key = s->loc_count; 
		int i = j-1;
    while(i>=0 && sorted_rings[i]->loc_count >= key){
      sorted_rings[i+1] = sorted_rings[i];
      i--;
    }
		sorted_rings[i+1] = s;
	}

  // if a ring has a macrocycle present, place at the front
  for (unsigned int j=1;j<r;j++){
    WLNRing *s = sorted_rings[j];
    unsigned int key = s->macro_return ? 1:0; 
		int i = j-1;
    while(i>=0 && (sorted_rings[i]->macro_return ? 1:0) <= key){
      sorted_rings[i+1] = sorted_rings[i];
      i--;
    }
		sorted_rings[i+1] = s;
	}


#if OPT_DEBUG
  for (unsigned int i=0;i<r;i++){
    fprintf(stderr,"  %d: ring-size:%d, locants: %d, multi-points: %d, bridges: %d\n",
        sorted_rings[i]->ranking,
        sorted_rings[i]->rsize,
        sorted_rings[i]->loc_count,
        sorted_rings[i]->multi_points,
        sorted_rings[i]->bridge_points
        ); 
  }

#endif

  std::map<WLNRing*, bool> seen_rings;
  std::set<WLNSymbol*> all_symbols; 
  for (unsigned int i=0;i<r;i++){
    WLNRing *ring = sorted_rings[i];
    WLNSymbol *node = ring->locants['A']; 
    if(!seen_rings[ring]){

      std::set<WLNSymbol*> reachable; 
      Reachable(node, reachable); 
      all_symbols.insert(reachable.begin(),reachable.end()); 

      for(std::set<WLNSymbol*>::iterator siter = reachable.begin();siter!=reachable.end();siter++){
        if( (*siter)->inRing )
          seen_rings[ (*siter)->inRing ] = true;
      }

      if(first_write){
        while(store.back() == '&')
          store.pop_back(); 

        store += " &";
      }

      CanonicalWLNRing(node, graph,0,store);
      first_write = true;
    }
  }

  if(first_write)
    store += " &"; 
  ChainOnlyCanonicalise(graph,all_symbols,store); 
  while(store.back() == '&' || store.back() == ' ')
    store.pop_back();
  return store; 
}

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
    wln_input = ptr; 

  unsigned int len = strlen(wln_input);

  WLNGraph wln_graph;
  BabelGraph obabel; 

  if(!ParseWLNString(ptr,wln_graph))
    return false;

    // needs to be this order to allow K to take the methyl groups
  if(!WLNKekulize(wln_graph))
    return Fatal(len,"Error: failed to kekulize mol");

  if(!ExpandWLNSymbols(wln_graph,len))
    return false;

  if(!obabel.ConvertFromWLN(mol,wln_graph,len))
    return false;

  obabel.NMOBSanitizeMol(mol);
  return true;
}




bool CanonicaliseWLN(const char *ptr, OBMol* mol)
{   
  if(!ptr){
    fprintf(stderr,"Error: could not read wln string pointer\n");
    return false;
  }
  else 
    wln_input = ptr; 

  WLNGraph wln_graph;
  BabelGraph obabel; 

  if(!ParseWLNString(ptr,wln_graph))
    return false;
  
  if(!WLNKekulize(wln_graph))
    return false; 

  // more minimal resolve step for certain groups, W removal
  unsigned int stop = wln_graph.symbol_count;
  stop = wln_graph.symbol_count;
  for (unsigned int i=0;i<stop;i++){
    WLNSymbol *sym = wln_graph.SYMBOLS[i];
    switch(sym->ch){
      case 'Y':
      case 'X':
      case 'K':
        if(!resolve_methyls(sym,wln_graph))
          return false;
        break;

      case 'W':
        if(sym->barr_n){
          sym->bond_array[0].order = 1;
          sym->bond_array[0].reverse->order = 1; 
        }
        if(sym->parr_n){
          sym->prev_array[0].order = 1;
          sym->prev_array[0].reverse->order = 1; 
        }
        break;
    }
  }
  
  std::string res; 
  // if no rings, choose a starting atom and flow from each, ions must be handled seperately
  if(!wln_graph.ring_count){
    std::set<WLNSymbol*> seen_set;
    ChainOnlyCanonicalise(wln_graph,seen_set,res); // bit more effecient 
  }
  else
    res = FullCanonicalise(wln_graph); 
  
  WritePostCharges(wln_graph, res); 
  std::cout << res << std::endl; 
  return true;
}

