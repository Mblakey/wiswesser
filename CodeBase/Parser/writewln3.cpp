

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector>
#include <stack>
#include <map>
#include <deque>
#include <iterator>
#include <sstream>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>


// --- macros ---
#define REASONABLE 1024

// --- inputs ---
const char *wln;
const char *dotfile;

// --- options ---
static bool opt_wln2dot = false;
static bool opt_allow = false;
static bool opt_debug = false;
static bool opt_convert = false;


// --- globals ---
struct WLNSymbol;
struct WLNRing; 

std::map<WLNSymbol *, unsigned int> index_lookup;
std::map<unsigned int, WLNSymbol*> symbol_lookup;
unsigned int glob_index = 0;


// --- pools --- 
std::vector<WLNSymbol*> symbol_mempool;
std::vector<WLNRing*>   ring_mempool;


enum WLNTYPE{ STANDARD = 0, LOCANT = 1, LINKER = 2, RING = 3, SPECIAL = 4};


// rule 2 - hierarchy - rules have diverged due to end terminator char, also use for locant setting from 14
std::map<unsigned char, unsigned int> char_hierarchy =
{
  {' ', 1}, {'-', 2}, {'/', 3}, 
  {'0', 4}, {'1', 5}, {'2', 6}, {'3', 7}, {'4', 8}, {'5', 9}, {'6', 10}, {'7', 11}, {'8', 12}, {'9', 13}, 
  {'A', 14}, {'B', 15}, {'C', 16}, {'D', 17}, {'E', 18}, {'F', 19}, {'G', 20}, {'H', 21}, {'I', 22}, {'J', 23},
  {'K', 24}, {'L', 25}, {'M', 26}, {'N', 27}, {'O', 28}, {'P', 29}, {'Q', 30}, {'R', 31}, {'S', 32}, {'T', 33}, 
  {'U', 34}, {'V', 35}, {'W', 36}, {'X', 37}, {'Y', 38}, {'Z', 40}, {'&', 41}
};


std::map<unsigned char,unsigned int> locant_symbols =
{
  {'A',1}, {'B',2}, {'C',3}, {'D',4}, {'E',5}, {'F',6}, {'G',7}, 
  {'H',8}, {'I',9}, {'J',10}, {'K',11}, {'L',12}, {'M',13}, 
  {'N',14}, {'O',15}, {'P',16}, {'Q',17}, {'R',18}, {'S',19}, 
  {'T',20}, {'U',21}, {'V',22}, {'W',23}, {'X',25}, {'Y',26}, {'Z',27}
};



// --- utilities ---

static bool isdigit_str(const std::string& s)
{
  for (char const &ch : s) {
    if (std::isdigit(ch) == 0) 
      return false;
  }
  return true;
}

static void get_notation(unsigned int s, unsigned int e, std::string &res){
  for (unsigned int i=s; i<=e; i++){
    res.push_back(wln[i]);
  }
}


static void Fatal(unsigned int pos){
  fprintf(stderr,"Fatal: %s\n",wln);
  fprintf(stderr,"       ");
  
  for(unsigned int i=0; i<pos;i++)
    fprintf(stderr," ");

  fprintf(stderr,"^\n");

  exit(1);
}

static void Reindex_lookups(){
  glob_index = 0; 
  for (WLNSymbol *node : symbol_mempool){
    index_lookup[node] = glob_index;
    symbol_lookup[glob_index] = node;
    glob_index++; 
  }
}



// should be all we need for a SCT XI connection table - obabel can handle coords
struct Atom{
  std::string symbol; 
  unsigned int atomic_num;

  int charge; 

  std::vector<Atom> bonded;
  std::vector<unsigned int> orders; 
};


struct AtomGraph{
  Atom *head; 
}; 





struct WLNSymbol
{

  unsigned char ch;
  unsigned int type; 

  unsigned int allowed_edges;
  unsigned int num_edges;

  WLNSymbol *previous;
  std::vector<WLNSymbol*>     children; // linked list of next terms chains
  std::vector<unsigned int>   orders;
  

  //if 'ch='*'
  std::string special; // string for element, or ring
  WLNRing* ring;

  // if default needed
  WLNSymbol()
  {
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
    
    previous = 0; 

    special = "";
    ring = 0; 
  }
  ~WLNSymbol(){};


  void set_edges(unsigned int edges){
    allowed_edges = edges; 
  }

  void set_type(unsigned int i){
    type = i; 
  }
  
  // resets the symbol 
  void reset(){
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
  }

  void add_special(unsigned int s, unsigned int e){
    for (unsigned int i=s; i<= e;i++)
      special.push_back(wln[i]);
  }

  void add_special(std::string str){
    special.append(str); 
  }

};

WLNSymbol *AllocateWLNSymbol(unsigned char ch){
  
  if(opt_debug)
    fprintf(stderr,"  allocating %c\n",ch);
  
  WLNSymbol *wln = new WLNSymbol;
  symbol_mempool.push_back(wln);
  wln->ch = ch;

  index_lookup[wln] = glob_index;
  symbol_lookup[glob_index] = wln;
  glob_index++; 
  return wln;
}

// these are expensive, but needed for some edge case notations
void DeallocateWLNSymbol(WLNSymbol *node){

  if(opt_debug)
    fprintf(stderr,"  manual deallocation: %c\n",node->ch);
  
  // find the node in the mem pool 
  unsigned int i = 0; 
  for (WLNSymbol *n : symbol_mempool){
    if (n == node)
      break;
    i++;
  }

  symbol_mempool.erase(symbol_mempool.begin() + i);
  delete node;
}


WLNSymbol* copy_symbol(WLNSymbol *src){
  
  WLNSymbol *copy = AllocateWLNSymbol(src->ch);
  copy->allowed_edges = src->allowed_edges;
  copy->num_edges = src->num_edges;
  
  for (unsigned int i=0; i<src->children.size();i++){
    copy->children.push_back(src->children[i]);
    copy->orders.push_back(src->orders[i]);
  }
  return copy;
}


/* struct to hold pointers for the wln ring */
struct WLNRing
{

  unsigned int size;
  bool aromatic;
  bool heterocyclic;

  std::vector<unsigned int>  ring_components; 
  std::map<unsigned char,unsigned char> fuse_points; // gives the fusing points for combined rings
  

  // nope this is actually very clever, use 1's for carbons
  // everything else gets a wln symbol, then converts to atoms


  std::map<WLNSymbol*,unsigned char> locants; 

  WLNRing()
  {
    size = 0;
    aromatic = false;
    heterocyclic = false;
  }
  ~WLNRing(){};


  unsigned int consume_ring_notation(std::string &block)
  {

    unsigned int local_size = 0; 
    unsigned int len = block.size() - 1;

    if (block.size() < 3){
      fprintf(stderr,"Error: not enough chars to build ring - %s\n",block.c_str());
      return false; 
    }

    if (block[0] == 'T')
      heterocyclic = true;
    else if (block[0] == 'L')
      heterocyclic = false;
    else{
      fprintf(stderr,"Error: first character in ring notation must be an L|T\n");
      return 0; 
    }

    if (block[len] != 'J'){
      fprintf(stderr,"Error: last character in ring notation must be J\n");
      return 0; 
    }

    if (block[len-1] == 'T'){
      aromatic = false;
    }
    else
      aromatic = true; 


    if (block[1] == ' '){
      // special ring types

    }
    else{
      // check how many locants are allowed 
      unsigned int it = 1;
      unsigned int rings = 0; 
      while (std::isdigit(block[it]) && block[it] != '\0')
      {
        unsigned int val = block[it] - '0';
        ring_components.push_back(val);

        local_size += val;
        rings++;
        it++;
      }

      if(rings > 1){
        // refactor size down 
        unsigned int term = rings - 2;
        unsigned int shared_atoms = rings + term;
        local_size = local_size - shared_atoms;
      }
      
      if (opt_debug)
        fprintf(stderr,"  evaluated ring to size %d\n",local_size);

      // create the pseudo ring

      /* process the substring*/
      process_interconnections(block.substr(it,block.length()));
      
    }
    
    return local_size; // returns size as a value
  }

  unsigned int process_interconnections(std::string block){

    /* 
    we assume this start after the ring blocks 
    i.e 'L66 AO TJ' -->' AO TJ'
    */ 

    // split the string on the spaces, and then process the blocks
    std::istringstream ss(block);
    std::string del;

    unsigned char locant = 'A'; // use the maps to assign locants
    WLNSymbol *assignment = 0; 

    while(getline(ss, del, ' ')){
     

      // process the locants as expected
      if(del.back() != 'J'){
        locant = del.front();
        //assignment = positional_symbols[locant];
      }

     
    }

    return 0;
  }



};


WLNRing *AllocateWLNRing(){
  WLNRing *wln_ring = new WLNRing;
  ring_mempool.push_back(wln_ring);
  return wln_ring;
}


// expensive but sometimes necessary for edge cases
void DeallocateWLNRing(WLNRing *ring){
  // find the ring in the mem pool 
  unsigned int i = 0; 
  for (WLNRing *r : ring_mempool){
    if (r == ring)
      break;
    i++;
  }

  ring_mempool.erase(ring_mempool.begin() + i);
  delete ring;
}



struct WLNGraph{

  WLNSymbol *root;

  WLNGraph() : root{(WLNSymbol *)0}{};
  ~WLNGraph(){
    for (WLNSymbol * node : symbol_mempool)
      delete node; 
    for (WLNRing * ring : ring_mempool)
      delete ring; 
  };


#ifdef DEV
  /* handles all inter ring defintions*/
  bool ParseInterRing(unsigned int start, unsigned int end, WLNRing *ring)
  {

    // locants are sequential if inline defined e.g AUO places O on B
    // start here should be where the cyclic values are ended
    // < so should not hit J

    bool pending_special = false;
    bool pending_locant = false;

    std::string special; 

    WLNSymbol *atom = 0;
    unsigned char cur_locant = '\0';

    for (unsigned int i = start; i < end; i++)
    {
      unsigned char ch = wln[i];
      switch (ch)
      {

      case 'A':
      case 'C':
      case 'D':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
      case 'J':
      case 'L':
      case 'Q':
      case 'R':
      case 'T':
      case 'X':
      case 'Y':
      case 'Z':
        if (pending_locant)
        {
          cur_locant = ch;
          atom = access_locant(ch, ring);
          if (!atom)
            return false;

          pending_locant = false;
        }
        else
        {
          fprintf(stderr, "Error: invalid definition in inter ring notation\n");
          return false;
        }
        break;

      case 'B':
      case 'K':
      case 'M':
      case 'N':
      case 'O':
      case 'P':
      case 'S':
        if (pending_locant)
        {
          cur_locant = ch;
          atom = access_locant(ch, ring);
          if (!atom)
            return false;

          pending_locant = false;
        }
        else
        {
          transform_symbol(atom, ch);
          atom = access_locant(cur_locant + 1, ring,false);
        }

        break;

      case 'U':
        if (pending_locant)
        {
          cur_locant = ch;
          atom = access_locant(ch, ring);
          if (!atom)
            return false;

          pending_locant = false;
        }
        else{
          atom = access_locant(cur_locant + 1, ring,false);
        }
        break;

      case 'V':
        if (pending_locant)
        {
          cur_locant = ch;
          atom = access_locant(ch, ring);
          if (!atom)
            return false;

          pending_locant = false;
        }
        else
        {
          WLNSymbol *oxy = AllocateWLNSymbol('O');
          if(!link_symbols(oxy,atom,1))
            return false;
          
          atom = access_locant(cur_locant + 1, ring,false);
        }

        break;

      case 'W':
        if (pending_locant)
        {
          cur_locant = ch;
          atom = access_locant(ch, ring);
          if (!atom)
            return false;

          pending_locant = false;
        }
        else
        {
          WLNSymbol *oxy_1 = AllocateWLNSymbol('O');
          WLNSymbol *oxy_2 = AllocateWLNSymbol('O');
          if(!link_symbols(oxy_1,atom,1) ||  !link_symbols(oxy_2,atom,1))
            return false;
          atom = access_locant(cur_locant + 1, ring,false);
        }
        break;

      case ' ':
        pending_locant = true;
        break;

      case '-': // allows inter-ring specific atoms

      default:
        fprintf(stderr, "Error: invalid symbol in inter ring notation - %c\n",ch);
        return false;
      }
    }

    return true;
  }





  /* one char so should be notation independent */
  WLNRing *consume_benzene(){
    
    WLNRing *ring = AllocateWLNRing();

    // 1) create a big ring
    WLNSymbol *rhead = AllocateWLNSymbol('C');
    ring->rhead = rhead; // set the rings head

    WLNSymbol *current = 0;
    WLNSymbol *prev = rhead;

    unsigned int locant = 0;
    for (unsigned int i = 1; i < 6; i++)
    {
      current = AllocateWLNSymbol('C');
      ring->locants[locant_symbols[locant++]] = current; // add the locants  
      add_aromatic(current,prev);
      prev = current;
    }
    add_aromatic(rhead,current);
    return ring; 
  }

  /* platform for launching ring notation build functions */
  WLNRing *consume_ring_notation(unsigned int start, unsigned int end)
  {

    bool handle_advanced = false;

    // 1) allocate the blank ring object
    WLNRing *wln_ring = AllocateWLNRing();

    // 2) minimum symbols for a ring notation is 3 - allows safe lookback

    if ((end - start) < 2)
    {
      fprintf(stderr, "Error: minimum chars for ring notation is 3 - found: %d\n", end - start);
      return (WLNRing *)0;
    }

    // 3) evaluate start character
    switch (wln[start])
    {

    case 'L':
      wln_ring->heterocyclic = false;
      break;

    case 'T':
      wln_ring->heterocyclic = true;
      break;

    default:
      fprintf(stderr, "Error: ring notation must start L|T ... not: %c\n", wln[start]);
      return (WLNRing *)0;
    }

    // 3) advanced vs standard notation test on second char

    switch (wln[start + 1])
    {

    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      break;

    case ' ':
      handle_advanced = true;
      break;

    default:
      fprintf(stderr, "Error: unknown second char in ring notation: %c\n", wln[start + 1]);
      return (WLNRing *)0;
    }

    // 4) aromatic testing on the second to last 'T'

    switch (wln[end - 1])
    {
    case 'T':
      wln_ring->aromatic = false;

      // we can move the end back in this case 'helps with locants'
      end += -1;
      break;

    default:
      wln_ring->aromatic = true;
    }

    if (handle_advanced)
    {
      // create the poly cyclic functions here
    }
    else
      if(!CreateStandardRing(start, end, wln_ring))
        return (WLNRing *)0; 

    return wln_ring;
  }
  
#endif

  /* should handle all bonding modes, adds child to parent->children
  'UU' bonding also added here */
  bool link_symbols(WLNSymbol *child, WLNSymbol *parent, unsigned int bond)
  {

    if (parent->ch == '*' && parent->ring){
      fprintf(stderr,"Error: trying to link a ring through standard notation, locants needed\n");
      return false; 
    }
   
    // if the child cannot handle the new valence
    if ( (child->num_edges + bond) > child->allowed_edges ){
      fprintf(stderr,"Error: wln character[%c] is exceeding allowed connections\n",child->ch);
      return false; 
    }

    // same for the parent
    if ( (parent->num_edges + bond) > parent->allowed_edges ){
      fprintf(stderr,"Error: wln character[%c] is exceeding allowed connections\n",parent->ch);
      return false; 
    }


    child->previous = parent; // keep the linked list so i can consume backwards rings

    child->num_edges  += bond;
    parent->num_edges += bond;

    parent->children.push_back(child);
    parent->orders.push_back(bond);

    return true;
  }


  WLNSymbol* define_element(std::vector<unsigned char> &special){
    
    // allocate a special wln
    WLNSymbol *created_wln = AllocateWLNSymbol('*');

    // some fancy switching

    switch(special[0]){

      case 'A':
        if(special[1] == 'C')
          created_wln->special = "Ac";
        else if (special[1] == 'G')
          created_wln->special = "Ag";
        else if (special[1] == 'L')
          created_wln->special = "Al";
        else if (special[1] == 'M')
          created_wln->special = "Am";
        else if (special[1] == 'R')
          created_wln->special = "Ar";
        else if (special[1] == 'S')
          created_wln->special = "As";
        else if (special[1] == 'T')
          created_wln->special = "At";
        else if (special[1] == 'U')
          created_wln->special = "Au";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'B':
        if(special[1] == 'A')
          created_wln->special = "Ba";
        else if (special[1] == 'E')
          created_wln->special = "Be";
        else if (special[1] == 'H')
          created_wln->special = "Bh";
        else if (special[1] == 'I')
          created_wln->special = "Bi";
        else if (special[1] == 'K')
          created_wln->special = "Bk";
        else if (special[1] == 'R')
          created_wln->special = "Br";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'C':
        if(special[1] == 'A')
          created_wln->special = "Ca";
        else if (special[1] == 'D')
          created_wln->special = "Cd";
        else if (special[1] == 'E')
          created_wln->special = "Ce";
        else if (special[1] == 'F')
          created_wln->special = "Cf";
        else if (special[1] == 'M')
          created_wln->special = "Cm";
        else if (special[1] == 'N')
          created_wln->special = "Cn";
        else if (special[1] == 'O')
          created_wln->special = "Co";
        else if (special[1] == 'R')
          created_wln->special = "Cr";
        else if (special[1] == 'S')
          created_wln->special = "Cs";
        else if (special[1] == 'U')
          created_wln->special = "Cu";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'D':
        if(special[1] == 'B')
          created_wln->special = "Db";
        else if (special[1] == 'S')
          created_wln->special = "Ds";
        else if (special[1] == 'Y')
          created_wln->special = "Dy";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'E':
        if(special[1] == 'R')
          created_wln->special = "Er";
        else if (special[1] == 'S')
          created_wln->special = "Es";
        else if (special[1] == 'U')
          created_wln->special = "Eu";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'F':
        if(special[1] == 'E')
          created_wln->special = "Fe";
        else if (special[1] == 'L')
          created_wln->special = "Fl";
        else if (special[1] == 'M')
          created_wln->special = "Fm";
        else if (special[1] == 'R')
          created_wln->special = "Fr";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'G':
        if(special[1] == 'A')
          created_wln->special = "Ga";
        else if (special[1] == 'D')
          created_wln->special = "Gd";
        else if (special[1] == 'E')
          created_wln->special = "Ge";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'H':
        if(special[1] == 'E')
          created_wln->special = "Ha";
        else if (special[1] == 'F')
          created_wln->special = "Hf";
        else if (special[1] == 'G')
          created_wln->special = "Hg";
        else if (special[1] == 'O')
          created_wln->special = "Ho";
        else if (special[1] == 'S')
          created_wln->special = "Hs";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;


      case 'I':
        if(special[1] == 'N')
          created_wln->special = "In";
        else if (special[1] == 'R')
          created_wln->special = "Ir";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'K':
        if(special[1] == 'R')
          created_wln->special = "Kr";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'L':
        if(special[1] == 'A')
          created_wln->special = "La";
        else if (special[1] == 'I')
          created_wln->special = "Li";
        else if (special[1] == 'R')
          created_wln->special = "Lr";
        else if (special[1] == 'U')
          created_wln->special = "Lu";
        else if (special[1] == 'V')
          created_wln->special = "Lv";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'M':
        if(special[1] == 'C')
          created_wln->special = "Mc";
        else if (special[1] == 'D')
          created_wln->special = "Md";
        else if (special[1] == 'G')
          created_wln->special = "Mg";
        else if (special[1] == 'N')
          created_wln->special = "Mn";
        else if (special[1] == 'O')
          created_wln->special = "Mo";
        else if (special[1] == 'T')
          created_wln->special = "Mt";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;


      case 'N':
        if(special[1] == 'A')
          created_wln->special = "Na";
        else if (special[1] == 'B')
          created_wln->special = "Nb";
        else if (special[1] == 'D')
          created_wln->special = "Nd";
        else if (special[1] == 'E')
          created_wln->special = "Ne";
        else if (special[1] == 'H')
          created_wln->special = "Nh";
        else if (special[1] == 'I')
          created_wln->special = "Ni";
        else if (special[1] == 'O')
          created_wln->special = "No";
        else if (special[1] == 'P')
          created_wln->special = "Np";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;


      case 'O':
        if(special[1] == 'G')
          created_wln->special = "Og";
        else if(special[1] == 'S')
          created_wln->special = "Os";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'P':
        if(special[1] == 'A')
          created_wln->special = "Pa";
        else if(special[1] == 'B')
          created_wln->special = "Pb";
        else if(special[1] == 'D')
          created_wln->special = "Pd";
        else if(special[1] == 'M')
          created_wln->special = "Pm";
        else if(special[1] == 'O')
          created_wln->special = "Po";
        else if(special[1] == 'R')
          created_wln->special = "Pr";
        else if(special[1] == 'T')
          created_wln->special = "Pt";
        else if(special[1] == 'U')
          created_wln->special = "Pu";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;
      
      case 'R':
        if(special[1] == 'A')
          created_wln->special = "Ra";
        else if(special[1] == 'B')
          created_wln->special = "Rb";
        else if(special[1] == 'E')
          created_wln->special = "Re";
        else if(special[1] == 'F')
          created_wln->special = "Rf";
        else if(special[1] == 'G')
          created_wln->special = "Rg";
        else if(special[1] == 'H')
          created_wln->special = "Rh";
        else if(special[1] == 'N')
          created_wln->special = "Rn";
        else if(special[1] == 'U')
          created_wln->special = "Ru";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;
      
      case 'S':
        if(special[1] == 'B')
          created_wln->special = "Sb";
        else if(special[1] == 'C')
          created_wln->special = "Sc";
        else if(special[1] == 'E')
          created_wln->special = "Se";
        else if(special[1] == 'I')
          created_wln->special = "Si";
        else if(special[1] == 'M')
          created_wln->special = "Sm";
        else if(special[1] == 'N')
          created_wln->special = "Sn";
        else if(special[1] == 'R')
          created_wln->special = "Sr";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;


      case 'T':
        if(special[1] == 'A')
          created_wln->special = "Ta";
        else if(special[1] == 'B')
          created_wln->special = "Tb";
        else if(special[1] == 'C')
          created_wln->special = "Tc";
        else if(special[1] == 'E')
          created_wln->special = "Te";
        else if(special[1] == 'H')
          created_wln->special = "Th";
        else if(special[1] == 'I')
          created_wln->special = "Ti";
        else if(special[1] == 'L')
          created_wln->special = "Tl";
        else if(special[1] == 'M')
          created_wln->special = "Tm";
        else if(special[1] == 'S')
          created_wln->special = "Ts";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;


      case 'X':
        if(special[1] == 'E')
          created_wln->special = "Xe";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'Y':
        if(special[1] == 'B')
          created_wln->special = "Yb";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      case 'Z':
        if(special[1] == 'N')
          created_wln->special = "Zn";
        else if (special[1] == 'R')
          created_wln->special = "Zr";
        else{
          fprintf(stderr,"Error: invalid element symbol in special definition\n");
          return (WLNSymbol*)0; 
        }
      break;

      default:
        fprintf(stderr,"Error: invalid character in special definition switch\n");
        return (WLNSymbol*)0; 
    }

    created_wln->allowed_edges = 8; // allow on octet default for these species. 

    return created_wln; 
  }


  WLNSymbol* return_open_branch(std::stack<WLNSymbol*> &branch_stack){

    if (branch_stack.empty())
      return (WLNSymbol*)0;

    WLNSymbol *top = branch_stack.top();

    // doesnt pop here, only '&' can pop the stack

    return top;
  }

  bool check_unbroken(unsigned int i){
    if(i > 1 && !(wln[i-1] == '&' && wln[i-2] == ' ')){
      fprintf(stderr,"Error: broken graph without ionic notation, check branches|locants and '&' count\n");
      return false;
    }
    
    return true;
  }

  WLNSymbol* pop_ringstack(unsigned int pops, std::stack <WLNSymbol*> &stack){

    if (pops >= stack.size()){
      fprintf(stderr,"Error: trying to pop too many rings check '&' count\n");
      return (WLNSymbol*)0; 
    }
    
    for (unsigned int i=0; i<pops;i++)
      stack.pop();
    
    return stack.top();
  }

  // this has a return clause in it and needs previous
  WLNSymbol* pop_branchstack(unsigned int pops, std::stack <WLNSymbol*> &stack, WLNSymbol *prev){

    if(!prev)
      fprintf(stderr,"Error: popping with no previous symbol\n");

    bool hard = false; 

    if(prev == stack.top())
      hard = true;

    if(opt_debug)
      fprintf(stderr,"  popping %d symbols down the stack: mode(%d) prev[%c]\n",pops, hard,prev->ch);

    if(hard){
      if(pops >= stack.size()){
        fprintf(stderr,"Error: to many stack pops - check '&' count\n");
        return 0; 
      }
      for (unsigned int i=0; i<pops;i++)
        stack.pop();      
    }
    else{
      if(pops > stack.size()){
        fprintf(stderr,"Error: to many stack pops - check '&' count\n");
        return 0; 
      }
      for (unsigned int i=1; i<pops;i++)
        stack.pop();      
    }
    return stack.top();
  }

  /* wraps popping for the linker and branch stacks */
  WLNSymbol *pop_standard_stacks(unsigned int pop_ticks, 
                        std::stack<WLNSymbol*> &branch_stack,
                        std::stack<WLNSymbol*> &linker_stack,
                        WLNSymbol* prev, unsigned int i)
  {
    WLNSymbol *ret = 0; 
    if(!branch_stack.empty())
      ret = pop_branchstack(pop_ticks,branch_stack,prev);
    else if(!linker_stack.empty())
      ret = pop_branchstack(pop_ticks,linker_stack,prev);
    else{
      fprintf(stderr,"Error: popping empty stacks - check '&' count\n");
      Fatal(i);
    }
      
    return ret;   
  }


  /* wraps the linking and graph checking functions */
  void create_bond( WLNSymbol *curr, WLNSymbol *prev, 
                    unsigned int bond_ticks, unsigned int i)
  {
    if(prev){
      if(!link_symbols(curr,prev,1 + bond_ticks))
        Fatal(i);
    }
    else{
      if(!check_unbroken(i))
        Fatal(i);
    }
  }

  /*  Wraps the creation of locant and bonding back ring assignment */
  void create_locant(WLNSymbol *curr,std::stack<WLNSymbol*> &ring_stack,unsigned int i){

    WLNSymbol *s_ring = 0; 
    unsigned char ch = wln[i];

    if(ring_stack.empty()){
      fprintf(stderr,"Error: no rings to assign locants to\n");
      Fatal(i);
    } 
    else
      s_ring = ring_stack.top();

    if(locant_symbols[ch] < s_ring->ring->size)
      s_ring->children.push_back(curr);
    else{
      fprintf(stderr,"Error: assigning locant greater than ring size\n");
      Fatal(i);
    }
  }

  /* a global segmentation using both rule sets - start merging */
  bool ParseWLNString(const char *wln, unsigned int len){


    std::stack <WLNSymbol*>   ring_stack;   // access through symbol
    std::stack <WLNSymbol*>   branch_stack; // between locants, clean branch stack
    std::stack <WLNSymbol*>   linker_stack; // used for branching ring systems 

    WLNSymbol *curr = 0;
    WLNSymbol *prev = 0;  

    bool pending_locant         = false;
    bool pending_special        = false;
    bool pending_closure        = false; 
    bool pending_inline_ring    = false;
    bool pending_spiro          = false;

    // allows consumption of notation after block parses
    unsigned int block_start =  0;
    unsigned int block_end    = 0;

    unsigned int pop_ticks  = 0; // '&' style popping
    unsigned int bond_ticks = 0; // 'U' style bonding

    for (unsigned int i=0; i<len; i++){
      unsigned char ch = wln[i];
      
      if(opt_debug)
        fprintf(stderr,"Parsing: %c\n",ch);
    

      switch (ch){

        case '0': // cannot be lone, must be an addition to another num
          if(pending_closure || pending_special){
            break;
          }
          if (i == 0)
            Fatal(i);
          else if(i > 0 && !std::isdigit(wln[i-1]))
            Fatal(i);
          else{
            curr = AllocateWLNSymbol(ch);
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
          if(pop_ticks){
            prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
            pop_ticks = 0; 
          }

          if(pending_closure || pending_special){
            break;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(2);

          create_bond(curr,prev,bond_ticks,i);

          prev = curr;
          break;


        // carbons

        case 'Y':
          
          if(pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);
            

            prev = curr; 
            pending_locant = false;
          }
          else{

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(3);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;

        case 'X':
          if(pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else {

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(4);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;


        // oxygens 

        case 'O':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else {

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(2);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;

        case 'Q':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else if(pending_closure || pending_special){
            continue; 
          }
          else{
            
            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(1);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = return_open_branch(branch_stack);
          }
          break;


        case 'V':
        case 'W':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else {

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(2);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr; 
          }
          break;



        // nitrogens 

        case 'N':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else {

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(3);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;

        case 'M':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else {

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0; 
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(2);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;

        case 'K':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else {

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(4);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;


        case 'Z':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else if(pending_closure || pending_special){
            continue; 
          }
          else{

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(1);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = return_open_branch(branch_stack);
          }
          break;


        // halogens - need to add rules for semi allowed hyper valence in ionions

        case 'E':
        case 'G':
        case 'F':
        case 'I':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else if(pending_closure || pending_special){
            continue; 
          }
          else{

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0; 
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(1);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            
            prev = return_open_branch(branch_stack);
          }
          break;

        
        // inorganics 

        case 'B':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else if(pending_closure || pending_special){
            continue; 
          }
          else{

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(3);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;

  
        case 'P':
        case 'S':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else if(pending_closure || pending_special){
            continue; 
          }
          else{

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol(ch);
            curr->set_type(STANDARD);
            curr->set_edges(6);

            branch_stack.push(curr);

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0;
            prev = curr;
          }
          break;


        // locants only?

        case 'A':
        case 'C':
        case 'D':
        case 'H':
          if (pending_closure || pending_special){
            break;
          }
          else if (pending_locant){
            
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else
            Fatal(i);

          break;



        // ring notation

        case 'J':
          if (pending_special){
            break;
          }
          if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }

          else if(pending_closure){
            block_end = i;

            curr = AllocateWLNSymbol('*');
            curr->set_type(RING);
            curr->ring = AllocateWLNRing();

            curr->add_special(block_start,block_end);
            block_start = 0; 
            block_end = 0; 

            curr->ring->size = curr->ring->consume_ring_notation(curr->special);

            ring_stack.push(curr);

            if(pending_spiro){
              prev->type = LINKER; // spiros are normal rings with dual linker notation
              prev->previous->type = LINKER; 
              pending_spiro = false; 
            }

            // does the incoming locant check
            if(prev){
              prev->children.push_back(curr);
              if(locant_symbols[prev->ch] > curr->ring->size){
                fprintf(stderr,"Error: attaching inline ring with out of bounds locant assignment\n");
                Fatal(i);
              }
            }

            bond_ticks = 0;
            prev = curr; 
            pending_closure = false;
          }
          else
            Fatal(i);

          break;

        case 'L':
        case 'T':
          if(pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else{
            
            if(i==0)
              pending_inline_ring = true; 
           
            if(!pending_inline_ring){
              fprintf(stderr,"Error: ring notation started without '-' denotion\n");
              Fatal(i);
            }
            else 
              pending_inline_ring = false;

  
            block_start = i; 
            pending_closure = true; 
          }
          break;

        case 'R':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT);
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else{

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0;
            }

            curr = AllocateWLNSymbol('*');
            curr->set_type(RING);
            curr->ring = AllocateWLNRing();

            curr->add_special("L6J");
            curr->ring->size = curr->ring->consume_ring_notation(curr->special);

            ring_stack.push(curr);

            curr->set_edges(1); // for inline R's they cannot have more than 1; 

            create_bond(curr,prev,bond_ticks,i);
            
            bond_ticks = 0;
            prev = curr; 
          }
          break;



        // bonding

        case 'U':
          if (pending_closure || pending_special){
            break;
          }
          else if(pending_locant){
            curr = AllocateWLNSymbol(ch);
            curr->set_type(LOCANT); 
            curr->set_edges(2); // locants always have two edges
            
            if(pending_inline_ring)
              create_bond(curr,prev,bond_ticks,i);
            else
              create_locant(curr,ring_stack,i);

            prev = curr; 
            pending_locant = false;
          }
          else
            bond_ticks++;

          break;
          

        
        // specials 


        case ' ':
          if (pending_closure){
            break;
          }
          if (pending_special){
            pending_special = false; 
            block_start = 0; 
            block_end = 0; 
          }

          // clear the branch stack betweek locants and ions
          while(!branch_stack.empty())
            branch_stack.pop();

          if(pop_ticks){
            prev = pop_ringstack(pop_ticks,ring_stack);
            if(!prev)
              Fatal(i);
            pop_ticks = 0;
          } 
           
          pending_locant = true;
          break;
          
        case '&': 
          if (pending_closure || pending_special){
            break;
          }

          if(pending_inline_ring){
            // spiro notation open
            pending_spiro = true; 
          }
          else if(pending_locant){
            // ionic species or spiro, reset the linkings
            prev = 0; 
            pending_locant = false;
          }
          else{
            pop_ticks++; // set the number of pops to do
          }
          break;


        case '-':
          if(!pending_inline_ring){
            pending_inline_ring = true;
            
            // send the linker into its own stack
            if(!branch_stack.empty() && (branch_stack.top()->num_edges < branch_stack.top()->allowed_edges) )
              linker_stack.push(branch_stack.top());
          }
            
          else if(pending_inline_ring){
            fprintf(stderr,"Error: only one pending ring can be active, check closures\n");
            Fatal(i);
          }
          else if(pending_special){

            if(pop_ticks){
              prev = pop_standard_stacks(pop_ticks,branch_stack,linker_stack,prev,i);
              pop_ticks = 0; 
            }

            block_end = i;
            curr = AllocateWLNSymbol('*');
            curr->add_special(block_start,block_end);

            block_start = 0;
            block_end = 0; 

            create_bond(curr,prev,bond_ticks,i);

            bond_ticks = 0; 
            prev = curr;
            pending_special = false;
          }
          else if (!pending_special)
            pending_special = true;
          
          break;

        case '/':
          prev = curr;
          curr = AllocateWLNSymbol(ch);
          break;

        default: 
          fprintf(stderr,"Error: unallowed character! - [A-Z][0-1][&-/' ']\n");
          Fatal(i);
        
      }


    }


    if(pending_closure){
      fprintf(stderr,"Error: expected 'J' to close ring\n");
      Fatal(len);
    }

    if(pending_locant){
      fprintf(stderr,"Error: expected locant to attach to ring\n");
      Fatal(len);
    }

    if(pending_inline_ring){
      fprintf(stderr,"Error: expected inline ring to be defined\n");
      Fatal(len);
    }

    if(pending_spiro){
      fprintf(stderr,"Error: expected sprio ring to be defined\n");
      Fatal(len);
    }
            
    return true; 
  }


  /* dump wln tree to a dotvis file */
  void WLNDumpToDot(FILE *fp)
  {

    std::map<WLNSymbol*, bool>::iterator sym_iter;
    

    fprintf(fp, "digraph WLNdigraph {\n");
    fprintf(fp, "  rankdir = LR;\n");
    for (WLNSymbol *node: symbol_mempool){

      fprintf(fp, "  %d", index_lookup[node]);
      if (node->ch == '*')
        fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
      else
        fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->ch);

      for (WLNSymbol *child : node->children)
      { 
        fprintf(fp, "  %d", index_lookup[node]);
        fprintf(fp, " -> ");
        fprintf(fp, "%d [arrowhead=none]\n", index_lookup[child]);
      }
    }
    fprintf(fp, "}\n");
  }
};





static void DisplayHelp(){
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

  wln = (const char *)0;
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
        wln = ptr;
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
  if (!wln)
  {
    fprintf(stderr, "Error: no wln string - nullptr\n");
    return 1;
  }


  WLNGraph wln_graph;

  wln_graph.ParseWLNString(wln,strlen(wln));
  Reindex_lookups();

  // create the wln dotfile
  if (opt_wln2dot)
  {
    FILE *fp = 0;
    fp = fopen("wln-graph.dot", "w");
    if (!fp){
      fprintf(stderr, "Error: could not open compiler dump file\n");
      fclose(fp);
      return 1;
    }
    else{
      wln_graph.WLNDumpToDot(fp);
      fclose(fp);
    }
      
  }

  return 0;
}