
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector>
#include <stack>
#include <map>
#include <deque>

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
std::map<WLNSymbol *, unsigned int> index_lookup;
std::map<unsigned int, WLNSymbol*> symbol_lookup;
unsigned int glob_index = 0;


// used for assigning charges with post notation
std::map<unsigned int, WLNSymbol*> pending_charge; 

// character type and instruction types
enum WLNType
{
  SINGLETON = 0,
  BRANCH = 1,
  LINKER = 2,
  TERMINATOR = 3,
  SPECIAL = 4
};

enum WLNCode
{
  ROOT = 0,
  STANDARD = 1,
  LOCANT = 2,
  CYCLIC = 3,
  BRIDGED = 4,
  SPIRO = 5,
  IONIC = 6
};

// follows connection table formats
enum WLNBond{
  SINGLE =    1, 
  DOUBLE =    2,
  TRIPLE =    3, 
  AROMATIC =  4
};




const char *code_hierarchy[] = {"ROOT", "STANDARD", "LOCANT", "CYCLIC", "BRIDGED", "SPIRO", "IONIC"};

// rule 2 - hierarchy - rules have diverged due to end terminator char, also use for locant setting from 14
std::map<unsigned char, unsigned int> char_hierarchy =
{
  {' ', 1}, {'-', 2}, {'/', 3}, {'0', 4}, {'1', 5}, {'2', 6}, {'3', 7}, {'4', 8}, {'5', 9}, {'6', 10}, {'7', 11}, {'8', 12}, {'9', 13}, {'A', 14}, {'B', 15}, {'C', 16}, {'D', 17}, {'E', 18}, {'F', 19}, {'G', 20}, {'H', 21}, {'I', 22}, {'J', 23}, {'K', 24}, {'L', 25}, {'M', 26}, {'N', 27}, {'O', 28}, {'P', 29}, {'Q', 30}, {'R', 31}, {'S', 32}, {'T', 33}, {'U', 34}, {'V', 35}, {'W', 36}, {'X', 37}, {'Y', 38}, {'Z', 40}, {'&', 41}
};

std::map<unsigned int, unsigned char> locant_symbols =
{
  {0, 'A'}, {1, 'B'}, {2, 'C'}, {3, 'D'}, {4, 'E'}, {5, 'F'}, {6, 'G'}, {7, 'H'}, {8, 'I'}, {9, 'J'}, {10, 'K'}, {11, 'L'}, {12, 'M'}, {13, 'N'}, {14, 'O'}, {15, 'P'}, {16, 'Q'}, {17, 'R'}, {18, 'S'}, {19, 'T'}, {20, 'U'}, {21, 'V'}, {22, 'W'}, {23, 'X'}, {24, 'Y'}, {25, 'Z'}
};


/*  assumes a bi-atomic fuse, max = 6*6 for bicyclic */
unsigned int calculate_ring_atoms(unsigned int rings, unsigned int max_atoms)
{

  unsigned int term = rings - 2;
  unsigned int shared_atoms = rings + term;

  return max_atoms - shared_atoms;
}


bool isdigit_str(const std::string& s)
{
  for (char const &ch : s) {
    if (std::isdigit(ch) == 0) 
      return false;
  }
  return true;
 }


struct WLNSymbol
{

  unsigned char ch;
  unsigned int type;

  int charge;

  unsigned int inc_bond; // can take values 1-3 for ['','U','UU']

  unsigned int allowed_edges;
  unsigned int num_edges;

  std::string special; // if ch='\0' then a special string is denoted e.g Mg
  // using string to maintain struct ownership

  WLNSymbol *prev; // should be a single term - wln symbol only has one incoming
  
  std::vector<WLNSymbol*> children; // linked list of next terms chains

  // if default needed
  WLNSymbol()
  {
    ch = '\0';
    type = 0;

    inc_bond = 0;
    charge = 0; 

    allowed_edges = 0;
    num_edges = 0;
    prev = 0;
  }

  bool init(unsigned char inp_char)
  {
    ch = inp_char;

    switch (ch)
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
      type = SINGLETON;
      allowed_edges = 4; // allows the U to be handled inline
      break;

    case 'A':
      type = SINGLETON;
      allowed_edges = 2;
      break;

    case 'B': // boron
      type = BRANCH;
      allowed_edges = 3;
      break;

    case 'C': // shortcut carbon atom
      type = BRANCH;
      allowed_edges = 4;
      break;

    case 'D':
      type = SINGLETON;
      allowed_edges = 2;
      break;

    case 'E': // halogens
    case 'F':
    case 'G':
    case 'I':
      type = BRANCH;
      allowed_edges = 3;
      break;

    case 'H': // closing hydrogen
      type = TERMINATOR;
      allowed_edges = 1;
      break;

    case 'J': // generic symbol for halogen
      type = BRANCH;
      allowed_edges = 3;
      break;

    case 'K':
      type = BRANCH;
      allowed_edges = 4;
      break;

    case 'L':
      type = LINKER;
      allowed_edges = 2;
      break;

    case 'M':
      type = BRANCH;
      allowed_edges = 2;
      break;

    case 'N':
      type = BRANCH;
      allowed_edges = 3;
      break;

    case 'O':
      type = SINGLETON;
      allowed_edges = 2;
      break;

    case 'P':
      type = BRANCH;
      allowed_edges = 5;
      break;

    case 'Q':
      type = TERMINATOR;
      allowed_edges = 1;
      break;

    case 'R':           // if handled in standard notation, only one inbound is allowed
      type = TERMINATOR;
      allowed_edges = 1;
      break;

    case 'S':
      type = BRANCH;
      allowed_edges = 6;
      break;

    case 'T':
    case 'U':
      type = LINKER;
      allowed_edges = 2;
      break;

    case 'V':
      type = SINGLETON;
      allowed_edges = 2;
      break;

    case 'W':
      type = LINKER;
      allowed_edges = 2;
      break;

    case 'X':
      type = BRANCH;
      allowed_edges = 4;
      break;

    case 'Y':
      type = BRANCH;
      allowed_edges = 3;
      break;

    case 'Z':
      type = TERMINATOR;
      allowed_edges = 1;
      break;

    case '&':
      type = TERMINATOR;
      allowed_edges = 1;
      break;

    case ' ':
    case '-':
    case '/':
      type = LINKER;
      allowed_edges = 2;
      break;

    case '\0':
      fprintf(stderr, "Error: end of string null char accessed!\n");
      return false;

    case '*':
      type = SPECIAL;
      break;

    default:
      fprintf(stderr, "Error: invalid wln symbol parsed: %c\n", ch);
      return false; 
    }

    prev = (WLNSymbol *)0;
    return true;
  }
};

/* struct to hold pointers for the wln ring - only for stack return */
struct WLNRing
{

  WLNSymbol *rhead;
  unsigned int ring_size;

  bool aromatic;
  bool heterocyclic;

  std::map<unsigned char, WLNSymbol *> locants;

  WLNRing()
  {
    rhead = 0;
    ring_size = 0;
    aromatic = false;
    heterocyclic = false;
  }

  void init()
  {
    rhead = 0;
    ring_size = 0;
    aromatic = false;
    heterocyclic = false;
  }

  void debug_map()
  {
    std::map<unsigned char, WLNSymbol *>::iterator iter;
    for (iter = locants.begin(); iter != locants.end(); iter++)
    {
      fprintf(stderr, "%p ---> %c\n", iter->second, iter->first);
    }
  }
};

struct WLNGraph{

  WLNSymbol *root;
  unsigned int wln_nodes = 0;
  unsigned int wln_rings = 0;
  std::vector<WLNSymbol *> symbol_mempool;
  std::vector<WLNRing *> ring_mempool;
  
  std::map<WLNSymbol*,bool> symbol_hide; 
  std::map<WLNRing *, WLNSymbol *> ring_access; // access the ring struct from locant A pointer

  WLNGraph() : root{(WLNSymbol *)0}, wln_nodes{0} {};
  ~WLNGraph()
  {
    for (WLNSymbol *allocedwln : symbol_mempool)
    { // if error free all the symbol_mempool --> stop leak
      delete allocedwln;
      allocedwln = 0;
    }

    for (WLNRing *allocedring : ring_mempool)
    { // if error free all the symbol_mempool --> stop leak
      delete allocedring;
      allocedring = 0;
    }
  }

  
  WLNSymbol *AllocateWLNSymbol(unsigned char ch)
  {
    wln_nodes++;
    WLNSymbol *wln = new WLNSymbol;
    if (wln->init(ch))
      symbol_mempool.push_back(wln);
    else
      return (WLNSymbol *)0;


    // add to globals --> needed for charge assignment
    index_lookup[wln] = glob_index;
    symbol_lookup[glob_index] = wln;
    glob_index++; 

    return wln;
  }

  // for a manual psuedo style deallocate
  void HideWLNSymbol(WLNSymbol *node){
    symbol_hide[node] = true;
  }

  WLNRing *AllocateWLNRing()
  {
    wln_rings++;
    WLNRing *wln_ring = new WLNRing;
    wln_ring->init();
    ring_mempool.push_back(wln_ring);
    return wln_ring;
  }

  void reset_indexes(){
    glob_index = 0; 
    for (WLNSymbol* node : symbol_mempool){
      if (symbol_hide[node])
        continue;

      index_lookup[node] = glob_index;
      symbol_lookup[glob_index] = node;
      glob_index++;
    }
  }

  /* re init a transformed symbol */
  WLNSymbol *transform_symbol(WLNSymbol *sym, unsigned char ch)
  {
    sym->init(ch);
    return sym;
  }

  
  bool copy_symbol_info(WLNSymbol *src, WLNSymbol *trg, bool SYM = false){
    
    if(SYM)
      trg->ch = src->ch;

    trg->type = src->type;
    trg->charge = src->charge;
    trg->inc_bond = src->inc_bond;
    trg->allowed_edges = src->allowed_edges;
    trg->num_edges = src->num_edges;
    trg->prev = src->prev;
    
 
    for(WLNSymbol *child: src->children){
      trg->children.push_back(child);
    }

    return true;
  }

  WLNSymbol *access_locant(unsigned char ch, WLNRing *ring, bool strict=true)
  {
    WLNSymbol *locant = 0;
    locant = ring->locants[ch];
    if (!locant)
    {
      if(strict)
        fprintf(stderr, "Error: invalid locant access - %c\n", ch);
      return 0;
    }
    return locant;
  }

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
          atom->inc_bond++;
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
          if(!add_symbol(oxy,atom,1))
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
          if(!add_symbol(oxy_1,atom,1) ||  !add_symbol(oxy_2,atom,1))
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

  /* inplace function, should only edit the WLNRing pointer */
  bool CreateStandardRing(unsigned int start, unsigned int end, WLNRing *ring)
  {

    if (!ring){
      fprintf(stderr, "Error: ring object incorrectly made!\n");
      return false;
    }

    unsigned int num_atoms = 0;
    unsigned int num_rings = 0;
    unsigned int digit_end = 0;

    std::vector<unsigned int> fuse_pattern;

    // 1) evaluate the number of rings

    unsigned int it = start + 1; // get the first num
    while (std::isdigit(wln[it]) && wln[it] != '\0')
    {
      unsigned int val = wln[it] - '0';
      num_atoms += val;
      num_rings++;
      fuse_pattern.push_back(val);
      it++;
    }

    digit_end = it;
    unsigned int ratoms = calculate_ring_atoms(num_rings, num_atoms);

    // 1) create a big ring
    WLNSymbol *rhead = AllocateWLNSymbol('C');
    ring->rhead = rhead; // set the rings head

    WLNSymbol *current = 0;
    WLNSymbol *prev = rhead;

    unsigned int locant = 0;
    for (unsigned int i = 1; i < ratoms; i++)
    {
      current = AllocateWLNSymbol('C');
      ring->locants[locant_symbols[locant++]] = current; // add the locants
      
      if(ring->aromatic)
        add_aromatic(current,prev);
      else
        add_symbol(current, prev,0);

      prev = current;
    }

    if(ring->aromatic)
      add_aromatic(rhead,current);
    else
      add_symbol(rhead, current,0);

    if (num_rings > 1)
    {
      // handle bicyclic fuse patterns here
    }

    // handle all inter atomic definitions here
    if(!ParseInterRing(digit_end, end, ring))
      return false; 
    
    
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
  

  /* used in ring notation only*/
  bool add_aromatic(WLNSymbol *child, WLNSymbol *parent){
    child->inc_bond   = 4; 
    child->num_edges  += 1;
    parent->num_edges += 1;
    parent->children.push_back(child);
    child->prev = parent; // keep the linked list functionality.
    return true;
  }

  /* should handle all bonding modes, adds child to parent->children
  'UU' bonding also added here */
  bool add_symbol(WLNSymbol *child, WLNSymbol *parent, unsigned int bond_ticks)
  {


    unsigned int bond_added = 1 + bond_ticks; // can be zero - local scope so reset outside

    // if the child cannot handle the new valence
    if ( (child->num_edges + bond_added) > child->allowed_edges ){
      fprintf(stderr,"Error: wln character[%c] is exceeding allowed connections\n",child->ch);
      return false; 
    }

    // same for the parent
    if ( (parent->num_edges + bond_added) > parent->allowed_edges ){
      fprintf(stderr,"Error: wln character[%c] is exceeding allowed connections\n",parent->ch);
      return false; 
    }

    
    // if these pass, we can add and change both num_edges 
    
    child->inc_bond   += bond_added; 

    child->num_edges  += bond_added;
    parent->num_edges += bond_added;

    parent->children.push_back(child);

    child->prev = parent; // keep the linked list functionality.

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

    /* this can perform the normal backtrack, excluding  '&' closure */
  WLNSymbol *backtrack_stack(std::stack<WLNSymbol *> &wln_stack)
  {
    WLNSymbol *tmp = 0;
    while (!wln_stack.empty())
    {
      tmp = wln_stack.top();
      if (tmp->type == BRANCH)
        return tmp;
      wln_stack.pop();
    }

    // if it goes all the way that string complete
    return (WLNSymbol *)0;
  }

  /* forces both the '&' closure and its parents branch to come off stack */
  WLNSymbol *force_closure(std::stack<WLNSymbol *> &wln_stack)
  {
    WLNSymbol *tmp = 0;
    unsigned int popped = 0;
    while (!wln_stack.empty())
    {
      tmp = wln_stack.top();
      if (tmp->type == BRANCH && popped > 1)
        return tmp;
      wln_stack.pop();
      popped++;
    }

    // if it goes all the way that string complete
    return (WLNSymbol *)0;
  }

  // can return head or tail
  WLNSymbol* consume_standard_notation2(unsigned int start, unsigned int end, bool tail = false){


    std::stack<WLNSymbol *> wln_stack;
    
    WLNSymbol *created_wln = AllocateWLNSymbol(wln[start]);
    WLNSymbol *prev = created_wln;
    WLNSymbol *root = created_wln;
    
    if (!created_wln)
      return (WLNSymbol *)0;
    
    if(created_wln->type == BRANCH)
      wln_stack.push(created_wln);

    bool open_special = false; 
    unsigned int bond_tick = 0;    
    std::vector<unsigned char> special;

    // stack rework, only push branching and have a condition for popping
    for (unsigned int i = start + 1; i <= end; i++){
      
      if(open_special && wln[i] != '-'){
        special.push_back(wln[i]);
        if(special.size() > 2){
          fprintf(stderr,"Error: invalid elemental notation in standard\n");
          return (WLNSymbol *)0;
        }
        continue;
      }

      if (wln[i] == 'U'){
        bond_tick++;
        continue; 
      }
      else if(wln[i] == '-'){
        if(!open_special){
          open_special = true;
          continue;
        }   
        if(open_special){
          created_wln = define_element(special);
          special.clear();
          open_special = false;
        }
      }
      else if(wln[i] == '&'){

        if (wln[i-1] == '&'){
          if(wln_stack.size() > 1){
            wln_stack.pop();
            prev = wln_stack.top();
          } else{ 
            fprintf(stderr,"Error: branching stack exhausted - extra '&' in notation\n");
            return (WLNSymbol *)0;
          }
        }
        else{
          if(!wln_stack.empty()){
            prev = wln_stack.top();
          }
          else{ 
            fprintf(stderr,"Error: branching stack exhausted - extra '&' in notation\n");
            return (WLNSymbol *)0;
          }
        }
        continue; 
      }
      else
        created_wln = AllocateWLNSymbol(wln[i]);

      if(created_wln->type == BRANCH)
        wln_stack.push(created_wln);

      // add the bond, and move prev across
      if (!add_symbol(created_wln, prev,bond_tick))
        return (WLNSymbol *)0;


      bond_tick = 0; // reset the bond counter;  

      // if bonded out then we pop off the stack due to bonding condition
      if(!wln_stack.empty() && prev == wln_stack.top()){
        if(prev->allowed_edges == prev->num_edges)
          wln_stack.pop();
      }

      // if a terminator we have to return to the stack
      if(!wln_stack.empty() && created_wln->type == TERMINATOR)
        prev = wln_stack.top();
      else
        prev = created_wln; 
    }
    
    // head or tail return

    if(tail)
      return created_wln;
    else
      return root; 

  }


#ifdef DEPRECATED
  /* consumes the standard blocks, uses stacks to handle branching */
  WLNSymbol *consume_standard_notation(unsigned int start, unsigned int end)
  {

    std::stack<WLNSymbol *> wln_stack;
    WLNSymbol *prev = 0;
    WLNSymbol *root = 0;

    WLNSymbol *created_wln = AllocateWLNSymbol(wln[start]);
    if (!created_wln)
      return (WLNSymbol *)0;
    wln_stack.push(created_wln);

    root = created_wln;

    bool open_special = false; 
    unsigned int bond_tick = 0;    
    std::vector<unsigned char> special;


    for (unsigned int i = start + 1; i <= end; i++)
    {

      // bonds
      if (wln[i] == 'U'){
        bond_tick++;
        continue; 
      }

      if(open_special && wln[i] != '-'){
        special.push_back(wln[i]);
        if(special.size() > 2){
          fprintf(stderr,"Error: invalid elemental notation in standard\n");
          return (WLNSymbol *)0;
        }
        continue;
      }

      // specials
      if(wln[i] == '-'){
        if(!open_special){
          open_special = true;
          continue;
        }   
        if(open_special){
          created_wln = define_element(special);
          special.clear();
          open_special = false;
        }
      }
      else
        created_wln = AllocateWLNSymbol(wln[i]);

      if (!created_wln)
        return (WLNSymbol*)0;

    
      if(wln_stack.empty()){
        fprintf(stderr,"Error: invalid branch notation - invalid tertiary access\n");
        return (WLNSymbol*)0;
      }
      else
        prev = wln_stack.top();

    
      wln_stack.push(created_wln); // push all of them

      if (!add_symbol(created_wln, prev,bond_tick))
        return (WLNSymbol *)0;

      bond_tick = 0; // reset the bond counter;  
     
      // options here depending on the forced closure of '&'
      if (created_wln->type == TERMINATOR)
      {
        if (created_wln->ch == '&' && prev->type == BRANCH)
          prev = force_closure(wln_stack);
        else
          prev = backtrack_stack(wln_stack);
      }
    }

    return root; // return start of the tree
  }

#endif


  bool create_chain(WLNSymbol *node, bool special=false){
    
    unsigned int atoms = 0; 
    if (special)
      atoms = std::stoi(node->special);
    else
      atoms = node->ch - '0';

    node = transform_symbol(node,'C'); // 1) transform the character

    // 2) create the chain
    WLNSymbol *prev = 0; 
    WLNSymbol *head = 0; 
    for (unsigned int k=0;k<atoms-1;k++){
      WLNSymbol *created = AllocateWLNSymbol('C');
      if(prev)
        add_symbol(prev,created,0);
      else
        head = created; 
      prev = created;
    }
    
    // 3) last prev will be the tail, copy over heads details
    copy_symbol_info(node,prev);
    
    // 4) remove children from node
    node->children.clear();
    // 5) bind the new chain
    node->children.push_back(head);

    return true;
  }

  // to be used on a mono substitued R benzene only
  // not in place uses hide mechanism. 
  bool create_benzene(WLNSymbol *node){
    
    symbol_hide[node] = true; 

    WLNSymbol *head = AllocateWLNSymbol('C');
    WLNSymbol *prev = head;

    for (unsigned int k=0;k<5;k++){
      WLNSymbol *created = AllocateWLNSymbol('C');
      add_aromatic(prev,created);
      prev = created;
    }

    add_aromatic(head,prev);
    prev->inc_bond = 1;

    

    // if it has a previous we bond back - if not, it has to either bond forward or to nothing
    if(node->prev)
      node->prev->children.push_back(prev);
    else if(!node->children.empty()){
      for (WLNSymbol *move : node->children)
        prev->children.push_back(move);
    }
    else
      prev->charge = -1; 
    

    return true;
  }

  /* search the mempool dfs style and find all concat points */
  bool ConcatNumerics(){

    std::stack<WLNSymbol*> node_stack; 
    std::map<WLNSymbol*,bool> visited; // avoid the loop issue i know is coming

    node_stack.push(symbol_mempool[0]);


    std::string chain; 
    std::vector<WLNSymbol*> streak; 

    WLNSymbol *head = 0; 
    WLNSymbol *tail = 0; 

    WLNSymbol *node = 0; 
    while(!node_stack.empty()){
      
      node = node_stack.top();
      node_stack.pop();
      visited[node] = true;

      if(std::isdigit(node->ch)){

        if(chain.empty())
          head = node;
        else  
          tail = node; 

        streak.push_back(node);
        chain.push_back(node->ch);
      }
        
      else{
        if (chain.size() > 1){
          WLNSymbol* chain_symbol = AllocateWLNSymbol('*');
          chain_symbol->special = chain;
          copy_symbol_info(tail,chain_symbol);
          if (head->prev){
            head->prev->children.push_back(chain_symbol);
          }
          for (WLNSymbol* n : streak)
            HideWLNSymbol(n);
        }

        head = 0; 
        tail = 0; 
        chain.clear();
        streak.clear();
      }

      for (WLNSymbol *child : node->children){
        if(child && !visited[child])
          node_stack.push(child);
      }

    }

    return true; 
  }

  // expands the graph to suite a pseudo smiles for conversion
  bool ExpandGraph(){
    
    unsigned int start_size = symbol_mempool.size(); // doesnt change
    for (unsigned int i=0; i<start_size;i++){
      WLNSymbol *node = symbol_mempool[i];
      if(symbol_hide[node])
        continue;

      switch (node->ch){

        case '1':
          node->ch = 'C';
          break; 

        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':{
          create_chain(node);
          break;
        }


        case 'E':
          node = transform_symbol(node,'*');
          node->special = "Br";
          break;

        case 'G':
          node = transform_symbol(node,'*');
          node->special = "Cl";
          break;
        

        case 'K':
          node = transform_symbol(node,'N');
          node->charge = 1;
          break;

        case 'M':{
          node = transform_symbol(node,'N');
          WLNSymbol *created = AllocateWLNSymbol('H');
          add_symbol(created,node,0);
          break;
        }

        case 'Z':{
          node = transform_symbol(node,'N');
          for (unsigned int k = 0; k < 2; k++){
            WLNSymbol *created = AllocateWLNSymbol('H');
            add_symbol(created,node,0);
          }
          break; 
        }

        case 'Q':{
          // expands to O-H
          node = transform_symbol(node,'O');
          WLNSymbol *created = AllocateWLNSymbol('H');
          add_symbol(created,node,0);
          break;
        }

        case 'V':{ 
          // expands to C=O
          node = transform_symbol(node,'C');
          WLNSymbol *created = AllocateWLNSymbol('O');
          add_symbol(created,node,1);
          break;
        }

        case 'Y':
        case 'X':{
          // expands to carbon, details should be kept about charge
          node->ch = 'C';
          break; 
        }

        case 'R':{
          create_benzene(node);
          break;
        }

       
        case 'W':
          fprintf(stderr,"Too handle!\n");
          break;


        case '&':
          HideWLNSymbol(node);
          break;


        case '*':
          if (isdigit_str(node->special))
            create_chain(node,true);
          
          break;


        // do nothings
        case 'F':
        case 'H':
        case 'I':
        case 'J':
        case 'L':
        case 'T':
        case 'U':
        case 'S':
          break;

        default:
          fprintf(stderr,"Error: unexpected char in graph expansion - %c\n",node->ch);
          break;
      }
    }

    return true; 
  }

  /*prints the WLN connection table in SCT XI format */
  void WLNConnectionTable(FILE *fp){

    /* 
    -- atom table ---
    |index| |type| |charge|

    --- bond table --- 
    |atom1| |atom2| |order| 
    */


    fprintf(fp, "---- atom table ----\n");
    fprintf(fp,"|index|\t|type|\t|charge|\n");
    for (WLNSymbol *node : symbol_mempool){
      
      if(symbol_hide[node])
        continue;

      if(node->ch == '*')
        fprintf(fp, "%d\t%s\t%d\n", index_lookup[node], node->special.c_str(),node->charge);
      else
        fprintf(fp, "%d\t%c\t%d\n", index_lookup[node], node->ch,node->charge);
      
    }
    fprintf(fp,"\n");
    

    fprintf(fp, "---- bond table ----\n");
    fprintf(fp,"|atom 1|\t|atom 2|\t|order|\n");
    for (WLNSymbol *node : symbol_mempool){
      if(symbol_hide[node])
        continue;

      for (WLNSymbol *child : node->children){
        if(symbol_hide[child])
          continue;

        if(!child->inc_bond){
          fprintf(stderr,"Error: undefined bond written into connection table\n");
          return ;
        }else
          fprintf(fp, "%d\t%d\t%d\n", index_lookup[node],index_lookup[child], child->inc_bond);
      }
    }    

    fprintf(fp,"\n");
  }

  /* dump wln tree to a dotvis file */
  void WLNDumpToDot(FILE *fp)
  {

    fprintf(fp, "digraph WLNdigraph {\n");
    fprintf(fp, "  rankdir = LR;\n");
    for (WLNSymbol *node : symbol_mempool)
    { 
      if(symbol_hide[node])
        continue;

      fprintf(fp, "  %d", index_lookup[node]);
      if (node->ch == '*')
        fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
      else
        fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->ch);

      for (WLNSymbol *child : node->children)
      { 
        if(symbol_hide[child])
          continue;

        if(child->inc_bond == SINGLE){
          fprintf(fp, "  %d", index_lookup[node]);
          fprintf(fp, " -> ");
          fprintf(fp, "%d [arrowhead=none]\n", index_lookup[child]);
        }
        else if (child->inc_bond == DOUBLE || child->inc_bond == TRIPLE){
          for (unsigned int i=0; i<child->inc_bond; i++){
            fprintf(fp, "  %d", index_lookup[node]);
            fprintf(fp, " -> ");
            fprintf(fp, "%d [arrowhead=none]\n", index_lookup[child]);
          }
        } 
        else if (child->inc_bond == AROMATIC){
          fprintf(fp, "  %d", index_lookup[node]);
          fprintf(fp, " -> ");
          fprintf(fp, "%d [arrowhead=none,color=blue]\n", index_lookup[child]);
        }
        else if (child->inc_bond == 0){
          fprintf(stderr,"Warning: plotting undefined bond, shown in red\n");
          fprintf(fp, "  %d", index_lookup[node]);
          fprintf(fp, " -> ");
          fprintf(fp, "%d [color=red]\n", index_lookup[child]);
        }
        
      }
    }
    fprintf(fp, "}\n");
  }
};

struct WLNInstruction
{

  unsigned int state;
  unsigned int start_ch;
  unsigned int end_ch;

  bool ring_linker;

  WLNInstruction *parent;
  std::vector<WLNInstruction *> next_instructions;

  void init_state(unsigned int c)
  {
    state = c;
    start_ch = 0;
    end_ch = 0;
    parent = (WLNInstruction *)0;
    ring_linker = false;
  }

  void add_start(unsigned int pos)
  {
    start_ch = pos;
  }

  void add_end(unsigned int pos)
  {
    end_ch = pos;
  }

  void add_prev(WLNInstruction *src)
  {
    parent = src;
  }

  void display()
  {
    if (state == ROOT)
      fprintf(stderr, "instruction: %10s\n", "ROOT");
    else if (state == LOCANT)
      fprintf(stderr, "instruction: %10s contains: %c\n", code_hierarchy[state], wln[start_ch]);
    else
    {
      fprintf(stderr, "instruction: %10s contains: ", code_hierarchy[state]);
      for (unsigned int i = start_ch; i <= end_ch; i++)
        fprintf(stderr, "%c", wln[i]);
      fprintf(stderr, "\n");
    }
  }
};

/* make a graph, split the string, call the subroutines - easy right? */
struct WLNParser
{

  WLNInstruction *root;
  unsigned int num_instructions;
  std::vector<WLNInstruction *> instruction_pool;

  WLNParser()
      : root{(WLNInstruction *)0}, num_instructions{0} {};
  ~WLNParser()
  {
    for (WLNInstruction *instruction : instruction_pool)
      delete instruction;
  }

  WLNInstruction *add_instruction(unsigned int int_code, unsigned int i)
  {
    WLNInstruction *instruction = new WLNInstruction;
    instruction->init_state(int_code);
    instruction->add_start(i);

    instruction_pool.push_back(instruction);
    num_instructions++;
    return instruction;
  }

  void display_instructions()
  {
    for (WLNInstruction *instruction : instruction_pool)
      instruction->display();
  }

  WLNSymbol *return_locant_symbol(unsigned int wln_pos, std::stack<WLNRing *> &ring_stack)
  {

    if (ring_stack.empty())
    {
      fprintf(stderr, "Error: accessing ring notation with zero rings\n");
      return 0;
    }

    WLNRing *ring = ring_stack.top();
    WLNSymbol *binder = ring->locants[wln[wln_pos]];

    if (!binder)
    {
      fprintf(stderr, "Error: accessed out of bound locant position\n");
      return 0;
    }

    return binder;
  }

  /* parse the wln string and create instruction set,
  I think its reasonable to have two parses for this */
  bool CreateWLNGraph(const char *wln, unsigned int len, WLNGraph &graph)
  {

    WLNInstruction *current = add_instruction(ROOT, 0);
    root = current;

    // these are now global for inline graph creation
    std::stack<WLNRing *> ring_stack;
    std::stack<WLNSymbol *> rlinker_stack; // for XR&R style ring definition

    bool pending_closure  = false;
    bool pending_locant   = false;
    bool pending_benzene  = false; // for substituted R notation

    WLNSymbol *binder = 0; // used to link the wln nodes to previous

    for (unsigned int i = 0; i < len; i++)
    {
      char ch = wln[i];
      switch (ch)
      {

      case 'L':
      case 'T':

        if (current->state == ROOT)
        {

          current = add_instruction(CYCLIC, i);

          pending_closure = true;
        }
        else if (current->state == STANDARD)
        {

          if (pending_locant)
          {
            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false; 
            
            pending_locant = false;
          }
        }
        else if (current->state == LOCANT)
        {

          current = add_instruction(CYCLIC, i);

          pending_closure = true;
        }
        else if (current->state == CYCLIC)
        {

          if (pending_locant)
          {

            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false; 

            pending_locant = false;
          }
        }

        break;

      case 'J': // pass
        if (current->state == STANDARD)
        {

          if (pending_locant)
          {
            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false; 

            pending_locant = false;
          }
        }
        else if (current->state == LOCANT || current->state == IONIC)
        {

          current = add_instruction(STANDARD, i);
        }
        else if (current->state == CYCLIC)
        {

          if (pending_closure)
          {
            current->add_end(i);
            // create the ring object

            WLNRing *ring = graph.consume_ring_notation(current->start_ch, current->end_ch);
            if(!ring)
              return false;
            
            ring_stack.push(ring);

            pending_closure = false;
          }
          else if (pending_locant)
          {

            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false; 

            pending_locant = false;
          }
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
      case 'S':
      case 'U':
      case 'V':
      case 'W':
      case 'X':
      case 'Y':
      case 'Z':
        if (current->state == ROOT)
        {
          current = add_instruction(STANDARD, i);
        }
        else if (current->state == STANDARD)
        {

          if (pending_locant)
          {

            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false;  

            pending_locant = false;
          }
        }
        else if (current->state == LOCANT || current->state == IONIC)
        {

          current = add_instruction(STANDARD, i);
        }
        else if (current->state == CYCLIC)
        {

          if (pending_locant)
          {

            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false; 

            pending_locant = false;
          }
        }

        break;


      case 'R': // we need to be able to place a pending system here for substituted benzene
        if (current->state == ROOT)
        { 
          if(i==0 && len == 1){
            current = add_instruction(STANDARD, i);
          }
          else
            pending_benzene = true;
        }
        else if (current->state == STANDARD)
        {

          if (pending_locant)
          {

            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false;  

            pending_locant = false;
          }
          else
            pending_benzene = true;
        }
        else if (current->state == LOCANT || current->state == IONIC)
        {
          pending_benzene = true;
        }
        else if (current->state == CYCLIC)
        {

          if (pending_locant)
          {

            current = add_instruction(LOCANT, i);
            current->add_end(i); // all locants terminate on one char

            binder = return_locant_symbol(i, ring_stack);
            if (!binder)
              return false; 

            pending_locant = false;
          }
        }

        break;

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

        if (current->state == ROOT)
        {
          if(pending_benzene){
            current = add_instruction(STANDARD, i-1);
            pending_benzene = false;
          }
          else
            current = add_instruction(STANDARD, i);
        }
        else if (current->state == LOCANT || current->state == IONIC)
        {

          current = add_instruction(STANDARD, i);
        }
        else if (current->state == CYCLIC)
        {

          if (wln[i - 1] == '&')
          {
            // attach a standard here, other rules handle ring,
            // match on no symbol for easy merge.

            unsigned int k = i - 1;
            unsigned int terms = 0;
            while (wln[k] == '&')
            {
              terms++;
              k--;
            }

            if (!current)
            {
              fprintf(stderr, "Error: no ring linker to return to via '&<x>-' - terminating parse\n");
              return false;
            }

            current = add_instruction(STANDARD, i);
          }
        }

        break;

      case ' ': // keep this simple, a '&' locant means ionic branch out
        if (current->state == ROOT){
          if(pending_benzene){
            current = add_instruction(CYCLIC,i-1);
            current->add_end(i-1); // the R is the cyclic notation
            WLNRing *ring = graph.consume_benzene();
            ring_stack.push(ring);
            pending_locant = true;
            pending_benzene = false;
          }
        }
        else if (current->state == STANDARD)
        {
          if(pending_benzene){

            current->add_end(i-2); // look back close
            WLNSymbol *tail = graph.consume_standard_notation2(current->start_ch,current->end_ch,true); // consume the chain before

            current = add_instruction(CYCLIC,i-1);
            current->add_end(i-1); // the R is the cyclic notation
            WLNRing *ring = graph.consume_benzene();
            ring_stack.push(ring);

            // bind the ring to the tail of the chain
            ring->rhead->children.push_back(tail);

            pending_locant = true;
            pending_benzene = false;
          }
          else{
          
            current->add_end(i - 1); // implied end of standard notation block

            // will need to perform certain checks here
            // for interplaying ring notation but build this slowly

            WLNSymbol *head = graph.consume_standard_notation2(current->start_ch, current->end_ch);
            if(!head)
              return false; 
            
            if (binder){
              binder->children.push_back(head);
              head->inc_bond = 1;
            }
              
            pending_locant = true;
          }
        }

        else if (current->state == LOCANT)
        {

          ;
        }

        else if (current->state == CYCLIC)
        {

          if (!pending_closure)
            pending_locant = true;
        }

        break;

      case '-':

        if (current->state == ROOT)
        {

          current = add_instruction(STANDARD, i);
        }
        else if (current->state == STANDARD)
        {

          ;
        }
        else if (current->state == LOCANT)
        {

          ;
        }
        else if (current->state == CYCLIC)
        {

          if (wln[i - 1] == '&')
          {
            // means a return to a branching linker - also pops off ring
            unsigned int k = i - 1;
            unsigned int terms = 0;
            while (wln[k] == '&')
            {
              terms++;
              k--;
            }

            if (!current)
            {
              fprintf(stderr, "Error: no ring linker to return to via '&<x>-' - terminating parse\n");
              return false;
            }
          }
        }

        else if (current->state == IONIC)
        {

          current = add_instruction(STANDARD, i);
        }

        break;

      case '&':
        if (current->state == STANDARD)
        {

          if (pending_locant)
          {
            current = add_instruction(IONIC, i); // ionic is always seperate in the graph
            current->add_end(i);

            // an ionic species means a complete blow through of the ring stack
            while (!ring_stack.empty())
              ring_stack.pop();

            pending_locant = false;
          }
        }
        else if (current->state == CYCLIC)
        {

          if (pending_locant)
          {
            current = add_instruction(IONIC, i); // ionic is always seperate in the graph
            current->add_end(i);

            // an ionic species means a complete blow through of the ring stack
            while (!ring_stack.empty())
              ring_stack.pop();

            pending_locant = false;
          }
        }
        break;

      default:
        fprintf(stderr, "Error: unrecognised symbol: %c\n", ch);
        return false;  
      }
    }


    if(!current){
      fprintf(stderr,"Error: no states could be assigned - broken notation\n");
      return false;
    }
    if (current->state == ROOT){
      fprintf(stderr,"Error: no states could be assigned - broken notation\n");
      return false;
    }


    // whatever the last instruction was, add len-1 as the end
    current->add_end(len - 1);

    // this also handles completely standard notation

    if (current->state == STANDARD)
    {
      WLNSymbol *head = graph.consume_standard_notation2(current->start_ch, current->end_ch);
      if (!head)
        return false; 

      if (binder){
        binder->children.push_back(head);
        head->inc_bond = 1;
      }
        
    }

    return true;
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
  WLNParser parser;

  // parse should exit 1 at errors
  if(!parser.CreateWLNGraph(wln, strlen(wln), wln_graph))
    return 1; 

  if (opt_debug)
    parser.display_instructions();

    
  if(opt_convert){
    wln_graph.ConcatNumerics();
    wln_graph.ExpandGraph();
    wln_graph.reset_indexes();
    
    if(opt_debug)
      wln_graph.WLNConnectionTable(stderr);
  }
  
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