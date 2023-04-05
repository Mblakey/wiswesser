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

#include <vector>
#include <stack>
#include <map>
#include <utility> // std::pair
#include <iterator>
#include <set>
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
#define INF 9999

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
struct WLNRing;

std::map<WLNSymbol *, unsigned int> index_lookup;
std::map<unsigned int, WLNSymbol *> symbol_lookup;

std::map<unsigned int,OpenBabel::OBAtom*> babel_atom_lookup;

unsigned int glob_index = 1; // babel starts from 1, keep consistent  

// --- pools ---
std::vector<WLNSymbol *>  symbol_mempool;
std::vector<WLNRing *>    ring_mempool;


// allows comp as string type
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

static bool isdigit_str(const std::string &s)
{
  for (char const &ch : s)
  {
    if (std::isdigit(ch) == 0)
      return false;
  }
  return true;
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

static void Reindex_lookups()
{
  glob_index = 1;
  for (WLNSymbol *node : symbol_mempool)
  {
    index_lookup[node] = glob_index;
    symbol_lookup[glob_index] = node;
    glob_index++;
  }
}



struct WLNSymbol
{
  std::string value; // allows easy modifiers
  std::string special; // string for element, or ring, if value = '*'

  // will deprecate special DEFO, and switch the first char of 
  // the value for graph expansion

  unsigned int type;
  unsigned int allowed_edges;
  unsigned int num_edges;

  WLNSymbol *previous;
  std::vector<WLNSymbol *> children; 
  std::vector<unsigned int> orders;

  // if default needed
  WLNSymbol()
  {
    value = "";
    allowed_edges = 0;
    num_edges = 0;
    previous = 0;
  }
  ~WLNSymbol(){};

  unsigned char get_ch(){
    if(value.empty()){
      fprintf(stderr,"Error: accessing empty value\n");
      return 0;
    }
    else
      return value[0];
  }

  void add_modifier(unsigned char mod){
    value.push_back(mod);
  }

  void set_edges(unsigned int edges)
  {
    allowed_edges = edges;
  }

  void set_type(unsigned int i)
  {
    type = i;
  }


  void add_special(unsigned int s, unsigned int e)
  {
    for (unsigned int i = s; i <= e; i++)
      special.push_back(wln[i]);
  }

  void add_special(std::string str)
  {
    special.append(str);
  }
};

WLNSymbol *AllocateWLNSymbol(unsigned char ch)
{

  WLNSymbol *wln = new WLNSymbol;
  symbol_mempool.push_back(wln);
  wln->value = ch;

  index_lookup[wln] = glob_index;
  symbol_lookup[glob_index] = wln;
  glob_index++;
  return wln;
}

// these are expensive, but needed for some edge case notations
void DeallocateWLNSymbol(WLNSymbol *node)
{
  if (opt_debug)
    fprintf(stderr, "  manual deallocation: %c\n", node->get_ch());

  // find the node in the mem pool
  unsigned int i = 0;
  for (WLNSymbol *n : symbol_mempool)
  {
    if (n == node)
      break;
    i++;
  }

  symbol_mempool.erase(symbol_mempool.begin() + i);
  delete node;
}

WLNSymbol *copy_symbol(WLNSymbol *src)
{

  WLNSymbol *copy = AllocateWLNSymbol(src->get_ch());
  copy->allowed_edges = src->allowed_edges;
  copy->num_edges = src->num_edges;

  for (unsigned int i = 0; i < src->children.size(); i++)
  {
    copy->children.push_back(src->children[i]);
    copy->orders.push_back(src->orders[i]);
  }
  return copy;
}

/* should handle all bonding modes, adds child to parent->children
'UU' bonding also added here - needs a bond specied! 1 for single */
bool link_symbols(WLNSymbol *child, WLNSymbol *parent, unsigned int bond, bool aromatic = false)
{

  if(!child || !parent){
    fprintf(stderr,"Error: attempting bond of non-existent symbols\n");
    return false;
  }

  // if the child cannot handle the new valence
  if ((child->num_edges + bond) > child->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->get_ch(),child->num_edges+bond, child->allowed_edges);
    return false;
  }
  // same for the parent
  if ((parent->num_edges + bond) > parent->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->get_ch(),parent->num_edges+bond, parent->allowed_edges);
    return false;
  }

  child->previous = parent; // keep the linked list so i can consume backwards rings
  child->num_edges += bond;
  parent->num_edges += bond;
  parent->children.push_back(child);
  
  if(aromatic)
    parent->orders.push_back(4);
  else
    parent->orders.push_back(bond);
  return true;
}

// use sparingly, loop check isnt ideal 
bool change_symbol_order(WLNSymbol *child, WLNSymbol* parent,unsigned int bond){

  if(!child || !parent){
    fprintf(stderr,"Error: attempting bond change of non-existent symbols\n");
    return false;
  }

  bool found = false;
  WLNSymbol *local_child = 0; 
  unsigned int i=0;
  for (i=0; i<parent->children.size();i++){
    local_child = parent->children[i];
    if(local_child == child){
      found = true;
      break;
    }
  }

  if(!found){
    fprintf(stderr,"Error: changing bond order of non-existent link\n");
    return false; 
  }

  // check can we increase the order, then access orders list with i. 

  unsigned int current_order = parent->orders[i]; 
  if(current_order == bond)
    return true; // save some work
  int diff = bond - current_order; // can be negative for a decrease in order
  // same checks
  if ((child->num_edges + diff) > child->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->get_ch(),child->num_edges+diff, child->allowed_edges);
    return false;
  }
  
  if ((parent->num_edges + diff) > parent->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->get_ch(),parent->num_edges+diff, parent->allowed_edges);
    return false;
  }

  child->num_edges += diff;
  parent->num_edges += diff;
  parent->orders[i] = bond;
  return true; 
}

bool increase_bond_order(WLNSymbol *child, WLNSymbol *parent){

  if(!child || !parent){
    fprintf(stderr,"Error: attempting bond increase of non-existent symbols\n");
    return false;
  }

  bool found = false;
  WLNSymbol *local_child = 0; 
  unsigned int i=0;
  for (i=0; i<parent->children.size();i++){
    local_child = parent->children[i];
    if(local_child == child){
      found = true;
      break;
    }
  }

  if(!found){
    fprintf(stderr,"Error: changing bond order of non-existent link\n");
    return false; 
  }

  unsigned int current_order = parent->orders[i];

  if(current_order == 4){
    fprintf(stderr,"Error: trying to increase order of aromatic bond - not possible\n");
    return false;
  }

  if ((child->num_edges + 1) > child->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->get_ch(),child->num_edges+1, child->allowed_edges);
    return false;
  }
  
  if ((parent->num_edges + 1) > parent->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->get_ch(),parent->num_edges+1, parent->allowed_edges);
    return false;
  }

  child->num_edges += 1;
  parent->num_edges += 1;
  parent->orders[i] += 1;

  return true; 
}



bool make_aromatic(WLNSymbol *child, WLNSymbol *parent, bool strict = true){

  if(!child || !parent){
    fprintf(stderr,"Error: attempting aromatic change of non-existent symbols, check locant range for rings or graph position for standard\n");
    return false;
  }

  if(child == parent){
    fprintf(stderr,"Warning: aromaticity change for equal WLNSymbol pointers\n");
    return true;
  }

  bool found = false;
  WLNSymbol *local_child = 0; 
  unsigned int i=0;
  for (i=0; i<parent->children.size();i++){
    local_child = parent->children[i];
    if(local_child == child){
      found = true;
      break;
    }
  }

  if(!found){
    if(strict)
      fprintf(stderr,"Error: attempting aromaticity of non-existent link\n");

    return false; 
  }


  unsigned int current_order = parent->orders[i];
  if (parent->orders[i] == 4)
    return true; // save some work

  // set aromatics
  switch(parent->get_ch()){

    case 'Y':
    case 'O':
    case 'M':
      parent->set_edges(2);
      break;
    
    case 'X':
    case 'C':
    case 'N':
      parent->set_edges(3);
      break;

    case 'K':
    case 'P':
    case 'S':
      parent->set_edges(4);
      break;

    case '*':
      parent->set_edges(8);
      break;
    
    default:
      fprintf(stderr,"Error: can not make %c symbol aromatic, please check definitions\n",parent->get_ch());
      return false;
  }

  // set aromatics
  switch(child->get_ch()){
    
    case 'Y':
    case 'O':
    case 'M':
      child->set_edges(2);
      break;

    case 'X':
    case 'C':
    case 'N':
      child->set_edges(3);
      break;

    case 'K':
    case 'P':
    case 'S':
      child->set_edges(4);
      break;

    case '*':
      child->set_edges(8);
      break; 
    
    default:
      fprintf(stderr,"Error: can not make %c symbol aromatic, please check definitions\n",child->get_ch());
      return false;
  }

  
  if (child->num_edges > child->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->get_ch(),child->num_edges, child->allowed_edges);
    return false;
  }

  if (parent->num_edges > parent->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->get_ch(),parent->num_edges, parent->allowed_edges);
    return false;
  }

  switch(current_order){
    case 1:
      break;
    case 2:
      // drop by 1
      child->num_edges += -1;
      parent->num_edges += -1;
      break; // equivilent order, do nothing
    case 3:
      // drop by 2
      child->num_edges += -2;
      parent->num_edges += -2;
      break;
    // already aromatic, do nothing
    case 4:
      break;
    
    default:
      fprintf(stderr,"Error: changing bond order of unknown bond type - %d\n",current_order);
  }


  parent->orders[i] = 4;
  return true; 
}


WLNSymbol* make_methyl(){

  WLNSymbol *carbon = AllocateWLNSymbol('C');
  carbon->set_edges(1);
  return carbon; 
}


/* resolve carbon methyl assumptions */
bool resolve_methyls(WLNSymbol *target){

  switch(target->get_ch()){

    case 'Y':
    case 'X':
    case 'K':
      while(target->num_edges < target->allowed_edges){
        WLNSymbol *methyl_head = make_methyl();
        link_symbols(methyl_head,target,1);
      }
      target->num_edges = target->allowed_edges;
      break;

    default:
      fprintf(stderr,"Error: resolving methyls performed on invalid symbol: %c\n",target->get_ch());
      return false;
  }

  return true;
}



WLNSymbol* define_element(std::string special){
    

  WLNSymbol *created_wln = AllocateWLNSymbol('*');
  created_wln->allowed_edges = INF; // allow anything for now;

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

  void print_distance(unsigned int *distance, unsigned int n){
    for (unsigned int i=0;i<n;i++){
      for (unsigned int j=0;j<n;j++){
        fprintf(stderr,"%d ",distance[i* n+j]); 
      }
      fprintf(stderr,"\n");
    }
  }

  /* creates distance matrix from connected grap */
  unsigned int *distance_matrix(unsigned int n){

    // distance is nxn matrix of a n vertices connected graph

    unsigned int *distance = (unsigned int*)malloc((n*n) *sizeof(unsigned int));

    // set the diagonal elements to zero dij = dii, other weights here are default 1
    for (unsigned int i=0;i<n;i++){
      for (unsigned int j=0;j<n;j++){
        if(i == j)
          distance[i* n+j] = 0; 
        else  
          distance[i* n+j] = INF; // n+1 for infinite 
      }
    }

    /*
    all other elements are the minimum 'steps' from one node to the next
    use Floyd-warshall to calculate distances, any 'infinites' here would indiciate a broken graph
    */

    // set the distance 1 pairs
    std::map<unsigned char,WLNSymbol*>::iterator map_iter;
    for (map_iter = locants.begin(); map_iter != locants.end(); map_iter++){
      WLNSymbol *current = map_iter->second; 
      unsigned int cur_int = locant_to_int(map_iter->first) - 1; // gives zero index for distance matrix
      for (WLNSymbol *child : current->children){
        unsigned int child_int = locant_to_int(locants_ch[child]) - 1;
        distance[cur_int* n+child_int] = 1;
        distance[child_int* n+cur_int] = 1;
      }
    } 

    for (unsigned int k=0;k<n;k++){
      for (unsigned int i=0;i<n;i++){
        for (unsigned int j=0;j<n;j++){
          if(distance[i* n+j] > distance[i* n+k] + distance[k* n+j])
            distance[i* n+j] = distance[i* n+k] + distance[k* n+j];
        }
      }
    }
    
    return distance; 
  }


  bool CreateMono(unsigned int local_size, bool aromatic){

    WLNSymbol *head = 0; 
    WLNSymbol *prev = 0;
    WLNSymbol *current = 0; 

    bool state = true;

    // assume already assigned locants
    for (unsigned int i=1;i<=local_size;i++){

      unsigned char loc = int_to_locant(i);

      if(!locants[loc]){
        current = assign_locant(loc,'C');
        current->set_edges(4);
      }
      else
        current = locants[loc];

      current->type = RING;

      if(aromatic)
        current->allowed_edges--; // take off 1 due to aromaticity.

      if (!head)
        head = current; 

      if(prev){
        if(!aromatic)
          state = link_symbols(current,prev,1);
        else
          state = link_symbols(current,prev,1,true);
      }
        
      prev = current;
    }

    if(!aromatic)
      state = link_symbols(head,prev,1);
    else
      state = link_symbols(head,prev,1,true);

    return state; 
  }

  /* creates poly rings, aromaticity is defined in reverse due to the nature of notation build */
  bool CreatePOLY(std::vector<std::pair<unsigned int,unsigned char>> &ring_assignments, std::vector<bool> &aromaticity){
     
    // perform the aromatic denotion check
    if (ring_assignments.size() != aromaticity.size()){
      fprintf(stderr,"Error: mismatch between number of rings and aromatic assignments\n");
      return false; 
    }

    unsigned int local_size = 0; 
    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i]; 

      if(local_size)
        local_size += component.first - 2;
      else
        local_size = component.first;
    }

    // create all the nodes in a large straight chain
    
    WLNSymbol *current = 0; 
    WLNSymbol *prev = 0; 
    for (unsigned int i=1;i<=local_size;i++){
      unsigned char loc = int_to_locant(i);
      if(!locants[loc]){
        current = assign_locant(loc,'C');
        current->set_edges(4);
      }
      else
        current = locants[loc];

      if(prev){
        if(!link_symbols(current,prev,1)){
          fprintf(stderr, "Error: inter-ring creating and bonding failed\n");
          return false; 
        }
      }
      prev = current;
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

      std::vector<unsigned char> ring_path;
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

          for (WLNSymbol *child : path->children){
            unsigned char child_loc = locants_ch[child];
            if(child_loc > highest_loc)
              highest_loc = child_loc;
          }    
          path = locants[highest_loc];
        }

        ring_path.push_back(locants_ch[path]); // add the last symbol
        bind_2 = highest_loc;
      }

      if(opt_debug)
        fprintf(stderr,"  fusing: %c <-- %c\n",bind_2,bind_1);
    
      if(!link_symbols(locants[bind_2],locants[bind_1],1)){
        fprintf(stderr,"Error: error in bonding locants together, check ring notation\n");
        return false;
      }

      if(aromatic){
        if(!AssignAromatics(ring_path))
          return false;
      }
        
      fuses++;
    }

    return true; 
  }


  /* interesting here that the multicyclic points are not used, i should perform a check for them though */
  bool CreatePERI(std::vector<std::pair<unsigned int,unsigned char>> &ring_assignments, 
                  std::vector<bool> &aromaticity,
                  std::vector<unsigned char> &multicyclic_locants,
                  unsigned char size_designator)
  {

    // perform the aromatic denotion check
    if (ring_assignments.size() != aromaticity.size()){
      fprintf(stderr,"Error: mismatch between number of rings and aromatic assignments\n");
      return false; 
    }

    // create a chain size of ring designator
    unsigned int local_size = locant_to_int(size_designator);
  
    // create all the nodes in a large straight chain
    WLNSymbol *current = 0; 
    WLNSymbol *prev = 0; 
    for (unsigned int i=1;i<=local_size;i++){
      unsigned char loc = int_to_locant(i);
      if(!locants[loc]){
        current = assign_locant(loc,'C');
        current->set_edges(4);
      }
      else
        current = locants[loc];

      if(prev){
        if(!link_symbols(current,prev,1)){
          fprintf(stderr, "Error: inter-ring creating and bonding failed\n");
          return false; 
        }
      }
      prev = current;
    }

    // calculate bindings and then traversals round the loops
    unsigned int comp_size = 0;
    unsigned char bind_1 = '\0';
    unsigned char bind_2 = '\0';
    unsigned int fuses = 0; 
    bool aromatic = false;

    unsigned int eval = 0; // for multicyclic
    unsigned int consumed = 0; 
   
    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i];
      comp_size = component.first;
      bind_1 = component.second;

      aromatic = aromaticity[i];
      WLNSymbol *path = locants[bind_1];

      std::vector<unsigned char> ring_path;

      // first pair can be calculated directly without a path travel
      if(!fuses){
        bind_2 = bind_1 + comp_size - 1; // includes start atom
        for (unsigned int i=0; i<comp_size;i++)
          ring_path.push_back(bind_1+i);
          
        consumed += comp_size;
      }
      else if (local_size - consumed < comp_size){ // ring wrap, which also needs to happen below
        
        // travel as much as we can, then add locants from bind_1 to get the binding position
        // from the size designator

        unsigned char highest_loc = '\0';  
        unsigned int track = 1; // include start
        while(!path->children.empty()){  
          track++;
          ring_path.push_back(locants_ch[path]);
          for (WLNSymbol *child : path->children){
            unsigned char child_loc = locants_ch[child];
            if(child_loc > highest_loc)
              highest_loc = child_loc;
          }    
          path = locants[highest_loc];
        }
        ring_path.push_back(locants_ch[path]); // add the last symbol


        // from the start locant, we move bind_1 diff spaces
        unsigned int diff = comp_size - track;

        bind_1 += diff; 
        bind_2 = size_designator; // created closure

        ring_path.push_back(bind_1);
      }
      else{
        
        if(path->num_edges == 3){ // standard multicylic points define a 3 ring share position, branching specials are handled differently

          unsigned char highest_loc = '\0';  
          for (unsigned int i=0;i<comp_size - 2; i++){  // 2 points are already defined
            ring_path.push_back(locants_ch[path]);

            for (WLNSymbol *child : path->children){
              unsigned char child_loc = locants_ch[child];
              if(child_loc > highest_loc)
                highest_loc = child_loc;
            }    
            path = locants[highest_loc];
          }
          ring_path.push_back(locants_ch[path]); // add the last symbol

          bind_2 = highest_loc;
          
          
          bind_1 += 1; // move to the next multicyclic position
          eval++;

          ring_path.push_back(bind_1); // add the last symbol
          consumed += comp_size - 3;
        }

        else{
          unsigned char highest_loc = '\0';
          for (unsigned int i=0;i<comp_size - 1; i++){
            ring_path.push_back(locants_ch[path]);

            for (WLNSymbol *child : path->children){
              unsigned char child_loc = locants_ch[child];
              if(child_loc > highest_loc)
                highest_loc = child_loc;
            }    
            path = locants[highest_loc];
          }

          ring_path.push_back(locants_ch[path]); // add the last symbol
          bind_2 = highest_loc;
          consumed += comp_size - 2;
        }          
      }

      if(opt_debug)
        fprintf(stderr,"  fusing: %c <-- %c\n",bind_2,bind_1);
      
      if(!link_symbols(locants[bind_2],locants[bind_1],1)){
        fprintf(stderr,"Error: error in bonding locants together, check ring notation\n");
        return false;
      }

      if(aromatic){
        if(!AssignAromatics(ring_path))
          return false;
      }

      fuses++;
    }


    return true; 
  }


  bool CreatePSDBRIDGE(std::vector<unsigned char> &fuses,
                       std::vector<unsigned int> &numerics,
                       unsigned int size)
  {
    unsigned int i = 1; 
    while (i < fuses.size()){
      if(opt_debug)
        fprintf(stderr,"  fusing: %c - %c\n",fuses[i-1],fuses[i]);
      i+=2;
    }
    return true; 
  }

  bool AssignAromatics(std::vector<unsigned char> &ring_path){
    
    return true;

    for (unsigned int i=1; i<ring_path.size();i+=1){

      unsigned char par = ring_path[i-1];
      unsigned char chi = ring_path[i];

      if(!make_aromatic(locants[chi],locants[par],false)){

        // try with the alternative
        if(!make_aromatic(locants[par],locants[chi])){
          fprintf(stderr,"Error: error in changing aromaticity - check ring notation [%c --> %c]\n",par,chi);
          return false;
        }
      }
    }

    if(!make_aromatic(locants[ring_path.back()],locants[ring_path.front()],false)){
      if(!make_aromatic(locants[ring_path.front()],locants[ring_path.back()])){
        fprintf(stderr,"Error: error in changing aromaticity - check ring notation [%c --> %c]\n",ring_path.front(),ring_path.back());
        return false; 
      }
    }


    return true;
  }

  unsigned char modifiy_locant(unsigned char positional_locant,unsigned int size_modifier){
    positional_locant += (size_modifier * 23);
    return positional_locant;
  }


  void FormWLNRing(std::string &block, unsigned int start){

    enum RingType{ MONO=0, POLY=1, PERI=2, BRIDGED=3, PSDBRIDGED = 4}; 
    const char* ring_strings[] = {"MONO","POLY","PERI","BRIDGED","PSDBRIDGED"};

    unsigned int ring_type = MONO;   // start in mono and climb up
    unsigned int end = 0;

    bool warned             = false;  // limit warning messages to console
    bool heterocyclic       = false;  // L|T designator can throw warnings


    // -- paths -- // int allows way more description in states
    unsigned int pending_multi        = 0; // 0 - closed, 1 - open multi notation, 2 - expect size denotation
    unsigned int pending_pseudo       = 0; 
    unsigned int pending_bridge       = 0;
    unsigned int pending_aromatics    = 0;
    unsigned int pending_special      = 0;
    
    unsigned int expected_locants     = 0;
    unsigned int size_modifier        = 0;       // multiple of 23 to move along for locant

    unsigned int  size_set              = 0; 
    
    unsigned char ring_size_specifier   = '\0';
    unsigned char positional_locant     = '\0'; 

    // allows stoi use
    std::string expanded_size; 
    std::string special; 

    std::vector<bool> aromaticity; 
    std::vector<std::pair<unsigned char, unsigned char>>  bond_increases; 

    std::vector<unsigned char> fuses; // read as pairs
    std::vector<unsigned char> bridge_locants;
    std::vector<unsigned char> multicyclic_locants;
    
    std::vector<std::pair<unsigned int, unsigned char>>  ring_components; 
   
    for (unsigned int i=0;i<block.size();i++){
      unsigned char ch = block[i];

      switch(ch){
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
          if(size_modifier && positional_locant){
            modifiy_locant(positional_locant,size_modifier);
            size_modifier = 0;
          }

          if (i > 1 && block[i-1] == ' '){
            pending_multi   = 1;
            expected_locants = ch - '0';
            break;
          }    
          else{
            if(positional_locant)
              ring_components.push_back({ch-'0',positional_locant});
            else
              ring_components.push_back({ch-'0','A'});

            positional_locant = '\0';
            break;
          }
          
            
        case '/':
          if(size_modifier && positional_locant){
            modifiy_locant(positional_locant,size_modifier);
            size_modifier = 0;
          }

          if(pending_special){
            fprintf(stderr,"Error: character %c is not allowed in '-<A><A>-' format where A is an uppercase letter\n",ch);
            Fatal(start+i);
          }

          expected_locants = 2; 
          pending_pseudo = true;
          ring_type = PSDBRIDGED; 
          break; 

        case '-':
          if(size_modifier && positional_locant){
            modifiy_locant(positional_locant,size_modifier);
            size_modifier = 0;
          }

          if (positional_locant){
            // opens up inter ring special definition
            if(!pending_special){
              pending_bridge = false;
              //expecting_component = 0;

              pending_special = true; 
            }
            else{
              locants[positional_locant] = define_element(special);
              locants[positional_locant]->type = RING;
              special.clear();
            }

          }
          break;

        // aromatics and locant expansion
        case '&':
          if(positional_locant)
            size_modifier++;

          break;

        case ' ':

          if(expected_locants){
            fprintf(stderr,"Error: %d more locants expected before space seperator\n",expected_locants);
            Fatal(start+i);
          }


          if(!size_set && ring_size_specifier)
            size_set = true;

          // resets any pendings and set states
          if(pending_multi){
            pending_multi     = false;
            if(ring_type < PERI)
              ring_type = PERI; 
          }
          else if (pending_bridge){
            if(ring_type < BRIDGED && positional_locant)
              ring_type = BRIDGED; 

            bridge_locants.push_back(positional_locant);
            pending_bridge = false; 
          }
          
          pending_pseudo    = false;
          //expecting_component = 0;
          positional_locant = '\0'; // hard reset the positional locant
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
          if(size_modifier && positional_locant){
            modifiy_locant(positional_locant,size_modifier);
            size_modifier = 0;
          }

          if(expected_locants){
            if(pending_multi){
              multicyclic_locants.push_back(ch);
              expected_locants--;
            }
            else if (pending_pseudo){
              fuses.push_back(ch);
              expected_locants--; 
            }
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }
            break;
          }

          else if(block[i-1] == ' '){
            // if(completed_multi && !ring_size_specifier){ // a size specifier is always needed
            //   ring_size_specifier = ch;
            //   size_set = 1; 
            // }
            positional_locant = ch;
            break;
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
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edges(5);
                positional_locant++; // allows inline defition continuation
                break;

              case 'Y':
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edges(3);
                positional_locant++; // allows inline defition continuation
                break;
              case 'N':
                if(!heterocyclic)
                  warned = true;
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edges(3);
                positional_locant++; // allows inline defition continuation
                break;

              case 'V':
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edges(2);
                positional_locant++; // allows inline defition continuation
                break;

              case 'M':
              case 'O':
                if(!heterocyclic)
                  warned = true;
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edges(2);
                positional_locant++; // allows inline defition continuation
                break;

              case 'X':
                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edges(4);
                positional_locant++; // allows inline defition continuation
                break;
              case 'K':
                if(!heterocyclic)
                  warned = true;

                new_locant = assign_locant(positional_locant,ch);
                new_locant->set_edges(4);
                positional_locant++; // allows inline defition continuation
                break;

              case 'U':
                if(opt_debug)
                  fprintf(stderr,"  increasing bond order from %c to %c by 1\n",positional_locant,positional_locant+1);

                bond_increases.push_back({positional_locant,positional_locant+1});
                break;

              
              default:
                fprintf(stderr,"Error: %c is not allowed as a atom assignment within ring notation\n",ch);
                Fatal(start+i);
            }            
            break;
          }
          else{
            positional_locant = ch;
            break;
          }
          

        // openers 
        case 'L':
          if(size_modifier && positional_locant){
            modifiy_locant(positional_locant,size_modifier);
            size_modifier = 0;
          }


          if(i==0){
           heterocyclic = false; 
           break;
          }
         
          else if(expected_locants){

            if(pending_multi){
              multicyclic_locants.push_back(ch);
              expected_locants--;
            }
            else if (pending_pseudo){
              fuses.push_back(ch);
              expected_locants--; 
            }
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }
            break;
          }
          else if (pending_special){
            special.push_back(ch);
            if(special.size() > 2){
              fprintf(stderr,"Error: special elemental notation must only have two characters\n");
              Fatal(start+i);
            } 
            break;
          }
          else if(block[i-1] == ' '){

            if(!ring_size_specifier){ // a size specifier is always needed
              ring_size_specifier = ch;
              size_set = 1; 
              positional_locant = ch;
            }
            else{
              positional_locant = ch;
              pending_bridge = true;
            }
            break;
          }
          else{
            positional_locant = ch;
            break;
          }

        case 'T':
          if(size_modifier && positional_locant){
            modifiy_locant(positional_locant,size_modifier);
            size_modifier = 0;
          }
          
          if(i==0){
            heterocyclic = true; 
            break; 
          }


          else if (pending_special){
            special.push_back(ch);
            if(special.size() > 2){
              fprintf(stderr,"Error: special elemental notation must only have two characters\n");
              Fatal(start+i);
            } 
            break;
          }
          if(expected_locants){

            if(pending_multi){
              multicyclic_locants.push_back(ch);
              expected_locants--;
            }
            else if (pending_pseudo){
              fuses.push_back(ch);
              expected_locants--; 
            }
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }
            break;
          }
          else if (pending_aromatics){
            aromaticity.push_back(false); // simple here
            break;
          }
          else if (positional_locant || std::isdigit(block[i-1])){
            pending_aromatics = true;
            aromaticity.push_back(false);
            positional_locant = ch;
            break;
          }
          else if(block[i-1] == ' '){

            if(!ring_size_specifier){ // a size specifier is always needed
              ring_size_specifier = ch;
              size_set = 1;
              positional_locant = ch;
            }
            else{
              positional_locant = ch;
              pending_bridge = true;
            }
            break;
          }
          else{
            positional_locant = ch;
            break;
          }

        //closure
        case 'J':
          end = i;
  
          if(expected_locants){
            if(pending_multi){
              multicyclic_locants.push_back(ch);
              expected_locants--;
            }
            else if (pending_pseudo){
              fuses.push_back(ch);
              expected_locants--; 
            }
            else{
              fprintf(stderr,"Error: unhandled locant rule\n");
              Fatal(start+i);
            }
            break;
          }
          else if (pending_special){
            special.push_back(ch);
            if(special.size() > 2){
              fprintf(stderr,"Error: special elemental notation must only have two characters\n");
              Fatal(start+i);
            } 
            break;
          }
          else if(block[i-1] == ' '){

            if(!ring_size_specifier){
              ring_size_specifier = ch; 
              size_set = 1;
              positional_locant = ch;
            }
            else{
              pending_bridge = true;
              positional_locant = ch;
            }
            break;
          }
          else{
            positional_locant = ch;
            break;
          }
            
        default:
          fprintf(stderr,"Error: unrecognised symbol in ring definition: %c\n",ch);
          Fatal(start + i);
      }
       
    }

    // set the ring type if not handled in parser
    if (ring_components.size() > 1 && ring_type < PERI)
      ring_type = POLY;
    

    // shorthand aromatic conditions
    if (aromaticity.size() == 1 && aromaticity[0] == false){
      while(aromaticity.size() < ring_components.size())
        aromaticity.push_back(false);
    }
    else if (aromaticity.empty()){
      while(aromaticity.size() < ring_components.size())
        aromaticity.push_back(true);
    }

    // reverse the aromaticity assignments, how the notation works
    std::reverse(aromaticity.begin(), aromaticity.end());


    if(warned)
      fprintf(stderr,"Warning: heterocyclic ring notation required for inter atom assignment, change starting 'L' to 'T'\n");
    
    if(ring_components.empty()){
      fprintf(stderr,"Error: error in reading ring components, check numerals in ring notation\n");
      Fatal(start+end);
    }

    // debug here
    if (opt_debug){
      
      fprintf(stderr,"  ring type: %s",ring_strings[ring_type]);

      fprintf(stderr,"  ring components: ");
      for (std::pair<unsigned int, unsigned char> comp : ring_components)
        fprintf(stderr,"%d(%c) ",comp.first,comp.second);
      fprintf(stderr,"\n");

      fprintf(stderr,"  aromaticity: ");
      for (bool aromatic : aromaticity)
        fprintf(stderr,"%d ",aromatic);
      fprintf(stderr,"\n");

      fprintf(stderr,"  multicyclic points: ");
      for (unsigned char loc : multicyclic_locants){
        fprintf(stderr,"%c ",loc == ' ' ? '_':loc);
      }
      fprintf(stderr,"\n");

      fprintf(stderr,"  bridge points: ");
      for (unsigned char loc : bridge_locants){
        fprintf(stderr,"%c ",loc == ' ' ? '_':loc);
      }
      fprintf(stderr,"\n");

      fprintf(stderr,"  hard fuses: ");
      for (unsigned int i=1;i<fuses.size();i+=2)
        fprintf(stderr,"(%c --> %c) ",fuses[i-1],fuses[i]);
      fprintf(stderr,"\n");

      
      fprintf(stderr,"  size denotion: %d\n",ring_size_specifier ? locant_to_int(ring_size_specifier) : 0);
      fprintf(stderr,"  heterocyclic: %s\n", heterocyclic ? "yes":"no");
    }

    
    bool state = true;
    switch(ring_type){
      case MONO:
        state = CreateMono(ring_components[0].first,aromaticity[0]);
        break;
      case POLY:
        state = CreatePOLY(ring_components,aromaticity);
        break;
      case PERI:
        state = CreatePERI(ring_components,aromaticity,multicyclic_locants,ring_size_specifier);
        break;
      case BRIDGED:
      case PSDBRIDGED:
        break;
    }

    if (!state)
      Fatal(start+end);

    for (std::pair<unsigned char, unsigned char> bond_pair : bond_increases)
      increase_bond_order(locants[bond_pair.second],locants[bond_pair.first]);
    
  }


};

WLNRing *AllocateWLNRing()
{
  WLNRing *wln_ring = new WLNRing;
  ring_mempool.push_back(wln_ring);
  return wln_ring;
}

// expensive but sometimes necessary for edge cases
void DeallocateWLNRing(WLNRing *ring)
{
  // find the ring in the mem pool
  unsigned int i = 0;
  for (WLNRing *r : ring_mempool)
  {
    if (r == ring)
      break;
    i++;
  }

  ring_mempool.erase(ring_mempool.begin() + i);
  delete ring;
}

struct WLNGraph
{
  
  WLNSymbol *root;

  WLNGraph() : root{(WLNSymbol *)0} {};
  ~WLNGraph()
  {
    for (WLNSymbol *node : symbol_mempool)
      delete node;
    for (WLNRing *ring : ring_mempool)
      delete ring;
  };

  bool expand_carbon_chain(WLNSymbol *head,unsigned int size){

    if (size > REASONABLE)
      fprintf(stderr,"Warning: making carbon chain over 1024 long, reasonable molecule?\n");
          
    head->value = 'C';
    head->set_edges(4);
    head->num_edges = 0;

    WLNSymbol *tmp = 0;
    unsigned int tmp_order = 0;
    // hold the bonds

    if(!head->children.empty()){
      tmp = head->children[0];
      tmp_order = head->orders[0];
      head->children.clear();
      head->orders.clear();
    }
          
    WLNSymbol *prev = head;
    for(unsigned int i=0;i<size-1;i++){
      WLNSymbol* carbon = AllocateWLNSymbol('C');
      carbon->set_edges(4); // allows hydrogen resolve
      link_symbols(carbon,prev,1);
      prev = carbon;
    } 

    if(tmp){
      prev->children.push_back(tmp);
      prev->orders.push_back(tmp_order);
    }

    return true;
  }

  /* must be performed before sending to obabel graph*/
  bool ExpandWLNGraph(){

 
    unsigned int stop = symbol_mempool.size();
    for (unsigned int i=0;i<stop;i++){
      WLNSymbol *sym = symbol_mempool[i];

      switch(sym->get_ch()){

        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if (!sym->special.empty())
            expand_carbon_chain(sym,std::stoi(sym->special));
          else
            expand_carbon_chain(sym,sym->get_ch() - '0');
          break;
        
        case 'Y':
        case 'X':
        case 'K':
          resolve_methyls(sym);
          break;

        default:
          break; // ignore
      }

    }

    Reindex_lookups();
    return true; 
  }


  WLNSymbol *return_open_branch(std::stack<WLNSymbol *> &branch_stack)
  {

    if (branch_stack.empty())
      return (WLNSymbol *)0;

    WLNSymbol *top = branch_stack.top();

    // doesnt pop here, only '&' can pop the stack

    return top;
  }

  bool check_unbroken(unsigned int i)
  {
    if (i > 1 && !(wln[i - 1] == '&' && wln[i - 2] == ' '))
    {
      fprintf(stderr, "Error: broken graph without ionic notation, check branches|locants and '&' count\n");
      return false;
    }

    return true;
  }

  WLNRing *pop_ringstack(unsigned int pops, std::stack<WLNRing *> &stack)
  {

    if (pops >= stack.size())
    {
      fprintf(stderr, "Error: trying to pop too many rings check '&' count\n");
      return (WLNRing *)0;
    }

    for (unsigned int i = 0; i < pops; i++)
      stack.pop();

    return stack.top();
  }

  // this has a return clause in it and needs previous
  WLNSymbol *pop_branchstack(unsigned int pops, std::stack<WLNSymbol *> &stack, WLNSymbol *prev)
  {

    if (!prev)
      fprintf(stderr, "Error: popping with no previous symbol\n");

    bool hard = false;

    if (prev == stack.top())
      hard = true;

    if (opt_debug)
      fprintf(stderr, "  popping %d symbols down the stack: mode(%d) prev[%c]\n", pops, hard, prev->get_ch());

    if (hard)
    {
      if (pops >= stack.size())
      {
        fprintf(stderr, "Error: to many stack pops - check '&' count\n");
        return 0;
      }
      for (unsigned int i = 0; i < pops; i++)
        stack.pop();
    }
    else
    {
      if (pops > stack.size())
      {
        fprintf(stderr, "Error: to many stack pops - check '&' count\n");
        return 0;
      }
      for (unsigned int i = 1; i < pops; i++)
        stack.pop();
    }
    return stack.top();
  }

  /* wraps popping for the linker and branch stacks */
  WLNSymbol *pop_standard_stacks(unsigned int pop_ticks,
                                 std::stack<WLNSymbol *> &branch_stack,
                                 std::stack<WLNSymbol *> &linker_stack,
                                 WLNSymbol *prev, unsigned int i)
  {
    WLNSymbol *ret = 0;
    if (!branch_stack.empty())
      ret = pop_branchstack(pop_ticks, branch_stack, prev);
    else if (!linker_stack.empty())
      ret = pop_branchstack(pop_ticks, linker_stack, prev);
    else
    {
      fprintf(stderr, "Error: popping empty stacks - check '&' count\n");
      Fatal(i);
    }

    return ret;
  }

  /* wraps the linking and graph checking functions */
  void create_bond(WLNSymbol *curr, WLNSymbol *prev,
                   unsigned int bond_ticks, unsigned int i)
  {
    if (prev)
    {
      if (!link_symbols(curr, prev, 1 + bond_ticks))
        Fatal(i);
    }
    else
    {
      if (!check_unbroken(i))
        Fatal(i);
    }
  }

  /*  Wraps the creation of locant and bonding back ring assignment */
  void create_locant(WLNSymbol *curr, std::stack<WLNRing *> &ring_stack, unsigned int i)
  {
    
    curr->type = LOCANT;

    WLNRing *s_ring = 0;
    unsigned char ch = wln[i];

    if (ring_stack.empty())
    {
      fprintf(stderr, "Error: no rings to assign locants to\n");
      Fatal(i);
    }
    else
      s_ring = ring_stack.top();

    if (s_ring->locants[ch]){
      if(!link_symbols(curr,s_ring->locants[ch],1))
        Fatal(i);
    } 
    else {
      fprintf(stderr, "Error: assigning locant outside of ring\n");
      Fatal(i);
    }

  }

  /* a global segmentation using both rule sets - start merging */
  bool ParseWLNString()
  {
    
    if (opt_debug)
      fprintf(stderr, "Parsing WLN notation:\n");

    unsigned int len = strlen(wln);

    std::stack<WLNRing *> ring_stack;   // access through symbol
    std::stack<WLNSymbol *> branch_stack; // between locants, clean branch stack
    std::stack<WLNSymbol *> linker_stack; // used for branching ring systems

    WLNSymbol *curr = 0;
    WLNSymbol *prev = 0;
    WLNRing   *ring = 0;

    bool pending_locant = false;
    bool pending_special = false;
    bool pending_closure = false;
    bool pending_inline_ring = false;
    bool pending_spiro = false;

    std::string special;

    // allows consumption of notation after block parses
    unsigned int block_start = 0;
    unsigned int block_end = 0;

    unsigned int pop_ticks = 0;  // '&' style popping
    unsigned int bond_ticks = 0; // 'U' style bonding

    for (unsigned int i = 0; i < len; i++)
    {
      unsigned char ch = wln[i];

      switch (ch)
      {

      case '0': // cannot be lone, must be an addition to another num
        if (pending_closure)
        {
          break;
        }
        if(prev && std::isdigit(prev->get_ch()))
          prev->special.push_back(ch);
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
        if (pending_closure)
          break;
        else if(pending_special){
          fprintf(stderr,"Error: character %c in special elemental definition are not allowed\n",ch);
          Fatal(i);
        }
        
        if (pop_ticks)
        {
          prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
          pop_ticks = 0;
        }

        if(!std::isdigit(prev->get_ch())){
          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(4);

          curr->special.push_back(ch); // prepare for a n > 10 chain

          create_bond(curr, prev, bond_ticks, i);
          prev = curr;
        }
        else 
          prev->special.push_back(ch);
        break;

      case 'Y':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {

          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(3);

          branch_stack.push(curr);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

      case 'X':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {

          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(4);

          branch_stack.push(curr);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

        // oxygens

      case 'O':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {

          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(2);

          branch_stack.push(curr);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

      case 'Q':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else if (pending_closure || pending_special)
        {
          continue;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(1);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = return_open_branch(branch_stack);
        }
        break;

      case 'V':
      case 'W':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {
          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(2);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

        // nitrogens

      case 'N':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(3);

          branch_stack.push(curr);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

      case 'M':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(2);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

      case 'K':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(4);

          branch_stack.push(curr);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

      case 'Z':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else if (pending_closure || pending_special)
        {
          continue;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(1);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = return_open_branch(branch_stack);
        }
        break;

        // halogens - need to add rules for semi allowed hyper valence in ionions

      case 'E':
      case 'G':
      case 'F':
      case 'I':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else if (pending_closure || pending_special)
        {
          continue;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(1);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;

          prev = return_open_branch(branch_stack);
        }
        break;

        // inorganics

      case 'B':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else if (pending_closure || pending_special)
        {
          continue;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(3);

          branch_stack.push(curr);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

      case 'P':
      case 'S':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else if (pending_closure || pending_special)
        {
          continue;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(6);

          branch_stack.push(curr);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
        }
        break;

        // locants only?

      case 'A':
      case 'C':
      case 'D':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {

          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

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
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {

          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else{
          // explicit hydrogens
          curr = AllocateWLNSymbol(ch);
          curr->set_type(STANDARD);
          curr->set_edges(1);

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          
          // this will automatically return a branch like a terminator
          prev = return_open_branch(branch_stack);
        }
        break;

        // ring notation

      case 'J':
        if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }

        else if (pending_closure)
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
            if (ring->locants[prev->get_ch()])
              create_bond(ring->locants[prev->get_ch()],prev,bond_ticks,i);
            else
            {
              fprintf(stderr, "Error: attaching inline ring with out of bounds locant assignment\n");
              Fatal(i);
            }
          }

          bond_ticks = 0;
          pending_closure = false;
        }
        else{
          fprintf(stderr,"Error: J notation used outside of a ring definition is undefined\n");
          Fatal(i);
        }
          

        break;

      case 'L':
      case 'T':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

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
          pending_closure = true;
        }
        break;

      case 'R':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          ring = AllocateWLNRing();

          std::string r_notation = "L6J";
          ring->FormWLNRing(r_notation,i);
          ring_stack.push(ring);

          if (prev)
            create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
        }
        break;

        // bonding

      case 'U':
        if(pending_closure)
          break;
        else if(pending_special){
          pending_inline_ring = false; // resets
          
          if(special.size() == 2){
            fprintf(stderr,"Error: special element definition must follow format '-<A><A>-' where A is an uppercase letter\n");
            Fatal(i);
          }
          else
            special.push_back(ch);
        }
        else if (pending_locant)
        {
          curr = AllocateWLNSymbol(ch);
          curr->set_type(LOCANT);
          curr->set_edges(2); // locants always have two edges

          if (pending_inline_ring)
            create_bond(curr, prev, bond_ticks, i);
          else
            create_locant(curr, ring_stack, i);

          prev = curr;
          pending_locant = false;
        }
        else
          bond_ticks++;

        break;

        // specials

      case ' ':
        if (pending_closure)
        {
          break;
        }
        if (pending_inline_ring)
        {
          pending_special = false; // these are set at the same time

          // send the linker into its own stack
          if (!branch_stack.empty() && (branch_stack.top()->num_edges < branch_stack.top()->allowed_edges))
            linker_stack.push(branch_stack.top());
        }

        // clear the branch stack betweek locants and ions
        while (!branch_stack.empty())
          branch_stack.pop();

        if (pop_ticks)
        {
          ring = pop_ringstack(pop_ticks, ring_stack);
          if (!prev)
            Fatal(i);
          pop_ticks = 0;
        }

        pending_locant = true;
        break;

      case '&':
        if (pending_closure)
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
        else
        {
          pop_ticks++; // set the number of pops to do
        }
        break;

      case '-':
        if (pending_closure)
          break;

        else if (pending_inline_ring)
        {
          fprintf(stderr, "Error: only one pending ring can be active, check closures\n");
          Fatal(i);
        }
        else if (pending_special)
        {

          if (pop_ticks)
          {
            prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
            pop_ticks = 0;
          }

          curr = define_element(special);
          special.clear();

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
          pending_special = false;
        }
        else{
          pending_inline_ring = true;
          pending_special = true;
        }
        break;



      case '/':
        if (pending_closure)
          break;
        else if (pending_special){
          fprintf(stderr,"Error: character %c in special elemental definition are not allowed\n",ch);
          Fatal(i);
        }
        prev = curr;
        curr = AllocateWLNSymbol(ch);
        break;

      default:
        fprintf(stderr, "Error: unallowed character! - [A-Z][0-1][&-/' ']\n");
        Fatal(i);
      }
    }

    if (pending_closure)
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


    Reindex_lookups();
    return true;
  }

  /* dump wln tree to a dotvis file */
  void WLNDumpToDot(FILE *fp)
  {  
    fprintf(fp, "digraph WLNdigraph {\n");
    fprintf(fp, "  rankdir = LR;\n");
    for (WLNSymbol *node : symbol_mempool)
    {
      fprintf(fp, "  %d", index_lookup[node]);
      if (node->get_ch() == '*')
        fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
      else if(node->type == LOCANT)
        fprintf(fp, "[shape=circle,label=\"%c\",color=blue];\n", node->get_ch());
      else if (node->type == RING)
        fprintf(fp, "[shape=circle,label=\"%c\",color=green];\n", node->get_ch());
      else{
        if(std::isdigit(node->get_ch())){
          if (!node->special.empty())
            fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
          else
            fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->get_ch());
        } 
        else
          fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->get_ch());
      }
        
      for(unsigned int i = 0; i<node->children.size(); i++){
        WLNSymbol *child = node->children[i];
        unsigned int bond_order = node->orders[i];
        
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


  OpenBabel::OBAtom* NMOBMolNewAtom(OpenBabel::OBMol* mol, unsigned int elem,unsigned int charge=0)
  {

    OpenBabel::OBAtom* result = mol->NewAtom();
    result->SetAtomicNum(elem);
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

  bool MBOBReduceSpecies(OpenBabel::OBMol* mol){

    OpenBabel::OBAtom* atom = 0;
    int reducing = 0; // reducing hydrogens count
    unsigned int valence = 0; 
    FOR_ATOMS_OF_MOL(atom,mol){
      // allows for wln notation implicits
      if(atom->GetImplicitHCount() > 0)
        continue;
      
      bool aromatic = false;
      if(atom->IsAromatic())
        aromatic = true;

      valence = atom->GetExplicitValence();

      switch(atom->GetAtomicNum()){
        
        case 6: // carbon
          reducing = 4;
          break;

        case 7:  // nitrogen
        case 15: // phosphorus
          reducing = 3;
          break;

        default:
          break;
      }

  
      if(aromatic)
        reducing--;

      if(reducing > 0){
        unsigned int defined_charge = atom->GetFormalCharge();
        if(defined_charge)
          reducing+=defined_charge;
        
        atom->SetImplicitHCount(reducing-valence);
      }

      
    }

    return true;
  }

  bool NMOBSanitizeMol(OpenBabel::OBMol* mol)
  {
    

    mol->SetAromaticPerceived(true);

    MBOBReduceSpecies(mol);

    if (!OBKekulize(mol))
        return false;
    mol->DeleteHydrogens();

    // WLN has no inherent stereochemistry, this can be a flag but should be off by default
    mol->SetChiralityPerceived(true);
    return true;
  }


  bool ConvertFromWLN(OpenBabel::OBMol* mol,WLNGraph &wln_graph){

    if(opt_debug)
      fprintf(stderr,"Converting wln to obabel mol object: \n");

    // set up atoms
    for (WLNSymbol *sym: symbol_mempool){

      if(sym->type != LOCANT){
        
        OpenBabel::OBAtom *atom = 0;

        unsigned int atomic_num = 0;
        unsigned int charge = 0; 
        unsigned int h_count = 0;

        switch(sym->get_ch()){

          case 'H':
            atomic_num = 1;
            h_count = 0;
            break; 

          case 'B':
            atomic_num = 5;
            break;

          case 'C':
          case 'X':
          case 'Y':
            atomic_num = 6; 
            break;

          case 'N':
            atomic_num = 7;
            break;

          case 'M':
            atomic_num = 7;
            h_count = 1;
            break;

          case 'Z':
            atomic_num = 7; 
            h_count = 2;
            break;

          case 'K':
            atomic_num = 7;
            charge = 1; 
            break;

          case 'O':
            atomic_num = 8;
            break;
          case 'Q':
            atomic_num = 8;
            h_count = 1;
            break;

          case 'F':
            atomic_num = 9;
            break;
          
          case 'P':
            atomic_num = 15;
            break;
          
          case 'S':
            atomic_num = 16;
            break;

          case 'G':
            atomic_num = 17;
            break;

          case 'E':
            atomic_num = 35;
            break;

          case 'I':
            atomic_num = 53;
            break;
        
          case '*':
            atomic_num = special_element_atm(sym->special);
            break;

          default:
            fprintf(stderr,"Error: unrecognised WLNSymbol* char in obabel mol build - %c\n",sym->get_ch());
            return false;
        }

        atom = NMOBMolNewAtom(mol,atomic_num,charge);
        atom->SetImplicitHCount(h_count);

        if(sym->type == RING)
          atom->SetInRing();

        babel_atom_lookup[index_lookup[sym]] = atom;
        if(opt_debug)
          fprintf(stderr,"  created: atom[%d] - atomic num(%d), charge(%d)\n",atom->GetIdx(),atomic_num,charge);
      }

    }

    // set bonds
    for(WLNSymbol *parent : symbol_mempool){

      if(parent->type == LOCANT)
        continue;

      unsigned int parent_id = index_lookup[parent];
      OpenBabel::OBAtom *par_atom = babel_atom_lookup[parent_id];

      for (unsigned int i=0;i<parent->children.size();i++){
        
        WLNSymbol *child = parent->children[i];
        unsigned int bond_order = 0; 
    
        // skip across locants
        if(child->type == LOCANT){
          bond_order = child->orders[0];
          child = child->children[0]; // skip to the start of the chain
        }
        else
          bond_order = parent->orders[i];
      
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

  if(!wln_graph.ParseWLNString()){
    fprintf(stderr,"Error: failed on wln graph formation\n");
    return false;
  }

  if(!wln_graph.ExpandWLNGraph()){
    fprintf(stderr,"Error: failed in expanding wln graph to SCT format\n");
    return false;
  }

  // create the wln dotfile
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


  if(!obabel.ConvertFromWLN(mol,wln_graph)){
    fprintf(stderr,"Error: failed on obabel mol object formation\n");
    return false;
  }

  if(!obabel.NMOBSanitizeMol(mol)){
    fprintf(stderr,"Error: failed on mol sanitize\n");
    return false; 
  }

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