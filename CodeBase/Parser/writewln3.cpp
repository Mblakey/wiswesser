

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
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

// --- macros ---
#define REASONABLE 1024
#define INF 9999

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
std::map<unsigned int, WLNSymbol *> symbol_lookup;
unsigned int glob_index = 0;

// --- pools ---
std::vector<WLNSymbol *> symbol_mempool;
std::vector<WLNRing *> ring_mempool;

enum WLNTYPE
{
  STANDARD = 0,
  LOCANT = 1,
  LINKER = 2,
  RING = 3,
  SPECIAL = 4
};

// rule 2 - hierarchy - rules have diverged due to end terminator char, also use for locant setting from 14
std::map<unsigned char, unsigned int> char_hierarchy =
    {
        {' ', 1}, {'-', 2}, {'/', 3}, {'0', 4}, {'1', 5}, {'2', 6}, {'3', 7}, {'4', 8}, {'5', 9}, {'6', 10}, {'7', 11}, {'8', 12}, {'9', 13}, {'A', 14}, {'B', 15}, {'C', 16}, {'D', 17}, {'E', 18}, {'F', 19}, {'G', 20}, {'H', 21}, {'I', 22}, {'J', 23}, {'K', 24}, {'L', 25}, {'M', 26}, {'N', 27}, {'O', 28}, {'P', 29}, {'Q', 30}, {'R', 31}, {'S', 32}, {'T', 33}, {'U', 34}, {'V', 35}, {'W', 36}, {'X', 37}, {'Y', 38}, {'Z', 40}, {'&', 41}};

std::map<unsigned char, unsigned int> locant_integer_map =
    {
        {'A', 1}, {'B', 2}, {'C', 3}, {'D', 4}, {'E', 5}, {'F', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'J', 10}, {'K', 11}, {'L', 12}, {'M', 13}, {'N', 14}, {'O', 15}, {'P', 16}, {'Q', 17}, {'R', 18}, {'S', 19}, {'T', 20}, {'U', 21}, {'V', 22}, {'W', 23}, {'X', 24}, {'Y', 25}, {'Z', 26}};

std::map<unsigned int, unsigned char> integer_locant_map =
    {
        {1, 'A'}, {2, 'B'}, {3, 'C'}, {4, 'D'}, {5, 'E'}, {6, 'F'}, {7, 'G'}, {8, 'H'}, {9, 'I'}, {10, 'J'}, {11, 'K'}, {12, 'L'}, {13, 'M'}, {14, 'N'}, {15, 'O'}, {16, 'P'}, {17, 'Q'}, {18, 'R'}, {19, 'S'}, {20, 'T'}, {21, 'U'}, {22, 'V'}, {23, 'W'}, {24, 'X'}, {25, 'Y'}, {26, 'Z'}};

// --- utilities ---

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
  glob_index = 0;
  for (WLNSymbol *node : symbol_mempool)
  {
    index_lookup[node] = glob_index;
    symbol_lookup[glob_index] = node;
    glob_index++;
  }
}

// should be all we need for a SCT XI connection table - obabel can handle coords
struct Atom
{
  std::string symbol;
  unsigned int atomic_num;

  int charge;

  std::vector<Atom> bonded;
  std::vector<unsigned int> orders;
};

struct AtomGraph
{
  Atom *head;
};

struct WLNSymbol
{

  unsigned char ch;
  unsigned int type;

  unsigned int allowed_edges;
  unsigned int num_edges;

  WLNSymbol *previous;
  std::vector<WLNSymbol *> children; // linked list of next terms chains
  std::vector<unsigned int> orders;

  // if 'ch='*'
  std::string special; // string for element, or ring
  
  // dont hold the ring anymore!

  // if default needed
  WLNSymbol()
  {
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
    previous = 0;
    special = "";
  }
  ~WLNSymbol(){};

  void set_edges(unsigned int edges)
  {
    allowed_edges = edges;
  }

  void set_type(unsigned int i)
  {
    type = i;
  }

  // resets the symbol
  void reset()
  {
    ch = '\0';
    allowed_edges = 0;
    num_edges = 0;
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
  wln->ch = ch;

  index_lookup[wln] = glob_index;
  symbol_lookup[glob_index] = wln;
  glob_index++;
  return wln;
}

// these are expensive, but needed for some edge case notations
void DeallocateWLNSymbol(WLNSymbol *node)
{

  if (opt_debug)
    fprintf(stderr, "  manual deallocation: %c\n", node->ch);

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

  WLNSymbol *copy = AllocateWLNSymbol(src->ch);
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
  // if the child cannot handle the new valence
  if ((child->num_edges + bond) > child->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges+bond, child->allowed_edges);
    return false;
  }
  // same for the parent
  if ((parent->num_edges + bond) > parent->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges+bond, parent->allowed_edges);
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
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges+diff, child->allowed_edges);
    return false;
  }
  
  if ((parent->num_edges + diff) > parent->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges+diff, parent->allowed_edges);
    return false;
  }

  child->num_edges += diff;
  parent->num_edges += diff;
  parent->orders[i] = bond;
  return true; 
}

bool make_aromatic(WLNSymbol *child, WLNSymbol *parent){

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
  if (parent->orders[i] == 4)
    return true; // save some work

  // set aromatics
  switch(parent->ch){
    case 'X':
    case 'C':
    case 'K':
      child->allowed_edges = 3;
      break;

    case 'Y':
    case 'N':
    case 'O':
      child->allowed_edges = 2;
      break;

    case 'P':
    case 'S':
      child->allowed_edges = 4;
      break;

    case '*':
      fprintf(stderr,"Error: aromaticity for specific elemental definitions currently unsupported\n");
      return false; 
    
    default:
      fprintf(stderr,"Error: can not make %c symbol aromatic, please check definitions\n");
      return false;
  }

  // set aromatics
  switch(child->ch){
    case 'X':
    case 'C':
    case 'K':
      child->allowed_edges = 3;
      break;

    case 'Y':
    case 'N':
    case 'O':
      child->allowed_edges = 2;
      break;

    case 'P':
    case 'S':
      child->allowed_edges = 4;
      break;

    case '*':
      fprintf(stderr,"Error: aromaticity for specific elemental definitions currently unsupported\n");
      return false; 
    
    default:
      fprintf(stderr,"Error: can not make %c symbol aromatic, please check definitions\n");
      return false;
  }

  

  if (child->num_edges > child->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", child->ch,child->num_edges, child->allowed_edges);
    return false;
  }

  if (parent->num_edges > parent->allowed_edges)
  {
    fprintf(stderr, "Error: wln character[%c] is exceeding allowed connections %d/%d\n", parent->ch,parent->num_edges, parent->allowed_edges);
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




/* struct to hold pointers for the wln ring */
struct WLNRing
{

  unsigned int size;
  bool heterocyclic;

  std::vector<unsigned int> rings;
  std::map<unsigned char, WLNSymbol *> locants;
  std::map<WLNSymbol*,unsigned char> locants_ch;


  // keep this simple

  WLNRing()
  {
    size = 0;
    heterocyclic = false;
  }
  ~WLNRing(){};

  // both lookups needed for QOL in ring building
  WLNSymbol* assign_locant(unsigned char loc,unsigned char type){
    WLNSymbol *locant = 0; 
    locant = AllocateWLNSymbol(type);
    locants[loc] = locant; 
    locants_ch[locant] = loc;
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
      unsigned int cur_int = locant_integer_map[map_iter->first] - 1; // gives zero index for distance matrix
      for (WLNSymbol *child : current->children){
        unsigned int child_int = locant_integer_map[locants_ch[child]] - 1;
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


  /* uses RP-PATH n3 algoritm from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2765087/ */
  bool SSRS(std::vector<std::vector<unsigned char>> &ring_subset, 
            std::vector<std::pair<unsigned int, unsigned char>> &ring_assignments)
  {

    unsigned int *distance = distance_matrix(locants.size());
    
    
    free(distance);
    
    return true;
  }


  bool CreateMono(unsigned int local_size, bool aromatic){

    WLNSymbol *head = 0; 
    WLNSymbol *prev = 0;
    WLNSymbol *current = 0; 

    bool state = true;

    size = local_size; // set for locant bonding outside of ring functions
    
    // assume already assigned locants
    for (unsigned int i=1;i<=local_size;i++){

      unsigned char loc = integer_locant_map[i];

      if(!locants[loc]){
        current = assign_locant(loc,'C');
        current->allowed_edges = 4;
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
  bool CreatePoly(std::vector<std::pair<unsigned int,unsigned char>> &ring_assignments, std::vector<bool> &aromaticity){
     
    // perform the aromatic denotion check
    if (ring_assignments.size() != aromaticity.size()){
      fprintf(stderr,"Error: mismatch between number of rings and aromatic assignments\n");
      return false; 
    }

    unsigned int prev_size = 1; // offset by 1 to get locant 'A' 
    unsigned int local_size = 0; 
    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i]; 

      if(local_size)
        local_size += component.first - 2;
      else
        local_size = component.first;

      prev_size = local_size;
    }

    // set the global size; 
    size = local_size; 

    // create all the nodes in a large straight chain with aromatics evaulated

    WLNSymbol *current = 0; 
    WLNSymbol *prev = 0; 
    for (unsigned int i=1;i<=size;i++){
      unsigned char loc = integer_locant_map[i];
      if(!locants[loc]){
        current = assign_locant(loc,'C');
        current->allowed_edges = 4;
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
    unsigned char bind_1 = '\0';
    unsigned char bind_2 = '\0';
    unsigned int fuses = 0; 

    // reverse the aromaticity assignments; 
    std::reverse(aromaticity.begin(), aromaticity.end());

    for (unsigned int i=0;i<ring_assignments.size();i++){
      std::pair<unsigned int, unsigned char> component = ring_assignments[i];
      bind_1 = component.second;
      bool aromatic = aromaticity[i];

      std::vector<unsigned char> ring_path;
      // first pair can be calculated directly without a path travel
      if(!fuses){
        bind_2 = bind_1 + component.first - 1; // includes start atom
        for (unsigned int i=0; i<component.first;i++)
          ring_path.push_back(bind_1+i);
      }
      else{
        //there needs to be a graph travel here taking the longest locant

        // 1. starting on bind_1, travel n-1 places through the maximal locant path, to calculate fuse
        
        // annoyingly n2 ... 
        WLNSymbol *path = locants[bind_1];
        unsigned char highest_loc = '\0';
        for (unsigned int i=0;i<component.first - 1; i++){
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

    
      if(!link_symbols(locants[bind_2],locants[bind_1],1)){
        fprintf(stderr,"Error: error in bonding locants together, check ring notation\n");
        return false;
      }

      if(aromatic){
        for (unsigned int i=1; i<ring_path.size();i+=1){
          unsigned char par = ring_path[i-1];
          unsigned char chi = ring_path[i];
          if(!make_aromatic(locants[chi],locants[par])){
            fprintf(stderr,"Error: error in changing aromaticity - check ring notation\n");
            return false; 
          }
        }
        if(!make_aromatic(locants[ring_path.back()],locants[ring_path.front()])){
          fprintf(stderr,"Error: error in changing aromaticity - check ring notation\n");
          return false; 
        }
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


  void FormWLNRing(std::string &block, unsigned int start){

    enum RingType{ MONO=0, POLY=1, PERI=2, BRIDGED=3, PSDBRIDGED = 4}; 
    unsigned int ring_type = MONO;   // start in mono and climb up
    unsigned int end = 0;

    bool warned = false; // limit warning messages to console
    bool heterocyclic       = false;  // L|T designator can throw warnings

    // -- stages --

    // based off certain rules, we can enforce ordering
    bool multi_completed      = false;
  
    // -- paths -- 
    bool pending_component  = false;
    bool pending_multi      = false; 
    bool pending_pseudo     = false; 
    bool pending_bridge     = false;
    bool pending_aromatics  = false;

    unsigned int expected_locants     = 0;
    unsigned char ring_size_specifier  = '\0';
    unsigned char positional_locant = '\0'; // this changes

    std::vector<bool> aromaticity; 

    std::vector<unsigned char> fuses; // read as pairs
    std::vector<unsigned char> bridge_locants;
    std::vector<unsigned char> multicyclic_locants;
    
    std::vector<std::pair<unsigned int, unsigned char>>  ring_components; 
   
   
    // things can now have multiple options, need to cover all variations
    // and move the pending bools around to give proper notation read. 

    WLNSymbol *wln_locant = 0; 

    for (unsigned int i=0;i<block.size();i++){
      unsigned char ch = block[i];

      switch(ch){
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          if (pending_component){
            if(!positional_locant)
              ring_components.push_back({ch - '0','A'});
            else{
              ring_components.push_back({ch - '0',positional_locant});
              positional_locant = '\0'; // reset the assigner
            }
            break;
          }
          else{
            pending_multi   = true;
            expected_locants = ch - '0';
            break;
          }
          
            
        case '/':
          expected_locants = 2; 
          pending_pseudo = true;
          ring_type = PSDBRIDGED; 
          break; 

        case '-':
          break;

        // aromatics
        case '&':
          pending_aromatics = true;
          aromaticity.push_back(true);
          break;

        case ' ':
          if(expected_locants){
            fprintf(stderr,"Error: %d more locants expected before space seperator\n",expected_locants);
            Fatal(start+i);
          }
          if(block[i-1] == ' '){
            fprintf(stderr,"Error: double spacing in ring notation is not allowed\n");
            Fatal(start+i);
          }


          // resets any pendings and set states
          if(pending_multi){
            pending_multi     = false;
            multi_completed   = true;

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
          pending_component = false;
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
        case 'R':
        case 'S':
        case 'U':
        case 'V':
        case 'W':
        case 'X':
        case 'Y':
        case 'Z':
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

            if(multi_completed && !ring_size_specifier){ // a size specifier is always needed
              ring_size_specifier = ch;
              positional_locant = ch;
            }
            else{
              positional_locant = ch;
              pending_component = true;
              pending_bridge = true;
            }
            break;
          }
          else if (positional_locant){
            pending_bridge = false;
            pending_component = false;

            if (opt_debug)
              fprintf(stderr,"  assigning WLNSymbol %c to position %c\n",ch,positional_locant);

            WLNSymbol *new_locant = 0; 

            switch(ch){
              case 'S':
              case 'P':
                if(!heterocyclic)
                  warned = true;
                new_locant = assign_locant(positional_locant,ch);
                new_locant->allowed_edges = 5;

                positional_locant++; // allows inline defition continuation
                break;

              case 'Y':
                new_locant = assign_locant(positional_locant,ch);
                new_locant->allowed_edges = 3;
                positional_locant++; // allows inline defition continuation
                break;
              case 'N':
                if(!heterocyclic)
                  warned = true;
                new_locant = assign_locant(positional_locant,ch);
                new_locant->allowed_edges = 3;
                positional_locant++; // allows inline defition continuation
                break;

              case 'V':
                new_locant = assign_locant(positional_locant,ch);
                new_locant->allowed_edges = 2;
                positional_locant++; // allows inline defition continuation
                break;

              case 'M':
              case 'O':
                if(!heterocyclic)
                  warned = true;
                new_locant = assign_locant(positional_locant,ch);
                new_locant->allowed_edges = 2;
                positional_locant++; // allows inline defition continuation
                break;

              case 'X':
                new_locant = assign_locant(positional_locant,ch);
                new_locant->allowed_edges = 4;
                positional_locant++; // allows inline defition continuation
                break;
              case 'K':
                if(!heterocyclic)
                  warned = true;

                new_locant = assign_locant(positional_locant,ch);
                new_locant->allowed_edges = 4;
                positional_locant++; // allows inline defition continuation
                break;

              case 'U':
                if(opt_debug)
                  fprintf(stderr,"  increasing bond order from %c to %c by 1\n",positional_locant,positional_locant+1);
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

          if(i==0){
           heterocyclic = false; 
           pending_component = true;
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
          else if(block[i-1] == ' '){

            if(multi_completed && !ring_size_specifier){ // a size specifier is always needed
              ring_size_specifier = ch;
              positional_locant = ch;
            }
            else{
              positional_locant = ch;
              pending_component = true;
              pending_bridge = true;
            }
            break;
          }
          else{
            positional_locant = ch;
            break;
          }

        case 'T':
          if(i==0){
            heterocyclic = true; 
            pending_component = true;
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
          else if (positional_locant && positional_locant == 'T'){
            pending_aromatics = true;
            aromaticity.push_back(false);
            positional_locant = 'T';
            break;
          }
          else if (i == block.size()-2){
            // this now must be the aromatic designator 
            if(opt_debug)
              fprintf(stderr,"  removing all aromaticity with singular T notation\n");

            pending_aromatics = true;
            for (unsigned int i=0;i<ring_components.size();i++)
              aromaticity.push_back(false);

            break;
          }
          else if(block[i-1] == ' '){

            if(multi_completed && !ring_size_specifier){ // a size specifier is always needed
              ring_size_specifier = ch;
              positional_locant = ch;
            }
            else{
              positional_locant = ch;
              pending_component = true;
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
          if(i==block.size() - 1){
            if(!pending_aromatics){
              for(unsigned int i=0;i<ring_components.size();i++)
                aromaticity.push_back(true);
            }
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
          else if(block[i-1] == ' '){

            if(multi_completed && !ring_size_specifier){ // a size specifier is always needed
              ring_size_specifier = ch;
              positional_locant = ch;
            }
            else{
              positional_locant = ch;
              pending_component = true;
              pending_bridge = true;
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

    // set the ring type

    if (ring_components.size() > 1 && ring_type < PERI)
      ring_type = POLY;
    
    

        
    // debug here

    if (opt_debug){
      
      switch(ring_type){
        case 0:
          fprintf(stderr,"  ring type: MONO\n");
          break;
        case 1:
          fprintf(stderr,"  ring type: POLY\n");
          break;
        case 2:
          fprintf(stderr,"  ring type: PERI\n");
          break;
        case 3:
          fprintf(stderr,"  ring type: BRIDGED\n");
          break;
        case 4:
          fprintf(stderr,"  ring type: PSDBRIDGED\n");
          break;
      }
      
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


      fprintf(stderr,"  size denotion: %c\n",ring_size_specifier);
      fprintf(stderr,"  heterocyclic: %s\n", heterocyclic ? "yes":"no");

    }

    if(warned)
      fprintf(stderr,"Warning: heterocyclic ring notation required for inter atom assignment, change starting 'L' to 'T'\n");
      
    
    bool state = true;
    switch(ring_type){
      case MONO:
        state = CreateMono(ring_components[0].first,aromaticity[0]);
        break;
      case POLY:
        state = CreatePoly(ring_components,aromaticity);
        break;
      case PERI:
      case BRIDGED:
      case PSDBRIDGED:
        break;
    }

    if (!state)
      Fatal(end);
    
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


  WLNSymbol *define_element(std::vector<unsigned char> &special)
  {

    // allocate a special wln
    WLNSymbol *created_wln = AllocateWLNSymbol('*');

    // some fancy switching

    switch (special[0])
    {

    case 'A':
      if (special[1] == 'C')
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
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'B':
      if (special[1] == 'A')
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
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'C':
      if (special[1] == 'A')
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
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'D':
      if (special[1] == 'B')
        created_wln->special = "Db";
      else if (special[1] == 'S')
        created_wln->special = "Ds";
      else if (special[1] == 'Y')
        created_wln->special = "Dy";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'E':
      if (special[1] == 'R')
        created_wln->special = "Er";
      else if (special[1] == 'S')
        created_wln->special = "Es";
      else if (special[1] == 'U')
        created_wln->special = "Eu";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'F':
      if (special[1] == 'E')
        created_wln->special = "Fe";
      else if (special[1] == 'L')
        created_wln->special = "Fl";
      else if (special[1] == 'M')
        created_wln->special = "Fm";
      else if (special[1] == 'R')
        created_wln->special = "Fr";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'G':
      if (special[1] == 'A')
        created_wln->special = "Ga";
      else if (special[1] == 'D')
        created_wln->special = "Gd";
      else if (special[1] == 'E')
        created_wln->special = "Ge";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'H':
      if (special[1] == 'E')
        created_wln->special = "Ha";
      else if (special[1] == 'F')
        created_wln->special = "Hf";
      else if (special[1] == 'G')
        created_wln->special = "Hg";
      else if (special[1] == 'O')
        created_wln->special = "Ho";
      else if (special[1] == 'S')
        created_wln->special = "Hs";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'I':
      if (special[1] == 'N')
        created_wln->special = "In";
      else if (special[1] == 'R')
        created_wln->special = "Ir";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'K':
      if (special[1] == 'R')
        created_wln->special = "Kr";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'L':
      if (special[1] == 'A')
        created_wln->special = "La";
      else if (special[1] == 'I')
        created_wln->special = "Li";
      else if (special[1] == 'R')
        created_wln->special = "Lr";
      else if (special[1] == 'U')
        created_wln->special = "Lu";
      else if (special[1] == 'V')
        created_wln->special = "Lv";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'M':
      if (special[1] == 'C')
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
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'N':
      if (special[1] == 'A')
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
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'O':
      if (special[1] == 'G')
        created_wln->special = "Og";
      else if (special[1] == 'S')
        created_wln->special = "Os";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'P':
      if (special[1] == 'A')
        created_wln->special = "Pa";
      else if (special[1] == 'B')
        created_wln->special = "Pb";
      else if (special[1] == 'D')
        created_wln->special = "Pd";
      else if (special[1] == 'M')
        created_wln->special = "Pm";
      else if (special[1] == 'O')
        created_wln->special = "Po";
      else if (special[1] == 'R')
        created_wln->special = "Pr";
      else if (special[1] == 'T')
        created_wln->special = "Pt";
      else if (special[1] == 'U')
        created_wln->special = "Pu";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'R':
      if (special[1] == 'A')
        created_wln->special = "Ra";
      else if (special[1] == 'B')
        created_wln->special = "Rb";
      else if (special[1] == 'E')
        created_wln->special = "Re";
      else if (special[1] == 'F')
        created_wln->special = "Rf";
      else if (special[1] == 'G')
        created_wln->special = "Rg";
      else if (special[1] == 'H')
        created_wln->special = "Rh";
      else if (special[1] == 'N')
        created_wln->special = "Rn";
      else if (special[1] == 'U')
        created_wln->special = "Ru";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'S':
      if (special[1] == 'B')
        created_wln->special = "Sb";
      else if (special[1] == 'C')
        created_wln->special = "Sc";
      else if (special[1] == 'E')
        created_wln->special = "Se";
      else if (special[1] == 'I')
        created_wln->special = "Si";
      else if (special[1] == 'M')
        created_wln->special = "Sm";
      else if (special[1] == 'N')
        created_wln->special = "Sn";
      else if (special[1] == 'R')
        created_wln->special = "Sr";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'T':
      if (special[1] == 'A')
        created_wln->special = "Ta";
      else if (special[1] == 'B')
        created_wln->special = "Tb";
      else if (special[1] == 'C')
        created_wln->special = "Tc";
      else if (special[1] == 'E')
        created_wln->special = "Te";
      else if (special[1] == 'H')
        created_wln->special = "Th";
      else if (special[1] == 'I')
        created_wln->special = "Ti";
      else if (special[1] == 'L')
        created_wln->special = "Tl";
      else if (special[1] == 'M')
        created_wln->special = "Tm";
      else if (special[1] == 'S')
        created_wln->special = "Ts";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'X':
      if (special[1] == 'E')
        created_wln->special = "Xe";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'Y':
      if (special[1] == 'B')
        created_wln->special = "Yb";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    case 'Z':
      if (special[1] == 'N')
        created_wln->special = "Zn";
      else if (special[1] == 'R')
        created_wln->special = "Zr";
      else
      {
        fprintf(stderr, "Error: invalid element symbol in special definition\n");
        return (WLNSymbol *)0;
      }
      break;

    default:
      fprintf(stderr, "Error: invalid character in special definition switch\n");
      return (WLNSymbol *)0;
    }

    created_wln->allowed_edges = 8; // allow on octet default for these species.

    return created_wln;
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
      fprintf(stderr, "  popping %d symbols down the stack: mode(%d) prev[%c]\n", pops, hard, prev->ch);

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
      fprintf(stderr, "Error: assigning locant greater than ring size - %d\n",s_ring->size);
      Fatal(i);
    }

  }

  /* a global segmentation using both rule sets - start merging */
  bool ParseWLNString(const char *wln, unsigned int len)
  {

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

    // allows consumption of notation after block parses
    unsigned int block_start = 0;
    unsigned int block_end = 0;

    unsigned int pop_ticks = 0;  // '&' style popping
    unsigned int bond_ticks = 0; // 'U' style bonding

    for (unsigned int i = 0; i < len; i++)
    {
      unsigned char ch = wln[i];

      if (opt_debug)
        fprintf(stderr, "Parsing: %c\n", ch);

      switch (ch)
      {

      case '0': // cannot be lone, must be an addition to another num
        if (pending_closure || pending_special)
        {
          break;
        }
        if (i == 0)
          Fatal(i);
        else if (i > 0 && !std::isdigit(wln[i - 1]))
          Fatal(i);
        else
        {
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
        if (pop_ticks)
        {
          prev = pop_standard_stacks(pop_ticks, branch_stack, linker_stack, prev, i);
          pop_ticks = 0;
        }

        if (pending_closure || pending_special)
        {
          break;
        }

        curr = AllocateWLNSymbol(ch);
        curr->set_type(STANDARD);
        curr->set_edges(3);

        create_bond(curr, prev, bond_ticks, i);

        prev = curr;
        break;

        // carbons

      case 'Y':

        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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

      case 'K':
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
      case 'H':
        if (pending_closure || pending_special)
        {
          break;
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
          Fatal(i);

        break;

        // ring notation

      case 'J':
        if (pending_special)
        {
          break;
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
            prev->type = LINKER; // spiros are normal rings with dual linker notation
            prev->previous->type = LINKER;
            pending_spiro = false;
          }

          // does the incoming locant check
          if (prev)
          {
            if (ring->locants[prev->ch])
              create_bond(ring->locants[prev->ch],prev,bond_ticks,i);
            else
            {
              fprintf(stderr, "Error: attaching inline ring with out of bounds locant assignment\n");
              Fatal(i);
            }
          }

          bond_ticks = 0;
          pending_closure = false;
        }
        else
          Fatal(i);

        break;

      case 'L':
      case 'T':
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_closure || pending_special)
        {
          break;
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
        if (pending_special)
        {
          pending_special = false;
          block_start = 0;
          block_end = 0;
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
        if (pending_closure || pending_special)
        {
          break;
        }

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
        if (!pending_inline_ring)
        {
          pending_inline_ring = true;

          // send the linker into its own stack
          if (!branch_stack.empty() && (branch_stack.top()->num_edges < branch_stack.top()->allowed_edges))
            linker_stack.push(branch_stack.top());
        }

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

          block_end = i;
          curr = AllocateWLNSymbol('*');
          curr->add_special(block_start, block_end);

          block_start = 0;
          block_end = 0;

          create_bond(curr, prev, bond_ticks, i);

          bond_ticks = 0;
          prev = curr;
          pending_special = false;
        }
        else if (!pending_special)
          pending_special = true;

        break;

      case '/':
        if (pending_closure || pending_special)
        {
          break;
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
      if (node->ch == '*')
        fprintf(fp, "[shape=circle,label=\"%s\"];\n", node->special.c_str());
      else if(node->type == LOCANT)
        fprintf(fp, "[shape=circle,label=\"%c\",color=blue];\n", node->ch);
      else if (node->type == RING)
        fprintf(fp, "[shape=circle,label=\"%c\",color=green];\n", node->ch);
      else if (node->type == LINKER)
        fprintf(fp, "[shape=circle,label=\"%c\",color=red];\n", node->ch);
      else
        fprintf(fp, "[shape=circle,label=\"%c\"];\n", node->ch);

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

  wln_graph.ParseWLNString(wln, strlen(wln));
  Reindex_lookups();

  // create the wln dotfile
  if (opt_wln2dot)
  {
    FILE *fp = 0;
    fp = fopen("wln-graph.dot", "w");
    if (!fp)
    {
      fprintf(stderr, "Error: could not open compiler dump file\n");
      fclose(fp);
      return 1;
    }
    else
    {
      wln_graph.WLNDumpToDot(fp);
      fclose(fp);
    }
  }

  return 0;
}