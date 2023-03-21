

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector>
#include <stack>
#include <map>
#include <utility> // std::pair
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
'UU' bonding also added here */
bool link_symbols(WLNSymbol *child, WLNSymbol *parent, unsigned int bond)
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
  parent->orders.push_back(bond);
  return true;
}


/* struct to hold pointers for the wln ring */
struct WLNRing
{

  unsigned int size;
  bool heterocyclic;

  std::vector<unsigned int> rings;
  std::map<unsigned char, WLNSymbol *> locants;


  // keep this simple

  WLNRing()
  {
    size = 0;
    heterocyclic = false;
  }
  ~WLNRing(){};

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


  unsigned int create_ring(unsigned int size, unsigned char start_locant){

    WLNSymbol *head = 0; 
    WLNSymbol *prev = 0;
    WLNSymbol *current = 0; 

    unsigned int locant_num = locant_integer_map[start_locant];
    
    for (unsigned int i=0;i<size;i++){
      current = AllocateWLNSymbol('C');
      current->type = RING;

      if(locants[integer_locant_map[locant_num]]){
        fprintf(stderr,"Error: overwriting locant in ring definition!\n");
        return 0; 
      }

      locants[integer_locant_map[locant_num]] = current;
      locant_num++; 

      if (!head)
        head = current; 

      if(prev)
        link_symbols(current,prev,0);

      prev = current;
    }

    link_symbols(head,prev,0);

    return locant_num; 
  }


  // creates a loop between two or more wln positions - locants are consecutive from the given char
  // ring size is implied with the two given atoms -- inputting 6 will create 4 new symbols

  // return last locant!
  unsigned int wrap_ring(std::vector<WLNSymbol*> &given_path, unsigned int size, unsigned char start_locant){
    
   
    unsigned int locant_num = locant_integer_map[start_locant + 1];
    
    WLNSymbol *prev = given_path.front();
    WLNSymbol *current = 0; 

    for (unsigned int i=0;i<size-given_path.size();i++){
      current = AllocateWLNSymbol('C');
      current->allowed_edges = 4;
      current->type = RING;

      if(locants[integer_locant_map[locant_num]]){
        fprintf(stderr,"Error: overwriting locant in ring definition!\n");
        return 0; 
      }

      locants[integer_locant_map[locant_num]] = current;
      locant_num++;

      if(prev)
        link_symbols(current,prev,0);

      prev = current; 
    }

    link_symbols(given_path.back(),prev,0);
     
    return locant_num - 1; // will always advance 1 extra past allocation
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
          break; 

        case '-':
          break;

        // aromatics
        case '&':
          pending_aromatics = true;
          if (positional_locant == 'T') // back search if we start this with T
            aromaticity.push_back(false);

          aromaticity.push_back(true);
          break;

        case ' ':
          if(expected_locants){
            fprintf(stderr,"Error: %d more locants expected before space seperator\n",expected_locants);
            Fatal(start+i);
          }

          // resets any pendings and set states
          if(pending_multi){
            pending_multi     = false;
            multi_completed   = true;
          }
          else if (pending_bridge){ 
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

            if(multi_completed){ // a size specifier is always needed
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

            if(!heterocyclic)
              warned = true;

            switch(ch){
              case 'S':
              case 'P':
                locants[positional_locant] = AllocateWLNSymbol(ch);
                locants[positional_locant]->allowed_edges = 5;
                positional_locant++; // allows inline defition continuation
                break;

              case 'Y':
              case 'N':
                locants[positional_locant] = AllocateWLNSymbol(ch);
                locants[positional_locant]->allowed_edges = 3;
                positional_locant++; // allows inline defition continuation
                break;

              case 'V':
              case 'M':
              case 'O':
                locants[positional_locant] = AllocateWLNSymbol(ch);
                locants[positional_locant]->allowed_edges = 2;
                positional_locant++; // allows inline defition continuation
                break;

              case 'X':
              case 'K':
                locants[positional_locant] = AllocateWLNSymbol(ch);
                locants[positional_locant]->allowed_edges = 4;
                positional_locant++; // allows inline defition continuation
                break;

              case 'U':
                if(opt_debug)
                  fprintf(stderr,"  increasing bond order from %c to %c by 1\n",positional_locant,positional_locant++);
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
          
          else{
            
            break;
          }

        case 'T':
          if(i==0){
            heterocyclic = true; 
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
              fprintf(stderr,"  removing all aromaticty with singular T notation\n");

            for (unsigned int i=0;i<ring_components.size();i++)
              aromaticity.push_back(false);
            
            break;
          }
          else{
            positional_locant = ch;
            break;
          }

        //closure
        case 'J':
          if(i==block.size() - 1){
            if(!pending_aromatics){
              for(unsigned int i=0;i<ring_components.size();i++)
                aromaticity.push_back(true);
            }
            break;
          }
            


          


        default:
          fprintf(stderr,"Error: unrecognised symbol in ring definition: %c\n",ch);
          Fatal(start + i);
      }
       
    }

    

    
    // debug here

    if (opt_debug){

      fprintf(stderr,"  ring components: ");
      for (std::pair<unsigned int, unsigned char> comp : ring_components)
        fprintf(stderr,"%d(%c) ",comp.first,comp.second);
      fprintf(stderr,"\n");

      fprintf(stderr,"  aromaticity: ");
      for (bool aromatic : aromaticity)
        fprintf(stderr,"%d ",aromatic);
      fprintf(stderr,"\n");

      fprintf(stderr,"  multicyclic points: ");
      for (unsigned char loc : multicyclic_locants)
        fprintf(stderr,"%c ",loc);
      fprintf(stderr,"\n");

      fprintf(stderr,"  bridge points: ");
      for (unsigned char loc : bridge_locants)
        fprintf(stderr,"%c ",loc);
      fprintf(stderr,"\n");

      fprintf(stderr,"  hard fuses: ");
      for (unsigned int i=1;i<fuses.size();i++)
        fprintf(stderr,"(%c --> %c) ",fuses[i-1],fuses[i]);
      fprintf(stderr,"\n");


      fprintf(stderr,"  size denotion: %c\n",ring_size_specifier);
      fprintf(stderr,"  heterocyclic: %s\n", heterocyclic ? "yes":"no");

    }

    if(warned)
      fprintf(stderr,"Warning: heterocyclic ring notation required for inter atom assignment, change starting 'L' to 'T'\n");
              
    
  }

  
  // creates and allocates the ring in WLNSymbol struct style
  bool create_symbol_ring(std::string block)
  {

    /*
    we assume this start after the ring blocks
    i.e 'L66 AO TJ' -->' AO TJ'
    */

    bool inter_ring = false;
    for (unsigned char ch : block)
    {
      if (ch == ' ')
      {
        inter_ring = true;
        break;
      }
    }

    if (inter_ring)
    {

      // split the string on the spaces, and then process the blocks
      std::istringstream ss(block);
      std::string del;

      unsigned int locant = 1; // use the maps to assign locants
      WLNSymbol *assignment = 0;

      unsigned int consumed = 0;
      while (getline(ss, del, ' '))
      {

        unsigned int assign_size = del.size();

        if (assign_size < 1)
          continue; // not sure i can help this?
        else if (assign_size == 1)
        {
          consumed += 2; // +1 for the space!

          // handles the lone 'J'
          if (consumed == block.size())
            break;

          // bridge defintions --> figure out a way to ignore the last 'J' as a 'J' bridge is intirely plausible
          if (opt_debug)
            fprintf(stderr, "  assigning bridge: %c\n", del[0]);

          //bridge_points.push_back(del[0]);
        }
        else if (assign_size > 1)
        {
          consumed += assign_size + 1; // +1 for the space!
          // standard atomic definitions

          unsigned int back = assign_size - 1;
          if (del[assign_size] == 'J')
          {
            back += -1;
            // if (!aromatic)
            //   back += -1;
          }

          // process the locants as expected in standard notation

          locant = locant_integer_map[del[0]];
          for (unsigned int i = 1; i < back; i++)
          {

            // need to add support for specials here

            if (del[i] != 'U')
            {
              if (opt_debug)
                fprintf(stderr, "  assigning symbol: locant(%c) --> symbol(%c)\n", integer_locant_map[locant], del[i]);

              locants[integer_locant_map[locant]] = AllocateWLNSymbol(del[i]);
              switch (del[i])
              {

              case 'B':
              case 'I':
                locants[integer_locant_map[locant]]->allowed_edges = 3;
                break;

              case 'N':
              case 'K':
                locants[integer_locant_map[locant]]->allowed_edges = 3;
                break;

              case 'O':
              case 'M':
              case 'V':
                locants[integer_locant_map[locant]]->allowed_edges = 2;
                break;

              case 'P':
              case 'S':
                locants[integer_locant_map[locant]]->allowed_edges = 5;
                break;
              }
            }

            locant++;
          }
        }
      }
    }

    // now we assign everything not coverered as 1's!
    // a nice property is that all locants are binded to eachother

    WLNSymbol *prev = 0;
    WLNSymbol *head = 0;
    for (unsigned int i = 1; i <= size; i++)
    {
      unsigned char loc = integer_locant_map[i];
      if (!locants[loc])
      {
        locants[loc] = AllocateWLNSymbol('C');
        locants[loc]->allowed_edges = 4;
      }

      locants[loc]->type = RING; 

      // if(aromatic)
      //   locants[loc]->allowed_edges += -1; // take off 1 when aromatic!

      if (!head)
        head = locants[loc];

      if (prev)
      {
        // if (aromatic){
        //   if(!link_aromatics(locants[loc], prev))
        //     return false;
        // }
        // else{
        //   if(!link_symbols(locants[loc], prev, 1))
        //     return false;
        // }
      }

      prev = locants[loc];
    }

    // if (aromatic){
    //   if(!link_aromatics(head, prev))
    //     return false;
    // }else{
    //   if(!link_symbols(head, prev, 1))
    //     return false;
    // }
    // fuse any points given
    // for (std::pair<unsigned char, unsigned char> fuse : fuse_points)
    // { 
    //   if (opt_debug)
    //     fprintf(stderr, "  fusing position %c to position %c\n",fuse.first,fuse.second);

    //   // if (aromatic){
    //   //   if(!link_aromatics(locants[fuse.first], locants[fuse.second]))
    //   //     return false;
    //   // }else{
    //   //   if(!link_symbols(locants[fuse.first], locants[fuse.second], 1))
    //   //     return false;
    //   // }
    // }

    return true;
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
      fprintf(stderr, "Error: assigning locant greater than ring size\n");
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
        curr->set_edges(2);

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

    std::map<WLNSymbol *, bool>::iterator sym_iter;

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