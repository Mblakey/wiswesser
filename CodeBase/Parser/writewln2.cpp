
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector>
#include <stack>
#include <map>
#include <deque>
#include <algorithm> // std::sort

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
static bool opt_valstrict = false; 
static bool opt_verbose = false;
static bool opt_canonical = false; 
static bool opt_returnwln = false;


// character type and instruction types
enum WLNType {SINGLETON = 0, BRANCH = 1, LINKER = 2, TERMINATOR = 3}; 
enum WLNCode {ROOT = 0, STANDARD = 1, LOCANT = 2, CYCLIC = 3, BRIDGED = 4, SPIRO = 5, IONIC = 6};

// rule 2 - hierarchy - rules have diverged due to end terminator char
std::map<unsigned char,unsigned int> char_hierarchy = 
{
  {' ',1}, {'-',2}, {'/',3}, 
  {'0',4},{'1',5},{'2',6},{'3',7},{'4',8},{'5',9},{'6',10},{'7',11},{'8',12},{'9',13},
  {'A',14},{'B',15},{'C',16},{'D',17},{'E',18},{'F',19},{'G',20},{'H',21},{'I',22},
  {'J',23},{'K',24},{'L',25},{'M',26},{'N',27},{'O',28},{'P',29},{'Q',30},{'R',31},
  {'S',32},{'T',33},{'U',34},{'V',35},{'W',36},{'X',37},{'Y',38},{'Z',40},{'&',41}
};

const char *code_hierarchy[] = {"ROOT","STANDARD", "LOCANT", "CYCLIC","BRIDGED", "SPIRO", "IONIC"};



struct WLNSymbol{

  unsigned char ch; 
  unsigned int type; 

  unsigned int allowed_edges;
  unsigned int num_edges;  

  WLNSymbol *prev; // should be a single term - wln symbol only has one incoming
  std::vector<WLNSymbol*> children; // linked list of next terms chains 

  bool init(unsigned char inp_char){
    ch = inp_char;

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
        type = SINGLETON; 
        allowed_edges = 2;
        break;
      
      case 'A':
        type = SINGLETON;
        allowed_edges = 2;
        break;

      case 'B':   // boron
        type = BRANCH;
        allowed_edges = 3; 
        break;
      
       case 'C':  // shortcut carbon atom
        type = BRANCH;
        allowed_edges = 4;
        break;
      
      case 'D': 
        type = SINGLETON; 
        allowed_edges = 2;
        break;
      
      case 'E':  // halogens
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

      case 'R':
        type = SINGLETON; 
        allowed_edges = 2;
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
        fprintf(stderr,"Error: end of string null char accessed!\n");
        return false;

      default:
        fprintf(stderr,"Error: invalid wln symbol parsed: %c\n",ch);
        return false; 
    }

    prev = (WLNSymbol*)0;
    return true; 
  } 

};


/* struct to hold pointers for the wln ring - only for stack return */ 
struct WLNRing{
  unsigned int ring_size;
  
  bool aromatic; 
  bool heterocyclic; 

  std::vector<WLNSymbol*> locants; 
  std::map<unsigned char, WLNSymbol*> lookup; 


  void init(){
    ring_size = 0; 
    aromatic = false; 
    heterocyclic = false;
  }
  
};

struct WLNGraph{

  WLNSymbol *root; 
  unsigned int wln_nodes = 0; 
  unsigned int wln_rings = 0; 
  std::vector<WLNSymbol*> symbol_mempool;
  std::vector<WLNRing*>   ring_mempool;

  std::map<WLNRing*,WLNSymbol*> ring_access; // access the ring struct from locant A pointer
  
  WLNGraph(): root{(WLNSymbol*)0},wln_nodes{0}{};
  ~WLNGraph(){
    for (WLNSymbol* allocedwln : symbol_mempool){  // if error free all the symbol_mempool --> stop leak 
      free(allocedwln); 
      allocedwln = 0;
    }
    for (WLNRing* allocedring : ring_mempool){  // if error free all the symbol_mempool --> stop leak 
      free(allocedring); 
      allocedring = 0;
    }
  }

  /* for std::sort on canonicalise */
  bool char_comp(const WLNSymbol* a, const WLNSymbol* b){
    return char_hierarchy[a->ch] > char_hierarchy[b->ch]; 
  }

  WLNSymbol* AllocateWLNSymbol(unsigned char ch){
    wln_nodes++; 
    WLNSymbol *wln = (WLNSymbol*)malloc( sizeof(WLNSymbol) ); 
    if(wln->init(ch)) 
      symbol_mempool.push_back(wln);
    else 
      return (WLNSymbol *)0;
    
    return wln; 
  }

  WLNRing* AllocateWLNRing(){
    wln_rings++; 
    WLNRing *wln_ring = (WLNRing*)malloc(sizeof(WLNRing));
    wln_ring->init();
    ring_mempool.push_back(wln_ring);
    return wln_ring;
  }


  /* e.g creates a 6-6 ring from 10 atoms - binds to symbol if given, returns locant A */
  WLNSymbol* create_ring( unsigned int atoms, 
                          std::vector<unsigned int> fuses, 
                          WLNSymbol *bind)
  {

    WLNRing *wln_ring = AllocateWLNRing();


    WLNSymbol *rhead = AllocateWLNSymbol('C'); 
    WLNSymbol *current = 0; 
    
    // 1) create a big ring
    WLNSymbol *prev = rhead; 
    for (unsigned int i=1; i<atoms; i++){
      current = AllocateWLNSymbol('C');
      add_symbol(current,prev);
    }

    // 2) loop back by adding to rhead; 
    add_symbol(rhead,current); 
    



    return rhead; 
  }

  bool handle_hypervalence(WLNSymbol *problem){

    // this will always lead to a positive ion species 
    switch(problem->ch){
      case 'M':   // tranforming this to a N is an easy solve
        if(opt_verbose)
          fprintf(stderr,"   transforming hypervalent M --> N\n");
        
        problem->ch = 'N';
        break; 

      case 'N':
        if(opt_verbose)
          fprintf(stderr,"   transforming hypervalent N --> K\n");
        
        problem->ch = 'N';
        break;


      case 'Y':  // can go to an X
        if (opt_verbose)
          fprintf(stderr,"   transforming hypervalent Y --> X\n");
        
        problem->ch = 'X';
        break; 

      default:
        if (opt_verbose) 
          fprintf(stderr,"Error: cannot handle hypervalent symbol: %c\n",problem->ch);
        return false;
    }
    return true; 
  }


  /* add src to the children vector of trg, handle hypervalent bonds if possible */
  bool add_symbol(WLNSymbol* src, WLNSymbol *trg){

    // handle exotic bonding - lookback 
    if (trg->ch == 'U'){
      if (trg->prev && trg->prev->ch == 'U')
        src->num_edges+=3;
      else
        src->num_edges+=2;
    }
    else
      src->num_edges++;

    if (src->num_edges > src->allowed_edges){
      if(!opt_valstrict){
        if(!handle_hypervalence(src))
          return false;
      }
      else {
        fprintf(stderr,"Error: (strict mode) hypervalence on WLN character %c\n",src->ch);
        return false; 
      }
    }
    
    if (trg->num_edges < trg->allowed_edges){
      trg->children.push_back(src);
      trg->num_edges++;
    }
    else {
      if(!opt_valstrict){
        if(!handle_hypervalence(trg))
          return false;
        else
          trg->children.push_back(src);
      }
      else {
        fprintf(stderr,"Error: (strict mode) hypervalence on WLN character %c\n",trg->ch);
        return false; 
      }
    }   

    return true;
  }

  /* this can perform the normal backtrack, excluding  '&' closure */ 
  WLNSymbol* backtrack_stack(std::stack<WLNSymbol*> &wln_stack){
    WLNSymbol *tmp = 0; 
    while(!wln_stack.empty()){
      tmp = wln_stack.top();
      if (tmp->type == BRANCH)  
        return tmp;
      wln_stack.pop();
    }
    
    // if it goes all the way that string complete
    return (WLNSymbol*)0;  
  }


  /* forces both the '&' closure and its parents branch to come off stack */
  WLNSymbol *force_closure(std::stack<WLNSymbol*> &wln_stack){
    WLNSymbol *tmp = 0; 
    unsigned int popped = 0;
    while(!wln_stack.empty()){
      tmp = wln_stack.top();
      if (tmp->type == BRANCH && popped > 1)  
        return tmp;
      wln_stack.pop();
      popped++; 
    }

    // if it goes all the way that string complete
    return (WLNSymbol*)0;
  }


  /* parses NON CYCLIC input string, mallocs graph nodes and sets up graph based on symbol read */
  WLNSymbol* ParseNonCyclic(const char *wln, unsigned int len){
    if(opt_verbose)
      fprintf(stderr,"   evaluating standard notation\n");

    std::stack <WLNSymbol*> wln_stack; 
    WLNSymbol *prev = 0;
    WLNSymbol *root = 0;

    WLNSymbol* created_wln = AllocateWLNSymbol(wln[0]);
    if (!created_wln)
      return (WLNSymbol*)0;
    wln_stack.push(created_wln);
    
    root = created_wln; 
    
    for (unsigned int i = 1; i<len; i++){
      
      created_wln = AllocateWLNSymbol(wln[i]);
      if (!created_wln)
        return (WLNSymbol*)0; 
      
      
      prev = wln_stack.top();
      wln_stack.push(created_wln); // push all of them

      if(!add_symbol(created_wln,prev))
        return (WLNSymbol*)0; 

      // options here depending on the forced closure of '&'
      if (created_wln->type == TERMINATOR){
        if (created_wln->ch == '&' && prev->type == BRANCH)
          prev = force_closure(wln_stack);
        else
          prev = backtrack_stack(wln_stack);
      }   
    }

    return root; // return start of the tree
  }



  /* reforms WLN string with DFS ordering */
  std::string ReformWLNString(WLNSymbol *root){
    std::string res; 

    std::stack<WLNSymbol*> wln_stack; 
    std::map<WLNSymbol*,bool> visit_map; 
    wln_stack.push(root);

    WLNSymbol *top = 0;
    while(!wln_stack.empty()){
      top = wln_stack.top(); 
      wln_stack.pop();
      visit_map[top] = true; 

      res.push_back(top->ch);

      for (WLNSymbol *child: top->children){
        if (!visit_map[child]){
          wln_stack.push(child);
        }
      }
    }

    return res; 
  }

#ifdef DEPRECATED
  /* uses a character sorting method to arrange the WLN symbols according to rule 2 BFS ensures a starting position */
  bool CanonicoliseWLN(WLNSymbol *root){

    fprintf(stderr,"CANONICOLISE NEEDS REWORKING -- SKIPPING\n");
    return true;

    std::deque<WLNSymbol*> wln_queue; 
    wln_queue.push_back(root);
    
    WLNSymbol *top = 0;
    while(!wln_queue.empty()){
      top = wln_queue.front();
      wln_queue.pop_front();

      if (top->children.size() > 1)
        std::sort(top->children.begin(),top->children.end(),char_comp);
      
      for (WLNSymbol* child : top->children){
        wln_queue.push_back(child);
      }
    }

    return true; 
  }

#endif

  /* dump wln tree to a dotvis file */
  void WLNDumpToDot(FILE *fp){

    // set up index map for node creation
    std::map <WLNSymbol*, unsigned int> index_map; 
    unsigned int glob_index = 0;
    for (WLNSymbol *node : symbol_mempool){
      index_map[node] = glob_index; 
      glob_index++;
    }
      
    fprintf(fp,"digraph WLNdigraph {\n");
    fprintf(fp,"  rankdir = LR;\n");
    for (WLNSymbol *node : symbol_mempool){
      fprintf(fp,"  %d",index_map[node]);
      fprintf(fp, "[shape=circle,label=\"%c\"];\n",node->ch);
      for (WLNSymbol *child : node->children){
        fprintf(fp,"  %d", index_map[node]);
        fprintf(fp," -> ");
        fprintf(fp,"%d\n",index_map[child]);
      }
    } 
    fprintf(fp,"}\n");
  }


};





/*  assumes a bi-atomic fuse, max = 6*6 for bicyclic */ 
unsigned int calculate_ring_atoms(unsigned int rings, unsigned int max_atoms){
  
  unsigned int term = rings - 2; 
  unsigned int shared_atoms = rings + term;

  return max_atoms - shared_atoms;
}


struct WLNInstruction{

  unsigned int state; 
  unsigned int start_ch; 
  unsigned int end_ch; 

  bool ring_linker; 
  
  WLNInstruction* parent; 
  std::vector<WLNInstruction*> next_instructions; 

  void init_state(unsigned int c){
    state = c;
    start_ch = 0; 
    end_ch = 0;
    parent = (WLNInstruction*)0; 
    ring_linker = false; 
  }

  void add_start(unsigned int pos){
    start_ch = pos; 
  }

  void add_end(unsigned int pos){
    end_ch = pos; 
  }

  void add_prev(WLNInstruction *src){
    parent = src; 
  }

  void display(){
    if(state == ROOT)
      fprintf(stderr,"instruction: %10s\n", "ROOT");
    else if (state == LOCANT)
      fprintf(stderr,"instruction: %10s contains: %c\n",code_hierarchy[state],wln[start_ch]);
    else{
      fprintf(stderr,"instruction: %10s contains: ",code_hierarchy[state]);
      for (unsigned int i=start_ch; i<=end_ch; i++)
        fprintf(stderr,"%c", wln[i]); 
      fprintf(stderr,"\n");
    }
  }

  /* handles single and bicyclic rings */
  WLNSymbol* construct_standard_ring()
  {

    if (state != CYCLIC){
      fprintf(stderr,"Error: constuct ring called on non-cyclic instruction\n");
      return (WLNSymbol*)0; 
    }

    unsigned int pos = 0; 
    unsigned char buffer[REASONABLE] = {'\0'};
    
    for (unsigned int i=start_ch; i<=end_ch; i++){
      buffer[pos++] = wln[i];
      if (pos == REASONABLE){
        fprintf(stderr,"Error: cyclic system greater than 1024 characters, limit hit\n");
        return (WLNSymbol*)0; 
      }
    }
      
    
    if(opt_verbose)
      fprintf(stderr,"constructing ring: %s\n",buffer);

    // globals
    bool ring_set     =   false; 
    bool heterocyclic =   false; 
    bool aromatic     =   false;

    // counters
    unsigned int ring_size = 0; 

    
    unsigned int it = 0; 
    while(buffer[it] != '\0'){

      // set the ring type
      if(it == 0){
        
        switch(buffer[it]){

          case 'L':
            heterocyclic = false;
            break;
          
          case 'T':
            heterocyclic = true;
            break; 

          default:
            fprintf(stderr,"Error: ring system starts with %c, must be L|T\n",buffer[it]);
            return (WLNSymbol*)0; 
        }
        // jump it along manually, allows some clever things
        it++;
        continue;
      }

      // setting the ring size
      if (!ring_set && (buffer[it] > '0' || buffer[it] < '9')){
        unsigned int fuses = 0; 
        while(std::isdigit(buffer[it]) && buffer[it] != '\0'){
          ring_size += buffer[it] - '0';
          fuses++; 
          it++; 
        }

        // refactor the ring size based on shared bonds
        ring_size = calculate_ring_atoms(fuses,ring_size);
        ring_set = true; 
        continue; 
      }



      it++; // prevent hanging on system build;
    }
    


    return (WLNSymbol*)0; 
  }

}; 


/* make a graph, split the string, call the subroutines - easy right? */
struct InstructionGraph{

  WLNInstruction *root; 
  unsigned int num_instructions; 
  std::vector<WLNInstruction*> instruction_pool; 

  InstructionGraph()
    :root{(WLNInstruction*)0},num_instructions{0}{};
  ~InstructionGraph(){
    for (WLNInstruction* instruction : instruction_pool)
      free(instruction);
  }

  WLNInstruction* add_instruction(unsigned int int_code, unsigned int i)
  {
    WLNInstruction *instruction = (WLNInstruction*)malloc(sizeof(WLNInstruction)); 
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

  /* gives linked list backtrack ability */
  void connect_instruction(WLNInstruction* parent, WLNInstruction *child)
  {
    parent->next_instructions.push_back(child);
    child->parent = parent; 
  }


  /* pops a number of rings off the stack */
  WLNInstruction* popdown_ringstack(std::stack<WLNInstruction*> &ring_stack, unsigned int terms){
    unsigned int popped = 0; 
    while(!ring_stack.empty()){
      ring_stack.pop();
      popped++; 

      if (terms == popped && !ring_stack.empty())
        return ring_stack.top();
    }
    return (WLNInstruction*)0;
  }


  WLNInstruction* backtrack_ringlinker(WLNInstruction* current)
  {
    // use the linked list
    while(current->parent){
      current = current->parent; 

      fprintf(stderr,"%s: %d\n",code_hierarchy[current->state],current->ring_linker);
      if(current->ring_linker)
        return current; 
    }

    return (WLNInstruction*)0;
  }


  /* parse the wln string and create instruction set,
  I think its reasonable to have two parses for this */
  bool CreateInstructionSet(const char *wln, unsigned int len){

    WLNInstruction *prev = 0; // use to add children for graph construct
    WLNInstruction *current = add_instruction(ROOT,0); 
    root = current;

    std::stack<WLNInstruction*> ring_stack; // keep

    bool pending_closure    = false; 
    bool pending_locant     = false; 
    bool pending_ring       = false;
  
  
    for (unsigned int i=0;i<len;i++){
      char ch = wln[i];
      switch(ch){
        
        case 'L':
        case 'T':
          pending_ring = false; // only stop pending once handled
          if(current->state == ROOT){
            prev = current;
            current = add_instruction(CYCLIC,i);
            ring_stack.push(current); // for back tracking if needed
            
            // update internal tracking
            pending_closure = true;

            connect_instruction(prev,current);
            break;
          }
          else if (current->state == STANDARD){
            
            if(pending_locant){
              current = add_instruction(LOCANT,i);
              current->add_end(i); // all locants terminate on one char
            
              pending_locant = false;

              if (!ring_stack.empty())
                ring_stack.top()->next_instructions.push_back(current);
              else{
                fprintf(stderr,"Error: no ring species to attach locant - terminating parse\n");
                return false;
              }  
            }
            
            break;              
          }
          else if(current->state == LOCANT){
            prev = current; 
            current = add_instruction(CYCLIC,i);
            ring_stack.push(current); // for back tracking if needed
          
            pending_closure = true;

            connect_instruction(prev,current);
            break;
          }
          else if (current->state == CYCLIC){
            
            if(pending_locant){
              prev = current; 

              current = add_instruction(LOCANT,i);
              current->add_end(i); // all locants terminate on one char
            
              pending_locant = false;

              connect_instruction(prev,current);
            }
            
            break;
          }
         

        case 'J': // pass 
          if (current->state == STANDARD){

            if(pending_locant){
              current = add_instruction(LOCANT,i);
              current->add_end(i); // all locants terminate on one char
              
              pending_locant = false;

              if (!ring_stack.empty())
                ring_stack.top()->next_instructions.push_back(current);
              else{
                // start some notation ending criteria
                fprintf(stderr,"Error: no ring species to attach locant - terminating parse\n");
                return false;
              }
            }

            break; 
          }
          else if (current->state == LOCANT || current->state == IONIC){
            prev = current; 
            current = add_instruction(STANDARD,i);
            connect_instruction(prev,current); 

            break;
          }
          else if (current->state == CYCLIC){
            
            if(pending_closure){
              current->add_end(i); // add the end and then update created
              pending_closure = false;  // allows internal atom positions to be passed
              
              // build the wln graph inline - greater flexibility for return backs
              
              // can do a ring check for poly - pericyclics (keep refactored)
              
              current->construct_standard_ring();
            }
            else if(pending_locant){
              prev = current; 

              current = add_instruction(LOCANT,i);
              current->add_end(i); // all locants terminate on one char
              
              pending_locant = false;

              connect_instruction(prev,current);
            }
            
            break;
          }
        
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
          pending_ring = false; // be very specific with pending ring

          if (current->state == STANDARD){
            
            if(pending_locant){
              prev = current; 
              current = add_instruction(LOCANT,i);
              current->add_end(i); // all locants terminate on one char
            
              pending_locant = false;
            }

            else if(pending_ring){
              connect_instruction(prev,current);
            }
            else{
              if (!ring_stack.empty())
                ring_stack.top()->next_instructions.push_back(current);
              else{
                fprintf(stderr,"Error: no ring species to attach locant - terminating parse\n");
                return false;
              }
            }
          
            break;
          }
          else if (current->state == LOCANT || current->state == IONIC){
            prev = current; 
            current = add_instruction(STANDARD,i);
            connect_instruction(prev,current); 
            break;
          } 
          else if (current->state == CYCLIC){

            if (pending_locant){
              prev = current; 

              current = add_instruction(LOCANT,i);
              current->add_end(i); // all locants terminate on one char
              
              pending_locant = false;

              connect_instruction(prev,current);
            }
            else{

              if(wln[i-1] == '&'){
                // attach a standard here, other rules handle ring, 
                // match on no symbol for easy merge.

                unsigned int k=i-1; 
                unsigned int terms = 0; 
                while(wln[k] == '&'){
                  terms++; 
                  k--;
                }
                // pop backs must be defined outer rule addition
                popdown_ringstack(ring_stack,terms); 

                current = backtrack_ringlinker(current);
                if(!current){
                  fprintf(stderr,"Error: no ring linker to return to via '&<x>-' - terminating parse\n");
                  return false; 
                }
                prev = current; 
                current = add_instruction(STANDARD,i);

                connect_instruction(prev,current);
              }

            }
            
            break;
          }
        

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
          pending_ring = false; // be very specific with pending ring

          if (current->state == ROOT){
            prev = current; 

            current = add_instruction(STANDARD,i);
             
            connect_instruction(prev,current);

            break;
          }
          else if (current->state == LOCANT || current->state == IONIC){
            prev = current; 
            current = add_instruction(STANDARD,i);
            connect_instruction(prev,current);
            
            break;
          }
          else if(current->state == CYCLIC){
            
            if(wln[i-1] == '&'){
              // attach a standard here, other rules handle ring, 
              // match on no symbol for easy merge.

              unsigned int k=i-1; 
              unsigned int terms = 0; 
              while(wln[k] == '&'){
                terms++; 
                k--;
              }
              // pop backs must be defined outer rule addition
              popdown_ringstack(ring_stack,terms); 

              current = backtrack_ringlinker(current);
              if(!current){
                fprintf(stderr,"Error: no ring linker to return to via '&<x>-' - terminating parse\n");
                return false; 
              }
              prev = current; 
              current = add_instruction(STANDARD,i);

              connect_instruction(prev,current);
            }
            
            break;
          }

        case ' ':  // keep this simple, a '&' locant means ionic branch out
          if (current->state == STANDARD){
            
            if (pending_ring){
              current->add_end(i-1);
              current->ring_linker = true; // use this for backtrack
        
              pending_locant = true;
            }
            else{
              current->add_end(i-1); // implied end of standard notation block

              // perform a check for '&' to look for ring burn condition
              if (wln[i-1] == '&'){

                fprintf(stderr,"char at pos: -%c-\n",wln[i]);
                // check how many rings are present behind
                unsigned int k=i-1; 
                unsigned int terms = 0; 
                while(wln[k] == '&'){
                  terms++; 
                  k--;
                }
                current = popdown_ringstack(ring_stack,terms);
                if (!current){
                  fprintf(stderr,"Error: notation contains too many '&', all rings popped - terminating parse\n");
                  return false; 
                }
              }
              pending_locant = true;
            }
            
            break;
          }
          
          else if (current->state == LOCANT){
            
            if(pending_ring){
              current = ring_stack.top();
              pending_locant = true; 
            }
            
            break;
          }

          else if (current->state == CYCLIC){
            
            if(!pending_closure)
              pending_locant = true;

            break;
          }
            
        
        case '-':
          pending_ring = false; // be very specific with pending ring

          if (current->state == ROOT){
            prev = current;
            current = add_instruction(STANDARD,i);
            connect_instruction(prev,current);
            break; 
          }
          else if(current->state == STANDARD){

            if (!ring_stack.empty())
              pending_ring = true; 
            
            break; 
          }
          else if (current->state == LOCANT){
            
            if (!ring_stack.empty())
              pending_ring = true;

            break;  
          }
          else if (current->state == CYCLIC){

            if(wln[i-1] == '&'){
              // means a return to a branching linker - also pops off ring
              unsigned int k=i-1; 
              unsigned int terms = 0; 
              while(wln[k] == '&'){
                terms++; 
                k--;
              }
              popdown_ringstack(ring_stack,terms); // we dont need current here;

              current = backtrack_ringlinker(current);
              if(!current){
                fprintf(stderr,"Error: no ring linker to return to via '&<x>-' - terminating parse\n");
                return false; 
              }
              prev = current; 
              current = add_instruction(STANDARD,i);
              connect_instruction(prev,current);

              pending_ring = true;               
            }

            break; 
          }

          else if (current->state == IONIC){
            prev = current; 
            current = add_instruction(STANDARD,i);
            connect_instruction(prev,current);

            break; 
          }

        case '&':
          if(current->state == STANDARD){
            
            if(pending_locant){
              current = add_instruction(IONIC,i); // ionic is always seperate in the graph
              current->add_end(i);

              // an ionic species means a complete blow through of the ring stack
              while(!ring_stack.empty())
                ring_stack.pop();

              pending_locant = false;
            }

            break;
          }
          else if (current->state == CYCLIC){
            
            
            if(pending_locant){
              current = add_instruction(IONIC,i); // ionic is always seperate in the graph
              current->add_end(i);

              // an ionic species means a complete blow through of the ring stack
              while(!ring_stack.empty())
                ring_stack.pop();

              pending_locant = false;
            }

            break;
          }
        
        default:
          fprintf(stderr,"Error: unrecognised symbol: %c\n",ch);
          return false; 
      }

    }

    // whatever the last instruction was, add len-1 as the end    
    current->add_end(len-1);
    return true; 
  }

  

  void DumpInstruction2Dot(FILE *fp, bool segment_string = false){
    // set up index map for node creation
    std::map <WLNInstruction*, unsigned int> index_map; 
    unsigned int glob_index = 0;
    for (WLNInstruction *node : instruction_pool){
      index_map[node] = glob_index; 
      glob_index++;
    }
      
    fprintf(fp,"digraph WLNdigraph {\n");
    fprintf(fp,"  rankdir = LR;\n");
    for (WLNInstruction *node : instruction_pool){
      fprintf(fp,"  %d",index_map[node]);
      if (segment_string){
        fprintf(fp, "[shape=circle,label=\"");
        for (unsigned int i=node->start_ch; i<=node->end_ch;i++)
          fprintf(fp,"%c",wln[i]);
        fprintf(fp,"\"];\n");
      }
      else
        fprintf(fp, "[shape=circle,label=\"%s\"];\n",code_hierarchy[node->state]);
      for (WLNInstruction *child : node->next_instructions){
        fprintf(fp,"  %d", index_map[node]);
        fprintf(fp," -> ");
        fprintf(fp,"%d\n",index_map[child]);
      }
    } 
    fprintf(fp,"}\n");
  }

};






static void DisplayUsage(){
  fprintf(stderr,"wln-writer <options> < input (escaped) >\n");
  fprintf(stderr,"<options>\n");
  fprintf(stderr,"  -v | --verbose                print messages to stdout\n");
  fprintf(stderr,"  -s | --strict                 fail on hypervalence, no symbol correction\n");
  fprintf(stderr,"  -c | --canonical              perform wln canonicalise procedure\n");
  fprintf(stderr,"  -r | --return-wln             return wln after altering procedure(s)\n");
  fprintf(stderr,"  --wln2dot                     dump wln trees to dot file\n");
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[]){

  const char *ptr=0; 
  int i,j; 

  wln     = (const char*)0;
  dotfile = (const char*)0;
  
  if (argc < 2)
   DisplayUsage();

  j=0; 
  for (i=1;i<argc;i++){

    ptr = argv[i];

    if (ptr[0]=='-' && ptr[1])
      switch (ptr[1]){

        case 'c':
          opt_canonical = true;
          break; 

        case 'r':
          opt_returnwln = true;
          break; 

        case 's':
          opt_valstrict = true; 
          break;

        case 'v':
          opt_verbose = true; 
          break; 

        case '-':
          if (!strcmp(ptr, "--wln2dot")){
            opt_wln2dot = true; 
            break;
          }
          else if (!strcmp(ptr, "--strict")){
            opt_valstrict = true; 
            break;
          }
          else if (!strcmp(ptr, "--verbose")){
            opt_verbose = true; 
            break;
          }
          else if (!strcmp(ptr, "--canonical")){
            opt_canonical = true; 
            break;
          }
          else if (!strcmp(ptr, "--return-wln")){
            opt_returnwln = true; 
            break;
          }

        default:
          fprintf(stderr,"Error: unrecognised input %s\n",ptr);
          break;
      }
    
    else switch(j++){
      case 0: wln = ptr; break;
      default: break;
    }
      
  }

  return;
}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
  if(!wln){
    fprintf(stderr,"Error: no wln string - nullptr\n");
    return 1;
  } 
  
  InstructionGraph parse_instructions;
  WLNGraph wln_graph;  

  if(!parse_instructions.CreateInstructionSet(wln,strlen(wln)))
    return 1; 
  
  if(opt_verbose)
    parse_instructions.display_instructions();

  
  // create the instruction graph
  if (opt_wln2dot){
    FILE *fp = 0;
    fp = fopen("instruction.dot","w");
    if(!fp){
      fprintf(stderr,"Error: coould not open compiler dump file\n");
      return 1;
    }
    else
      parse_instructions.DumpInstruction2Dot(fp,false);
  }

  

  return 0;
}