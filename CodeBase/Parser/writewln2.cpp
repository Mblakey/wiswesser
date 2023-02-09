
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector>
#include <stack>
#include <map>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

// forward dec
struct WLNSymbol; 

// --- inputs ---  
const char *wln; 
const char *dotfile; 


// --- options --- 
bool opt_wln2dot = false;
bool opt_valstrict = false; 
bool opt_verbose = false; 



// --- memory --- 
std::vector<WLNSymbol*> mempool; 

static void empty_mempool(){
  for (WLNSymbol* allocedwln : mempool){  // if error free all the mempool --> stop leak 
    free(allocedwln); 
    allocedwln = 0;
  }
}


enum WLNType {SINGLETON = 0, BRANCH = 1, LINKER = 2, TERMINATOR = 3}; 

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
      case ' ':
        type = TERMINATOR;
        allowed_edges = 1;
        break;

      case '-':
      case '/': 
        type = LINKER; 
        allowed_edges = 2;
        break;

      default:
        fprintf(stderr,"Error: invalid wln symbol parsed: %c\n",ch);
        return false; 
    }

    prev = (WLNSymbol*)0;
    return true; 
  } 

};


WLNSymbol* AllocateWLNSymbol(unsigned char ch){

  WLNSymbol *wln = (WLNSymbol*)malloc( sizeof(WLNSymbol) ); 
  if(wln->init(ch)) 
    mempool.push_back(wln);
  else 
    return (WLNSymbol *)0;
  
  return wln; 
}


bool handle_hypervalence(WLNSymbol *problem){

  // this will always lead to a positive ion species 
  switch(problem->ch){


    case 'M':   // tranforming this to a N is an easy solve
      if(opt_verbose)
        fprintf(stderr,"Status: transforming hypervalent M --> N\n");
      
      problem->ch = 'N';
      break; 

    case 'Y':  // can go to an X
      if (opt_verbose)
        fprintf(stderr,"Status: transforming hypervalent Y --> X\n");
      
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
  fprintf(stderr,"Error: returning nullptr from stack backtrack\n");
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
  fprintf(stderr,"Error: returning nullptr from force closure\n");
  return (WLNSymbol*)0;
}


/* parses NON CYCLIC input string, mallocs graph nodes and sets up graph based on symbol read */
WLNSymbol* ParseNonCyclic(const char *wln, unsigned int len){
  
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

    add_symbol(created_wln,prev);

    // options here depending on the forced closure of '&'
    if (created_wln->type == TERMINATOR){
      if (created_wln->ch == '&' && prev->type == BRANCH)
        prev = force_closure(wln_stack);
      else
        prev = backtrack_stack(wln_stack);
    }   
  }

  // for aid with tree reordering, a '&' can be used to close all wln notation
  
  prev = created_wln; // last created 
  created_wln = AllocateWLNSymbol('&');
  add_symbol(created_wln,prev);

  return root; // return start of the tree
}


/* dump wln tree to a dotvis file */
void WLNDumpToDot(FILE *fp){

  // set up index map for node creation
  std::map <WLNSymbol*, unsigned int> index_map; 
  unsigned int glob_index = 0;
  for (WLNSymbol *node : mempool){
    index_map[node] = glob_index; 
    glob_index++;
  }
    
  fprintf(fp,"digraph WLNdigraph {\n");
  fprintf(fp,"  rankdir = LR;\n");
  for (WLNSymbol *node : mempool){
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

static void DisplayUsage(){
  fprintf(stderr,"wln-writer <input> (escaped)\n");
  fprintf(stderr,"<options>\n");
  fprintf(stderr,"  -v | --verbose                print messages to stdout\n");
  fprintf(stderr,"  -s | --strict                 fail on hypervalence, no symbol correction\n");
  fprintf(stderr,"  --wln2dot <dotfile.dot>       dump wln tree to dot file\n");
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

        case 's':
          opt_valstrict = true; 
          break;

        case 'v':
          opt_verbose = true; 
          break; 

        case '-':
          if (!strcmp(ptr, "--wln2dot")){
            opt_wln2dot = true; 
            if (i == argc - 1){
              fprintf(stderr,"Error: --wlndot requires a <file>.dot as next arguement\n");
              DisplayUsage();
            }
            i++; // increment argv
            if (argv[i][0] != '-')
              dotfile = argv[i];
            else{
              fprintf(stderr,"Error: --wlndot requires a <file>.dot as next arguement\n");
              DisplayUsage();
            } 
            break;
          }

          if (!strcmp(ptr, "--strict")){
            opt_valstrict = true; 
            break;
          }
          
          if (!strcmp(ptr, "--verbose")){
            opt_verbose = true; 
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
      
    // end argc loop
  }

  return;
}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
  fprintf(stderr,"Parsing: %s\n",wln);

  WLNSymbol *root = ParseNonCyclic(wln, strlen(wln));
  if (!root){
    empty_mempool();
    return 1; 
  }

  if (opt_wln2dot){
    FILE *fp = 0; 
    fp = fopen(dotfile, "w");
    if (!fp)
      fprintf(stderr,"Error: could not write %s as .dot file - skipping\n",dotfile);
    else
      WLNDumpToDot(fp);
  }
  
  empty_mempool();
  return 0;
}