

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


struct WLNSymbol;  
struct WLNGraph; 

const char *wln; 
std::vector<WLNSymbol*> mempool; 
enum WLNType { UNRESOLVED = 0, SINGLETON = 1, BRANCH = 2, LINKER = 3, TERMINATOR = 4}; 


struct WLNGraph{
  WLNSymbol *head; // linked list beg
  WLNSymbol *tail; // linked list end
  unsigned int count; 

  WLNGraph(){
    head = (WLNSymbol*)0;
    tail = (WLNSymbol*)0;
    count = 0; 
  }

}; 

static void empty_mempool(){
  for (WLNSymbol* allocedwln : mempool){  // if error free all the mempool --> stop leak 
    free(allocedwln); 
    allocedwln = 0;
  }
}

struct WLNSymbol{

  unsigned char ch; 
  unsigned int type; 
  unsigned int num_childs;  

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
        num_childs = 1;
        break;
      
      case 'A':
        type = SINGLETON;
        num_childs = 1;
        break;

      case 'B':   // boron
        type = BRANCH;
        num_childs = 2; 
        break;
      
       case 'C':  // shortcut carbon atom
        type = BRANCH;
        num_childs = 3;
        break;
      
      case 'D': 
        type = SINGLETON; 
        num_childs = 1;
        break;
      
      case 'E':  // halogens
      case 'F':
      case 'G': 
      case 'I': 
        type = BRANCH; 
        num_childs = 2; 
        break;
      
      case 'H': // closing hydrogen
        type = TERMINATOR; 
        num_childs = 0; 
        break;

      case 'J': // generic symbol for halogen 
        type = BRANCH; 
        num_childs = 2; 
        break;
      
      case 'K': 
        type = BRANCH; 
        num_childs = 3; 
        break;

      case 'L':
        type = LINKER;
        num_childs = 1; 
        break;
      
      case 'M':
        type = BRANCH; 
        num_childs = 1; 
        break;

      case 'N': 
        type = BRANCH; 
        num_childs = 2;
        break; 

      case 'O':
        type = BRANCH; 
        num_childs = 1;
        break;

      case 'P':
        type = BRANCH; 
        num_childs = 4;
        break;

      case 'Q': 
        type = TERMINATOR; 
        num_childs = 0;
        break;

      case 'R':
        type = SINGLETON; 
        num_childs = 0; 
        break;

      case 'S': 
        type = BRANCH;
        num_childs = 5;
        break;

      case 'T':
      case 'U': 
        type = LINKER;
        num_childs = 1;
        break; 

      case 'V': 
        type = SINGLETON; 
        num_childs = 1;
        break;

      case 'W':
        type = LINKER; 
        num_childs = 1;
        break;

      case 'X':
        type = BRANCH; 
        num_childs = 4;
        break;

      case 'Y':
        type = BRANCH; 
        num_childs = 3;
        break;

      case 'Z': 
        type = TERMINATOR;
        num_childs = 0; 
        break;

      case '&':
        type = TERMINATOR;
        num_childs = 0; 
        break;

      case '-':
      case '/': 
        type = LINKER; 
        num_childs = 1; 
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


/* create the functions used for adding the WLN tree */



/* add src to the children vector of trg */
void add_symbol(WLNSymbol* src, WLNSymbol *trg){

#ifdef DEBUGWLN
  if (trg->children.size() < trg->num_childs){
    fprintf(stderr,"adding symbol %c to bonds of %c\n",src->ch, trg->ch);
    trg->children.push_back(src);
  }
  else {
    fprintf(stderr,"Warning: allowing hypervalence on WLN character %c\n",trg->ch);
    fprintf(stderr,"adding symbol %c to bonds of %c\n",src->ch, trg->ch);
    trg->children.push_back(src);
  }   
#else
  if (trg->children.size() < trg->num_childs){
    trg->children.push_back(src);
  }
  else {
    fprintf(stderr,"Warning: allowing hypervalence on WLN character %c\n",trg->ch);
    trg->children.push_back(src);
  }   
#endif

}



/* parses NON CYCLIC input string, mallocs graph nodes and sets up graph based on symbol read */
bool ParseNonCyclic(const char *wln, unsigned int len, WLNGraph *symbol_tree){
  
  WLNSymbol *prev = 0; // hold previous term for single linear adding 
  std::stack <WLNSymbol*> wln_stack; 

  for (unsigned int i = 0; i<len; i++){

    WLNSymbol* created_wln = AllocateWLNSymbol(wln[i]);
    if (!created_wln)
      return false; 

    if(!symbol_tree->head){
      symbol_tree->head = created_wln; // set up the first access point
      if (created_wln->type == BRANCH)
        wln_stack.push(created_wln);
      prev = created_wln; 
    }
    else {
      add_symbol(created_wln,prev);
      if (created_wln->type == TERMINATOR && !wln_stack.empty())
        prev = wln_stack.top();
      else if (created_wln->type == BRANCH){
        wln_stack.push(created_wln);
        prev = created_wln;
      }else
        prev = created_wln; 
    }
  }
  return true;
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
  fprintf(stderr,"--wln2dot <dotfile.dot>       dump wln tree to dot file\n")
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[]){

  const char *ptr=0; 
  int i,j; 

  wln = (const char*)0; 
  
  if (argc < 2)
   DisplayUsage();

  j=0; 
  for (i=1;i<argc;i++){

    ptr = argv[i];

    if (ptr[0]=='-' && ptr[1])
      fprintf(stderr,"Error: writer only takes in single strings, option detected!\n");
    
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

  WLNGraph *symbol_tree = new WLNGraph();  // defaults set

  
  if(!ParseNonCyclic(wln, strlen(wln), symbol_tree)){
    empty_mempool();
    delete symbol_tree; 
    return 1; 
  }

  WLNDumpToDot(stdout);

  empty_mempool();
  delete symbol_tree; 
  return 0;
}