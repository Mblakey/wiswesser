

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector>
#include <stack>

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

  enum WLNType { UNRESOLVED = 0, CARBON = 1, ATOM = 2, FRAGMENT = 3, LINKER = 4, LOCANT = 5}; 

  unsigned char ch; 
  unsigned int type; 
  unsigned int max_next_size;  

  WLNSymbol *prev; // should be a single term - wln symbol only has one incoming
  std::vector<WLNSymbol*> next; // linked list of next terms chains 

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
        type = CARBON; 
        max_next_size = 1;
        break;
      
      case 'A':
        type = LOCANT; 
        max_next_size = 1;
        break;

      case 'B':   // boron
        type = ATOM;
        max_next_size = 2;
        break;
      
       case 'C':  // shortcut carbon atom
        type = ATOM;
        max_next_size = 3;
        break;
      
      case 'D': 
        type = LOCANT;
        max_next_size = 1; 
        break;
      
      case 'E':  // halogens
      case 'F':
      case 'G': 
      case 'I': 
        type = ATOM; 
        max_next_size = 2;
        break;
      
      case 'H': // closing hydrogen
        type = ATOM; 
        max_next_size = 0;
        break;

      case 'J': // generic symbol for halogen 
        type = ATOM; 
        max_next_size = 2;
        break;
      
      case 'K': 
        type = ATOM; 
        break;

      case 'L':
        type = LINKER;
        max_next_size = 1; 
        break;
      
      case 'M':
        type = ATOM; 
        max_next_size = 1;
        break;

      case 'N': 
        type = ATOM; 
        max_next_size = 2;
        break; 

      case 'O':
        type = ATOM; 
        max_next_size = 1;
        break;

      case 'P':
        type = ATOM; 
        max_next_size = 4;
        break;

      case 'Q': 
        type = FRAGMENT; 
        max_next_size = 1;
        break;

      case 'R':
        type = FRAGMENT; 
        max_next_size = 0; 
        break;

      case 'S': 
        type = ATOM;
        max_next_size = 5;
        break;

      case 'T':
      case 'U': 
        type = LINKER;
        max_next_size = 1;
        break; 

      case 'V': 
        type = FRAGMENT; 
        max_next_size = 1;
        break;

      case 'W':
        type = LINKER; 
        max_next_size = 1;
        break;

      case 'X':
        type = LINKER; 
        max_next_size = 4;
        break;

      case 'Y':
        type = LINKER; 
        max_next_size = 3;
        break;

      case 'Z': 
        type = ATOM;
        max_next_size = 0; 
        break;

      case '&':
      case '-':
      case '/': 
        type = LINKER; 
        max_next_size = 1; 
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



/* add src to the next vector of trg */
void add_symbol(WLNSymbol* src, WLNSymbol *trg){

#ifdef DEBUGWLN
  if (trg->next.size() < trg->max_next_size){
    fprintf(stderr,"adding symbol %c to bonds of %c\n",src->ch, trg->ch);
    trg->next.push_back(src);
  }
  else {
    fprintf(stderr,"Warning: allowing hypervalence on WLN character %c\n",trg->ch);
    trg->next.push_back(src);
  }   

#else
  if (trg->next.size() < trg->max_next_size){
    trg->next.push_back(src);
  }
  else {
    fprintf(stderr,"Warning: allowing hypervalence on WLN character %c\n",trg->ch);
    trg->next.push_back(src);
  }   
#endif

}



/* parses input string, mallocs graph nodes and sets up graph based on symbol */
bool ParseWLN(const char *wln, unsigned int len, WLNGraph *symbol_tree){


  std::stack <WLNSymbol*> wln_stack; 
  WLNSymbol *top = 0; 
  for (unsigned int i = 0; i<len; i++){

    WLNSymbol* created_wln = AllocateWLNSymbol(wln[i]);
    if (!created_wln)
      return false; 

    if(!symbol_tree->head){
      symbol_tree->head = created_wln; // set up the first access point
      wln_stack.push(created_wln);
    }
    else {
      
      top = wln_stack.top(); 
      
      // we always add into the graph node of top, and we move top around instead

      add_symbol(created_wln,top);

      if (top->next.size() < top->max_next_size)
        wln_stack.pop();

      wln_stack.push(created_wln);


    }
    
  }

  return true;
}




static void DisplayUsage(){
  fprintf(stderr,"wln-writer <input> (escaped)\n");
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

  
  if(!ParseWLN(wln, strlen(wln), symbol_tree)){
    empty_mempool();
    delete symbol_tree; 
    return 1; 
  }


  

  empty_mempool();
  delete symbol_tree; 
  return 0;
}