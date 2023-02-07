

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
        break;
      
      case 'A':
        type = LOCANT; 
        break;

      case 'B':
      case 'C':
        type = ATOM; 
        break;
      
      case 'D': 
        type = LOCANT; 
        break;
      
      case 'E':
      case 'F':
      case 'G': 
      case 'H': 
      case 'I': 
        type = ATOM; 
        break;

      case 'J':
        type = FRAGMENT; // generic symbol for halogen 
        break;
      
      case 'K': 
        type = ATOM; 
        break;

      case 'L':
        type = LINKER; 
        break;
      
      case 'M':
      case 'N': 
      case 'O':
      case 'P':
        type = ATOM; 
        break;

      case 'Q': 
      case 'R':
        type = FRAGMENT; 
        break;

      case 'S': 
        type = ATOM;
        break;

      case 'T':
      case 'U': 
        type = LINKER;
        break; 

      case 'V': 
        type = FRAGMENT; 
        break;

      
      case 'W':
      case 'X':
      case 'Y':
        type = LINKER; 
        break;

      case 'Z': 
        type = FRAGMENT; 
        break;

      case '&':
      case '-':
      case '/': 
        type = LINKER; 
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


/* parses input string, mallocs graph nodes and sets up graph based on symbol */
bool ParseWLN(const char *wln, unsigned int len, WLNGraph *symbol_tree){

  WLNSymbol *prev_seen = 0;
  for (unsigned int i = 0; i<len; i++){
    WLNSymbol* created_wln = AllocateWLNSymbol(wln[i]);
    if (!created_wln)
      return false; 

    if(!symbol_tree->head){
      symbol_tree->head = created_wln; // set up the first access point
      prev_seen = created_wln;
      continue;
    }

    created_wln->prev = prev_seen; // update previous (this is simple)


    
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