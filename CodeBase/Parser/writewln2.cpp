

#include <stdlib.h>
#include <stdio.h>

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


const char *wln; 


// Carbon decomp can be done at init time 


struct WLNSymbol{

  enum type { CARBON = 0, ATOM = 1, FRAGMENT = 2, LINKER = 3}; 
  unsigned char ch; 

  WLNSymbol *dn; // should be a single term
  WLNSymbol *ac; // linked list of across chains 







}; 



static void DisplayUsage(){
  fprintf(stderr,"wln-writer <input> (escaped)\n");
  exit(1);
}


static unsigned int ProcessCommandLine(int argc, char *argv[]){

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

}



int main(int argc, char *argv[]){


  ProcessCommandLine(argc, argv);

  fprintf(stderr,"wln: %s\n",wln);

  return 0;
}