
#include <stdlib.h>
#include <stdio.h> 

#include "rfsm.h"
#include "wlndfa.h"
#include "wlnzip.h"

const char *input;

static void DisplayUsage()
{
  fprintf(stderr, "wlnzip <options> <input> > <out>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -c          compress input\n");
  fprintf(stderr, "  -d          decompress input\n");
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  input = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          exit(1); 
      }
    }
    else{
      switch(j++){
        case 0:
          input = ptr; 
          break;
        default:
          fprintf(stderr,"Error: multiple files not currently supported\n");
          exit(1);
      }
    }
  }

  if(!input){
    fprintf(stderr,"Error: no input file given\n");
    exit(1);
  }

  return;
}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
    
  FSMAutomata *wlnmodel = CreateWLNDFA(REASONABLE*2,REASONABLE*4); // build the model 

  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept){
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\0');  
    }
  }
  wlnmodel->AssignEqualProbs(); // initalise the order -1 model
  

  bool ending = true;
  // perhaps need to add the escape sequence as a transition. 
  std::string bitstream; 
  if(!WLNPPMCompressBuffer(input, wlnmodel, bitstream,ending))
    return 1; 
   
  if(ending)
    WLNPPMDecompressBuffer(bitstream, wlnmodel);
  
  delete wlnmodel;
  return 0;
}
