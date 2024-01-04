#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <map>

#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"

unsigned saved_bytes = 0;
unsigned int opt_mode = 0;
unsigned int opt_verbose = false;
const char *input;


bool train_on_file(FILE *ifp, FSMAutomata *wlnmodel){
  
  unsigned char ch = 0;
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;
  
  while(fread(&ch, sizeof(unsigned char), 1, ifp)){
    for(edge=curr->transitions;edge;edge=edge->nxt){
      if(edge->ch == ch){
        edge->c++; 
        curr = edge->dwn; 
        break;
      }
    }
  }

  // write 32 bit numbers in index order 
  for(unsigned int i=0;i< wlnmodel->num_edges;i++){
    unsigned int edge_freq = wlnmodel->edges[i]->c;
    fwrite(&edge_freq,sizeof(unsigned int),1,stdout);
  }

  return true;
}


static void DisplayUsage()
{
  fprintf(stderr, "wlntrain\n");
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
        case 'a':
          opt_mode = 1;
          break;

        case 'h':
          opt_mode = 2;
          break;
        
        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
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

  return;
}


int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);

  FSMAutomata *wlnmodel = CreateWLNDFA(); // build the model 

  // model arithmetic coder
  wlnmodel->AddTransition(wlnmodel->root,wlnmodel->root,'\0');  
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept)
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
  }

  unsigned int singles = 0;
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->transitions && !wlnmodel->states[i]->transitions->nxt)
    {
      singles++;
    }
  } 
 
  fprintf(stderr,"wln has %d singles\n",singles);

  delete wlnmodel;
  return 0;
}