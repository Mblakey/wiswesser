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
  fprintf(stderr, "wlntrain <input> <type> <out>\n");
  fprintf(stderr,"types:\n");
  fprintf(stderr,"-a    create train file for arthimetic coder (wlncompress)\n");
  fprintf(stderr,"-h    create train file for huffman coder    (wlncompress2)\n");
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

  if(!input){
    fprintf(stderr,"Error: no input file given\n");
    DisplayUsage();
  }

  if(!opt_mode){
    fprintf(stderr,"Error: no choice for type of training file selected\n");
    DisplayUsage();
  }

  return;
}


int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);

  FSMAutomata *wlnmodel = CreateWLNDFA(); // build the model 

  // make the root an EOF 
  wlnmodel->AddTransition(wlnmodel->root,wlnmodel->root,'\0');

  if(opt_mode == 1){
  
    // to every accept add the newline character pointing back to the root
    for(unsigned int i=0;i<wlnmodel->num_states;i++){
      if(wlnmodel->states[i]->accept)
        wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
    }
  }
  else if (opt_mode == 2){
    // add a null byte to each accept only
    for(unsigned int i=0;i<wlnmodel->num_states;i++){
      if(wlnmodel->states[i]->accept)
        wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->states[i],'\0');
    }
  }

  // add 1 to avoid the zero frequency problem
  for(unsigned int i=0;i< wlnmodel->num_edges;i++)
    wlnmodel->edges[i]->c = 1;


  FILE *fp = 0; 
  fp = fopen(input,"rb");
  if(fp){
    train_on_file(fp,wlnmodel);
    fclose(fp);
  }
  else{
    fprintf(stderr,"Error: could not open file at %s\n",input);
    return 1;
  }

  delete wlnmodel;
  return 0;
}