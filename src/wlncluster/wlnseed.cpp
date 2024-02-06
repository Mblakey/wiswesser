

#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <vector>

#include "wlndfa.h"
#include "rfsm.h"

std::vector<const char*> train_files; 


bool seed_from_file(FILE *ifp, FSMAutomata *wlnmodel){
  unsigned char ch = 0;
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;
  
  while(fread(&ch, sizeof(unsigned char), 1, ifp)){
    for(edge=curr->transitions;edge;edge=edge->nxt){
      if(edge->ch == ch){
        if(edge->c < UINT32_MAX)
          edge->c++; 

        curr = edge->dwn; 
        break;
      }
    }
  }

  return true;
}


bool WriteEdgeCounts(FSMAutomata *wlnmodel){
  for(unsigned int i=0;i<wlnmodel->num_edges;i++){
    fprintf(stdout,"%d\n",wlnmodel->edges[i]->c);
  }
  return true;
}


static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;
  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        default:
          fprintf(stderr, "Error: unrecognised option %s\n", ptr);
          exit(1);
      }
    }
    else{
      switch(j++){
        default:
          train_files.push_back(ptr); 
          break;
      }
    }
  }

  return;
}




int main(int argc, char *argv[])
{

  ProcessCommandLine(argc, argv);

  FSMAutomata *wlnmodel = CreateWLNDFA(REASONABLE,REASONABLE,false); // build the model 

  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept)
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
  }

  unsigned int f = 0;
  if(!train_files.empty()){
    for (const char *trainfile : train_files){
      f++;
      FILE *tfp = fopen(trainfile,"r");
      if(!tfp){
        fprintf(stderr,"Error: could not open train file %d - skipping\n",f);
        continue;
      }
      else{
        seed_from_file(tfp,wlnmodel);
        fclose(tfp);
        tfp = 0;
      }
    }
  }
  else {
    fprintf(stderr,"Error: no files provided\n"); 
    return 1;
  }
  
  WriteEdgeCounts(wlnmodel);
  delete wlnmodel;
  return 0;
}
