#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <random>

#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"

unsigned int gen_length = 5;
unsigned int gen_count = 10; 
unsigned int opt_verbose = false;

const char *trainfile;


bool train_on_file(FILE *ifp, FSMAutomata *wlnmodel){
  
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

  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    FSMState *s = wlnmodel->states[i];
    unsigned int count = 0;
    FSMEdge *e = 0; 
    for(e=s->transitions;e;e=e->nxt)
      count+= e->c; // get total counts

    for(e=s->transitions;e;e=e->nxt)
      e->p = (double)e->c/(double)count;
  }
  return true;
}

bool StatisticalGenerate(FSMAutomata *wlnmodel){
  unsigned int count; 
  unsigned int length; 

  std::random_device rd;
  std::mt19937 gen(rd());

  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  while(count < gen_count){

    std::vector<FSMEdge*> e = {};
    std::vector<double> p = {}; // vector of probabilities 
    for(edge=state->transitions;edge;edge=edge->nxt){
      e.push_back(edge);
      p.push_back(edge->p);
    }

    std::discrete_distribution<> d(p.begin(), p.end());
    unsigned int chosen = d(gen); 
    if(e[chosen]->ch == '\n'){

      if(gen_length <= length){ // only accepts will have new line, ensures proper molecule
        fputc('\n',stdout);
        count++;
        length = 0;
        state = wlnmodel->root;
      }
      else{
        // choose something else
        while(e[chosen]->ch == '\n')
          chosen = d(gen); 
      }
    }

    fputc(e[chosen]->ch,stdout);
    length++;

    state = e[chosen]->dwn;
  }


  return true;
}


static void DisplayUsage()
{
  fprintf(stderr, "wlngen <options> <trainfile>\n");
  fprintf(stderr,"options:\n");
  fprintf(stderr,"-l <int>    set length for generation (default 5)\n");
  fprintf(stderr,"-n <int>    set target count for generation (default 10)\n");
  exit(1);
}


static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  trainfile = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {
    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        case 'n':
          if(i > argc-1){
            fprintf(stderr,"Error: count flag must be followed by a number\n");
            DisplayUsage();
          }
          else if(!sscanf(argv[i+1],"%d",&gen_count)){
            fprintf(stderr,"Error: count flag must be followed by a number\n");
            DisplayUsage();
          }
          i++;
          break;

        case 'l':
          if(i > argc-1){
            fprintf(stderr,"Error: length flag must be followed by a number\n");
            DisplayUsage();
          }
          else if(!sscanf(argv[i+1],"%d",&gen_length)){
            fprintf(stderr,"Error: length flag must be followed by a number\n");
            DisplayUsage();
          }
          i++;
          break;
        
        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
    else{
      switch(j++){
        case 0:
          trainfile = ptr; 
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

  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept)
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
  }

  wlnmodel->AssignEqualProbs(); // default and prevents div 0

  if(trainfile){
    FILE *tfp = fopen(trainfile,"r");
    if(!tfp){
      fprintf(stderr,"Error: could not open train file\n");
      return 1;
    }

    train_on_file(tfp,wlnmodel);
    fclose(tfp);
  }


  StatisticalGenerate(wlnmodel);

  delete wlnmodel;
  return 0;
}