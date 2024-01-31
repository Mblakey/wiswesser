#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <random>

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/kekulize.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/groupcontrib.h>

#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"
#include "parser.h"
#include "descriptor.h"

unsigned int gen_length = 5;
unsigned int gen_count = 10; 
unsigned int opt_verbose = false;

const char *trainfile;


using namespace OpenBabel;


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


bool Validate(const char * wln_str, OBMol *mol){
  fprintf(stderr,"testing: %s\n",wln_str);
  if(!ReadWLN(wln_str,mol))
    return false;
  return true;
}

double RewardFunction(const char *wln_strm, OBMol *mol){
  OBDescriptor* pDesc = OBDescriptor::FindType("logP");
  if(pDesc)
    return pDesc->Predict(mol); 
  else
    return 0.0;
};


/* 
--- notes ---

Epsilon-greedy policy - either take an action based on the qtable
values, (exploitation), or the random FSM values, (exploration).

I think really the Q table needs to be on the edges per state, since
two transitions could go to the same state, but have wildly different
outcomes to the final string. 

logic might follow: mealey Q-learning machine? hmmm

dummy scoring:

+1 = valid WLN
+2 = unique
+5 = in target range

From a good Q table, rescale the probabilities, some might be zero which
means transitions are removed? Generate compounds


Mode collapse - not a black box model so there will be methods to handle 
                this
*/


/* watch for overflows on the unsigned int, floats MAY be better 
for incremental rewards that level off */
void debug_Qtable(FSMAutomata *wlnmodel, std::map<FSMEdge*,unsigned int> &Qtable){
  for(unsigned int i=0;i<wlnmodel->num_edges;i++){
    FSMEdge *e = wlnmodel->edges[i]; 
    fprintf(stderr,"%d - %d\n",e->id,Qtable[e]);
  }
}

/* uses Q learning to generate compounds from the language FSM
as a markov decision process, WLN is small enough that with 20
characters a large scope of chemical space can be covered. */
bool QGenerateWLN(FSMAutomata *wlnmodel){
  unsigned int count = 0; 
  unsigned int length = 0; 

  // set up a qtable for each state
  std::map<FSMEdge*,unsigned int> Qtable; 
  for(unsigned int i=0;i<wlnmodel->num_edges;i++)
    Qtable[wlnmodel->edges[i]] = 0;

  
  std::random_device rd;
  std::mt19937 gen(rd());

  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  std::string wlnstr; 
  std::vector<FSMEdge*> path;
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
        length = 0;
        state = wlnmodel->root;

        // in beta, the faster i make readWLN the better this is
        
        OBMol mol;
        if(Validate(wlnstr.c_str(),&mol)){
          count++;
          //fprintf(stdout,"%s = %f\n",wlnstr.c_str(),RewardFunction(wlnstr.c_str(),&mol));
        
          // go back through all the edges and give them the +1 score
          for(FSMEdge *e:path){
            Qtable[e]++;
          }
        
        }
    
        
        wlnstr.clear();
        path.clear(); // can assign learning here
      }
      else{
        // choose something else
        while(e[chosen]->ch == '\n')
          chosen = d(gen); 
      }
    }
    else{
      path.push_back(e[chosen]);
      wlnstr += e[chosen]->ch;
      length++;
    }
   

    state = e[chosen]->dwn;
  }

  //debug_Qtable(wlnmodel,Qtable);
  return true;
}


static void DisplayUsage()
{
  fprintf(stderr, "wlngen <options> <trainfile>\n");
  fprintf(stderr,"options:\n");
  fprintf(stderr,"-l <int>    set length for generation (default 5)\n");
  fprintf(stderr,"-n <int>    set target count for generation (default 10)\n");
  fprintf(stderr,"-v          timing and debugging statements to console\n");
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
        
        case 'v':
          opt_verbose = true;
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

  FSMAutomata *wlnmodel = CreateWLNDFA(REASONABLE,REASONABLE,false); // build the model 

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

  QGenerateWLN(wlnmodel);
  delete wlnmodel;
  return 0;
}