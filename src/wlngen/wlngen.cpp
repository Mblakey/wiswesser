#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <set>
#include <unordered_map> // traditional hash map
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
#include "wlndfa.h"
#include "parser.h"

using namespace OpenBabel;

int length = 5;
int count = 10; 
int episodes = 5;

int dt = 0;
double molwt = 0.0;
double logp = 0.0;

double epsilon = 0.5;
double learning_rate = 0.5;
double discount_rate = 0.85;
double decay_rate = 0.005;

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


/* 
--- notes ---

Previous attempts at the classic Q learning approach does not work, rewards 
are way to sparse and time delaying also doesnt work. The main one is the greedy here
is not suitable since mulitple transitions can take each state. 

Each HER and time delayed measures to bring this into line.

Epsilon-greedy policy - either take an action based on the qtable
values, (exploitation), or the random FSM values, (exploration).

dummy scoring:

+1 = valid WLN
+2 = unique
+5 = in target range

From a good Q table, rescale the probabilities, some might be zero which
means transitions are removed? Generate compounds
*/


bool Validate(const char *wln_str, OBMol *mol){
  if(!ReadWLN(wln_str,mol))
    return false;
  return true;
}


// https://open-babel.readthedocs.io/en/latest/Descriptors/descriptors.html


double LogP(const char *wln_strm, OBMol *mol){
  OBDescriptor* pDesc = OBDescriptor::FindType("logP");
  if(pDesc)
    return pDesc->Predict(mol); 
  else
    return 0.0;
}

double MolWt(const char *wln_strm, OBMol *mol){
  OBDescriptor* pDesc = OBDescriptor::FindType("MW");
  if(pDesc)
    return pDesc->Predict(mol); 
  else
    return 0.0;
}

#if LEARNING
FSMEdge *EpsilonGreedy(FSMState *curr, double epsilon, std::mt19937 &rgen){

  std::uniform_real_distribution<> dis(0, 1);
  double choice = dis(rgen); 

  FSMEdge *e = 0;
  if(choice > epsilon){
    
    // Explotation, take the best Q score
    FSMEdge *best = curr->transitions;
    for(e = curr->transitions;e;e=e->nxt){
      if(e->c > best->c)
        best = e; 
    }

    return best; 
  }
  else{
    std::vector<FSMEdge*> edge_vec = {};
    std::vector<double> prob_vec = {}; // vector of probabilities 
    for(e=curr->transitions;e;e=e->nxt){
      edge_vec.push_back(e);
      prob_vec.push_back(e->p);
    }

    std::discrete_distribution<> d(prob_vec.begin(), prob_vec.end());
    unsigned int chosen = d(rgen);
    return edge_vec[chosen];
  }
}
#endif 

FSMEdge *RandomEdge(FSMState *curr, std::mt19937 &rgen){
  std::uniform_real_distribution<> dis(0, 1);
  std::vector<FSMEdge*> edge_vec = {};
  std::vector<double> prob_vec = {}; // vector of probabilities 
  
  FSMEdge *e = 0;
  for(e=curr->transitions;e;e=e->nxt){
    edge_vec.push_back(e);
    prob_vec.push_back(e->p);
  }

  std::discrete_distribution<> d(prob_vec.begin(), prob_vec.end());
  unsigned int chosen = d(rgen);
  return edge_vec[chosen];
}

/* uses Q learning to generate compounds from the language FSM
as a markov decision process, WLN is small enough that with 20
characters a large scope of chemical space can be covered. */
bool QGenerateWLN(FSMAutomata *wlnmodel){
  int hits = 0; 
  int misses = 0;
  int duplicates = 0;
  int out_range = 0;
 
  // from the compression data, we can use edge->c to be the q-scoring

  std::random_device rd;
  std::mt19937 gen(rd());

  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  std::string wlnstr; 
  std::set<FSMEdge*> path; // avoid the duplicate path increases
  std::unordered_map<std::string, bool> unique;

  int strlength = 0;
  while(hits < count){

    edge = RandomEdge(state,gen);

    if(edge->ch == '\n'){

      if(strlength >= length){ // only accepts will have new line, ensures proper molecule
        strlength = 0;
        state = wlnmodel->root;

        // in beta, the faster i make ReadWLN the better this is
        
        OBMol mol;
        unsigned int score = 0;
        if(Validate(wlnstr.c_str(),&mol)){
          score+= 1;

          if(!unique[wlnstr]){
            score += 1;
            unique[wlnstr] = true;

            double lp = 0.0;
            double mw = 0.0;
            switch(dt){
              case 0:
                hits++;
                fprintf(stderr,"%s\n",wlnstr.c_str());
                break;

              case 1:
                lp = LogP(wlnstr.c_str(),&mol);
                if (lp >= logp-0.5 && lp <= logp+0.5){
                  score+= 3;
                  hits++;
                  fprintf(stderr,"%s - %f\n",wlnstr.c_str(),lp);
                }
                else
                  out_range++;
                break;

              case 2: 
                mw = MolWt(wlnstr.c_str(),&mol);
                if (mw >= molwt-50 && mw <= molwt+50){
                  score+= 3;
                  hits++;
                  fprintf(stderr,"%s - %f\n",wlnstr.c_str(),mw);
                }
                else
                  out_range++;
                break;
            }
          }
          else 
            duplicates++;
        
          // go back through all the edges and give them the +1 score
          if(score){
            path.insert(edge);
            for(FSMEdge *pe:path){
              if(pe->c < UINT32_MAX)
                pe->c+=score;
            }
          }
        }
        else
          misses++;
    
        wlnstr.clear();
        path.clear(); // can assign learning here
      }
      else{
        // choose something else
        while(edge->ch == '\n')
          edge = RandomEdge(state,gen);

        path.insert(edge);
      }
    }
    else{
      path.insert(edge);
      wlnstr += edge->ch;
      strlength++;
    }
   
    state = edge->dwn;
  }

  fprintf(stderr,"%d hits, %d misses, %d duplicates, %d out of target range\n",hits,misses,duplicates,out_range);
  return true;
}


bool prefix(const char *pre, const char *str)
{
  return strncmp(pre, str, strlen(pre)) == 0;
}

static void DisplayUsage()
{
  fprintf(stderr, "wlngen <options> <trainfile>\n");
  fprintf(stderr,"options:\n");
  fprintf(stderr,"-l|--length=<int>      set length for generation        (default 5)\n");
  fprintf(stderr,"-c|--count=<int>       set target count for generation  (default 10)\n");
  
  fprintf(stderr,"\ntuning:\n");
  fprintf(stderr,"-r|--runs=<int>        set learning episodes            (default 5)\n");
  fprintf(stderr,"-e|--epsilon=<double>  set epsilon hyperparameter       (default 0.5)\n");
  fprintf(stderr,"-d|--decay=<double>    set decay rate hyperparameter    (default 0.005)\n");
  fprintf(stderr,"-a|--alpha=<double>    set learning rate hyperparameter (default 0.5)\n");
  fprintf(stderr,"-g|--gamma=<double>    set discount rate hyperparameter (default 0.85)\n");
  
  fprintf(stderr,"\ngeneral:\n");
  fprintf(stderr,"-p|--print             show all set hyperparameters and exit\n");
  fprintf(stderr,"-h|--help              show this help menu and exit\n");
  
  fprintf(stderr,"\ndescriptors:\n");
  fprintf(stderr,"--logp=<double>        set logp  target value, range is +/- 0.5 from this value\n");
  fprintf(stderr,"--molwt=<double>       set molwt target value, range is +/- 50  from this value\n");
  exit(1);
}


static void DisplayParameters(){
  fprintf(stderr,"----------------------------\n");
  fprintf(stderr,"target count:      %d\n",count);
  fprintf(stderr,"target length:     %d\n",length);
  
  fprintf(stderr,"\nepisodes:          %d\n",episodes);
  fprintf(stderr,"learning rate:     %f\n",learning_rate);
  fprintf(stderr,"discount rate:     %f\n",discount_rate);
  fprintf(stderr,"epsilon:           %f\n",epsilon);
  fprintf(stderr,"decay rate:        %f\n",decay_rate);
  
  fprintf(stderr,"\nlogp target:       %f\n",logp);
  fprintf(stderr,"molwt target:      %f\n",molwt);
  switch(dt){
    case 0:
      fprintf(stderr,"dt mode:           %d (no descriptor)\n",dt);
      break;
    case 1:
      fprintf(stderr,"dt mode:           %d (logp)\n",dt);
      break;
    case 2:
      fprintf(stderr,"dt mode:           %d (molwt)\n",dt);
      break;
  }
  fprintf(stderr,"----------------------------\n");
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  unsigned char ch = 0;
  int i,j,l;
  j = 0;
  for (i = 1; i < argc; i++)
  {
    ptr = argv[i];
    ch = *ptr; 
    if (ptr[0] == '-' && ptr[1]){
      l = 0;
      switch (ptr[1]){

        case 'p':
          DisplayParameters();
          break;
         

        case 'h':
          DisplayUsage();
          break;


        case 'r':
          while(ch != '=' && ch){
            ch = *(++ptr);
            l++;
          }

          if(!(*ptr) || l != 2){
            fprintf(stderr,"Error: incorrect flag format\n");
            DisplayUsage();
          }
          else{
            ptr++;
            episodes = atoi(ptr);
            if(episodes < 0){
              fprintf(stderr,"Error: runs must be a positive integer\n");
              DisplayUsage();
            }
          }
          break;


        case 'c':
          while(ch != '=' && ch){
            ch = *(++ptr);
            l++;
          }

          if(!(*ptr) || l != 2){
            fprintf(stderr,"Error: incorrect flag format\n");
            DisplayUsage();
          }
          else{
            ptr++;
            count = atoi(ptr);
            if(count < 0){
              fprintf(stderr,"Error: count must be a positive integer\n");
              DisplayUsage();
            }
          }
          break;

        case 'l':
          while(ch != '=' && ch){
            ch = *(++ptr);
            l++;
          }

          if(!(*ptr) || l != 2){
            fprintf(stderr,"Error: incorrect flag format\n");
            DisplayUsage();
          }
          else{
            ptr++;
            length = atoi(ptr);
            if(length < 0){
              fprintf(stderr,"Error: length must be a positive integer\n");
              DisplayUsage();
            }
          }
          break;

        case 'e':
          while(ch != '=' && ch){
            ch = *(++ptr);
            l++;
          }

          if(!(*ptr) || l != 2){
            fprintf(stderr,"Error: incorrect flag format\n");
            DisplayUsage();
          }
          else{
            ptr++;
            char* endptr; 
            epsilon = strtod(ptr, &endptr); 
            if (epsilon < 0 || epsilon > 1){
              fprintf(stderr,"Error: range for epsilon is [0,1]\n");
              DisplayUsage();
            }
          }
          break;

        case 'd':
          while(ch != '=' && ch){
            ch = *(++ptr);
            l++;
          }

          if(!(*ptr) || l != 2){
            fprintf(stderr,"Error: incorrect flag format\n");
            DisplayUsage();
          }
          else{
            ptr++;
            char* endptr; 
            decay_rate = strtod(ptr, &endptr); 
            if (decay_rate < 0 || decay_rate > 1){
              fprintf(stderr,"Error: range for decay rate is [0,1]\n");
              DisplayUsage();
            }
          }
          break;

        case 'a':
          while(ch != '=' && ch){
            ch = *(++ptr);
            l++;
          }

          if(!(*ptr) || l != 2){
            fprintf(stderr,"Error: incorrect flag format\n");
            DisplayUsage();
          }
          else{
            ptr++;
            char* endptr; 
            learning_rate = strtod(ptr, &endptr); 
            if (learning_rate < 0 || learning_rate > 1){
              fprintf(stderr,"Error: range for learning rate is [0,1]\n");
              DisplayUsage();
            }
          }
          break;

        case 'g':
          while(ch != '=' && ch){
            ch = *(++ptr);
            l++;
          }

          if(!(*ptr) || l != 2){
            fprintf(stderr,"Error: incorrect flag format\n");
            DisplayUsage();
          }
          else{
            ptr++;
            char* endptr; 
            discount_rate = strtod(ptr, &endptr); 
            if (discount_rate < 0 || discount_rate > 1){
              fprintf(stderr,"Error: range for discount rate is [0,1]\n");
              DisplayUsage();
            }
          }
          break;


        case '-':
          ptr++;
          if(prefix("-logp",ptr)){
            if(dt){
              fprintf(stderr,"Error: targeting two descriptors is currently unsupported\n");
              DisplayUsage();
            }

            ch = *ptr;
            while(ch != '=' && ch)
              ch = *(++ptr); 

            ptr++;
            char* endptr; 
            logp = strtod(ptr, &endptr); 
            dt = 1;
          }
          else if(prefix("-molwt",ptr)){
            if(dt){
              fprintf(stderr,"Error: targeting two descriptors is currently unsupported\n");
              DisplayUsage();
            }

            ch = *ptr;
            while(ch != '=' && ch)
              ch = *(++ptr); 

            ptr++;
            char* endptr; 
            molwt = strtod(ptr, &endptr); 
            dt = 2;
          }
          else if (prefix("-print",ptr)){
            DisplayParameters();
          }
          else if (prefix("-runs",ptr)){
            while(ch != '=' && ch)
              ch = *(++ptr);

            if(!(*ptr)){
              fprintf(stderr,"Error: incorrect flag format\n");
              DisplayUsage();
            }
            else{
              ptr++;
              episodes = atoi(ptr);
              if(episodes < 0){
                fprintf(stderr,"Error: count must be a positive integer\n");
                DisplayUsage();
              }
            }
          }
          else if (prefix("-epsilon",ptr)){
            while(ch != '=' && ch)
              ch = *(++ptr);
            
            if(!(*ptr)){
              fprintf(stderr,"Error: incorrect flag format\n");
              DisplayUsage();
            }
            else{
              ptr++;
              char* endptr; 
              epsilon = strtod(ptr, &endptr); 
              if (epsilon < 0 || epsilon > 1){
                fprintf(stderr,"Error: range for epsilon is [0,1]\n");
                DisplayUsage();
              }
            }
          }
          else if (prefix("-decay",ptr)){
            while(ch != '=' && ch)
              ch = *(++ptr);
            
            if(!(*ptr)){
              fprintf(stderr,"Error: incorrect flag format\n");
              DisplayUsage();
            }
            else{
              ptr++;
              char* endptr; 
              decay_rate = strtod(ptr, &endptr); 
              if (decay_rate < 0 || decay_rate > 1){
                fprintf(stderr,"Error: range for decay rate is [0,1]\n");
                DisplayUsage();
              }
            }
          }
          else if (prefix("-alpha",ptr)){
            while(ch != '=' && ch)
              ch = *(++ptr);
            
            if(!(*ptr)){
              fprintf(stderr,"Error: incorrect flag format\n");
              DisplayUsage();
            }
            else{
              ptr++;
              char* endptr; 
              learning_rate = strtod(ptr, &endptr); 
              if (learning_rate < 0 || learning_rate > 1){
                fprintf(stderr,"Error: range for learning_rate is [0,1]\n");
                DisplayUsage();
              }
            }
          }
          else if (prefix("-gamma",ptr)){
            while(ch != '=' && ch)
              ch = *(++ptr);
            
            if(!(*ptr)){
              fprintf(stderr,"Error: incorrect flag format\n");
              DisplayUsage();
            }
            else{
              ptr++;
              char* endptr; 
              discount_rate = strtod(ptr, &endptr); 
              if (discount_rate < 0 || discount_rate > 1){
                fprintf(stderr,"Error: range for discount rate is [0,1]\n");
                DisplayUsage();
              }
            }
          }
          else if (prefix("-count",ptr)){
            while(ch != '=' && ch)
              ch = *(++ptr);

            if(!(*ptr)){
              fprintf(stderr,"Error: incorrect flag format\n");
              DisplayUsage();
            }
            else{
              ptr++;
              count = atoi(ptr);
              if(count < 0){
                fprintf(stderr,"Error: count must be a positive integer\n");
                DisplayUsage();
              }
            }
          }
          else if (prefix("-length",ptr)){
            while(ch != '=' && ch)
              ch = *(++ptr);

            if(!(*ptr)){
              fprintf(stderr,"Error: incorrect flag format\n");
              DisplayUsage();
            }
            else{
              ptr++;
              length = atoi(ptr);
              if(length < 0){
                fprintf(stderr,"Error: length must be a positive integer\n");
                DisplayUsage();
              }
            }
          }
          else{
            fprintf(stderr,"Error: incorrect input -%s\n",ptr);
            DisplayUsage();
          }
          break;

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
    else{
      switch(j++){
        default:
          train_files.push_back(ptr); 
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

  QGenerateWLN(wlnmodel);
  delete wlnmodel;
  return 0;
}
