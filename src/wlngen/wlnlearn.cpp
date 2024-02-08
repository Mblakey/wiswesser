/*
 first iteration of WLN learn, the idea is to use sparse rewards to inorder to 
 reduce the number of misses, this MAY (likely) lead to something like mode collapse 
 where the same inputs are given each time, this is good, at least we are trending to
 an learning answer. 
 */

#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

#define GEN_DEBUG 1
#define COUNT 5000
#define LOGP 2.5

int length = 5;
int episodes = 5;

double start_epsilon = 0.5;
double learning_rate = 0.5;
double decay_rate = 0.005;

bool opt_strings = false;

std::vector<const char*> train_files;


bool prefix(const char *pre, const char *str)
{
  return strncmp(pre, str, strlen(pre)) == 0;
}

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

bool seed_from_string(const char *str, FSMAutomata *wlnmodel){
  unsigned char ch = *str;
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;
  
  while(ch){
    for(edge=curr->transitions;edge;edge=edge->nxt){
      if(edge->ch == ch){
        if(edge->c < UINT32_MAX)
          edge->c+=100; // arbituary but seeds the starting string more than files 

        curr = edge->dwn; 
        break;
      }
    }
    ch = *(++str); 
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


FSMEdge *write_seed(const char *seed, FSMAutomata *wlnmodel, std::string &buffer ){
  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;
  unsigned char ch = *seed;
  while(ch){
    buffer += ch;
    for(edge=state->transitions;edge;edge=edge->nxt){
      if(edge->ch == ch){
        state = edge->dwn;
        break;
      }
    }
    ch = *(++seed);
  }
  return edge; 
}




double LogP(OBMol *mol){
  OBDescriptor* pDesc = OBDescriptor::FindType("logP");
  if(pDesc)
    return pDesc->Predict(mol); 
  else
    return 0.0;
}



/* reward function, if the WLN is valid back score */
bool Validate(const char *wln_str, OBMol *mol){
  if(!ReadWLN(wln_str,mol))
    return false;
  return true;
}

unsigned int ScoreFunction(const char*wln_str){
  unsigned int score = 0;
  OBMol mol;

  if(!Validate(wln_str,&mol))
    return score;
  else
    score +=1;
 
  double logp = LogP(&mol);
  if(logp < (LOGP + 0.5) && logp > (LOGP-0.5))
    score += 3;
  
  return score;
}

 
double DecayEpsilon(double epsilon_N0, double decay_rate, unsigned int iteration){
  double new_epsilon = epsilon_N0 * exp(-decay_rate*(iteration*10));
  if (new_epsilon > 0.1)
    return new_epsilon; 
  else
    return 0.1;
}

/* return an edge based on completely random distribution - highest exploration */
FSMEdge *RandomEdge(FSMState *curr, std::mt19937 &rgen){
  FSMEdge *e = 0; 
  unsigned int total = 0; 
  std::vector<FSMEdge*> edge_vec = {};
  std::vector<double> prob_vec = {}; // vector of probabilities 
  for(e=curr->transitions;e;e=e->nxt){
    edge_vec.push_back(e);
    total++; 
  }

  for (unsigned int i=0; i<total;i++)
    prob_vec.push_back(1/(double)total);
  
  std::discrete_distribution<> d(prob_vec.begin(), prob_vec.end());
  unsigned int chosen = d(rgen);
  return edge_vec[chosen];
}


/* return an edge based on the learnt probabilities - still high exploration */
FSMEdge *LikelyEdge(FSMState *curr, std::mt19937 &rgen){
  FSMEdge *e = 0; 
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


FSMEdge *ChooseEdge(FSMState *curr, double epsilon,std::mt19937 &rgen){
  std::uniform_real_distribution<> dis(0, 1);
  double choice = dis(rgen);
  if(choice > epsilon)
    return LikelyEdge(curr,rgen);
  else
    return RandomEdge(curr, rgen);
}


/* used to update the Q values in the chain if a successful hit is made */
void BellManEquation(FSMEdge *curr, unsigned int score){
 // curr->p += 0.02;
 // return;
  curr->p = ((1-learning_rate)*curr->p) + (learning_rate *score);
  return;
}


void NormaliseState(FSMState* state){
  FSMEdge *e = 0;
  double sum = 0.0; // should be a safe normalise value
  
  for (e=state->transitions;e;e=e->nxt)
    sum += e->p;   
  
//  double check = 0.0;
  for (e=state->transitions;e;e=e->nxt) {
    e->p = e->p/sum;
//    check += e->p;
//    fprintf(stderr,"%f ",e->p);
  }
  
//  fprintf(stderr,"\ncheck: %f\n",check);
  return;
}

void RewardPath(std::set<FSMEdge*> &path,unsigned int score){
  for(FSMEdge *e : path)
    BellManEquation(e, score);

  for(FSMEdge *e: path)
    NormaliseState(e->dwn);
  return;
}


/* uses Q learning to generate compounds from the language FSM
as a markov decision process, WLN is small enough that with 20
characters a large scope of chemical space can be covered. */
void QLearnWLN(  FSMAutomata *wlnmodel, double epsilon){
  int hits = 0; 
  int misses = 0;
  
  std::random_device rd;
  std::mt19937 rgen(rd());

  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  std::string wlnstr; 
  std::set<FSMEdge*> path; // avoid the duplicate path increases
  std::unordered_map<std::string, bool> unique;

  // if there is a target, start the string with that target
  
  if(opt_strings)
    edge = write_seed(train_files[0],wlnmodel,wlnstr);

  int strlength = 0;
  while(hits < COUNT){
    edge = ChooseEdge(state,epsilon,rgen); 
    
    if(edge->ch == '\n'){

      if(strlength >= length){ // only accepts will have new line, ensures proper molecule
        strlength = 0;
        state = wlnmodel->root;
        
        unsigned int score = 0; 
        score = ScoreFunction(wlnstr.c_str());   
        
        if(score){
          //fprintf(stderr,"%s\n", wlnstr.c_str());
          path.insert(edge);
          RewardPath(path,score);
          hits++;
        }
        else
         misses++;
        
        // move the machine to ensure that subunit is found
      
        wlnstr.clear();
        path.clear(); 
        if(opt_strings)
          edge = write_seed(train_files[0], wlnmodel, wlnstr);
          
      }
      else{
        // choose something else
        while(edge->ch == '\n')
          edge = ChooseEdge(state,epsilon,rgen);

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

  fprintf(stderr,"epsilon: %f, misses: %d\n", epsilon,misses);
}

void GenerateWLN(FSMAutomata *wlnmodel){
  int hits = 0; 

  std::random_device rd;
  std::mt19937 rgen(rd());

  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  std::string wlnstr; 
  std::unordered_map<std::string, bool> unique;

  if(opt_strings)
    edge = write_seed(train_files[0], wlnmodel, wlnstr);

  int strlength = 0;
  while(hits < COUNT){
    edge = LikelyEdge(state,rgen); 
    
    if(edge->ch == '\n'){

      if(strlength >= length){ // only accepts will have new line, ensures proper molecule
        strlength = 0;
        state = wlnmodel->root;
        
        OBMol mol; 
        if(Validate(wlnstr.c_str(), &mol)){
          fprintf(stderr,"%s\n", wlnstr.c_str());
          hits++;
        }
 
        wlnstr.clear();
        if(opt_strings)
          edge = write_seed(train_files[0], wlnmodel, wlnstr);
      }
      else{
          // choose something else
        while(edge->ch == '\n')
          edge = LikelyEdge(state,rgen);
      }
    }
    else{
      wlnstr += edge->ch;
      strlength++;
    }
   
    state = edge->dwn;
  }

}

/* compare old Q learning factors with new current, need 2 copies of the FSM*/
bool RunEpisodes(FSMAutomata *wlnmodel){
  // Get the initial Q-learning values for update loop. 
  // ~ 10 seconds per 50K molecules on ARM64 M1.   
  double epsilon = start_epsilon;
  for (unsigned int i=0;i<(unsigned int)episodes;i++){
    QLearnWLN(wlnmodel,epsilon);
    epsilon = DecayEpsilon(start_epsilon, decay_rate,i);
  }

  return true;
}
static void DisplayUsage(){
  fprintf(stderr, "wlngen <options> <trainfile>\n");
  fprintf(stderr,"options:\n");
  fprintf(stderr,"-l|--length=<int>      set length for generation        (default 5)\n");
  
  fprintf(stderr,"\ntuning:\n");
  fprintf(stderr,"-r|--runs=<int>        set learning episodes            (default 5)\n");
  fprintf(stderr,"-e|--epsilon=<double>  set epsilon hyperparameter       (default 0.5)\n");
  fprintf(stderr,"-d|--decay=<double>    set decay rate hyperparameter    (default 0.005)\n");
  fprintf(stderr,"-a|--alpha=<double>    set learning rate hyperparameter (default 0.5)\n");
  
  fprintf(stderr,"\ngeneral:\n");
  fprintf(stderr,"-p|--print             show all set hyperparameters and exit\n");
  fprintf(stderr,"-h|--help              show this help menu and exit\n");
  fprintf(stderr,"-s|--string            treat trainfile as input string\n");
  
  exit(1);
}


static void DisplayParameters(){
  fprintf(stderr,"----------------------------\n");
  fprintf(stderr,"target count:      %d\n",COUNT);
  fprintf(stderr,"target length:     %d\n",length);
  
  fprintf(stderr,"\nepisodes:          %d\n",episodes);
  fprintf(stderr,"learning rate:     %f\n",learning_rate);
  fprintf(stderr,"epsilon:           %f\n",start_epsilon);
  fprintf(stderr,"decay rate:        %f\n",decay_rate);
  fprintf(stderr,"strings:           %d\n", opt_strings); 
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
            start_epsilon = strtod(ptr, &endptr); 
            if (start_epsilon < 0 || start_epsilon > 1){
              fprintf(stderr,"Error: range for epsilon is [0,1]\n");
              DisplayUsage();
            }
          }
          break;

        case 's':
          opt_strings = true;
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

       case '-':
          ptr++;
          if (prefix("-print",ptr)){
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
              start_epsilon = strtod(ptr, &endptr); 
              if (start_epsilon < 0 || start_epsilon > 1){
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
      if(opt_strings){
        seed_from_string(trainfile, wlnmodel);
      }
      else{
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
  }

  RunEpisodes(wlnmodel);
  GenerateWLN(wlnmodel);
  delete wlnmodel;
  return 0;
}
