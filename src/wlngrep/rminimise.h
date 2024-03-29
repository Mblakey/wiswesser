/*##############################################################

Minimise and prove language equivalence for DFA

###############################################################*/

#ifndef REG_MINIMISE
#define REG_MINIMISE

#include <stdlib.h>
#include <stdio.h>

#include <set>
#include <map>
#include <stack>
#include <vector>

#include "rfsm.h"
#include "rtransitions.h"


struct StateRank{
  FSMState* state;
  unsigned int p; // which partition its currently in.
};


void SortOnID(StateRank **arr,unsigned int len){
  for (unsigned int j=1;j< len;j++){
		unsigned int key = arr[j]->state->id;
    StateRank *ptr = arr[j];
		int i = j-1;
		while(i >= 0 && arr[i]->state->id > key){
			arr[i+1] = arr[i];
			i--;
		}
		arr[i+1] = ptr;
	}
}

void SortOnPartition(StateRank **arr,unsigned int len){
	for (unsigned int j=1;j< len;j++){
		unsigned int key = arr[j]->p;
    StateRank *ptr = arr[j];
		int i = j-1;
		while(i >= 0 && arr[i]->p > key){
			arr[i+1] = arr[i];
			i--;
		}
		arr[i+1] = ptr;
	}
}

/* least significant is the state ID num, then partition is major */
void OrderPartition(StateRank **arr,unsigned int len){
  SortOnID(arr,len);
  SortOnPartition(arr,len);
}
  
  
void PrintPartitions(StateRank **arr,unsigned int len, bool partition_nums){
  int last_seen = -1;
  for (unsigned int i=0;i<len;i++){

    if(!arr[i])
      continue;
    int partition = arr[i]->p;
    if(last_seen != partition){

      if(last_seen != -1)
        fprintf(stderr,"},");
      
      fprintf(stderr,"{ ");
      last_seen = partition;
    }
    fprintf(stderr,"%d ",arr[i]->state->id);
  }
  fprintf(stderr,"}\n");

  if(partition_nums){
    int last_seen = -1;
    for (unsigned int i=0;i<len;i++){

      if(!arr[i])
        continue;
      int partition = arr[i]->p;
      if(last_seen != partition){

        if(last_seen != -1)
          fprintf(stderr,"},");
        
        fprintf(stderr,"{ ");
        last_seen = partition;
      }
      fprintf(stderr,"%d ",arr[i]->p);
    }
    fprintf(stderr,"}\n");
  }
}


/* Two states p and q are distinguishable in partition Pk for any input symbol ‘a’,
if δ (p, a) and δ (q, a) are in different sets in partition Pk-1. 
- make all comparsions, then move marked states? P-1 is the old partition before any movment */
bool Distinguishable(FSMState *p,FSMState *q,StateRank **FSMPartition, FSMAutomata *dfa){
  for (unsigned char ch = 0; ch<255;ch++){
    if(dfa->alphabet[ch]){
      FSMState *fp = SingletonTransition(p,ch);
      FSMState *fq = SingletonTransition(q,ch);

      if(fp==fq)
        continue;

      if(fp && fq){ // both must transition, if one doesnt there could be an optimisation here

        unsigned int fp_partition = 0;
        unsigned int fq_partition = 0;
        
        for(unsigned int i=0;i<dfa->num_states;i++){
          if(FSMPartition[i]->state == fp)
            fp_partition = FSMPartition[i]->p;
            
          if(FSMPartition[i]->state == fq)
            fq_partition = FSMPartition[i]->p;
        }

        if(fp_partition != fq_partition){
          //fprintf(stderr,"  %d and %d are distinguishable: fp(%d) is %d, fq(%d) is %d\n",p->id,q->id,fp->id,fp_partition,fq->id,fq_partition);
          return true;
        }
      }
      else if(fp || fq){
        //fprintf(stderr,"  %d and %d are distinguishable: different outputs\n",p->id,q->id);
        return true; 
      }
    }
  }
  return false;
}

void CopyPartition(StateRank **src, StateRank **trg, unsigned int len){
  for(unsigned int i=0;i < len;i++){
    trg[i]->state = src[i]->state;
    trg[i]->p = src[i]->p;
  }
}

/* Partition Refinement algorithm using equivalance theorum */
void PartitionRefinement(StateRank **FSMPartition, FSMAutomata *dfa){

  unsigned int current_partitions = 2; 
  StateRank **FSMPartitionPREV = (StateRank**)malloc(sizeof(StateRank*) * dfa->num_states);  // P-1 therefore a copy needed 
  for(unsigned int i=0;i < dfa->num_states;i++)
    FSMPartitionPREV[i] = (StateRank*)malloc(sizeof(StateRank));
    
  CopyPartition(FSMPartition,FSMPartitionPREV,dfa->num_states); // inits the pointers

  bool WorkDone = true;
  while(WorkDone){

    WorkDone = false;
    for(unsigned int i=0;i<dfa->num_states;i++){
      StateRank *s              =   FSMPartition[i];
      FSMState  *x              =   s->state;
      unsigned int j            =   i+1; 
      unsigned int x_partition  =   s->p; // starting partition

      // we want to test every forward state in the partition
      while(j<dfa->num_states && FSMPartition[j]->p == x_partition){
        StateRank *r = FSMPartition[j];
        FSMState  *y =  r->state;

        if(Distinguishable(x,y,FSMPartitionPREV,dfa)){
          r->p = current_partitions; // if distinguishable, we move it to a new set 
          WorkDone = true;
        } 
        j++;
      }

      if(WorkDone){
        current_partitions++; // new partition must of been made
        OrderPartition(FSMPartition,dfa->num_states);
        CopyPartition(FSMPartition,FSMPartitionPREV,dfa->num_states);
        break;
      }
    }
  }

  // release all memory
  for(unsigned int i=0;i < dfa->num_states;i++)
    free(FSMPartitionPREV[i]);
  
  free(FSMPartitionPREV);
  return;
}

/* create the minised DFA */
FSMAutomata* CreateMinimalDFA(StateRank **FSMPartition, FSMAutomata *dfa){

  FSMAutomata *minimal = new FSMAutomata(dfa->max_states,dfa->max_edges); // will sort out sizing etc 

  FSMState *min_state = 0;
  std::map<FSMState*,FSMState*> new_states; 
  std::map<unsigned int,bool> partition_created;

  // init the new states 
  for(unsigned int i=0;i<dfa->num_states;i++){
    StateRank *s = FSMPartition[i];
    if(!partition_created[s->p]){
      min_state = minimal->AddState(s->state->accept);
      partition_created[s->p] = true;
    }
      
    if(s->state == dfa->root)
      minimal->root = min_state;
      
    new_states[s->state] = min_state;
  }

  FSMEdge *edge = 0;
  for(unsigned int i=0;i<dfa->num_states;i++){
    StateRank *s = FSMPartition[i];
    FSMState* src = s->state;
    for(edge=src->transitions;edge;edge=edge->nxt){
      FSMState *trg = edge->dwn;
      minimal->AddTransition(new_states[src],new_states[trg],edge->ch);
    }
  }

  minimal->RemoveUnreachables();
  return minimal;
}

/* use hopcrofts partitions with pre-init lists to minimise DFA, uses Radix style partition sort */
FSMAutomata *MinimiseDFA(FSMAutomata *dfa){
  if(dfa->type != DFA){
    fprintf(stderr,"Error: calling minimise on a FSM other than a DFA is undefined\n");
    return 0;
  }

  if(!dfa->num_accepts){
    fprintf(stderr,"Error: minimising DFA without any accept states is undefined");
    return 0;
  }

  dfa->RemoveUnreachables(); // also reindexes


  FSMAutomata *optimal = 0;
  StateRank **FSMPartition = (StateRank**)malloc(sizeof(StateRank*) * dfa->num_states); // ensures the set properties, pigeon hole principle
  for(unsigned int i=0;i < dfa->num_states;i++){
    FSMPartition[i] = (StateRank*)malloc(sizeof(StateRank));
    FSMPartition[i]->state = dfa->states[i];
    FSMPartition[i]->p = 0; // ahh okay here as well
  }

  // place all accepts in one partition, then refine
  for(unsigned int i=0;i < dfa->num_states;i++){
    if(FSMPartition[i]->state->accept)
      FSMPartition[i]->p = 1;
  }

  OrderPartition(FSMPartition,dfa->num_states); // group the partitions together
  PartitionRefinement(FSMPartition,dfa);
  optimal = CreateMinimalDFA(FSMPartition,dfa);

  // release all memory
  for(unsigned int i=0;i < dfa->num_states;i++)
    free(FSMPartition[i]);
  
  free(FSMPartition);

  optimal->Categorize();
  if(optimal->type != DFA){
    fprintf(stderr,"Error in DFA minimisation\n");
    delete optimal;
    return 0;
  }
  else
    return optimal;
}



#endif 