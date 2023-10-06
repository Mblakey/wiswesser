/*##############################################################

Merging operations for the regular automaton

###############################################################*/

#ifndef REG_MERGE
#define REG_MERGE

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <stack>
#include <vector>
#include <map>

#include "rfsm.h"
#include "rminimise.h"


FSMAutomata *Copy(FSMAutomata *fsm){
  
  FSMState *state = 0;
  FSMEdge *edge = 0;
  FSMAutomata *copy = new FSMAutomata(fsm->max_states,fsm->max_edges); 

  std::stack<FSMState*> stack; 
  std::map<FSMState*,bool> visited; 
  std::map<FSMState*,FSMState*> new_states; 

  stack.push(fsm->root);
  while(!stack.empty()){
    state = stack.top();
    visited[state] = true;
    stack.pop();

    FSMState *new_parent = 0;
    if(!new_states[state]){
      new_parent = copy->AddState(state->accept);
      new_states[state] = new_parent; 
    }
    else
      new_parent = new_states[state];
      
    for (edge=state->transitions;edge;edge=edge->nxt){
      FSMState *new_child = 0;
      if(edge->dwn){
        if(!new_states[edge->dwn]){
          new_child = copy->AddState(edge->dwn->accept);
          new_states[edge->dwn] = new_child; 
        }
        else
          new_child = new_states[edge->dwn];
          
        copy->AddTransition(new_parent,new_child,edge->ch);

        if(!visited[edge->dwn])
          stack.push(edge->dwn);
      }
    }

  }

  return copy; 
}

/* merges n FSM into an eNFA that is a union 
of all the input languages */
FSMAutomata* MergeUnion(std::vector<FSMAutomata*> fsms){
  FSMAutomata *enfa = new FSMAutomata(REASONABLE,REASONABLE); 

  FSMState *q = 0; 
  FSMEdge *e = 0;

  std::map<FSMState*,FSMState*> new_states; 
  FSMState *nfa_root = enfa->AddState(); 
  FSMState *nfa_accept = enfa->AddState(true); // init the root and final accept
  
  for (FSMAutomata *fsm : fsms){

    fsm->RemoveUnreachables();
    fsm->Reindex();
     
    for (unsigned int i=0;i<fsm->num_states;i++){
      FSMState *new_state = enfa->AddState(); 
      q = fsm->states[i];
      new_states[q] = new_state; 

      if (q == fsm->root)
        enfa->AddTransition(nfa_root,new_state,0);
      if(q->accept)
        enfa->AddTransition(new_state,nfa_accept,0);
    }

    // add all the edges in with the id offset 
    for (unsigned int i=0;i<fsm->num_states;i++){
      q = fsm->states[i];
      FSMState *nfa_src = new_states[q]; 
      
      for(e=q->transitions;e;e=e->nxt){
        FSMState *nfa_trg = new_states[e->dwn]; 
        enfa->AddTransition(nfa_src,nfa_trg,e->ch);
      }
    }

  }

  return enfa; 
}

/* merges all fsms into one new object, note
this will have multiple start states, and cannot be used in
matching - use for language equvilance*/
FSMAutomata* MergeParallel(std::vector<FSMAutomata*> fsms){
  FSMAutomata *parallel = new FSMAutomata(REASONABLE,REASONABLE); 

  FSMState *q = 0; 
  FSMEdge *e = 0;

  std::map<FSMState*,FSMState*> new_states;

  for (FSMAutomata *fsm : fsms){
    
    fsm->RemoveUnreachables();
    fsm->Reindex();
     
    for (unsigned int i=0;i<fsm->num_states;i++){
      FSMState *q = fsm->states[i];
      FSMState *new_state = parallel->AddState();
      new_states[q] = new_state;  
    } 
      
    // add all the edges in with the id offset 
    for (unsigned int i=0;i<fsm->num_states;i++){
      q = fsm->states[i];
      FSMState *parallel_src = new_states[q]; 
      if(q->accept)
        parallel->MakeAccept(parallel_src);
        
      for(e=q->transitions;e;e=e->nxt){
        FSMState *parallel_trg = new_states[e->dwn]; 
        parallel->AddTransition(parallel_src,parallel_trg,e->ch);
      }
    }
  }
	
  return parallel; 
}




#endif 