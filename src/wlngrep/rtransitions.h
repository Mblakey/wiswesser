/*##############################################################

Traditional Transition, extended transition functions

###############################################################*/


#ifndef REG_TRANSITIONS
#define REG_TRANSITIONS

#include <stdlib.h>
#include <stdio.h>

#include <set>
#include <utility>

#include "rfsm.h"

FSMState *SingletonTransition(FSMState *state,unsigned char ch){
  FSMEdge *edge = 0;
  for (edge=state->transitions;edge;edge=edge->nxt){
    if (edge->ch == ch)
      return edge->dwn;
  }
  return 0;
}

FSMState *ExtendedSingletonTransition(FSMState *state,const char *w){
  FSMState *curr = state; 
  unsigned char ch = *w; 
  while(ch){
    curr = SingletonTransition(curr,ch);
    if(!curr){
      fprintf(stderr,"Error: Extended transition from state is not possible\n");
      return 0;
    }
    ch = *(++w);
  }

  return curr;
}

bool SetTransition(std::set<FSMState*> states,unsigned char ch,std::set<FSMState*> &transition_states){
  FSMEdge *edge = 0;

  std::set<FSMState*>::iterator it; 
  for (it = states.begin(); it != states.end();it++){
    for (edge=(*it)->transitions;edge;edge=edge->nxt){
      if (edge->ch ==ch)
        states.insert(edge->dwn);
    }
  }
  
  return true;
}


#endif