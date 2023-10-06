/*##############################################################

Convert between the regular class of languages representations 
DFA, NFA, eNFA and regular expressions 

Subset algorithm is contained here along with epsilion transitions

###############################################################*/


#ifndef CONVERT_H
#define CONVERT_H

#include <stdio.h>
#include <stdlib.h>

#include <set>
#include <stack>
#include <queue> 
#include <map>

#include "rfsm.h"
#include "rminimise.h"

unsigned int max(unsigned int a, unsigned int b){
  return a > b ? a:b; 
}


std::set<FSMState*> TransitionFunction(std::set<FSMState*> S,unsigned char ch){
  std::set<FSMState*> valid_transitions;
  std::set<FSMState*>::iterator it;
  for (it = S.begin(); it != S.end(); ++it){
    FSMState *q = (*it); // starting state
    FSMEdge* e = 0;
    for (e = q->transitions;e;e=e->nxt){
      if(e->ch == ch)
        valid_transitions.insert(e->dwn);
    }
  }
  return valid_transitions;
}


/* perform epsilon closures, the edge sets have everything we need
- pass by reference allows easy set unions
returns bool based on whether eclosure states should all be accepts */
bool EpsilonClosure(std::set<FSMState*> S,std::set<FSMState*> &eclosure){

  bool accept_closure = false; 
  std::set<FSMState*>::iterator it;

  for (it = S.begin(); it != S.end(); ++it){

    FSMState *q = (*it); // starting state
    eclosure.insert(q); // e closure always includes itself 
    if(q->accept)
      accept_closure = true;

    std::stack<FSMState*> queue; 
    std::map<FSMState*,bool> visited; 
    queue.push(q);

    FSMState *top = 0;
    while(!queue.empty()){
      top = queue.top();
      queue.pop();

      visited[top] = true;

      FSMEdge* e = 0;
      for (e = top->transitions;e;e=e->nxt){
        if(!e->ch && !visited[e->dwn]){
          eclosure.insert(e->dwn);
          if(e->dwn->accept)
            accept_closure = true; 
  
          queue.push(e->dwn);
        } 
      }
    }
  }

  return accept_closure;
}


FSMAutomata *eNFAtoNFA(FSMAutomata *machine){
  FSMAutomata *nfa = new FSMAutomata(machine->max_states,machine->max_edges);

  FSMState *curr = 0; 
  FSMState* new_state = 0;
  bool accept_set = false;

  std::map<FSMState*,FSMState*> new_states;
  std::set<FSMState*> first_eclosure;   
  std::set<FSMState*> second_eclosure;  
  std::set<FSMState*> vtransitions; 
  
  std::set<FSMState*>::iterator it; 
  std::map<FSMState*,bool> seen;

  std::queue<FSMState*> reachable; 
  reachable.push(machine->root);

  while(!reachable.empty()){

    curr = reachable.front();
    reachable.pop(); // already evaluated;
    seen[curr] = true;

    accept_set = EpsilonClosure({curr},first_eclosure); // eclosure regardless of character

    if(!new_states[curr]){
      new_state = nfa->AddState();
      new_states[curr] = new_state; 
    }
    else
      new_state = new_states[curr];
    
    if(accept_set)
      nfa->MakeAccept(new_state);

    // closure per character in potential alphabet. 
    for (unsigned int ch=1;ch<255;ch++){
      if (!machine->alphabet[ch])
        continue;

      vtransitions = TransitionFunction(first_eclosure,ch); 
      if(vtransitions.empty())
        continue;

      accept_set = EpsilonClosure(vtransitions,second_eclosure);

      // update reachable states from the set closure
      for (it = second_eclosure.begin();it != second_eclosure.end();++it){
        FSMState* new_child = 0;
        if(!seen[(*it)]){
          reachable.push(*it);
          seen[(*it)] = true;
          new_child = nfa->AddState();
          new_states[(*it)] = new_child; 
        }
        else
          new_child = new_states[(*it)];

        if(accept_set)
          nfa->MakeAccept(new_child);

        nfa->AddTransition(new_state,new_child,ch);
      }

      second_eclosure.clear();
    }

    first_eclosure.clear();
  }
  
  nfa->RemoveUnreachables();
  return nfa;
}

/* uses subset contruction to turn an NFA into DFA */
FSMAutomata *NFAtoDFA(FSMAutomata *machine){
  FSMAutomata *dfa = new FSMAutomata(machine->max_states,machine->max_edges);

  FSMState* new_state = 0;
  FSMState* trans_state = 0;
  std::set<FSMState*> curr; 
  std::set<FSMState*> vtransitions; 
  std::map< std::set<FSMState*>,bool> seen_sets; 
  std::map<std::set<FSMState*>,FSMState*> set_to_singleton; 
  
  std::set<FSMState*>::iterator it;
  
  std::queue<std::set<FSMState*>> reachable; 
  reachable.push({machine->root});

  while(!reachable.empty()){

    curr = reachable.front();
    reachable.pop(); // already evaluated;

    if(!seen_sets[curr]){
      new_state = dfa->AddState();
      set_to_singleton[curr] = new_state;
      seen_sets[curr] = true;
    }
    else
      new_state = set_to_singleton[curr];
      
    for (unsigned int ch=1;ch<255;ch++){
      if (!machine->alphabet[ch])
        continue;
    
      vtransitions = TransitionFunction(curr,ch);
      if(vtransitions.empty())
        continue;
      
      if(!seen_sets[vtransitions]){
        reachable.push(vtransitions);
        trans_state = dfa->AddState();
        set_to_singleton[vtransitions] = trans_state;

        for (it = vtransitions.begin();it != vtransitions.end();++it){
          if((*it)->accept){
            dfa->MakeAccept(trans_state);
            break;
          }
        }

        seen_sets[vtransitions] = true;
      }
      else
        trans_state = set_to_singleton[vtransitions];

      dfa->AddTransition(new_state,trans_state,ch);
    }
  
  }

  dfa->RemoveUnreachables();
  return dfa;
}


/* a direct eNFA to DFA with Epsilon closures, saves computation on sets
compared to eNFA->NFA->DFA */
FSMAutomata *eNFAtoDFA(FSMAutomata *machine){
  FSMAutomata *dfa = new FSMAutomata(machine->max_states,machine->max_edges);

  bool accept_set = false;
  FSMState* new_state = 0;
  FSMState* trans_state = 0;
  std::set<FSMState*> curr; 
  std::set<FSMState*> vtransitions; 
  std::set<FSMState*> first_eclosure;  
  std::set<FSMState*> second_eclosure;  
  std::map< std::set<FSMState*>,bool> seen_sets; 
  std::map< std::set<FSMState*>,FSMState*> set_to_singleton; 

  std::set<FSMState*>::iterator it;
  
  std::queue<std::set<FSMState*>> reachable; 
  reachable.push({machine->root});

  while(!reachable.empty()){
    curr = reachable.front();
    reachable.pop(); // already evaluated;

    accept_set = EpsilonClosure({curr},first_eclosure); // eclosure regardless of character

    if(!seen_sets[first_eclosure]){
      new_state = dfa->AddState(accept_set);

      set_to_singleton[first_eclosure] = new_state;
      seen_sets[first_eclosure] = true;
    }
    else
      new_state = set_to_singleton[first_eclosure];


    for (unsigned int ch=1;ch<255;ch++){
      if (!machine->alphabet[ch])
        continue;
    
      vtransitions = TransitionFunction(first_eclosure,ch);
      accept_set = EpsilonClosure(vtransitions,second_eclosure); 
      if(second_eclosure.empty())
        continue;

      if(!seen_sets[second_eclosure]){
        reachable.push(second_eclosure);
        trans_state = dfa->AddState(accept_set);
        set_to_singleton[second_eclosure] = trans_state;
        seen_sets[second_eclosure] = true;
      }
      else
        trans_state = set_to_singleton[second_eclosure];
      
      dfa->AddTransition(new_state,trans_state,ch);
      second_eclosure.clear();
    }

    first_eclosure.clear(); 
  }

  dfa->RemoveUnreachables();
  return dfa; 
}



/* Creates DFA from eNFA or NFA, memory is not handled, remember to delete machines */
FSMAutomata *ConvertToDFA(FSMAutomata *machine){
  
  FSMAutomata *dfa = 0;
  machine->Categorize();

  if(machine->type == NFA)
    dfa = NFAtoDFA(machine);
  else if (machine->type == eNFA)
    dfa = eNFAtoDFA(machine);
  else{
    fprintf(stderr,"Warning: FSM already DFA returning ptr\n");
    return machine;
  }
  
  dfa->Categorize();
  if(dfa->type != DFA){
    fprintf(stderr,"Error: conversion to DFA FAILED\n");
    return 0;
  }

  return dfa; 
}



#endif