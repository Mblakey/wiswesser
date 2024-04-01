#ifndef READ_DOT_H
#define READ_DOT_H

#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <stdio.h>

#include "rfsm.h"
#include "read_file.h" 


bool ParseStates(const char *inp, FSMAutomata *fsm){
  FILE *fp = fopen(inp,"r"); 
  if (!fp){
    fprintf(stderr, "Error: could not open file\n"); 
    return 0;
  }
  char buffer [1024] = {0}; 
  while(ReadLineFromFile(fp, buffer, 1024,0)){
    int first_state = 0; 
    char label_type[64] = {0}; 
    if (sscanf(buffer, "%d [shape= %s,label=\"\"]", &first_state, label_type) == 2){
      if(!strncmp(label_type, "doublecircle",12)){
        // make an accept
        FSMState *s = fsm->AddState(true); 
        s->id = first_state; 
      }
      else{
        // make a normal
        FSMState *s = fsm->AddState(false); 
        s->id = first_state; 
      }
    }
  }
  fclose(fp);
  return true; 
}

bool ParseEdges(const char *inp, FSMAutomata *fsm){
  FILE *fp = fopen(inp,"r"); 
  if (!fp){
    fprintf(stderr, "Error: could not open file\n"); 
    return 0;
  }


  char buffer [1024] = {0}; 
  while(ReadLineFromFile(fp, buffer, 1024,0)){
    int first_state = 0; 
    int second_state = 0;
    char transition = 0; 
    if (sscanf(buffer, "%d -> %d [label=\"%c\"]", &first_state, &second_state, &transition) == 3) {
        FSMState *f = 0; 
        FSMState *l = 0; 
        for(unsigned int i=0;i<fsm->num_states;i++){
          if(fsm->states[i]->id == first_state)
            f = fsm->states[i]; 

          if(fsm->states[i]->id == second_state)
            l = fsm->states[i]; 
        }
        
        fsm->AddTransition(f, l, transition);
        fprintf(stderr,"assigning: %d - %c-> %d\n",first_state,transition,second_state); 
    }
  }
  fclose(fp);
  return true; 
}

FSMAutomata *FSMFromDotFile(const char *inp){
  FSMAutomata *fsm =  new FSMAutomata(8512, 50000); 
  ParseStates(inp,fsm); 
  ParseEdges(inp, fsm); 
  return fsm; 
}


#endif 

