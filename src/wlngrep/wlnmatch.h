/*##############################################################

Text match simulation for the automata of regular languaes, 
DFA, NFA and eNFA simultion. 

###############################################################*/



#ifndef REG_MATCH_H
#define REG_MATCH_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rfsm.h"

#define BUFF_SIZE 2048


struct MatchOptions{

};

void display_line(char *line){
  fprintf(stdout, "%s\n",line);
}

void display_highlighted_line(char *line, unsigned int spos, unsigned int epos){
  for(unsigned int i=0;i<BUFF_SIZE;i++){
    if(!line[i])
      break;
    if(i<spos)
      fprintf(stdout,"%c",line[i]);
    else if(i<epos){
      fprintf(stdout,"\033[1;31m");
      fprintf(stdout,"%c",line[i]);
      fprintf(stdout,"\033[0m");
    }
    else
      fprintf(stdout,"%c",line[i]);
  }

  fprintf(stdout,"\n");
}

void display_highlighted_match(char *line, unsigned int spos, unsigned int epos){
  
  for(unsigned int i=0;i<BUFF_SIZE;i++){
    if(!line[i] || line[i] == '\n')
      break;

    if(i>=spos && i<epos){
      fprintf(stdout,"\033[1;31m");
      fprintf(stdout,"%c",line[i]);
      fprintf(stdout,"\033[0m");
    }
  }
  fprintf(stdout,"\n");
}

void display_match(char *line, unsigned int spos, unsigned int epos){
  
  for(unsigned int i=0;i<BUFF_SIZE;i++){
    if(!line[i] || line[i] == '\n')
      break;

    if(i>=spos && i<epos)
      fprintf(stdout,"%c",line[i]);
  }
  fprintf(stdout,"\n");
}


/* matches the longest possible word using DFA
- 1 matches only, 2 - exact match only, 3- return all matches */
unsigned int DFAGreedyMatchLine(const char *inp, FSMAutomata *dfa, bool highlight, unsigned int opt_match_option=0, bool count=false){
  
  char line[BUFF_SIZE] = {0};
  strcpy(line,inp);
  unsigned int len = strlen(line);

  FSMState *state = dfa->root;

  int spos = -1; 
  int apos = -1;

  unsigned int n = 0;
  unsigned int match = 0; // 1 if true 
  unsigned char inp_char = *inp;

  while(n <= len){

    if((inp_char && inp_char != '\n') && state->access[inp_char]){
      state = state->access[inp_char];
      if(spos == -1)
        spos = n;

      if(state->accept)
        apos = n;
    }
    else{

      if(opt_match_option == 2){
        if(spos == 0 && (!inp_char || inp_char == '\n')){
          if(count)
            match++;
          else if(highlight)
            display_highlighted_line(line,0,n);
          else
            display_line(line);
        }
      }
      else if(apos >= 0 && spos >= 0 && spos <= apos){ // turn off single letter match here
        // failed in non_accepting state, what was the last accept state we saw
        if(count)
          match++;
        else if(highlight){
          if(opt_match_option==1)
            display_highlighted_match(line,spos,apos+1);
          else
            display_highlighted_line(line,spos,apos+1);
        }
        else{
          if(opt_match_option==1)
            display_match(line,spos,apos+1);
          else
            display_line(line);
        }

        apos = -1;
      }

      // resets the machine
      if(dfa->root->access[inp_char]){
        state = dfa->root->access[inp_char];
        spos = n;
        if(state->accept)
          apos = n;
      }
      else{
        state = dfa->root;
        spos = -1; // resets the start character
      }
    }

    n++;
    inp_char = *(++inp);      
  }
  
  return match;
}


  /* takes in input word and creates an epsilon bonded NFA
  - used for grep style word searching */
  bool ReadWord(const char *inp, FSMAutomata *nfa){
    if(!nfa->root)
      nfa->root = nfa->AddState(0);

    unsigned char ch = *inp;
    FSMState *prev = nfa->AddState();
    nfa->AddTransition(nfa->root,prev,0);
    FSMState *q = 0;
    while(ch){
      q = nfa->AddState();
      nfa->AddTransition(prev,q,ch);
      prev = q;
      ch = *(++inp);
    }

    nfa->MakeAccept(q);
    return true;
  }

#endif