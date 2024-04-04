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
#define SINGLE_CHAR 0

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



// If the previous character was not locant 
void StackAmpersands(unsigned char ch,std::stack<unsigned char> &amp_stack){
  unsigned int closures = 0; 
  switch(ch){
    case 'Y':
    case 'B':
    case 'N':
      closures = 2; 
      break;

    case 'X':
    case 'K':
      closures = 3;
      break;

    case 'S':
      closures = 5;
      break;

    case 'P':
      closures = 4;
      break;
  }

  for(unsigned int i=0;i<closures;i++)
    amp_stack.push('&');

  return; 
}

bool isTerminator(unsigned char ch){
  switch(ch){
    case 'Q':
    case 'Z':
    case 'E':
    case 'F':
    case 'G':
    case 'I':
      return true;
  }
  
  return false; 
}

// If not on locant, pop stack
bool PopAmpersand(std::stack<unsigned char> &amp_stack){
  if(amp_stack.empty())
    return false;

  amp_stack.pop(); 
  return true; 
}

/* matches the longest possible word using DFA
- 1 matches only, 2 - exact match only, 3- return all matches */
unsigned int DFAGreedyMatchLine(const char *inp, FSMAutomata *dfa, bool highlight, bool invert,unsigned int opt_match_option=0, bool count=false){
  
  enum match_mode {WHOLE_LINE=0,MATCH_ONLY=1,EXACT=2};
  
  char line[BUFF_SIZE] = {0};
  strcpy(line,inp);
  unsigned int len = strlen(line);

  FSMState *state = dfa->root;

  int spos = -1; 
  int apos = -1;
  
  unsigned int l = 0;
  unsigned int n = 0;
  unsigned int match = 0; // 1 if true 
  unsigned char inp_char = *inp;
 
  bool reading_ring = false;
  bool expecting_locant = false;
  unsigned char locant = 0;

  std::stack<unsigned char> ampersand_stack; 
  ampersand_stack.push('&'); // all notation is allowed one closure. 

  while(n <= len){
    
    if(inp_char && state->access['*']) // accept all of barrie walkers changes
      inp_char = '*'; 

    if(inp_char && state->access[inp_char]){
      
      if(reading_ring){
        if(!expecting_locant && inp_char == 'J')
          reading_ring = false;
        // else do nothing
      }
      else{
        
        if(inp_char == ' '){
          expecting_locant = true; 
        }
        else if(expecting_locant){
          if(!locant){
            if(inp_char == '&'){
              // ion condition
              expecting_locant = false;
              while(!ampersand_stack.empty())
                ampersand_stack.pop(); 
            }
            else 
              locant = inp_char;
          }
          else if(inp_char != '&'){
            locant = 0;
            expecting_locant = false;  
          }
        }
        else if(inp_char == '&' && !PopAmpersand(ampersand_stack)){
          goto return_match;  
        }
      }

      state = state->access[inp_char];
      if(spos == -1)
        spos = n;

      if(state->accept)
        apos = n;
      
      l++;
    }
    else{
return_match:
      if(opt_match_option == EXACT){
        if(invert){
          if(n != len){
            if(count)
              match++;
            else
              display_line(line);
          }
        }
        else if(spos == 0 && !inp_char && state->accept && l > 1){
          if(count)
            match++;
          else if(highlight)
            display_highlighted_line(line,0,n);
          else
            display_line(line);
        }
        return match; 
      }
      else if(apos >= 0 && spos >= 0 && spos <= apos && l > 1){ // turn off single letter match here
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
          if(opt_match_option==MATCH_ONLY)
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

      l = 0;
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
