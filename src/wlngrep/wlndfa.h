#ifndef WLN_DFA_H
#define WLN_DFA_H

#include "rfsm.h"
#include "rconvert.h"

/* will purge and delete to merge into the main fsm */
FSMState * InsertAcyclic(FSMAutomata *acyclic){
  FSMState *root = acyclic->AddState(false);

  FSMState *func_group = acyclic->AddState(true); 
  FSMState *digits = acyclic->AddState(true);
  
  FSMState *branch = acyclic->AddState(true); // '&' characters
  FSMState *double_bond = acyclic->AddState(false); // 'U'
  FSMState *triple_bond = acyclic->AddState(false);  // 'UU'

  FSMState *open_dash = acyclic->AddState(false); // 
  FSMState *close_dash = acyclic->AddState(true);

  FSMState *element_a = acyclic->AddState(false);
  FSMState *element_b = acyclic->AddState(false);
  FSMState *hypervalent = acyclic->AddState(false);

  // set up digits state, cant be zero on entrance
  for(unsigned char ch = '0';ch <= '9';ch++)
    acyclic->AddTransition(digits,digits,ch);

  // set up special elements 
  acyclic->AddTransition(root,open_dash,'-');
  acyclic->AddTransition(func_group,open_dash,'-');
  acyclic->AddTransition(digits,open_dash,'-');

  // allow any combo for element definition, at the moment
  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    acyclic->AddTransition(open_dash,element_a,ch);
    acyclic->AddTransition(element_a,element_b,ch);
  }

  // add in single char hypervalence 
  acyclic->AddTransition(open_dash,hypervalent,'P');
  acyclic->AddTransition(open_dash,hypervalent,'S');
  acyclic->AddTransition(open_dash,hypervalent,'E');
  acyclic->AddTransition(open_dash,hypervalent,'F');
  acyclic->AddTransition(open_dash,hypervalent,'G');
  acyclic->AddTransition(open_dash,hypervalent,'I');
  acyclic->AddTransition(open_dash,hypervalent,'E');

  // close both dash blocks
  acyclic->AddTransition(hypervalent,close_dash,'-');
  acyclic->AddTransition(element_b,close_dash,'-');

  // could be element followed by element
  acyclic->AddTransition(close_dash,open_dash,'-');


  // link double and triple bonds
  acyclic->AddTransition(double_bond,triple_bond,'U');


  // same for open dash
  acyclic->AddTransition(double_bond,open_dash,'-');
  acyclic->AddTransition(triple_bond,open_dash,'-');

  // allow the double start to come from symbols
  acyclic->AddTransition(func_group,double_bond,'U');
  acyclic->AddTransition(digits,double_bond,'U');
  acyclic->AddTransition(close_dash,double_bond,'U');

  
  // add unlimited branching to and from any symbol 
  acyclic->AddTransition(func_group,branch,'&');
  acyclic->AddTransition(digits,branch,'&');
  acyclic->AddTransition(close_dash,branch,'&');
  
  // allow repeats
  acyclic->AddTransition(branch,branch,'&');
  acyclic->AddTransition(branch,open_dash,'&');

  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    switch(ch){
      case 'L':
      case 'T':
      case 'D':
      case 'J':
      case 'A':
      case ' ':
      case '-':
      case '&':
      case '/':
      case 'U':
      //case 'H': 
      case 'R': // seperate out benzene here
        break;
      
      default:
        acyclic->AddTransition(root,func_group,ch);
        acyclic->AddTransition(func_group,func_group,ch);
        acyclic->AddTransition(digits,func_group,ch);

        acyclic->AddTransition(close_dash,func_group,ch);

        acyclic->AddTransition(double_bond,func_group,ch);
        acyclic->AddTransition(triple_bond,func_group,ch);

        acyclic->AddTransition(branch,func_group,ch);
        break;
    }
  }

  // same for digits
  for(unsigned char ch = '1';ch <= '9';ch++){
    acyclic->AddTransition(root,digits,ch);

    acyclic->AddTransition(func_group,digits,ch);

    acyclic->AddTransition(close_dash,digits,ch);

    acyclic->AddTransition(double_bond,digits,ch);
    acyclic->AddTransition(triple_bond,digits,ch);

    acyclic->AddTransition(branch,digits,ch);
  }

  acyclic->AddTransition(branch,double_bond,'U');
   
  return root;
}


/* just handles the internal ring structure */
FSMState *InsertCyclic(FSMAutomata *cyclic){

  FSMState *root = cyclic->AddState();

  /* ######### CYCLIC DEFINITIONS ######### */
  
  FSMState *open_ring = cyclic->AddState(false);
  FSMState *close_ring = cyclic->AddState(true);

  cyclic->AddTransition(root,open_ring,'L');
  cyclic->AddTransition(root,open_ring,'T');
  cyclic->AddTransition(root,open_ring,'D');


  // handle all the ring digits 

  FSMState *ring_digits = cyclic->AddState(false); // not an accept now
  for(unsigned char ch = '0';ch <= '9';ch++)
    cyclic->AddTransition(ring_digits,ring_digits,ch); // zero here allows the matching if metallocens


  for(unsigned char ch = '1';ch <= '9';ch++)
    cyclic->AddTransition(open_ring,ring_digits,ch);

  FSMState *big_ring_dash_open = cyclic->AddState(false);
  FSMState *big_ring_dash_close = cyclic->AddState(false);
  FSMState *big_ring_digits = cyclic->AddState(false);

  cyclic->AddTransition(open_ring,big_ring_dash_open,'-'); // L-


  for(unsigned char ch = '1';ch <= '9';ch++)
    cyclic->AddTransition(big_ring_dash_open,big_ring_digits,ch); // L-6
  
  for(unsigned char ch = '0';ch <= '9';ch++)
    cyclic->AddTransition(big_ring_digits,big_ring_digits,ch); // L-666... 

  cyclic->AddTransition(big_ring_digits,big_ring_dash_close,'-'); // L-666-
  
  for(unsigned char ch = '1';ch <= '9';ch++)
    cyclic->AddTransition(big_ring_dash_close,ring_digits,ch); // L-666-6
  
  cyclic->AddTransition(ring_digits,big_ring_dash_open,'-'); // L6-
  cyclic->AddTransition(big_ring_dash_close,big_ring_dash_open,'-'); // L-6--


  // POLY CYCLICS RING NODES

  FSMState *digit_space = cyclic->AddState(false);
  FSMState *digit_locant = cyclic->AddState(false);

  cyclic->AddTransition(digit_locant,digit_locant,'&'); //  L E&6
  cyclic->AddTransition(digit_locant,digit_locant,'-'); //  L E&6

  cyclic->AddTransition(digit_locant,digit_space,' '); //  forward bridge notation

  cyclic->AddTransition(open_ring,digit_space,' '); // L' '

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    cyclic->AddTransition(digit_space,digit_locant,ch); // 'L' 'A

  for(unsigned char ch = '1';ch <= '9';ch++)
    cyclic->AddTransition(digit_locant,ring_digits,ch); // L' 'A6..

  cyclic->AddTransition(digit_locant,big_ring_dash_open,'-'); // L' 'A-6-..

  cyclic->AddTransition(ring_digits,digit_space,' '); // L666' 'A
  cyclic->AddTransition(big_ring_dash_close,digit_space,' '); // L66-6-' 'A

  // MULTI CYCLICS NOTATION

  FSMState *multi_space = cyclic->AddState(false);
  FSMState *multi_digit = cyclic->AddState(false);
  FSMState *multi_locants = cyclic->AddState(false);

  cyclic->AddTransition(digit_locant,multi_space,' '); //  forward bridge notation
  
  cyclic->AddTransition(ring_digits,multi_space,' '); // L6' '
  cyclic->AddTransition(big_ring_dash_close,multi_space,' '); // L-6-' '
  
  for(unsigned char ch = '1';ch <= '9';ch++)
    cyclic->AddTransition(multi_space,multi_digit,ch); //  L6' '2

  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    cyclic->AddTransition(multi_digit,multi_locants,ch); //  L6' '2A
    cyclic->AddTransition(multi_locants,multi_locants,ch); //  L6' '2AAAA...
  }

  FSMState *multi_size_space = cyclic->AddState(false);
  FSMState *multi_size = cyclic->AddState(false);

  cyclic->AddTransition(multi_locants,multi_locants,'&'); //  L E&6
  cyclic->AddTransition(multi_locants,multi_locants,'-'); //  L E-6 // broken locant 

  cyclic->AddTransition(multi_locants,multi_size_space,' '); //  L6' '2AAAA' '
  

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    cyclic->AddTransition(multi_size_space,multi_size,ch); //  L6' '2AAAA...P    

  // give size expansion ability
  cyclic->AddTransition(multi_size,multi_size,'&'); //  L E&6

  // Branching rings notation 
  FSMState *pair_slash = cyclic->AddState(false);
  FSMState *pair_loc_a = cyclic->AddState(false); 
  FSMState *pair_loc_b = cyclic->AddState(false);  
  
  cyclic->AddTransition(ring_digits,pair_slash,'/'); //  L6/
  cyclic->AddTransition(big_ring_dash_close,pair_slash,'/'); //L-6-/

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    cyclic->AddTransition(pair_slash,pair_loc_a,ch); //  L6' '2AAAA...P    

  cyclic->AddTransition(pair_loc_a,pair_loc_a,'&'); //  L E&6
  cyclic->AddTransition(pair_loc_a,pair_loc_a,'-'); //  L E&6

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    cyclic->AddTransition(pair_loc_a,pair_loc_b,ch); //  L6' '2AAAA...P

  cyclic->AddTransition(pair_loc_b,pair_loc_b,'&'); //  L E&6
  cyclic->AddTransition(pair_loc_b,pair_loc_b,'-'); //  L E&6

  cyclic->AddTransition(pair_loc_b,pair_slash,'/');
  cyclic->AddTransition(pair_loc_b,digit_space,' ');
  cyclic->AddTransition(pair_loc_b,multi_space,' ');

  // RING HETERO ATOMS

  // handle all the ring heteroatomic assignments, (must come from digits)  
  FSMState *hetero_space =  cyclic->AddState(false);
  FSMState *hetero_locant = cyclic->AddState(false);
  FSMState *hetero_atom =   cyclic->AddState(false);

  FSMState *hetero_open_dash =   cyclic->AddState(false);
  FSMState *hetero_close_dash =   cyclic->AddState(false);
  FSMState *hetero_element_a =   cyclic->AddState(false);
  FSMState *hetero_element_b =   cyclic->AddState(false);
  FSMState *hetero_hypervalent =   cyclic->AddState(false);


  cyclic->AddTransition(ring_digits,hetero_open_dash,'-'); // L6P
  cyclic->AddTransition(big_ring_dash_close,hetero_open_dash,'-'); // L-6-P

  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    cyclic->AddTransition(hetero_open_dash,hetero_element_a,ch);
    cyclic->AddTransition(hetero_element_a,hetero_element_b,ch);
  }


  cyclic->AddTransition(hetero_hypervalent,hetero_close_dash,'-');
  cyclic->AddTransition(hetero_element_b,hetero_close_dash,'-');

  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    switch(ch){
      case 'L':
      case 'T':
      case 'D':
      case 'J':
      case 'A':
      case ' ':
      case '-':
      //case '&':
      case '/':
      // case 'U':
      // case 'H': // this are allowed here!
        break;
      
      default:
        cyclic->AddTransition(hetero_locant,hetero_atom,ch);
        cyclic->AddTransition(hetero_atom,hetero_atom,ch); // L6PPP...
        cyclic->AddTransition(ring_digits,hetero_atom,ch); // L6P
        cyclic->AddTransition(big_ring_dash_close,hetero_atom,ch); // L-6-P
        cyclic->AddTransition(hetero_close_dash,hetero_atom,ch);
        break; // direct connections
    }
  }

  cyclic->AddTransition(hetero_locant,hetero_open_dash,'-'); 

  cyclic->AddTransition(hetero_close_dash,hetero_open_dash,'-'); 
  cyclic->AddTransition(hetero_close_dash,hetero_space,' ');
 
  
  cyclic->AddTransition(hetero_atom,hetero_space,' '); // L6P A...
  cyclic->AddTransition(hetero_atom,hetero_open_dash,'-'); // L6P A...

  cyclic->AddTransition(ring_digits,hetero_space,' '); // L6 AP...
  cyclic->AddTransition(big_ring_dash_close,hetero_space,' '); // L-6-' '.

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    cyclic->AddTransition(hetero_space,hetero_locant,ch); // L6' 'A..

  cyclic->AddTransition(hetero_locant,hetero_locant,'&'); // L6P A&...
  

  // multi atom attachment 
  cyclic->AddTransition(multi_size,hetero_space,' '); //  L6' '2AAAA' 'PJ

  // allow briding and locant expansion
  cyclic->AddTransition(hetero_locant,hetero_space,' '); 


  // AROMATICS
  FSMState *aromatics = cyclic->AddState(false);
  cyclic->AddTransition(aromatics,aromatics,'&');
  cyclic->AddTransition(aromatics,aromatics,'T'); // L ... &T&T&TJ

  cyclic->AddTransition(big_ring_dash_close,aromatics,'&');
  cyclic->AddTransition(big_ring_dash_close,aromatics,'T');

  cyclic->AddTransition(ring_digits,aromatics,'&');
  cyclic->AddTransition(ring_digits,aromatics,'T');

  cyclic->AddTransition(hetero_close_dash,aromatics,'&');
  cyclic->AddTransition(hetero_close_dash,aromatics,'T');

  cyclic->AddTransition(hetero_atom,aromatics,'&');
  cyclic->AddTransition(hetero_atom,aromatics,'T');

  // multi atom attachment 
  cyclic->AddTransition(multi_size,aromatics,' '); //  L6' '2AAAA' 'PJ
  cyclic->AddTransition(multi_size,aromatics,'T');
  cyclic->AddTransition(multi_size,aromatics,'&');


  // ALL CLOSURES
  cyclic->AddTransition(ring_digits,close_ring,'J');
  cyclic->AddTransition(multi_size,close_ring,'J'); //  L6' '2AAAA' 'PJ
  cyclic->AddTransition(hetero_atom,close_ring,'J'); // L6PJ
  cyclic->AddTransition(hetero_close_dash,close_ring,'J');
  cyclic->AddTransition(aromatics,close_ring,'J');
  cyclic->AddTransition(big_ring_dash_close,close_ring,'J'); // L-6-J


  cyclic->AddTransition(close_ring,close_ring,'&'); // instant closures

  return root; 
}

/* uses merges to be more specific on ring locant block */
void BuildWLNFSM2(FSMAutomata *wln){

  wln->AddState(); // create a root that points with epsiolon transitions to each block

  
  FSMState *cyclic_root = InsertCyclic(wln);
  wln->AddTransition(wln->root,cyclic_root,0); // e-transition


  // only one possible
  FSMState *cycle_accept = 0;
  for(unsigned int i=0;i<wln->num_states;i++){
    if(wln->states[i]->accept)
      cycle_accept = wln->states[i];
  }

  // handle locants, spiro and inline rings
  FSMState *locant_open = wln->AddState(false); 
  FSMState *locant_char = wln->AddState(false); 

  wln->AddTransition(locant_char,locant_char,'&'); // expansion

  FSMState *acyclic_ring_root = InsertAcyclic(wln);

  FSMState *inline_open = wln->AddState(false); 
  FSMState *inline_locant = wln->AddState(false); 
  wln->AddTransition(inline_open,inline_locant,' '); // make use of the epsilons here


  FSMState *spiro_open = wln->AddState(false); 
  FSMState *spiro_confirm = wln->AddState(false); 
  FSMState *spiro_locant = wln->AddState(false); 

  // out of line U- bonding 
  FSMState *out_double = wln->AddState(false); 

  FSMState *benzene = wln->AddState(true);

  for(unsigned int i=0;i<wln->num_states;i++){
    if(wln->states[i]->accept){
      wln->AddTransition(wln->states[i],locant_open,' ');
  
      if(wln->states[i] != cycle_accept){
        wln->AddTransition(wln->states[i],inline_open,'-');
        wln->AddTransition(wln->states[i],out_double,'U');
      }
    }
  }

  wln->AddTransition(out_double,inline_open,'-');
  
  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    wln->AddTransition(locant_open,locant_char,ch);
    wln->AddTransition(inline_locant,cyclic_root,ch); 
    wln->AddTransition(spiro_locant,cyclic_root,ch); 
  }

  

  // pi bonding
  wln->AddTransition(locant_open,locant_char,'0');
  wln->AddTransition(inline_locant,cyclic_root,'0'); 
  wln->AddTransition(spiro_locant,cyclic_root,'0'); 
  
  // make use of the epsilons here
  wln->AddTransition(locant_char,cyclic_root,0); 
  wln->AddTransition(locant_char,cyclic_root,'U');
  wln->AddTransition(locant_char,out_double,'U');
  wln->AddTransition(locant_char,acyclic_ring_root,0); 
  wln->AddTransition(locant_char,acyclic_ring_root,'U'); 
  
  wln->AddTransition(locant_char,inline_open,'-'); 

  wln->AddTransition(locant_char,spiro_open,'-'); 
  wln->AddTransition(spiro_open,spiro_confirm,'&'); 
  wln->AddTransition(spiro_confirm,spiro_locant,' '); 


  FSMState *non_cycles = InsertAcyclic(wln);
  wln->AddTransition(wln->root,non_cycles,0); // e-transition
  wln->AddTransition(wln->root,benzene,'R'); 
  wln->AddTransition(benzene,acyclic_ring_root,0);

  // handle benzene shorthand here
  wln->AddTransition(benzene,benzene,'&');
  wln->AddTransition(benzene,benzene,'R');
  wln->AddTransition(benzene,locant_open,' ');
  wln->AddTransition(benzene,inline_open,'-');

  wln->AddTransition(locant_char,benzene,'R');

  // ions are just repeats
  FSMState *ion = wln->AddState(false);
  // point all accepts at ion space
  for(unsigned int i=0;i<wln->num_states;i++){
    if(wln->states[i]->accept){
      wln->AddTransition(wln->states[i],ion,' ');
      wln->AddTransition(wln->states[i],benzene,'R');// this will jump out of acyclic
    }
  }
  wln->AddTransition(ion,wln->root,'&'); 
  
  return;
}


FSMAutomata * CreateWLNDFA(unsigned int node_size, unsigned int edge_size){
  FSMAutomata *wln = new FSMAutomata(node_size,edge_size);
  FSMAutomata *wlnDFA = 0;
  FSMAutomata *wlnMinimal = 0;

  BuildWLNFSM2(wln);
  wlnDFA = ConvertToDFA(wln);
  wlnMinimal = MinimiseDFA(wlnDFA);
  
  if(wlnMinimal)
    wlnMinimal->InitJumpTable();

  delete wln;
  delete wlnDFA;
  return wlnMinimal;
}


#endif