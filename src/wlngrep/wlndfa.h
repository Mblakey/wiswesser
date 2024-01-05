#ifndef WLN_DFA_H
#define WLN_DFA_H

#include "rfsm.h"
#include "rconvert.h"

/* will purge and delete to merge into the main fsm */
FSMAutomata *CreateAcyclic(){
  FSMAutomata *acyclic = new FSMAutomata(REASONABLE,REASONABLE); 

  FSMState *root = acyclic->AddState(true);

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
   

  
  return acyclic;
}

/* for ease of building, this will build an NFA, which will subset down to a DFA */
void BuildWLNFSM(FSMAutomata *wlnNFA){

  FSMState *root = wlnNFA->root;
  if(!root){
    fprintf(stderr,"Error: root not initialised\n");
    return;
  }

  FSMState *first_allowed = wlnNFA->AddState(true); 
  FSMState *digits = wlnNFA->AddState(true);

  // set up digits state, cant be zero on entrance
  for(unsigned char ch = '1';ch <= '9';ch++){
    wlnNFA->AddTransition(root,digits,ch);
    wlnNFA->AddTransition(first_allowed,digits,ch);
  }
    
  // 0 allowed on repeating digit state
  for(unsigned char ch = '0';ch <= '9';ch++)
    wlnNFA->AddTransition(digits,digits,ch);

  // set up first set of transitions 
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
        wlnNFA->AddTransition(root,first_allowed,ch);
        wlnNFA->AddTransition(first_allowed,first_allowed,ch);
        wlnNFA->AddTransition(digits,first_allowed,ch);
        break;
    }
  }


  // allow branching '& here
  FSMState *branch = wlnNFA->AddState(true);
  wlnNFA->AddTransition(first_allowed,branch,'&'); 
  wlnNFA->AddTransition(branch,branch,'&'); // allow the syntax of repeated pops 
  wlnNFA->AddTransition(digits,branch,'&'); 

  // return after branching 
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
      case 'H':
      case 'R': // seperate out benzene here
        break;
      
      default:
        wlnNFA->AddTransition(branch,first_allowed,ch);
        break;
    }
  }

  // cant be zero on entrance
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(branch,digits,ch);
  

  // what can double bond? these are not accepts
  FSMState *db_only = wlnNFA->AddState(false);

  wlnNFA->AddTransition(db_only,db_only,'U');
  wlnNFA->AddTransition(first_allowed,db_only,'U');
  wlnNFA->AddTransition(digits,db_only,'U');
  wlnNFA->AddTransition(branch,db_only,'U');

  // double bonds
  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    switch(ch){
      case 'L':
      case 'T':
      case 'D':
      case 'J':
      case 'A':
      case ' ':
      case '-':
      case '/':
      case 'U':
      case 'H':
      case 'C':
      case 'R':
        break;
      
      default:
        wlnNFA->AddTransition(db_only,first_allowed,ch);
        break;
    }
  }

  // cant be zero on entrance
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(db_only,digits,ch);
  

  // Dash element specification for elements and hypervalence
  FSMState *element_dash_start = wlnNFA->AddState(false);
  FSMState *element_dash_end = wlnNFA->AddState(true);

  wlnNFA->AddTransition(root,element_dash_start,'-');
  wlnNFA->AddTransition(first_allowed,element_dash_start,'-');
  wlnNFA->AddTransition(db_only,element_dash_start,'-');
  wlnNFA->AddTransition(digits,element_dash_start,'-');
  wlnNFA->AddTransition(branch,element_dash_start,'-');

  FSMState *char_1 = wlnNFA->AddState(false);
  FSMState *char_2 = wlnNFA->AddState(false);

  // allow any letter
  for(unsigned char ch = 'A'; ch <= 'Z';ch++){
    wlnNFA->AddTransition(element_dash_start,char_1,ch);
    wlnNFA->AddTransition(char_1,char_2,ch);
  }
  wlnNFA->AddTransition(char_2,element_dash_end,'-');
  
  FSMState *hypo_char = wlnNFA->AddState(false);
  // hyper valent single elements
  wlnNFA->AddTransition(element_dash_start,hypo_char,'P');
  wlnNFA->AddTransition(element_dash_start,hypo_char,'S');
  wlnNFA->AddTransition(element_dash_start,hypo_char,'E');
  wlnNFA->AddTransition(element_dash_start,hypo_char,'F');
  wlnNFA->AddTransition(element_dash_start,hypo_char,'G');
  wlnNFA->AddTransition(element_dash_start,hypo_char,'I');
  wlnNFA->AddTransition(element_dash_start,hypo_char,'E');
  wlnNFA->AddTransition(hypo_char,element_dash_end,'-');
  
  // connect dash end to various states
  // return after branching 
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
      case 'R': // seperate out benzene here
      //case 'H': H is allowed here
        break;
      
      default:
        wlnNFA->AddTransition(element_dash_end,first_allowed,ch);
        break;
    }
  }

  // cant be zero on entrance
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(element_dash_end,digits,ch);

  wlnNFA->AddTransition(element_dash_end,db_only,'U');
  wlnNFA->AddTransition(element_dash_end,branch,'&');


  // IONS
  FSMState *ion_space = wlnNFA->AddState(false);
  FSMState *ion_ampersand = wlnNFA->AddState(false);
  wlnNFA->AddTransition(ion_space,ion_ampersand,'&');

  // ions cannot come from root or off a double bond
  wlnNFA->AddTransition(first_allowed,ion_space,' ');
  wlnNFA->AddTransition(digits,ion_space,' ');
  wlnNFA->AddTransition(branch,ion_space,' ');
  wlnNFA->AddTransition(element_dash_end,ion_space,' ');

  // connect the ampersand to various
  // return after branching 
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
      case 'H':
        break;
      
      default:
        wlnNFA->AddTransition(ion_ampersand,first_allowed,ch);
        break;
    }
  }

  // cant be zero on entrance
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(ion_ampersand,digits,ch);

  wlnNFA->AddTransition(ion_ampersand,element_dash_start,'-');

  // expand ions to gather charge post notation 
  FSMState *charge_start = wlnNFA->AddState(false);
  FSMState *charge_end = wlnNFA->AddState(true);
  FSMState *charge_slash = wlnNFA->AddState(false);

  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(ion_ampersand,charge_start,ch);

  for(unsigned char ch = '0';ch <= '9';ch++)
    wlnNFA->AddTransition(charge_start,charge_start,ch);

  wlnNFA->AddTransition(charge_start,charge_slash,'/');

  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(charge_slash,charge_end,ch);

  for(unsigned char ch = '0';ch <= '9';ch++)
    wlnNFA->AddTransition(charge_end,charge_end,ch);

  wlnNFA->AddTransition(charge_end, ion_space,' ');


  /* ######### BENZENE DEFINITIONS ######### */

  FSMState *benzene = wlnNFA->AddState(true);
  wlnNFA->AddTransition(root,benzene,'R');
  wlnNFA->AddTransition(benzene,benzene,'R');

  wlnNFA->AddTransition(first_allowed,benzene,'R');
  wlnNFA->AddTransition(digits,benzene,'R');
  wlnNFA->AddTransition(db_only,benzene,'R');
  wlnNFA->AddTransition(element_dash_end,benzene,'R');
  wlnNFA->AddTransition(branch,benzene,'R');

  // benzene can be used inline 
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
      case 'H':
        break;
      
      default:
        wlnNFA->AddTransition(benzene,first_allowed,ch);
        break;
    }
  }

  // cant be zero on entrance
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(benzene,digits,ch);
  
  wlnNFA->AddTransition(benzene,branch,'&');
  wlnNFA->AddTransition(benzene,element_dash_start,'-');
  wlnNFA->AddTransition(benzene,db_only,'U');
  wlnNFA->AddTransition(benzene,ion_space,' ');

  // locants 
  FSMState *locant_space = wlnNFA->AddState(false);
  FSMState *locant_ch    = wlnNFA->AddState(true);

  wlnNFA->AddTransition(benzene,locant_space,' ');

  // allow any transition
  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    wlnNFA->AddTransition(locant_space,locant_ch,ch);

  wlnNFA->AddTransition(locant_space,locant_ch,'0'); // metallocene

  // point the locant first allowed, but include U for outer ring unsaturations 
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
      case 'H':
        break;
      
      default:
        wlnNFA->AddTransition(locant_ch,first_allowed,ch);
        break;
    }
  }

  // cant be zero on entrance
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(locant_ch,digits,ch);

  wlnNFA->AddTransition(locant_ch,element_dash_start,'-');
  wlnNFA->AddTransition(locant_ch,branch,'&');
  wlnNFA->AddTransition(locant_ch,db_only,'U');

  wlnNFA->AddTransition(first_allowed,locant_space,' ');
  wlnNFA->AddTransition(digits,locant_space,' ');
  wlnNFA->AddTransition(branch,locant_space,' ');
  wlnNFA->AddTransition(element_dash_end,locant_space,' ');

  /* ######### CYCLIC DEFINITIONS ######### */
  
  FSMState *open_ring = wlnNFA->AddState(false);
  FSMState *close_ring = wlnNFA->AddState(true);

  wlnNFA->AddTransition(root,open_ring,'L');
  wlnNFA->AddTransition(root,open_ring,'T');

  // all locants ...
  wlnNFA->AddTransition(close_ring,locant_space,' ');
  
  // all ring ions
  wlnNFA->AddTransition(close_ring,ion_space,' ');
  wlnNFA->AddTransition(ion_ampersand,open_ring,'L');
  wlnNFA->AddTransition(ion_ampersand,open_ring,'T');

  // immediate closures 
  wlnNFA->AddTransition(close_ring,close_ring,'&');

  // handle all the ring digits 

  FSMState *ring_digits = wlnNFA->AddState(false); // not an accept now
  for(unsigned char ch = '0';ch <= '9';ch++)
    wlnNFA->AddTransition(ring_digits,ring_digits,ch); // zero here allows the matching if metallocens


  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(open_ring,ring_digits,ch);

  wlnNFA->AddTransition(ring_digits,close_ring,'J');

  FSMState *big_ring_dash_open = wlnNFA->AddState(false);
  FSMState *big_ring_dash_close = wlnNFA->AddState(false);
  FSMState *big_ring_digits = wlnNFA->AddState(false);

  wlnNFA->AddTransition(open_ring,big_ring_dash_open,'-'); // L-


  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(big_ring_dash_open,big_ring_digits,ch); // L-6
  
  for(unsigned char ch = '0';ch <= '9';ch++)
    wlnNFA->AddTransition(big_ring_digits,big_ring_digits,ch); // L-666... 

  wlnNFA->AddTransition(big_ring_digits,big_ring_dash_close,'-'); // L-666-
  
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(big_ring_dash_close,ring_digits,ch); // L-666-6
  
  wlnNFA->AddTransition(ring_digits,big_ring_dash_open,'-'); // L6-
  wlnNFA->AddTransition(big_ring_dash_close,big_ring_dash_open,'-'); // L-6--

  wlnNFA->AddTransition(big_ring_dash_close,close_ring,'J'); // L-6-J


  // POLY CYCLICS RING NODES

  FSMState *digit_space = wlnNFA->AddState(false);
  FSMState *digit_locant = wlnNFA->AddState(false);

  wlnNFA->AddTransition(digit_locant,digit_locant,'&'); //  L E&6
  wlnNFA->AddTransition(digit_locant,digit_locant,'-'); //  L E&6

  wlnNFA->AddTransition(digit_locant,digit_space,' '); //  forward bridge notation

  wlnNFA->AddTransition(open_ring,digit_space,' '); // L' '

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    wlnNFA->AddTransition(digit_space,digit_locant,ch); // 'L' 'A

  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(digit_locant,ring_digits,ch); // L' 'A6..

  wlnNFA->AddTransition(digit_locant,big_ring_dash_open,'-'); // L' 'A-6-..

  wlnNFA->AddTransition(ring_digits,digit_space,' '); // L666' 'A
  wlnNFA->AddTransition(big_ring_dash_close,digit_space,' '); // L66-6-' 'A

  // MULTI CYCLICS NOTATION

  FSMState *multi_space = wlnNFA->AddState(false);
  FSMState *multi_digit = wlnNFA->AddState(false);
  FSMState *multi_locants = wlnNFA->AddState(false);

  wlnNFA->AddTransition(digit_locant,multi_space,' '); //  forward bridge notation
  
  wlnNFA->AddTransition(ring_digits,multi_space,' '); // L6' '
  wlnNFA->AddTransition(big_ring_dash_close,multi_space,' '); // L-6-' '
  
  for(unsigned char ch = '1';ch <= '9';ch++)
    wlnNFA->AddTransition(multi_space,multi_digit,ch); //  L6' '2

  for(unsigned char ch = 'A';ch <= 'Z';ch++){
    wlnNFA->AddTransition(multi_digit,multi_locants,ch); //  L6' '2A
    wlnNFA->AddTransition(multi_locants,multi_locants,ch); //  L6' '2AAAA...
  }

  FSMState *multi_size_space = wlnNFA->AddState(false);
  FSMState *multi_size = wlnNFA->AddState(false);

  wlnNFA->AddTransition(multi_locants,multi_locants,'&'); //  L E&6
  wlnNFA->AddTransition(multi_locants,multi_locants,'-'); //  L E-6 // broken locant 

  wlnNFA->AddTransition(multi_locants,multi_size_space,' '); //  L6' '2AAAA' '
  

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    wlnNFA->AddTransition(multi_size_space,multi_size,ch); //  L6' '2AAAA...P    

  // give size expansion ability
  wlnNFA->AddTransition(multi_size,multi_size,'&'); //  L E&6
  wlnNFA->AddTransition(multi_size,close_ring,'J'); //  L6' '2AAAA' 'PJ

  // Branching rings notation 
  FSMState *pair_slash = wlnNFA->AddState(false);
  FSMState *pair_loc_a = wlnNFA->AddState(false); 
  FSMState *pair_loc_b = wlnNFA->AddState(false);  
  
  wlnNFA->AddTransition(ring_digits,pair_slash,'/'); //  L6/
  wlnNFA->AddTransition(big_ring_dash_close,pair_slash,'/'); //L-6-/

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    wlnNFA->AddTransition(pair_slash,pair_loc_a,ch); //  L6' '2AAAA...P    

  wlnNFA->AddTransition(pair_loc_a,pair_loc_a,'&'); //  L E&6
  wlnNFA->AddTransition(pair_loc_a,pair_loc_a,'-'); //  L E&6

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    wlnNFA->AddTransition(pair_loc_a,pair_loc_b,ch); //  L6' '2AAAA...P

  wlnNFA->AddTransition(pair_loc_b,pair_loc_b,'&'); //  L E&6
  wlnNFA->AddTransition(pair_loc_b,pair_loc_b,'-'); //  L E&6

  wlnNFA->AddTransition(pair_loc_b,pair_slash,'/');
  wlnNFA->AddTransition(pair_loc_b,digit_space,' ');
  wlnNFA->AddTransition(pair_loc_b,multi_space,' ');

  // RING HETERO ATOMS

  // handle all the ring heteroatomic assignments, (must come from digits)  
  FSMState *hetero_space =  wlnNFA->AddState(false);
  FSMState *hetero_locant = wlnNFA->AddState(false);
  FSMState *hetero_atom =   wlnNFA->AddState(false);

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
        wlnNFA->AddTransition(hetero_locant,hetero_atom,ch);
        wlnNFA->AddTransition(hetero_atom,hetero_atom,ch); // L6PPP...
        wlnNFA->AddTransition(ring_digits,hetero_atom,ch); // L6P
        wlnNFA->AddTransition(big_ring_dash_close,hetero_atom,ch); // L-6-P
        break; // direct connections
    }
  }

  wlnNFA->AddTransition(hetero_atom,close_ring,'J'); // L6PJ
  wlnNFA->AddTransition(hetero_atom,hetero_space,' '); // L6P A...

  wlnNFA->AddTransition(ring_digits,hetero_space,' '); // L6 AP...
  wlnNFA->AddTransition(big_ring_dash_close,hetero_space,' '); // L-6-' '.

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    wlnNFA->AddTransition(hetero_space,hetero_locant,ch); // L6' 'A..

  wlnNFA->AddTransition(hetero_locant,hetero_locant,'&'); // L6P A&...
  

  // multi atom attachment 
  wlnNFA->AddTransition(multi_size,hetero_space,' '); //  L6' '2AAAA' 'PJ

  // allow briding and locant expansion
  wlnNFA->AddTransition(hetero_locant,hetero_space,' '); 


  // AROMATICS
  FSMState *aromatics = wlnNFA->AddState(false);
  wlnNFA->AddTransition(aromatics,aromatics,'&');
  wlnNFA->AddTransition(aromatics,aromatics,'T'); // L ... &T&T&TJ
  wlnNFA->AddTransition(aromatics,close_ring,'J');

  wlnNFA->AddTransition(big_ring_dash_close,aromatics,'&');
  wlnNFA->AddTransition(big_ring_dash_close,aromatics,'T');

  wlnNFA->AddTransition(ring_digits,aromatics,'&');
  wlnNFA->AddTransition(ring_digits,aromatics,'T');

  wlnNFA->AddTransition(hetero_atom,aromatics,'&');
  wlnNFA->AddTransition(hetero_atom,aromatics,'T');

  // multi atom attachment 
  wlnNFA->AddTransition(multi_size,aromatics,' '); //  L6' '2AAAA' 'PJ
  wlnNFA->AddTransition(multi_size,aromatics,'T');
  wlnNFA->AddTransition(multi_size,aromatics,'&');

  // Recursive Ring defintions

  FSMState *inline_ring = wlnNFA->AddState(false);
  FSMState *inline_space = wlnNFA->AddState(false);
  FSMState *inline_locant = wlnNFA->AddState(false);

  wlnNFA->AddTransition(inline_ring,inline_ring,'&'); // allow spiro
  
  wlnNFA->AddTransition(first_allowed,inline_ring,'-');
  wlnNFA->AddTransition(element_dash_end,inline_ring,'-');
  wlnNFA->AddTransition(digits,inline_ring,'-');
  wlnNFA->AddTransition(branch,inline_ring,'-');
  wlnNFA->AddTransition(db_only,inline_ring,'-');

  wlnNFA->AddTransition(locant_ch,inline_ring,'-');

  wlnNFA->AddTransition(inline_ring,inline_space,' ');

  for(unsigned char ch = 'A';ch <= 'Z';ch++)
    wlnNFA->AddTransition(inline_space,inline_locant,ch); 
  wlnNFA->AddTransition(inline_space,inline_locant,'0'); // metallocene 
  
  wlnNFA->AddTransition(inline_locant,inline_locant,'&');
  
  wlnNFA->AddTransition(inline_locant,open_ring,'L');
  wlnNFA->AddTransition(inline_locant,open_ring,'T');

}


FSMAutomata * CreateWLNDFA(){
  FSMAutomata *wln = new FSMAutomata(REASONABLE,REASONABLE);
  FSMAutomata *wlnDFA = 0;
  FSMAutomata *wlnMinimal = 0;

  wln->AddState(); 
  BuildWLNFSM(wln); 

  wlnDFA = ConvertToDFA(wln);

  wlnMinimal = MinimiseDFA(wlnDFA);
  
  if(wlnMinimal)
    wlnMinimal->InitJumpTable();


  FSMAutomata *test = CreateAcyclic();
  FSMAutomata *testDFA = ConvertToDFA(test);
  FSMAutomata *testMinimal = MinimiseDFA(testDFA);

  testMinimal->InitJumpTable();
  return testMinimal;

  delete wln;
  delete wlnDFA;
  return wlnMinimal;
}


#endif