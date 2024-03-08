#include <cstdint>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <inttypes.h>

#include <string>

#include "rfsm.h"
#include "wlndfa.h"
#include "ppm.h"

const char *input;
unsigned int decode_zero_added = 0; 

#define NGRAM 4

#define FIXED_POINT_FRACTIONAL_BITS 16
#define FIXED_POINT_I_MAX INT32_MAX
#define FIXED_POINT_I_MIN INT32_MIN
#define FIXED_POINT_C INT32_C

/* use this for probability scaling */ 
unsigned int double_to_fixed(double input) {
  input *= FIXED_POINT_C(1) << FIXED_POINT_FRACTIONAL_BITS;
  input += (input < 0) ? -0.5 : 0.5;
  if (input >= (FIXED_POINT_I_MAX/2 + 1)*2.0) {
    return FIXED_POINT_I_MAX;
  }

  if (input - FIXED_POINT_I_MIN <= -1) {
    return FIXED_POINT_I_MIN;
  }

  return (unsigned int)input;
}

/* This also adds inf 1s to output */
void fread_string(unsigned char *ch, std::string &bitstream){
  (*ch) = 0; // set bits to zero
  unsigned int i=0;
  for(i=0;i<8;i++){
    if(i>=bitstream.size()){
      if(decode_zero_added)
        (*ch) ^= (1 << (7-i));
      else
        decode_zero_added = 1;
    }
    else if(bitstream[i] == 1)
      (*ch) ^= (1 << (7-i));
  }
  if(!bitstream.empty())
    bitstream.erase(0,i); 
}



bool WLNPPMCompressBuffer(const char *str, FSMAutomata *wlnmodel, std::string &bitstream, bool add_terminal){
  
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \n";

  wlnmodel->AssignEqualProbs(); // you moron
  
  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  unsigned int read = 0; 

  unsigned int low = 0; 
  unsigned int high = UINT32_MAX; // set all the bits to 11111...n 
  unsigned int underflow_bits = 0;
  
  unsigned int seen_context = 0;
  unsigned char lookback[NGRAM+1] = {0}; 

  unsigned char ch = *str;
  bool stop = false;
  
  Node *root = AllocateTreeNode('0', 0); 
  root->c = 1; 

  // the FSM is a fallback model, if PPM does not work, not intergrated with PPM itself. 
  // our order-1 is already incredibley good

  for(;;){

    // set an ngram look back for forward trie
    if(seen_context < NGRAM)
      lookback[seen_context++] = ch;
    else{      
      // escaping with no context is pointless
      for(unsigned int i=0;i<NGRAM-1;i++)
        lookback[i] = lookback[i+1]; 
      lookback[NGRAM-1] = ch; 
 
    }

    unsigned int T = 0;
    unsigned int edges = 0; 
    for(edge=state->transitions;edge;edge=edge->nxt){
      T+= edge->c;
      edges++; 
    }
    
    double test_prob = 0.0;
    for(unsigned int a=0;a<strlen(wln);a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, 'A', seen_context-1, strlen(wln)); 
      test_prob += prob; 
    }

    fprintf(stderr,"sum: %f\n",test_prob); 
    
    unsigned int Cc = 0; 
    unsigned int Cn = 0; 
    bool found = false;
    for(edge = state->transitions;edge;edge=edge->nxt){
      if(ch == edge->ch){
        Cn += edge->c;
        state = edge->dwn; // move the state;
        found = 1; 
        break;
      }
      else
        Cc += edge->c;
    }
    Cn += Cc;
    
    if(!found){
      fprintf(stderr,"Error: invalid wln notation\n");
      return false;
    }


// ################################################

    // standard arithmetic coder 32 bit int. 
    uint64_t range = ((uint64_t)high+1)-(uint64_t)low;
    uint64_t new_low = (uint64_t)low + (uint64_t)floor((range*Cc)/T); 
    uint64_t new_high = (uint64_t)low + (uint64_t)floor((range*Cn)/T);  

    // truncate down
    low = new_low;
    high = new_high;

    unsigned char lb = low & (1 << 31) ? 1:0;
    unsigned char hb = high & (1 << 31) ? 1:0;
    unsigned char lb2 = low & (1 << 30) ? 1:0;
    unsigned char hb2 = high & (1 << 30) ? 1:0;
    unsigned char ubit = lb ? 0:1;

    if(lb == hb){
      while(lb == hb){
        bitstream += lb;
  
        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;
        lb = low & (1 << 31) ? 1:0;
        hb = high & (1 << 31) ? 1:0;

        if(underflow_bits){
          for(unsigned int i=0;i<underflow_bits;i++)
            bitstream += ubit; 
          underflow_bits = 0;
        }
      }
    }
    else if (lb2 && !hb2){      
      low <<= 1;
      high <<= 1; 

      low ^= (1 << 31);
      high ^= (1 << 31);
      high ^= 1;
      
      underflow_bits++;
    }

// #################################################################################

    // PPM trie update happens after an encoding so the decompressor can keep up

    fprintf(stderr,"%s\n",lookback); 
    // update context trie
    if(!BuildContextTree(root, (const char*)lookback, seen_context)){
      fprintf(stderr,"Error: failure in building n-gram trie\n");
      return false;
    }

// #################################################################################

    if(stop)
      break;
    else
     ch = *(++str);
    
    if(!ch){
      if(!add_terminal)
        break;
      else
        stop = true;     
      ch = '\0';
    }



    read++; 
  }
  
  if(add_terminal)
    read++; 

  FILE *fp = fopen("./debug.dot","w");
  WriteDotFile(root, fp);
  fclose(fp);

  RReleaseTree(root);

  
  fprintf(stderr,"bit stream: %d/%d\n",bitstream.size(), read*8); 
  return true; 
}

bool WLNPPMDecompressBuffer(std::string &bitstream, FSMAutomata *wlnmodel){
  
  wlnmodel->AssignEqualProbs(); // you moron

  FSMState *state = wlnmodel->root;
  FSMEdge *edge = 0 ;
  
  unsigned int low = 0; 
  unsigned int high = UINT32_MAX;

  unsigned int enc_pos = 0; 
  unsigned int encoded = 0;
  unsigned char ch = 0;

  for(unsigned int i=0;i<4;i++){ // read 32 max into encoded
    fread_string(&ch, bitstream);
    for(int j=7;j>=0;j--){
      if(ch & (1 << j))  
        encoded ^= (1 << (31-enc_pos) );      
      
      enc_pos++;  
    }
  }
  
  // pre-load next char if relevent
  fread_string(&ch,bitstream);
  enc_pos = 0;

  for(;;){
    
    unsigned int T = 0; 
    unsigned int Cc = 0;
    unsigned int Cn = 0;
    for(edge=state->transitions;edge;edge=edge->nxt)
      T += edge->c; 
    
    uint64_t range = ((uint64_t)high+1)-(uint64_t)low;
    uint64_t scaled_sym = floor((T*(uint64_t)(encoded-low+1)-1)/range); // potential -1 herei
    
    for(edge=state->transitions;edge;edge=edge->nxt){
      Cn += edge->c;
      if(scaled_sym >= Cc && scaled_sym < Cn){  
        if(!edge->ch)
          return true;
        else
          fputc(edge->ch,stdout);
        
        state = edge->dwn;
        break;
      }
      else
        Cc += edge->c;
    }

    uint64_t new_low = (uint64_t)low + (uint64_t)floor((range*Cc)/T); 
    uint64_t new_high = (uint64_t)low + (uint64_t)floor((range*Cn)/T);  // should there be a minus 1 here?
                                                                          
    low = new_low;
    high = new_high;

    unsigned char lb = low & (1 << 31) ? 1:0;
    unsigned char hb = high & (1 << 31) ? 1:0;
    unsigned char lb2 = low & (1 << 30) ? 1:0;
    unsigned char hb2 = high & (1 << 30) ? 1:0;

    if(lb == hb){
      while(lb == hb){
  
        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;

        encoded <<= 1;
        if(ch & (1 << 7))
          encoded ^= 1;

        ch <<= 1;
        enc_pos++;

        if(enc_pos == 8){ // read the next block
          fread_string(&ch, bitstream);
          enc_pos = 0;
        } 

        lb = low & (1 << 31) ? 1:0;
        hb = high & (1 << 31) ? 1:0;
      }
    }
    else if (lb2 && !hb2){      
      unsigned int p = 0;
      unsigned int encoded_shift = 0; 
      for(int j=31;j>=0;j--){
        if(j != 30){
          if (encoded & (1 << j))
            encoded_shift ^= (1 << (31-p) );
          p++;
        }
      }
      
      if(ch & (1 << 7))
        encoded_shift ^= 1;

      ch <<= 1;
      enc_pos++;

      if(enc_pos == 8){ // read the next block
        fread_string(&ch, bitstream);
        enc_pos = 0;
      } 
      
      encoded = encoded_shift; // set the bit to spliced

      low <<= 1;
      high <<= 1; 

      low ^= (1 << 31);
      high ^= (1 << 31);
      high ^= 1;
    }
  }
  
  return true; 
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  input = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          exit(1); 
      }
    }
    else{
      switch(j++){
        case 0:
          input = ptr; 
          break;
        default:
          fprintf(stderr,"Error: multiple files not currently supported\n");
          exit(1);
      }
    }
  }

  if(!input){
    fprintf(stderr,"Error: no input file given\n");
    exit(1);
  }

  return;
}
int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
    
  double test = 0.2500;
  fprintf(stderr,"%f = %d\n",test,double_to_fixed(test)); 

  FSMAutomata *wlnmodel = CreateWLNDFA(REASONABLE*2,REASONABLE*4); // build the model 

  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept){
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\0');  
    }
  }
  wlnmodel->AssignEqualProbs(); // initalise the order -1 model
  

  bool ending = false;
  // perhaps need to add the escape sequence as a transition. 
  std::string bitstream; 
  if(!WLNPPMCompressBuffer(input, wlnmodel, bitstream,ending))
    return 1; 
   
  if(ending)
    WLNPPMDecompressBuffer(bitstream, wlnmodel);
  
  delete wlnmodel;
  return 0;
}
