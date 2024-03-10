#include <cstdint>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>

#include "rfsm.h"
#include "ppm.h"

#include "wlnzip.h"

unsigned int decode_zero_added = 0; 

#define NGRAM 5
#define PPM 1
#define ALPHABET 42
#define TERMINATE 'x'

#define FIXED_POINT_FRACTIONAL_BITS 16
#define FIXED_POINT_I_MAX INT16_MAX
#define FIXED_POINT_I_MIN INT16_MIN
#define FIXED_POINT_C INT16_C

/* use this for probability scaling */ 
unsigned short double_to_fixed(double input) {
  input *= FIXED_POINT_C(1) << FIXED_POINT_FRACTIONAL_BITS;
  input += (input < 0) ? -0.5 : 0.5;
  if (input >= (FIXED_POINT_I_MAX/2 + 1)*2.0) {
    return FIXED_POINT_I_MAX;
  }

  if (input - FIXED_POINT_I_MIN <= -1) {
    return FIXED_POINT_I_MIN;
  }
  
  return (unsigned short int)input;
}


/* linked list bitstream, just for single string encoding
 * purposes 
*/
void ReadStream(BitStream *stream){
  BitStream *bit = stream;

  while(bit){
    fprintf(stderr,"%d",bit->b ? 1:0);
    bit = bit->nxt; 
  } 
  fprintf(stderr,"\n");
}

void DeleteStream(BitStream *stream){
  BitStream *bit = stream; 
  while(bit){
    BitStream *r = bit; 
    bit = bit->nxt; 
    free(r); 
  }
}

void Append(unsigned char b, BitStream *stream){
  BitStream *nb = (BitStream*)malloc(sizeof(BitStream)); 
  nb->b = b; 
  nb->nxt = 0; 

  BitStream *bit = stream; 
  while(bit->nxt)
    bit = bit->nxt; 

  bit->nxt = nb;   
}


void OutputStream(std::string &stream){
  unsigned int  pos = 0; 
  unsigned char ch  = 0; 
  for(unsigned int i=0;i<stream.size();i++){
    if(stream[i])
      ch ^= (1 << (7-pos));
    
    pos++; 
    if(pos == 8){
      fputc(ch, stdout); 
      pos = 0; 
      ch = 0; 
    }
  }
  
  while(pos < 8){
    if(decode_zero_added)
     ch ^= (1 << (7-pos));
    else
      decode_zero_added = true;
    pos++;

    fprintf(stderr,"adding ones\n"); 
  }

  fputc(ch, stdout);
}

BitStream* WLNPPMCompressBuffer(const char *str, FSMAutomata *wlnmodel, unsigned char escape_mode,bool add_terminal){ 
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx"; // TERMINATE = 'x'
  wlnmodel->AssignEqualProbs();
  
  BitStream *bitstream = (BitStream*)malloc(sizeof(BitStream)); 
  bitstream->b = 0; 
  bitstream->nxt = 0; 


  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  unsigned int read = 0; 

  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX; // set all the bits to 11111...n 
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

    bool found = 0; 
    unsigned int T = 0;
    unsigned int Cc = 0; 
    unsigned int Cn = 0; 

#if PPM
    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_mode, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob); 
      T += fp_prob;
    }

    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_mode, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob); 
    
      if(ch == wln[a]){
        Cn += fp_prob;
        found = 1; 
        break;
      }
      else
        Cc += fp_prob; 
    }
    Cn += Cc; 
#else 
    unsigned int edges = 0; 
    for(edge=state->transitions;edge;edge=edge->nxt){
      T+= edge->c;
      edges++; 
    }
 
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
      fprintf(stderr,"Error: invalid wln notation - %c\n",ch);
      return false;
    }
#endif


    // standard arithmetic coder 16 bit int precision. 
    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int new_low = (unsigned int)low + (unsigned int)floor((range*Cc)/T); 
    unsigned int new_high = (unsigned int)low + (unsigned int)floor((range*Cn)/T)-1; // yes!  

    // truncate down
    low = new_low;
    high = new_high;

    unsigned char lb = low & (1 << 15) ? 1:0;
    unsigned char hb = high & (1 << 15) ? 1:0;
    unsigned char lb2 = low & (1 << 14) ? 1:0;
    unsigned char hb2 = high & (1 << 14) ? 1:0;
    unsigned char ubit = lb ? 0:1;

    if(lb == hb){
      while(lb == hb){
        Append(lb,bitstream);
  
        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;
        lb = low & (1 << 15) ? 1:0;
        hb = high & (1 << 15) ? 1:0;

        if(underflow_bits){
          for(unsigned int i=0;i<underflow_bits;i++)
            Append(ubit, bitstream);  
          underflow_bits = 0;
        }
      }
    }
    else if (lb2 && !hb2){      
      low <<= 1;
      high <<= 1; 

      low ^= (1 << 15);
      high ^= (1 << 15);
      high ^= 1;
      
      underflow_bits++;
    }


// #################################################################################

// PPM trie update happens after an encoding so the decompressor can keep up
// update context trie

    if(seen_context < NGRAM)
      lookback[seen_context++] = ch;
    else{      
      // escaping with no context is pointless
      for(unsigned int i=0;i<NGRAM-1;i++)
        lookback[i] = lookback[i+1]; 
      lookback[NGRAM-1] = ch; 
    }

    if(!BuildContextTree(root, (const char*)lookback, seen_context)){
      fprintf(stderr,"Error: failure in building n-gram trie\n");
      return 0;
    }

 // #################################################################################

    if(stop)
      break;
    else
      ch = *(++str);

    if(!ch){
      if(!add_terminal)
        break;
      else{
        stop = true;
        ch = TERMINATE; // special character for WLN ending
      }
    }

    read++; 
  }
  

  if(add_terminal)
    read++; 

  RReleaseTree(root);
  return bitstream; 
}



bool WLNPPMDecompressBuffer(BitStream *bitstream, FSMAutomata *wlnmodel, unsigned char escape_type){
  
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx";
  wlnmodel->AssignEqualProbs(); // you moron

  FSMState *state = wlnmodel->root;
  FSMEdge *edge = 0 ;

  BitStream *bit = bitstream->nxt; // get first bit
  
  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX;
  unsigned short int encoded = 0;
  unsigned int enc_pos = 0; 

  unsigned int seen_context = 0;
  unsigned char lookback[NGRAM+1] = {0}; 

  Node *root = AllocateTreeNode('0', 0); 
  root->c = 1; 

  while(bit && enc_pos < 16){ // read 16 max into encoded
    if(bit->b)
      encoded ^= (1 << (15-enc_pos) );

    enc_pos++;
    bit = bit->nxt;
  }
  
  // pre-load next char if relevent
  enc_pos = 0;
  for(;;){
    
    unsigned int T = 0; 
    unsigned int Cc = 0;
    unsigned int Cn = 0;

#if PPM
    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_type, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob); 
      T+= fp_prob;
    }

    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range); 
    
    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_type, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob); 
    
      Cn += fp_prob; 
      if(scaled_sym >= Cc && scaled_sym < Cn){
        if(wln[a] == TERMINATE)
          return true;
        else{
          // add the char to the buffer, and then build the trie
          fputc(wln[a],stdout);
          
          // build the trie immediately
          if(seen_context < NGRAM)
            lookback[seen_context++] = wln[a];
          else{   
            for(unsigned int i=0;i<NGRAM-1;i++)
              lookback[i] = lookback[i+1]; 
            lookback[NGRAM-1] = wln[a];  
          }
          
          if(!BuildContextTree(root, (const char*)lookback, seen_context)){
            fprintf(stderr,"Error: failure in building n-gram trie\n");
            return false;
          }

        }
        break;
      }
      else
        Cc += fp_prob; 
    }
   
#else
    for(edge=state->transitions;edge;edge=edge->nxt)
      T += edge->c; 
    
    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range); // potential -1 herei
    
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
#endif

    unsigned int new_low = (unsigned int)low + (unsigned int)floor((range*Cc)/T); 
    unsigned int new_high = (unsigned int)low + (unsigned int)floor((range*Cn)/T)-1;  // should there be a minus 1 here?
    
    // truncate down
    low = new_low;
    high = new_high;

    unsigned char lb = low & (1 << 15) ? 1:0;
    unsigned char hb = high & (1 << 15) ? 1:0;
    unsigned char lb2 = low & (1 << 14) ? 1:0;
    unsigned char hb2 = high & (1 << 14) ? 1:0;

    if(lb == hb){
      while(lb == hb){
  
        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;

        encoded <<= 1;

        if(bit){
          if(bit->b)
            encoded ^= 1; 
          bit = bit->nxt;
        }              
        else if(!decode_zero_added)
          decode_zero_added = 1;
        else
          encoded ^= 1;

        lb = low & (1 << 15) ? 1:0;
        hb = high & (1 << 15) ? 1:0;
      }
    }
    else if (lb2 && !hb2){      
      enc_pos = 0;
      unsigned short int encoded_shift = 0; 
      for(int j=15;j>=0;j--){
        if(j != 14){
          if (encoded & (1 << j))
            encoded_shift ^= (1 << (15-enc_pos) );
          enc_pos++;
        }
      }


      if(bit){
        if(bit->b)
          encoded_shift ^= 1; 

        bit = bit->nxt;
      }
      else if(!decode_zero_added)
        decode_zero_added = 1;
      else
        encoded_shift ^= 1; // inf ones if needed to end the message

      encoded = encoded_shift; // set the bit to spliced

      low <<= 1;
      high <<= 1; 

      low ^= (1 << 15);
      high ^= (1 << 15);
      high ^= 1;
    }
  }
  

  DeleteStream(bitstream); 
  RReleaseTree(root); 
  return true; 
}



bool WLNPPMCompressFile(FILE *ifp, FSMAutomata *wlnmodel, unsigned char escape_type){
  std::string bitstream; 

  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx"; // TERMINATE = 'x'
  wlnmodel->AssignEqualProbs(); // you moron
  
  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  unsigned int read = 0; 

  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX; // set all the bits to 11111...n 
  unsigned int underflow_bits = 0;
  
  unsigned int seen_context = 0;
  unsigned char lookback[NGRAM+1] = {0}; 

  unsigned char ch = 0;
  bool stop = false;
  
  Node *root = AllocateTreeNode('0', 0); 
  root->c = 1; 

  if(!fread(&ch,sizeof(unsigned char),1,ifp)){
    fprintf(stderr,"Error: no data in file\n"); 
    return false;
  }

  for(;;){

    bool found = 0; 
    unsigned int T = 0;
    unsigned int Cc = 0; 
    unsigned int Cn = 0; 

#if PPM
    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_type, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob); 
      T+= fp_prob;
    }

    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_type, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob); 
    
      if(ch == wln[a]){
        Cn += fp_prob;
        found = 1; 
        break;
      }
      else
        Cc += fp_prob; 
    }
    Cn += Cc; 
#else 
    unsigned int edges = 0; 
    for(edge=state->transitions;edge;edge=edge->nxt){
      T+= edge->c;
      edges++; 
    }
 
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
      fprintf(stderr,"Error: invalid wln notation - char: %c, accept?: %d\n",ch, state->accept);
      return false;
    }
#endif

// ################################################

    // standard arithmetic coder 32 bit int. 
    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int new_low = (unsigned int)low + (unsigned int)floor((range*Cc)/T); 
    unsigned int new_high = (unsigned int)low + (unsigned int)floor((range*Cn)/T)-1;  

    // truncate down
    low = new_low;
    high = new_high;

    unsigned char lb = low & (1 << 15) ? 1:0;
    unsigned char hb = high & (1 << 15) ? 1:0;
    unsigned char lb2 = low & (1 << 14) ? 1:0;
    unsigned char hb2 = high & (1 << 14) ? 1:0;
    unsigned char ubit = lb ? 0:1;

    if(lb == hb){
      while(lb == hb){
        bitstream += lb;
  
        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;
        lb = low & (1 << 15) ? 1:0;
        hb = high & (1 << 15) ? 1:0;

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

      low ^= (1 << 15);
      high ^= (1 << 15);
      high ^= 1;
      
      underflow_bits++;
    }


// #################################################################################

// the order of these build --> add arguements seem to lead to better or worse compression?

// PPM trie update happens after an encoding so the decompressor can keep up
// update context trie
    if(seen_context < NGRAM)
      lookback[seen_context++] = ch;
    else{      
      // escaping with no context is pointless
      for(unsigned int i=0;i<NGRAM-1;i++)
        lookback[i] = lookback[i+1]; 
      lookback[NGRAM-1] = ch; 
    }

    if(!BuildContextTree(root, (const char*)lookback, seen_context)){
      fprintf(stderr,"Error: failure in building n-gram trie\n");
      return false;
    }

 // #################################################################################

    if(stop)
      break;
    else if(!fread(&ch,sizeof(unsigned char),1,ifp)){
      stop = true;
      ch = TERMINATE; // special character for WLN ending
    }

    read++; 
  }
  
  OutputStream(bitstream); 

  fprintf(stderr,"bit stream: %d/%d = %f\n",bitstream.size(), read*8, (read*8)/(double)bitstream.size()); 
  RReleaseTree(root);  
  return true;
}




bool WLNPPMDecompressFile(FILE *ifp, FSMAutomata *wlnmodel, unsigned char escape_type){
  
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx";
  wlnmodel->AssignEqualProbs(); // you moron

  FSMState *state = wlnmodel->root;
  FSMEdge *edge = 0 ;
  
  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX;
  unsigned short int encoded = 0;

  unsigned int enc_pos = 0; 
  unsigned char ch = 0;

  unsigned int seen_context = 0;
  unsigned char lookback[NGRAM+1] = {0}; 

  Node *root = AllocateTreeNode('0', 0); 
  root->c = 1; 

  for(unsigned int i=0;i<2;i++){ // read 16 max into encoded
    fread(&ch,sizeof(unsigned char),1,ifp);
    for(int j=7;j>=0;j--){
      if(ch & (1 << j))  
        encoded ^= (1 << (15-enc_pos) );      
      
      enc_pos++;  
    }
  }
  
  // pre-load next char
  fread(&ch,sizeof(unsigned char),1,ifp);
  enc_pos = 0;
  for(;;){
    
    unsigned int T = 0; 
    unsigned int Cc = 0;
    unsigned int Cn = 0;

#if PPM
    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_type, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob); 
      T+= fp_prob;
    }

    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range);
    
    for(unsigned int a=0;a<ALPHABET;a++){
      double prob = PredictPPM((const char*)lookback, wln[a], root, escape_type, seen_context, ALPHABET); 
      unsigned short int fp_prob = double_to_fixed(prob);     
      Cn += fp_prob; 
      
      if(scaled_sym >= Cc && scaled_sym < Cn){
        if(wln[a] == TERMINATE)
          return true;
        else{
          // add the char to the buffer, and then build the trie
          fputc(wln[a],stdout);
          
          // build the trie immediately
          if(seen_context < NGRAM)
            lookback[seen_context++] = wln[a];
          else{   
            for(unsigned int i=0;i<NGRAM-1;i++)
              lookback[i] = lookback[i+1]; 
            lookback[NGRAM-1] = wln[a];  
          }
          
          if(!BuildContextTree(root, (const char*)lookback, seen_context)){
            fprintf(stderr,"Error: failure in building n-gram trie\n");
            return false;
          }

        }
        break;
      }
      else
        Cc += fp_prob; 
    }
   
#else
    for(edge=state->transitions;edge;edge=edge->nxt)
      T += edge->c; 
    
    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range); // potential -1 herei
    
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
#endif

    unsigned int new_low = (unsigned int)low + (unsigned int)floor((range*Cc)/T); 
    unsigned int new_high = (unsigned int)low + (unsigned int)floor((range*Cn)/T)-1;  // should there be a minus 1 here?
                                                                          
    low = new_low;
    high = new_high;

    unsigned char lb = low & (1 << 15) ? 1:0;
    unsigned char hb = high & (1 << 15) ? 1:0;
    unsigned char lb2 = low & (1 << 14) ? 1:0;
    unsigned char hb2 = high & (1 << 14) ? 1:0;

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
          if(!fread(&ch,sizeof(unsigned char),1,ifp)){
            ch = UINT8_MAX;
            
            fprintf(stderr,"upper?\n");
            exit(1); 

            if(!decode_zero_added){
              ch ^= (1 << 7);
              decode_zero_added = true;
            }
          }
          enc_pos = 0;
        } 

        lb = low & (1 << 15) ? 1:0;
        hb = high & (1 << 15) ? 1:0;
      }
    }
    else if (lb2 && !hb2){      
      unsigned int p = 0;
      unsigned short int encoded_shift = 0; 
      for(int j=15;j>=0;j--){
        if(j != 14){
          if (encoded & (1 << j))
            encoded_shift ^= (1 << (15-p) );
          p++;
        }
      }
      
      if(ch & (1 << 7))
        encoded_shift ^= 1;

      ch <<= 1;
      enc_pos++;

      if(enc_pos == 8){ // read the next block
        if(!fread(&ch,sizeof(unsigned char),1,ifp)){
          ch = UINT8_MAX;

          fprintf(stderr,"upper?\n"); 
          exit(1); 

          if(!decode_zero_added){
            ch ^= (1 << 7);
            decode_zero_added = true;
            fprintf(stderr,"here?\n");
            exit(1);
          }
        }
        enc_pos = 0;
      } 
      
      encoded = encoded_shift; // set the bit to spliced

      low <<= 1;
      high <<= 1; 

      low ^= (1 << 15);
      high ^= (1 << 15);
      high ^= 1;
    }
  }
  
  RReleaseTree(root); 
  return true; 
}
