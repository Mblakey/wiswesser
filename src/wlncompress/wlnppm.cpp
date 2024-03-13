#include <cstdint>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "rfsm.h"
#include "ctree.h"
#include "wlnzip.h"

#define NGRAM 5
#define PPM 1
#define ALPHABET 42
#define TERMINATE 'x'
#define UPDATE_EXCLUSION 0
#define ASCII_EXCLUDES 1
#define ESCAPE 'C'
/* linked list bitstream, just for single string encoding
 * purposes 
*/

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

/* updates the lookback array, returns the node for the longest 
 * context ready for next iteration */
Node* UpdateCurrentContext(Node *root, unsigned char *lookback, unsigned int seen_context){
  Node * curr_context = root;
  Edge *cedge = 0; 
  for(unsigned int i=1;i<seen_context;i++){
    for(cedge = curr_context->leaves;cedge;cedge=cedge->nxt){
      if(cedge->dwn->ch == lookback[i]){
        curr_context = cedge->dwn;
        break;
      }
    }
  }
  return curr_context;
}


BitStream* WLNPPMCompressBuffer(const char *str, FSMAutomata *wlnmodel){ 
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx"; // TERMINATE = 'x'
  wlnmodel->AssignEqualProbs();
  
  BitStream *bitstream = (BitStream*)malloc(sizeof(BitStream)); 
  bitstream->b = 0; 
  bitstream->nxt = 0; 

  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;

  unsigned int read_bytes = 0; 
  unsigned int out_bits = 0; 

  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX; // set all the bits to 11111...n 
  unsigned int underflow_bits = 0;
  
  unsigned int seen_context = 0;
  unsigned char lookback[NGRAM+1] = {0}; 
  bool ascii_exclude[255] = {false};  // for exclusion mechanism
  unsigned int excluded = 0; 

  unsigned char ch = *str;
  read_bytes++; 

  bool stop = false;
   
  Node *root = AllocateTreeNode('0', 0); // place to return to
  Node *curr_context = 0; 
  Edge *cedge = 0; 
  root->c = 1; 

  // the FSM is a fallback model, if PPM does not work, not intergrated with PPM itself. 
  // our order-1 is already incredibley good

  for(;;){
    unsigned int T = 0;
    unsigned int Cc = 0; 
    unsigned int Cn = 0; 

#if PPM
    bool encoding_escape = false; // this stops ch pointer incrementing
    unsigned int e_o = 1; 
    
    // order- -1 model, escape is not needed here 
    if(!curr_context){
      T = ALPHABET-excluded; 
      for (unsigned int i=0; i<ALPHABET; i++){
        if(!ascii_exclude[wln[i]]){
          if(wln[i] == ch){
            Cn += 1;
            break;
          }
          else
            Cc += 1;
        }
      }
      Cn+= Cc;

      memset(ascii_exclude, 0, 255); // reset the exclusions
      excluded = 0; 
    }
    else{
      bool found = false;
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
        if(!ascii_exclude[cedge->dwn->ch])
          T += cedge->dwn->c;
      }
      
      // methods for escape calculation go here
      T+= e_o; // add the escape frequency in
      
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
        if(!ascii_exclude[cedge->dwn->ch]){
          
          if(cedge->dwn->ch == ch){
            Cn += cedge->dwn->c;
            found = true; 
            break;
          }
          Cc+= cedge->dwn->c;

        }
      }

      Cn+=Cc; 
      
      if(!found){
        encoding_escape = true;
        Cn += e_o; // escape probability to high range
#if ASCII_EXCLUDES
      // exclude the characters
        for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
          if(!ascii_exclude[cedge->dwn->ch]){
            ascii_exclude[cedge->dwn->ch]= true;
            excluded++;
          }
        }
#endif
      }
      
      // encoded char will be something just not in the trie, test binded candidancy
    }
#else
    bool found = 0; 
    for(edge=state->transitions;edge;edge=edge->nxt)
      T += edge->c;
 
    for(edge = state->transitions;edge;edge=edge->nxt){
      if(ch == edge->ch){
        Cn += edge->c;
        state = edge->dwn; // move the state;
        
        edge->c++;
        if(edge->c == 64)
          edge->c >>= 1;
        
        found = 1; 
        break;
      }
      else
        Cc += edge->c;
    }
    Cn += Cc;
    
    if(!found){
      fprintf(stderr,"Error: invalid wln notation - %c\n",ch);
      return (BitStream*)0;
    }
#endif

    // standard arithmetic coder 16 bit int precision. 
    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int new_low = (unsigned int)low + (unsigned int)floor((range*Cc)/T); 
    unsigned int new_high = (unsigned int)low + (unsigned int)floor((range*Cn)/T)-1; // yes!  

    // truncate down
    low = new_low;
    high = new_high;

    for(;;){
      
      unsigned char lb = low & (1 << 15) ? 1:0;
      unsigned char hb = high & (1 << 15) ? 1:0;
      unsigned char lb2 = low & (1 << 14) ? 1:0;
      unsigned char hb2 = high & (1 << 14) ? 1:0;

      if(lb == hb){
        Append(lb,bitstream);
        out_bits++;

        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;

        for(unsigned int i=0;i<underflow_bits;i++){
          Append(!lb, bitstream);  
          out_bits++;
        }
        underflow_bits = 0;
      }
      else if (lb2 && !hb2){      

        high <<= 1; 
        high ^= (1 << 15);
        high ^= 1;

        low <<= 1;
        low &= (1 << 15)-1;
        
        underflow_bits++;
      }
      else 
        break; 
    }


// #################################################################################
// PPM trie update happens after an encoding so the decompressor can keep up
#if PPM
    if(!encoding_escape){
      
      if(seen_context < NGRAM)
        lookback[seen_context++] = ch;
      else{      
        // escaping with no context is pointless
        for(unsigned int i=0;i<NGRAM-1;i++)
          lookback[i] = lookback[i+1]; 
        lookback[NGRAM-1] = ch; 
      }

      BuildContextTree(root, (const char*)lookback, seen_context,UPDATE_EXCLUSION); 
      memset(ascii_exclude,0,255);
      excluded = 0; 
    
      curr_context = UpdateCurrentContext(root,lookback,seen_context);    
    }
    else
      curr_context = curr_context->vine; // move to the lower context
  
// #################################################################################
    if(!encoding_escape){
      if(stop)
        break; 
      else{
        ch = *(++str);  
        read_bytes++; 
      }
    }
#else
    if(stop)
      break;
    else{
      ch = *(++str);  
      read_bytes++; 
    }
#endif

    if(!ch){
      stop = true;
      ch = TERMINATE; // special character for WLN ending
    }
  }
     
//  WriteDotFile(root, stdout); 
  
  Append(0, bitstream);
  Append(1, bitstream);
  out_bits+=2; 

  fprintf(stderr,"%d/%d bits = %f\n",out_bits,read_bytes*8, (read_bytes*8)/(double)out_bits); 
  RReleaseTree(root);
  return bitstream; 
}



bool WLNPPMDecompressBuffer(BitStream *bitstream, FSMAutomata *wlnmodel){
   
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx";
  wlnmodel->AssignEqualProbs(); // you moron

  FSMState *state = wlnmodel->root;
  FSMEdge *edge = 0 ;

  BitStream *bit = bitstream->nxt; // get first bit
   
  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX;
  unsigned short int encoded = UINT16_MAX;
  unsigned int shift_pos = 0; 

  unsigned int seen_context = 0;   
  unsigned char lookback[NGRAM+1] = {0}; 

  bool ascii_exclude[255] = {false};  // for exclusion mechanism
  unsigned int excluded = 0; 

  Node *root = AllocateTreeNode('0', 0); // place to return to
  Node *curr_context = 0; 
  Edge *cedge = 0; 
  root->c = 1; 

  while(bit->nxt && shift_pos < 16){ // read 16 max into encoded
    if(!bit->b)
      encoded ^= (1 << (15-shift_pos) );

    shift_pos++;
    bit = bit->nxt;
  }

  // pre-load next char if relevent
  for(;;){
     
    unsigned int T = 0; 
    unsigned int Cc = 0;
    unsigned int Cn = 0;
    
    unsigned int e_o = 1; 

#if PPM
    if(!curr_context){
      T = ALPHABET-excluded; 
    }
    else{
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
        if(!ascii_exclude[cedge->dwn->ch])
          T+= cedge->dwn->c;
      }
      T+= e_o; 
    }

    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range); 
    
    if(!curr_context){
      for(unsigned int a=0;a<ALPHABET;a++){
        if(!ascii_exclude[wln[a]]){
          Cn += 1; 
          if(scaled_sym >= Cc && scaled_sym < Cn){

            if(wln[a] == TERMINATE)
              return true;
            else{
              fputc(wln[a],stdout);
              if(seen_context < NGRAM)
                lookback[seen_context++] = wln[a];
              else{   
                for(unsigned int i=0;i<NGRAM-1;i++)
                  lookback[i] = lookback[i+1]; 
                lookback[NGRAM-1] = wln[a];  
              }
             
              BuildContextTree(root, (const char*)lookback, seen_context,UPDATE_EXCLUSION);
              memset(ascii_exclude,0,255);
              excluded = 0; 

              curr_context = UpdateCurrentContext(root,lookback,seen_context);    
              break;
            }
          }
          else
            Cc += 1; 
        }
      }
    }
    else{ 
      bool found = false;
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){   
        if(!ascii_exclude[cedge->dwn->ch]){
          Cn += cedge->dwn->c; 
          if(scaled_sym >= Cc && scaled_sym < Cn){   
            found = true;
            fputc(cedge->dwn->ch,stdout);

            if(seen_context < NGRAM)
              lookback[seen_context++] = cedge->dwn->ch;
            else{   
              for(unsigned int i=0;i<NGRAM-1;i++)
                lookback[i] = lookback[i+1]; 
              lookback[NGRAM-1] = cedge->dwn->ch;  
            }
             
            BuildContextTree(root, (const char*)lookback, seen_context,UPDATE_EXCLUSION);
            memset(ascii_exclude,0,255);
            excluded = 0; 
             
            curr_context = UpdateCurrentContext(root,lookback,seen_context);    
            break;
          }
          else 
            Cc += cedge->dwn->c; 
        }
      }
      
      if(!found){

#if ASCII_EXCLUDES
        // must be an escape character, add ascii exclusion and move back a vine
        for(cedge = curr_context->leaves;cedge;cedge=cedge->nxt){
          if(!ascii_exclude[cedge->dwn->ch]){
            ascii_exclude[cedge->dwn->ch] = true;
            excluded++; 
          }
        }
#endif
        curr_context = curr_context->vine;
        Cn += e_o; 
      }
    }
#else
    for(edge=state->transitions;edge;edge=edge->nxt)
      T += edge->c; 
    
    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range); // potential -1 herei
    
    for(edge=state->transitions;edge;edge=edge->nxt){
      Cn += edge->c;
      if(scaled_sym >= Cc && scaled_sym < Cn){  
        if(edge->ch == TERMINATE)
          return true;
        else
          fputc(edge->ch,stdout);
        
        edge->c++;
        if(edge->c == 64)
          edge->c >>= 1;

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

    for(;;){
      unsigned char lb = low & (1 << 15) ? 1:0;
      unsigned char hb = high & (1 << 15) ? 1:0;
      unsigned char lb2 = low & (1 << 14) ? 1:0;
      unsigned char hb2 = high & (1 << 14) ? 1:0;

      if(lb == hb){

        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;

        encoded <<= 1;
        encoded |= bit->b; 
        if(bit->nxt)
          bit = bit->nxt;

        lb = low & (1 << 15) ? 1:0;
        hb = high & (1 << 15) ? 1:0;
      }
      else if (lb2 && !hb2){      
        unsigned short int msb = encoded >> 15;
        unsigned short int rest = encoded & 0x3fff; // bits 0 to 14
        encoded = (msb<<15)|(rest<<1);
        encoded |= bit->b;
        if(bit->nxt)
          bit = bit->nxt;
    

        high <<= 1; 
        high ^= (1 << 15);
        high ^= 1;
        
        low <<= 1;
        low &= (1 << 15) -1;
      }
      else 
        break;
    }
  }

  RReleaseTree(root); 
  return true; 
}



/* used for outputting to the stream */
void append_bit(unsigned char bit,unsigned char *ch){
  bit = bit?1:0;
  (*ch) <<= 1;
  (*ch) |= bit; 
}

void transfer_bit(unsigned char *ch, unsigned short int *encoded){
  (*encoded) |= ((*ch) >> 7);
  (*ch) <<= 1; 
}

void print_bits(unsigned char *ch){
  for(int i=7;i>=0;i--){
    if( (*ch) & (1<<i) )
      fputc('1',stderr);
    else
      fputc('0',stderr);
  }
}



bool WLNPPMCompressFile(FILE *ifp, FSMAutomata *wlnmodel){
    
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx"; // TERMINATE = 'x', EXCLUSION = 'e'
  wlnmodel->AssignEqualProbs(); // you moron
   
  FSMState *state = wlnmodel->root; 
  FSMEdge *edge = 0;
  
  unsigned int read_bytes = 0; 
  unsigned int out_bits = 0; 

  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX; // set all the bits to 11111...n 
  unsigned int underflow_bits = 0;

  unsigned int seen_context = 0;
  unsigned char lookback[NGRAM+1] = {0}; 
  bool ascii_exclude[256] = {false};  // for exclusion mechanism
  unsigned int excluded = 0; 

  bool stop = false;
   
  Node *root = AllocateTreeNode('.', 0); // place to return to
  Node *curr_context = 0; 
  Edge *cedge = 0; 
  root->c = 1; 
  
  unsigned int stream_pos = 0; 
  unsigned char stream = 0;
  
  unsigned char ch = 0; 
  read_bytes++;
  if(!fread(&ch,sizeof(unsigned char),1,ifp)){
    fprintf(stderr,"Error: no data in file\n"); 
    return false;
  }
  

  for(;;){
    
    bool found = 0; 
    unsigned int T = 0;
    unsigned int Cc = 0; 
    unsigned int Cn = 0; 
    unsigned int e_o = 1; 

#if PPM
    bool encoding_escape = 0; // this stops ch pointer incrementing
    // order- -1 model, escape is not needed here 
    if(!curr_context){
      T = ALPHABET-excluded; 
      for (unsigned int i=0; i<ALPHABET; i++){
        if(!ascii_exclude[wln[i]]){
          if(wln[i] == ch){
            Cn += 1;
            break;
          }
          else
            Cc += 1;
        }
      }
      Cn+= Cc;

      memset(ascii_exclude, 0, 255); // reset the exclusions
      excluded = 0; 
    }
    else{
      bool found = false;
      unsigned int scontexts = 0; 
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
        if(!ascii_exclude[cedge->dwn->ch]){
          T += cedge->dwn->c;
          scontexts++; 
        }
      }
       
      if(ESCAPE == 'C')
        e_o = scontexts ? scontexts:1;
      // methods for escape calculation go here

      T+= e_o; // add the escape frequency in
      
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
        if(!ascii_exclude[cedge->dwn->ch]){
          if(cedge->dwn->ch == ch){
            Cn += cedge->dwn->c;
            found = true; 
            break;
          }
          Cc+= cedge->dwn->c;
        }
      }
      Cn+=Cc; 
      
      if(!found){
        encoding_escape = true;
        Cn += e_o; // escape probability to high range
              
#if ASCII_EXCLUDES
      // exclude the characters
        for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
          if(!ascii_exclude[cedge->dwn->ch]){
            ascii_exclude[cedge->dwn->ch]= true;
            excluded++;
          }
        }
#endif
      }
    }
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

    for(;;){
      
      unsigned char lb = low & (1 << 15) ? 1:0;
      unsigned char hb = high & (1 << 15) ? 1:0;
      unsigned char lb2 = low & (1 << 14) ? 1:0;
      unsigned char hb2 = high & (1 << 14) ? 1:0;

      if(lb == hb){
        append_bit(lb, &stream);
        stream_pos++;
        if(stream_pos == 8){
          fputc(stream,stdout);
          stream = 0;
          stream_pos = 0;
        }
        out_bits++; 

        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;

        for(unsigned int i=0;i<underflow_bits;i++){
          append_bit(!lb, &stream);          
          stream_pos++;
          if(stream_pos == 8){
            fputc(stream,stdout);
            stream = 0; 
            stream_pos = 0;
          }
          out_bits++; 
        }

        underflow_bits = 0;
      }    
      else if (lb2 && !hb2){      
        
        high <<= 1; 
        high |= (1 << 15);
        high |= 1;
        
        low <<= 1;
        low &= (1 << 15)-1;

        underflow_bits++;
      }
      else 
        break;
    }
      
// #################################################################################
#if PPM
    if(!encoding_escape){

      if(seen_context < NGRAM)
        lookback[seen_context++] = ch;
      else{      
        for(unsigned int i=0;i<NGRAM-1;i++)
          lookback[i] = lookback[i+1]; 
        lookback[NGRAM-1] = ch; 
      }

      BuildContextTree(root, (const char*)lookback, seen_context,UPDATE_EXCLUSION); 
      memset(ascii_exclude,0,255);
      excluded = 0; 
    
      curr_context = UpdateCurrentContext(root,lookback,seen_context);    
    }
    else
      curr_context = curr_context->vine; // move to the lower context
  
// #################################################################################
    if(!encoding_escape){
      if(stop)
        break; 
      else{
        if(!fread(&ch,sizeof(unsigned char),1,ifp)){
          ch = TERMINATE;
          stop = true;
        }
        else
          read_bytes++; 
      }
    }
#else
    if(stop)
      break;
    else{
      if(!fread(&ch,sizeof(unsigned char),1,ifp)){
        ch = TERMINATE; 
        stop = true;
      }
      else
        read_bytes++; 
    }
#endif
  }

  append_bit(0, &stream); 
  out_bits++; 
  stream_pos++;
  if(stream_pos == 8){
    fputc(stream,stdout);
    stream_pos = 0;
    stream = 0;
  }

  while(stream_pos < 8){
    append_bit(1, &stream);
    stream_pos++;
    out_bits++;
  }
  fputc(stream,stdout);
  
  fprintf(stderr,"%d/%d bits = %f\n",out_bits,read_bytes*8, (read_bytes*8)/(double)out_bits); 
  RReleaseTree(root);  
  return true;
}



bool WLNPPMDecompressFile(FILE *ifp, FSMAutomata *wlnmodel){

  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789&/- \nx";
  wlnmodel->AssignEqualProbs(); // you moron

  FSMState *state = wlnmodel->root;
  FSMEdge *edge = 0 ;
  
  unsigned short int low = 0; 
  unsigned short int high = UINT16_MAX;
  unsigned short int encoded = UINT16_MAX;

  unsigned int shift_pos = 0; 
  unsigned char ch = 0;
  
  unsigned int seen_context = 0;   
  unsigned char lookback[NGRAM+1] = {0}; 
  bool ascii_exclude[255] = {false};  // for exclusion mechanism
  unsigned int excluded = 0; 

  Node *root = AllocateTreeNode('0', 0); // place to return to
  Node *curr_context = 0; 
  Edge *cedge = 0; 
  root->c = 1; 

  for(unsigned int i=0;i<2;i++){ // read 16 max into encoded
    if(!fread(&ch,sizeof(unsigned char),1,ifp))
      ch = UINT8_MAX;
    
    for(int j=7;j>=0;j--){
      if( (ch & (1 << j))==0 )  
        encoded ^= (1 << (15-shift_pos) );      
      
      shift_pos++;  
    }
  }
   
  // pre-load next char, ready for transfer
  shift_pos = 0;
  if(!fread(&ch,sizeof(unsigned char),1,ifp))
    ch = UINT8_MAX; 

  for(;;){
    
    unsigned int T = 0; 
    unsigned int Cc = 0;
    unsigned int Cn = 0;
    unsigned int e_o = 1; 

#if PPM
    if(!curr_context){
      T = ALPHABET-excluded; 
    }
    else{
      unsigned int scontexts = 0;
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){
        if(!ascii_exclude[cedge->dwn->ch]){
          T+= cedge->dwn->c;
          scontexts++; 
        }
      }

      if(ESCAPE == 'C')
        e_o = scontexts ? scontexts:1;
      T+= e_o; 
    }

    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range); 
    
    if(!curr_context){
      for(unsigned int a=0;a<ALPHABET;a++){
        if(!ascii_exclude[wln[a]]){
          Cn += 1; 
          if(scaled_sym >= Cc && scaled_sym < Cn){
            if(wln[a] == TERMINATE)
              return true;
            else{
              fputc(wln[a],stdout);
              if(seen_context < NGRAM)
                lookback[seen_context++] = wln[a];
              else{   
                for(unsigned int i=0;i<NGRAM-1;i++)
                  lookback[i] = lookback[i+1]; 
                lookback[NGRAM-1] = wln[a];  
              }
              BuildContextTree(root, (const char*)lookback, seen_context,UPDATE_EXCLUSION);
              memset(ascii_exclude,0,255);
              excluded = 0; 

              curr_context = UpdateCurrentContext(root,lookback,seen_context);    
              break;
            }
          }
          else
            Cc += 1; 
        }
      }
    }
    else{ 
      bool found = false;
      for(cedge=curr_context->leaves;cedge;cedge=cedge->nxt){   
        if(!ascii_exclude[cedge->dwn->ch]){
          Cn += cedge->dwn->c; 
          if(scaled_sym >= Cc && scaled_sym < Cn){   
            found = true;
            fputc(cedge->dwn->ch,stdout);
            if(seen_context < NGRAM)
              lookback[seen_context++] = cedge->dwn->ch;
            else{   
              for(unsigned int i=0;i<NGRAM-1;i++)
                lookback[i] = lookback[i+1]; 
              lookback[NGRAM-1] = cedge->dwn->ch;  
            }
             
            BuildContextTree(root, (const char*)lookback, seen_context,UPDATE_EXCLUSION);
            memset(ascii_exclude,0,255);
            excluded = 0;

            curr_context = UpdateCurrentContext(root,lookback,seen_context);    
            break;
          }
          else 
            Cc += cedge->dwn->c; 
        }
      }
      
      if(!found){
        // must be an escape character, add ascii exclusion and move back a vine
#if ASCII_EXCLUDES
        for(cedge = curr_context->leaves;cedge;cedge=cedge->nxt){
          if(!ascii_exclude[cedge->dwn->ch]){
            ascii_exclude[cedge->dwn->ch] = true;
            excluded++;
          }
        }
#endif
        curr_context = curr_context->vine;
        Cn += e_o; 
      }
    }
#else
    for(edge=state->transitions;edge;edge=edge->nxt)
      T += edge->c; 
    
    unsigned int range = ((unsigned int)high+1)-(unsigned int)low;
    unsigned int scaled_sym = floor((T*(unsigned int)(encoded-low+1)-1)/range); // potential -1 herei
    
    for(edge=state->transitions;edge;edge=edge->nxt){
      Cn += edge->c;
      if(scaled_sym >= Cc && scaled_sym < Cn){  
        if(edge->ch == TERMINATE)
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
    unsigned int new_high = (unsigned int)low + (unsigned int)floor((range*Cn)/T)-1;  
                                                                          
    low = new_low;
    high = new_high;

    for(;;){
      
      unsigned char lb = low & (1 << 15) ? 1:0;
      unsigned char hb = high & (1 << 15) ? 1:0;
      unsigned char lb2 = low & (1 << 14) ? 1:0;
      unsigned char hb2 = high & (1 << 14) ? 1:0;

      if(lb == hb){
 
        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;

        encoded <<= 1;
        transfer_bit(&ch, &encoded);
        shift_pos++;
        if(shift_pos == 8){
          if(!fread(&ch,sizeof(unsigned char),1,ifp))
            ch = UINT8_MAX;
          shift_pos = 0;
        }

        // move a bit from ch to encoded.
      }
      else if (lb2 && !hb2){      
        unsigned short int msb = encoded >> 15;
        unsigned short int rest = encoded & 0x3fff; // bits 0 to 14
        encoded = (msb<<15)|(rest<<1); // remember to OR the bit with char. 
        transfer_bit(&ch, &encoded); 
        shift_pos++; 
        if(shift_pos == 8){
          if(!fread(&ch,sizeof(unsigned char),1,ifp))
            ch = UINT8_MAX; // 01 has been assumed to be encoded from compressor
          shift_pos = 0; 
        }

        high <<= 1;
        high |= (1 << 15);
        high |= 1;

        low <<= 1;
        low &= (1 << 15)-1;
      }
      else
        break;
    }
  }
 
  RReleaseTree(root); 
  return true; 
}
