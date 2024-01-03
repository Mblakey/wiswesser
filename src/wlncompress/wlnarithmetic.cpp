#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <map>

#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"

unsigned int opt_mode = 0;
unsigned int opt_verbose = false;

const char *input;
const char *trainfile; 



void print_bits(void *ptr, int size, int offset=0) 
{
  if(size>8){
    fprintf(stderr,"Error: print functions maxes at 64 bits\n");
    return;
  }
  long long *ch = (long long*)ptr;
  int size_bits = size * 8;
  for(int i = size_bits-1-offset; i>=0; i--)
    fprintf(stderr,"%lld", (*ch >> i) & 1) ;
  
}

void stream_to_bytes(std::string &stream){
  unsigned int char_pos = 0;
  unsigned char out = 0;
  for(unsigned int i=0;i<stream.size();i++){
    if(stream[i])
      out ^= (1 << 7-char_pos);
  
    if(char_pos == 7){
      char_pos = 0;
      fwrite(&out,sizeof(unsigned char),1,stdout);
      out = 0;
    }
    else
      char_pos++;
  }

  unsigned int zero_added = 0; 
  while(char_pos < 8){
    if(zero_added)
      out ^= (1 << 7-char_pos);
    
    zero_added = 1;
    char_pos++; 
  }
  fwrite(&out,sizeof(unsigned char),1,stdout);

  out = 128; // 0111...
  if(!zero_added)
   fwrite(&out,sizeof(unsigned char),1,stdout);
}


bool encode_file(FILE *ifp, FSMAutomata *wlnmodel){
  
  unsigned int lines = 1;
  unsigned char ch = 0;
  fseek(ifp,0,SEEK_END);
  unsigned int file_len = ftell(ifp); // bytes
  fseek(ifp,0,SEEK_SET);

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  unsigned int low = 0; 
  unsigned int high = UINT32_MAX;
  unsigned int underflow_bits = 0;

  std::string cstream;

  for(unsigned int i=0;i<file_len+1;i++){
    if(i<file_len)
      fread(&ch, sizeof(unsigned char), 1, ifp);
    else 
      ch = 0;
    
    if(ch == '\n')
      lines++;

    unsigned int T = 0;
    for(edge=curr->transitions;edge;edge=edge->nxt)
      T += edge->c;
    
    bool found = 0;
    unsigned int Cc = 0;
    unsigned int Cn = 0;
    for(edge=curr->transitions;edge;edge=edge->nxt){
      if(edge->ch == ch){
        Cn += edge->c;
        found = 1;
        curr = edge->dwn;
        break; 
      }
      else
        Cc += edge->c;
    }
    Cn += Cc; 

    if(!found){
      fprintf(stderr,"Error: invalid wln syntax - please remove line: %d\n",lines);
      return false;
    }

    uint64_t range = ((uint64_t)high+1)-(uint64_t)low;
    uint64_t new_low = (uint64_t)low + (uint64_t)floor((range*Cc)/T); 
    uint64_t new_high = (uint64_t)low + (uint64_t)floor((range*Cn)/T); // potentially off by -1 here

    low = new_low; 
    high = new_high; // truncate

    unsigned char lb = low & (1 << 31) ? 1:0;
    unsigned char hb = high & (1 << 31) ? 1:0;
    unsigned char lb2 = low & (1 << 30) ? 1:0;
    unsigned char hb2 = high & (1 << 30) ? 1:0;
    unsigned char ubit = lb ? 0:1;

    if(lb == hb){
      while(lb == hb){
        cstream += lb; 
  
        low <<= 1; // shift in the zero 
        high <<= 1; // shift in zero then set to 1.
        high ^= 1;
        lb = low & (1 << 31) ? 1:0;
        hb = high & (1 << 31) ? 1:0;

        if(underflow_bits){
          for(unsigned int i=0;i<underflow_bits;i++)
            cstream += ubit;
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
  }

  if(opt_verbose){
    fprintf(stderr,"%d to %d bits: %f compression ratio\n",
            file_len*8,cstream.size(),
            (double)(file_len*8)/(double)cstream.size() );
  }

  stream_to_bytes(cstream);
  return true;
}



bool decode_file(FILE *ifp, FSMAutomata *wlnmodel){
  
  fseek(ifp,0,SEEK_END);
  unsigned int file_len = ftell(ifp); // bytes
  fseek(ifp,0,SEEK_SET);

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  unsigned int low = 0; 
  unsigned int high = UINT32_MAX;
  
  // init encoded to 32 bits
  unsigned int zero_added = 0;
  unsigned int enc_pos = 0;
  unsigned int encoded = 0; 

  unsigned char ch = 0;
  unsigned int i=0;

  while(i < 4){ // read 32 bits
    fread(&ch, sizeof(unsigned char), 1, ifp);
    for(int j = 7;j>=0;j--){
      if(ch & (1 << j))
        encoded ^= (1 << 31-enc_pos);
      
      enc_pos++;  
    }
    i++;
  }

  // saved the next char, will read its MSB to shift in
  fread(&ch, sizeof(unsigned char), 1, ifp);
  i++;

  enc_pos = 0;
  unsigned int safety = 0;
  for(;;){
    // safety++;
    // if(safety == 1000)
    //   return true;

    uint64_t T = 0;
    uint64_t Cc = 0;
    uint64_t Cn = 0;
    for(edge=curr->transitions;edge;edge=edge->nxt)
      T += edge->c;
    
  
    uint64_t range = ((uint64_t)high+1)-(uint64_t)low;
    uint64_t scaled_sym = floor((T*(uint64_t)(encoded-low+1)-1)/range); // potential -1 here


    for(edge=curr->transitions;edge;edge=edge->nxt){
      Cn += edge->c;
      if(scaled_sym >= Cc && scaled_sym < Cn){
        if(!edge->ch)
          return true; 
        else
          fwrite(&edge->ch,sizeof(unsigned char),1,stdout);

        curr = edge->dwn;
        break;
      }
      else
        Cc = Cn;
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
          if(i<file_len)
            fread(&ch, sizeof(unsigned char), 1, ifp);
          else
            ch = UINT8_MAX;
          
          i++;
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
            encoded_shift ^= (1 << 31-p);
          p++;
        }
      }
      
      if(ch & (1 << 7))
        encoded_shift ^= 1;

      ch <<= 1;
      enc_pos++;

      if(enc_pos == 8){ // read the next block
        if(i<file_len)
          fread(&ch, sizeof(unsigned char), 1, ifp);
        else
          ch = UINT8_MAX;

        i++;
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


static void DisplayUsage()
{
  fprintf(stderr, "wlncompress <options> <input> > <out>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -c          compress input\n");
  fprintf(stderr, "  -d          decompress input\n");
  fprintf(stderr, "  -t          add an optional train file for edge frequencies (see wlntrain)\n");
  fprintf(stderr, "  -v          verbose debugging statements on\n");
  exit(1);
}


static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  input = (const char *)0;
  trainfile = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        
        case 'c': 
          opt_mode = 1; 
          break;

        case 'd':
          opt_mode = 2; 
          break;

        case 'v':
          opt_verbose = true;
          break;

        case 't':
          if(i < argc - 1){
            i++;
            trainfile = argv[i];
          }
          else{
            fprintf(stderr,"Error: -t must be followed with a file\n");
            DisplayUsage();
          }
          break;

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
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
    DisplayUsage();
  }

  if(!opt_mode){
    fprintf(stderr,"Error: please choose -c or -d for file\n");
    DisplayUsage();
  }

  return;
}


int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);

  FSMAutomata *wlnmodel = CreateWLNDFA(); // build the model 
  wlnmodel->MakeAccept(wlnmodel->root);
  wlnmodel->AddTransition(wlnmodel->root,wlnmodel->root,'\0');

  // to every accept, add the null character, and the newline character pointing back to the root
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept)
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
  }


  if(!trainfile){
    if(opt_verbose)
      fprintf(stderr,"Warning: using order-0 probabilities\n");

    wlnmodel->AssignEqualProbs();
  }
  else{
    FILE *tfp = 0; 
    tfp = fopen(trainfile,"rb");
    if(!tfp){
      fprintf(stderr,"Error: cannot open train file\n");
      return 1;
    }

    unsigned int i=0;
    unsigned int freq = 0; 
    while(fread(&freq,sizeof(unsigned int),1,tfp))
      wlnmodel->edges[i++]->c = freq;
    
    for(unsigned int i=0;i<wlnmodel->num_edges;i++){
      if(!wlnmodel->edges[i]->c){
        fprintf(stderr,"Warning - null count for edge %d, using 1\n",i);
        wlnmodel->edges[i]->c = 1;
      }
    }

    if(opt_verbose)
      fprintf(stderr,"Warning: train file read, using probabilities\n");

    fclose(tfp);
  }

  FILE *fp = 0; 
  fp = fopen(input,"rb");
  if(fp){
    if(opt_mode == 1)
      encode_file(fp,wlnmodel);
    else if (opt_mode == 2)
      decode_file(fp,wlnmodel);

    fclose(fp);
  }
  else{
    fprintf(stderr,"Error: could not open file at %s\n",input);
    return 1;
  }


  delete wlnmodel;
  return 0;
}