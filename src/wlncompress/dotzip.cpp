#include <stdlib.h>
#include <stdio.h> 

#include "readdot.h"
#include "rfsm.h"
#include "wlndfa.h"
#include "wlnzip.h"

const char *input;
const char *dotfile;
unsigned int mode = 0; 

#define DEFLATE 1 

static void DisplayUsage()
{
  fprintf(stderr, "dotzip <options> <input> <savefile> > <out>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -c   compress input\n");
  fprintf(stderr, "  -d   decompress input\n");
  fprintf(stderr, "  -s   string input compress (debugging)\n"); 
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  input = (const char *)0;
  dotfile = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        case 'c':
          mode = 1;
          break;
        case 'd':
          mode = 2;
          break;
        case 's':
          mode = 3; 
          break; 

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
        case 1:
          dotfile = ptr; 
          break;
        default:
          fprintf(stderr,"Error: multiple files not currently supported\n");
          exit(1);
      }
    }
  }
  
  if(!input || !dotfile){
    fprintf(stderr,"Error: no input file given\n");
    DisplayUsage();
  }

  if(!mode){
    fprintf(stderr,"Error: select compress/decompress mode\n");
    DisplayUsage(); 
  }

  return;
}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
  
  FILE *fp = 0; 
  FSMAutomata *dotmodel = FSMFromDotFile(dotfile); 

  for(unsigned int i=0;i<dotmodel->num_states;i++){
    if(dotmodel->states[i]->accept){
      dotmodel->AddTransition(dotmodel->states[i],dotmodel->root,'\n');
      dotmodel->AddTransition(dotmodel->states[i],dotmodel->root,127);
    }
  }
  
  dotmodel->AddTransition(dotmodel->root,dotmodel->root,127);  // x is a chosen terminal here, can change

#if DEFLATE // experimental, for comparison only
  if(mode == 1){
    // ppm compress file
    fp = fopen(input, "rb"); 
    if(!fp){
      fprintf(stderr,"Error: could not open file\n"); 
      return 1;
    }
    
    if(!WLNdeflate(fp, dotmodel)){
      fprintf(stderr,"Error: failed to compress file\n"); 
      return 1;
    }

    fclose(fp); 
  }
  else if(mode == 2){
    // ppm decompress file
    fp = fopen(input, "rb"); 
    if(!fp){
      fprintf(stderr,"Error: could not open file\n"); 
      return 1;
    }

    if(!WLNinflate(fp, dotmodel)){
      fprintf(stderr,"Error: failed to compress file\n"); 
      return 1;
    }

    fclose(fp); 
  }
  else{
    fprintf(stderr,"NOP: string deflate not avaliable\n"); 
    return 1; 
  }
#else
  if(mode == 1){
    // ppm compress file
    fp = fopen(input, "rb"); 
    if(!fp){
      fprintf(stderr,"Error: could not open file\n"); 
      return 1;
    }
    
    if(!WLNPPMCompressFile(fp, dotmodel)){
      fprintf(stderr,"Error: failed to compress file\n"); 
      return 1;
    }

    fclose(fp); 
  }
  else if(mode == 2){
    // ppm decompress file
    fp = fopen(input, "rb"); 
    if(!fp){
      fprintf(stderr,"Error: could not open file\n"); 
      return 1;
    }

    if(!WLNPPMDecompressFile(fp, dotmodel)){
      fprintf(stderr,"Error: failed to compress file\n"); 
      return 1;
    }

    fclose(fp); 
  }
  else if (mode == 3){
     
    BitStream *bitstream = WLNPPMCompressBuffer(input, dotmodel);
    if(!bitstream)
      return 1; 
   
    if(!WLNPPMDecompressBuffer(bitstream, dotmodel))
      return 1; 
    
    fprintf(stdout,"\n"); 
    DeleteStream(bitstream); 
  }
#endif

  delete dotmodel;
  return 0;
}
