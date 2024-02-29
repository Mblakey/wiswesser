

#include <stdlib.h>
#include <stdio.h>

#include "fingerprint.h"

const char *str1; 

static void DisplayUsage(){
  fprintf(stderr, "wlndesc <string>\n");
  exit(1); 
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i,j;

  str1 = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        
        case 'h':
          DisplayUsage();
          break;

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
    else{
      switch(j++){
        case 0:
          str1 = ptr;  
          break;
            
          
        default:
          fprintf(stderr,"Error: descriptor debugging takes in a single arguement\n");
          exit(1);
      }
    }
  }

  if(!str1){
    fprintf(stderr,"Error: no inputs given\n");
    DisplayUsage();
  } 

  return;
}

int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
  WLNDescriptors(str1); 
  return 0; 
}

