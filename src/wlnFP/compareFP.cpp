
#include <stdlib.h>
#include <stdio.h> 
#include <inttypes.h>


#include "fingerprint.h"

bool opt_verbose = false; 
unsigned int opt_mode = 0;
const char *str1;
const char *str2; 

static void DisplayUsage(){
  fprintf(stderr, "wlnfp <string> <string>\n");
  exit(1); 
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i,j;

  str1 = (const char *)0;
  str2 = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){

        case 'h':
          DisplayUsage();
          break;

        case 'v':
          opt_verbose = true;
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

        case 1:
          str2 = ptr; 
          break;
            
        default:
          fprintf(stderr,"Error: n-wise comparisons not currently supported\n");
          exit(1);
      }
    }
  }

  if(!str1||!str2){
    fprintf(stderr,"Error: no inputs given\n");
    DisplayUsage();
  } 

  return;
}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
  
  double obfp = OBabelTanimoto(str1, str2);
  fprintf(stderr,"ObabelFP MACCS: %f\n", obfp);

  u_int8_t *fp1 = WLNFingerprint(str1); 
  u_int8_t *fp2 = WLNFingerprint(str2);
  double wlnfp = WLNFPTanimoto(fp1, fp2); 
  fprintf(stderr,"wlnFP: %f\n", wlnfp);
  free(fp1);
  free(fp2); 
  

  fp1 = WLNBitScreen(str1); 
  fp2 = WLNBitScreen(str2);
  double wlnbs = WLNBSTanimoto(fp1, fp2); 
  fprintf(stderr,"wlnBS: %f\n", wlnbs);
  free(fp1);
  free(fp2); 
  
  double wlnlingo = LingoTanimoto(str1, str2); 
  fprintf(stderr,"WLNlingo: %f\n", wlnlingo);

  return 0; 
}
