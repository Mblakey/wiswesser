#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "fingerprint.h"

double WLNFPTanimoto(uint16_t *fp1, uint16_t *fp2){
  unsigned int AnB = 0; 
  unsigned int A = 0; 
  unsigned int B = 0;
  for(unsigned int i=0;i<FPSIZE;i++){
    uint16_t a = fp1[i];
    uint16_t b = fp2[i];
    for(int i=15;i>=0;i--){
      if( ((a >> i) & 1) == ((b >> i) & 1) )
        AnB++; 
    }
    
    A += 16;
    B += 16; 
  }

  if(VERBOSE){
    fprintf(stderr,"\nA|B = %d\nAnB = %d\n", A,AnB); 
  }

  return (double)(AnB)/(double)(A+B-AnB);
}
