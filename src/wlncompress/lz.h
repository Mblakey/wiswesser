
#ifndef LZ_H
#define LZ_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define LZBUCKETS 30

/* ######################################################################################### */

// again not the most efficient but means we wont get lost
typedef struct{
  unsigned int lstart;
  unsigned int dstart; 

  unsigned int lbits;
  unsigned int dbits;
} LLBucket; 



/* we can keep the original DEFLATE specification for distances, but for lengths we choose our own */
LLBucket ** init_buckets(){
  LLBucket **buckets = (LLBucket**)malloc(sizeof(LLBucket*)*LZBUCKETS); // use for instant look up
  memset(buckets,0,sizeof(LLBucket*));

  for(unsigned int i=0;i<LZBUCKETS;i++)
    buckets[i] = (LLBucket*)malloc(sizeof(LLBucket));

  buckets[0]->dstart = 1;
  buckets[0]->dbits = 0;

  buckets[1]->dstart = 2;
  buckets[1]->dbits = 0;

  buckets[2]->dstart = 3;
  buckets[2]->dbits = 0;

  buckets[3]->dstart = 4;
  buckets[3]->dbits = 0;

  buckets[4]->dstart = 5;
  buckets[4]->dbits = 1;

  buckets[5]->dstart = 7;
  buckets[5]->dbits = 1;

  buckets[6]->dstart = 9;
  buckets[6]->dbits = 2;
  
  buckets[7]->dstart = 13;
  buckets[7]->dbits = 2;

  buckets[8]->dstart = 17;
  buckets[8]->dbits = 3;

  buckets[9]->dstart = 25;
  buckets[9]->dbits = 3;

  buckets[10]->dstart = 33;
  buckets[10]->dbits = 4;
  
  buckets[11]->dstart = 49;
  buckets[11]->dbits = 4;

  buckets[12]->dstart = 65;
  buckets[12]->dbits = 5;

  buckets[13]->dstart = 97;
  buckets[13]->dbits = 5;
  
  buckets[14]->dstart = 129;
  buckets[14]->dbits = 6;

  buckets[15]->dstart = 193;
  buckets[15]->dbits = 6;

  
  buckets[16]->dstart = 257;
  buckets[17]->dstart = 385;
  buckets[18]->dstart = 513;
  buckets[19]->dstart = 769;
  buckets[20]->dstart = 1025;
  buckets[21]->dstart = 1537;
  buckets[22]->dstart = 2049;
  buckets[23]->dstart = 3073;
  buckets[24]->dstart = 4097;
  buckets[25]->dstart = 6145;
  buckets[26]->dstart = 8193;
  buckets[27]->dstart = 12289;
  buckets[28]->dstart = 16385;
  buckets[29]->dstart = 24577;

  return buckets;
}

void PurgeBuckets(LLBucket **buckets){
  for(unsigned int i=0;i<LZBUCKETS;i++){
    if(buckets[i])
      free(buckets[i]);
  }
  free(buckets);
}



#endif