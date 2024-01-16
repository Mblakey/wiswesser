
#ifndef LZ_H
#define LZ_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define LZBUCKETS 30
#define WINDOW 290
#define BACKREFERENCE 32768
#define BUFFSIZE WINDOW+BACKREFERENCE


/* 
we can keep the original DEFLATE specification for distances, but for lengths we,
keep the original characters and compact. 
*/

typedef struct{
  unsigned char symbol;
  unsigned int lstart;
  unsigned int dstart; 

  unsigned int lbits;
  unsigned int dbits;
} LLBucket; 


LLBucket ** init_buckets(){
  LLBucket **buckets = (LLBucket**)malloc(sizeof(LLBucket*)*LZBUCKETS); // use for instant look up
  memset(buckets,0,sizeof(LLBucket*));

  for(unsigned int i=0;i<LZBUCKETS;i++){
    buckets[i] = (LLBucket*)malloc(sizeof(LLBucket));
    buckets[i]->symbol = 'a' + i;
  }
  
  buckets[0]->dstart = 1;
  buckets[0]->dbits = 0;
  buckets[0]->lstart = 3;
  buckets[0]->lbits = 0;

  buckets[1]->dstart = 2;
  buckets[1]->dbits = 0;
  buckets[1]->lstart = 4;
  buckets[1]->lbits = 0;

  buckets[2]->dstart = 3;
  buckets[2]->dbits = 0;
  buckets[2]->lstart = 5;
  buckets[2]->lbits = 0;

  buckets[3]->dstart = 4;
  buckets[3]->dbits = 0;
  buckets[3]->lstart = 6;
  buckets[3]->lbits = 0;

  buckets[4]->dstart = 5;
  buckets[4]->dbits = 1;
  buckets[4]->lstart = 7;
  buckets[4]->lbits = 0;

  buckets[5]->dstart = 7;
  buckets[5]->dbits = 1;
  buckets[5]->lstart = 8;
  buckets[5]->lbits = 0;

  buckets[6]->dstart = 9;
  buckets[6]->dbits = 2;
  buckets[6]->lstart = 9;
  buckets[6]->lbits = 0;
  
  buckets[7]->dstart = 13;
  buckets[7]->dbits = 2;
  buckets[7]->lstart = 10;
  buckets[7]->lbits = 0;

  buckets[8]->dstart = 17;
  buckets[8]->dbits = 3;
  buckets[8]->lstart = 11;
  buckets[8]->lbits = 1;

  buckets[9]->dstart = 25;
  buckets[9]->dbits = 3;
  buckets[9]->lstart = 13;
  buckets[9]->lbits = 1;

  buckets[10]->dstart = 33;
  buckets[10]->dbits = 4;
  buckets[10]->lstart = 15;
  buckets[10]->lbits = 1;
  
  buckets[11]->dstart = 49;
  buckets[11]->dbits = 4;
  buckets[11]->lstart = 17;
  buckets[11]->lbits = 1;

  buckets[12]->dstart = 65;
  buckets[12]->dbits = 5;
  buckets[12]->lstart = 19;
  buckets[12]->lbits = 2;

  buckets[13]->dstart = 97;
  buckets[13]->dbits = 5;
  buckets[13]->lstart = 23;
  buckets[13]->lbits = 2;
  
  buckets[14]->dstart = 129;
  buckets[14]->dbits = 6;
  buckets[14]->lstart = 27;
  buckets[14]->lbits = 2;

  buckets[15]->dstart = 193;
  buckets[15]->dbits = 6;
  buckets[15]->lstart = 31;
  buckets[15]->lbits = 2;

  buckets[16]->dstart = 257;
  buckets[16]->dbits = 7;
  buckets[16]->lstart = 35;
  buckets[16]->lbits = 3;

  buckets[17]->dstart = 385;
  buckets[17]->dbits = 7;
  buckets[17]->lstart = 43;
  buckets[17]->lbits = 3;

  buckets[18]->dstart = 513;
  buckets[18]->dbits = 8;
  buckets[18]->lstart = 51;
  buckets[18]->lbits = 3;

  buckets[19]->dstart = 769;
  buckets[19]->dbits = 8;
  buckets[19]->lstart = 59;
  buckets[19]->lbits = 3;

  buckets[20]->dstart = 1025;
  buckets[20]->dbits = 9;
  buckets[20]->lstart = 67;
  buckets[20]->lbits = 4;

  buckets[21]->dstart = 1537;
  buckets[21]->dbits = 9;
  buckets[21]->lstart = 83;
  buckets[21]->lbits = 4;

  buckets[22]->dstart = 2049;
  buckets[22]->dbits = 10;
  buckets[22]->lstart = 99;
  buckets[22]->lbits = 4;

  buckets[23]->dstart = 3073;
  buckets[23]->dbits = 10;
  buckets[23]->lstart = 115;
  buckets[23]->lbits = 4;

  buckets[24]->dstart = 4097;
  buckets[24]->dbits = 11;
  buckets[24]->lstart = 131;
  buckets[24]->lbits = 5;

  buckets[25]->dstart = 6145;
  buckets[25]->dbits = 11;
  buckets[25]->lstart = 163;
  buckets[25]->lbits = 5;

  buckets[26]->dstart = 8193;
  buckets[26]->dbits = 12;
  buckets[26]->lstart = 195;
  buckets[26]->lbits = 5;

  buckets[27]->dstart = 12289;
  buckets[27]->dbits = 12;
  buckets[27]->lstart = 227;
  buckets[27]->lbits = 5;

  buckets[28]->dstart = 16385;
  buckets[28]->dbits = 13;
  buckets[28]->lstart = 259;
  buckets[28]->lbits = 5;

  buckets[29]->dstart = 24577;
  buckets[29]->dbits = 13;
  buckets[29]->lstart = 290;
  buckets[29]->lbits = 0;

  return buckets;
}

void free_buckets(LLBucket **buckets){
  for(unsigned int i=0;i<LZBUCKETS;i++)
    free(buckets[i]);
  free(buckets);
}

LLBucket *length_bucket(unsigned int length, LLBucket **buckets){
  LLBucket *lb = 0;
  for(unsigned int i=0;i<LZBUCKETS-1;i++){
    if(!lb && (length >= buckets[i]->lstart && length < buckets[i+1]->lstart)  )
      lb = buckets[i];
  }

  if(!lb)
    lb = buckets[LZBUCKETS-1];

  return lb;
}

LLBucket *distance_bucket(unsigned int distance, LLBucket **buckets){
  LLBucket *db = 0;
  for(unsigned int i=0;i<LZBUCKETS-1;i++){
    if(!db && (distance >= buckets[i]->dstart && distance < buckets[i+1]->dstart))
      db = buckets[i];
  }

  if(!db)
    db = buckets[LZBUCKETS-1];

  return db;
}


#endif