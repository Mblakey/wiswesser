
#ifndef LZ_H
#define LZ_H

#define LZBUCKETS 30
#define WINDOW 290
#define BACKREFERENCE 32768
#define BUFFSIZE WINDOW+BACKREFERENCE

typedef struct{
  unsigned char symbol;
  unsigned int lstart;
  unsigned int dstart; 
  unsigned int lbits;
  unsigned int dbits;
} LLBucket; 

LLBucket ** init_buckets();
void free_buckets(LLBucket **buckets);
LLBucket *length_bucket(unsigned int length, LLBucket **buckets);
LLBucket *distance_bucket(unsigned int distance, LLBucket **buckets);

#endif
