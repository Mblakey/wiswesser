
#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <stdlib.h>
#include <string> 

#define FPSIZE 42
#define SCREENSIZE 6

typedef struct LingoEntry{
  std::string str; 
  unsigned int n; 
}LingoEntry; 

typedef struct LingoTable {
  LingoEntry **ltable;
  unsigned int size; 
}LingoTable;

u_int8_t* WLNFingerprint(const char* string);
u_int8_t *WLNBitScreen(const char *string);
bool WLNDescriptors(const char *string); 
LingoTable* WLNLingo(const char *str, unsigned int len, unsigned int lingo);

double WLNFPTanimoto(u_int8_t *fp1, u_int8_t *fp2, bool bitscreen); 
double OBabelTanimoto(const char *str1, const char *str2);

#endif 

