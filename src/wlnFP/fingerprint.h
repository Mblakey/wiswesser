
#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <stdlib.h>

#define FPSIZE 42
#define SCREENSIZE 6

u_int8_t* WLNFingerprint(const char* string);
u_int8_t *WLNBitScreen(const char *string);
bool WLNDescriptors(const char *string); 

double WLNFPTanimoto(u_int8_t *fp1, u_int8_t *fp2, bool bitscreen); 
double OBabelTanimoto(const char *str1, const char *str2);

#endif 

