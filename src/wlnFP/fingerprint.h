
#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <stdlib.h>

#define FPSIZE 32

u_int8_t* WLNFingerprint(const char* string);
bool WLNDescriptors(const char *string); 

double WLNFPTanimoto(u_int8_t *fp1, u_int8_t *fp2); 
double OBabelTanimoto(const char *str1, const char *str2);

#endif 

