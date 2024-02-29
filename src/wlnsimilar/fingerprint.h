
#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <inttypes.h>

#define VERBOSE 1 
#define FPSIZE 32

uint16_t* WLNFingerprint(const char* string);
bool WLNDescriptors(const char *string); 

double WLNFPTanimoto(uint16_t *fp1, uint16_t *fp2); 


#endif 

