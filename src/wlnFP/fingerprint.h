
#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <stdlib.h>
#include <string> 
#include <set>

#define FPSIZE 42
#define SCREENSIZE 6
#define LINGO 3

u_int8_t* WLNFingerprint(const char* string);
u_int8_t *WLNBitScreen(const char *string);
bool WLNDescriptors(const char *string); 

std::set<std::string> WLNLingo(const char *str, unsigned int len);

unsigned int Intersection(std::set<std::string> &v1, std::set<std::string> &v2);
unsigned int Union(std::set<std::string> &v1, std::set<std::string> &v2);

double WLNFPTanimoto(u_int8_t *fp1, u_int8_t *fp2); 
double WLNBSTanimoto(u_int8_t *fp1, u_int8_t *fp2); 
double OBabelTanimoto(const char *str1, const char *str2);
double LingoTanimoto(const char *str1, const char *str2);

#endif 

