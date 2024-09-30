#ifndef PARSER_H
#define PARSER_H

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

typedef unsigned char       u8; 
typedef unsigned short      u16; 
typedef unsigned int        u32; 
typedef unsigned long long  u64; 

bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol);
bool WriteWLN(std::string &buffer, OpenBabel::OBMol* mol);
bool CanonicaliseWLN(const char *ptr, OpenBabel::OBMol* mol);


int C_ReadWLN(const char *ptr, OpenBabel::OBMol* mol); 
int C_ReadWLNFile(FILE *fp, OpenBabel::OBMol* mol, OpenBabel::OBConversion *conv); 

int WriteWLNFile(FILE *fp, OpenBabel::OBMol* mol, OpenBabel::OBConversion *conv);
#endif 
