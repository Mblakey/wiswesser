#ifndef PARSER_H
#define PARSER_H

#include <openbabel/mol.h>

typedef unsigned char       u8; 
typedef unsigned short      u16; 
typedef unsigned int        u32; 
typedef unsigned long long  u64; 

bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol);
bool WriteWLN(std::string &buffer,OpenBabel::OBMol* mol, bool modern);
bool CanonicaliseWLN(const char *ptr, OpenBabel::OBMol* mol);


int C_ReadWLN(const char *ptr, OpenBabel::OBMol* mol); 
bool NMReadWLN(const char *ptr, OpenBabel::OBMol* mol);

#endif 
