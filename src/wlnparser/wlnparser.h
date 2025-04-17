#ifndef WLN_PARSER_H
#define WLN_PARSER_H

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <inttypes.h>
#include <stdbool.h>

bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol);
bool WriteWLN(std::string &buffer, OpenBabel::OBMol* mol);
bool CanonicaliseWLN(const char *ptr, OpenBabel::OBMol* mol);

int C_ReadWLN(const char *ptr, OpenBabel::OBMol* mol); 
int C_ReadWLNFile(FILE *fp, OpenBabel::OBMol* mol, OpenBabel::OBConversion *conv); 

int WriteWLNFile(FILE *fp, OpenBabel::OBMol* mol, OpenBabel::OBConversion *conv);
#endif 
