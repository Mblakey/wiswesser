#ifndef WLN_PARSER_H
#define WLN_PARSER_H

#include <inttypes.h>
#include <stdbool.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

// bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol);
// bool WriteWLN(std::string &buffer, OpenBabel::OBMol* mol);
// bool CanonicaliseWLN(const char *ptr, OpenBabel::OBMol* mol);

bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol); 
bool ReadWLNFile(FILE *fp, OpenBabel::OBMol* mol, OpenBabel::OBConversion *conv); 

int WriteWLNFile(FILE *fp, OpenBabel::OBMol* mol, OpenBabel::OBConversion *conv);
#endif 
