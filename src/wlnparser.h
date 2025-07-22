#ifndef OBABEL_WLN_PARSER_H
#define OBABEL_WLN_PARSER_H

#include <inttypes.h>
#include <stdbool.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

bool ReadWLN(const char *ptr, OpenBabel::OBMol* mol); 
bool WriteWLN(char *ptr, OpenBabel::OBMol* mol); 

#endif 
