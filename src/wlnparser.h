#ifndef OBABEL_WLN_PARSER_H
#define OBABEL_WLN_PARSER_H

#include <inttypes.h>
#include <stdbool.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

bool ReadWLN(const char *buffer, OpenBabel::OBMol* mol); 
bool WriteWLN(char *buffer, unsigned int nbuffer, OpenBabel::OBMol* mol); 

#endif 
