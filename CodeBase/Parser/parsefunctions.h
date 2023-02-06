#ifndef READERLIB_WLN_NEXTMOVE_H
#define READERLIB_WLN_NEXTMOVE_H

//
// Created by Michael Blakey on 08/02/2022.
//

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/kekulize.h>
#include <openbabel/obconversion.h>

#include <vector>
#include <numeric>

bool NMReadWLN(const char *ptr, OpenBabel::OBMol* mol);
bool MBWLNMurckoScaffold(const char *ptr, OpenBabel::OBMol* mol);
bool MBWLNRGroupDecomp(const char *ptr, OpenBabel::OBMol* mol);
bool MBWLNGraphScaffold(const char *ptr, OpenBabel::OBMol* mol);


// conversion functions
void ConvertWLN(const char *input_string, const char *format, OpenBabel::OBMol *Mol);
const char* ConvertSMI(const char *input_string, const char *format); 


const char*  MurckoScaffold(const char *test_string, const char *format);
const char*  GraphScaffold(const char *test_string, const char *format);
const char*  RGroupDecomp(const char *test_string, const char *format);


bool MBWriterWLN(OpenBabel::OBMol *mol, std::string &buffer);

#endif //READERLIB_WLN_NEXTMOVE_H