#ifndef PARSER_H
#define PARSER_H

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/kekulize.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/graphsym.h>

using namespace OpenBabel; 

bool ReadWLN(const char *ptr, OBMol* mol);
bool WriteWLN(std::string &buffer, OBMol* mol);
bool WriteModernWLN(std::string &buffer, OBMol* mol);

#endif 