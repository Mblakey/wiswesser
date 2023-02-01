//
// Created by Michael Blakey on 12/07/2022.
//

#ifndef WISSWESSER_SOUTHAMPTON_WLN_H
#define WISSWESSER_SOUTHAMPTON_WLN_H

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <iterator>

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

bool MBWriterWLN(OpenBabel::OBMol *mol, std::string &buffer);

#endif //WISSWESSER_SOUTHAMPTON_WLN_H
