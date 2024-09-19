#ifndef WLN_API
#define WLN_API

// hate the cpp, but easy to work with python
#include <stdio.h>
#include <string>
#include <vector>

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/obconversion.h>
#include <openbabel/babelconfig.h>

class WLNParser {
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;
  int conv_type; 

  int SetConvMethod(const char *type)
  {

  }

  std::string ReadWLNString(std::string &inp)
  {


  }


  std::vector<std::string> ReadWLNFile(FILE *fp)
  {

  }
  
  std::string WriteWLNString(std::string *inp)
  {

  }
  
  std::vector<std::string> WriteWLNFile(FILE *fp)
  {

  }

};

#endif 
