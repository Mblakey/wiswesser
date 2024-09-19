#ifndef WLN_API
#define WLN_API

#include <stdio.h>

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

  int ReadWLNString(const char *ptr)
  {


  }


  int ReadWLNFile(FILE *fp)
  {

  }
  
  int WriteWLNString(FILE *fp)
  {

  }
  
  int WriteWLNFile(FILE *fp)
  {

  }

  int WLNGrepBuffer(const char *ptr)
  {

  }

  int WLNGrepFile(FILE *fp)
  {

  }

};

#endif 
