

#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <iostream>

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

const char *smiles; 


static void DisplayUsage()
{
  fprintf(stderr, "can_babel <smiles>\n");
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i, j;

  smiles = (const char *)0;

  if (argc < 2)
    DisplayUsage();

  j = 0;
  for (i = 1; i < argc; i++)
  {
    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      fprintf(stderr,"Error: no options for obabel canonical\n");
      DisplayUsage();
    }
    else{
      switch(j){
        case 0:
          smiles = ptr;
          break;
      }
      j++;
    }
  }

  return;
}

int main(int argc, char *argv[]){
  ProcessCommandLine(argc,argv);
  std::string res; 

  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;

  conv.SetInFormat("smi");
  conv.SetOutFormat("can"); 

  conv.ReadString(&mol,smiles);
  res = conv.WriteString(&mol);
  std::cout << res;
  return 0; 
}