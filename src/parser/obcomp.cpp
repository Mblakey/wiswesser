

#include <stdlib.h>
#include <stdio.h>

#include <string>
#include <iostream>
#include <limits>

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
#include <openbabel/fingerprint.h>

#include <openbabel/stereo/stereo.h> // remove stereo

const char *smiles_1; 
const char *smiles_2; 


bool AreSame(double a, double b, const double EPSILON)
{
  return fabs(a - b) < EPSILON;
}


static void DisplayUsage()
{
  fprintf(stderr, "can_babel <smiles> <smiles>\n");
  exit(1);
}


static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i, j;

  smiles_1 = (const char *)0;

  if (argc < 3){
    fprintf(stderr,"Error: not enough args\n");
    DisplayUsage();
  }
    

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
          smiles_1 = ptr;
          break;
        case 1:
          smiles_2 = ptr;
          break;
      }
      j++;
    }
  }

  return;
}

using namespace OpenBabel;

int main(int argc, char *argv[]){
  ProcessCommandLine(argc,argv);
  
  const double epsilon = std::numeric_limits<double>::epsilon(); 

  OpenBabel::OBMol mol_1;
  OpenBabel::OBMol mol_2;

  OpenBabel::OBConversion conv;

  conv.SetInFormat("smi");

  conv.ReadString(&mol_1,smiles_1);
  conv.ReadString(&mol_2,smiles_2);
  
  // removes all the stereo
  mol_1.DeleteData(27);
  mol_2.DeleteData(27);

  std::vector<unsigned int> first_fp;
  std::vector<unsigned int> second_fp;

  OpenBabel::OBFingerprint* fp = OpenBabel::OBFingerprint::FindFingerprint("");

  fp->GetFingerprint(&mol_1, first_fp);
  fp->GetFingerprint(&mol_2, second_fp);

  double tanimoto = OBFingerprint::Tanimoto(first_fp,second_fp);

  // single atom fingerprints cause issues
  if(mol_1.NumAtoms() == 1 && mol_2.NumAtoms() == 1){
    if(mol_1.GetAtom(1)->GetAtomicNum() == mol_1.GetAtom(1)->GetAtomicNum())
      fprintf(stdout,"%d",1);
    else
      fprintf(stdout,"%d",0);

    return 0;
  }
  
  if(AreSame(tanimoto,1.00,epsilon))
    fprintf(stdout,"%d",1);
  else
    fprintf(stdout,"%d",0);

  return 0;
}