#include <stdlib.h>
#include <stdio.h>

#include <vector>
#include <string>

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

#include "fingerprint.h"
#include "parser.h"


double WLNFPTanimoto(u_int8_t *fp1, u_int8_t *fp2){
  unsigned int AnB = 0; 
  unsigned int A = 0; 
  unsigned int B = 0;
  for(unsigned int i=0;i<FPSIZE;i++){
    uint16_t a = fp1[i];
    uint16_t b = fp2[i];
    for(int i=7;i>=0;i--){
      if( ((a >> i) & 1) == ((b >> i) & 1) )
        AnB++; 
    }
    
    A += 8;
    B += 8; 
  }

  return (double)(AnB)/(double)(A+B-AnB);
}

double OBabelTanimoto(const char *str1, const char *str2){
  
  std::string first_smiles;
  std::string second_smiles; 
  std::vector<unsigned int> first_fp;
  std::vector<unsigned int> second_fp;

  OpenBabel::OBFingerprint* fp = OpenBabel::OBFingerprint::FindFingerprint("MACCS");
  
  OpenBabel::OBMol mol_1;
  OpenBabel::OBMol mol_2; 
  OBConversion conv; 
  conv.SetOutFormat("smi");


  if(!ReadWLN(str1, &mol_1))
    return 0.0;
  if(!ReadWLN(str2, &mol_2))
    return 0.0;

  first_smiles = conv.WriteString(&mol_1);
  second_smiles = conv.WriteString(&mol_2);
  
  fprintf(stderr,"1: %s", first_smiles.c_str());
  fprintf(stderr,"2: %s", second_smiles.c_str());

  fp->GetFingerprint(&mol_1, first_fp);
  fp->GetFingerprint(&mol_2, second_fp);

  double tanimoto = OBFingerprint::Tanimoto(first_fp,second_fp);
  return tanimoto; 
}





