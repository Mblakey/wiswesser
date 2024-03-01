#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <set>

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

using namespace OpenBabel; 

double WLNFPTanimoto(u_int8_t *fp1, u_int8_t *fp2){
  unsigned int AnB = 0; 
  unsigned int A = 0; 
  unsigned int B = 0;

    for(unsigned int i=0;i<FPSIZE;i++){
      uint8_t a = fp1[i];
      uint8_t b = fp2[i];
      
      if(a==b)
        AnB++;
      
      A += 1;
      B += 1; 
    }
  
  return (double)(AnB)/(double)(A+B-AnB);
}


double WLNBSTanimoto(u_int8_t *fp1, u_int8_t *fp2){
  unsigned int AnB = 0; 
  unsigned int A = 0; 
  unsigned int B = 0;

  for(unsigned int i=0;i<SCREENSIZE;i++){
    uint16_t a = fp1[i];
    uint16_t b = fp2[i];
    for(int i=7;i>=0;i--){
      if( ( a & (1 << i)) == ( b & (1 << i)) )
        AnB++; 
        
    }
    A += 8;
    B += 8; 
  }

  return (double)(AnB)/(double)(A+B-AnB);
}

double OBabelTanimoto(const char *str1, const char *str2){
  
  std::vector<unsigned int> first_fp;
  std::vector<unsigned int> second_fp;

  OpenBabel::OBFingerprint* fp = OpenBabel::OBFingerprint::FindFingerprint("MACCS");
  
  OBMol mol1;
  OBMol mol2; 
  OBConversion conv; 
  conv.SetInFormat("smi");
    
  conv.ReadString(&mol1, str1);
  conv.ReadString(&mol2, str2);
  
  fp->GetFingerprint(&mol1, first_fp);
  fp->GetFingerprint(&mol2, second_fp);

  double tanimoto = OBFingerprint::Tanimoto(first_fp,second_fp);
  return tanimoto; 
}


double LingoTanimoto(const char *str1, const char *str2){

  std::set<std::string> L1 = WLNLingo(str1,strlen(str1)); 
  std::set<std::string> L2 = WLNLingo(str2,strlen(str2)); 
  
  unsigned int AnB = Intersection(L1, L2); 
  unsigned int AuB = Union(L1, L2); 
  return AnB/(double)AuB; 
}



