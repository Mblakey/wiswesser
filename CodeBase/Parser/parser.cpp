//
// Created by Michael Blakey on 04/07/2022.
//

#include <iostream>
#include "southampton-wln.h"
#include "wln-nextmove.h"
#include <string>
#include <cstring>
#include <algorithm>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <iterator>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <iterator>

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>


#define BUFFER_SIZE 8*4096
const char *filename; 
const char *wln;




static bool ReadLineFromFile(FILE *fp, char *buffer, unsigned int n){
  char *end = buffer+n;
  char *ptr;
  int ch;

  ptr = buffer;
  do {
    ch = getc_unlocked(fp); // this increments fp
    if (ch == '\n') {
      *ptr++ = '\0';
      *ptr = '\0';
      return true;
    }
    if (ch == '\f') {
      *ptr++ = '\0';
      *ptr = '\0';
      return true;
    }
    if (ch == '\r') {
      *ptr++ = '\0';
      *ptr = '\0';
      ch = getc_unlocked(fp);
      if (ch != '\n') {
        if (ch == -1)
          return false;
        ungetc(ch,fp);
      }
      return true;
    }
    if (ch == -1) {
      *ptr++ = '\0';
      *ptr = '\0';
      return ptr-buffer > 1;
    }
    *ptr++ = ch;
  } while (ptr < end);
  *ptr = 0;
  
  fprintf(stderr, "Warning: line too long!\n");
  return false;
}


static void ProcessCommandLine(int argc, char *argv[]){

  const char *ptr=0; 
  int i,j; 

  filename = (const char*)0;
  wln = (const char*)0;

  

  j=0; 
  for (i=1;i<argc;i++){
    ptr = argv[i];

    if (ptr[0]=='-' && ptr[1]){
      // bracket flag
      if (!strcmp(ptr,"-r") || !strcmp(ptr,"--read")){ 
        // expects a folder following a space
        if (i == argc-1){
          fprintf(stderr,"ERROR: Provide a valid wln after -r\n");
          exit(1);
        }
        else{
          ptr = argv[i+1];
          if (!ptr[0] || ptr[0] == '-'){
            fprintf(stderr,"ERROR: Provide a valid wln after -r\n");
            exit(1);
          }else {
            wln = ptr;
            i++;
          } 
        } 
      }
        
    }
    else switch(j++){
      case 0: filename = ptr; break; 
      default:
        fprintf(stderr,"No options\n");
        exit(1);
    }
  }    
}


// returns the number of matches in a file from a file pointer
unsigned int WLNReadFilePointer(FILE *ifp){
  fprintf(stderr,"matching on disc file\n");
  // arbituary values; 
  char buffer[BUFFER_SIZE];
  unsigned int match_count = 0; 
  // this now processes a file and returns the match_count

  while (ReadLineFromFile(ifp,buffer,BUFFER_SIZE-1)) {
    std::string str = std::string(buffer);
    std::transform(str.begin(), str.end(),str.begin(), ::toupper);
    const char *smiles = WLNToSmiles(str.c_str(),"smi");
    if (strcmp(smiles,"NULL") != 0){
      fprintf(stdout,"%s\t%s\t%d\n",str.c_str(),smiles,strlen(str.c_str()));
      match_count++;
    }
  }

  fprintf(stdout,"Valid WLN: %d\n",match_count);
  return match_count;
}

int main(int argc, char *argv[]){
  ProcessCommandLine(argc,argv);

  if (wln && *wln && strcmp(wln,"-")) {
    const char *smiles = WLNToSmiles(wln,"smi");
    fprintf(stderr,"%s    %s\n",wln,smiles);
    if (smiles){
      delete [] smiles; 
      smiles = 0;
    }
      
    return 0;
  }


  FILE *ifp =0;
  if (filename && *filename && strcmp(filename,"-")) {
    ifp = fopen(filename,"r");
    if (!ifp) {
      fprintf(stderr,"Error: Cannot read input file: %s\n",filename);
      return 1;
    }
    WLNReadFilePointer(ifp);
    return 0;
  } 



  return 0;
}
