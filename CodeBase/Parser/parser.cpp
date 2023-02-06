//
// Created by Michael Blakey on 04/07/2022.
//

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <string>

#include "parsefunctions.h"

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>


#define BUFFER_SIZE 8*4096
const char *filename; 
const char *wln;

const char *inpformat; 
const char *outformat; 


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


static bool ReadFormat(const char *ptr){

  fprintf(stderr,"format being read is: %s\n",ptr);ß

}


static void DisplayUsage(){
  fprintf(stderr,"wiswesser -i<format> -o<format> <input>\n");
}



static unsigned int ProcessCommandLine(int argc, char *argv[]){

  const char *ptr=0; 
  int i,j; 

  inpformat = (const char*)0; 
  outformat = (const char*)0; 

  if (argc < 2)
    DisplayUsage();

  j=0; 
  for (i=1;i<argc;i++){

    ptr = argv[i];

    if (ptr[0]=='-' && ptr[1]){
      switch (ptr[1]){


        case 'i':
          



        default:
          fprintf(stderr,"Error: Unrecognised letter option - %c\n",ptr[1]);
          break;
      }
      


    }
      


    // end argc loop
  }

  // check cfx file
  const char *cfx_ext = strrchr(cfxname, '.');
  if (!cfx_ext) {
    fprintf(stderr,"ERROR: Could not recognise .cfx(2) file extension %s\n",cfxname);
    exit(1);
  }
  else{
    if(strcmp(cfx_ext+1,"cfx2")==0)
      twolevel = 1;
  }

  if (j < 2)
    DisplayUsage();
  
  if (isDirectory(inpname))
    dflag = 1;

  if (rflag && !dflag){
    fprintf(stderr, "ERROR: Recursive search cannot be enabled for a non-dir\n");
    DisplayUsage();
  }

  // stdin linux
  if (!isatty(1))
    hflag = false; 

  // precedence flags
  if (oflag && hflag){
    fprintf(stderr,"WARNING: -o and -h are incompatible, proceeding with -o\n");
    hflag = 0;
  }

  if (cflag && hflag){
    fprintf(stderr,"WARNING: -c and -h are incompatible, proceeding with -c\n");
    hflag = 0;
  }

  // incompatable flags
  if (!twolevel && dfsflag){
    fprintf(stderr,"ERROR: -dfs is not a valid mode for one-level matching\n");
    DisplayUsage();
  }
}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc,argv);


#ifdef WORK

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


#endif



  return 0;
}
