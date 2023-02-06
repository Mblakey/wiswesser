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

const char *inpname; 
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


// access global strings for type def 
static bool ParseChemicalNotation(){

  fprintf(stderr,"input is: %s\n",inpname);
  

}


static bool ReadInpFormat(const char *ptr){
  ptr+=2; 
  
  if (!strcmp(ptr,"wln") || !strcmp(ptr,"WLN")){
    fprintf(stderr,"Runtime: setting wln as input format\n");
    inpformat = "wln"; 
    return true;
  }
  else if (!strcmp(ptr,"smi") || !strcmp(ptr,"smiles") || !strcmp(ptr,"SMI")){
    fprintf(stderr,"Runtime: setting smiles as input format\n");
    inpformat = "smi"; 
    return true;
  }
  else {
    fprintf(stderr,"Runtime: unrecognised format entered for input format - %s\n",ptr);
    return false;
  }


  return false;
}


static bool ReadOutFormat(const char *ptr){
  ptr+=2; 
  
  if (!strcmp(ptr,"wln") || !strcmp(ptr,"WLN")){
    fprintf(stderr,"Runtime: setting wln as output format\n");
    outformat = "wln"; 
    return true;
  }
  else if (!strcmp(ptr,"smi") || !strcmp(ptr,"smiles") || !strcmp(ptr,"SMI")){
    fprintf(stderr,"Runtime: setting smiles as output format\n");
    outformat = "smi"; 
    return true;
  }
  else {
    fprintf(stderr,"Runtime: unrecognised format entered for output format - %s\n",ptr);
    return false;
  }

  return false;
}


static void DisplayUsage(){
  fprintf(stderr,"wiswesser -i<format> -o<format> <input>\n");
  exit(1);
}



static unsigned int ProcessCommandLine(int argc, char *argv[]){

  const char *ptr=0; 
  int i,j; 

  inpformat = (const char*)0; 
  outformat = (const char*)0;
  inpname = (const char*)0; 

  if (argc < 2)
   DisplayUsage();

  j=0; 
  for (i=1;i<argc;i++){

    ptr = argv[i];

    if (ptr[0]=='-' && ptr[1]){
      switch (ptr[1]){
        case 'i':
          if(!ReadInpFormat(ptr))
            DisplayUsage();
          break;
        case 'o':
          if(!ReadOutFormat(ptr))
            DisplayUsage();
          break;
        default:
          fprintf(stderr,"Error: Unrecognised letter option - %c\n",ptr[1]);
          break;
      }
      
    }
    else switch(j++){
      case 0: inpname = ptr; break;
      default: break;
    }
      
    // end argc loop
  }

  if (!inpname){
    fprintf(stderr,"Error: no input string | file given for parsing\n");
    DisplayUsage();
  }

}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc,argv);

  
  // file check in here
  ParseChemicalNotation(inpname);


  return 0;
}
