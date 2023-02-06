//
// Created by Michael Blakey on 04/07/2022.
//

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <string>
#include <vector>

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


std::vector<const char *> file_queue;
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

// access global strings for type def 
static bool ParseChemicalNotation(){

  //fprintf(stderr,"input is: %s\n",inpname);

  // smiles input block
  if (!strcmp(inpformat,"smi")){

    //for (const char * str : file_queue)
      //ConvertSMI(str, outformat);

  }

  // wln input block
  else if (!strcmp(inpformat,"wln")){

    for (const char * str : file_queue){
      OpenBabel::OBMol *Mol = new OpenBabel::OBMol; 
      ConvertWLN(str, outformat, Mol);
      delete Mol;
      Mol = 0;
    }
      
  }

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
      case 0: file_queue.push_back(ptr); break;
      default: break;
    }
      
    // end argc loop
  }

  if (file_queue.empty()){
    fprintf(stderr,"Error: no input string(s) | file(sß) given for parsing\n");
    DisplayUsage();
  }

}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc,argv);

  
  // file check in here
  ParseChemicalNotation();


  return 0;
}
