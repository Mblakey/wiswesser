#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <map>
#include <vector> // can optimise this out later
#include <string> // just for prototype

unsigned int opt_mode = 0;
unsigned int opt_verbose = false;
const char *input;

bool ReadLineFromFile(FILE *fp, char *buffer, unsigned int n, bool add_nl=true){
  char *end = buffer+n;
  char *ptr;
  int ch;

  ptr = buffer;
  do {
    ch = getc_unlocked(fp); // this increments fp
    if (ch == '\n') {
      if (add_nl)
        *ptr++ = '\n'; // if i want the newline or not
      *ptr = '\0';
      return true;
    }
    if (ch == '\f') {
      *ptr++ = '\n';
      *ptr = '\0';
      return true;
    }
    if (ch == '\r') {
      *ptr++ = '\n';
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
      *ptr++ = '\n';
      *ptr = '\0';
      return ptr-buffer > 1;
    }
    *ptr++ = ch;
  } while (ptr < end);
  *ptr = 0;
  
  fprintf(stderr, "Warning: line too long!\n");
  return false;
}

static void DisplayUsage()
{
  fprintf(stderr, "compresswln <options> <input file> > <out>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -c          compress input file\n");
  fprintf(stderr, "  -d          decompress input file\n");
  fprintf(stderr, "  -v          verbose debugging statements on\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- WLN Compression ---\n");
  fprintf(stderr, "This exec writes a wln file into a 6 bit representation, and can perform\n"
                  "various compression schemes which are selected in options.\n"
                  "This is part of michaels PhD investigations into compressing chemical strings.\n\n");
  DisplayUsage();
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  input = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        
        case 'c': 
          opt_mode = 1; 
          break;

        case 'd':
          opt_mode = 2; 
          break;

        case 'v':
          opt_verbose = true;
          break;

        case 'h':
          DisplayHelp();

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
    else{
      switch(j++){
        case 0:
          input = ptr; 
          break;
        default:
          fprintf(stderr,"Error: multiple files not currently supported\n");
          exit(1);
      }
    }
  }

  if(!input){
    fprintf(stderr,"Error: no input file given\n");
    DisplayUsage();
  }

  return;
}


void initialise_maps(  std::map<unsigned char,unsigned int> &encode, 
                      std::map<unsigned int, unsigned char> &decode)
{
  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 -/&";
  unsigned int j=1;
  for (unsigned int i=0;i<40;i++){
    if(opt_verbose)
      fprintf(stderr,"%c --> %d\n",wln[i],j);
  
    encode[wln[i]] = j;
    decode[j] = wln[i];
    j++;
  }
}





int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);

  std::map<unsigned char,unsigned int> encode;
  std::map<unsigned int, unsigned char> decode;

  initialise_maps(encode,decode);

  


  return 0;
}