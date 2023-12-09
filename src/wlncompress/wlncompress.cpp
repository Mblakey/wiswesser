

#include <stdlib.h>
#include <stdio.h>


unsigned int opt_mode = 0;
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
  fprintf(stderr, "  -c                    compress input file\n");
  fprintf(stderr, "  -d                    decompress input file\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- WLN Compression ---\n\n");
  fprintf(stderr, " This exec writes a wln file into a 6 bit representation, and can perform\n"
                  " various compression schemes which are selected in options\n"
                  " this is part of michaels PhD investigations into compressing chemical strings\n");
  DisplayUsage();
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i;

  input = (const char *)0;

  if (argc < 2)
    DisplayUsage();

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

        case 'h':
          DisplayHelp();


       


        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
  }

  if(!input){
    fprintf(stderr,"Error: no input string entered\n");
    DisplayUsage();
  }

  return;
}


int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);


  return 0;
}