

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>


#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"

const char *filename;
unsigned int lines_parsed = 0; 

unsigned int opt_dump         = 0;
unsigned int opt_match_option = 0; // 0 - return whole line, 1 - return matches only, 2 - exact match only, 3- invert exact match
unsigned int opt_count        = 0;
unsigned int opt_string_file  = 0;

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


static bool MatchFile(FILE *fp,FSMAutomata *machine){

  unsigned int matches = 0; 
  char *buffer = (char*)malloc(sizeof(char) * BUFF_SIZE+1);
  memset(buffer,0,BUFF_SIZE+1);

  while(ReadLineFromFile(fp,buffer,BUFF_SIZE,false)){
    lines_parsed++;
    matches += DFAGreedyMatchLine(buffer,machine,isatty(1),opt_match_option, opt_count);
  }
  
  if(opt_count)
    fprintf(stderr,"%d matches\n",matches);

  free(buffer);
  return true;
}

static void DisplayUsage()
{
  fprintf(stderr, "usage: wlngrep <options> <file>\n");
  fprintf(stderr, "options:\n");
  fprintf(stderr, "-c|--only-count        return number of matches instead of string\n");
  fprintf(stderr, "-d|--dump              dump resultant machine to dot file\n");
  fprintf(stderr, "-o|--only-match        print only the matched parts of line\n");
  fprintf(stderr, "-s|--string            interpret <file> as a string to match\n");
  fprintf(stderr, "-x|--exact-match       return string if whole line matches\n");
  fprintf(stderr, "-v|--invert-match      return string if whole line does not match\n");
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i, j;

  filename = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];

    if (ptr[0] == '-' && ptr[1])
      switch (ptr[1])
      {

        case 'c':
          opt_count = 1;
          break;

        case 'd':
          opt_dump = 1; 
          break;

        case 'o':
          opt_match_option = 1; 
          break;

        case 's':
          opt_string_file = 1;
          break;

        case 'x':
          opt_match_option = 2;
          break;

        case 'v':
          opt_match_option = 3;
          break;

        case '-':
          if(!strcmp(ptr,"--only-count"))
            opt_count = 1;
          else if(!strcmp(ptr,"--dump"))
            opt_dump = 1;
          else if(!strcmp(ptr,"--only-match"))
            opt_match_option = 1;
          else if(!strcmp(ptr,"--exact-match"))
            opt_match_option = 2;
          else if(!strcmp(ptr,"--invert-match"))
            opt_match_option = 2;
          else if(!strcmp(ptr,"--string"))
            opt_string_file = 1;

          break;

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }

    else
      switch (j++)
      {
      case 0:
        filename = ptr;
        break;
      }
  }

  if(!filename){
    fprintf(stderr,"Error: not enough args\n");
    DisplayUsage();
  }

  return;
}

int main(int argc, char* argv[])
{
  ProcessCommandLine(argc,argv); 

  FSMAutomata *wlnDFA = CreateWLNDFA(REASONABLE,REASONABLE);
  if(!wlnDFA || wlnDFA->type != DFA)
    return 1;

  if(opt_dump){
    fprintf(stderr,"machines dumped, exiting\n");
    delete wlnDFA;
    return 0;
  }

  if(!opt_string_file){
    FILE *fp = fopen(filename,"r");
    if(!fp){
      fprintf(stderr,"Error: unable to open file at: %s\n",filename);
      return 1; 
    }

    MatchFile(fp,wlnDFA);
    fprintf(stderr,"%d lines parsed\n",lines_parsed);

    fclose(fp);
  }
  else{
    unsigned int matches = DFAGreedyMatchLine(filename,wlnDFA,isatty(0),opt_match_option,opt_count);
    if(opt_count)
      fprintf(stderr,"%d matches\n",matches);
  }

  delete wlnDFA; 
  return 0;
}
