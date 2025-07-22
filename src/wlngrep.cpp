#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "wlnfsm.h" // generated file

// 0 - return whole line, 
// 1 - return matches only, 
// 2 - exact match only, 
// 3- invert exact match

FILE *fp;
unsigned char latty = 0;

unsigned int opt_match_option = 0; 
unsigned int opt_count        = 0;
unsigned int opt_invert_match  = 0;

unsigned int matches = 0; 
unsigned int lines_parsed = 0; 


unsigned char 
readline(FILE *fp, char *buffer, unsigned int n, char add_nl){
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
      return 1;
    }
    if (ch == '\f') {
      *ptr++ = '\n';
      *ptr = '\0';
      return 1;
    }

    if (ch == '\r') {
      *ptr++ = '\n';
      *ptr = '\0';
      ch = getc_unlocked(fp);
      if (ch != '\n') {
        if (ch == -1)
          return 0;
        ungetc(ch,fp);
      }
      return 1;
    }
    if (ch == -1) {
      *ptr++ = '\n';
      *ptr = '\0';
      return ptr-buffer > 1;
    }
    *ptr++ = ch;
  } while (ptr < end);
  *ptr = 0;
  
  fprintf(stderr, "Error: line too long for buffer - %d\n", n);
  return 0;
}


static void 
MatchBuffer(char *buffer, 
            char *match_map,
            struct fsm_state *fsm) 
{
  fsm_state *head = &fsm[1]; 
  const unsigned int len = strlen(buffer); 

  for (unsigned int i=0; i<len; i++) {
    unsigned char ch = buffer[i];
    if (head->jmp[ch] != 0) for (unsigned int j=i;j<=len;j++) {
      ch = buffer[j];
      unsigned int jmpid = head->jmp[ch];
      if (!jmpid) {
        fprintf(stderr, "failed? - %d\n", head->final); 
        if (!head->final) 
          head = &fsm[1]; // reset to root
        else {
          i = j; // greedy match
          fprintf(stderr, "matched\n");
        }
      }
      else {
        head = &fsm[jmpid]; 
        fprintf(stderr, "moved on %c\n", ch); 
      }
    }
  }
}


static bool 
process_file(FILE *fp, struct fsm_state *fsm)
{
  char buffer[4096];  
  char match_map[4096];  
  while (readline(fp, buffer, 4096, false)){
    lines_parsed++;
    MatchBuffer(buffer, match_map, fsm);  
  }
  
  if(opt_count)
    fprintf(stderr,"%d matches\n",matches);
  return true;
}


static void 
display_usage()
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

static void 
process_cml(int argc, char *argv[])
{
  const char *ptr = 0;
  fp = NULL; 
  int i, j;
  j = 0;
  for (i = 1; i < argc; i++) {
    ptr = argv[i];
    if (ptr[0] == '-' && !ptr[1]) {
      fp = stdin; 
      j++;
    }
    if (ptr[0] == '-' && ptr[1]) switch (ptr[1]) {
      case 'c': opt_count = 1; break;
      case 'o': opt_match_option = 1; break;
      case 'x': opt_match_option = 2; break;
      case 'v': opt_invert_match = 1; break;

      case '-':
        if (!strcmp(ptr,"--only-count"))
          opt_count = 1;
        else if (!strcmp(ptr,"--only-match"))
          opt_match_option = 1;
        else if (!strcmp(ptr,"--exact-match"))
          opt_match_option = 2;
        else if (!strcmp(ptr,"--invert-match"))
          opt_invert_match = true;
        break;

      default:
        fprintf(stderr, "Error: unrecognised input %s\n", ptr);
        display_usage();
    }
    else switch (j++) {
      case 0: 
        fp = fopen(ptr, "r");  
        if (!fp) {
          fprintf(stderr, "Error: could not open file at %s\n", ptr); 
          display_usage(); 
        }
        break; 
    }
  }

  if (!fp) 
    fp = stdin; 
  return;
}


int 
main(int argc, char* argv[])
{
  process_cml(argc,argv); 
  fsm_state *wlnfsm = generate_WLN_fsm();  
  
  process_file(fp, wlnfsm);  
  if (fp != stdin) 
    fclose(fp); 
  free(wlnfsm); 
  return 0;
}
