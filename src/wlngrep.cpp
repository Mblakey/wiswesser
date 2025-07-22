#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "readline.h"

#define RESET  "\e[0;0m"
#define RED    "\e[0;31m"

// 0 - return whole line, 
// 1 - return matches only, 
// 2 - exact match only, 

#define LINE_MATCH  0
#define MATCH_ONLY  1
#define EXACT_MATCH 2

FILE *fp;
unsigned char latty = 0;

unsigned int opt_match_option;  
unsigned int opt_invert_match; 
unsigned int opt_count; 

unsigned int matches = 0; 
unsigned int lines_parsed = 0; 


struct fsm_state {
  unsigned char final; 
  unsigned short jmp[256];
};


static void 
fsm_state_make_final(fsm_state *s) {
  s->final = 1; 
}


static void 
fsm_state_add_transition(fsm_state *s, 
                         unsigned short trg_id, 
                         unsigned char ch)
{
  s->jmp[ch] = trg_id;
}


static void 
fsm_state_add_transition_range(fsm_state *s, 
                               unsigned short trg_id, 
                               unsigned char s_ch, 
                               unsigned char t_ch)
{
  for (unsigned char ch=s_ch; ch <= t_ch;ch++) 
    s->jmp[ch] = trg_id;
}


static struct fsm_state*
wlnmatcher_alloc() {
  struct fsm_state *fsm = (struct fsm_state*)malloc(sizeof(struct fsm_state) * 300);
  memset(fsm, 0, sizeof(struct fsm_state) * 300); 

#define INIT   1
#define BRANCH 2

#define RING_OPEN 3
#define RING_SSSR 4
#define RING_CLOSE 5

#define BRANCH_WC  6; // after a cycle logic is expanded

  fsm_state *init       = &fsm[INIT]; 
  fsm_state *branch     = &fsm[BRANCH]; 
  fsm_state *ring_open  = &fsm[RING_OPEN]; 
  fsm_state *ring_SSSR  = &fsm[RING_SSSR]; 
  fsm_state *ring_close = &fsm[RING_CLOSE]; 
  
  fsm_state_make_final(branch); 
  fsm_state_make_final(ring_close); 

  fsm_state_add_transition_range(init, BRANCH, '0', '9'); 
  fsm_state_add_transition(init, BRANCH, 'B'); 
  fsm_state_add_transition(init, BRANCH, 'C'); 

  fsm_state_add_transition_range(init, BRANCH, 'E', 'I'); 
  fsm_state_add_transition(init, BRANCH, 'K'); 
  fsm_state_add_transition_range(init, BRANCH, 'M', 'Q'); 
  fsm_state_add_transition(init, BRANCH, 'S'); 
  fsm_state_add_transition_range(init, BRANCH, 'V', 'Z'); 

  fsm_state_add_transition(init, RING_OPEN, 'L'); 
  fsm_state_add_transition(init, RING_OPEN, 'T');

  fsm_state_add_transition_range(branch, BRANCH, '0', '9'); 
  fsm_state_add_transition(branch, BRANCH, 'B'); 
  fsm_state_add_transition(branch, BRANCH, 'C'); 

  fsm_state_add_transition_range(branch, BRANCH, 'E', 'I'); 
  fsm_state_add_transition(branch, BRANCH, 'K'); 
  fsm_state_add_transition_range(branch, BRANCH, 'M', 'Q'); 
  fsm_state_add_transition(branch, BRANCH, 'S'); 
  fsm_state_add_transition(branch, BRANCH, 'U'); 
  fsm_state_add_transition_range(branch, BRANCH, 'V', 'Z'); 


  fsm_state_add_transition(branch, BRANCH, '&'); 

  fsm_state_add_transition_range(ring_open, RING_SSSR, '0', '9');
  fsm_state_add_transition_range(ring_SSSR, RING_SSSR, '0', '9');
  fsm_state_add_transition(ring_SSSR, RING_CLOSE, 'J');
  return fsm;
}


static bool
match_buffer(char *buffer, 
             char *match_map,
             unsigned int len,
             struct fsm_state *fsm) 
{
  bool any_match = false;
  fsm_state *root = &fsm[1]; 
  fsm_state *head = root; 

  for (unsigned int i=0; i<len; i++) {
    unsigned char ch = buffer[i];
    if (head->jmp[ch] != 0) for (unsigned int j=i;j<=len;j++) {
      ch = buffer[j];
      unsigned int jmpid = head->jmp[ch];
      if (!jmpid) {
        if (!head->final) {
          memset(&match_map[i], 0, len-i); 
        }
        else {
          i = j; // greedy match
          any_match = true;
        }
        head = root; 
        break;
      }
      else {
        head = &fsm[jmpid];
        match_map[j] = 1;
      }
    }
  }
  return any_match; 
}


static void 
print_buffer(char *buffer, 
             char *match_map,
             unsigned int len, 
             unsigned char inverse) 
{
  char matching = 0; 

  if (!inverse) for (unsigned int i=0; i<len; i++) {
    if (match_map[i]) { 
      if (!matching) {
        matching = 1;
        if (latty) printf(RED);
      }
    }
    else if (matching) {
      matching = 0; 
      if (latty) printf(RESET);
    }
    fputc(buffer[i], stdout); 
  }
  else for (unsigned int i=0; i<len; i++) {
    if (!match_map[i]) { 
      if (!matching) {
        matching = 1;
        if (latty) printf(RED);
      }
    }
    else if (matching) {
      matching = 0; 
      if (latty) printf(RESET);
    }
    fputc(buffer[i], stdout); 
  }
  fputc('\n', stdout); 
}


static void 
print_matches(char *buffer, 
              char *match_map,
              unsigned int len, 
              unsigned char inverse) 
{
  char matching = 0;
  if (!inverse) for (unsigned int i=0; i<len; i++) {
    if (match_map[i]) {
      if (!matching)
        matching = 1;
      fputc(buffer[i], stdout); 
    }
    else if (matching) {
      matching = 0;
      fputc('\n', stdout); 
    }
  }
  else for (unsigned int i=0; i<len; i++) {
    if (!match_map[i]) {
      if (!matching)
        matching = 1;
      fputc(buffer[i], stdout); 
    }
    else if (matching) {
      matching = 0;
      fputc('\n', stdout); 
    }
  }
}


static bool 
process_file(FILE *fp, struct fsm_state *fsm)
{
  char buffer[4096];  
  char match_map[4096];  
  unsigned int len; 
  while (readline(fp, buffer, 4096, false)){
    lines_parsed++;
    len = strlen(buffer);
    memset(match_map, 0, 4096);
    if (match_buffer(buffer, match_map, len, fsm)) {
      
      switch (opt_match_option) {
        case LINE_MATCH: print_buffer(buffer, match_map, len, opt_invert_match); break;
        case MATCH_ONLY: print_matches(buffer, match_map, len, opt_invert_match); break; 
        case EXACT_MATCH: break; 
      }

    }
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
  opt_count = 0; 
  opt_invert_match = 0; 
  opt_match_option = LINE_MATCH; 

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
      case 'o': opt_match_option = MATCH_ONLY; break;
      case 'x': opt_match_option = EXACT_MATCH; break;
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

  if (isatty(STDOUT_FILENO))
    latty = 0xff;
  else 
    latty = 0x0; 

  struct fsm_state *wlnfsm = wlnmatcher_alloc(); 
  process_file(fp, wlnfsm);  
  if (fp != stdin) 
    fclose(fp); 
  free(wlnfsm); 
  return 0;
}
