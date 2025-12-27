#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define RESET  "\e[0;0m"
#define RED    "\e[0;31m"

// 0 - return whole line, 
// 1 - return matches only, 
// 2 - exact match only, 

#define LINE_MATCH  0
#define MATCH_ONLY  1
#define EXACT_MATCH 2

unsigned char latty = 0;

int opt_match_option;  
bool opt_invert_match; 
bool opt_count; 

unsigned int matches = 0; 
unsigned int lines_parsed = 0; 

FILE *ifp; 

struct fsm_state {
  bool final; 
  unsigned short jmp[256];
};  


static void make_accept(struct fsm_state *blk, unsigned int id)
{
  blk[id].final = true; 
}


static void add_transition(struct fsm_state *blk, 
                           unsigned int dst_id, 
                           unsigned int src_id, 
                           unsigned char ch)
{
  blk[dst_id].jmp[ch] = src_id; 
}


static void add_non_cyclic_transitions(struct fsm_state *blk, 
                                       unsigned int dst_id, 
                                       unsigned int src_id)
{
  blk[dst_id].jmp['B'] = src_id; 
  blk[dst_id].jmp['C'] = src_id; 
  blk[dst_id].jmp['E'] = src_id; 
  blk[dst_id].jmp['F'] = src_id; 
  blk[dst_id].jmp['G'] = src_id; 
  blk[dst_id].jmp['H'] = src_id; 
  blk[dst_id].jmp['I'] = src_id; 
  blk[dst_id].jmp['K'] = src_id; 
  blk[dst_id].jmp['M'] = src_id; 
  blk[dst_id].jmp['N'] = src_id; 
  blk[dst_id].jmp['O'] = src_id; 
  blk[dst_id].jmp['P'] = src_id; 
  blk[dst_id].jmp['Q'] = src_id; 
  blk[dst_id].jmp['R'] = src_id; 
  blk[dst_id].jmp['S'] = src_id; 
  blk[dst_id].jmp['U'] = src_id; 
  blk[dst_id].jmp['V'] = src_id; 
  blk[dst_id].jmp['W'] = src_id; 
  blk[dst_id].jmp['X'] = src_id; 
  blk[dst_id].jmp['Y'] = src_id; 
  blk[dst_id].jmp['Z'] = src_id; 


  for (unsigned char ch='0'; ch <= '9'; ch++)
    blk[dst_id].jmp[ch] = src_id; 
}


static struct fsm_state* create_wlnfsm() 
{
  struct fsm_state *blk = (struct fsm_state*)malloc(sizeof(struct fsm_state)*64);
  unsigned int nstates = 0; 
  const unsigned int head = nstates++;  
  const unsigned int non_cyclic = nstates++; 
  const unsigned int cyclic = nstates++; 
  
  make_accept(blk, non_cyclic); 
  add_non_cyclic_transitions(blk, head, non_cyclic); 
  add_non_cyclic_transitions(blk, non_cyclic, non_cyclic); 
  add_transition(blk, non_cyclic, non_cyclic, '&'); 
  
  const unsigned int non_cyclic_dash_open  = nstates++; 
  const unsigned int non_cyclic_dash_elem  = nstates++; 
  const unsigned int non_cyclic_dash_close = nstates++; 

  /* generic elemental characters */
  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) {
    add_transition(blk, non_cyclic_dash_elem, non_cyclic_dash_open, ch); 
    add_transition(blk, non_cyclic_dash_close, non_cyclic_dash_elem, ch); 
  }
  
  add_transition(blk, non_cyclic_dash_open, non_cyclic, '-'); 
  add_transition(blk, non_cyclic, non_cyclic_dash_close, '-'); 

  return blk; 
}


static unsigned char readline(char *buffer, 
                              unsigned int n, 
                              FILE *fp)
{
  char *end = buffer+n;
  char *ptr;
  int ch;

  ptr = buffer;
  do {
    ch = getc_unlocked(fp); // this increments fp
    if (ch == '\n') {
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



static bool match_buffer(struct fsm_state *wlnfsm,
                         char *buffer, 
                         char *match_map,
                         unsigned int len)
{
  bool any_match = false;
  fsm_state *root = &wlnfsm[0]; 
  fsm_state *head = root; 

  for (unsigned int i=0; i<len; i++) {
    unsigned char ch = buffer[i];
    if (head->jmp[ch] != 0) for (unsigned int j=i;j<=len;j++) {
      ch = buffer[j];
      unsigned int jmpid = head->jmp[ch];
      if (!jmpid) {
        if (!head->final) 
          memset(&match_map[i], 0, len-i); 
        else {
          i = j; // greedy match
          any_match = true;
        }
        head = root; 
        break;
      }
      else {
        head = &wlnfsm[jmpid];
        match_map[j] = 1;
      }
    }
  }
  return any_match; 
}


static void print_buffer(char *buffer, 
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


static bool process_file(struct fsm_state *wlnfsm, FILE *fp)
{
  char buffer[4096];  
  char match_map[4096];  
  unsigned int len; 
  while (readline(buffer, 4096, fp)) {
    lines_parsed++;
    len = strlen(buffer);
    memset(match_map, 0, 4096);
    if (match_buffer(wlnfsm, buffer, match_map, len)) {
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


static void display_usage()
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


static void process_cml(int argc, char *argv[])
{
  const char *ptr;
  opt_count = 0; 
  opt_invert_match = 0; 
  opt_match_option = LINE_MATCH; 

  int i, j;
  j = 0;
  for (i = 1; i < argc; i++) {
    ptr = argv[i];
    if (ptr[0] == '-' && !ptr[1]) {
      ifp = stdin; 
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
        ifp = fopen(ptr, "r");  
        if (!ifp) {
          fprintf(stderr, "Error: could not open file at %s\n", ptr); 
          display_usage(); 
        }
        break; 
    }
  }
  
  if (!j)
    ifp = stdin;
  return;
}


int main(int argc, char* argv[])
{
  process_cml(argc,argv); 

  if (isatty(STDOUT_FILENO))
    latty = 0xff;
  else 
    latty = 0x0; 

  struct fsm_state *wlnfsm = create_wlnfsm(); 
  process_file(wlnfsm, ifp);  
  
  free(wlnfsm); 
  if (ifp != stdin) 
    fclose(ifp); 
  return 0;
}
