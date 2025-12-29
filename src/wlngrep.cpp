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
#define COUNT_ONLY  3

int opt_match_option;  
unsigned char latty;
unsigned long nmatches; 

FILE *ifp; 

struct fsm_state {
  bool final; 
  unsigned short jmp[256];
};  


static void make_accept(struct fsm_state *blk, unsigned int id)
{
  blk[id].final = true; 
}


static inline bool is_accept(struct fsm_state *state)
{
  return state->final; 
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
  memset(blk, 0, sizeof(struct fsm_state)*64); 

  unsigned int nstates = 0; 
  const unsigned int neg_state = nstates++; /* allows the memset for false match */
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


static void handle_match(char *buffer, int match_end, int consumed)
{
  switch (opt_match_option)  {
    case LINE_MATCH:
      if (latty) printf(RED); 
      fwrite(buffer, match_end, 1, stdout); 
      if (latty) printf(RESET); 

      if (consumed > match_end)
        fwrite(buffer+match_end, consumed-match_end, 1, stdout); 
      break; 

    case MATCH_ONLY: 
      fwrite(buffer, match_end, 1, stdout); 
      fputc('\n', stdout); 
      break; 

    case COUNT_ONLY: nmatches++; break; 
  }
}


/* uses goto to keep logic clean, nothing to be scared of */
static bool process_file(struct fsm_state *wlnfsm, FILE *fp)
{
  const size_t bufsize = 16384; 
  const size_t wlnbuf_size = 4096; 

  char buffer[bufsize]; 
  char wln_buf[wlnbuf_size]; 
  unsigned int wlnlen = 0; 

  fsm_state *root = &wlnfsm[1]; 
  fsm_state *head = root; 
  unsigned int jmpid; 
  unsigned int match_end = 0; 

  nmatches = 0; 
  
  while (true) {
jmp_fetch_data:
    size_t bytes = fread(buffer, 1, bufsize, fp); 
    if (!bytes) {
      if (match_end > 0) 
        handle_match(wln_buf, match_end+1, wlnlen);  
      return true;
    }

    unsigned int i = 0; 

    if (head != root)
      goto jmp_matching; 

jmp_not_matching:
    for (; i < bytes; i++) {
      unsigned char ch = buffer[i]; 
      jmpid = head->jmp[ch]; 
      if (jmpid) {
        wlnlen    = 0; 
        match_end = 0; 
        head = &wlnfsm[jmpid]; 
        wln_buf[wlnlen++] = ch; 
        i++; // skip to next char
        goto jmp_matching;
      }
      else if (opt_match_option == LINE_MATCH) 
        fputc(ch, stdout);
    }
    goto jmp_fetch_data; 

jmp_matching:
    for (; i < bytes; i++) { 
      unsigned char ch = buffer[i]; 
      jmpid = head->jmp[ch]; 
      if (!jmpid) {
        head = root; 
        if (match_end > 0) 
          handle_match(wln_buf, match_end+1, wlnlen);  
        match_end = 0; 
        goto jmp_not_matching;
      }
      else if (wlnlen >= wlnbuf_size) {
        fprintf(stderr, "Warning: WLN string too long!\n");
        goto jmp_not_matching;
      }
      else {
        head = &wlnfsm[jmpid]; 
        wln_buf[wlnlen] = ch; 
        if (is_accept(head))
          match_end = wlnlen; 
        wlnlen++; 
      }
    }
    goto jmp_fetch_data; 
  }
  
  return true;
}


static void display_usage()
{
  fprintf(stderr, "usage: wlngrep <options> <file>\n");
  fprintf(stderr, "options:\n");
  fprintf(stderr, "-c|--only-count        return number of matches instead of string\n");
  fprintf(stderr, "-o|--only-match        print only the matched parts of line\n");
  exit(1);
}


static void process_cml(int argc, char *argv[])
{
  const char *ptr;
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
      case 'c': opt_match_option = COUNT_ONLY; break;
      case 'o': opt_match_option = MATCH_ONLY; break;

      case '-':
        if (!strcmp(ptr,"--only-count"))
          opt_match_option = COUNT_ONLY;
        else if (!strcmp(ptr,"--only-match"))
          opt_match_option = MATCH_ONLY;
        else {
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          display_usage();
        }
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
  latty = isatty(STDOUT_FILENO) ? 0xff:0x0; 
  process_cml(argc,argv); 

  struct fsm_state *wlnfsm = create_wlnfsm(); 

  process_file(wlnfsm, ifp);  
  if (opt_match_option == COUNT_ONLY)
    printf("%lu matches\n", nmatches); 
  
  free(wlnfsm); 
  if (ifp != stdin) 
    fclose(ifp); 
  return 0;
}


