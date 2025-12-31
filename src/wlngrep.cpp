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

#define NO_JMP 0xffff

int opt_match_option;  
unsigned char latty;
unsigned long nmatches; 

FILE *ifp; 

struct fsm_state {
  bool final; 
  unsigned short jmp[256]; 
};  

unsigned short nstates;
struct fsm_state wlnfsm[64]; 

static unsigned short state_init()
{
  wlnfsm[nstates].final = false; 
  memset(wlnfsm[nstates].jmp, 0xff, sizeof(unsigned short)*256); 
  return nstates++; 
}

#define is_accept(x) wlnfsm[x].final
#define make_accept(x) wlnfsm[x].final = true
#define make_jmp(dst,src,ch) wlnfsm[src].jmp[ch] = dst
#define state_jmp(x,ch) wlnfsm[x].jmp[ch]


static void make_symbol_transitions(unsigned int dst_id, unsigned int src_id)
{
	make_jmp(dst_id, src_id, 'B');

#ifdef WLN_C_MATCH
  // The C symbol in WLN is rare, and adding C as a potential 
  // match token directly conflicts with SMILES. e.g most normal 
  // organic SMILES will match due to C chains. Omitting this token
  // here as default, can be recompiled to match C symbols.
	
  make_jmp(dst_id, src_id, 'C');
#endif
	
  make_jmp(dst_id, src_id, 'E');
	make_jmp(dst_id, src_id, 'F');
	make_jmp(dst_id, src_id, 'G');
	make_jmp(dst_id, src_id, 'H');
	make_jmp(dst_id, src_id, 'I');
	make_jmp(dst_id, src_id, 'K');
	make_jmp(dst_id, src_id, 'M');
	make_jmp(dst_id, src_id, 'N');
	make_jmp(dst_id, src_id, 'O');
	make_jmp(dst_id, src_id, 'P');
	make_jmp(dst_id, src_id, 'Q');
	make_jmp(dst_id, src_id, 'R');
	make_jmp(dst_id, src_id, 'S');
	make_jmp(dst_id, src_id, 'U');
	make_jmp(dst_id, src_id, 'V');
	make_jmp(dst_id, src_id, 'W');
	make_jmp(dst_id, src_id, 'X');
	make_jmp(dst_id, src_id, 'Y');
	make_jmp(dst_id, src_id, 'Z');
}


static void init_wlnfsm() 
{
  nstates = 0; 
  const unsigned int head        = state_init(); 

  const unsigned int chain_open  = state_init(); 
  const unsigned int locant_space = state_init(); 
  
  make_accept(chain_open); 
  make_symbol_transitions(chain_open, head); 
  make_symbol_transitions(chain_open, chain_open); 

  for (unsigned char ch='0'; ch <= '9'; ch++) {
    make_jmp(chain_open, head, ch); 
    make_jmp(chain_open, chain_open, ch); 
  }

  make_jmp(chain_open, chain_open, '&'); 
  
  const unsigned int dash_open  = state_init();
  const unsigned int dash_elem  = state_init();
  const unsigned int dash_close = state_init();

  /* generic elemental characters */
  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) {
    make_jmp(dash_elem, dash_open, ch); 
    make_jmp(dash_close, dash_elem, ch); 
  }
  
  make_jmp(dash_open, head, '-'); 
  make_jmp(dash_open, chain_open, '-'); 
  make_jmp(chain_open, dash_close, '-'); 

  /* cyclic internal block states */

  const unsigned int ring_open   = state_init(); 
  const unsigned int ring_close  = state_init(); 
  const unsigned int ring_digits = state_init(); 

  make_accept(ring_close); 

  make_jmp(ring_open, head, 'L');  
  make_jmp(ring_open, head, 'T');  

  for (unsigned char ch = '0'; ch <= '9'; ch++) {
    make_jmp(ring_digits, ring_open, ch);
    make_jmp(ring_digits, ring_digits, ch);
  }

  make_jmp(ring_close, ring_digits, 'J'); 

  const unsigned int digit_space  = state_init(); 
  const unsigned int digit_locant = state_init(); 
  
  /* ring digit locants and bridge atoms */
  make_jmp(digit_space, ring_open, ' '); 
  make_jmp(digit_space, ring_digits, ' '); 
  
  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) 
    make_jmp(digit_locant, digit_space, ch); 
  make_jmp(digit_locant, digit_locant, '&'); 

  make_jmp(digit_space, digit_locant, ' '); 

  for (unsigned char ch = '0'; ch <= '9'; ch++) 
    make_jmp(ring_digits,digit_locant, ch); 
  
  make_jmp(ring_close, digit_locant, 'J'); 

  /* multicyclic ring block defintion */
  // n, locants, space, size

  const unsigned int multi_digit    = state_init(); 
  const unsigned int multi_locants  = state_init(); 
  const unsigned int multi_break    = state_init(); 
  const unsigned int multi_size     = state_init(); 
  for (unsigned char ch = '0'; ch <= '9'; ch++) {
    make_jmp(multi_digit, digit_space, ch); 
    make_jmp(multi_digit, multi_digit, ch); 
  }

  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) {
    make_jmp(multi_locants, multi_digit, ch); 
    make_jmp(multi_locants, multi_locants, ch); 
  }
  make_jmp(multi_locants, multi_locants, '&'); 
  make_jmp(multi_break, multi_locants, ' '); 

  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) 
    make_jmp(multi_size, multi_break, ch); 
  make_jmp(multi_size, multi_size, '&'); 
  
  make_jmp(ring_close, multi_size, 'J'); 

  /* hetero element definitions */
  const unsigned int ring_hetero        = state_init(); 
  const unsigned int hetero_dash_open   = state_init(); 
  const unsigned int hetero_dash_elem_a = state_init(); 
  const unsigned int hetero_dash_elem_b = state_init(); 
  const unsigned int hetero_dash_close  = state_init(); 
  const unsigned int hetero_space       = state_init(); 
  const unsigned int hetero_locant      = state_init(); 

  make_symbol_transitions(ring_hetero, ring_digits); 
  make_symbol_transitions(ring_hetero, digit_locant); 
  make_symbol_transitions(ring_hetero, multi_size); 
  make_symbol_transitions(ring_hetero, ring_hetero); 
  make_symbol_transitions(ring_hetero, hetero_locant); 
  make_symbol_transitions(ring_hetero, hetero_dash_close); 

  make_jmp(hetero_dash_open, ring_digits, '-'); 
  make_jmp(hetero_dash_open, digit_locant, '-'); 
  make_jmp(hetero_dash_open, multi_size, '-'); 
  make_jmp(hetero_dash_open, ring_hetero, '-'); 
  make_jmp(hetero_dash_open, hetero_locant, '-'); 

  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) {
    make_jmp(hetero_dash_elem_a, hetero_dash_open, ch); 
    make_jmp(hetero_dash_elem_b, hetero_dash_elem_a, ch); 
  }
  make_jmp(hetero_dash_close, hetero_dash_elem_b, '-'); 

  make_jmp(hetero_space, multi_size, ' '); 
  make_jmp(hetero_space, ring_hetero, ' '); 

  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) 
    make_jmp(hetero_locant, hetero_space, ch); 

  make_jmp(ring_close, hetero_dash_close, 'J'); 
  make_jmp(hetero_space, hetero_dash_close, ' '); 
  make_jmp(ring_close, ring_hetero, 'J'); 

  make_jmp(ring_close, ring_close, '&'); 
  
  /* aromaticity assignments */
  const unsigned int aromacity  = state_init(); 
  make_jmp(ring_close, aromacity, 'J'); 
  make_jmp(aromacity, aromacity, 'T'); 
  make_jmp(aromacity, aromacity, '&'); 
  
  make_jmp(aromacity, ring_digits, 'T'); 
  make_jmp(aromacity, ring_digits, '&'); 

  make_jmp(aromacity, digit_locant, 'T'); 
  make_jmp(aromacity, digit_locant, '&'); 

  make_jmp(aromacity, multi_size, 'T'); 
  make_jmp(aromacity, multi_size, '-'); 

  make_jmp(aromacity, ring_hetero, 'T'); 
  make_jmp(aromacity, ring_hetero, '&'); 

  make_jmp(aromacity, hetero_dash_close, 'T'); 
  make_jmp(aromacity, hetero_dash_close, '&'); 

  /* ring locant jump into specific blocks */ 
  const unsigned int rgroup_locant = state_init(); 

  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) 
    make_jmp(rgroup_locant, locant_space, ch); 
  make_jmp(rgroup_locant, rgroup_locant, '&'); 

  make_symbol_transitions(chain_open, rgroup_locant); 
  for (unsigned char ch='0'; ch <= '9'; ch++) 
    make_jmp(chain_open, rgroup_locant, ch); 
  make_jmp(dash_open, rgroup_locant, '-'); 
  
  make_jmp(locant_space, ring_close, ' '); 
  make_jmp(locant_space, rgroup_locant, ' '); 
  make_jmp(locant_space, chain_open, ' '); 

  /* ions */
  make_jmp(head, locant_space, '&'); 

  /* inline ring defintions */
  const unsigned int inline_ring_space  = state_init(); 
  const unsigned int inline_ring_locant = state_init(); 
  const unsigned int inline_spiro       = state_init(); 
  
  make_jmp(inline_ring_space, dash_open, ' '); 
  for (unsigned char ch = 'A'; ch <= 'Z'; ch++) 
    make_jmp(inline_ring_locant, inline_ring_space, ch); 
  make_jmp(inline_ring_locant, inline_ring_locant, '&'); 

  make_jmp(ring_open, inline_ring_locant, 'L');  
  make_jmp(ring_open, inline_ring_locant, 'T');  

  /* inline spiro definition */
  make_jmp(inline_spiro, dash_open, '&');
  make_jmp(inline_ring_space, inline_spiro, ' '); 

  fprintf(stderr, "%u states\n", nstates); 
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
static bool process_file(FILE *fp)
{
  const size_t bufsize = 16384; 
  const size_t wlnbuf_size = 4096; 

  char buffer[bufsize]; 
  char wln_buf[wlnbuf_size]; 
  unsigned int wlnlen = 0; 
  
  const unsigned short root_id = 0; 
  unsigned short curr_id = root_id; 
  
  unsigned int jmp_id; 
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

    if (curr_id)
      goto jmp_matching; 

jmp_not_matching:
    for (; i < bytes; i++) {
      unsigned char ch = buffer[i]; 
      jmp_id = state_jmp(curr_id, ch);  
      if (jmp_id != NO_JMP) {
        wlnlen    = 0; 
        match_end = 0; 
        curr_id = jmp_id; 
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
      jmp_id = state_jmp(curr_id, ch);  
      if (jmp_id == NO_JMP) {
        curr_id = root_id; 
        if (match_end > 0) 
          handle_match(wln_buf, match_end+1, wlnlen);  
        else if (opt_match_option == LINE_MATCH)
          fwrite(wln_buf, wlnlen, 1, stdout); 
        match_end = 0; 
        goto jmp_not_matching;
      }
      else if (wlnlen >= wlnbuf_size) {
        fprintf(stderr, "Warning: WLN string too long!\n");
        goto jmp_not_matching;
      }
      else {
        curr_id = jmp_id; 
        wln_buf[wlnlen] = ch; 
        match_end = is_accept(curr_id) ? wlnlen : match_end; 
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

  init_wlnfsm(); 

  process_file(ifp);  
  if (opt_match_option == COUNT_ONLY)
    printf("%lu matches\n", nmatches); 
  
  if (ifp != stdin) 
    fclose(ifp); 
  return 0;
}


