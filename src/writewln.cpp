#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef USING_OPENBABEL
#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/kekulize.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#endif

#include "wlnparser.h"

FILE *fp; 
const char *format; 


static unsigned char 
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
process_file(FILE *fp)
{
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;

  conv.SetInFormat(format);
  
  char wlnout[1024] = {0}; 

  char buffer[4096]; 
  while (readline(fp, buffer, 4096, 0)) {
    conv.ReadString(&mol, buffer); 
    if (WriteWLN(wlnout, 1024, &mol))
      fprintf(stdout, "%s\n", wlnout); 
    else 
      fprintf(stdout, "NULL\n"); 
    mol.Clear(); 
  }
}


static void 
display_usage()
{
  fprintf(stderr, "writewln <options> -i<format> <input>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -i                    choose input format (-ismi, -iinchi, -ican, -imol)\n");
  exit(1);
}


static void 
process_cml(int argc, char *argv[])
{
  const char *ptr = 0;
  int i;
  unsigned int j = 0;
  
  fp = NULL; 
  format = (const char *)0;

  for (i = 1; i < argc; i++) {
    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]) switch (ptr[1]) {
      case 'i': format = &ptr[2]; break; 
      default:
        fprintf(stderr, "Error: unrecognised input %s\n", ptr);
        display_usage();
    }
    else switch (j++) {
      case 0:
        if (ptr[0] == '-' && !ptr[1])
          fp = stdin;
        else {
          fp = fopen(ptr, "r"); 
          if (!fp) {
            fprintf(stderr, "Error: could not open file at %s\n", ptr); 
            display_usage(); 
          }
        }
        break;
    }
  }

  if (!format) {
    fprintf(stderr,"Error: no input format selected\n");
    display_usage();
  }

  if (!fp) 
    fp = stdin; 
  return;
}


int 
main(int argc, char *argv[])
{
  process_cml(argc, argv);
  process_file(fp);  
  if (fp != stdin) fclose(fp); 
  return 0;
}


