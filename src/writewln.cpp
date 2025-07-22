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
#include "readline.h"

FILE *fp; 
const char *format; 

static void 
process_file(FILE *fp)
{
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;

  conv.SetInFormat(format);
  
  char wlnout[1024]; 
  char buffer[4096]; 
  
  while (readline(fp, buffer, 4096, 0)) {
    if (WriteWLN(wlnout,&mol))
      fprintf(stdout, "%s\n", wlnout); 
    else 
      fprintf(stdout, "NULL\n"); 
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
        fp = fopen(ptr, "r"); 
        if (!fp) {
          fprintf(stderr, "Error: could not open file at %s\n", ptr); 
          display_usage(); 
        }
        break;
    }
  }

  if (!format) {
    fprintf(stderr,"Error: no input format selected\n");
    display_usage();
  }

  if (!j) fp = stdin; 
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


