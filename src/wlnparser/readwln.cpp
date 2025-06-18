/*********************************************************************
Author : Michael Blakey

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdbool.h>

#include <iostream>

#include "wlnparser.h"

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>

const char *cli_inp;
const char *format; 
bool opt_string; 


static void 
display_usage()
{
  fprintf(stderr, "--- wisswesser notation parser ---\n");
  fprintf(stderr, "This parser reads and evaluates wiswesser\n"
                  "line notation (wln), the parser is native C\n"
                  "with a plug in function to OpenBabel\n\n");
  fprintf(stderr, "readwln <options> -o<format> [-s <input (escaped)> | <infile>]\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, " -h                   show the help for executable usage\n");
  fprintf(stderr, " -o                   choose output format (-osmi, -oinchi, -okey, -ocan)\n");
  exit(1);
}


static void 
process_cml(int argc, char *argv[])
{
  const char *ptr = 0;
  int i;
  int j = 0;

  cli_inp = (const char *)0;
  format = (const char *)0;
  opt_string = false; 

  if (argc < 2)
    display_usage();

  for (i = 1; i < argc; i++) {
    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]) {
      switch (ptr[1]) {
        case 'h':
          display_usage(); 

        case 's':
          if (i >= argc - 1) {
            fprintf(stderr, "Error: -s flag must be followed by an escaped argument\n");
            display_usage();
          }
          cli_inp = argv[++i];  
          j = 1; 
          opt_string = true; 
          break; 

        case 'o':
          if (!strcmp(ptr, "-osmi"))
          {
            format = "smi";
            break;
          }
          else if (!strcmp(ptr, "-oinchi"))
          {
            format = "inchi";
            break;
          }
          else if(!strcmp(ptr,"-okey"))
          {
            format  = "inchikey";
            break;
          }
          else if(!strcmp(ptr,"-owln"))
          {
            format  = "WLN";
            break;
          }
          else {
            fprintf(stderr,"Error: unrecognised format, choose between ['smi','inchi','key']\n");
            display_usage();
          } 
        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          display_usage();
        }
    }
    else switch (j) {
      case 0:
        cli_inp = ptr;
        break;
      default:
        fprintf(stderr,"Error: reader input already set - %s\n", cli_inp);
        display_usage();
    }
    j++;
  }

  if (!format) {
    fprintf(stderr,"Error: no output format selected\n");
    display_usage();
  }
  if (!cli_inp) {
    fprintf(stderr,"Error: no input entered\n");
    display_usage();
  }
}

int 
main(int argc, char *argv[])
{
  process_cml(argc, argv);

  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;
  
  conv.SetOutFormat(format);
  conv.AddOption("h",OpenBabel::OBConversion::OUTOPTIONS);
  
  if (opt_string) {
    if (!ReadWLN(cli_inp, &mol))
      return 1;
    else
      conv.Write(&mol ,&std::cout);
  }
  else {
    FILE *fp = fopen(cli_inp, "r"); 
    if (!fp) {
      fprintf(stderr,"Error: file could not be opened\n"); 
      return 1;
    }
    else {
      if (!ReadWLNFile(fp, &mol, &conv)) {
        fclose(fp); 
        return 1; 
      }
      fclose(fp); 
    }
  }
  return 0;
}


