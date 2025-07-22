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
#include <stdbool.h>

#include <iostream>

#include "wlnparser.h"

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>

#define WLN_MAX_INPUT 64

unsigned int ninput;
const char *input[WLN_MAX_INPUT];
const char *format; 

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
display_usage()
{
  fprintf(stderr, "--- wisswesser notation parser ---\n");
  fprintf(stderr, "This parser reads and evaluates wiswesser\n"
                  "line notation (wln), the parser is C\n"
                  "with a C++ plug in function to OpenBabel\n\n");
  fprintf(stderr, "readwln <options> -o<format> [<input (escaped)> | <infile>]\n");
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
  
  ninput = 0; 
  format = (const char *)0;

  if (argc < 2)
    display_usage();

  for (i = 1; i < argc; i++) {
    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]) {
      switch (ptr[1]) {
        case 'h':
          display_usage(); 
        case 'o': format = ptr+2; break;
        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          display_usage();
        }
    }
    else {
      if (ninput == WLN_MAX_INPUT) {
        fprintf(stderr, "Error: max input value reached - %u\n", WLN_MAX_INPUT); 
        exit(1); 
      }
      else 
        input[ninput++] = ptr; 
    } 
  }

  if (!format) {
    fprintf(stderr,"Error: no output format selected\n");
    display_usage();
  }

  if (!ninput) // stdin default
    input[ninput++] = "-";  
}


int 
main(int argc, char *argv[])
{
  process_cml(argc, argv);

  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;
  
  conv.SetOutFormat(format);
  conv.AddOption("h",OpenBabel::OBConversion::OUTOPTIONS);
  
  char buffer[1024]; 
  for (int i = 0; i < ninput; i++) {
    const char *inp = input[i]; 
    if (strcmp(inp, "-") == 0) while (readline(stdin, buffer, 1024, 0)) {
      if (ReadWLN(buffer, &mol))
        conv.Write(&mol ,&std::cout);
      else
       fprintf(stdout, "NULL\n");
      mol.Clear(); 
    }
    else {
      FILE *fp = fopen(inp, "r"); 
      if (fp) { 
        while (readline(fp, buffer, 1024, 0)) {
          if (ReadWLN(buffer, &mol))
            conv.Write(&mol ,&std::cout);
          else
           fprintf(stdout, "NULL\n");
          mol.Clear(); 
        }
        fclose(fp); 
      }
    }
  }
  return 0;
}


