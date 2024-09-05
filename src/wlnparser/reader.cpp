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
#include <iostream>

#include "parser.h"

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/obconversion.h>
#include <openbabel/babelconfig.h>

const char *cli_inp;
const char *format; 
bool opt_old = false;

static void DisplayUsage()
{
  fprintf(stderr, "readwln <options> -o<format> <input (escaped)>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, " -h                   show the help for executable usage\n");
  fprintf(stderr, " -o                   choose output format (-osmi, -oinchi, -okey, -ocan, -owln *)\n");
  fprintf(stderr, "                      * selecting -owln will return the shortest possible wln string\n");
  fprintf(stderr, " --old                use the old wln parser (nextmove software)\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser reads and evaluates wiswesser\n"
                  " line notation (wln), the parser is native C\n"
                  " with a plug in function to OpenBabel\n"
        );
  fprintf(stderr, " Input is expected to either be a WLN string\n"
                  " or a file of WLN strings seperated with newline\n"
                  " characters. Detection is done by checking for a\n"
                  " file extension seperator . \n"
        );

  DisplayUsage();
}


static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i;
  unsigned int j = 0;

  cli_inp = (const char *)0;
  format = (const char *)0;

  if (argc < 2)
    DisplayUsage();

  for (i = 1; i < argc; i++)
  {
    ptr = argv[i];

    if (ptr[0] == '-' && ptr[1]) {

      if (ptr[1] >= 'A' && ptr[1] <= 'Z' && !j) {
        cli_inp = ptr;
        j++; 
      }
    
      else {
        switch (ptr[1]) {

        case 'h':
          DisplayHelp();

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
          else if (!strcmp(ptr, "-ocan"))
          {
            format = "can";
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
            fprintf(stderr,"Error: unrecognised format, choose between ['smi','inchi','can','key','wln']\n");
            DisplayUsage();
          } 
        
        case '-':
          if (!strcmp(ptr, "--old")) {
            opt_old = true;
            break;
          }
          
        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
        }
      }
    }
    else {
      switch (j) {
        case 0:
          cli_inp = ptr;
          break;
        
        default:
          fprintf(stderr,"Error: wln string already set - %s\n",cli_inp);
          DisplayUsage();
      }
      j++;
    }
  }

  if (!format) {
    fprintf(stderr,"Error: no output format selected\n");
    DisplayUsage();
  }

  if (!cli_inp) {
    fprintf(stderr,"Error: no input entered\n");
    DisplayUsage();
  }

  return;
}

int main(int argc, char *argv[])
{
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv;
  
  ProcessCommandLine(argc, argv);
  
  if(opt_old){
    if(!NMReadWLN(cli_inp,&mol))
      return 1;
  }
  else if(!strcmp(format, "WLN")){
    if(!CanonicaliseWLN(cli_inp,&mol))
      return 1;
  }
  else{
    if(!C_ReadWLN(cli_inp,&mol))
      return 1;
    conv.AddOption("h",OpenBabel::OBConversion::OUTOPTIONS);
    conv.SetOutFormat(format);
    conv.Write(&mol,&std::cout);
  }
  return 0;
}
