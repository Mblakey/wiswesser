

#include <stdlib.h>
#include <stdio.h>

#include "parser.h"

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
#include <openbabel/graphsym.h>

const char *cli_inp;
const char *format; 


static void DisplayUsage()
{
  fprintf(stderr, "readwln <options> -o<format> -s <input (escaped)>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, " -h                   show the help for executable usage\n");
  fprintf(stderr, " -o                   choose output format (-osmi, -oinchi, -ocan)\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser reads and evaluates wiswesser\n"
                  " line notation (wln), the parser is native\n"
                  " and will can return either a reformatted string*\n"
                  " *if rules do not parse exactly, and the connection\n"
                  " table which can be used in other libraries\n");
  DisplayUsage();
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i;

  cli_inp = (const char *)0;
  format = (const char *)0;

  if (argc < 2)
    DisplayUsage();

  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];

    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1])
      {

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
        else{
          fprintf(stderr,"Error: unrecognised format, choose between ['smi','inchi','can']\n");
          DisplayUsage();
        } 

      case 's':
        if(i+1 >= argc){
          fprintf(stderr,"Error: must add string after -s\n");
          DisplayUsage();
        }
        else{
          cli_inp = argv[i+1];
          i++;
        }
        break;

      default:
        fprintf(stderr, "Error: unrecognised input %s\n", ptr);
        DisplayUsage();
      }
    }
  }

  if(!format){
    fprintf(stderr,"Error: no output format selected\n");
    DisplayUsage();
  }

  if(!cli_inp){
    fprintf(stderr,"Error: no input string entered\n");
    DisplayUsage();
  }

  return;
}

int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);
  
  std::string res;
  OBMol mol;
  if(!ReadWLN(cli_inp,&mol))
    return 1;
  
  OBConversion conv;
  conv.AddOption("h",OBConversion::OUTOPTIONS);
  conv.SetOutFormat(format);

  res = conv.WriteString(&mol);
  std::cout << res;
  return 0;
}