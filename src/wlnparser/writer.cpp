#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

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

#include "parser.h"

using namespace OpenBabel; 

#define REASONABLE 1024

const char *cli_inp;
const char *format; 

bool opt_modern = false;

static void DisplayUsage()
{
  fprintf(stderr, "writewln <options> -i<format> <input (escaped)>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -h                    show the help for executable usage\n");
  fprintf(stderr, "  -i                    choose input format (-ismi, -iinchi, -ican, -imol)\n");
  fprintf(stderr, "  -m                    write mwln (modern) strings (part of michaels PhD work\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- wisswesser notation parser ---\n\n");
  fprintf(stderr, " This parser writes to wiswesser\n"
                  " line notation (wln) from smiles/inchi, the parser is built on OpenBabels\n"
                  " toolkit and will return the WLN string for the given input.\n"
                  " Note: Rule 2 canonicalisation is not implemented here, ordering of the input atoms\n"
                  "       will effect the resulting WLN, a canconicalisation algorithm to obey such\n"
                  "       rules is in development\n");
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
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){

        case 'h':
          DisplayHelp();

        case 'i':
          if (!strcmp(ptr, "-ismi"))
          {
            format = "smi";
            break;
          }
          else if (!strcmp(ptr, "-iinchi"))
          {
            format = "inchi";
            break;
          }
          else if (!strcmp(ptr, "-ican"))
          {
            format = "can";
            break;
          }
          else if (!strcmp(ptr, "-imol"))
          {
            format = "mol";
            break;
          }
          else{
            fprintf(stderr,"Error: unrecognised format, choose between ['smi','inchi','can','mol']\n");
            DisplayUsage();
          }

        case 'm':
          opt_modern = true;
          break;

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
    else{
      switch (j)
      {
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

  if(!format){
    fprintf(stderr,"Error: no input format selected\n");
    DisplayUsage();
  }

  if(!cli_inp){
    fprintf(stderr,"Error: no input string entered\n");
    DisplayUsage();
  }


  if(opt_modern)
    fprintf(stderr,"Warning: modern wln functions not fully functional\n");

  return;
}

int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);
  
  bool res;
  OBMol mol;
  OBConversion conv;

  conv.SetInFormat(format);  
  if(!strcmp(format,"mol"))
    res = conv.ReadFile(&mol, cli_inp);
  else
    res = conv.ReadString(&mol,cli_inp);

  if(!res){
    fprintf(stderr,"Error: failed to read incoming notation\n");
    return 1; 
  }

  std::string buffer;
  buffer.reserve(1000);
  if(!WriteWLN(buffer,&mol,opt_modern))
    return 1;
  
  std::cout << buffer << std::endl;
  return 0;
}
