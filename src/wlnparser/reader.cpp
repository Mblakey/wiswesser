

#include <cstring>
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

      if(ptr[1] >= 'A' && ptr[1] <= 'Z' && !j){
        cli_inp = ptr;
        j++; 
      }
    
      else{
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
          else{
            fprintf(stderr,"Error: unrecognised format, choose between ['smi','inchi','can','key','wln']\n");
            DisplayUsage();
          } 
        
        case '-':
          if(!strcmp(ptr, "--old")){
            opt_old = true;
            break;
          }
          

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
        }
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

  if(opt_old){
    if(!NMReadWLN(cli_inp,&mol))
      return 1;
  }
  else if(!strcmp(format, "WLN")){
    if(!CanonicaliseWLN(cli_inp,&mol))
      return 1;
  }
  else{
    if(!ReadWLN(cli_inp,&mol))
      return 1;
    OBConversion conv;
    conv.AddOption("h",OBConversion::OUTOPTIONS);
    conv.SetOutFormat(format);

    res = conv.WriteString(&mol);
    std::cout << res;
  }
  return 0;
}
