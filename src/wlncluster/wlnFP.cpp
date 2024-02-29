
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>

bool opt_verbose = false; 
const char *input; 


static void DisplayUsage(){
  fprintf(stderr, "wlnfp <string>\n");
  exit(1); 
}

static void ProcessCommandLine(int argc, char *argv[])
{

  const char *ptr = 0;
  int i,j;

  input = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        
        case 'h':
          DisplayUsage();
          break;

        case 'v':
          opt_verbose = true;
          break;

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
    else{
      switch(j++){
        case 0:
          input = ptr; 
          break;
            
        default:
          fprintf(stderr,"Error: multiple files not currently supported\n");
          exit(1);
      }
    }
  }

  if(!input){
    fprintf(stderr,"Error: no input given\n");
    DisplayUsage();
  }  
  return;
}


double WLNFingerprint(const char* string){

  
  // inorganics
  unsigned int Bsymbol = 0; // boron
  unsigned int Psymbol = 0; // phosphorous 
  unsigned int Ssymbol = 0; // sulphur 
  

  // nitrogens
  unsigned int Ksymbol = 0; // 4 nitro
  unsigned int Msymbol = 0; // 2 nitro
  unsigned int Nsymbol = 0; // 3 nitro 
  unsigned int Zsymbol = 0; // terminal nitro

  // carbons
  unsigned int Ysymbol =  0; // 3 carbonyl 
  unsigned int Xsymbol =  0; 

  // oxygens 
  unsigned int Osymbol = 0; 
  unsigned int Qsymbol = 0; 

  // halogens
  unsigned int Esymbol = 0; // terminal bromine
  unsigned int Fsymbol = 0; // terminal flourine 
  unsigned int Gsymbol = 0; // terminal chlorine 
  unsigned int Hsymbol = 0; // terminal hydrogen
  unsigned int Isymbol = 0; // terminal Iodide

  
  // Functional
  unsigned int Vsymbol = 0; // C=O
  unsigned int Wsymbol = 0; // X(=O)(=O)
  unsigned int Rsymbol = 0; // c1cccc1


  // patterns 
  unsigned int CarbonChains = 0; 
  unsigned int BondUnsaturations = 0; 
  unsigned AtomOther = 0; 

  // Cycles - Build this as we go  
  unsigned int Scaffolds = 0; // L/T - J notation 
  unsigned int Subcycles = 0; // Linear time ring perception 
    
  bool pending_locant = false;
  bool pending_J_closure = false;
  bool reading_chain = false; 
  bool reading_dash = false; 

  unsigned char ch = *string; 
  while(ch)
  { 
    switch (ch)
    {

    case '0': // cannot be lone, must be an addition to another num
      if(pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      break;

    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
      if(pending_J_closure)
          break; 
      else if (reading_dash)
        break;
      
      reading_chain = true; 
      break;
    

    case 'Y':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
       Ysymbol++; 
      break;

    case 'X':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }     
      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false; 
      else
        Xsymbol++; 
      break;


    case 'O':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Osymbol++; 
      break;

    case 'Q':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }     

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Qsymbol++;
      break;

    case 'V':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Vsymbol++; 
      break;

    case 'W':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Wsymbol++; 
      break;

      // nitrogens

    case 'N':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Nsymbol++; 
      break;

    case 'M':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)  
        pending_locant = false;
      else
        Msymbol++; 
      break;

    case 'K':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Ksymbol++; 
      break;

    case 'Z':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
       Zsymbol++;
      break;


    case 'E':
     if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Esymbol++; 
      break; 


    case 'G':
     if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Gsymbol++;  
      break;


    case 'F':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Fsymbol++; 
      break;

    case 'I':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break; 
      else if (pending_locant)
        pending_locant = false;
      else
        Isymbol++;
      
      break;

      // inorganics

    case 'B':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)  
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Bsymbol++;
      break;

    case 'P':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Psymbol++;

      break; 


    case 'S':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      
      else if (pending_locant)
        pending_locant = false;
      else
        Ssymbol++; 

      break;

    // multiply bonded carbon, therefore must be at least a double bond
    case 'C':
     if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      } 

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false; 
      
      // special carbon case 
      break;


    case 'A':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
          
      break;
        
    // this can start a chelating ring compound, so has the same block as 'L\T'
    case 'D':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
    
      break;
        
        
    // hydrogens explicit

    case 'H':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Hsymbol++; 
      break;


    case 'J':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (reading_dash)
        break;
      else if(pending_locant)
        pending_locant = false; 
      else if(pending_J_closure){
        pending_J_closure = false; 
        Scaffolds++; 
      }
      break;

    case 'L':
    case 'T':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
       pending_J_closure = true;      
      break;

    case 'R':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        Rsymbol++; 
      break;

      // bonding

    case 'U':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        BondUnsaturations++;

      // unsaturate here
      break;

      // specials


    case ' ':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if(reading_dash){
        reading_dash = false; 
        pending_locant = true; 
      }
      else if (pending_J_closure){
        break;
      }
      else 
        pending_locant = true; 
      
      break;


    case '&': // this could be ignored, unless doing number of branches as a fingerprint metric
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      
      break;


    case '-':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }

      if (pending_J_closure)
        break;

      if(reading_dash){
        reading_dash = false;
        AtomOther++; 
      }
      else
        reading_dash = true; 
      break;
    
    
    case '/':
      if(reading_chain){
        reading_chain = false; 
        CarbonChains++; 
      }
      break;

    default:
      return fprintf(stderr,"Error: unallowed character! - alphabet: [A-Z][0-1][&-/' ']");
      exit(1);
    }

    ch = *(++string); 
  }

  if(reading_chain)
    CarbonChains++; 

  if(opt_verbose){

    fprintf(stderr,"CarbonChains: %d\n", CarbonChains); 
    fprintf(stderr,"X symbols: %d\n", Xsymbol); 
    fprintf(stderr,"Y symbols: %d\n", Ysymbol);

    fprintf(stderr,"K symbols: %d\n", Ksymbol); 
    fprintf(stderr,"M symbols: %d\n", Msymbol); 
    fprintf(stderr,"N symbols: %d\n", Nsymbol); 
    
    fprintf(stderr,"O symbols: %d\n", Osymbol); 
    fprintf(stderr,"Q symbols: %d\n", Qsymbol); 
    
    fprintf(stderr,"P symbols: %d\n", Psymbol); 
    fprintf(stderr,"S symbols: %d\n", Ssymbol); 
    fprintf(stderr,"B symbols: %d\n", Bsymbol); 

    fprintf(stderr,"V symbols: %d\n", Vsymbol); 
    fprintf(stderr,"W symbols: %d\n", Wsymbol); 
    fprintf(stderr,"R symbols: %d\n", Rsymbol); 
 
    fprintf(stderr,"E symbols: %d\n", Esymbol); 
    fprintf(stderr,"F symbols: %d\n", Fsymbol); 
    fprintf(stderr,"G symbols: %d\n", Gsymbol); 
    fprintf(stderr,"H symbols: %d\n", Hsymbol); 
    fprintf(stderr,"I symbols: %d\n", Isymbol); 

    fprintf(stderr,"Unsaturations: %d\n", BondUnsaturations); 
    fprintf(stderr,"Other Atoms: %d\n", AtomOther); 
    
    fprintf(stderr,"Scaffolds: %d\n", Scaffolds); 
    fprintf(stderr,"Subcycles: %d\n", Subcycles); 
  }

  return 0; 
}



int main(int argc, char *argv[]){

  ProcessCommandLine(argc, argv);
  WLNFingerprint(input); 

  return 0; 
}
