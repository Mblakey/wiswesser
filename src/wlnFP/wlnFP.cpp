#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fingerprint.h"
  
typedef struct{ 
 
// inorganics
  u_int16_t Bsymbol; // boron
  u_int16_t Psymbol; // phosphorous 
  u_int16_t Ssymbol; // sulphur 

// nitrogens
  u_int16_t Ksymbol; // 4 nitro
  u_int16_t Msymbol; // 2 nitro
  u_int16_t Nsymbol; // 3 nitro 
  u_int16_t Zsymbol; // terminal nitro

  // carbons
  u_int16_t Ysymbol; // 3 carbonyl 
  u_int16_t Xsymbol; 

  // oxygens 
  u_int16_t Osymbol; 
  u_int16_t Qsymbol; 

  // halogens
  u_int16_t Esymbol; // terminal bromine
  u_int16_t Fsymbol; // terminal flourine 
  u_int16_t Gsymbol; // terminal chlorine 
  u_int16_t Hsymbol; // terminal hydrogen
  u_int16_t Isymbol; // terminal Iodide


  // Functional
  u_int16_t Vsymbol; // C=O
  u_int16_t Wsymbol; // X(=O)(=O)
  u_int16_t Rsymbol; // c1cccc1

  // patterns 
  u_int16_t CarbonChains; 
  u_int16_t CarbonAtoms; // the total atoms counted by chains
  u_int16_t BondUnsaturations; 
  u_int16_t AtomOther; 

  // Cycles  
  u_int16_t ScaffoldAtoms; 
  u_int16_t HeteroScaffolds;      // L/T - J notation 
  u_int16_t CarbonScaffolds;      // L/T - J notation 

  u_int16_t AromSubcycles;  // Linear time ring perception 
  u_int16_t AlipSubcycles;  // Linear time ring perception 
  
  u_int16_t MultiCyclics;   // Linear time ring internals 
  u_int16_t BridgeAtoms;    // Linear time ring internals 
  u_int16_t SpiroPoints;    // Linear time ring internals 
  
}Descriptors;



void init_descriptors(Descriptors *desc){

  desc->CarbonAtoms = 0;
  desc->CarbonChains = 0;
  desc->Xsymbol = 0; 
  desc->Ysymbol = 0; 

  desc->Ksymbol = 0;
  desc->Nsymbol = 0; 
  desc->Msymbol = 0; 
  desc->Zsymbol = 0; 
  
  desc->Osymbol = 0;
  desc->Qsymbol = 0;
  desc->Vsymbol = 0; 
  desc->Wsymbol = 0; 

  desc->Bsymbol = 0;
  desc->Ssymbol = 0;
  desc->Psymbol = 0;
  desc->AtomOther = 0; 

  desc->Esymbol = 0; 
  desc->Fsymbol = 0; 
  desc->Gsymbol = 0; 
  desc->Hsymbol = 0; 
  desc->Isymbol = 0;


  desc->Rsymbol = 0;
  desc->BondUnsaturations = 0; 
  
  desc->ScaffoldAtoms = 0; 

  desc->HeteroScaffolds = 0;
  desc->CarbonScaffolds = 0;
  
  desc->AromSubcycles = 0; 
  desc->AlipSubcycles = 0;
  desc->MultiCyclics = 0; 
  desc->BridgeAtoms = 0;
  desc->SpiroPoints = 0; 
}


void debug_descriptors(Descriptors *desc){

  fprintf(stderr,"\nCarbon Atoms: %d\n", desc->CarbonAtoms); 
  fprintf(stderr,"Alkyl Chains: %d\n", desc->CarbonChains); 
  fprintf(stderr,"X symbols: %d\n", desc->Xsymbol); 
  fprintf(stderr,"Y symbols: %d\n", desc->Ysymbol);

  fprintf(stderr,"K symbols: %d\n", desc->Ksymbol); 
  fprintf(stderr,"M symbols: %d\n", desc->Msymbol); 
  fprintf(stderr,"N symbols: %d\n", desc->Nsymbol); 
  
  fprintf(stderr,"O symbols: %d\n", desc->Osymbol); 
  fprintf(stderr,"Q symbols: %d\n", desc->Qsymbol); 
  
  fprintf(stderr,"P symbols: %d\n", desc->Psymbol); 
  fprintf(stderr,"S symbols: %d\n", desc->Ssymbol); 
  fprintf(stderr,"B symbols: %d\n", desc->Bsymbol); 

  fprintf(stderr,"V symbols: %d\n", desc->Vsymbol); 
  fprintf(stderr,"W symbols: %d\n", desc->Wsymbol); 
  fprintf(stderr,"R symbols: %d\n", desc->Rsymbol); 

  fprintf(stderr,"E symbols: %d\n", desc->Esymbol); 
  fprintf(stderr,"F symbols: %d\n", desc->Fsymbol); 
  fprintf(stderr,"G symbols: %d\n", desc->Gsymbol); 
  fprintf(stderr,"H symbols: %d\n", desc->Hsymbol); 
  fprintf(stderr,"I symbols: %d\n", desc->Isymbol); 

  fprintf(stderr,"Unsaturations: %d\n", desc->BondUnsaturations); 
  fprintf(stderr,"Other Atoms:   %d\n", desc->AtomOther); 
  
  fprintf(stderr,"Scaffold Atoms:          %d\n", desc->ScaffoldAtoms); 
  fprintf(stderr,"Carbon Scaffolds:        %d\n", desc->CarbonScaffolds); 
  fprintf(stderr,"Hetero Scaffolds:        %d\n", desc->HeteroScaffolds); 
  fprintf(stderr,"Aliphatic Subcycles:     %d\n", desc->AlipSubcycles); 
  fprintf(stderr,"Aromatic  Subcycles:     %d\n", desc->AromSubcycles);
  fprintf(stderr,"Multicyclic Ring Points: %d\n", desc->MultiCyclics); 
  fprintf(stderr,"Spiro Points:            %d\n", desc->SpiroPoints);
  fprintf(stderr,"Ring Bridges:            %d\n", desc->BridgeAtoms);
}



unsigned int locant_to_int(unsigned char ch){
  if(ch < 'A'){
    fprintf(stderr,"Error: %c is not a valid locant\n",ch); 
    return 0; 
  }

  return (ch - 'A') +1; 
}



bool WLNRingParse(const char *cpy, unsigned int s, unsigned int e, Descriptors *desc){
  
  bool expecting_locant = false;
  bool expecting_size = false; 
  bool reading_dash = false; 
  unsigned int multi_skips = 0; 
  unsigned int pseudo_skips = 0; 

  unsigned int read_size = 0; 
  unsigned int total_cycles = 0; 
  unsigned int ali_cycles = 0; 
  unsigned int arom_cycles = 0; 

  unsigned int subcycles[64] = {0}; 
  unsigned char locant_read = 0; 
  
  char dash_capture[3] = {0};
  unsigned int dash_pos = 0; 


  for(unsigned int i=s;i<e;i++){
    unsigned char ch = cpy[i]; 

    if(multi_skips){
      multi_skips--;
      continue;
    }

    if(pseudo_skips){
      pseudo_skips--;
      continue;
    }

    switch(ch){
      
      case '0': 
      case '1':
      case '2':
      case '3':
      case '4': 
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        if(expecting_locant){
          desc->MultiCyclics += ch - '0';
          expecting_size = true; 
          expecting_locant = false;
          multi_skips = ch - '0';
        }
        else if(reading_dash){

          if(dash_pos > 2){
            fprintf(stderr,"Error: overflowing ring size buffer\n"); 
            return false; 
          }
          
          dash_capture[dash_pos++] = ch; 
        }
        else{
          subcycles[total_cycles++] = ch - '0'; 
          locant_read = 0;
        }
        break;
      
      case 'B':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch;
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch);
        }
        else {
          desc->Bsymbol++;
          locant_read = 0; 
        }
        break; 
      
      case 'E':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Esymbol++;
          locant_read = 0; 
        }
        break;

      case 'F':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Fsymbol++; 
          locant_read = 0;
        }
        break;

      case 'G':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Gsymbol++; 
          locant_read = 0; 
        }
        break;

      case 'H':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Hsymbol++; 
          locant_read = 0; 
        }
        break;

      case 'I':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Isymbol++; 
          locant_read = 0; 
        }
        break;

      case 'A':
      case 'C':
      case 'D':
      case 'J':
      case 'Z':
      case 'Q':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          fprintf(stderr,"Error: locant only character read as atom\n"); 
          return false;
        }
        break;
 

      case 'K':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Ksymbol++; 
          locant_read = 0; 
        }
        break;

     case 'M':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Msymbol++; 
          locant_read = 0; 
        }
        break;

     case 'N':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Nsymbol++; 
          locant_read = 0; 
        }
        break;

     case 'O':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Osymbol++; 
          locant_read = 0;
        }
        break;

      case 'P':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch);
        }
        else{
          desc->Psymbol++; 
          locant_read = 0; 
        }
        break;


     case 'S':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch); 
        }
        else{
          desc->Ssymbol++; 
          locant_read = 0; 
        }
        break;


      case 'L':
        if(reading_dash)
          break;
        else if(i==s)
          desc->CarbonScaffolds++;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size){
          read_size = locant_to_int(ch);
          locant_read = 0;
        }
        break;


      case 'T':
        if(reading_dash)
          break;
        else if(i==s)
          desc->HeteroScaffolds++;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size)
          read_size = locant_to_int(ch); 
        else{
          if(locant_read){
            desc->BridgeAtoms++;
            locant_read = 0; 
          }

          ali_cycles++;  
        }
        break;

      case 'U':
        if(reading_dash)
          break;
        else if (expecting_locant){
          expecting_locant = false;
          locant_read = ch; 
        }
        else if (expecting_size && !read_size)
          read_size = locant_to_int(ch); 
        else{
          desc->BondUnsaturations++; 
          locant_read = 0; 
        }
        break;

      case '&':
        if (expecting_size && read_size){
          read_size += 23;
          arom_cycles++; 
          locant_read = 0; 
        }
        break;

      case '-':
        if(expecting_size){
          expecting_size = false;
        }
        else if (reading_dash){
          if(dash_pos > 0){
            subcycles[total_cycles++] = atoi(dash_capture);
            memset(dash_capture,0,3); 
            dash_pos = 0; 
          }
          else {
            desc->AtomOther++; 
          }
          reading_dash = false; 
        }
        else {
          reading_dash = true; 
        }
        break;

      case ' ':
        if(reading_dash){
          reading_dash = false;
          i++; // sort of cheating but why not, skips the double bond assignment
          locant_read = 0;
        }
        if(expecting_size && !read_size){
          break;
        }
        else if (locant_read){
          desc->BridgeAtoms++; 
          locant_read = 0;
          expecting_locant = true; 
        }
        else{
          expecting_locant = true;
          expecting_size = false;
        }
        break;

      case '/':
        pseudo_skips = 2;
        break;

      default: 
        fprintf(stderr,"Error: unallowed character! - (%c) alphabet: [A-Z][0-1][&-/' ']\n", ch);
        return false;  
    }

  }
  

  // -- aromatic and aliphatic subcycles 
  if(ali_cycles + arom_cycles == total_cycles){
    desc->AromSubcycles += arom_cycles;
    desc->AlipSubcycles += ali_cycles; 
  }
  else if (ali_cycles == 1)
    desc->AlipSubcycles += total_cycles; 
  else if (arom_cycles == 0)
    desc->AromSubcycles += total_cycles; 
  else{
    fprintf(stderr,"Error: subcycles read do match the aromaticity characters within the scaffold\n"); 
    return false; 
  }
  
  // -- ring sizes
  if(read_size)
    desc->ScaffoldAtoms += read_size;
  else{
    // use calc
    read_size+= subcycles[0];
    for(unsigned int i=1;i<total_cycles;i++){
      read_size += (subcycles[i] - 2); 
    }
    desc->ScaffoldAtoms += read_size;
  }

  return true; 
}

bool WLNParse(const char* string, Descriptors *desc){
  
  bool pending_locant = false;
  bool pending_J_closure = false; // use as a ring flag 
  bool reading_chain = false; 
  bool reading_dash = false; 

  u_int16_t chain_pos = 0; // track carbon counts
  char chain[3] = {0}; 

  const char *cpy = string; 
  unsigned char ch = *string; 
  unsigned int p = 0; 
  unsigned int r_start = 0; 

  while(ch)
  { 
    switch (ch)
    {

    case '0': // cannot be lone, must be an addition to another num
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
      if(chain_pos<3)
        chain[chain_pos++] = ch; 
      else{
        fprintf(stderr,"Error: overflowing carbon chain\n");
        return false; 
      }
      break;
    
    case 'Y':
      if(reading_chain){
        reading_chain = false;
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
       desc->Ysymbol++; 
      break;

    case 'X':
      if(reading_chain){
        reading_chain = false;
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 

        desc->CarbonChains++; 
      }     
      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false; 
      else
        desc->Xsymbol++; 
      break;


    case 'O':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Osymbol++; 
      break;

    case 'Q':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }     

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Qsymbol++;
      break;

    case 'V':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Vsymbol++; 
      break;

    case 'W':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Wsymbol++; 
      break;

      // nitrogens

    case 'N':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Nsymbol++; 
      break;

    case 'M':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)  
        pending_locant = false;
      else
        desc->Msymbol++; 
      break;

    case 'K':
      if(reading_chain){
        reading_chain = false;
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Ksymbol++; 
      break;

    case 'Z':
      if(reading_chain){
        reading_chain = false;
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
       desc->Zsymbol++;
      break;


    case 'E':
     if(reading_chain){
        reading_chain = false;
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Esymbol++; 
      break; 


    case 'G':
     if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Gsymbol++;  
      break;


    case 'F':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Fsymbol++; 
      break;

    case 'I':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break; 
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Isymbol++;
      
      break;

      // inorganics

    case 'B':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)  
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Bsymbol++;
      break;

    case 'P':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Psymbol++;

      break; 


    case 'S':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Ssymbol++; 

      break;

    // multiply bonded carbon, therefore must be at least a double bond
    case 'C':
     if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
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
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
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
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
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
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Hsymbol++; 
      break;


    case 'J':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (reading_dash)
        break;
      else if(pending_locant)
        pending_locant = false; 
      else if(pending_J_closure){
        pending_J_closure = false; 
        if(!WLNRingParse(cpy, r_start, p,desc))
          return (uint16_t*)0;
        r_start = 0; 
      }
      break;

    case 'L':
    case 'T':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else{
        pending_J_closure = true;      
        r_start = p; 
      }
      break;

    case 'R':
      if(reading_chain){
        reading_chain = false;
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->Rsymbol++; 
      break;

      // bonding

    case 'U':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;
      else if (reading_dash)
        break;
      else if (pending_locant)
        pending_locant = false;
      else
        desc->BondUnsaturations++;

      // unsaturate here
      break;

      // specials


    case ' ':
      if(reading_chain){
        reading_chain = false;
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
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
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }
      
      if (pending_J_closure)
        break;
      else if (reading_dash){
        desc->SpiroPoints++;
        reading_dash = false;
      }
      break;


    case '-':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }

      if (pending_J_closure)
        break;

      if(reading_dash){
        reading_dash = false;
        desc->AtomOther++; 
      }
      else
        reading_dash = true; 
      break;
    
    
    case '/':
      if(reading_chain){
        reading_chain = false; 
        desc->CarbonAtoms += atoi(chain); 
        memset(chain, 0, 3);
        chain_pos = 0; 
        desc->CarbonChains++; 
      }
      break;

    default:
      fprintf(stderr,"Error: unallowed character! - alphabet: [A-Z][0-1][&-/' ']");
      return false; 
    }

    ch = *(++string);
    p++;
  }

  if(reading_chain){
    desc->CarbonChains++;
    desc->CarbonAtoms += atoi(chain); 
  }

  return true; 
}

u_int16_t *WLNFingerprint(const char *string){
  Descriptors *desc = (Descriptors*)malloc(sizeof(Descriptors));
  init_descriptors(desc);

  if(!WLNParse(string, desc)){
    free(desc); 
    return 0;
  }
  u_int16_t *FP = (u_int16_t*)malloc(sizeof(u_int16_t)* FPSIZE);
  memset(FP, 0, sizeof(u_int16_t) * FPSIZE); 

  // -- Carbons -- 
  FP[0] = desc->CarbonAtoms; 
  FP[1] = desc->CarbonChains;
  FP[2] = desc->Xsymbol;
  FP[3] = desc->Ysymbol; 

  // -- Nitrogens
  FP[4] = desc->Ksymbol;
  FP[5] = desc->Nsymbol; 
  FP[6] = desc->Msymbol; 
  FP[7] = desc->Zsymbol; 

  // -- Oxygens 
  FP[8] = desc->Osymbol; 
  FP[9] = desc->Qsymbol; 
  FP[10] = desc->Vsymbol;
  FP[11] = desc->Wsymbol; 

  // -- Inorganics
  FP[12] = desc->Bsymbol;
  FP[13] = desc->Ssymbol;
  FP[14] = desc->Psymbol; 
  FP[15] = desc->AtomOther; 

  // -- Halogens
  FP[16] = desc->Esymbol;
  FP[17] = desc->Fsymbol;
  FP[18] = desc->Gsymbol;
  FP[19] = desc->Hsymbol;
  FP[20] = desc->Isymbol;

  // -- Bonding
  FP[21] = desc->BondUnsaturations; 

  // -- Cycles
  FP[22] = desc->ScaffoldAtoms; 
  FP[23] = desc->CarbonScaffolds;
  FP[24] = desc->HeteroScaffolds; 
  FP[25] = desc->AromSubcycles; 
  FP[26] = desc->AlipSubcycles; 
  FP[27] = desc->MultiCyclics; 
  FP[28] = desc->BridgeAtoms; 
  free(desc); 
  return FP; 
}


bool WLNDescriptors(const char *string){
  Descriptors *desc = (Descriptors*)malloc(sizeof(Descriptors));
  init_descriptors(desc);
  if(!WLNParse(string, desc)){
    free(desc);
    return false;
  }
  else 
    debug_descriptors(desc);

  free(desc); 
  return true;
}

