#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <type_traits>

#include "fingerprint.h"
  
typedef struct{ 
 
// inorganics
  u_int8_t Bsymbol; // boron
  u_int8_t Psymbol; // phosphorous 
  u_int8_t Ssymbol; // sulphur 

// nitrogens
  u_int8_t Ksymbol; // 4 nitro
  u_int8_t Msymbol; // 2 nitro
  u_int8_t Nsymbol; // 3 nitro 
  u_int8_t Zsymbol; // terminal nitro

  // carbons
  u_int8_t Ysymbol; // 3 carbonyl 
  u_int8_t Xsymbol; 

  // oxygens 
  u_int8_t Osymbol; 
  u_int8_t Qsymbol; 

  // halogens
  u_int8_t Esymbol; // terminal bromine
  u_int8_t Fsymbol; // terminal flourine 
  u_int8_t Gsymbol; // terminal chlorine 
  u_int8_t Hsymbol; // terminal hydrogen
  u_int8_t Isymbol; // terminal Iodide


  // Functional
  u_int8_t Vsymbol; // C=O
  u_int8_t Wsymbol; // X(=O)(=O)
  u_int8_t Rsymbol; // c1cccc1

  // patterns 
  u_int8_t CarbonChains; 
  u_int8_t CarbonAtoms; // the total atoms counted by chains
  u_int8_t BondUnsaturations; 
  u_int8_t AtomOther; 

  // Cycles  
  u_int8_t RingAtoms; 
  u_int8_t HeteroScaffolds;      // L/T - J notation 
  u_int8_t CarbonScaffolds;      // L/T - J notation 

  u_int8_t Arom3cycles;
  u_int8_t Arom4cycles; 
  u_int8_t Arom5cycles; 
  u_int8_t Arom6cycles; 
  u_int8_t Arom7cycles; 
  u_int8_t Arom8cycles;
  u_int8_t Arom9cycles;
  u_int8_t AromBigCycle;
 
  u_int8_t Alip3cycles;
  u_int8_t Alip4cycles; 
  u_int8_t Alip5cycles; 
  u_int8_t Alip6cycles; 
  u_int8_t Alip7cycles; 
  u_int8_t Alip8cycles;
  u_int8_t Alip9cycles;
  u_int8_t AlipBigCycle;

  u_int8_t MultiCyclics;   // Linear time ring internals 
  u_int8_t BridgeAtoms;    // Linear time ring internals 
  u_int8_t SpiroPoints;    // Linear time ring internals 
  
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
  
  desc->RingAtoms = 0; 

  desc->HeteroScaffolds = 0;
  desc->CarbonScaffolds = 0;
  
  desc->Arom3cycles = 0;
  desc->Arom4cycles = 0; 
  desc->Arom5cycles = 0; 
  desc->Arom6cycles = 0; 
  desc->Arom7cycles = 0; 
  desc->Arom8cycles = 0;
  desc->Arom9cycles = 0;
  desc->AromBigCycle = 0;
 
  desc->Alip3cycles = 0;
  desc->Alip4cycles = 0; 
  desc->Alip5cycles = 0; 
  desc->Alip6cycles = 0; 
  desc->Alip7cycles = 0; 
  desc->Alip8cycles = 0;
  desc->Alip9cycles = 0;
  desc->AlipBigCycle = 0;

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
  
  fprintf(stderr,"Ring Atoms:             %d\n", desc->RingAtoms); 
  fprintf(stderr,"Carbon Scaffolds:       %d\n", desc->CarbonScaffolds); 
  fprintf(stderr,"Hetero Scaffolds:       %d\n", desc->HeteroScaffolds); 

  
  fprintf(stderr,"Arom3:                  %d\n", desc->Arom3cycles); 
  fprintf(stderr,"Arom4:                  %d\n", desc->Arom4cycles); 
  fprintf(stderr,"Arom5:                  %d\n", desc->Arom5cycles); 
  fprintf(stderr,"Arom6:                  %d\n", desc->Arom6cycles); 
  fprintf(stderr,"Arom7:                  %d\n", desc->Arom7cycles); 
  fprintf(stderr,"Arom8:                  %d\n", desc->Arom8cycles); 
  fprintf(stderr,"Arom9:                  %d\n", desc->Arom9cycles); 
  fprintf(stderr,"AromBig:                %d\n", desc->AromBigCycle); 

  fprintf(stderr,"Alip3:                  %d\n", desc->Alip3cycles); 
  fprintf(stderr,"Alip4:                  %d\n", desc->Alip4cycles); 
  fprintf(stderr,"Alip5:                  %d\n", desc->Alip5cycles); 
  fprintf(stderr,"Alip6:                  %d\n", desc->Alip6cycles); 
  fprintf(stderr,"Alip7:                  %d\n", desc->Alip7cycles); 
  fprintf(stderr,"Alip8:                  %d\n", desc->Alip8cycles); 
  fprintf(stderr,"Alip9:                  %d\n", desc->Alip9cycles); 
  fprintf(stderr,"AlipBig:                %d\n", desc->AlipBigCycle); 


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
  unsigned int cycles_read = 0; 
  
  bool subcycle_type[64] = {0}; 
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
          subcycle_type[cycles_read++] = false; 
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
          locant_read = 0; 
          expecting_size = false;
        }
        else
         subcycle_type[cycles_read++] = true;  
        
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
    
  if(total_cycles != cycles_read){
    
    if(cycles_read == 0){
      for(unsigned int i=0;i<total_cycles;i++)
        subcycle_type[i] = true;
    }
    else if (cycles_read == 1 && subcycle_type[0] == false){
      for(unsigned int i=0;i<total_cycles;i++)
        subcycle_type[i] = false;
    }
    else{
      fprintf(stderr,"Error: aromaticity assignments do no match the number of rings\n"); 
      return false; 
    }
  }
  
  // -- aromatic and aliphatic subcycles 
  for(unsigned int i=0;i<total_cycles;i++){
    unsigned int subcycle_size = subcycles[i];
    bool arom = subcycle_type[i];
    switch (subcycle_size) { 
      case 3:
        if(arom)
          desc->Arom3cycles++;
        else
         desc->Alip3cycles++;
        break;

      case 4:
        if(arom)
          desc->Arom4cycles++;
        else
         desc->Alip4cycles++;
        break;

      case 5:
        if(arom)
          desc->Arom5cycles++;
        else
         desc->Alip5cycles++;
        break;

      case 6:
        if(arom)
          desc->Arom6cycles++;
        else
         desc->Alip6cycles++;
        break;

      case 7:
        if(arom)
          desc->Arom7cycles++;
        else
         desc->Alip7cycles++;
        break;

      case 8:
        if(arom)
          desc->Arom8cycles++;
        else
         desc->Alip8cycles++;
        break;

      case 9:
        if(arom)
          desc->Arom9cycles++;
        else
         desc->Alip9cycles++;
        break;


      default:
        if(arom)
          desc->AromBigCycle++;
        else
         desc->AlipBigCycle++;
        break; 
    }
  }

  
  // -- ring sizes
  if(read_size)
    desc->RingAtoms += read_size;
  else{
    // use calc
    read_size+= subcycles[0];
    for(unsigned int i=1;i<total_cycles;i++){
      read_size += (subcycles[i] - 2); 
    }
    desc->RingAtoms += read_size;
  }

  return true; 
}

bool WLNParse(const char* string, Descriptors *desc){
  
  bool pending_locant = false;
  bool pending_J_closure = false; // use as a ring flag 
  bool reading_chain = false; 
  bool reading_dash = false;
  bool dash_numerical = false; 

  u_int8_t chain_pos = 0; // track carbon counts
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
      else if (reading_dash){
        dash_numerical = true; 
        break;
      }
      
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
          return (u_int8_t*)0;
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

      if (pending_J_closure){
        if(p - r_start == 1)
          pending_J_closure = false;
        break;
      }
      else if(reading_dash && !dash_numerical){
        reading_dash = false;
        dash_numerical = false; 
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

u_int8_t *WLNFingerprint(const char *string){
  Descriptors *desc = (Descriptors*)malloc(sizeof(Descriptors));
  init_descriptors(desc);

  if(!WLNParse(string, desc)){
    free(desc); 
    return 0;
  }
  u_int8_t *FP = (u_int8_t*)malloc(sizeof(u_int8_t)* FPSIZE);
  memset(FP, 0, sizeof(u_int8_t) * FPSIZE); 

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
  FP[22] = desc->RingAtoms; 
  FP[23] = desc->CarbonScaffolds;
  FP[24] = desc->HeteroScaffolds; 
  
  FP[25] = desc->Arom3cycles; 
  FP[26] = desc->Arom4cycles; 
  FP[27] = desc->Arom5cycles; 
  FP[28] = desc->Arom6cycles; 
  FP[29] = desc->Arom8cycles; 
  FP[30] = desc->Arom9cycles; 
  FP[31] = desc->AromBigCycle; 


  FP[32] = desc->Alip3cycles; 
  FP[33] = desc->Alip4cycles; 
  FP[34] = desc->Alip5cycles; 
  FP[35] = desc->Alip6cycles; 
  FP[36] = desc->Alip7cycles; 
  FP[37] = desc->Alip8cycles; 
  FP[38] = desc->Alip9cycles; 
  FP[39] = desc->AlipBigCycle; 

  FP[40] = desc->MultiCyclics; 
  FP[41] = desc->BridgeAtoms; 
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

