
#ifndef WLNZIP_H
#define WLNZIP_H

#include <stdlib.h>
#include <stdio.h>

#include <string>

#include "rfsm.h"

typedef struct  BitStream{
  unsigned char b; 
  struct BitStream *nxt; 
}BitStream; 

bool WLNdeflate(FILE *ifp, FSMAutomata *wlnmodel); 
bool WLNinflate(FILE *ifp, FSMAutomata *wlnmodel); 

BitStream* WLNPPMCompressBuffer(const char *str, FSMAutomata *wlnmodel, unsigned char escape_type,bool add_terminal); 
bool WLNPPMDecompressBuffer(BitStream *bitstream, FSMAutomata *wlnmodel, unsigned char escape_type); 


bool WLNPPMCompressFile(FILE *ifp, FSMAutomata *wlnmodel,   unsigned char escape_type); 
bool WLNPPMDecompressFile(FILE *ifp, FSMAutomata *wlnmodel, unsigned char escape_type); 


#endif
