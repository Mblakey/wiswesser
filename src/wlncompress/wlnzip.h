
#ifndef WLNZIP_H
#define WLNZIP_H

#include <stdlib.h>
#include <stdio.h>

#include <string>

#include "rfsm.h"

/* debugging bitstream struct, not for release application */
typedef struct  BitStream{
  unsigned char b; 
  struct BitStream *nxt; 
}BitStream;


void ReadStream(BitStream *stream);
void DeleteStream(BitStream *stream);
void Append(unsigned char b, BitStream *stream);

bool WLNdeflate(FILE *ifp, FSMAutomata *wlnmodel); 
bool WLNinflate(FILE *ifp, FSMAutomata *wlnmodel); 

BitStream* WLNPPMCompressBuffer(const char *str, FSMAutomata *wlnmodel, unsigned char escape_type); 
bool WLNPPMDecompressBuffer(BitStream *bitstream, FSMAutomata *wlnmodel, unsigned char escape_type); 


bool WLNPPMCompressFile(FILE *ifp, FSMAutomata *wlnmodel,   unsigned char escape_type); 
bool WLNPPMDecompressFile(FILE *ifp, FSMAutomata *wlnmodel, unsigned char escape_type); 


#endif
