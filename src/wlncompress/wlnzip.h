#ifndef WLNZIP_H
#define WLNZIP_H

#include <stdio.h>

#include "rfsm.h"
#include "lz.h"

bool WLNENCODE(FILE *ifp, FSMAutomata *wlnmodel);
bool WLNDECODE(FILE *ifp, FSMAutomata *wlnmodel);
unsigned int EncodedBits(const char*str, FSMAutomata *wlnmodel);


void LeftShift(unsigned char *arr, unsigned int len, unsigned int n); 
void stream_to_bytes(std::vector<unsigned char> &stream);

void FlushMachine(FSMAutomata *wlnfsm);
unsigned int ScoreBackReference(  unsigned int length, unsigned int distance, 
                                  FSMState*curr, LLBucket **buckets);


#endif
