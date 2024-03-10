
#ifndef WLNZIP_H
#define WLNZIP_H

#include <stdlib.h>
#include <stdio.h>

#include <string>

#include "rfsm.h"

bool WLNdeflate(FILE *ifp, FSMAutomata *wlnmodel); 
bool WLNinflate(FILE *ifp, FSMAutomata *wlnmodel); 

bool WLNPPMCompressBuffer(const char *str, FSMAutomata *wlnmodel, std::string &bitstream, bool add_terminal); 
bool WLNPPMDecompressBuffer(std::string &bitstream, FSMAutomata *wlnmodel); 
bool WLNPPMCompressFile(FILE *ifp, FSMAutomata *wlnmodel, std::string &bitstream); 


#endif
