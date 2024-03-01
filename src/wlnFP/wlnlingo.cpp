
#include <stdlib.h>
#include <map>
#include <string>

#include "fingerprint.h"

LingoTable *WLNLingo(const char *str, unsigned int len, unsigned int lingo){
  LingoTable *LTable = (LingoTable*)malloc(sizeof(LingoTable)); 
  std::map<std::string, unsigned int>  hash_map; 
  std::string lingo_str; 
  for(unsigned int i=0;i<len-lingo;i++){
    for(unsigned int j=0;j<lingo;j++)
      lingo_str += str[i+j];
    
    hash_map[lingo_str]++;
    lingo_str.clear(); 
  }

  LTable->ltable = (LingoEntry**)malloc(sizeof(LingoEntry) * hash_map.size()); 
  LTable->size = hash_map.size(); 
    
  unsigned int i = 0;
  for(std::map<std::string,unsigned int>::iterator miter = hash_map.begin(); miter != hash_map.end(); miter++){
    LTable->ltable[i] = (LingoEntry*)malloc(sizeof(LingoEntry));
    LTable->ltable[i]->str = miter->first; 
    LTable->ltable[i]->n = miter->second; 
    i++;
  }

  return LTable;  
}


/* as given in Rogers Paper */ 


