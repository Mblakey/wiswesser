
#include <stdio.h>
#include <string.h>

#include "align.h"

unsigned int WLNAlignment(const char *a, const char *b){
  if (strcmp(a, b)) {
    AlignStruct stralign(strlen(a)+1, strlen(b)+1);
    stralign.Init();
    ResultStruct result(strlen(a),strlen(b)); 
    stralign.NeedlemanWuncsh(a,b,result);
  
    result.WriteInstructionsVerbose(stderr); 
    return result.num_changes; 
  }
  else
    return 0;
}

