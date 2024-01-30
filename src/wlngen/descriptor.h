#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

#include <stdlib.h>
#include <stdio.h>


unsigned int HeavyAtomCount(const char*str){
  unsigned char ch = *str; 

  unsigned int hacount = 0;

  bool in_ring = false;
  unsigned char last_loc = 0;


  while(ch){


    if(!in_ring){
      switch (ch){
        case 'O':
        case 'P':
        case 'S':
        case 'F':
        case 'G':
        case 'I':
          hacount++;
          break;
      }
    }

    ch = *(++str);
  }


  return hacount;
}






#endif 