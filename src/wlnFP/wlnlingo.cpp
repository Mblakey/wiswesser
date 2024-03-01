
#include <stdlib.h>
#include <string>
#include <set>
#include <vector>
#include <algorithm>

#include "fingerprint.h"

std::set<std::string> WLNLingo(const char *str, unsigned int len){
  std::set<std::string> lset; 
  std::string lingo_str; 
  for(unsigned int i=0;i<len-WLINGO;i++){
    for(unsigned int j=0;j<WLINGO;j++)
      lingo_str += str[i+j];
    
    lset.insert(lingo_str); 
    lingo_str.clear(); 
  }

  return lset;   
}




unsigned int Intersection(std::set<std::string> &v1, std::set<std::string> &v2){ 
  std::vector<std::string> v_intersection; 
  std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(),
                          std::back_inserter(v_intersection));

  return v_intersection.size(); 
}

unsigned int Union(std::set<std::string> &v1, std::set<std::string> &v2){ 
  v1.insert(v2.begin(),v2.end());
  return v1.size(); 
}





