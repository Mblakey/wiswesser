#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <map>
#include <vector> // can optimise this out later
#include <string> // just for prototype


unsigned saved_bytes = 0;
unsigned int opt_mode = 0;
unsigned int opt_verbose = false;
const char *input;


void print_bits(unsigned char val) {
  for (int i = 7; i >= 0; i--)
    fprintf(stderr,"%d", (val & (1 << i)) ? 1:0);
  fprintf(stderr,"\n");
}

void write_6bits(unsigned char val, std::string &buffer) {
  for (int i = 7; i >= 2; i--)
    buffer += ((val & (1 << i)) ? '1':'0');
}


unsigned char * encode_string(  const char *wln, 
                                std::map<unsigned char,unsigned int> &encode)
{
  
  // calculate the number of padding bits we need. 
  unsigned int sbits = 8;
  unsigned int bits = 6;  // write the null character as a block of 6
  unsigned int padding = 0; 
  std::string bitstring; 

  unsigned int i=0;
  while(wln[i] != 0){
    if(!encode[wln[i]]){
      fprintf(stderr,"Error: character %c is not in the wln alphabet\n",wln[i]);
      return 0;
    }
    
    i++;
    bits += 6; 
    sbits+= 8;
  }

  while((bits + padding) % 8 != 0)
    padding++; 
  

  unsigned char *encoded_str = (unsigned char*)malloc(sizeof(unsigned char) * (bits+padding)/8);
  memset(encoded_str,0,(bits+padding)/8);

  // parse the character string and replace with encode number, left shift the bits by 2
  // and then read the first 6 bits to create a binary 'string' using chars

  i = 0;
  while(wln[i] != 0){
    unsigned int encoding = encode[wln[i++]];
    write_6bits(encoding << 2,bitstring);
  }

  write_6bits(0,bitstring); // write the null character and add padding bits
  for(unsigned j=0;j<padding;j++)
    bitstring += '0';

  unsigned int p = 0;
  for(unsigned int j=0;j<(bits+padding);j+=8){
    i = 0;
    char buffer[8] = {0}; 

    for(unsigned int k=j; k<j+8;k++)
      buffer[i++] = bitstring[k];
    
    encoded_str[p++] = (unsigned char)strtol(buffer,0,2);
  }

  saved_bytes += (sbits - (bits+padding))/8;
  return encoded_str;
} 


unsigned char* decode_string( const char *encoded_string, 
                              std::map<unsigned int,unsigned char> &decode){

  if(!encoded_string){
    fprintf(stderr,"Error: decoding null string\n");
    return 0;
  }
  
  // easier to reverse the process into the binary bitstring and read from there
  std::string bitstring; 
  unsigned int i=0;
  while(encoded_string[i] != 0){
    for (int k = 7; k >= 0; k--)
      bitstring += (encoded_string[i] & (1 << k)) ? '1':'0';
    i++;
  }

  unsigned int bits = bitstring.size();
  while(bits % 6 != 0)
    bits--; 

  unsigned char *decoded_str = (unsigned char*)malloc(sizeof(unsigned char) * (bits)/6);
  memset(decoded_str,0,(bits)/6);

  unsigned int d = 0;
  for(unsigned int i=0; i <bitstring.size();i+=6){
    unsigned int p = 2;
    char buffer[8] = {0};
    buffer[0] = '0';
    buffer[1] = '0'; // effectively reverse the leftshift

    for(unsigned int k=i;k<i+6;k++)
      buffer[p++] = bitstring[k];
      
    decoded_str[d++] = decode[(unsigned int)strtol(buffer,0,2)];
  }
  
  return decoded_str;
}


bool ReadLineFromFile(FILE *fp, char *buffer, unsigned int n, bool add_nl=true){
  char *end = buffer+n;
  char *ptr;
  int ch;

  ptr = buffer;
  do {
    ch = getc_unlocked(fp); // this increments fp
    if (ch == '\n') {
      if (add_nl)
        *ptr++ = '\n'; // if i want the newline or not
      *ptr = '\0';
      return true;
    }
    if (ch == '\f') {
      *ptr++ = '\n';
      *ptr = '\0';
      return true;
    }
    if (ch == '\r') {
      *ptr++ = '\n';
      *ptr = '\0';
      ch = getc_unlocked(fp);
      if (ch != '\n') {
        if (ch == -1)
          return false;
        ungetc(ch,fp);
      }
      return true;
    }
    if (ch == -1) {
      *ptr++ = '\n';
      *ptr = '\0';
      return ptr-buffer > 1;
    }
    *ptr++ = ch;
  } while (ptr < end);
  *ptr = 0;
  
  fprintf(stderr, "Warning: line too long!\n");
  return false;
}

static void DisplayUsage()
{
  fprintf(stderr, "compresswln <options> <input> > <out>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -c          compress input\n");
  fprintf(stderr, "  -d          decompress input\n");
  fprintf(stderr, "  -v          verbose debugging statements on\n");
  exit(1);
}

static void DisplayHelp()
{
  fprintf(stderr, "\n--- WLN Compression ---\n");
  fprintf(stderr, "This exec writes a wln file into a 6 bit representation, and can perform\n"
                  "various compression schemes which are selected in options.\n"
                  "This is part of michaels PhD investigations into compressing chemical strings.\n\n");
  DisplayUsage();
}

static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  input = (const char *)0;

  j = 0;
  for (i = 1; i < argc; i++)
  {

    ptr = argv[i];
    if (ptr[0] == '-' && ptr[1]){
      switch (ptr[1]){
        
        case 'c': 
          opt_mode = 1; 
          break;

        case 'd':
          opt_mode = 2; 
          break;

        case 'v':
          opt_verbose = true;
          break;

        case 'h':
          DisplayHelp();

        default:
          fprintf(stderr, "Error: unrecognised input %s\n", ptr);
          DisplayUsage();
      }
    }
    else{
      switch(j++){
        case 0:
          input = ptr; 
          break;
        default:
          fprintf(stderr,"Error: multiple files not currently supported\n");
          exit(1);
      }
    }
  }

  if(!input){
    fprintf(stderr,"Error: no input file given\n");
    DisplayUsage();
  }

  if(!opt_mode){
    fprintf(stderr,"Error: please choose -c or -d for file\n");
    DisplayUsage();
  }

  return;
}


int main(int argc, char *argv[])
{
  ProcessCommandLine(argc, argv);

  const char *wln = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 -/&";

  std::map<unsigned char,unsigned int> encode;
  std::map<unsigned int, unsigned char> decode;

  unsigned int j=1;
  for (unsigned int i=0;i<40;i++){
    encode[wln[i]] = j;
    decode[j] = wln[i];
    j++;
  }

  FILE *fp = 0; 
  fp = fopen(input,"rb");
  if(!fp){
    if(opt_mode == 1){
      unsigned char *encoded_str = encode_string(input,encode);
      if(!encoded_str)
        return 1;
      
      unsigned char *decoded_str = decode_string((const char*)encoded_str,decode);

      if(strcmp(input,(const char*)decoded_str) == 1)
        fprintf(stderr,"Error: decoding round trip failed - %s\n",decoded_str);
      else
        fprintf(stdout,"%s\n",encoded_str);
      
      free(encoded_str);
      free(decoded_str);
    }
    else if(opt_mode == 2)
      fprintf(stderr,"Error: decoding from terminal string input will ignore special character values\n");
  }
  else{
    fclose(fp);
  }

  if(opt_verbose)
    fprintf(stderr,"saved %d bytes\n",saved_bytes);

  return 0;
}