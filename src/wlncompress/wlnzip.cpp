#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include <vector>

#include "rfsm.h"
#include "huffman.h"
#include "lz.h"
#include "wlnzip.h"

#define CSIZE 64

// general functions
void LeftShift(unsigned char *arr, unsigned int len, unsigned int n)   
{
  memmove(&arr[0], &arr[n], (len-n)*sizeof(unsigned char));
  memset(&arr[len-n], 0, n * sizeof(unsigned char));
}

/* there will be some optimisations here, specially to get rid of vectors */
void stream_to_bytes(std::vector<unsigned char> &stream){
  unsigned int char_pos = 0;
  unsigned char out = 0;

  for(unsigned int i=0;i<stream.size();i++){
    if(stream[i])
      out ^= (1 << (7-char_pos) );

    char_pos++;
    if(char_pos == 8){
      char_pos = 0;
      fputc(out,stdout);
      out = 0;
    }
  }

  // always get the last one
  if(char_pos && out)
    fputc(out,stdout);
}

/* ######################################################################################### */


void FlushMachine(FSMAutomata *wlnfsm){
  for(unsigned int i=0;i<wlnfsm->num_edges;i++)
    wlnfsm->edges[i]->c = 1;
}

/* run the virtual FSM and score how many bits a backreference will save, take the highest */
unsigned int ScoreBackReference(  unsigned int length, unsigned int distance, 
                                  FSMState*curr, LLBucket **buckets)
{
  LLBucket *lb = 0;
  LLBucket *db = 0;

  for(unsigned int i=0;i<LZBUCKETS;i++){
    if(!lb && length <= buckets[i]->lstart)
      lb = buckets[i];

    if(!db && distance <= buckets[i]->dstart)
      db = buckets[i];
  }

  if(!lb)
    lb = buckets[LZBUCKETS-1];
  if(!db)
    db = buckets[LZBUCKETS-1];
  
  // the bits saved would be the codes + offset, vs if just huffman coding them normally
  // might be a good way to approximate this

  // approximate that each code will be 6 in length for fast calcs 
  unsigned int approx_lz = (lb->lbits + db->dbits) + 8 + 8;  
  unsigned int approx_normal = length * 4;
  
  int saved = approx_normal - approx_lz;
  if(saved < 0)
    return 0;
  else 
    return saved;
}


bool ParseFilePtr(FILE *ifp, FSMAutomata *wlnmodel){
  unsigned char ch = 0;
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;
  
  while(fread(&ch, sizeof(unsigned char), 1, ifp)){
    for(edge=curr->transitions;edge;edge=edge->nxt){
      if(edge->ch == ch){
        if(edge->c < UINT32_MAX)
          edge->c++; 
        curr = edge->dwn; 
        break;
      }
    }
  }

  rewind(ifp); // set the file pointer back to normal
  return true;
}



bool WLNENCODE(FILE *ifp, FSMAutomata *wlnmodel){

  unsigned char *buffer = (unsigned char*)malloc(sizeof(unsigned char)*BUFFSIZE);   
  memset(buffer,0,BUFFSIZE);

  bool reading_data = true;
  unsigned char ch = 0;
  unsigned int molecules = 0;

  Node *htree = 0;
  PQueue *priority_queue = (PQueue*)malloc(sizeof(PQueue)); 
  init_heap(priority_queue,512); // safe value for WLN


  // set up the distance tree
  for(unsigned int c = 0;c<LZBUCKETS;c++){
    Node *n = AllocateNode(('a'+c),1); // use 'a' as offset, only adds 1-2 KB overhead
    insert_term(n,priority_queue);
  }

  Node *lz_tree = ConstructHuffmanTree(priority_queue);
  LLBucket **buckets = init_buckets();

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  unsigned char code[CSIZE] = {0};
  std::vector<unsigned char> bitstream;


  // PARSE FILE FOR SMALL TREES
 
  wlnmodel->AssignEqualProbs();
  // ParseFilePtr(ifp, wlnmodel);
  
  // fill up the forward window
  for(unsigned int i=BACKREFERENCE;i<BUFFSIZE;i++){
    if(!fread(&ch,sizeof(unsigned char),1,ifp))
      reading_data = false;
    else
      buffer[i] = ch;
    
    if(!reading_data)
      break;
  } 

  if(!buffer[BACKREFERENCE]){
    fprintf(stderr,"Error: no data!\n");
    free(buffer);
    return false;
  }


  while(buffer[BACKREFERENCE]){
    unsigned int distance = 0;
    unsigned int length   = 0;
    unsigned int best_length = 0;
    unsigned int best_distance = 0;
    unsigned int saved = 0; // potential for overflow here, another optimisation point

    for(unsigned int i=0;i<BUFFSIZE;i++){
      if(i >= BACKREFERENCE && !length)
        break;
      
      if(buffer[i] == buffer[BACKREFERENCE+length]){
        length++;
        if(!distance)
          distance = BACKREFERENCE - i;
      }
      else if(length > 2){
        
        // if good backreference, store the positions in best, this is an approximation
        unsigned int lsaved = ScoreBackReference(length,distance,curr,buckets);
        if(lsaved > saved){
          best_distance = distance;
          best_length = length;
          saved = lsaved; 
        }

        distance = 0;
        length = 0;
      }
      else{
        // reset
        distance = 0;
        length = 0;
      }
    }

    if(best_length && best_distance){
      // here i need to encode, best length, bits little endian, best distance, bits little endian
      // move the FSM through the shifts and increase the adaptive count for those edges

      // potential increase the adaptive count per special symbol - optimisation once working
      bitstream.push_back(0);
      bitstream.push_back(0);

      LLBucket *lb = length_bucket(best_length,buckets);
      LLBucket *db = distance_bucket(best_distance,buckets);
      
      unsigned int loffset = best_length - lb->lstart;
      unsigned int doffset = best_distance - db->dstart;
      
      // write code for length bucket
      unsigned int clen = WriteHuffmanCode(lz_tree,lb->symbol,code);
      if(!clen){
        fprintf(stderr,"Error: length code generation\n");
        return false;
      }

      for(unsigned int c = 0;c<clen;c++)
        bitstream.push_back(code[c]);
      memset(code,0,CSIZE);

      // write the bits expected for the length offset in little endian order 
      for(unsigned int j=0;j<lb->lbits;j++){
        if( loffset & (1 << j) )
          bitstream.push_back(1); 
        else
          bitstream.push_back(0);
      }

      // repeat for the distance symbol
      clen = WriteHuffmanCode(lz_tree,db->symbol,code);
      if(!clen){
        fprintf(stderr,"Error: distance code generation\n");
        return false;
      }
        

      for(unsigned int c = 0;c<clen;c++)
        bitstream.push_back(code[c]);
      memset(code,0,CSIZE);

      // write the bits expected for the length offset in little endian order 
      for(unsigned int j=0;j<db->dbits;j++){
        if( doffset & (1 << j) )
          bitstream.push_back(1); 
        else
          bitstream.push_back(0);
      }

      // shift everything down and move the machine in lockstep
      for(unsigned int j=0;j<best_length;j++){
        for(edge=curr->transitions;edge;edge=edge->nxt){
          if(edge->ch == buffer[BACKREFERENCE]){
            curr = edge->dwn;
            break;
          }
        }

        if(buffer[BACKREFERENCE] == '\n')
          molecules++;        

        LeftShift(buffer,BUFFSIZE,1);
        if(reading_data){
          if(!fread(&ch,sizeof(unsigned char),1,ifp))
            reading_data = 0;
          else
            buffer[BUFFSIZE-1] = ch;
        }
      }
    }
    else{
     
      for(edge=curr->transitions;edge;edge=edge->nxt){
        if(edge->c > 0){
          Node *n = AllocateNode(edge->ch,edge->c);
          insert_term(n,priority_queue);
        }
      }

      htree = ConstructHuffmanTree(priority_queue);
      ReserveCode("00",'*',htree);
      if(!htree){
        fprintf(stderr,"Huffman tree allocation fault\n");
        return false;
      }

      if(priority_queue->size){
        fprintf(stderr,"Error: queue is not being fully dumped on tree creation\n");
        return false;
      }


      unsigned int clen = WriteHuffmanCode(htree,buffer[BACKREFERENCE],code);
      if(!clen){
        fprintf(stderr,"Error: literal code generation - line %d\n",molecules);
        return false;
      }

      for(unsigned int c = 0;c<clen;c++)
        bitstream.push_back(code[c]);
      
      memset(code,0,CSIZE);

      for(edge=curr->transitions;edge;edge=edge->nxt){
        if(edge->ch == buffer[BACKREFERENCE]){
          curr = edge->dwn;
          edge->c++;
          break;
        }
      }

      if(buffer[BACKREFERENCE] == '\n')
        molecules++;

      LeftShift(buffer,BUFFSIZE,1);
      if(reading_data){
        if(!fread(&ch,sizeof(unsigned char),1,ifp))
          reading_data = 0;
        else
          buffer[BUFFSIZE-1] = ch;
      }


      free_huffmantree(htree); 

    }

    // gets deleted no matter what as well
  }


  // get the last ones still in the stream
  stream_to_bytes(bitstream);

  free_huffmantree(lz_tree);
  free_buckets(buckets);
  free_heap(priority_queue);
  free(buffer);

  fprintf(stderr,"%d molecules compressed\n",molecules);
  return true;
}


bool WLNDECODE(FILE *ifp, FSMAutomata *wlnmodel){
  
  unsigned char *buffer = (unsigned char*)malloc(sizeof(unsigned char)*BUFFSIZE); 
  memset(buffer,0,BUFFSIZE);

  unsigned char ch = 0; 
  unsigned int reading_bits = 0;
  unsigned int length = 0;
  unsigned int distance = 0;

  unsigned int offset = 0;
  unsigned int offpos = 0;
  
  Node *htree = 0; 
  Node *tree_root = 0;

  PQueue *priority_queue = (PQueue*)malloc(sizeof(PQueue)); 
  init_heap(priority_queue,512); // safe value for WLN

  // set up the distance tree
  for(unsigned int c = 0;c<LZBUCKETS;c++){
    Node *n = AllocateNode(('a'+c),1); // use 'a' as offset, only adds 1-2 KB overhead
    insert_term(n,priority_queue);
  }

  Node *lz_tree = ConstructHuffmanTree(priority_queue);
  LLBucket **buckets = init_buckets();

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  // init the first huffman tree for state, change on switch now
  for(edge=curr->transitions;edge;edge=edge->nxt){
    Node *n = AllocateNode(edge->ch,edge->c);
    insert_term(n,priority_queue);
  }

  tree_root = ConstructHuffmanTree(priority_queue);
  ReserveCode("00",'*',tree_root);
  htree = tree_root;


  wlnmodel->AssignEqualProbs();

  while(fread(&ch,sizeof(unsigned char),1,ifp)){
    for(int i=7;i>=0;i--){
      unsigned char bit = (ch & (1 << i))?1:0;


      if(reading_bits){

        if(bit)
          offset += (1 << offpos);

        reading_bits--;
        offpos++;

        if(!reading_bits){

          if(!distance)
            length += offset;
          else{
            distance += offset; 

            for(unsigned int i=0;i<length;i++)
              buffer[i+BACKREFERENCE] = buffer[BACKREFERENCE-distance+i];
             
            while(buffer[BACKREFERENCE]){
              for(edge=curr->transitions;edge;edge=edge->nxt){
                if(edge->ch == buffer[BACKREFERENCE]){
                  curr = edge->dwn;
                  //edge->c++;
                  break;
                }
              }

              fputc(buffer[BACKREFERENCE],stdout);
              LeftShift(buffer,BUFFSIZE,1);
            }

    
            // create the new queue
            for(edge=curr->transitions;edge;edge=edge->nxt){
              Node *n = AllocateNode(edge->ch,edge->c);
              insert_term(n,priority_queue);
            }

            // create the new tree
            tree_root = ConstructHuffmanTree(priority_queue);
            ReserveCode("00",'*',tree_root);
            htree = tree_root;

            if(!htree){
              fprintf(stderr,"Huffman tree allocation fault\n");
              return false;
            }

            distance = 0;
            length = 0;
          }

          offpos = 0;
          offset = 0;
        }
      }
      else{

        if(bit)
          htree = htree->r;
        else
          htree = htree->l;

        if(!htree){
          fprintf(stderr,"Error: dead traversal\n");
          return false;
        }
        else if(htree->ch == '*'){
          // tree swap marker
          htree = lz_tree;
        }
        else if(htree->ch >= 'a' && htree->ch < ('a'+LZBUCKETS)){
          
          if(!length){
            length = buckets[htree->ch - 'a']->lstart;
            reading_bits = buckets[htree->ch - 'a']->lbits;
            htree = lz_tree;
          }
          else if (!distance){
            distance = buckets[htree->ch - 'a']->dstart;
            reading_bits = buckets[htree->ch - 'a']->dbits;

            if(!reading_bits){

              for(unsigned int i=0;i<length;i++)
                buffer[i+BACKREFERENCE] = buffer[BACKREFERENCE-distance+i];
             
              while(buffer[BACKREFERENCE]){
                for(edge=curr->transitions;edge;edge=edge->nxt){
                  if(edge->ch == buffer[BACKREFERENCE]){
                    curr = edge->dwn;
                    //edge->c++;
                    break;
                  }
                }

                fputc(buffer[BACKREFERENCE],stdout);
                LeftShift(buffer,BUFFSIZE,1);
              }

  
              // create the new queue
              for(edge=curr->transitions;edge;edge=edge->nxt){
                Node *n = AllocateNode(edge->ch,edge->c);
                insert_term(n,priority_queue);
              }

              // create the new tree
              tree_root = ConstructHuffmanTree(priority_queue);
              ReserveCode("00",'*',tree_root);
              htree = tree_root;

              if(!htree){
                fprintf(stderr,"Huffman tree allocation fault\n");
                return false;
              }

              length = 0;
              distance = 0;
              offset = 0;
              offpos = 0;
            }
          }
        }
        else if (htree->ch){

          LeftShift(buffer,BACKREFERENCE,1);
          fputc(htree->ch,stdout);
          buffer[BACKREFERENCE-1] = htree->ch;

          // move the machine to the next state
          for(edge=curr->transitions;edge;edge=edge->nxt){
            if(edge->ch == htree->ch){
              curr = edge->dwn;
              edge->c++;
              break; 
            }
          }

          // delete the huffman tree 
          free_huffmantree(tree_root);

          // create the new queue
          for(edge=curr->transitions;edge;edge=edge->nxt){
            Node *n = AllocateNode(edge->ch,edge->c);
            insert_term(n,priority_queue);
          }

          // create the new tree
          tree_root = ConstructHuffmanTree(priority_queue);
          ReserveCode("00",'*',tree_root);
          htree = tree_root;

          if(!htree){
            fprintf(stderr,"Huffman tree allocation fault\n");
            return false;
          }

        }

      }
    }

  }
  
  free_heap(priority_queue);
  free(buffer);
  return true;
}

/* this function will NOT reset the edge->c values to zero on start */
unsigned int EncodedBits(const char*str, FSMAutomata *wlnmodel){

  unsigned char *buffer = (unsigned char*)malloc(sizeof(unsigned char)*BUFFSIZE); 
  memset(buffer,0,BUFFSIZE);
  
  bool reading_data = true;
  unsigned char ch = *str;

  Node *htree = 0;
  PQueue *priority_queue = (PQueue*)malloc(sizeof(PQueue)); 
  init_heap(priority_queue,512); // safe value for WLN


  // set up the distance tree
  for(unsigned int c = 0;c<LZBUCKETS;c++){
    Node *n = AllocateNode(('a'+c),1); // use 'a' as offset, only adds 1-2 KB overhead
    insert_term(n,priority_queue);
  }

  Node *lz_tree = ConstructHuffmanTree(priority_queue);
  LLBucket **buckets = init_buckets();

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  unsigned char code[CSIZE] = {0};
  std::vector<unsigned char> bitstream;
  
  // fill up the forward window
  for(unsigned int i=BACKREFERENCE;i<BUFFSIZE;i++){
    if(!ch)
      reading_data = false;
    else
      buffer[i] = ch;
    
    if(!reading_data)
      break;
    
    ch = *(++str);
  } 

  if(!buffer[BACKREFERENCE]){
    fprintf(stderr,"Error: no data!\n");
    free(buffer);
    return false;
  }


  while(buffer[BACKREFERENCE]){
    unsigned int distance = 0;
    unsigned int length   = 0;
    unsigned int best_length = 0;
    unsigned int best_distance = 0;
    unsigned int saved = 0; // potential for overflow here, another optimisation point

    for(unsigned int i=0;i<BUFFSIZE;i++){
      if(i >= BACKREFERENCE && !length)
        break;
      
      if(buffer[i] == buffer[BACKREFERENCE+length]){
        length++;
        if(!distance)
          distance = BACKREFERENCE - i;
      }
      else if(length > 2){
        
        // if good backreference, store the positions in best, this is an approximation
        unsigned int lsaved = ScoreBackReference(length,distance,curr,buckets);
        if(lsaved > saved){
          best_distance = distance;
          best_length = length;
          saved = lsaved; 
        }

        distance = 0;
        length = 0;
      }
      else{
        // reset
        distance = 0;
        length = 0;
      }
    }


    if(best_length && best_distance){
      // here i need to encode, best length, bits little endian, best distance, bits little endian
      // move the FSM through the shifts and increase the adaptive count for those edges
      
      // potential increase the adaptive count per special symbol - optimisation once working
      bitstream.push_back(0);
      bitstream.push_back(0);

      LLBucket *lb = length_bucket(best_length,buckets);
      LLBucket *db = distance_bucket(best_distance,buckets);
      
      unsigned int loffset = best_length - lb->lstart;
      unsigned int doffset = best_distance - db->dstart;
      
      // write code for length bucket
      unsigned int clen = WriteHuffmanCode(lz_tree,lb->symbol,code);
      if(!clen){
        fprintf(stderr,"Error: length code generation\n");
        return false;
      }

      for(unsigned int c = 0;c<clen;c++)
        bitstream.push_back(code[c]);
      memset(code,0,CSIZE);

      // write the bits expected for the length offset in little endian order 
      for(unsigned int j=0;j<lb->lbits;j++){
        if( loffset & (1 << j) )
          bitstream.push_back(1); 
        else
          bitstream.push_back(0);
      }

      // repeat for the distance symbol
      clen = WriteHuffmanCode(lz_tree,db->symbol,code);
      if(!clen){
        fprintf(stderr,"Error: distance code generation\n");
        return false;
      }
        

      for(unsigned int c = 0;c<clen;c++)
        bitstream.push_back(code[c]);
      memset(code,0,CSIZE);

      // write the bits expected for the length offset in little endian order 
      for(unsigned int j=0;j<db->dbits;j++){
        if( doffset & (1 << j) )
          bitstream.push_back(1); 
        else
          bitstream.push_back(0);
      }

      // shift everything down and move the machine in lockstep
      for(unsigned int j=0;j<best_length;j++){
        for(edge=curr->transitions;edge;edge=edge->nxt){
          if(edge->ch == buffer[BACKREFERENCE]){
            curr = edge->dwn;
            break;
          }
        }   

        LeftShift(buffer,BUFFSIZE,1);
        if(reading_data){
          if(!ch)
            reading_data = 0;
          else
            buffer[BUFFSIZE-1] = ch;

          ch = *(++str);
        }
      }
    }
    else{
      
      /* huffman tree gets made no matter what */
      for(edge=curr->transitions;edge;edge=edge->nxt){
        Node *n = AllocateNode(edge->ch,edge->c);
        insert_term(n,priority_queue);
      }

      htree = ConstructHuffmanTree(priority_queue);
      ReserveCode("00",'*',htree);
      if(!htree){
        fprintf(stderr,"Huffman tree allocation fault\n");
        return false;
      }

      if(priority_queue->size){
        fprintf(stderr,"Error: queue is not being fully dumped on tree creation\n");
        return false;
      }


      unsigned int clen = WriteHuffmanCode(htree,buffer[BACKREFERENCE],code);
      if(!clen)
        return false;
      
      for(unsigned int c = 0;c<clen;c++)
        bitstream.push_back(code[c]);
      
      memset(code,0,CSIZE);

      for(edge=curr->transitions;edge;edge=edge->nxt){
        if(edge->ch == buffer[BACKREFERENCE]){
          curr = edge->dwn;
          edge->c++;
          break;
        }
      }

      LeftShift(buffer,BUFFSIZE,1);
      if(reading_data){
        if(!ch)
          reading_data = 0;
        else
          buffer[BUFFSIZE-1] = ch;
        ch = *(++str);
      }

      free_huffmantree(htree); 
    }

  }
 

  // get the last ones still in the stream

  free_huffmantree(lz_tree);
  free_buckets(buckets);
  free_heap(priority_queue);
  free(buffer);

  return bitstream.size();
}


/* ######################################################################################### */
