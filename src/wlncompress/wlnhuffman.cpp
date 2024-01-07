#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <stack>
#include <vector>
#include <map>

#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"

#define RCODE_BYTES 2

unsigned int opt_mode = 0;
unsigned int opt_verbose = false;

const char *input;

/* priority queue and huffman tree code */
struct Node;
struct PQueue;

struct Node{
  unsigned int freq;
  unsigned char ch; 
  Node *l; 
  Node *r; 
  Node *p;
};

struct PQueue{
  Node **arr; 
  unsigned int cap;
  unsigned int size;  
}; 

Node *AllocateNode(unsigned char ch, unsigned int f){
  Node *n = (Node*)malloc(sizeof(Node)); 
  n->freq = f;
  n->ch = ch; 
  n->l = 0; 
  n->r = 0;
  n->p = 0;
  return n; 
}

bool init_heap(PQueue *heap, unsigned int cap){
  heap->arr = (Node**)malloc(sizeof(Node*) * cap); 
  if(!heap->arr){
    fprintf(stderr,"Error: could not allocate memory\n");
    return false;
  }

  for (unsigned int i=0;i<cap;i++)
    heap->arr[i] = 0;
  
  heap->size = 0; 
  heap->cap = cap;
  return true;   
}


void shift_right(unsigned int low,unsigned int high, PQueue *heap){
  unsigned int root = low; 
  while(2*root+1 <= high){
    unsigned int leftIdx = 2*root+1;
    unsigned int rightIdx = 2*root+2; 
    unsigned int swapIdx = root;  

    if (leftIdx <= high && heap->arr[leftIdx]->freq < heap->arr[swapIdx]->freq )
      swapIdx = leftIdx;
    if (rightIdx <= high && heap->arr[rightIdx]->freq < heap->arr[swapIdx]->freq)
      swapIdx = rightIdx;

    if(swapIdx != root){
      Node* tmp = heap->arr[root];
      heap->arr[root] = heap->arr[swapIdx];
      heap->arr[swapIdx] = tmp;
      root = swapIdx; 
    }
    else
      break;
  }
}

void heapify(unsigned int low, unsigned int high, PQueue *heap)
{
  if(!heap->size)
    return;

  int midIdx = 0;
  if(!high)
    midIdx = 0;
  else
    midIdx = (high - low -1)/2;

  while (midIdx >= 0)
  {
    shift_right(midIdx, high,heap);
    midIdx--;
  }
}

bool blind_insert(Node *term, PQueue *heap){
  if(heap->size+1 > heap->cap){
    fprintf(stderr,"Error: maxing array capacity\n");
    return false; 
  }

  heap->arr[heap->size++] = term;
  return true; 
}

bool insert_term(Node *term, PQueue *heap){
  if(!heap->size)
    return blind_insert(term, heap); 
  else{
    if(blind_insert(term, heap)){
      heapify(0,heap->size-1,heap);
      return true;
    }else
      return false; 
  }
}


Node *pop_front(PQueue *priority_queue){
  if(!priority_queue->size){
    fprintf(stderr,"Error: popping empty heap\n");
    return 0;
  }
  
  Node *f = priority_queue->arr[0];
  priority_queue->arr[0] = 0; 

  unsigned int pos = 0;
  for(unsigned int i=0;i<priority_queue->size;i++){
    Node *s = priority_queue->arr[i]; 
    if(s){
      priority_queue->arr[i] = NULL;
      priority_queue->arr[pos++] = s;
    }
  }

  priority_queue->size = pos;   
  heapify(0,pos-1,priority_queue);

  return f;
}


Node *ConstructHuffmanTree(PQueue *priority_queue){
  if(!priority_queue->size){
    fprintf(stderr,"Error: constructing huffman tree from empty queue\n");
    return 0;
  }
  else if(priority_queue->size == 1){
    Node *root = AllocateNode('*',0);
    Node *f = pop_front(priority_queue);
    f->p = root; 
    root->l = f;
    return root;
  }
  else{
    while(priority_queue->size > 1){
      // get the first 2 and remove
      Node *f = pop_front(priority_queue);
      Node *s = pop_front(priority_queue);

      Node *sum = AllocateNode('*',f->freq + s->freq);
      sum->l = f; 
      sum->r = s; 

      f->p = sum;
      s->p = sum;
      insert_term(sum,priority_queue);
    }

    Node *root = pop_front(priority_queue);
    if(root && root->ch != '*')
      fprintf(stderr,"Error: building huffman trees, root has ch - %c(%d)\n",root->ch,root->ch);
  
    return root; 
  }
}

void DeleteHuffmanTree(Node *root){
  if(!root)
    return;

  Node *top = 0; 
  std::stack<Node*> stack;
  stack.push(root); 
  while(!stack.empty()){
    top = stack.top();
    stack.pop();
    if(top->l)
      stack.push(top->l);
    if(top->r)
      stack.push(top->r);

    free(top);
    top = 0;
  }
}

/* builds the code in reverse and writes to stream */
bool WriteHuffmanCode(Node *root,unsigned char ch, std::string &cstream){
  
  Node *ch_node = 0;
  Node *top = 0;
  std::stack<Node*> stack; 
  stack.push(root);
  while(!stack.empty()){
    top = stack.top();
    stack.pop();

    if(top->ch == ch){
      ch_node = top;
      break;
    }

    if(top->l)
      stack.push(top->l);
    
    if(top->r)
      stack.push(top->r);
  }

  if(!ch_node){
    fprintf(stderr,"Error: could not find %c (%d) in states huffman tree\n",ch,ch);
    return false;
  }

  unsigned int clen = 0;
  char code[128] = {0};

  Node *curr = ch_node; 
  Node *prev = ch_node; 
  while(curr->p){
    prev = curr;
    curr = curr->p; 

    if(prev == curr->l)
      code[clen++] = 0;
    else if (prev == curr->r)
      code[clen++] = 1; 
  }

  for(int i=clen-1;i>=0;i--)
    cstream += code[i];

  return true;
}


unsigned int count_bytes(char *buffer, unsigned int n){
  for(unsigned int i=0;i<n;i++){
    if(!buffer[i])
      return i;
  }
  return n;
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


void debug_buffer(unsigned char*buffer){
  for(unsigned int i=0;i<BUFF_SIZE;i++){
    if(!buffer[i])
      break;
    fprintf(stderr,"%d,",buffer[i]);
  }

  fprintf(stderr,"\n");
}

/* huffman style */
void stream_to_bytes(std::string &stream){
  unsigned int char_pos = 0;
  unsigned char out = 0;
  for(unsigned int i=0;i<stream.size();i++){
    if(stream[i])
      out ^= (1 << (7-char_pos) );

    char_pos++;
    if(char_pos == 8){
      char_pos = 0;
      fwrite(&out,sizeof(unsigned char),1,stdout);
      out = 0;
    }

  }
}

// convert the table number into the characters needed in the fsm, using a 2-byte scheme
// as the minimum ring entry needed is 3 chars, does a lot of work in string manipulation
void uint_to_chars(unsigned int val, unsigned char *buffer){
  unsigned int char_pos = 0;
  unsigned int char_bit = 7;

  unsigned char ch = 0;
  for(int i=31;i>=0;i--){
    bool bit = val & (1 << i) ? 1:0;
    if(bit)
      ch ^= (1 << char_bit);
    
    char_bit--;
    if(i % 8==0 || i==0){
      if(ch)
        buffer[char_pos++] = ch;
      
      ch = 0;    
      char_bit = 7;
    }
  }
}

bool encode_file(FILE *ifp, FSMAutomata *wlnmodel, std::map<FSMState*,PQueue*> &queue_lookup){

  std::map<std::string,bool> seen_before;

  unsigned char ch = 0;
  unsigned int bytes_read = 0;

  Node *htree = 0;
  PQueue *priority_queue = 0;
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  std::string cstream;
  unsigned int stream_bits = 0;

  // get the ring close, little forced but works quickly
  FSMState *ring_close = 0;
  const char *find_close = "L6J";
  while( (*find_close) ){
    for(edge=curr->transitions;edge;edge=edge->nxt){
      if(edge->ch == *find_close)
        curr = edge->dwn; 
    }
    find_close++;
  }
  ring_close = curr; 


  FSMState *ring_dict_entry = 0;
   
  unsigned int ring_tsize = 256;        
  std::vector<std::string> ring_table; // this can get very big, so use ring_tsize to limit
  ring_table.reserve(UINT16_MAX); // -1 as we never want 0 0 0 0
  
  // need these to be global
  std::string ring_fragment;
  unsigned char ring_code[RCODE_BYTES] = {0}; // 2 byte codes + entery 'r'
  
  curr = wlnmodel->root;
  unsigned char buffer[BUFF_SIZE] = {0};
  while(ReadLineFromFile(ifp,(char*)buffer,BUFF_SIZE,true)){
    if(opt_verbose)
      bytes_read += count_bytes((char*)buffer,BUFF_SIZE);

    for(unsigned int i=0;i<BUFF_SIZE;i++){
      ch = buffer[i];
      if(!ch){
        memset(buffer,0,BUFF_SIZE);
        break;
      }

#define RING_DICT 1
#if RING_DICT 
      if(curr == wlnmodel->root && (ch == 'L' || ch == 'T')){
        // we either jump through a table entry already made, or we create a transition to be encoded
        
        // letter r for rings, 'f' for fragments
        if(!ring_dict_entry){
          ring_dict_entry = wlnmodel->AddState(false);
          edge = wlnmodel->AddTransition(wlnmodel->root,ring_dict_entry,'r'); 

          // add a priority queue pointer here
          if(!queue_lookup[ring_dict_entry]){
            PQueue *nqueue = (PQueue*)malloc(sizeof(PQueue)); 
            init_heap(nqueue,256); // safe value for alphabet 
            queue_lookup[ring_dict_entry] = nqueue;
          }

        }
        
        // yank out the ring fragment
        FSMState *rsearch = curr;
        FSMEdge *redge = 0; 

        unsigned int j = i;
        while(rsearch != ring_close){
          for(redge = rsearch->transitions;redge;redge=redge->nxt){
            if(redge->ch == buffer[j]){
              rsearch = redge->dwn;
              ring_fragment += redge->ch;
              j++;
              break;
            }
          }
        }

        // if found, replace section with dictionary char
        unsigned int tpos = 0;
        bool found = false; // just to keep the code readable.
        for(tpos=0;tpos<ring_tsize;tpos++){
          if(ring_fragment == ring_table[tpos]){
            found = true;    
            break;
          }
        }

        if(found){

          // now that the minimum len we get is 3 bytes, ('r' + code) 
          //we can simply overwrite, mark and sweep

          uint_to_chars(tpos,ring_code);
          for(unsigned int p=i;p<j;p++)
            buffer[p] = 0; // null while yanked section
          
          buffer[i] = 'r'; // r opener

          for(unsigned int k=0;k<RCODE_BYTES;k++) // write the code into the stream
            buffer[i+k+1] = ring_code[k];

          // mark sweep any nulls out of the string
          unsigned int idx = 0;
          for(unsigned int k=0;k<BUFF_SIZE;k++){
            unsigned char bchar = buffer[k];
            buffer[k] = 0;
            if(bchar)
              buffer[idx++] = bchar;
          }

          FSMState *cstate = ring_dict_entry; // pass the 'r', already made
          FSMEdge *cedge = 0;

          // have to use full 2 byte codes, to keep dfa property. 
          // * need to point the last one at the ring closure

          for(cedge=cstate->transitions;cedge;cedge=cedge->nxt){
            if(cedge->ch == ring_code[0]){
              cstate = cedge->dwn;
              break;
            }
          }

          // create the jump state
          if(cstate == ring_dict_entry){
            cstate = wlnmodel->AddState();
            cedge = wlnmodel->AddTransition(ring_dict_entry,cstate,ring_code[0]);

            PQueue *nqueue = (PQueue*)malloc(sizeof(PQueue)); 
            init_heap(nqueue,256); // safe value for alphabet 
            queue_lookup[cstate] = nqueue; 
          }
          
          // add the last code back to the ring closure. 
          cedge = wlnmodel->AddTransition(cstate,ring_close,ring_code[1]);
      
          ring_fragment.clear();
          ch = buffer[i]; // set the buffer to the first code letter, should be 'r' if code
        }
     
        // if not found means not in table, yet, 
        // decompresser needs this to transition once before adding        
      }
      else if((ring_tsize < UINT16_MAX) &&  curr == ring_close && ring_fragment.size() > 1){
        
        // if we dont create transitions here, any single seen fragments only take up table space, not transitions 
        // create the code in the FSM, post read, decompresser is lock step
        uint_to_chars(ring_tsize,ring_code);
        while(!ring_code[0] || !ring_code[1])
          uint_to_chars(ring_tsize++,ring_code);
        
        ring_table[ring_tsize++] = ring_fragment;
        ring_fragment.clear();
      }
#endif

      // construct tree based on C values
      priority_queue = queue_lookup[curr]; 
      for(edge=curr->transitions;edge;edge=edge->nxt){
        Node *n = AllocateNode(edge->ch,edge->c);
        insert_term(n,priority_queue);
      }

      htree = ConstructHuffmanTree(priority_queue);
      if(!htree){
        fprintf(stderr,"Huffman tree allocation fault\n");
        return false;
      }

      for(edge=curr->transitions;edge;edge=edge->nxt){
        if(edge->ch == ch){
          curr = edge->dwn;
          edge->c++;
          break;
        }
      }

      if(!WriteHuffmanCode(htree,ch,cstream))
        return false;
      

      DeleteHuffmanTree(htree); 

      // keep the memory in check, every 32 bytes clear and start again
      if(cstream.size() == 256){
        stream_to_bytes(cstream);
        cstream.clear();
        stream_bits+=256;
      }
    }
  }


  // // write a byte from the root of the machine indicating EOF
  priority_queue = queue_lookup[wlnmodel->root];
  for(edge=wlnmodel->root->transitions;edge;edge=edge->nxt){
    Node *n = AllocateNode(edge->ch,edge->c);
    insert_term(n,priority_queue);
  }

  htree = ConstructHuffmanTree(priority_queue);
  WriteHuffmanCode(htree,0,cstream);
  DeleteHuffmanTree(htree);

  stream_to_bytes(cstream);
  stream_bits+=cstream.size();

  if(opt_verbose){
    fprintf(stderr,"ring table size: %d\n",ring_tsize);
    fprintf(stderr,"%d to %d bits: %f compression ratio\n",
            bytes_read*8,stream_bits,
            (double)(bytes_read*8)/stream_bits);
  }

  return true;
}


bool decode_file(FILE *ifp, FSMAutomata *wlnmodel,std::map<FSMState*,PQueue*> &queue_lookup){
  
  Node *tree_root = 0;
  Node *htree = 0; 
  PQueue *priority_queue = 0;
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  unsigned char bit_char = 0;
  unsigned int bit_pos = 0;
  fread(&bit_char,sizeof(unsigned char),1,ifp); // read first bit

  for(;;){

    unsigned char ch_read = 0;
    priority_queue = queue_lookup[curr]; 
    for(edge=curr->transitions;edge;edge=edge->nxt){
      Node *n = AllocateNode(edge->ch,edge->c);
      insert_term(n,priority_queue);
    }

    tree_root = ConstructHuffmanTree(priority_queue);
    htree = tree_root;
    if(!htree){
      fprintf(stderr,"Huffman tree allocation fault\n");
      return false;
    }

    while(!ch_read){
      unsigned char bit = bit_char & (1 << (7 - bit_pos) ) ? 1:0;
      bit_pos++; 

      // fetch next char from stream
      if (bit_pos == 8){
        if(!fread(&bit_char,sizeof(unsigned char),1,ifp)){
          DeleteHuffmanTree(tree_root);
          return true;
        }
          
        bit_pos = 0;
      } 

      if(bit)
        htree = htree->r;
      else
        htree = htree->l;

      if(!htree){
        fprintf(stderr,"Error: dead traversal\n");
        return false;
      }
    
      if(htree->ch != '*'){
        ch_read = htree->ch; 

        if(!ch_read){
          if(curr == wlnmodel->root){
            DeleteHuffmanTree(tree_root);
            return true; 
          }
          else{
            fprintf(stderr,"Error: reading null byte not at fsm root\n");
            DeleteHuffmanTree(tree_root);
            return false;
          }
        }
        else{
          fputc(ch_read,stdout);
          for(edge=curr->transitions;edge;edge=edge->nxt){
            if(edge->ch == ch_read){
              curr = edge->dwn;
              edge->c++;
              break; 
            }
          }

        }
         // should free all?
        DeleteHuffmanTree(tree_root);
        tree_root = 0;
      }

     
    }
  }


  return true;
}


static void DisplayUsage()
{
  fprintf(stderr, "wlnhuffman <options> <input> > <out>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -c          compress input\n");
  fprintf(stderr, "  -d          decompress input\n");
  fprintf(stderr, "  -v          verbose debugging statements on\n");
  exit(1);
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

  FSMAutomata *wlnmodel = CreateWLNDFA(REASONABLE*2,REASONABLE*4); // build the model 

  // minic arithmetic 
  wlnmodel->AddTransition(wlnmodel->root,wlnmodel->root,'\0');  
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept)
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->root,'\n');
  }

  wlnmodel->AssignEqualProbs();

  std::map<FSMState*,PQueue*> queue_lookup;   
  // create a queue for each state, and adaptively build as we go. 
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    FSMState *s = wlnmodel->states[i];
    PQueue *priority_queue = (PQueue*)malloc(sizeof(PQueue)); 
    init_heap(priority_queue,256); // safe value for ASCII values
    queue_lookup[s] = priority_queue;
  }


  FILE *fp = 0; 
  fp = fopen(input,"rb");
  if(fp){
    if(opt_mode == 1)
      encode_file(fp,wlnmodel,queue_lookup);
    else if (opt_mode == 2)
      decode_file(fp,wlnmodel,queue_lookup);

    fclose(fp);
  }
  else{
    fprintf(stderr,"Error: could not open file at %s\n",input);
    return 1;
  }


  // safe with memory
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    FSMState *s = wlnmodel->states[i];
    PQueue *priority_queue = queue_lookup[s];
    if(priority_queue){
      while(priority_queue->size){
        free(priority_queue->arr[priority_queue->size-1]);
        priority_queue->size--;
      }

      free(priority_queue->arr);
      free(priority_queue);
      priority_queue = 0;
    }
  }

  delete wlnmodel;
  return 0;
}