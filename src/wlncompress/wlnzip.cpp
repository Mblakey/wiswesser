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
    Node *root = AllocateNode(0,0);
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

      Node *sum = AllocateNode(0,f->freq + s->freq);
      sum->l = f; 
      sum->r = s; 

      f->p = sum;
      s->p = sum;
      insert_term(sum,priority_queue);
    }

    Node *root = pop_front(priority_queue);
    if(root && root->ch != 0)
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
      fputc(out,stdout);
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

  unsigned int look_ahead = 256;
  unsigned int look_behind = 1024;  

  unsigned char ahead[BUFF_SIZE] = {0};
  unsigned char behind[BUFF_SIZE] = {0};
  unsigned char store[look_ahead] = {0}; // handle data overlaps
  
  // read in big chunks for a good look ahead
  if(!fread(ahead,sizeof(unsigned char),BUFF_SIZE,ifp))
    fprintf(stderr,"Error: no data to fill buffer\n");

  // keep this constant with 

 


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

  unsigned int buff_pos = 0;
  unsigned char buffer[BUFF_SIZE] = {0};
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
          fprintf(stdout,"%s",buffer);
          buff_pos = 0;
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
    
      if(htree->ch != 0){
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
          buffer[buff_pos++] = ch_read;
          if(ch_read == '\n'){
            fprintf(stdout,"%s",buffer);
            memset(buffer,0,BUFF_SIZE);
            buff_pos = 0;
          }

          // transition
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