#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <stack>
#include <vector>

#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"

#define CSIZE 64
#define LZBUCKETS 30

unsigned int opt_mode = 0;
unsigned int opt_verbose = false;

const unsigned int window = 258;
const unsigned int backreference = 32768;
const unsigned int buff_size = window+backreference;

const char *input;

/* ######################################################################################### */

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


void delete_heap(PQueue *priority_queue){
  while(priority_queue->size){
    free(priority_queue->arr[priority_queue->size-1]);
    priority_queue->size--;
  }
    
  free(priority_queue->arr);
  free(priority_queue);
  priority_queue = 0;
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
unsigned int WriteHuffmanCode(Node *root,unsigned char ch, unsigned char *code){
  
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
    return 0;
  }

  unsigned int clen = 0;

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

  // reverse the code
  if(clen > 1){
    unsigned char *ptr1, *ptr2, temp; 
    ptr1 = code; 
    ptr2 = code+(clen-1); // point to last element off set

    while(ptr1 < ptr2){
      temp = *ptr1;
      *ptr1 = *ptr2; 
      *ptr2 = temp;
      
      ptr1++;
      ptr2--;
    }
  }  

  return clen;
}

/* ######################################################################################### */


/* ######################################################################################### */

// again not the most efficient but means we wont get lost
typedef struct{
  unsigned int lstart;
  unsigned int dstart; 

  unsigned int lbits;
  unsigned int dbits;
} LLBucket; 



/* we can keep the original DEFLATE specification for distances, but for lengths we choose our own */
LLBucket ** init_buckets(){
  LLBucket **buckets = (LLBucket**)malloc(sizeof(LLBucket*)*LZBUCKETS); // use for instant look up
  memset(buckets,0,sizeof(LLBucket*));

  for(unsigned int i=0;i<LZBUCKETS;i++)
    buckets[i] = (LLBucket*)malloc(sizeof(LLBucket));

  buckets[0]->dstart = 1;
  buckets[0]->dbits = 0;

  buckets[1]->dstart = 2;
  buckets[1]->dbits = 0;

  buckets[2]->dstart = 3;
  buckets[2]->dbits = 0;

  buckets[3]->dstart = 4;
  buckets[3]->dbits = 0;

  buckets[4]->dstart = 5;
  buckets[4]->dbits = 1;

  buckets[5]->dstart = 7;
  buckets[5]->dbits = 1;

  buckets[6]->dstart = 9;
  buckets[6]->dbits = 2;
  
  buckets[7]->dstart = 13;
  buckets[7]->dbits = 2;

  buckets[8]->dstart = 17;
  buckets[8]->dbits = 3;

  buckets[9]->dstart = 25;
  buckets[9]->dbits = 3;


  buckets[10]->dstart = 33;
  buckets[11]->dstart = 49;
  buckets[12]->dstart = 65;
  buckets[13]->dstart = 97;
  buckets[14]->dstart = 129;
  buckets[15]->dstart = 193;
  buckets[16]->dstart = 257;
  buckets[17]->dstart = 385;
  buckets[18]->dstart = 513;
  buckets[19]->dstart = 769;
  buckets[20]->dstart = 1025;
  buckets[21]->dstart = 1537;
  buckets[22]->dstart = 2049;
  buckets[23]->dstart = 3073;
  buckets[24]->dstart = 4097;
  buckets[25]->dstart = 6145;
  buckets[26]->dstart = 8193;
  buckets[27]->dstart = 12289;
  buckets[28]->dstart = 16385;
  buckets[29]->dstart = 24577;

  return buckets;
}

void PurgeBuckets(LLBucket **buckets){
  for(unsigned int i=0;i<LZBUCKETS;i++){
    if(buckets[i])
      free(buckets[i]);
  }
  free(buckets);
}


/* ######################################################################################### */



/* ######################################################################################### */


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
  unsigned int bits = 0;

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


/* run the virtual FSM and score how many bits a backreference will save, take the highest */
unsigned int ScoreBackReference(unsigned int length, unsigned int distance, FSMState*curr){




}

bool WLNENCODE(FILE *ifp, FSMAutomata *wlnmodel){

  unsigned char *buffer = (unsigned char*)malloc(sizeof(unsigned char)*buff_size); 
  memset(buffer,0,buff_size);

  bool reading_data = true;
  bool reference = false;
  unsigned char ch = 0;
  unsigned int bytes_read = 0;
  unsigned int fpos = backreference;
  unsigned int best_length = 0;
  unsigned int best_distance = 0;
  
  Node *htree = 0;
  PQueue *priority_queue = (PQueue*)malloc(sizeof(PQueue)); 
  init_heap(priority_queue,512); // safe value for WLN

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  unsigned char code[CSIZE] = {0};
  std::vector<unsigned char> bitstream;
  
  // fill up the forward window
  while(fpos < buff_size){
    if(!fread(&ch,sizeof(unsigned char),1,ifp))
      reading_data = 0;

    buffer[fpos++] = ch;
    if(!reading_data)
      break;
  } 

  if(!buffer[backreference]){
    fprintf(stderr,"Error: no data!\n");
    free(buffer);
    return false;
  }


  while(buffer[backreference]){
    reference = false;
    unsigned int distance = 0;
    unsigned int length   = 0;

    for(unsigned int i=0;i<buff_size;i++){
      if(i >= backreference && !length)
        break;
      
      if(buffer[i] == buffer[backreference+length]){
        length++;
        if(!distance)
          distance = backreference - i;
      }
      else if(length > 2){
        //reference = true;

        // // if good backreference, break here


      

        best_distance = distance;
        best_length = length; 

        distance = 0;
        length = 0;
      }
      else{
        // reset
        distance = 0;
        length = 0;
      }
    }

    if(reference){

    
      
    
    
    }
    else{

      // write the character code
      // construct tree based on C values
      for(edge=curr->transitions;edge;edge=edge->nxt){
        Node *n = AllocateNode(edge->ch,edge->c);
        insert_term(n,priority_queue);
      }

      htree = ConstructHuffmanTree(priority_queue);
      if(!htree){
        fprintf(stderr,"Huffman tree allocation fault\n");
        return false;
      }

      if(priority_queue->size){
        fprintf(stderr,"Error: queue is not being fully dumped on tree creation\n");
        return false;
      }

      unsigned int clen = WriteHuffmanCode(htree,buffer[backreference],code);
      
      for(unsigned int c = 0;c<clen;c++){
        bitstream.push_back(code[c]);
        
        // if(bitstream.size() == 258){
        //   stream_to_bytes(bitstream);
        //   bitstream.clear();
        // }

      }
      
      memset(code,0,CSIZE);

      DeleteHuffmanTree(htree); 

      for(edge=curr->transitions;edge;edge=edge->nxt){
        if(edge->ch == buffer[backreference]){
          curr = edge->dwn;
          edge->c++;
          break;
        }
      }

      LeftShift(buffer,buff_size,1);
      if(reading_data){
        if(!fread(&ch,sizeof(unsigned char),1,ifp))
          reading_data = 0;
        else
          buffer[buff_size-1] = ch;
      }
    }
  }

  // get the last ones still in the stream
  stream_to_bytes(bitstream);

  // need to add a end of file marker. 

  delete_heap(priority_queue);
  free(buffer);
  return true;
}


bool WLNDECODE(FILE *ifp, FSMAutomata *wlnmodel){
  
  unsigned char *buffer = (unsigned char*)malloc(sizeof(unsigned char)*buff_size); 
  memset(buffer,0,buff_size);

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

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  // init the first huffman tree for state, change on switch now
  for(edge=curr->transitions;edge;edge=edge->nxt){
    Node *n = AllocateNode(edge->ch,edge->c);
    insert_term(n,priority_queue);
  }

  tree_root = ConstructHuffmanTree(priority_queue);
  htree = tree_root;

  while(fread(&ch,sizeof(unsigned char),1,ifp)){
    for(int i=7;i>=0;i--){
      unsigned char bit = (ch & (1 << i))?1:0;


      if(reading_bits){




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
        else if(htree->ch){
          
          LeftShift(buffer,backreference,1);
          fputc(htree->ch,stdout);
          buffer[backreference-1] = htree->ch;

          // move the machine to the next state
          for(edge=curr->transitions;edge;edge=edge->nxt){
            if(edge->ch == htree->ch){
              curr = edge->dwn;
              edge->c++;
              break; 
            }
          }

          // delete the huffman tree 
          DeleteHuffmanTree(tree_root);

          // create the new queue
          for(edge=curr->transitions;edge;edge=edge->nxt){
            Node *n = AllocateNode(edge->ch,edge->c);
            insert_term(n,priority_queue);
          }

          // create the new tree
          tree_root = ConstructHuffmanTree(priority_queue);
          htree = tree_root;

          if(!htree){
            fprintf(stderr,"Huffman tree allocation fault\n");
            return false;
          }

          if(priority_queue->size){
            fprintf(stderr,"Error: queue is not being fully dumped on tree creation\n");
            return false;
          }
        }

      }
    }

  }
  
  delete_heap(priority_queue);
  free(buffer);
  return true;
}


/* ######################################################################################### */


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




  FILE *fp = 0; 
  fp = fopen(input,"rb");
  if(fp){
    if(opt_mode == 1)
      WLNENCODE(fp,wlnmodel);
    else if (opt_mode == 2)
      WLNDECODE(fp,wlnmodel);

    fclose(fp);
  }
  else{
    fprintf(stderr,"Error: could not open file at %s\n",input);
    return 1;
  }

  delete wlnmodel;
  return 0;
}

/* ######################################################################################### */
