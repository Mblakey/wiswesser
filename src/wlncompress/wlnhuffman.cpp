#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <stack>
#include <map>

#include "rfsm.h"
#include "wlnmatch.h"
#include "wlndfa.h"

unsigned int opt_mode = 0;
unsigned int opt_verbose = false;

const char *input;
const char *trainfile; 

/* priority queue and huffman tree code */
struct Node;
struct PQueue;

struct Node{
  unsigned freq;
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

void debug_heap(PQueue *heap){
  for (unsigned int i=0;i<heap->size;i++)
    fprintf(stderr,"%c, ",heap->arr[i] ? '1':'0');

  fprintf(stderr,"\n");
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
  if(priority_queue->size == 1){
    Node *root = AllocateNode('*',0);
    Node *f = pop_front(priority_queue);
    f->p = root; 
    root->l = f;
    return root;
  }
  

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

void DeleteHuffmanTree(Node *root){
  if(root->l)
    DeleteHuffmanTree(root->l);
  if(root->r)
    DeleteHuffmanTree(root->r);
  free(root);
}

/* builds the code in reverse and writes to stream */
void WriteHuffmanCode(Node *root,unsigned char ch, std::string &cstream){
  
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
    fprintf(stderr,"Error: could not find %c in states huffman tree\n",ch);
    return;
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


bool encode_file(FILE *ifp, FSMAutomata *wlnmodel, std::map<FSMState*,Node*> &huffman_lookup){
  
  unsigned char ch = 0;
  fseek(ifp,0,SEEK_END);
  unsigned int file_len = ftell(ifp); // bytes
  fseek(ifp,0,SEEK_SET);

  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  std::string cstream;

  while(fread(&ch,sizeof(unsigned char),1,ifp)){

    Node *htree = huffman_lookup[curr]; 
    if(!htree){
      fprintf(stderr,"Huffman tree is uninitialised\n");
      return false;
    }

    if(ch == '\n'){
      WriteHuffmanCode(htree,0,cstream);
      curr = wlnmodel->root; 
    }
    else{
      WriteHuffmanCode(htree,ch,cstream);
      for(edge=curr->transitions;edge;edge=edge->nxt){
        if(edge->ch == ch){
          curr = edge->dwn;
          break;
        }
      }
    }
  }

  // // write a byte from the root of the machine indicating EOF
  Node *rtree = huffman_lookup[wlnmodel->root];
  WriteHuffmanCode(rtree,0,cstream);

  if(opt_verbose){
    fprintf(stderr,"%d to %d bits: %f compression ratio\n",
            file_len*8,(unsigned int)cstream.size(),
            (double)(file_len*8)/(double)cstream.size() );
  }

  stream_to_bytes(cstream);
  return true;
}



bool decode_file(FILE *ifp, FSMAutomata *wlnmodel,std::map<FSMState*,Node*> &huffman_lookup){
  
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  unsigned char bit_char = 0;
  unsigned int bit_pos = 0;
  fread(&bit_char,sizeof(unsigned char),1,ifp);

  for(;;){

    unsigned char ch_read = 0;
    Node *htree = huffman_lookup[curr]; 
    if(!htree){
      fprintf(stderr,"Error: dead huffman tree\n");
      return false;
    }

    while(!ch_read){
      unsigned char bit = bit_char & (1 << (7 - bit_pos) ) ? 1:0;
      bit_pos++; 

      // fetch next char from stream
      if (bit_pos == 8){
        if(!fread(&bit_char,sizeof(unsigned char),1,ifp))
          return true;
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
          if(curr == wlnmodel->root)
            return true; 
          else{
            fputc('\n',stdout);
            curr = wlnmodel->root; 
            ch_read = 1;
          }
        }
        else{
          fputc(ch_read,stdout);
          for(edge=curr->transitions;edge;edge=edge->nxt){
            if(edge->ch == ch_read){
              curr = edge->dwn;
              break; 
            }
          }

        }
      }

    }
  }

  return true;
}


static void DisplayUsage()
{
  fprintf(stderr, "wlncompress2 <options> <input> > <out>\n");
  fprintf(stderr, "<options>\n");
  fprintf(stderr, "  -c          compress input\n");
  fprintf(stderr, "  -d          decompress input\n");
  fprintf(stderr, "  -t          add an optional train file for edge frequencies (see wlntrain)\n");
  fprintf(stderr, "  -v          verbose debugging statements on\n");
  exit(1);
}


static void ProcessCommandLine(int argc, char *argv[])
{
  const char *ptr = 0;
  int i,j;

  input = (const char *)0;
  trainfile = (const char *)0;

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

        case 't':
          if(i < argc - 1){
            i++;
            trainfile = argv[i];
          }
          else{
            fprintf(stderr,"Error: -t must be followed with a file\n");
            DisplayUsage();
          }
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

  FSMAutomata *wlnmodel = CreateWLNDFA(); // build the model 
  wlnmodel->MakeAccept(wlnmodel->root);
  
  // to every accept, add the null character used to reset the machine
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    if(wlnmodel->states[i]->accept)
      wlnmodel->AddTransition(wlnmodel->states[i],wlnmodel->states[i],'\0');
  }

  if(!trainfile){
    if(opt_verbose)
      fprintf(stderr,"Warning: using order-0 probabilities\n");

    wlnmodel->AssignEqualProbs();
  }
  else{
    FILE *tfp = 0; 
    tfp = fopen(trainfile,"rb");
    if(!tfp){
      fprintf(stderr,"Error: cannot open train file\n");
      return 1;
    }

    unsigned int i=0;
    unsigned int freq = 0; 
    while(fread(&freq,sizeof(unsigned int),1,tfp))
      wlnmodel->edges[i++]->c = freq;
    
    for(unsigned int i=0;i<wlnmodel->num_edges;i++){
      if(!wlnmodel->edges[i]->c){
        fprintf(stderr,"Warning - null count for edge %d, using 1\n",i);
        wlnmodel->edges[i]->c = 1;
      }
    }

    if(opt_verbose)
      fprintf(stderr,"Warning: train file read, using probabilities\n");

    fclose(tfp);
  }

  // create a huffman tree for each state
  FSMState *s = 0; 
  FSMEdge *e = 0;
  std::map<FSMState*,Node*> huffman_lookup; 
  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    
    PQueue *priority_queue = (PQueue*)malloc(sizeof(PQueue)); 
    init_heap(priority_queue,256); // safe value for alphabet 
    
    s = wlnmodel->states[i];
    for(e=s->transitions;e;e=e->nxt){
      Node *n = AllocateNode(e->ch,e->c);
      insert_term(n,priority_queue);
    }

    Node *htree = ConstructHuffmanTree(priority_queue);
    huffman_lookup[s] = htree;

    free(priority_queue->arr);
    free(priority_queue);
    priority_queue = 0;
  }


  FILE *fp = 0; 
  fp = fopen(input,"rb");
  if(fp){
    if(opt_mode == 1)
      encode_file(fp,wlnmodel,huffman_lookup);
    else if (opt_mode == 2)
      decode_file(fp,wlnmodel,huffman_lookup);

    fclose(fp);
  }
  else{
    fprintf(stderr,"Error: could not open file at %s\n",input);
    return 1;
  }

  for(unsigned int i=0;i<wlnmodel->num_states;i++){
    s = wlnmodel->states[i];
    DeleteHuffmanTree(huffman_lookup[s]);
  }

  delete wlnmodel;
  return 0;
}