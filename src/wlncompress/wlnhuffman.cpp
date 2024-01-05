#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <string>
#include <stack>
#include <deque>
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


/* idea can we look back 5 states e.g 5 characters from an accept, and get back to an accept 
accept restriction is to eliminate weird half locant */

bool encode_file(FILE *ifp, FSMAutomata *wlnmodel, std::map<FSMState*,PQueue*> &queue_lookup){
  
  unsigned char ch = 0;
  unsigned int bytes_read = 0;

  Node *htree = 0;
  PQueue *priority_queue = 0;
  FSMState *curr = wlnmodel->root;
  FSMEdge *edge = 0;

  std::string cstream;


  unsigned int link_size = 0;
  std::string test;
  std::deque<FSMState*> accept_links; 

  unsigned char pattern = 91; // 1 above Z
  std::map<std::string,unsigned char> pattern_map;
  std::map<unsigned char,std::string> get_pattern;
  
  // read line allows the look ahead
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

      if(link_size == 5){
        if(curr->accept && (accept_links[0]->accept || accept_links[0] == wlnmodel->root) ){
          //fprintf(stderr,"%s could be a pattern\n",test.c_str());

          unsigned char pattern_ch = pattern_map[test];
          if(!pattern_ch){
            pattern_ch = pattern;
            get_pattern[pattern] = test;
            pattern_map[test] = pattern++;
          }

#define TEST 1
#if TEST
          // duplicate edges ignored
          FSMEdge *dict_edge = wlnmodel->AddTransition(accept_links.front(),accept_links.back(),pattern_ch);
          dict_edge->c = 1;
#endif
        }


  #if TEST
        for(unsigned int i=0;i<accept_links.size()-1;i++){
          FSMState *s = accept_links[i];
          FSMEdge *e = s->transitions;
          for(e=s->transitions;e;e=e->nxt){
            if(e->ch == test[i]){
              if(accept_links[i+1] != e->dwn){
                fprintf(stderr,"no go! failing on %c(%d)\n",test[i],test[i]);
                return false;
              }
            }
          }
        }
  #endif

        test.clear();
        accept_links.clear();
        link_size = 0;
      }
      else if (ch == '\n'){
        test.clear();
        accept_links.clear();
        link_size = 0;
      }
      else{
        accept_links.push_back(curr);
        test+=ch;
        link_size++;
      }

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

      FSMEdge *standard = 0;
      FSMEdge *jump = 0;
      unsigned int j = i;
      for(edge=curr->transitions;edge;edge=edge->nxt){
      
        if(edge->ch == ch)
          standard = edge;
        else if(edge->ch > 'Z'){

          std::string dict_str = get_pattern[edge->ch];

          fprintf(stderr,"potentially looking for %s\n",dict_str.c_str());
          fprintf(stderr,"buffer: %s",(buffer) );
          fprintf(stderr,"remaining: %s\n",(buffer+i) );
          
          unsigned int search = 0;
          for(j=i;j<BUFF_SIZE;j++){
            if(!buffer[j] || (dict_str[search] != buffer[j]) )
              break;
            
          
            search++;
            if(search == 5)
              break;
          }

          if(search==5){
            ch = edge->ch;
            jump = edge;
          }
        
        }
      }

      if(jump){
        curr = jump->dwn;
        jump->c++;
        
        i = j; 
        test.clear();
        accept_links.clear();
        link_size = 0;
      }
      else if(standard){
        curr = standard->dwn;
        standard->c++;
      }


      WriteHuffmanCode(htree,ch,cstream);
      DeleteHuffmanTree(htree); 
    }
  }


  fprintf(stderr,"%d patterns\n",pattern-'Z');
  // // write a byte from the root of the machine indicating EOF
  priority_queue = queue_lookup[wlnmodel->root];
  for(edge=wlnmodel->root->transitions;edge;edge=edge->nxt){
    Node *n = AllocateNode(edge->ch,edge->c);
    insert_term(n,priority_queue);
  }

  htree = ConstructHuffmanTree(priority_queue);
  WriteHuffmanCode(htree,0,cstream);
  DeleteHuffmanTree(htree); 

  if(opt_verbose){
    fprintf(stderr,"%d to %d bits: %f compression ratio\n",
            bytes_read*8,(unsigned int)cstream.size(),
            (double)(bytes_read*8)/(double)cstream.size() );
  }

  stream_to_bytes(cstream);
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

  FSMAutomata *wlnmodel = CreateWLNDFA(); // build the model 

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
    init_heap(priority_queue,256); // safe value for alphabet 
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

    while(priority_queue->size){
      free(priority_queue->arr[priority_queue->size-1]);
      priority_queue->size--;
    }

    free(priority_queue->arr);
    free(priority_queue);
    priority_queue = 0;
  }

  delete wlnmodel;
  return 0;
}