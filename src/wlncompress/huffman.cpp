
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <stack>
#include "huffman.h"

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


void free_heap(PQueue *priority_queue){
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

void free_huffmantree(Node *root){
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

/* uses tree splicing to reverse a specific pattern */
void ReserveCode(const char*code,unsigned char sym,Node* tree_root){
  unsigned char ch = *code; 
  Node *htree = tree_root;

  if(!htree)
    return;

  unsigned int clen = 0;
  while(ch){

    if(ch == '0'){
      if(!htree->l){
        htree->l = AllocateNode(htree->ch,0);
        htree->l->p = htree;
        htree->ch = 0;
      }
      
      htree = htree->l;
    }
      
    else if (ch == '1'){
      if(!htree->r){
        htree->r = AllocateNode(htree->ch,0);
        htree->r->p = htree;
        htree->ch = 0;
      }

      htree = htree->r;
    }
    
    clen++;
    ch = *(++code);
  }

  if(clen<1)
    fprintf(stderr,"Error: reserving single bit code undefined\n");
  

  // now there are some splicing conditions
  Node *splice_parent = htree->p; 

  unsigned int p = 0;
  if(htree == splice_parent->l){
    p = 1;
    if(!splice_parent->r){
      splice_parent->r = htree;
      splice_parent->l = AllocateNode(sym,1);
      splice_parent->l->p = splice_parent;
      return;
    }
  }
  else if(htree == splice_parent->r){
    p = 2;
    if(!splice_parent->l){
      splice_parent->l = htree;
      splice_parent->r = AllocateNode(sym,1);
      splice_parent->r->p = splice_parent;
      return;
    }
  }
   
  Node *branch = AllocateNode(0,0);
  htree->p = 0;

  // move the splice locations to the new node
  branch->l = splice_parent->l;
  branch->r = splice_parent->r; 

  // remove from the parent
  splice_parent->l = 0;
  splice_parent->r = 0;

  // if they exist, add their parents 
  if(branch->l)
    branch->l->p = branch;
    
  if(branch->r)
    branch->r->p = branch;

  // add branches parent as the splice node.
  branch->p = splice_parent;
  if(p==1){
    splice_parent->r = branch;
    splice_parent->l = AllocateNode(sym,1);
    splice_parent->l->p = splice_parent;
  }
  else if(p==2){
    splice_parent->l = branch;
    splice_parent->r = AllocateNode(sym,1);
    splice_parent->r->p = splice_parent;  
  }
  return;
}

