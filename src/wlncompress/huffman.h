
#ifndef HUFFMAN_H
#define HUFFMAN_H

typedef struct Node{
  unsigned int freq;
  unsigned char ch; 
  struct Node *l; 
  struct Node *r; 
  struct Node *p;
} Node;

typedef struct{
  Node **arr; 
  unsigned int cap;
  unsigned int size;  
} PQueue; 

Node *AllocateNode(unsigned char ch, unsigned int f);
bool init_heap(PQueue *heap, unsigned int cap);
void free_heap(PQueue *priority_queue);

void shift_right(unsigned int low,unsigned int high, PQueue *heap);
void heapify(unsigned int low, unsigned int high, PQueue *heap);

bool blind_insert(Node *term, PQueue *heap);
bool insert_term(Node *term, PQueue *heap);

Node *pop_front(PQueue *priority_queue);

Node *ConstructHuffmanTree(PQueue *priority_queue);
void free_huffmantree(Node *root);

unsigned int WriteHuffmanCode(Node *root,unsigned char ch, unsigned char *code);
void ReserveCode(const char*code,unsigned char sym,Node* tree_root);

#endif
