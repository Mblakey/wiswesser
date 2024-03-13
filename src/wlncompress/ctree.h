#ifndef CTREE_H
#define CTREE_H

#include <stdio.h>

typedef struct Edge Edge;
typedef struct Node Node;

struct Edge{
  Node *dwn; 
  Edge *nxt; 
};

struct Node{
  unsigned int id;
  unsigned char ch; 
  unsigned int c; 
  Edge *leaves;
  Node *vine;

  Node *prev; // debugging only, remove when algorithm working 
 };


Node *AllocateTreeNode(unsigned char ch, unsigned int id);
Edge *AllocateTreeEdge(Node *p, Node *c);

void RReleaseTree(Node *root);
void WriteDotFile(Node *root, FILE *stream);

void BuildContextTree(Node *root,const char *str, unsigned int context_len,bool update_exclusion);
void RunbackContext(Node *node);


#endif 
