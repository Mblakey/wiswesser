#ifndef TPPM_H
#define TPPM_H

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

void BuildContextTree(Node *root,const char *str, unsigned int context_len);
bool BuildContextTreeUpdateExclusion(Node *root,const char *str, unsigned int context_len);
void RunbackContext(Node *node);

double PredictPPMExclusion(const char *message, unsigned ch_pred, Node *tree, unsigned char mode, unsigned int context_len);
double PredictPPMLazyExclusion(const char *message, unsigned ch_pred, Node *tree, unsigned char mode, unsigned int context_len);


#endif 
