#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ctree.h"

unsigned int debug_id = 1; 

Node *AllocateTreeNode(unsigned char ch, unsigned int id){
  Node* n = (Node*)malloc(sizeof(Node));
  n->c = 0; 
  n->ch = ch;
  n->leaves = 0;
  n->id = id;
  
  n->vine = 0;

  n->prev = 0; // debugging only, remove when algorithm tested
  return n; 
}

Edge *AllocateTreeEdge(Node *p, Node *c){
  if(!p||!c)
    fprintf(stderr,"Error: dead nodes\n");

  Edge* e = (Edge*)malloc(sizeof(Edge));
  e->dwn = c; 
  e->nxt = 0;
  if(!p->leaves)
    p->leaves = e;
  else{
    Edge *q = p->leaves;
    while(q->nxt){
      q = q->nxt;
    }
    q->nxt = e;
  }
  
  c->prev = p; 
  return e; 
}


void RReleaseTree(Node *root){
  if(!root){
    fprintf(stderr,"Error: dead node\n"); 
    return; 
  }

  Edge *e = root->leaves;
  Edge *p = 0; 

  while(e){
    p = e;
    RReleaseTree(e->dwn);
    e = e->nxt;
    free(p); 
  }

  free(root); 
}

void RdotTraverse(Node *n, FILE *fp){
  if(n){
    fprintf(fp,"\t%d [label=\"%c (%d)\"];\n", n->id, n->ch, n->c);
    Edge *e = 0; 
    for(e = n->leaves;e;e=e->nxt){
      fprintf(fp,"\t%d -> %d;\n",n->id,e->dwn->id); 
      RdotTraverse(e->dwn, fp);
    }
  
//    if(n->vine)  
//      fprintf(stdout,"\t%d -> %d [color=\"green\"];\n",n->id,n->vine->id); 
  }
}

void WriteDotFile(Node *root, FILE *fp){
  fprintf(fp,"digraph ContextTree{\n");
  RdotTraverse(root,fp); 
  fprintf(fp,"}\n");
}



Node *search_tree(const char *str, Node*root, unsigned int k){ 
  unsigned int i=0; 
  Node *v = root;
  Edge *e = 0;
  unsigned char ch = *str; 
  while(ch && i < k ){
    bool found = false;
    for(e=v->leaves;e;e=e->nxt){
      if(e->dwn->ch ==ch){
        v = e->dwn;
        found = true;
        break;
      }
    }   
    if(!found)
      return 0;

    i++;
    ch = *(++str); 
  }
  return v; 
}

/* Builds a forward context tree, character by character, this will return the last
 * node seen/created, meaning decrementing through the vines to the root node is a way
 * to quickly parse contexts, will return the longest node, use prev if context_len == top
 * */
void BuildContextTree(Node *root, const char *str,unsigned int context_len,bool update_exclusion){

  Node *t = root;
  unsigned int j = 0;
  Node *prev = 0;

  while(j < context_len){
    bool found = false;
    t = root; 
    for(unsigned int k = j; k < context_len;k++){
      found=false;
      for(Edge *e=t->leaves;e;e=e->nxt){
        if(e->dwn->ch == str[k]){
          t = e->dwn;
          found = true;
        }
      }
    
      if(!found){
        Node *n = AllocateTreeNode(str[k], debug_id++);
        AllocateTreeEdge(t, n);      
        if(t==root){
          n->vine = t;
          t->c++; 
        }

        t = n;
        t->c=1;
      }

    }

    if(prev)
      prev->vine = t; 
    prev = t; 
    
    if(found){
      
      if(update_exclusion){
        if(j==0)
          t->c++;
      }
      else
        t->c++;

#define BASIC_SCALE 1
#if BASIC_SCALE
      if(t->c == 64)
        t->c = 16;
#endif
    }

    j++;
  }

  return;
}


void RunbackContext(Node *node){
  unsigned int bpos = 0; 
  unsigned char buffer[32] = {0};
  Node *n= node;
  while(n){
    buffer[bpos++] = n->ch; 
    n = n->prev; 
  }

  for(int i = bpos-1;i>=0;i--)
    fputc(buffer[i],stderr);
  fputc('\n',stderr); 
}


double PredictPPM(  const char *message, unsigned ch_pred, Node *tree, unsigned char mode, 
                    unsigned int context_len, unsigned int avaliable_chars){
 
  if(!tree){
    fprintf(stderr,"Error: dead root\n");
    return 0; 
  }

  // with the vines we should only need one traversal, find the highest order message
  Node *t = tree; 
  Edge *e = 0; 

  unsigned int char_occurance = 0;
  
  unsigned int order = context_len; 
  long double prob = 0.0;
  long double weight = 1.0; 

  for(unsigned int i=0;i<context_len;i++){ 
    for(e=t->leaves;e;e=e->nxt){
      if(e->dwn->ch == message[i])
        t = e->dwn;
    }
  }

  while(t){
    unsigned int edges = 0; 
    unsigned int Co = 0; // t->c-1
    for(e=t->leaves;e;e=e->nxt){
      edges++; 
      if(e->dwn->ch == ch_pred){
        char_occurance = e->dwn->c;
      }
      Co += e->dwn->c; 
    }

    unsigned int unique_occurance = edges; // i think so! 
    
    double e_o = 1.0; 
    if(mode=='A')
      e_o = 1/(double)(Co+1);
    else if (mode == 'B' && unique_occurance && Co > 1)
      e_o = unique_occurance/(double)(Co); 
    else if(mode == 'C' && unique_occurance)
      e_o = unique_occurance/(double)(Co+unique_occurance);
    
//    fprintf(stderr,"order %d:: q_o: %d, c_o: %d, C_o: %d, e_o: %f\n",order,edges, char_occurance, Co,e_o);   
    if(char_occurance){
      double WoPo = 0.0;  
      if(mode == 'A')
        WoPo = (char_occurance)/(double)(Co+1);
      else if (mode == 'B' && char_occurance > 1 && Co)      
        WoPo = (char_occurance-1)/(double)(Co);
      else if (mode == 'C')
        WoPo = (char_occurance)/(double)(Co+unique_occurance); 

      WoPo *= weight; 
      prob += WoPo; 
    }
    
    weight *= e_o; 
    t = t->vine;

    order--; 
  }

  // order-1 model, equal probs
  prob+= weight * 1/(double)avaliable_chars; 
  return prob; 
}





/* Predict PPM with exclusions, function takes in memory buffer pointer, and returns
 * number of values placed in the array before a probability was found. A zero value
 * indicates an escape should be encoded. 
 *
 * How the escape character range is calculated is moved into the arithmetic coder, where 
 * T is incremented, and counts changed to reflect the inclusion of the escape character
*/
unsigned int PredictPPMExclusion(   const char *message, unsigned ch_pred, 
                                    Node *tree, unsigned char mode, unsigned int context_len,
                                    unsigned int *frequency_buffer)
{ 
  Node *t = tree; 
  Edge *e = 0; 
  unsigned int char_occurance = 0;
  unsigned char ch = *message; 

  unsigned int order = 0;
  unsigned int freq_pos = 0; 
  
  bool ascii_exclude[255] = {false}; 
  
  while(ch){ 
    order++; 
    for(e=t->leaves;e;e=e->nxt){
      if(e->dwn->ch == ch)
        t = e->dwn;
    }
    ch = *(++message); 
  }
  
  unsigned int excluded = 0; 
  while(t){
    unsigned int edges = 0; 
    unsigned int Co = 0;
    for(e=t->leaves;e;e=e->nxt){
      if(e->dwn->ch == ch_pred)
        char_occurance = e->dwn->c;
      
      if(!ascii_exclude[e->dwn->ch]){
        Co += e->dwn->c;
        edges++; 
      }
    }
    
    // some stuff for mode 'B' here!
    if(!char_occurance || (mode == 'B' && char_occurance == 1)){
      for(e=t->leaves;e;e=e->nxt){
        if(!ascii_exclude[e->dwn->ch]){
          if(mode=='B' && e->dwn->c <= 1)
            continue; // only skip characters for B when they hit the conditions

          ascii_exclude[e->dwn->ch]= true;
          excluded++;
        }
      }
    }

//    fprintf(stderr,"order  %d:: q_o: %d, c_o: %d, C_o: %d, e_o: %f, weight: %f\n",order,edges, char_occurance, Co,e_o,weight);  

    frequency_buffer[freq_pos++] = 0; // an escape should be encoded.
    
    order--;
    t = t->vine;
  }
  
  double prob = 0.0;
  double weight = 0.0; 
 // fprintf(stderr,"order -1:: excluded: %d,  weight: %f\n",excluded,weight);  
  prob += weight * 1/(double)(5-excluded); 
  return prob; 
}



double PredictPPMLazyExclusion(const char *message, unsigned ch_pred, Node *tree, unsigned char mode, unsigned int context_len){
 
  if(!tree){
    fprintf(stderr,"Error: dead root\n");
    return 0; 
  }

  int i= strlen(message);
  unsigned int count = 0;
  while(i>=0 && count < context_len){ 
    i--;
    count++;
  }  
  i++; 
  message += i;  // move the pointer to context prediction <addres> + ds:0x[context_len] 

  // with the vines we should only need one traversal, find the highest order message
  Node *t = tree; 
  Edge *e = 0; 
  unsigned int char_occurance = 0;
  unsigned char ch = *message; 

  unsigned int order = 0;
  double prob = 0.0;
  double weight = 1.0; 
  
  while(ch){ 
    order++; 
    for(e=t->leaves;e;e=e->nxt){
      if(e->dwn->ch == ch)
        t = e->dwn;
    }
    ch = *(++message); 
  }
 

  while(t){
    unsigned int edges = 0; 
    unsigned int Co = 0; //t->c-1;

    for(e=t->leaves;e;e=e->nxt){
      edges++;
      if(e->dwn->ch == ch_pred)
        char_occurance = e->dwn->c;

      Co += e->dwn->c; 
    }
    
    double e_o = 1.0;
    if(mode=='A')
      e_o = 1/(double)(Co+1);
    else if (mode == 'B' && edges && Co)
      e_o = edges/(double)(Co); 
    else if(mode == 'C' && edges)
      e_o = edges/(double)(Co+edges);

//    fprintf(stderr,"order  %d:: q_o: %d, c_o: %d, C_o: %d, e_o: %f, weight: %f\n",order,edges, char_occurance, Co,e_o,weight);  

    if(char_occurance){
      double WoPo = 0.0;  
      if(mode == 'A'){
        WoPo = (char_occurance)/(double)(Co+1);
        return WoPo * weight;
      }
      else if (mode == 'B' && char_occurance > 1){      
        WoPo = (char_occurance-1)/(double)(Co);
        return WoPo * weight;
      }
      else if (mode == 'C'){
        WoPo = (char_occurance)/(double)(Co+edges); 
        return WoPo * weight;
      }
    }
    
    weight *= e_o; 
    order--;
    t = t->vine;
  }
  
  prob += weight * 1/(double)(5); 
  return prob; 
}
