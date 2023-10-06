/*##############################################################

Define the structs needed for automaton based representation 
of regular languages, handles DFA, NFA, eNFA

###############################################################*/



#ifndef REG_FSM_H
#define REG_FSM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <set>
#include <stack>
#include <map>

#define REALLOC 512  // reallocate 512 more states or edges at a time
#define REASONABLE 1024 // basic init 

enum FSMType{DFA=0,NFA=1,eNFA=2};

struct FSMState;
struct FSMEdge; 

struct FSMState{
  bool accept;
  unsigned int id; 
	FSMEdge *transitions;
  FSMState *access[255]; // instant access array for matching

  FSMState():accept{0},id{0},transitions{0}{
    for (unsigned int i=0;i<255;i++)
      access[i] = 0; 
  }
};

struct FSMEdge{
  unsigned int id;
	unsigned char ch; // 0 will be a epsilon transition 
	FSMState *dwn;
	FSMEdge *nxt;

  FSMEdge(){
    id = 0;
    ch = 0;
    dwn = 0;
    nxt = 0;
  };
};

struct FSMAutomata{

	unsigned char type;
	FSMState** states;
  FSMEdge**  edges;
	FSMState* root; 

  bool alphabet [255];

  unsigned int num_states;
  unsigned int num_edges;
  unsigned int num_accepts;

  unsigned int max_states; 
  unsigned int max_edges; 

	FSMAutomata(unsigned int node_size, unsigned int edge_size){
    states =  (FSMState**)malloc(sizeof(FSMState*) * node_size);
    edges =   (FSMEdge**)malloc(sizeof(FSMEdge*) * edge_size);
    
    for (unsigned int i=0;i<node_size;i++)
      states[i] = 0;
    for (unsigned int i=0;i<edge_size;i++)
      edges[i]  = 0;
    
    max_states = node_size;
    max_edges = edge_size;

    root = 0;
    type = 0;

    num_states = 0;
    num_edges = 0;
    num_accepts = 0;

    for(unsigned char ch = 0;ch<255;ch++)
      alphabet[ch] = 0;
  };

	~FSMAutomata(){
		DeleteFSM();
	}

  /* only need to call explicity if trying to save memory*/
  bool DeleteFSM(){
    for (unsigned int i=0;i<max_states;i++)
			delete states[i];
    
    for (unsigned int i=0;i<max_edges;i++)
      delete edges[i];
    
    free(states); 
    free(edges);
    return true;
  }

  bool ReallocateStateSpace(){

    unsigned int old_max = max_states;
    max_states += REALLOC; // add 512
    FSMState **new_space = (FSMState**)realloc(states,sizeof(FSMState*) *max_states);
    if(!new_space){
      fprintf(stderr,"Error: reallocation of states memory failed\n");
      return false;
    }
    else if(new_space != states) // if the memory chunk has moved
      states = new_space;

    // null the new space - safety
    for(unsigned int i = old_max+1;i<max_states;i++)
      states[i] = 0;

    return true;
  }

  bool ReallocateEdgeSpace(){
    unsigned int old_max = max_edges;
    max_edges += REALLOC; // add 512
    FSMEdge** new_space = (FSMEdge**)realloc(edges,sizeof(FSMEdge*) *max_edges);
    if(!new_space){
      fprintf(stderr,"Error: reallocation of edge memory failed\n");
      return false;
    }
    else if(new_space != edges) // if the memory chunk has moved
      edges = new_space;

    // null the new space - safety
    for(unsigned int i =old_max+1;i<max_edges;i++)
      edges[i] = 0;

    if(!edges)
      return false;
    else
      return true;
  }



	FSMState* AddState(bool accept=false){
    		
    FSMState *state = 0;
    if(num_states == max_states && !ReallocateStateSpace())
      return 0;
    

    state = new FSMState;
    if(!root)
      root = state; 

    state->id =  num_states++;
    states[state->id] = state;
    state->accept = 0; 
    state->transitions = 0;

    // null potential jump table
    for (unsigned int i=0;i<255;i++)
      state->access[i] = 0;

    if(accept)
      state = MakeAccept(state);
    
    return state;
	}

  /* must be called on the state in the memory pool */
  bool RemoveState(FSMState *state){

    unsigned int h_edges = 0; 
    FSMEdge * hold_edges[255];

    FSMEdge *e = 0; 
    for (e = state->transitions; e; e = e->nxt)
      hold_edges[h_edges++] = e;

    for (unsigned int i=0;i<h_edges;i++)
      RemoveTransition(state,hold_edges[i]);
  
      
    // find any instance of the state in a previous edge 
    std::vector<FSMEdge*> eremove;
    eremove.reserve(255);
    for(unsigned int i=0;i<num_states;i++){
      FSMState *par = states[i];
      if(par){
        for(e = par->transitions;e;e=e->nxt){
          if(e->dwn == state)
            eremove.push_back(e);
        }
    
        for(FSMEdge *er : eremove)
          RemoveTransition(par,er);

        eremove.clear();
      }
    }    

    if(state->accept)
      num_accepts--; 
    
    num_states--;
    states[state->id] = NULL; // null the global
    delete state;
    return true; 
  }


	FSMEdge* AddTransition(FSMState *src, FSMState *trg, unsigned char ch){
    if(!src|| !trg){
      fprintf(stderr,"Error: attempting transition on dead pointers\n");
      return 0;
    }

    // if the edge is already there
    FSMEdge *e = 0;
    for (e = src->transitions;e;e = e->nxt){
      if(e->dwn == trg && e->ch == ch){
        //fprintf(stderr,"Warning: Adding duplicate transition\n");
        return e;
      }
       
    }

    FSMEdge *edge = 0; 
   
    if(num_edges == max_edges && !ReallocateEdgeSpace())
      return 0;
    
    edge = new FSMEdge;
    edge->id = num_edges++; 
    edges[edge->id] = edge;
		edge->ch = ch;
    edge->dwn = trg; 
    alphabet[ch] = true; // store for look up on transitions

    e = src->transitions;
    if(!e)
      src->transitions = edge;
    else{
      while(e->nxt)
        e = e->nxt;
      
      e->nxt = edge;
    }

		return edge;
	}

  // need to find in the pools to null all safely
  bool RemoveTransition(FSMState* state,FSMEdge* edge){

    bool found = false;
    FSMEdge *e = state->transitions; 
    FSMEdge *prev = 0; 
    
    while(e){
      if(e == edge){
        found = true; 
        break; 
      }
      prev = e;
      e = e->nxt;
    }

    if(found){
      if(prev)
        prev->nxt = e->nxt; 

      num_edges--; 
      edges [edge->id]= NULL;
      delete edge;
      return true; 
    }
    else
      return false; 
  }

  FSMState *MakeAccept(FSMState *state){
    if(!state){
      fprintf(stderr,"Error: attempting accept on nullptr\n");
      return 0;
    }
    
    state->accept = true;
    num_accepts++; 
    return state;
  }


  /* allows a character access to avalible states
  - DFA only */
  bool InitJumpTable(){
    if(type != DFA)
      return false; 

    FSMEdge *edge = 0;
    FSMState *state = 0;
    for (unsigned int i=0;i<max_states;i++){
      state = states[i];
      if(state){
        for (edge=state->transitions;edge;edge=edge->nxt)
          state->access[edge->ch] = edge->dwn;
      }
    }
    return true;
  }

  /* Travels the states and categorises the FSM into 
  DFA,NFA,e-NFA*/
  void Categorize(unsigned int verbose=0){
    for (unsigned int i=0; i<num_states;i++){
      FSMState *state = states[i];
      FSMEdge *edge = 0;
      bool seen[255] = {false};
      for(edge = state->transitions;edge;edge=edge->nxt){
        if(seen[edge->ch]){
          
          if(verbose)
            fprintf(stderr,"%c transition seen twice on state: %d\n",edge->ch,state->id);

          if(type<NFA)
            type = NFA; 
        }
        
        seen[edge->ch] = true;
        if(!edge->ch)
          type = eNFA; 
      }
    }

    if(verbose){
      switch(type){
        case DFA:
          fprintf(stderr,"fsm type: DFA\n");
          break;
        case NFA:
          fprintf(stderr,"fsm type: NFA\n");
          break;
        
        case eNFA:
        fprintf(stderr,"fsm type: e-NFA\n");
          break;
      }
    }
  }

  // dead non accepting states
  bool DeadState(FSMState* state){
    if(state->accept)
      return false;
    else{
      FSMEdge *e = 0;
      for(e=state->transitions;e;e=e->nxt){
        if(e->dwn != state)
          return false;
      }
      return true;
    }
  }


  /* in-line modification of the fsm, removing
  all unreachable and dead states */
  bool RemoveUnreachables(){
    std::map<FSMState*,bool> visited; 
    std::map<FSMState*,bool> dead_states; 
    std::stack<FSMState*>     stack;
    FSMState *top             = 0; 
    FSMEdge *e                = 0;

    stack.push(root);
    while(!stack.empty()){
      top = stack.top();
      visited[top] = true;
      stack.pop();

      if(DeadState(top))
        dead_states[top] = true;

      for (e=top->transitions;e;e=e->nxt){
        if(!visited[e->dwn]){
          stack.push(e->dwn);
          visited[e->dwn] = true;
        }
      }
    }

    for (unsigned int i=0;i<num_states;i++){
      if(states[i] && (!visited[states[i]] || dead_states[states[i]] ))
        RemoveState(states[i]);  
    }

    Reindex();
    return true;
  }

  void ReOrderRoot(){
    if(root == states[0])
      return;

    for (unsigned int i=0;i< max_states;i++){
      if(states[i] == root){
        states[i] = states[0];
        states[0] = root;
        return;
      }
    }
  }

  /* reindex the whole fsm - make sure root is always at zero */
  void Reindex(){

    ReOrderRoot();

    unsigned int state_index  = 0; 
    unsigned int edge_index   = 0; 

    for (unsigned int i=0;i< max_states;i++){
      FSMState *s =  states[i];
      FSMEdge  *e =  edges[i];
      if(s){
        states[i] = NULL;
        s->id = state_index;
        states[state_index++] = s;
      }
      if(e){
        edges[i] = NULL;
        e->id = edge_index;
        edges[edge_index++] = e;
      }
    }

    num_states = state_index;
    num_edges = edge_index; 
  }



  bool DumpFSM(const char *filename){
    fprintf(stderr,"--- FSM Dump ---\n");
    fprintf(stderr, "states:  %d\n"
                    "edges:   %d\n"
                    "accepts: %d\n",num_states,num_edges,num_accepts);

    FILE *fp = 0;
    fp = fopen(filename,"w");
    if(!fp){
      fprintf(stderr,"Error: could not create file pointer\n");
      return false;
    }

    fprintf(fp, "digraph FSMdigraph {\n");
    fprintf(fp, "  rankdir = LR;\n");
    for (unsigned int i=0; i<num_states;i++)
    {

      FSMState *node = states[i];
      if(!node)
        continue;

      fprintf(fp, "  %d", node->id);

      if (node == root){
        if(node->accept)
          fprintf(fp, "[shape=doublecircle, fillcolor=red,style=filled];\n");
        else
          fprintf(fp, "[shape=circle, fillcolor=red,style=filled];\n");
      }
      else if (node->accept)
        fprintf(fp, "[shape=doublecircle];\n");
      else
        fprintf(fp, "[shape=circle];\n");
    
      FSMEdge *edge = 0;
      for (edge = node->transitions;edge;edge = edge->nxt){
        FSMState *child = edge->dwn;
        fprintf(fp, "  %d -> %d [label=\"%c\"]\n",node->id,child->id,edge->ch ? edge->ch:'*');
      }
    }

    fprintf(fp, "}\n");
    fclose(fp);

    return true;
  }

  void DebugFSM(){
    fprintf(stderr,"--- FSM pointer dump ---\n"); 
    fprintf(stderr,"states[ %d ]:\n",num_states);
    for(unsigned int i=0;i<num_states;i++)
      fprintf(stderr,"%d: %p\n",i,states[i]);
    fprintf(stderr,"\n");

    fprintf(stderr,"edges[ %d ]:\n",num_edges);
    for(unsigned int i=0;i<num_edges;i++){
      if(edges[i])
        fprintf(stderr,"%d: %p --> %p\n",i,edges[i],edges[i]->dwn);
      else
        fprintf(stderr,"%d: %p\n",i,edges[i]);
    }
      
    fprintf(stderr,"\n");

  }
};

#endif