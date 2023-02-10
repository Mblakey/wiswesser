
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector>
#include <stack>
#include <map>
#include <deque>
#include <algorithm> // std::sort

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

// --- forward declarations --- (functions can have indirect recursion)
struct WLNSymbol; 
WLNSymbol* ParseNonCyclic(const char *wln, unsigned int len); 
WLNSymbol* ParseCyclic(const char *wln, unsigned int len);

// --- inputs ---  
const char *wln; 
const char *dotfile; 
unsigned int pos;   // global track of where we are in the wln string
unsigned int glob_len;   // total length of the wln string

// --- options --- 
static bool opt_wln2dot = false;
static bool opt_valstrict = false; 
static bool opt_verbose = false;
static bool opt_canonical = false; 
static bool opt_returnwln = false;

// macros 
#define REASONABLE 1024

// pending notation
bool pending_ring = false; 

// --- memory --- 
std::vector<WLNSymbol*> mempool; 

static void empty_mempool(){
  for (WLNSymbol* allocedwln : mempool){  // if error free all the mempool --> stop leak 
    free(allocedwln); 
    allocedwln = 0;
  }
}


// character type and pending state types
enum WLNType {SINGLETON = 0, BRANCH = 1, LINKER = 2, TERMINATOR = 3}; 
enum WLNState{PASS = 0, NONCYCLIC = 1, CYCLIC = 2, POLYCYCLIC=3, PERICYCLIC=4, BRIDGED=5, SPIRO=6};


// rule 2 - hierarchy - rules have diverged due to end terminator char
std::map<unsigned char,unsigned int> char_hierarchy = 
{
  {' ',1}, {'-',2}, {'/',3}, 
  {'0',4},{'1',5},{'2',6},{'3',7},{'4',8},{'5',9},{'6',10},{'7',11},{'8',12},{'9',13},
  {'A',14},{'B',15},{'C',16},{'D',17},{'E',18},{'F',19},{'G',20},{'H',21},{'I',22},
  {'J',23},{'K',24},{'L',25},{'M',26},{'N',27},{'O',28},{'P',29},{'Q',30},{'R',31},
  {'S',32},{'T',33},{'U',34},{'V',35},{'W',36},{'X',37},{'Y',38},{'Z',40},{'&',41}
};

// quick validity checker 
bool ValidCharParse(const char *wln, unsigned int len){
  for (unsigned int i=0;i<len;i++){
    if (!char_hierarchy[wln[i]]){
      fprintf(stderr,"Error: invalid wln character %c in string\n",wln[i]);
      return false;
    }
  }
  return true; 
}


struct WLNSymbol{

  unsigned char ch; 
  unsigned int type; 

  unsigned int allowed_edges;
  unsigned int num_edges;  

  WLNSymbol *prev; // should be a single term - wln symbol only has one incoming
  std::vector<WLNSymbol*> children; // linked list of next terms chains 

  bool init(unsigned char inp_char){
    ch = inp_char;

    switch(ch){
      
      case '0':
      case '1':
      case '2':
      case '3': 
      case '4':
      case '5': 
      case '6': 
      case '7': 
      case '8':
      case '9':
        type = SINGLETON; 
        allowed_edges = 2;
        break;
      
      case 'A':
        type = SINGLETON;
        allowed_edges = 2;
        break;

      case 'B':   // boron
        type = BRANCH;
        allowed_edges = 3; 
        break;
      
       case 'C':  // shortcut carbon atom
        type = BRANCH;
        allowed_edges = 4;
        break;
      
      case 'D': 
        type = SINGLETON; 
        allowed_edges = 2;
        break;
      
      case 'E':  // halogens
      case 'F':
      case 'G': 
      case 'I': 
        type = BRANCH; 
        allowed_edges = 3;
        break;
      
      case 'H': // closing hydrogen
        type = TERMINATOR; 
        allowed_edges = 1;
        break;

      case 'J': // generic symbol for halogen 
        type = BRANCH; 
        allowed_edges = 3;
        break;
      
      case 'K': 
        type = BRANCH; 
        allowed_edges = 4;
        break;

      case 'L':
        type = LINKER;
        allowed_edges = 2;
        break;
      
      case 'M':
        type = BRANCH; 
        allowed_edges = 2; 
        break;

      case 'N': 
        type = BRANCH; 
        allowed_edges = 3;
        break; 

      case 'O':
        type = SINGLETON; 
        allowed_edges = 2;
        break;

      case 'P':
        type = BRANCH; 
        allowed_edges = 5;
        break;

      case 'Q': 
        type = TERMINATOR; 
        allowed_edges = 1;
        break;

      case 'R':
        type = SINGLETON; 
        allowed_edges = 2;
        break;

      case 'S': 
        type = BRANCH;
        allowed_edges = 6;
        break;

      case 'T':
      case 'U': 
        type = LINKER;
        allowed_edges = 2;
        break; 

      case 'V': 
        type = SINGLETON; 
        allowed_edges = 2;
        break;

      case 'W':
        type = LINKER; 
        allowed_edges = 2;
        break;

      case 'X':
        type = BRANCH; 
        allowed_edges = 4;
        break;

      case 'Y':
        type = BRANCH; 
        allowed_edges = 3;
        break;

      case 'Z': 
        type = TERMINATOR;
        allowed_edges = 1;
        break;

      case '&':
        type = TERMINATOR;
        allowed_edges = 1;
        break;

      case ' ':
      case '-':
      case '/': 
        type = LINKER; 
        allowed_edges = 2;
        break;

      case '\0':
        fprintf(stderr,"Error: end of string null char accessed!\n");
        return false;

      default:
        fprintf(stderr,"Error: invalid wln symbol parsed: %c\n",ch);
        return false; 
    }

    prev = (WLNSymbol*)0;
    return true; 
  } 

};

/* for std::sort on canonicalise */
bool char_comp(const WLNSymbol* a, const WLNSymbol* b){
  return char_hierarchy[a->ch] > char_hierarchy[b->ch]; 
}


WLNSymbol* AllocateWLNSymbol(unsigned char ch){

  WLNSymbol *wln = (WLNSymbol*)malloc( sizeof(WLNSymbol) ); 
  if(wln->init(ch)) 
    mempool.push_back(wln);
  else 
    return (WLNSymbol *)0;
  
  return wln; 
}


bool handle_hypervalence(WLNSymbol *problem){

  // this will always lead to a positive ion species 
  switch(problem->ch){


    case 'M':   // tranforming this to a N is an easy solve
      if(opt_verbose)
        fprintf(stderr,"   transforming hypervalent M --> N\n");
      
      problem->ch = 'N';
      break; 

    case 'N':
      if(opt_verbose)
        fprintf(stderr,"   transforming hypervalent N --> K\n");
      
      problem->ch = 'N';
      break;


    case 'Y':  // can go to an X
      if (opt_verbose)
        fprintf(stderr,"   transforming hypervalent Y --> X\n");
      
      problem->ch = 'X';
      break; 

    default:
      if (opt_verbose) 
        fprintf(stderr,"Error: cannot handle hypervalent symbol: %c\n",problem->ch);
      return false;
  }

  return true; 
}


/* add src to the children vector of trg, handle hypervalent bonds if possible */
bool add_symbol(WLNSymbol* src, WLNSymbol *trg){

  // handle exotic bonding - lookback 
  if (trg->ch == 'U'){
    if (trg->prev && trg->prev->ch == 'U')
      src->num_edges+=3;
    else
      src->num_edges+=2;
  }
  else
    src->num_edges++;

  if (src->num_edges > src->allowed_edges){
    if(!opt_valstrict){
      if(!handle_hypervalence(src))
        return false;
    }
    else {
      fprintf(stderr,"Error: (strict mode) hypervalence on WLN character %c\n",src->ch);
      return false; 
    }
  }
  
  if (trg->num_edges < trg->allowed_edges){
    trg->children.push_back(src);
    trg->num_edges++;
  }
  else {
    if(!opt_valstrict){
      if(!handle_hypervalence(trg))
        return false;
      else
        trg->children.push_back(src);
    }
    else {
      fprintf(stderr,"Error: (strict mode) hypervalence on WLN character %c\n",trg->ch);
      return false; 
    }
  }   

  return true;
}


/* handles special symbol combinations to enter main subroutines */
unsigned int pending_states(unsigned char ch){

  switch(ch){

    case '-':
      if (pending_ring)
        pending_ring = false; 
      else 
        pending_ring = true;
      break;

    case ' ':
      if(pending_ring){
        if (opt_verbose)
          fprintf(stderr,"   ring system detected in branch\n");
        pending_ring = false; // we handle it here
        return CYCLIC;
      }
      break;


    default:
      if(pending_ring)
        pending_ring = false; 
      break; 
  }
  
  return 0;
}

/* this can perform the normal backtrack, excluding  '&' closure */ 
WLNSymbol* backtrack_stack(std::stack<WLNSymbol*> &wln_stack){
  WLNSymbol *tmp = 0; 
  while(!wln_stack.empty()){
    tmp = wln_stack.top();
    if (tmp->type == BRANCH)  
      return tmp;
    wln_stack.pop();
  }
  
  // if it goes all the way that string complete
  return (WLNSymbol*)0;  
}

/* forces both the '&' closure and its parents branch to come off stack */
WLNSymbol *force_closure(std::stack<WLNSymbol*> &wln_stack){
  WLNSymbol *tmp = 0; 
  unsigned int popped = 0;
  while(!wln_stack.empty()){
    tmp = wln_stack.top();
    if (tmp->type == BRANCH && popped > 1)  
      return tmp;
    wln_stack.pop();
    popped++; 
  }

  // if it goes all the way that string complete
  return (WLNSymbol*)0;
}

/* parses NON CYCLIC input string, mallocs graph nodes and sets up graph based on symbol read */
WLNSymbol* ParseNonCyclic(const char *wln, unsigned int len){
  
  std::stack <WLNSymbol*> wln_stack; 
  WLNSymbol *prev = 0;
  WLNSymbol *root = 0;

  WLNSymbol* created_wln = AllocateWLNSymbol(wln[0]);
  if (!created_wln)
    return (WLNSymbol*)0;
  wln_stack.push(created_wln);
  
  root = created_wln; 
  
  for (unsigned int i = 1; i<len; i++){
    
    // this can detect based on the char, and then change created_wln as per given special
    unsigned int special = pending_states(wln[i]);
    switch(special){

      case WLNState::CYCLIC: 
        created_wln = ParseCyclic(wln+=i,len);
        break;



      default:
        created_wln = AllocateWLNSymbol(wln[i]);
        break;
    }

    if (!created_wln)
      return (WLNSymbol*)0; 
    
    
    prev = wln_stack.top();
    wln_stack.push(created_wln); // push all of them

    if(!add_symbol(created_wln,prev))
      return (WLNSymbol*)0; 

    // options here depending on the forced closure of '&'
    if (created_wln->type == TERMINATOR){
      if (created_wln->ch == '&' && prev->type == BRANCH)
        prev = force_closure(wln_stack);
      else
        prev = backtrack_stack(wln_stack);
    }   
  }

  // for aid with tree reordering, a '&' can be used to close all wln notation
  
  prev = created_wln; // last created 
  created_wln = AllocateWLNSymbol('&');
  prev->children.push_back(created_wln);

  return root; // return start of the tree
}




WLNSymbol* parse_locant(unsigned int locant_start, unsigned int locant_end){
  WLNSymbol *branch_root = 0;

  unsigned char substr[REASONABLE] = {'\0'}; 
  unsigned int arr_len = locant_end - locant_start;

  if (arr_len > REASONABLE){
    fprintf(stderr,"Error: branch in ring system exceeds 1024 characters - termination\n");
    return (WLNSymbol*)0;
  }

  unsigned int access = 0;  // we dont include the locant symbol - shift by 1
  for (unsigned int tmp = locant_start+1; tmp<locant_end+1; tmp++)
    substr[access++] = wln[tmp];


  if (opt_verbose)
    fprintf(stderr,"   bonding %s to locant %c\n",substr,wln[locant_start]);

  branch_root = ParseNonCyclic((const char*)substr, arr_len);
  return branch_root;
} 


/* parse the 'first' CYCLIC wln species */
WLNSymbol* ParseCyclic(const char *wln, unsigned int len){

  WLNSymbol *prev = 0;
  WLNSymbol *root = 0;
  WLNSymbol *jsymbol = 0;
  WLNSymbol* created_wln = AllocateWLNSymbol(wln[0]);
  root = prev = created_wln; 

  unsigned int j_pos = 0; 
  for(unsigned int i = 1; i < len; i++){

    created_wln = AllocateWLNSymbol(wln[i]); 
    prev->children.push_back(created_wln);
    prev = created_wln; 
    if(wln[i] == 'J'){
      j_pos = i;
      jsymbol = created_wln;
      break; 
    }
  }

  // looks for J ring closure 
  if(!j_pos || !jsymbol){ 
    fprintf(stderr,"Error: ring system not closed with a J\n");
    return (WLNSymbol *)0;
  }

  // look for immediate ring exit <-- important for joined cyclics  'L6TJ&'
  if (j_pos+1 < len && wln[j_pos+1] == '&'){
    if(opt_verbose)
      fprintf(stderr,"   forced '&' ring closure detected\n");
    return root;
  }
     
  // current i position on J, therefore add 1, should be on space, add 2 get the first locant letter
  unsigned int locant_start = j_pos+2; 
  unsigned int locant_end = 0; 
  for (unsigned int i=j_pos+2;i<len;i++){
    
    if(wln[i] == ' '){
      locant_end = i-1; // zero indexed
        
      WLNSymbol *branch_root = parse_locant(locant_start,locant_end);
      if (!branch_root){
        fprintf(stderr,"Error: could not parse locant - return nullptr\n");
        return (WLNSymbol *)0;
      }

      // create a locant node and bind branch
      WLNSymbol *locant_node = AllocateWLNSymbol(wln[locant_start]);
      locant_node->children.push_back(branch_root);

      // bind the locant to the 'J' value of the ring --> evaluation to smiles later
      jsymbol->children.push_back(locant_node);
      locant_start = i+1;
    }
    else if (i == len -1){     // on notation end, last locant

      locant_end = i; // zero indexed
      WLNSymbol *branch_root = parse_locant(locant_start,locant_end);
      if (!branch_root){
        fprintf(stderr,"Error: could not parse locant - return nullptr\n");
        return (WLNSymbol *)0;
      }
        

      // create a locant node and bind branch
      WLNSymbol *locant_node = AllocateWLNSymbol(wln[locant_start]);
      locant_node->children.push_back(branch_root);

      // bind the locant to the 'J' value of the ring --> evaluation to smiles later
      jsymbol->children.push_back(locant_node);
    }
    
  }
  
  return root; 
}




/* uses a character sorting method to arrange the WLN symbols according to rule 2
BFS ensures a starting position */
bool CanonicoliseWLN(WLNSymbol *root){
  std::deque<WLNSymbol*> wln_queue; 
  wln_queue.push_back(root);
  
  WLNSymbol *top = 0;
  while(!wln_queue.empty()){
    top = wln_queue.front();
    wln_queue.pop_front();

    if (top->children.size() > 1)
      std::sort(top->children.begin(),top->children.end(),char_comp);
    
    for (WLNSymbol* child : top->children){
      wln_queue.push_back(child);
    }
  }

  return true; 
}


/* reforms WLN string with DFS ordering */
std::string ReformWLNString(WLNSymbol *root){
  std::string res; 

  std::stack<WLNSymbol*> wln_stack; 
  std::map<WLNSymbol*,bool> visit_map; 
  wln_stack.push(root);

  WLNSymbol *top = 0;
  while(!wln_stack.empty()){
    top = wln_stack.top(); 
    wln_stack.pop();
    visit_map[top] = true; 

    res.push_back(top->ch);

    for (WLNSymbol *child: top->children){
      if (!visit_map[child]){
        wln_stack.push(child);
      }
    }
  }

  return res; 
}

/* dump wln tree to a dotvis file */
void WLNDumpToDot(FILE *fp){

  // set up index map for node creation
  std::map <WLNSymbol*, unsigned int> index_map; 
  unsigned int glob_index = 0;
  for (WLNSymbol *node : mempool){
    index_map[node] = glob_index; 
    glob_index++;
  }
    
  fprintf(fp,"digraph WLNdigraph {\n");
  fprintf(fp,"  rankdir = LR;\n");
  for (WLNSymbol *node : mempool){
    fprintf(fp,"  %d",index_map[node]);
    fprintf(fp, "[shape=circle,label=\"%c\"];\n",node->ch);
    for (WLNSymbol *child : node->children){
      fprintf(fp,"  %d", index_map[node]);
      fprintf(fp," -> ");
      fprintf(fp,"%d\n",index_map[child]);
    }
  } 
  fprintf(fp,"}\n");
}

static void DisplayUsage(){
  fprintf(stderr,"wln-writer <options> < input (escaped) >\n");
  fprintf(stderr,"<options>\n");
  fprintf(stderr,"  -v | --verbose                print messages to stdout\n");
  fprintf(stderr,"  -s | --strict                 fail on hypervalence, no symbol correction\n");
  fprintf(stderr,"  -c | --canonical              perform wln canonicalise procedure\n");
  fprintf(stderr,"  -r | --return-wln             return wln after altering procedure(s)\n");
  fprintf(stderr,"  --wln2dot <dotfile.dot>       dump wln tree to dot file\n");
  exit(1);
}

static void ProcessCommandLine(int argc, char *argv[]){

  const char *ptr=0; 
  int i,j; 

  wln     = (const char*)0;
  dotfile = (const char*)0;
  
  if (argc < 2)
   DisplayUsage();

  j=0; 
  for (i=1;i<argc;i++){

    ptr = argv[i];

    if (ptr[0]=='-' && ptr[1])
      switch (ptr[1]){

        case 'c':
          opt_canonical = true;
          break; 

        case 'r':
          opt_returnwln = true;
          break; 

        case 's':
          opt_valstrict = true; 
          break;

        case 'v':
          opt_verbose = true; 
          break; 

        case '-':
          if (!strcmp(ptr, "--wln2dot")){
            opt_wln2dot = true; 
            if (i == argc - 1){
              fprintf(stderr,"Error: --wlndot requires a <file>.dot as next arguement\n");
              DisplayUsage();
            }
            i++; // increment argv
            if (argv[i][0] != '-')
              dotfile = argv[i];
            else{
              fprintf(stderr,"Error: --wlndot requires a <file>.dot as next arguement\n");
              DisplayUsage();
            } 
            break;
          }
          else if (!strcmp(ptr, "--strict")){
            opt_valstrict = true; 
            break;
          }
          else if (!strcmp(ptr, "--verbose")){
            opt_verbose = true; 
            break;
          }
          else if (!strcmp(ptr, "--canonical")){
            opt_canonical = true; 
            break;
          }
          else if (!strcmp(ptr, "--return-wln")){
            opt_returnwln = true; 
            break;
          }

        default:
          fprintf(stderr,"Error: unrecognised input %s\n",ptr);
          break;
      }
    
    else switch(j++){
      case 0: wln = ptr; break;
      default: break;
    }
      
  }

  return;
}


int main(int argc, char *argv[]){
  ProcessCommandLine(argc, argv);
  if(!wln)
    return 1;
  else{
    pos = 0; 
    glob_len = strlen(wln);
  }

  if(!ValidCharParse(wln, glob_len)) // quick check the wln formula 
    return 1;
  
  WLNSymbol *root = 0;

  if(opt_verbose)
    fprintf(stderr,"-- parsing input: %s\n",wln);

  if(wln[0] == 'L' || wln[0] == 'T')
    root = ParseCyclic(wln, glob_len);
  else
    root = ParseNonCyclic(wln, glob_len);
  
  if (!root){
    if(opt_verbose)
      fprintf(stderr,"   failed parse\n");
    empty_mempool();
    return 1; 
  }

  if(opt_verbose)
    fprintf(stderr,"\n");

  if (opt_canonical){
    if(opt_verbose)
      fprintf(stderr,"-- canonicaling wln...\n");

    CanonicoliseWLN(root);

    if(opt_verbose)
      fprintf(stderr,"\n");
  }

  if (opt_wln2dot){
    if(opt_verbose)
      fprintf(stderr,"-- dumping wln to dot file...\n");
    FILE *fp = 0; 
    fp = fopen(dotfile, "w");
    if (!fp)
      fprintf(stderr,"Error: could not write %s as .dot file - skipping\n",dotfile);
    else
      WLNDumpToDot(fp);

    if(opt_verbose)
      fprintf(stderr,"\n");
  }

  if(opt_returnwln){
    if(opt_verbose)
      fprintf(stderr,"-- reforming wln string...\n");
    
    std::string res = ReformWLNString(root);
    
    if(opt_verbose){
      fprintf(stderr,"   %s\n",res.c_str());
      fprintf(stderr,"\n");
    }
    else
      fprintf(stderr,"%s\n",res.c_str()); 
  }
  
  empty_mempool();
  return 0;
}