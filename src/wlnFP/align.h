#ifndef _ALIGN_H_
#define _ALIGN_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stack>
#include <string>

const unsigned int EDIT_UPPER_BOUND = 32; 

static int max(int a, int b) {
  return a < b ? b : a;
}

static int max(int a, int b, int c) {
  return max(max(a,b),c);
}


static int SWScoreChar(char a, char b){
  // minus one for a missing char, will handle in the SW loop
  if (a==b)
    return 1;  
  else
    return -1; 
}

static int NWScoreChar(char a, char b){
  // minus one for a missing char, will handle in the SW loop
  if (a==b)
    return +1;  
  else
    return -1; 
}


struct Edit{
  char ch;
  char mut_ch; 

  // must be static for switch
  static const unsigned char ADD = 1; 
  static const unsigned char DEL = 2;
  static const unsigned char MUT = 3;

  // 1 = add, 2 = del, 3 = Mut
  unsigned char instruct; 
  Edit()
    :ch{'\0'},mut_ch{'\0'},instruct{0}{};

  Edit(unsigned char _instruct, char _ch)
    :ch{_ch},mut_ch{'\0'},instruct{_instruct}{};

  Edit(unsigned char _instruct, char _ch, char _mut)
    :ch{_ch},mut_ch{_mut},instruct{_instruct}{};


  // copy hunting 
  //Edit(Edit const& copy){PRINT_COPY("Edit");}

    void SetEditValues(unsigned int _instruct, char _ch){
      instruct = _instruct; 
      ch = _ch; 
    }

    void SetEditValues(unsigned char _instruct, char _ch, char _mut){
      instruct = _instruct; 
      ch = _ch; 
      mut_ch = _mut; 
    }


  bool SendToFileVerbose(FILE *fp){
    if(!fp){
      fprintf(stderr, "ERROR: Invalid *fp for Edit::SendToFileVerbose()\n");
      return false;
    }
    switch(instruct){
      case ADD:
        fprintf(fp,"add(%c)",ch);
        return true;
      case DEL:
        fprintf(fp,"del(%c)",ch);
        return true;
      case MUT: 
        fprintf(fp,"mut(%c,%c)",ch,mut_ch);
        return true;
      default:
        fprintf(stderr, "ERROR: Invalid switch instruction for Edit::SendToFileVerbose()\n");
        return false;
    }
  }

  bool SendToFile(FILE *fp){
    if(!fp){
      fprintf(stderr, "ERROR: Invalid *fp for Edit::SendToFile()\n");
      return false;
    }
    switch(instruct){
      case ADD:{
        fprintf(fp, "%d%c", ADD,ch);
        return true;
      }
      case DEL:{
        fprintf(fp, "%d%c", DEL,ch);
        return true;
      }
      case MUT:{
        fprintf(fp, "%d%c%c", MUT,ch,mut_ch);
        return true;
      }
      default:
        fprintf(stderr, "ERROR: Invalid switch instruction for Edit::SendToFile()\n");
        return false;
    }
  }

};

// result should only get created once. and then passed down
// so pointers should be consistent 
struct ResultStruct{
  std::string a; 
  std::string b; 
  unsigned int num_changes=0; 
  Edit *pedit = 0; 

  ResultStruct(){};

  ResultStruct(unsigned int a_len, unsigned int b_len){
    a.reserve(a_len+1);
    b.reserve(b_len+1);
  }

  ~ResultStruct(){
    if (pedit){
      delete [] pedit;
      pedit = 0;
    }
  }

  bool display(FILE *fp){
    if (!fp){
      fprintf(stderr, "ERROR: Invalid *fp for AlignStruct::display()\n");
      return false;
    }
    fprintf(fp,"%s\t%s",a.c_str(),b.c_str());
    return true;
  }

  bool WriteInstructionsVerbose(FILE *fp){
    if (!fp){
      fprintf(stderr, "ERROR: Invalid *fp for AlignStruct::WriteInstructionsVerbose()\n");
      return false;
    }

    display(fp);
    fputc('\t',fp);
    fputc('[',fp);
    for (unsigned int i=0;i<num_changes;i++){
      pedit[i].SendToFileVerbose(fp);
      if (i != num_changes-1)
        fputs(" ,",fp);
    }
    fputc(']',fp);
    fputc('\n',fp);
    return true; 
  }

  bool WriteInstructionsOnly(FILE *fp){
    if (!fp){
      fprintf(stderr, "ERROR: Invalid *fp for AlignStruct::WriteInstructionsOnly()\n");
      return false;
    }

    for (unsigned int i=0;i<num_changes;i++){
      pedit[i].SendToFile(fp);
      if (i != num_changes-1)
        fputs("\t",fp);
    }
    fputc('\n',fp);
    return true; 
  }
};


struct TracePosition{
  int index; 
  int score;
  int x;
  int y;
  
  bool travel_left;
  bool travel_up;
  bool travel_up_left; 

  // no undefined behaviour <-- use new in malloc array
  TracePosition()
      :index{0},score{0},x{0},y{0},travel_left{false},travel_up{false},travel_up_left{false}{};

  TracePosition(int _index,int _score,int _x, int _y)
      :index{_index},score{_score},x{_x},y{_y},travel_left{false},travel_up{false},travel_up_left{false}{};

  void display(){
    printf("idx: %i, score: %i (%i,%i)\n",index,score, x,y);
  }

  bool SetTraceCoordinates(int _index,int _score,int _x, int _y){
    index = _index;
    score = _score;
    x = _x; 
    y = _y;
    return true; 
  }
};


 
struct AlignStruct{

  int n; 
  int m; 
  int *matrix; 
  TracePosition *ptrace = 0;
  unsigned int max_ptrace=0;

  AlignStruct(int _n, int _m)
        :n{_n},m{_m}{}; 

  ~AlignStruct(){
    //fprintf(stderr,"Destructor Called\n")]
    if (matrix)
      free(matrix);
    if (ptrace)
      delete [] ptrace;
    ptrace = 0;
    matrix = 0; 
  };

  bool Init(){
    //fprintf(stderr,"%u,%u %u %llu\n", n,m, (n*m),(n*m) *sizeof(int));
    matrix = (int*)malloc((n*m) *sizeof(int));
    if (!matrix){
      fprintf(stderr, "ERROR: Cannot allocate memory SAMatrix::Init()\n");
      return false; 
    }
    memset(matrix,0,(n*m) *sizeof(int)); 
    

    // worst case for the matrix plus 1 extra for safety, zero index stops array looping
    // using new as this calls the the default constructor
    ptrace = new TracePosition [m+n+1]; 
  

    if (!ptrace){
      fprintf(stderr, "ERROR: Invalid *ptrace SAMatrix::Init()\n");
      return false; 
    }

    return true; 
  }

  bool DisplayMatrix(const char *a, const char *b){
    if (!matrix)
      return false; 
    printf("|  ");
    for (int j=0; j<n ;j++)
      printf("%2c ", j ? a[j-1] : ' ');
    printf("|\n");
    for (int i=0; i<m ;i++){
      printf("|%2c", i ? b[i-1] : ' ');
      for (int j=0; j<n ;j++)
        printf("%2i ", matrix[i*n+j]);
      printf("|\n");
    }
    printf("\n");
    return true;
  }


  bool DisplayPtrace(){

    unsigned int iter=1;
    while(ptrace[iter].index != 0){
      ptrace[iter].display();
      iter++;
    }
    return true;
  }

  bool Set(int value, int x, int y){
    if (x < 0 || x >= n || y < 0 || y >= m){
      fprintf(stderr, "ERROR: Index out of Bounds, Max: (%i,%i), input(%i,%i) SAMatrix::Set()\n",n-1,m-1,x,y);
      //abort();
      return false; 
    }

    matrix[(y*n)+x] = value; 
    return true; 
  }

  int Get(int x, int y){
    if (x < 0 || x >= n || y < 0 || y >= m){
      fprintf(stderr, "ERROR: Index out of Bounds, Max: (%i,%i), input(%i,%i) SAMatrix::Get()\n",n-1,m-1,x,y);
      //abort();
      return 0; 
    }
    return matrix[(y*n)+x];
  }

  bool SmithWaterman(const char *a, const char *b, ResultStruct &result){
    
    int gap_penalty = -2;

     // keep track of highest score position
    int highest_score=0;
    int high_x = 0;
    int high_y = 0;

    int left=0;
    int up=0;
    int up_left=0;

    for (int i=1; i<n; i++) {
      for (int j=1; j<m; j++) {       

        left = gap_penalty + Get(i-1,j);   
        up =  gap_penalty + Get(i,j-1);   
        up_left = Get(i-1,j-1) + SWScoreChar(a[i-1],b[j-1]);  

        int score = max(up,left,up_left);
        if (score < 0)
          score=0;
                
        if (score >= highest_score){
          highest_score = score;
          high_x = i;
          high_y = j;
        }
        Set( score ,i,j);
      }
    }

    // matrix is filled, performing traceback
    std::stack<TracePosition> branch_stack;
    SWTraceBackIteration(highest_score, high_x,high_y,branch_stack);
    if(!AlignStrings(a,b,result))
      return false;
    else
      return true;
  }


  // MULTIPLE  EQUIVILENT ALIGNMENTS ARE OFF
  bool SWTraceBackIteration(  int highest_score, int x_coord, int y_coord, 
                              std::stack<TracePosition> &branch_stack){

    int iter = 0; 
    do { 
      
      if (!branch_stack.empty()){
        
        highest_score = branch_stack.top().score;
        x_coord = branch_stack.top().x;
        y_coord = branch_stack.top().y;
        iter = branch_stack.top().index + 1; // will be the next index from the stack track
        
       
        // all diagonals taken by default, no need check here
        if (branch_stack.top().travel_up){
          highest_score = Get(x_coord,y_coord-1);
          y_coord = y_coord-1;
          branch_stack.top().travel_up = false;
          if (!branch_stack.top().travel_up && !branch_stack.top().travel_left)
            branch_stack.pop();
        }

        else if (branch_stack.top().travel_left){
          highest_score = Get(x_coord-1,y_coord);
          x_coord = x_coord-1;
          branch_stack.top().travel_left = false;
          if (!branch_stack.top().travel_up && !branch_stack.top().travel_left)
            branch_stack.pop();
        }
        
        // there is guarteed to be a branch for 3 zero values on the padding, we dont want that
        if (highest_score != 0)
          ptrace[iter].SetTraceCoordinates(iter,highest_score,x_coord,y_coord);
      }
      
      while(highest_score !=0){
        unsigned int counter=0;
        iter++;

        ptrace[iter].SetTraceCoordinates(iter,highest_score,x_coord,y_coord);
      
        int left = Get(x_coord-1,y_coord);
        int up = Get(x_coord,y_coord-1);
        int up_left = Get(x_coord-1,y_coord-1);
        highest_score = max(up,left,up_left);    
        
        // allows equal paths to be known before setting coordinates
        if (highest_score == up_left){
          ptrace[iter].travel_up_left = true;
          counter++;
        }
        if (highest_score == up){
          ptrace[iter].travel_up = true;
          counter++;
        }
        if (highest_score == left){
          ptrace[iter].travel_left = true;
          counter++;
        }

        // prioritise diagonal and assignment coordinates
        if (ptrace[iter].travel_up_left){
          x_coord--;
          y_coord--;
          ptrace[iter].travel_up_left = false;
        }
        else if(ptrace[iter].travel_up){
          y_coord--;
          ptrace[iter].travel_up = false;
        }
        else{
          x_coord--;
          ptrace[iter].travel_left = false;
        }
      
        // TURN ON FOR MULTIPLE ALIGNMENT
        // if (counter > 1)
        //   branch_stack.push(state);
      }
      
    }while(!branch_stack.empty());
    max_ptrace = iter;
    return true; 
  }


  bool NeedlemanWuncsh(const char *a, const char *b, ResultStruct &result){

    //[i*n+j]
    // inits the matrix to negative decending values based off the gap penalty

    int gap_penalty = -1;
    for (int j=0;j<n;j++)
      matrix[(0*n)+j] = 0 + (j*gap_penalty);

    for (int i=0;i<m;i++)
      matrix[(i*n)+0] = 0 + (i*gap_penalty);
    

    // no need to keep track of highest score position here
    int left=0;
    int up=0;
    int up_left=0;

    for (int i=1; i<n; i++) {
      for (int j=1; j<m; j++) {       

        left = gap_penalty + Get(i-1,j);   
        up =  gap_penalty + Get(i,j-1);   
        up_left = Get(i-1,j-1) + NWScoreChar(a[i-1],b[j-1]);  

        int score = max(up,left,up_left);
      
        Set(score ,i,j);
      }
    }

    // would need to malloc to have this as an init struct member
    std::stack<TracePosition> branch_stack;
    NWTraceBackIteration(branch_stack);
    
    if(!AlignStrings(a,b,result))
      return false;
    else
      return true; 
  }

    // MULTIPLE EQUIVILENT ALIGNMENTS ARE OFF
  bool NWTraceBackIteration(std::stack<TracePosition> &branch_stack){
    unsigned int iter = 0; 
    
    // init loop to last value in matrix
    int x_coord = n-1;
    int y_coord = m-1;
    int highest_score = Get(x_coord,y_coord);
    int left=0;
    int up=0;
    int up_left=0;

    // first do while to handle the branch cases
    do{

      // place the stack operations here, same as SW

      // ...
      // ... 
      // ... 

      // need a safety net for traversing the top or side of the matrix, 
      // as that is possible but completes the path for NW
      while(x_coord + y_coord != 0){

        iter++;
        unsigned int counter=0;
        // create the state and add to the trace
       
        ptrace[iter].SetTraceCoordinates(iter,highest_score,x_coord,y_coord);

        // matrix edge safety conditions
        if(!x_coord){
          up = Get(x_coord,y_coord-1);
          highest_score = up;
        }
        else if(!y_coord){
          left = Get(x_coord-1,y_coord);
          highest_score = left; 
        }
        else{
          left = Get(x_coord-1,y_coord);
          up = Get(x_coord,y_coord-1);
          up_left = Get(x_coord-1,y_coord-1);
          highest_score = max(up,left,up_left);
        } 
        
        // allows equal paths to be known before setting coordinates
        if (highest_score == up_left){
          ptrace[iter].travel_up_left = true;
          counter++;
        }
        if (highest_score == up){
          ptrace[iter].travel_up = true;
          counter++;
        }
        if (highest_score == left){
          ptrace[iter].travel_left = true;
          counter++;
        }

        // prioritise diagonal and assignment coordinates
        if (ptrace[iter].travel_up_left){
          x_coord--;
          y_coord--;
          ptrace[iter].travel_up_left = false;
        }
        else if(ptrace[iter].travel_up){
          y_coord--;
          ptrace[iter].travel_up = false;
        }
        else{
          x_coord--;
          ptrace[iter].travel_left = false;
        }

        
        // TURN ON FOR MULTIPLE ALIGNMENT
        //if (counter > 1)
        //  branch_stack.push(state);
        
      }

    }while(!branch_stack.empty());
    max_ptrace = iter;
    return true; 
  }


  bool AlignStrings(const char *a, const char *b, ResultStruct &result){
    
    unsigned int changes = 0;
    unsigned int iter = max_ptrace;

    // init only for if the start for SW is non-zero
    char a_ch = '-'; 
    char b_ch = '-'; 
    int current_x = ptrace[iter].x; 
    int current_y = ptrace[iter].y; 
    int prev_x = 0;
    int prev_y = 0;

    // destructor handles memory
    result.pedit = new Edit [EDIT_UPPER_BOUND]; 

    while(ptrace[iter].index != 0){

      if (changes > EDIT_UPPER_BOUND)
        return false;
      
      Edit edit; 
      current_x = ptrace[iter].x;
      current_y = ptrace[iter].y;

      if (current_x ==0){
        a_ch = '-';
        b_ch = b[current_y-1];
        result.pedit[changes].SetEditValues(1,b_ch);
        changes++;
      }
      else if (current_y ==0){
        a_ch = a[current_x-1];
        b_ch = '-';
        result.pedit[changes].SetEditValues(2,a_ch);
        changes++;
      }
      else{

        // set up a diagonal condition <-- diagonal is a kept character
        if (current_x == prev_x+1 && current_y == prev_y+1){
          a_ch = a[current_x-1];
          b_ch = b[current_y-1];

          if (a_ch != b_ch){
            result.pedit[changes].SetEditValues(3,a_ch,b_ch);
            changes++;
          }
        }
          
        else{
          if (current_x == prev_x+1){
            a_ch = a[current_x-1];
            b_ch = '-';
            result.pedit[changes].SetEditValues(2,a_ch);
            changes++;
          }
          else if (current_y == prev_y+1){
            a_ch = '-';
            b_ch = b[current_y-1];
            result.pedit[changes].SetEditValues(1,b_ch);
            changes++;
          }
        }
      }
      result.a.push_back(a_ch);
      result.b.push_back(b_ch);

      prev_x = current_x;
      prev_y = current_y; 
      iter--;
    }
    result.num_changes = changes;
    return true;
  }

};



unsigned int WLNAlignment(const char *a, const char *b);

#endif //_ALIGN_H_
