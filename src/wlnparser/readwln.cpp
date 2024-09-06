/*********************************************************************
 
Author : Michael Blakey
Description: WLN reader C file, labelled as cpp to link with cpp toolkits.  

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/babelconfig.h>

#include "parser.h"

#define MAX_DEGREE 8
#define REASONABLE 32

#define DEBUG_FUNCS 1

#define SPACE_READ  0x01 
#define DIGIT_READ  0x02 
#define RING_READ   0x04

#define CARBON  6
#define NITRO   7
#define OXYGEN  8

typedef struct symbol_t symbol_t; 

/* returns atomic number of element packing */
u16 get_atomic_num(u8 high, u8 low){
  switch (high){
    case 'A':
      switch (low) {
        case 'C':
          return 89;
        case 'G':
          return 47;
        case 'L':
          return 13;
        case 'M':
          return 95;
        case 'R':
          return 18;
        case 'S':
          return 33;
        case 'T':
          return 85;
        case 'U':
          return 79;
      }
      break; 

    case 'B':
      switch (low) {
        case 0: 
          return 5;
        case 'A':
          return 56;
        case 'E':
          return 4;
        case 'H':
          return 107;
        case 'I':
          return 83;
        case 'K':
          return 97;
        case 'R':
          return 35;
      }
      break; 

    case 'C':
      switch (low) {
        case 0:
          return 6; 
        case 'A':
          return 20;
        case 'D':
          return 48;
        case 'E':
          return 58;
        case 'F':
          return 98;
        case 'M':
          return 96;
        case 'N':
          return 112;
        case 'O':
          return 27;
        case 'R':
          return 24;
        case 'S':
          return 55;
        case 'U':
          return 29;
      }
      break; 

    case 'D':
      switch (low) {
        case 'B':
          return 105;
        case 'S':
          return 110; 
        case 'Y':
          return 66; 
      }
      break;

    case 'E':
      switch (low) {
        case 0:
          return 35; 
        case 'R':
          return 68; 
        case 'S':
          return 99; 
        case 'U':
          return 63; 
      }

    case 'F':
      switch (low) {
        case 0:
          return 9; 
        case 'E':
          return 26;
        case 'L':
          return 114;
        case 'M':
          return 100;
        case 'R':
          return 87;
      }
      break;

    case 'G':
      switch (low) {
        case 0: 
          return 17;
        case 'A':
          return 31;
        case 'D':
          return 64;
        case 'E':
          return 32;
      }
      break;

    case 'H':
      switch (low) {
        case 'E':
          return 2;
        case 'F':
          return 72;
        case 'G':
          return 80;
        case 'O':
          return 67;
        case 'S':
          return 108;
      }
      break;

    case 'I':
      switch (low) {
        case 0: 
          return 53;
        case 'N':
          return 49;
        case 'R':
          return 77;
      }
      break;

    case 'K':
      switch (low) {
        case 0:  
          return 7;
        case 'R':
          return 36;
        case 'A':
          return 19;
      }
      break;

    case 'L':
      switch (low) {
        case 'A':
          return 57;
        case 'I':
          return 3;
        case 'R':
          return 103;
        case 'U':
          return 71;
        case 'V':
          return 116;
      }
      break;

    case 'M':
      switch (low) {
        case 0:
          return 7; 
        case 'C':
          return 115;
        case 'D':
          return 101;
        case 'G':
          return 12;
        case 'N':
          return 25;
        case 'O':
          return 42;
        case 'T':
          return 109;
      }
      break;

    case 'N':
      switch (low) {
        case 0:
          return 7; 
        case 'A':
          return 11;
        case 'B':
          return 41;
        case 'D':
          return 60;
        case 'E':
          return 10;
        case 'H':
          return 113;
        case 'I':
          return 28;
        case 'O':
          return 102;
        case 'P':
          return 93;
      }
      break; 


    case 'O':
      switch (low) {
        case 0:
          return 8; 
        case 'G':
          return 118;
        case 'S':
          return 76;
      }
      break;

    case 'P':
      switch (low) {
        case 0:
          return 15; 
        case 'A':
          return 91;
        case 'B':
          return 82;
        case 'D':
          return 46;
        case 'M':
          return 61;
        case 'O':
          return 84;
        case 'R':
          return 59;
        case 'T':
          return 78;
        case 'U':
          return 94;
      }
      break;
    
    case 'Q':
      return 8; 

    case 'R':
      switch (low) {
        case 'A':
          return 88;
        case 'B':
          return 37;
        case 'E':
          return 75;
        case 'F':
          return 104;
        case 'G':
          return 111;
        case 'H':
          return 45;
        case 'N':
          return 86;
        case 'U':
          return 44;
      }
      break;

    case 'S':
      switch (low) {
        case 0: 
          return 16; 
        case 'B':
          return 51;
        case 'C':
          return 21;
        case 'E':
          return 34;
        case 'G':
          return 106;
        case 'I':
          return 14;
        case 'M':
          return 62;
        case 'N':
          return 50;
        case 'R':
          return 38;
      }
      break;

    case 'T':
      switch (low) {
        case 'A':
          return 73;
        case 'B':
          return 65;
        case 'C':
          return 43;
        case 'E':
          return 52;
        case 'H':
          return 90;
        case 'I':
          return 22;
        case 'L':
          return 81;
        case 'M':
          return 69;
        case 'S':
          return 117;
      }
      break;

    case 'U':
      if(low == 'R')
        return 92;
      break;

    case 'V':
      if(low == 'A')
        return 23;
      break;
  
    case 'W':
      if(low == 'T')
        return 74;
      break; 

    case 'X':
      if (!low)
        return 6; 
      else if (low == 'E')
        return 54;
      break;

    case 'Y':
      if (!low)
        return 6; 
      else if (low == 'T')
        return 39;
      else if (low == 'B')
        return 70;
      break;

    case 'Z':
      if (!low)
        return 7; 
      else if (low == 'N')
        return 30;
      else if (low == 'R')
        return 40;
      break;
  }
  return 0;
}


typedef struct {
  symbol_t* c; 
  u8 order; 
} edge_t; 

struct symbol_t {
  u8 atomic_num;   
  u8 valence_pack;  // [  max u4    ][ curr     u4 ]  [0-8][0-8] 
  u8 n_bonds; 
  edge_t bonds[MAX_DEGREE]; // directional (most memory usage)
};

typedef struct {
  symbol_t *s; 
  u8 hloc; // used for pathsolver 
  u8 r_pack; // [ (of) 2b ][ arom 1b ][ bridging 1b ][ dangling u4 ] 
  edge_t *off_path[2]; // 0 = -. 1 = '&' 
} locant; 

typedef struct {
  u8  size; 
  locant path[1]; // malloc sizeof(locant) *size-1 + (1 byte for size)
} ring_t;  

/* memory pool - handle and reuse allocations */
typedef struct {
  u16 s_num; 
  u16 s_max; 
  symbol_t  *symbols; 
} graph_t; 


static void alloc_graph_t(graph_t *g, const size_t size)
{
  g->s_num    = 0; 
  g->s_max    = size; 
  g->symbols = (symbol_t*)malloc(sizeof(symbol_t) * size);  
  memset(g->symbols,0,sizeof(symbol_t) * size); 
}

static void realloc_graph_t(graph_t *g, const size_t size)
{
  g->symbols = (symbol_t*)realloc(g->symbols,sizeof(symbol_t) * (g->s_max + size) );  
  memset(g->symbols+g->s_max,0,sizeof(symbol_t) * size); 
  g->s_max += size; 
}

static void free_graph_t(graph_t *g)
{
  free(g->symbols); 
  memset(g,0, (sizeof(u16) * 4) + sizeof(symbol_t*)); 
}

/* 
 * overwrite and add have to use the same func signiture to avoid branch 
 * prediction when using fn_ptr lookup
 * */
static symbol_t* overwrite_symbol(graph_t *g, symbol_t *s, const u16 id, const u8 lim_valence)
{
  s->atomic_num   = id; 
  s->n_bonds      = 0;
  s->valence_pack = lim_valence; 
  s->valence_pack <<= 4;
  memset(s->bonds,0,sizeof(edge_t) * MAX_DEGREE); 
  return s; 
}

static symbol_t* add_symbol(graph_t *g, symbol_t *c, const u16 id, const u8 lim_valence)
{
  if(g->s_num == g->s_max)
    realloc_graph_t(g,g->s_max + REASONABLE); 
  
  symbol_t *s = &g->symbols[g->s_num++];
  s->atomic_num   = id; 
  s->n_bonds      = 0;
  s->valence_pack = lim_valence; 
  s->valence_pack <<= 4; 
  return s; 
}

static __always_inline edge_t* next_virtual_edge(symbol_t *p)
{
  edge_t *e = &p->bonds[p->n_bonds++]; 
  e->order += (e->order==0); 
  return e; 
}

/* if the edge is vitual, spawn it in by modifying the packing
 * for the child, packing for parent is modified on virtual spawn
 * but only needs checked when a real edge is made. 
 */
static edge_t* set_edge(edge_t *e, symbol_t *p, symbol_t *c)
{ 
  e->c = c; 
  c->valence_pack += e->order; 

  // TODO - nibble bit trick is definitely possible
  if ((p->valence_pack & 0x0F) >= (p->valence_pack >> 4) || 
      (c->valence_pack & 0x0F) >= (c->valence_pack >> 4)) 
  {
    fprintf(stderr,"Error: symbol reached WLN allowed valence - %d/%d & %d/%d\n",
            p->valence_pack & 0x0F, p->valence_pack >> 4,
            c->valence_pack & 0x0F, c->valence_pack >> 4); 
    return (edge_t*)0; 
  }
  return e; 
}


typedef struct  {
  u8 r_loc; // use to calculate size + add in off-branch positions
  u8 r_size;
  u8 arom; 
} r_assignment; 


#if 0
static int path_solverIII(graph_t *g, ring_t *r, 
                          r_assignment *SSSR, u8 SSSR_ptr, 
                          u8 state_pseudo) 
{
  /*
   * PathsolverIII - solving WLN Hamiltonian Paths 
   *
   * 2. When lookback locants /XX are used, path property is 
   *    broken and a flood fill required. 
   */
   
  edge_t *e; 
  locant *l;
  locant *start, *end; 
  u8 steps;
  r_assignment *subcycle; 


  // create the initial chain
  l = &r->path[0]; 
  if (!l->s) {
    l->s = add_symbol(g, 6, 4); 
    l->hloc = 1; 
    l->r_pack += 0x2; 
  }

  for (u16 i=1; i<r->size-1; i++) {
    l = &r->path[i]; 
    if (!l->s) { // edges might be made from unsaturates
      l->s = add_symbol(g, 6, 4); 
      e = next_virtual_edge(r->path[i-1].s); 
      e = set_edge(e, r->path[i-1].s, l->s); 
    }
    r->path[i-1].hloc = i; 
  }

  if (!state_pseudo) {
    // 1. WLN locants maximise path traversal, which is a flipped
    //    version of mimising the function sum as specified. 
    for (u16 i=0; i<SSSR_ptr; i++) {
      subcycle = &SSSR[i];   
      steps    = subcycle->r_size; 
      start    = &r->path[subcycle->r_loc]; 
      
      for (u16 s=0; s<steps; s++)
        end = &r->path[start->hloc]; 
      
      e = next_virtual_edge(start->s); 
      e = set_edge(e, start->s, end->s); 
      start->hloc = end - &r->path[0]; 
    }

  }
  else { 



  }


  return 1;   
}
#endif

static ring_t* parse_cyclic(const char *s_ptr, const char *e_ptr, graph_t *g) 
{
  ring_t *ring; 

  u8 locant_ch  = 0; 
  u8 arom_count = 0; 

  u8 SSSR_ptr  = 0; 
  r_assignment SSSR[REASONABLE]; // this really is upper bound sensible for WLN

  u8 dash_ptr = 0; 
  unsigned char dash_chars[3]; // last byte is mainly for overflow 
                               
  u16 state = 0; // bit field:                      
                 // [][][][][][][][][][][][][][dash][space][rings read]

  // note: The WLN ring system will always be <= than the sum of the ring assignment
  //       values, therefore alloc'ing the sum wastes the least amount of space and 
  //       avoids a hash map for bridging/pseudo bridging flags. 
  
  // ring_t *ring = (ring_t*)malloc((sizeof(locant)*REASONABLE-1)+1);  // realloc if needed
  
  // note: WLN has a lot of ambiguty due to the limited char set, therefore
  //       you have to reduce possible states in order to correctly handle. 
  //       Reading arom assignments first allows determined state of '& symbols. 

  while (*e_ptr == '&' || *e_ptr == 'T') {
    arom_count++; 
    e_ptr--; 
  }
  e_ptr++; 

  for (u16 i=0; i<arom_count; i++) 
    SSSR[i].arom = (*(e_ptr+i) == '&'); // & symbols can now only be a position expansion


  unsigned char ch = 1; 
  while (s_ptr != e_ptr) {
    ch = *(s_ptr++); 
    switch (ch) {
      case '0':
        break;

      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        SSSR[SSSR_ptr].r_size = ch - '0'; 
        SSSR[SSSR_ptr++].r_loc = ((locant_ch == 0) * 'A') + locant_ch; // branchless check
        locant_ch = 0; 
        break; 

      // locant only symbols within ring
      case 'A':
      case 'C':
      case 'D':
      case 'J':
      case 'L':
      case 'R':
        if (1) {
          locant_ch = ch; 
        }
        break; 



      case 'T':
        break;

      case '&':
        if (locant_ch)
          locant_ch += 23; 
        else {
          fprintf(stderr, "Error: & expansion used without previous locant\n"); 
          return (ring_t*)0; 
        }
        break; 

      case ' ':
        break; 

    }

  }
  
  return ring; 
}

/* assumes the head node contains only virtual bonds */
static symbol_t *add_alkyl_chain(graph_t *g, symbol_t *p, int size)
{
  edge_t *e;
  symbol_t *c=p;  
  for (u16 i=0; i<size; i++) {
    p->valence_pack++; 
    e = next_virtual_edge(c); 
    c = add_symbol(g, c, CARBON, 4);
    e = set_edge(e, p, c); 
    p = c; 
  }
  return c;
}

static void default_methyls(graph_t *g, symbol_t *c, const u8 n)
{
  edge_t *e; 
  symbol_t *m; 
  for (u8 i=(c->valence_pack & 0x0F); i<n; i++) {
    c->valence_pack++; 
    m = add_symbol(g, m, CARBON, 4); 
    e = next_virtual_edge(c); 
    e = set_edge(e, c, m); 
  }
  c->n_bonds = 0;  
}


/*
 * -- Parse WLN Notation --
 *
 */
static int parse_wln(const char *ptr, graph_t *g)
{
  edge_t   *e=0; 
  symbol_t *c=0;
  symbol_t *p=0;
  ring_t   *r=0;

  u8 locant_ch  = 0; 
  u8 alkyl_len  = 0; 
  u8 ring_chars = 0; 
  
  u8 state = 0; // bit field: ordering allows quick U checks
                // [0][0][0][0] [0][ring skip][digit][space]

  u16 stack_ptr = 0; 
  struct stack_frame {
    union f_addr{ // safter than void* cast. 
      symbol_t *s; 
      ring_t   *r; 
    } addr; 
    signed char ref; // -1 for (ring_t*) else (symbol_t*) 
  } stack[REASONABLE]; 
  
  // avoids a branch on each symbol case
  symbol_t* (*sym_fnPtr[2])(graph_t*, symbol_t*, const u16, const u8) = {&overwrite_symbol, &add_symbol}; 

  unsigned char ch = 1; 
  while(ch) {
    ch = *(ptr++);
    switch (ch) {
      case '0':
        if (state < DIGIT_READ) {
          fprintf(stderr,"Error: zero numeral without prefix digits\n"); 
          return 0; 
        }
        else
          alkyl_len *= 10; 
        break;

      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        if (state == RING_READ)
          ring_chars++; 
        else {
          alkyl_len *= 10; 
          alkyl_len += ch - '0'; 
          if (*ptr < '0' || *ptr > '9') {  // ptr is a +1 lookahead

            if (e) {
              p->valence_pack += (e->c == 0); 
              c = sym_fnPtr[e->c==0](g, e->c, CARBON, 4); 
              e = set_edge(e, p, c); 
            }
            else 
              c = add_symbol(g, 0, CARBON, 4); 
            
            c = add_alkyl_chain(g, c, alkyl_len-1);  
            e = next_virtual_edge(c); 
            alkyl_len = 0; 
            p = c; 
          }
        }
        break;
      
      case 'A':
      case 'B':
      case 'C':
      case 'D':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
        break; 
     
      case 'J':
        if (state == RING_READ) {
        
          if (*ptr == '&' || *ptr == ' ' || *ptr == 0) {
            // J can be used inside ring notation, requires lookahead 
            // condition 
            
            // note: the ptr passed in does not include the starting L/T or ending J
            r = parse_cyclic(ptr-ring_chars, ptr-2, g);
            if (!r) 
              return 0; 
            else {
              stack[stack_ptr].addr.r = r; 
              stack[stack_ptr++].ref = -1; 
              state &= ~(0x04);
              ring_chars = 0; 
            }
          }
          else
            ring_chars++; 
        }

        break; 

      
      case 'L':
      case 'T':
        if (state == RING_READ) {
          ring_chars++; 
        }
        break; 
     
#if 0
      /* nitrogen symbols */
      case 'N':
        if (state == RING_READ) 
          ring_chars++; 
        else {
          c = add_symbol(g, 7, 3); 
          stack[stack_ptr].addr.s = c; 
          stack[stack_ptr++].ref = 2;
          if (p) {
            e = set_edge(p, c); 
            state &= ~(0x30); 
          }
          p = c; 
        }
        break;
#endif

      case 'Y':
        if (state == RING_READ) 
          ring_chars++; 
        else {
          
          if (e) {
            p->valence_pack += (e->c == 0); 
            c = sym_fnPtr[(e->c==0)](g, e->c, CARBON, 4); 
            e = set_edge(e, p, c); 
          }
          else 
            c = add_symbol(g, c, CARBON, 4); 
          
          default_methyls(g, c, 3);  

          stack[stack_ptr].addr.s = c; 
          stack[stack_ptr].ref  = 2;
          stack_ptr++; 
          
          e = next_virtual_edge(c); 
          p = c; 
        }
        break;
      

      case ' ':
        if (state == RING_READ)
          ring_chars++; 
        break; 

      case '&':
        if (state == RING_READ) 
          ring_chars++; 
        else {
          if (!stack_ptr) {
            fprintf(stderr,"Error: backtracking stack is empty - too many &\n");
            return 0; 
          }
          else if (stack[stack_ptr-1].ref == -1) {
            // ring closures  
            free(stack[--stack_ptr].addr.r);
            stack[stack_ptr].addr.r = 0; 
            stack[stack_ptr].ref = 0; 
          }
          else {
            // stack logic for implied closures
            while (stack_ptr) { 
              c = stack[stack_ptr-1].addr.s; 
              if (stack[stack_ptr-1].ref == 0) 
                stack_ptr--; 
              else
                break;
            }
            
            if (!stack_ptr) {
              fprintf(stderr,"Error: backtrack stack is empty - too many &\n");
              return 0; 
            }
            else { 
              stack[stack_ptr-1].ref--; 
              p = c;
            }
          }
        }
        break;

      case 'U':
        if (state == 0x04) 
          ring_chars++;
        else if (e) {
          e->order += 1; 
          p->valence_pack++; 
        }
        else {
          fprintf(stderr,"Error: unsaturation called without previous bond\n");
          return 0;
        }
        break; 

      case 0:
        return 1; // allows the +1 look-ahead 

      default:
        fprintf(stderr,"Error: invalid character - %c (%d) for WLN notation\n",ch,ch);
        return 0; 
    }
  }
  
  return 1; 
}


OpenBabel::OBAtom* ob_add_atom(OpenBabel::OBMol* mol, u16 elem, char charge, char hcount)
{
  OpenBabel::OBAtom* result = mol->NewAtom();
  if(!result)
    return 0;

  result->SetAtomicNum(elem);
  result->SetFormalCharge(charge);
  if(hcount >= 0)
    result->SetImplicitHCount(hcount);
  return result;
}


OpenBabel::OBBond* ob_add_bond(OpenBabel::OBMol* mol, OpenBabel::OBAtom* s, OpenBabel::OBAtom* e, u8 order)
{
  OpenBabel::OBBond* bptr = 0; 
  if (!s || !e) {
    fprintf(stderr,"Error: could not find atoms in bond, bond creation impossible\n");
    return bptr;
  }

  if (!mol->AddBond(s->GetIdx(), e->GetIdx(), order)) {
    fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n",s->GetIdx(),e->GetIdx());
    return bptr;
  }
      
  bptr = mol->GetBond(mol->NumBonds() - 1);
  return bptr;
}


int ob_convert_wln_graph(OpenBabel::OBMol *mol, graph_t *g) {  
  for (u16 i=0; i<g->s_num; i++) {
    symbol_t *node = &g->symbols[i]; 
    switch (node->atomic_num) {
      case 6: // carbons
        ob_add_atom(mol, node->atomic_num, 0, 4 - (node->valence_pack & 0x0F)); 
        break; 
      
      case 7: // nitrogens
        ob_add_atom(mol, node->atomic_num, 0, 3 - (node->valence_pack & 0x0F)); 
        break; 
        
      default:
        break;
        //ob_add_atom(mol, atomic_num, 0, (node->max_valence-node->valence)); 
    }
  }
  
  /* will have a 1-1 indexing */
  for (u16 i=0; i<g->s_num; i++) {
    symbol_t *node = &g->symbols[i];
    for (u16 j=0; j<MAX_DEGREE; j++) {
      edge_t *e = &node->bonds[j];
      if (e->c) {
        u16 beg = node - g->symbols; 
        u16 end = e->c - g->symbols;
        ob_add_bond(mol, mol->GetAtom(beg+1), mol->GetAtom(end+1), e->order); 
      }
    }
  }
  return 1; 
}


int C_ReadWLN(const char *ptr, OpenBabel::OBMol* mol)
{   
  graph_t wln_graph; 
  graph_t *g = &wln_graph; 
  alloc_graph_t(g, REASONABLE); 
  
  if (!parse_wln(ptr, g))
    return 0; 

  // graph_to_dotfile(stderr, g); 
  ob_convert_wln_graph(mol,g); 

  free_graph_t(g); 
  return 1; 
}

