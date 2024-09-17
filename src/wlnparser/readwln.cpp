/*********************************************************************
 
Author : Michael Blakey
Description: WLN reader - write out SMILES etc from WLN

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

#define SPACE_READ  0x01 
#define DIGIT_READ  0x02 
#define DASH_READ   0x04

#define RING_READ   0x08
#define BIND_READ   0x10

#define SSSR_READ   0x04

// Error codes 
#define ERR_NONE   0 // success
#define ERR_ABORT  1 
#define ERR_MEMORY 2 

// Element "magic numbers"
#define DUMMY   0
#define CARBON  6
#define NITRO   7
#define OXYGEN  8

typedef struct symbol_t symbol_t; 

static u8 error(const char *message) 
{
  fprintf(stderr,"%s\n", message); 
  return ERR_ABORT;  
}

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


// 9 bytes
typedef struct {
  symbol_t* c; 
  u8 order; 
} edge_t; 

// 3 + (8*9) = 75 bytes per symbol
// 80 = 16*5, multiple of 16 stack aligned
struct symbol_t {
  u8 atomic_num;   
  u8 valence_pack;  // [  max u4    ][ curr     u4 ]  [0-8][0-8] 
  u8 n_bonds; 
  edge_t bonds[MAX_DEGREE]; // directional 
};

typedef struct {
  symbol_t *s; 
  u8 hloc; // used for pathsolver 
  u8 r_pack; // [ (of) 2b ][ arom 1b ][ bridging 1b ][ dangling u4 ] 
  edge_t *off_path[2]; // 0 = -. 1 = '&' 
} locant; 

typedef struct {
  u8  size; 
  locant path[1]; // malloc sizeof(locant) * (size-1) + (1 byte for size)
} ring_t;  


typedef struct {
  u16 s_num; 
  u16 s_max; 

  u8 stack_ptr;  // WLN DFS style branch and ring stack
  struct {
    void *addr; 
    signed char ref; // -1 for (ring_t*) else (symbol_t*) 
  } stack[32];  

  symbol_t  *symbols; 
} graph_t; 


static void gt_alloc(graph_t *g, const size_t size)
{
  g->s_num      = 0; 
  g->s_max      = size; 
  g->stack_ptr  = 0; 
  g->symbols    = (symbol_t*)malloc(sizeof(symbol_t) * size);  
  memset(g->symbols, 0, sizeof(symbol_t) * size); 
}

// seperated in order to reuse the symbols memory on file read
static void gt_stack_flush(graph_t *g) 
{
  u16 stack_ptr = g->stack_ptr; 
  for (u16 i=0;i<stack_ptr;i++) {
    if (g->stack[i].ref < 0) {
      free(g->stack[i].addr);
      g->stack[i].addr = 0; 
      g->stack[i].ref = 0;
    }
  }
}

static void gt_clear(graph_t *g)
{
  gt_stack_flush(g); 
  g->stack_ptr = 0; 
  g->s_num     = 0; 
}

static void gt_free(graph_t *g)
{
  gt_stack_flush(g); 
  g->stack_ptr = 0; 
  g->s_max     = 0; 
  g->s_num     = 0; 
  free(g->symbols); 
}


static symbol_t* next_symbol(graph_t *g, edge_t *e, const u16 id, const u8 lim_valence)
{
  if (g->s_num == g->s_max) {
    fprintf(stderr,"Warning: symbol limit reached (%d), reallocating...\n", g->s_max); 
    return (symbol_t*)0; 
  }
  
  symbol_t *s = 0; 
  
  if (!e->c)
    s = &g->symbols[g->s_num++];
  else 
    s = e->c; 

  s->atomic_num   = id; 
  s->n_bonds      = 0;
  s->valence_pack = lim_valence; 
  s->valence_pack <<= 4;
  memset(s->bonds,0,sizeof(edge_t) * MAX_DEGREE); 
  return s; 
}

static symbol_t* new_symbol(graph_t *g, const u16 id, const u8 lim_valence)
{
  if (g->s_num == g->s_max) {
    fprintf(stderr,"Warning: symbol limit reached, reallocating...\n"); 
    return (symbol_t*)0; 
  }
  
  symbol_t *s = &g->symbols[g->s_num++];

  s->atomic_num   = id; 
  s->n_bonds      = 0;
  s->valence_pack = lim_valence; 
  s->valence_pack <<= 4;
  memset(s->bonds,0,sizeof(edge_t) * MAX_DEGREE); 
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
static edge_t* set_virtual_edge(edge_t *e, symbol_t *p, symbol_t *c)
{
  p->valence_pack += (e->c == 0); // unsaturations modify this directly
  e->c = c; 
  c->valence_pack += e->order; 

  // TODO - nibble bit trick is definitely possible
  if ((p->valence_pack & 0x0F) > (p->valence_pack >> 4) || 
      (c->valence_pack & 0x0F) > (c->valence_pack >> 4)) 
  {
    fprintf(stderr,"Error: symbol reached WLN allowed valence - %d/%d & %d/%d\n",
            p->valence_pack & 0x0F, p->valence_pack >> 4,
            c->valence_pack & 0x0F, c->valence_pack >> 4); 
    return (edge_t*)0; 
  }
  return e; 
}


static void read_stack_frame(symbol_t **p, edge_t **e, ring_t **r, graph_t *g)
{
  u16 stack_top = g->stack_ptr - 1; 
  if (g->stack[stack_top].ref > 0) {
    *p = (symbol_t*)g->stack[stack_top].addr; 
    *e = next_virtual_edge(*p); 
  }
  else {
    *p = 0;
    *e = 0; 
    *r = (ring_t*)g->stack[stack_top].addr; 
  }
}


typedef struct  {
  u8 r_loc; // use to calculate size + add in off-branch positions
  u8 r_size;
  u8 arom; 
} r_assignment; 


static ring_t* path_solverIII(graph_t *g, ring_t *r, 
                              r_assignment *SSSR, u8 SSSR_ptr, 
                              u8 state_pseudo) 
{
  u8 steps;
  edge_t *e; 
  locant *c, *p=0;
  locant *start, *end; 
  r_assignment *subcycle; 

  for (u16 i=0; i<r->size; i++) {
    c = &r->path[i]; 
    if (!c->s) { 
      c->r_pack = 0x1 + (!i | (i==r->size-1)); 
      c->s = new_symbol(g, CARBON, 4); 

      if (p) {
        e = next_virtual_edge(p->s); 
        e = set_virtual_edge(e, p->s, c->s); 
      }
    }
    c->hloc = i+1; // pointing to the value in front 
    p = c; 
  }
 
  /*
   * PathsolverIII Algorithm:
   *
   * Named after the original attempts at WLN hamiltonian paths
   * from lynch et al. PathsolverIII iterates a given hamiltonian 
   * path by using the "allowed connections" property. Please
   * refer to my thesis for more details. 
   *
   * In short - ring bonds can have a maximum of 3 connections
   *            unless specified as bridging (-1) or expanded (+1)
   *
   * The path is maximised at each step which mirrors the minimisation
   * of the fusion sum as mentioned in the manuals. 
   */

  if (!state_pseudo) {

    for (u16 i=0; i<SSSR_ptr; i++) {
      subcycle = &SSSR[i];   
      steps    = subcycle->r_size; 
      start    = &r->path[subcycle->r_loc]; 
      end      = start; 

      for (u16 s=0; s<steps-1; s++)
        end = &r->path[end->hloc]; 
      
      e = next_virtual_edge(start->s); 
      e = set_virtual_edge(e, start->s, end->s);

      start->hloc = end - &r->path[0]; 
    }
  }
  else { 
    // 2. when pseudo bonds are specified, a flood fill is needed to solve
    // the path. Significantly slower therefore prefer the upper branch


  }

  return r;   
}

static ring_t* parse_cyclic(const char *s_ptr, const char *e_ptr, graph_t *g) 
{
  // note: The WLN ring system will always be <= than the sum of 
  //       the ring assignment values, therefore alloc'ing the sum 
  //       wastes the least amount of space and avoids a hash map for 
  //       bridging/pseudo bridging flags. 
  
  ring_t *ring; 

  u8 locant_ch  = 0; 
  u8 arom_count = 0; 
  u8 upper_r_size = 0; 

  u8 SSSR_ptr  = 0; 
  r_assignment SSSR[32]; // this really is sensible for WLN

  u8 dash_ptr = 0; 
  unsigned char dash_chars[3]; // last byte is mainly for overflow 
                               
  u8 state = 0; // bit field:                      
                 // [][][][][dash][SSSR][digit][space]

  
  // note: WLN has a lot of ambiguty due to the limited char set, therefore
  //       you have to reduce possible states in order to correctly handle. 
  //       Reading arom assignments first allows determined state of '& symbols. 
  while (*e_ptr == '&' || *e_ptr == 'T') {
    arom_count++; 
    e_ptr--; 
  }
  e_ptr++; 

  // & symbols can now only be a position expansion
  for (u16 i=0; i<arom_count; i++) 
    SSSR[i].arom = (*(e_ptr+i) == '&'); 

  unsigned char ch = 1; 
  state = SSSR_READ; 
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
        if (state == SSSR_READ){
          upper_r_size += ch - '0'; 
          SSSR[SSSR_ptr].r_size = ch - '0'; 
          SSSR[SSSR_ptr++].r_loc = locant_ch;  
          locant_ch = 0; 
        }
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
  
  // simplest rings
  if (state == SSSR_READ) {
    ring = (ring_t*)malloc(sizeof(ring_t) + sizeof(locant)*(upper_r_size-1)); 
    memset(ring, 0, sizeof(ring_t) + sizeof(locant)*(upper_r_size-1)); 
  }
 
  // do a mem allocation check here

  ring->size = upper_r_size; 
  return path_solverIII(g, ring, SSSR, SSSR_ptr, 0);     
}


static void default_methyls(graph_t *g, symbol_t *c, const u8 n)
{
  edge_t *e; 
  symbol_t *m; 
  for (u8 i=(c->valence_pack & 0x0F); i<n; i++) {
    e = next_virtual_edge(c); 
    m = next_symbol(g, e, CARBON, 4); 
    e = set_virtual_edge(e, c, m); 
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
  u8 ring_chars = 0; 
  u8 atom_num   = 0; 
  u16 alkyl_len = 0; 
  
  u8 state = 0; // bit field: 
                // [0][0][0][inline state][ring skip][dash][digit][space]
                //
                // bind_prev indicates either a inline ring or spiro that
                // must be bound to the previous chain
  
  u8 dash_ptr = 0; 
  unsigned char dash_chars[3]; // last byte is mainly for overflow 
 
  // init conditions, make one dummy atom, and one bond - work of the virtual bond
  // idea entirely *--> grow...
  
  c = new_symbol(g, DUMMY, 1); 
  e = next_virtual_edge(c); 
  e->order = DUMMY; 
  p = c; 
  
  // charges are assigned through string indexing - tf: need a strlen() 
    
  unsigned char ch; 
  unsigned char ch_nxt; 
  u16 len = strlen(ptr); 

  for (u16 sp=0; sp<len; sp++) {
    ch     = ptr[sp]; 
    ch_nxt = ptr[sp+1]; // one lookahead is defined behaviour 
    switch (ch) {
      case '0':
        if (!(state & DIGIT_READ))
          return error("Error: zero numeral without prefix digits"); 
        else {
          alkyl_len *= 10; 
          if (ch_nxt < '0' || ch_nxt > '9') {  
                       
            for (u16 i=0; i<alkyl_len; i++) {
              c = next_symbol(g, e, CARBON, 4);
              if (!c)
                return ERR_MEMORY; 
              else
                e = set_virtual_edge(e, p, c); 

              if (!e)
                return ERR_ABORT; 
              else {
                p = c;
                e = next_virtual_edge(p); 
              }
            }

            alkyl_len = 0; 
          }
        }
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
        if (state & RING_READ)
          ring_chars++; 
        else {
          alkyl_len *= 10; 
          alkyl_len += ch - '0'; 
          if (ch_nxt < '0' || ch_nxt > '9') { 
            
            for (u16 i=0; i<alkyl_len; i++) {
              c = next_symbol(g, e, CARBON, 4);
              if (!c)
                return ERR_MEMORY; 
              else
                e = set_virtual_edge(e, p, c); 

              if (!e)
                return ERR_ABORT; 

              p = c;
              e = next_virtual_edge(p); 
            }

            alkyl_len = 0; 
          }
          else
            state |= DIGIT_READ; 
        }
        break;
      
      case 'A':
      case 'B':
        if (state & RING_READ)
          ring_chars++; 
        else if (state & DASH_READ) {
          if (dash_ptr == 3) 
            return error("Error: elemental code can only have 2 character symbols"); 
          else
            dash_chars[dash_ptr++] = ch; 
        }
        else if (state & SPACE_READ)  {
          // switch on packing - state popcnt() >= hit bits. 
          locant_ch = ch - 'A'; 
          state &= ~(SPACE_READ); 
          
          if (state & BIND_READ)
            break; 
          else if (r && locant_ch < r->size) {
            c = r->path[locant_ch].s; 
            e = next_virtual_edge(c); 
            p = c; 
          } 
          else 
            return error("Error: out of bounds locant access"); 
        }
        break; 


      case 'C':
      case 'D':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
        if (state & RING_READ) 
          ring_chars++; 
        else if (state & DASH_READ) {
          if (dash_ptr == 3)
            return error("Error: elemental code can only have 2 character symbols"); 
          else
            dash_chars[dash_ptr++] = ch; 
        }
        break; 
     
      case 'J':
        if (state & RING_READ) {
        
          if (*ptr == '&' || *ptr == ' ' || *ptr == 0) {
            // J can be used inside ring notation, requires lookahead 
            // condition 
            
            // note (2 magic number): The ptr passed in does not include the 
            //                        starting L/T or ending J
            
            r = parse_cyclic(ptr-ring_chars, ptr-2, g);
            if (!r) 
              return 0; 
            else {
              g->stack[g->stack_ptr].addr = r; 
              g->stack[g->stack_ptr++].ref = -1; 
              state &= ~(RING_READ);
              ring_chars = 0; 
            }

            if (state == BIND_READ){
              // virutal edge should be dangling at this frame. 

              // this needs wrapping to avoid seg fault on a user error
              c = r->path[locant_ch].s;

              e = set_virtual_edge(e, p, c);  
              state &= ~(BIND_READ); 
            }
          }
          else
            ring_chars++; 
        }
        else if (state & DASH_READ) {
          if (dash_ptr == 3) 
            return error("Error: elemental code can only have 2 character symbols"); 
          else
            dash_chars[dash_ptr++] = ch; 
        }

        break; 

      
      case 'L':
      case 'T':
        if (state & RING_READ) 
          ring_chars++; 
        else if (state & DASH_READ) {
          if (dash_ptr == 3) 
            return error("Error: elemental code can only have 2 character symbols"); 
          else
            dash_chars[dash_ptr++] = ch; 
        }
        else {
          state |= RING_READ; 
          ring_chars++; 
        }
        break; 
     
      /* nitrogen symbols */
      case 'N':
        if (state & RING_READ) 
          ring_chars++; 
        else if (state & DASH_READ) {
          if (dash_ptr == 3) 
            return error("Error: elemental code can only have 2 character symbols"); 
          else
            dash_chars[dash_ptr++] = ch; 
        }
        else {
          c = next_symbol(g, e, NITRO, 4);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            g->stack[g->stack_ptr].addr = c; 
            g->stack[g->stack_ptr].ref  = 2;
            g->stack_ptr++; 

            p = c;
            e = next_virtual_edge(c); 
          }
        }
        break;

      case 'X':
        if (state & RING_READ) 
          ring_chars++; 
        else if (state & DASH_READ) {
          if (dash_ptr == 3) 
            return error("Error: elemental code can only have 2 character symbols"); 
          else
            dash_chars[dash_ptr++] = ch; 
        }
        else {
          c = next_symbol(g, e, CARBON, 4);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            default_methyls(g, c, 4);  

            g->stack[g->stack_ptr].addr = c; 
            g->stack[g->stack_ptr].ref  = 3;
            g->stack_ptr++; 

            p = c; 
            e = next_virtual_edge(c); 
          }
        }
        break;

      case 'Y':
        if (state & RING_READ) 
          ring_chars++; 
        else {
          c = next_symbol(g, e, CARBON, 4);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            default_methyls(g, c, 3);  

            g->stack[g->stack_ptr].addr = c; 
            g->stack[g->stack_ptr].ref  = 2;
            g->stack_ptr++; 

            p = c; 
            e = next_virtual_edge(c); 
          }
        }
        break;
      
      case 'Z':
        if (state & RING_READ) 
          ring_chars++; 
        else {
          c = next_symbol(g, e, NITRO, 1);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
          
            // terminating symbol
            if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                read_stack_frame(&p, &e, &r, g); 
            }
            else {
              p = c; 
              e = next_virtual_edge(p); 
            }
          }
        }
        break;

      case 'U':
        if (state & RING_READ) 
          ring_chars++;
        else if (e) {
          e->order += 1; 
          p->valence_pack++; 
        }
        else
          return error("Error: unsaturation called without previous bond"); 
        break; 
      
      case '-':
        if (state & RING_READ)
          ring_chars++; 
        else if (state & DASH_READ) {
          
          atom_num = get_atomic_num(dash_chars[0],dash_chars[1]); 
          if (!atom_num) 
            return error("Error: invalid element two character code"); 

          
          memset(dash_chars,0,3); 
          state &= ~DASH_READ; 
        }
        else {
          if (ch_nxt == ' ') 
            state |= BIND_READ; 
          else {
            dash_ptr = 0; 
            state |= DASH_READ; 
          }
        }
        break; 

      case ' ':
        if (state & RING_READ)
          ring_chars++; 
        else 
          state |= SPACE_READ; 
        break; 

      case '&':
        if (state & RING_READ) 
          ring_chars++; 
        else if (state & SPACE_READ) {
          // either a new ion, or a charge assignment 
          // - create a new dummy, handle charge with
          // alkyl values

          c = new_symbol(g, DUMMY, 1); 
          e = next_virtual_edge(c); 
          e->order = DUMMY; 
          p = c; 

          gt_stack_flush(g); 
          g->stack_ptr = 0; 

          state &= ~SPACE_READ; 
        }
        else if (!g->stack_ptr) 
          return error("Error: empty stack - too many &?"); 
        else { 

          if (g->stack[g->stack_ptr-1].ref < 0) {
            // ring closures  
            free(g->stack[--g->stack_ptr].addr);
            g->stack[g->stack_ptr].addr = 0; 
            g->stack[g->stack_ptr].ref = 0;
          }
          else {
            // branch closures
            // note: methyl contractions will have a live virtual 
            // bond, therefore can be used checked with AND. 

            c = (symbol_t*)g->stack[g->stack_ptr-1].addr; 
            g->stack_ptr -= (p == c) & (c->bonds[c->n_bonds].c == 0); 
            g->stack[g->stack_ptr-1].ref--; 
            g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 
          }

          // reseting block
          if (!g->stack_ptr) 
            return error("Error: empty stack - too many &?"); 
          else 
            read_stack_frame(&p, &e, &r, g); 
        }
        break;

      default:
        return error("Error: invalid character read for WLN notation"); 
    }
  }
  
  return ERR_NONE; 
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
    fprintf(stderr,"Error: could not find atoms in bond, bond creation impossible s: %p, e: %p\n",s,e);
    return bptr;
  }

  if (!mol->AddBond(s->GetIdx(), e->GetIdx(), order)) {
    fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n",s->GetIdx(),e->GetIdx());
    return bptr;
  }
      
  bptr = mol->GetBond(mol->NumBonds() - 1);
  return bptr;
}


/*
 * A dummy atom is created at new instance positions, siginificantly simplifies 
 * and speeds up the parsing code - removes a lot of branch checks, requires a mapping
 * to build mol
 */
int ob_convert_wln_graph(OpenBabel::OBMol *mol, graph_t *g) {  
  OpenBabel::OBAtom *atom; 
  OpenBabel::OBAtom **amapping = (OpenBabel::OBAtom **)alloca(sizeof(OpenBabel::OBAtom*) * g->s_num); 
  for (u16 i=1; i<g->s_num; i++) {
    symbol_t *node = &g->symbols[i]; 
    
    // transition metals complicate this somewhat
    
    switch (node->atomic_num) {
      case CARBON: 
        atom = ob_add_atom(mol, node->atomic_num, 0, 4 - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
        break; 

      case NITRO:
        atom = ob_add_atom(mol, node->atomic_num, 0, 3 - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
        break; 
      
      case DUMMY: // used to simplify grow code
        break; 

      default:
        atom = ob_add_atom(mol, node->atomic_num, 0, ((node->valence_pack & 0xF0) >> 4) - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
    }
  }
  
  for (u16 i=1; i<g->s_num; i++) {
    symbol_t *node = &g->symbols[i];
    if (node->atomic_num != DUMMY) {
      for (u16 j=0; j<MAX_DEGREE; j++) {
        edge_t *e = &node->bonds[j];
        if (e->c) {
          u16 beg = i; 
          u16 end = e->c - g->symbols;
          ob_add_bond(mol, amapping[beg], amapping[end], e->order); 
        }
      }
    }
  }
  return 1; 
}


int C_ReadWLN(const char *ptr, OpenBabel::OBMol* mol)
{   
  int RET_CODE = 0; 
  u16 st_pool_size = 128; 
  graph_t wln_graph; 
  graph_t *g = &wln_graph; 
  gt_alloc(g, st_pool_size); 

  // allows recovery and resetting for super large WLN strings
  for (;;) {

    RET_CODE = parse_wln(ptr, g); 
    switch (RET_CODE) {
      
      case ERR_MEMORY: 
        if (st_pool_size == 1024) {
          fprintf(stderr,"Error: WLN string specifies > 1024 atoms\n"); 
          return 0; 
        }
        else {
          st_pool_size *= 2; 
          gt_free(g); 
          gt_alloc(g, st_pool_size); 
        }
        break; 

      case ERR_ABORT:
        gt_free(g); 
        return 0; 
      
      default:
        ob_convert_wln_graph(mol,g); 
        gt_free(g); 
        return 1; 
    }
  }
}

