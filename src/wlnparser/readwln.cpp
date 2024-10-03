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

/* TODO: 
 * 
*/

#define OPENBABEL 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if OPENBABEL
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>

#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#endif 

#include "parser.h"

#define DEBUG 1 // debug log - lvls: 0 - none, 1 - minimal, 2 - all

#define MAX_DEGREE 8

// common bit fields
#define SPACE_READ  0x01 
#define DIGIT_READ  0x02 
#define DASH_READ   0x04

// standard parse fields
#define RING_READ   0x08
#define BIND_READ   0x10
#define CHARGE_READ 0x20

// ring parse fields
#define SSSR_READ   0x08
#define MULTI_READ  0x10
#define PSEUDO_READ 0x20

// Error codes 
#define ERR_NONE   0 // success
#define ERR_ABORT  1 
#define ERR_MEMORY 2 

// Element "magic numbers"
#define DUM   0
#define BOR   5
#define CAR   6
#define NIT   7
#define OXY   8
#define FLU   9
#define PHO   15
#define SUL   16
#define CHL   17
#define BRO   35
#define IOD   53

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

// 4 + (8*9) = 76 bytes per symbol
// 80 = 16*5, multiple of 16 stack aligned
struct symbol_t {
  u8 atomic_num;   
  u8 charge;   
  u8 valence_pack;  // [  max u4    ][ curr     u4 ]  [0-8][0-8] 
  u8 n_bonds; 
  edge_t bonds[MAX_DEGREE]; // directional 
};

typedef struct locant {
  symbol_t *s; 
  u8 hloc; // used for pathsolver 
  u8 r_pack; // [ (of) 2b ][ arom 1b ][ bridging 1b ][ dangling u4 ] 
  locant *off_path[2]; // 0 = -. 1 = '&' 
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
  symbol_t  **idx_symbols; 
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
  s->charge       = 0; 
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
  s->charge       = 0; 
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

static ring_t* new_ring(const size_t size) 
{
  ring_t *ring = 0; 
  ring = (ring_t*)malloc(sizeof(ring_t) + sizeof(locant)*(size-1)); 
  memset(ring, 0, sizeof(ring_t) + sizeof(locant)*(size-1)); 
  return ring; 
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


static ring_t* pathsolverIII_fast(graph_t *g, ring_t *r, 
                                   r_assignment *SSSR, u8 SSSR_ptr) 
{
  u8 steps, s_pos;
  edge_t *e; 
  locant *c, *p=0;
  locant *start, *end; 
  r_assignment *subcycle; 

  for (u16 i=0; i<r->size; i++) {
    c = &r->path[i]; 
    if (!c->s) 
      c->s = new_symbol(g, CAR, 4); 

    if (p) {
      e = next_virtual_edge(p->s); 
      e = set_virtual_edge(e, p->s, c->s); 
    }

    c->r_pack++; 
    c->hloc = i+1; // point to the next locant 
    p = c; 
  }

  // end chain movement
  r->path[0].r_pack++; 
  r->path[r->size-1].r_pack++; 
  r->path[r->size-1].hloc = r->size-1; 

  /*
   * PathsolverIII FAST Algorithm (Michael Blakey):
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
   *
  */
  
  for (u16 i=0; i<SSSR_ptr; i++) {
    subcycle = &SSSR[i];   
    steps    = subcycle->r_size; 
    s_pos    = subcycle->r_loc; 
    start    = &r->path[s_pos]; 
    end      = start; 
  
    for (u16 s=0; s<steps-1; s++)
      end = &r->path[end->hloc]; 

    // if used max times in ring, shift along path
    while ((start->r_pack & 0x0F) == 0 && s_pos < r->size)
      start = &r->path[++s_pos];  

#if DEBUG 
    fprintf(stderr,"%d: %c --> %c\n",steps,start - &r->path[0] + 'A',end - &r->path[0] + 'A'); 
#endif
    
    start->r_pack--; 
    end->r_pack--; 

    e = next_virtual_edge(start->s); 
    e = set_virtual_edge(e, start->s, end->s);

    start->hloc = end - &r->path[0]; 
  }

  return r;   
}

static ring_t* pathsolverIII(graph_t *g, ring_t *r, 
                             r_assignment *SSSR, u8 SSSR_ptr, 
                             u8 *connection_table)
{
  u8 steps, s_pos;
  edge_t *e; 
  locant *c, *p=0;
  locant *start, *end; 
  r_assignment *subcycle; 

  for (u16 i=0; i<r->size; i++) {
    c = &r->path[i]; 
    if (!c->s) 
      c->s = new_symbol(g, CAR, 4); 

    if (p) {
      e = next_virtual_edge(p->s); 
      e = set_virtual_edge(e, p->s, c->s); 
      connection_table[i * r->size +(i-1)] = 1; 
      connection_table[(i-1) * r->size + i] = 1; 
    }

    connection_table[i*r->size+i] = 1; // bond to itself?

    c->r_pack++; 
    c->hloc = i+1; // point to the next locant 
    p = c; 
  }

  // end chain movement
  r->path[0].r_pack++; 
  r->path[r->size-1].r_pack++; 
  r->path[r->size-1].hloc = r->size-1; 

  /*
   * PathsolverIII Algorithm (Michael Blakey):
   *
   * Named after the original attempts at WLN hamiltonian paths
   * from lynch et al. PathsolverIII iterates a given hamiltonian 
   * path by using the "allowed connections" property. Please
   * refer to my thesis for more details. 
   *
   * Pseudo locants break the iterative walk, and a flood fill is required to 
   * find the maximal path through the ring system, this can be optimised somewhat
   * by using a priority queue on the walks, and bonding the pseudo positions
   * during the notation parse
   *
   * connection table allows the floodfill to be done without another data structure, plus a 
   * easy pass through for pseudo locants defined in the ring parse 
  */
  


  return r;  
}




static ring_t* parse_cyclic(const char *ptr, const u16 s, u16 e, graph_t *g) 
{
  
  symbol_t *c; 
  ring_t *ring; 

  u8 locant_ch     = 0; 
  u8 arom_count    = 0;
  u8 upper_r_size  = 0;  
  // note: The WLN ring system will always be <= than the sum of 
  //       the ring assignment values, therefore alloc'ing the sum 
  //       wastes the least amount of space and avoids a hash map for 
  //       bridging/pseudo bridging flags. 
  

  u8 SSSR_ptr  = 0; 
  r_assignment SSSR[32]; // this really is sensible for WLN

  u8 buff_ptr = 0; 
  unsigned char buffer[3]; 
    
  u8 *connection_table = 0; // stack allocates on pseudo read

  u8 state = 0; // bit field:                      
                 // [][][][pseudo][SSSR][dash][digit][space]

  unsigned char ch; 
  
  // note: WLN has a lot of ambiguty due to the limited char set, therefore
  //       you have to reduce possible states in order to correctly handle. 
  //       Reading arom assignments first allows determined state of '&' symbols. 

#ifdef AROM
  while (ptr[e] == '&' || ptr[e] == 'T') {
    SSSR[arom_count++].arom = (ptr[e] == '&'); // will require reversal 
    e--; 
  }
  e++; 
#endif

  state = SSSR_READ; 
  for (u16 sp=s; sp<e; sp++){
    ch     = ptr[sp]; 
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
        if (state & SPACE_READ) {
          // multicyclic block start
          // move forward space to skip redundant block - should appear on space, 
          // if not, error.
          //
          
          sp += ch - '0' + 1; 
          if (sp >= e || ptr[sp] != ' ') {
            fprintf(stderr, "Error: invalid format for multicyclic ring\n"); 
            return (ring_t*)0; 
          }
          else {
            state |= MULTI_READ; 
            state &= ~SSSR_READ; // turn off SSSR if present
            state &= ~SPACE_READ; 
          }
        }
        else if (state & SSSR_READ){
          upper_r_size += ch - '0'; 
          SSSR[SSSR_ptr].r_size = ch - '0'; 
          SSSR[SSSR_ptr++].r_loc = locant_ch;  
          locant_ch = 0; 
        }
        break; 
      
      
      case 'L':
      case 'M':
      case 'S':
      case 'A':
      case 'C':
      case 'D':
      case 'F':
      case 'J':
      case 'R':
      case 'T':
        if (state & PSEUDO_READ) {
          if (buff_ptr == 2) {
            fprintf(stderr,"Error: more than two pseudo bonds specified\n"); 
            return (ring_t*)0;
          }
          else 
            buffer[buff_ptr++] = ch; 
        }
        else if (state & SPACE_READ) {
          locant_ch = ch - 'A';
          state &= ~SPACE_READ; 
        } 
        break;
    
      case 'G':
      case 'N':
        if (state & PSEUDO_READ) {
          if (buff_ptr == 2) {
            fprintf(stderr,"Error: more than two pseudo bonds specified\n"); 
            return (ring_t*)0;
          }
          else 
            buffer[buff_ptr++] = ch; 
        }
        else if (state & SPACE_READ) {
          locant_ch = ch - 'A';
          state &= ~SPACE_READ; 
        } 
        else if (state & MULTI_READ) {
          // must be the incoming ring size
          ring->size = ch - 'A' + 1; 
        }
        else if (state & SSSR_READ) {
          // end the SSSR read, start element reading
          if (locant_ch > upper_r_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            return (ring_t*)0; 
          }
          else {
            ring = new_ring(upper_r_size); 
            c = ring->path[locant_ch].s = new_symbol(g, NIT, 3); 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; // allows sequetial assignment 
          }

          state &= ~SSSR_READ; 
        }
        else {
          // elemental assignment 
          if (locant_ch > upper_r_size) {
            fprintf(stderr,"Error: out of bounds locant access - %d\n",locant_ch); 
            return (ring_t*)0; 
          }
          else {
            c = ring->path[locant_ch].s = new_symbol(g, NIT, 3); 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; // allows sequetial assignment 
          }
        }
        break; 

      case '&':
        if (state & MULTI_READ)
          ring->size += 23; 
        else if (state & PSEUDO_READ)
          buffer[buff_ptr-1] += 23; 
        else if (locant_ch)
          locant_ch += 23; 
        else {
          fprintf(stderr, "Error: & expansion used without previous locant\n"); 
          return (ring_t*)0; 
        }
        break; 

      case '/':
        if (state & SSSR_READ) {
          ring = new_ring(upper_r_size); 
          state &= ~(SSSR_READ); 
          buff_ptr = 0; 
          state |= PSEUDO_READ; 
        }
        else if (state & PSEUDO_READ) {
          if (buff_ptr != 2) {
            fprintf(stderr,"Error: pseudo locants must come in pairs\n"); 
            return (ring_t*)0; 
          }


        }
        else {
          buff_ptr = 0; 
          state |= PSEUDO_READ; 
        }
        break; 

      case ' ':
        if (state & PSEUDO_READ) {
          // make the pseudo bond immediately avaliable, 
          if (buff_ptr != 2) {
            fprintf(stderr,"Error: pseudo locants must come in pairs\n"); 
            return (ring_t*)0; 
          }
          else {
#if DEBUG 
            fprintf(stderr,"pseudo: %c --> %c\n", buffer[0], buffer[1]); 
#endif 
            
            buff_ptr = 0; 
            memset(buffer,0,2); 

            state &= ~(PSEUDO_READ & SSSR_READ);  // no further rings can come from pseudo bonds
          }
        }

        state |= SPACE_READ; 
        locant_ch = 0; 
        break;

      default:
        break; 
    }
  }
  
  // simplest rings
  if ((state & SSSR_READ) | (state & MULTI_READ)) {
    ring = new_ring(upper_r_size); 
    ring->size = upper_r_size; 
  }
 
  
  if (connection_table)
    return pathsolverIII(g, ring, SSSR, SSSR_ptr, connection_table);
  else 
    return pathsolverIII_fast(g, ring, SSSR, SSSR_ptr);     
}


static u8 default_methyls(graph_t *g, symbol_t *c, const u8 n)
{
  edge_t *e; 
  symbol_t *m; 
  for (u8 i=(c->valence_pack & 0x0F); i<n; i++) {
    e = next_virtual_edge(c); 
    m = next_symbol(g, e, CAR, 4); 
    if (!m)
      return ERR_MEMORY; 
    else 
      e = set_virtual_edge(e, c, m); 
  }

  c->n_bonds = 0;  
  return ERR_NONE; 
}

static u8 add_oxy(graph_t *g, symbol_t *p)
{
  edge_t *e = next_virtual_edge(p); 
  symbol_t *c = next_symbol(g, e, OXY, 2); 

  if (!c)
    return ERR_MEMORY; 
  else {
    e->order++; 
    p->valence_pack++; 
    set_virtual_edge(e, p, c); 
  }
  return ERR_NONE; 
}

/* placeholder for charged alterations */ 
static u8 add_dioxy(graph_t *g, symbol_t *p)
{
  u8 ret = 0; 
  ret = add_oxy(g, p); 
  if (ret != ERR_NONE)
    return ret; 
  else 
    return add_oxy(g, p); 
}; 

/*
 * -- Parse WLN Notation --
 *
 */
static int parse_wln(const char *ptr, const u16 len, graph_t *g)
{
  edge_t   *e=0; 
  symbol_t *c=0;
  symbol_t *p=0;
  ring_t   *r=0;

  u8 locant_ch  = 0; 
  u8 ring_chars = 0; 
  u8 atom_num   = 0; 
  u16 digit_n = 0; 
  
  u8 state = 0; // bit field: 
                // [0][0][charge][inline ring][ring skip][dash][digit][space]
                //
                // bind_prev indicates either a inline ring or spiro that
                // must be bound to the previous chain
  
  u8 dash_ptr = 0; 
  unsigned char dash_chars[3]; // last byte is mainly for overflow 
 
  // init conditions, make one dummy atom, and one bond - work of the virtual bond
  // idea entirely *--> grow...
  
  c = new_symbol(g, DUM, 1); 
  e = next_virtual_edge(c); 
  e->order = DUM; 
  p = c; 
  g->idx_symbols[0] = p;  // overwrite dummy when 0 pos is requested
  
  // charges are assigned through string indexing - tf: need a strlen() 
    
  unsigned char ch; 
  unsigned char ch_nxt; 

  for (u16 sp=0; sp<len; sp++) {
    ch     = ptr[sp]; 
    ch_nxt = ptr[sp+1]; // one lookahead is defined behaviour 
    
    if (state & RING_READ) {
      if (ch == 'J' && ch_nxt < '0') {
        // J can be used inside ring notation, requires lookahead 
        // condition 
        
        // note: The ptr passed in does not include the 
        //       starting L/T or ending J (<sp) 
        
        r = parse_cyclic(ptr, sp-ring_chars+1, sp, g);
        if (!r) 
          return ERR_ABORT; 
        else {
          g->stack[g->stack_ptr].addr = r; 
          g->stack[g->stack_ptr++].ref = -1; 
          state &= ~(RING_READ);
          ring_chars = 0; 
        }

        if (state == BIND_READ){
          // virtual edge should be dangling at this frame. 
          if (locant_ch > r->size) {
            fprintf(stderr,"Error: out of bounds locant access"); 
            return ERR_ABORT; 
          }
          else 
            c = r->path[locant_ch].s;

          e = set_virtual_edge(e, p, c);  
          state &= ~(BIND_READ); 
        }
      }
      else
        ring_chars++; 
    }
    else if (state & DASH_READ && ch != '-') {
      if (dash_ptr == 3) 
        return error("Error: elemental code can only have 2 character symbols"); 
      else
        dash_chars[dash_ptr++] = ch; 
    }
    else if (state & SPACE_READ) {
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
    else {
      switch (ch) {

        case '0':
          digit_n *= 10; 
          if (ch_nxt == '/') {
            // positive charge assignment
            if (digit_n > len || !g->idx_symbols[digit_n]) {
              fprintf(stderr,"Error: charge assignment out of bounds\n");
              return ERR_ABORT; 
            }
            else {
              g->idx_symbols[digit_n]->charge++; 
              state |= CHARGE_READ; 
              digit_n = 0; 
              sp++; // skip the slash, only ever used again in R notation
            }
          } 
          else if ( state & (CHARGE_READ | DIGIT_READ) 
                    && (ch_nxt < '0' || ch_nxt > '9')) 
          {  
                       
            if (state & CHARGE_READ) {
              // negative charge assignment
              if (digit_n > len || !g->idx_symbols[digit_n]) {
                fprintf(stderr,"Error: charge assignment out of bounds\n");
                return ERR_ABORT; 
              }
              else {
                g->idx_symbols[digit_n]->charge--; 
                state &= ~CHARGE_READ; 
                digit_n = 0; 
              }

            }
            else {
              // create the alkyl chain
              for (u16 i=0; i<digit_n; i++) {
                c = next_symbol(g, e, CAR, 4);
                if (!c)
                  return ERR_MEMORY; 
                else
                  e = set_virtual_edge(e, p, c); 

                if (!e)
                  return ERR_ABORT; 

                p = c;
                e = next_virtual_edge(p); 
              }

              g->idx_symbols[sp+1] = p; 
              digit_n = 0; 
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
          digit_n *= 10; 
          digit_n += ch - '0'; 
          
          if (ch_nxt == '/') {
            // positive charge assignment
            if (digit_n > len || !g->idx_symbols[digit_n]) {
              fprintf(stderr,"Error: charge assignment out of bounds\n");
              return ERR_ABORT; 
            }
            else {
              g->idx_symbols[digit_n]->charge++; 
              state |= CHARGE_READ; 
              digit_n = 0; 
              sp++; // skip the slash, only ever used again in R notation
            }
          } 
          else if (ch_nxt < '0' || ch_nxt > '9') { 
            
            if (state & CHARGE_READ) {
              // negative charge assignment
              if (digit_n > len || !g->idx_symbols[digit_n]) {
                fprintf(stderr,"Error: charge assignment out of bounds\n");
                return ERR_ABORT; 
              }
              else {
                g->idx_symbols[digit_n]->charge--; 
                state &= ~CHARGE_READ; 
                digit_n = 0; 
              }

            }
            else {
              // create the alkyl chain
              for (u16 i=0; i<digit_n; i++) {
                c = next_symbol(g, e, CAR, 4);
                if (!c)
                  return ERR_MEMORY; 
                else
                  e = set_virtual_edge(e, p, c); 

                if (!e)
                  return ERR_ABORT; 

                p = c;
                e = next_virtual_edge(p); 
              }

              g->idx_symbols[sp+1] = p; 
              digit_n = 0; 
            }
          }
          else
            state |= DIGIT_READ; 
          break;
        
        case 'A':
        case 'J':
          return error("Error: non-atomic symbol used in chain"); 
        
        case 'B':
          c = next_symbol(g, e, BOR, 3);
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
            g->idx_symbols[sp+1] = p; 
            e = next_virtual_edge(c); 
          }
          break; 
        
        case 'C':
          fprintf(stderr,"Full unsaturate needs handling\n"); 
          break; 

        // extrememly rare open chelate notation
        case 'D':
          fprintf(stderr,"chelate ring open needs handling\n"); 
          break; 

        case 'E':
          c = next_symbol(g, e, BRO, 1);
          if (!c)
            return ERR_MEMORY; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

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
          break;

        case 'F':
          c = next_symbol(g, e, FLU, 1);
          if (!c)
            return ERR_MEMORY; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

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
          break;

        case 'G':
          c = next_symbol(g, e, CHL, 1);
          if (!c)
            return ERR_MEMORY; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

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
          break;

        case 'H':
          fprintf(stderr,"Explicit hydrogen needs implementation\n"); 
          break; 

        case 'I':
          c = next_symbol(g, e, IOD, 1);
          if (!c)
            return ERR_MEMORY; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

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
          break;
       
        case 'K':
          c = next_symbol(g, e, NIT, 4);
          c->charge++; 
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
            g->idx_symbols[sp+1] = c; 
            e = next_virtual_edge(c); 
          }
          break;

        case 'L':
        case 'T':
          state |= RING_READ; 
          ring_chars++; 
          break; 
        
        case 'M':
          c = next_symbol(g, e, NIT, 3);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            p = c;
            g->idx_symbols[sp+1] = p; 
            e = next_virtual_edge(c); 
          }
          break;
       
        /* nitrogen symbols */
        case 'N':
          c = next_symbol(g, e, NIT, 3);
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
            g->idx_symbols[sp+1] = p; 
            e = next_virtual_edge(c); 
          }
          break;
        
        case 'O':
          c = next_symbol(g, e, OXY, 2);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            p = c;
            g->idx_symbols[sp+1] = p; 
            e = next_virtual_edge(c); 
          }
          break;


        case 'Q':
          c = next_symbol(g, e, OXY, 1);
          if (!c)
            return ERR_MEMORY; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

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
          break;
        
        case 'U':
          if (e) {
            e->order++; 
            p->valence_pack++; 
          }
          else
            return error("Error: unsaturation called without previous bond"); 
          break; 

        case 'V':
          c = next_symbol(g, e, CAR, 4);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            if (add_oxy(g, c) == ERR_MEMORY)
              return ERR_MEMORY; 
            else {
              p = c; 
              g->idx_symbols[sp+1] = c; 
              e = next_virtual_edge(c); 
            }
          }
          break;
        
        // if not previous, create dummy carbon
        case 'W':
          c = next_symbol(g, e, CAR, 4);
          if (!c)
            return ERR_MEMORY; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            
            if (add_dioxy(g, c) == ERR_MEMORY)
              return ERR_MEMORY; 
            else {
              ; 
            }
          }
          break;

          break; 

        case 'X':
          c = next_symbol(g, e, CAR, 4);
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
            g->idx_symbols[sp+1] = c; 
            e = next_virtual_edge(c); 
          }
          break;

        case 'Y':
          c = next_symbol(g, e, CAR, 4);
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
          break;
        
        case 'Z':
          c = next_symbol(g, e, NIT, 1);
          if (!c)
            return ERR_MEMORY; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

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
          break;

        
        case '-':
          if (state & DASH_READ) {
            
            atom_num = get_atomic_num(dash_chars[0],dash_chars[1]); 
            if (!atom_num) 
              return error("Error: invalid element two character code"); 
            else {

            }
            
            memset(dash_chars,0,3); 
            state &= ~DASH_READ; 
            dash_ptr = 0; 
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
          state |= SPACE_READ; 
          break; 

        case '&':
          if (!g->stack_ptr) 
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

        case '/':
          return error("Error: slash seen outside of ring - multipliers currently unsupported");
        
        case '\n':
          break; 

        default:
          return error("Error: invalid character read for WLN notation"); 
      }
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
      case CAR: 
        atom = ob_add_atom(mol, node->atomic_num, node->charge, 4 - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
        break; 

      case NIT:
        atom = ob_add_atom(mol, node->atomic_num, node->charge, 3 - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
        break; 
      
      case OXY:
        atom = ob_add_atom(mol, node->atomic_num, node->charge, 2 - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
        break; 

      case FLU:
      case CHL:
      case BRO:
      case IOD:
        atom = ob_add_atom(mol, node->atomic_num, node->charge, 1 - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
        break; 

      
      case DUM: // used to simplify grow code
        break; 

      default:
        atom = ob_add_atom(mol, node->atomic_num, node->charge, ((node->valence_pack & 0xF0) >> 4) - (node->valence_pack & 0x0F) ); 
        amapping[i] = atom; 
    }
  }
  
  for (u16 i=1; i<g->s_num; i++) {
    symbol_t *node = &g->symbols[i];
    if (node->atomic_num != DUM) {
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


    // WLN has no inherent stereochemistry, this can be a flag but should be off by default
  mol->SetChiralityPerceived(true);
  mol->SetAromaticPerceived(false);

  mol->DeleteHydrogens();

  return 1; 
}


// openbabel format reader function
int C_ReadWLN(const char *ptr, OpenBabel::OBMol* mol)
{   
  int ret = 0; 
  u16 st_pool_size = 128; 
  graph_t wln_graph; 
  graph_t *g = &wln_graph; 
  const u16 len = strlen(ptr); 
  
  gt_alloc(g, st_pool_size); 
  g->idx_symbols = (symbol_t**)malloc(sizeof(symbol_t*) * len+1); 
  memset(g->idx_symbols, 0, sizeof(symbol_t*) * len+1); 

  // allows recovery and resetting for large molecules (ideally never invoked)  
  ret = parse_wln(ptr, len, g); 
  while (ret == ERR_MEMORY && st_pool_size < 1024) {
    st_pool_size *= 2; 
    gt_free(g); 
    gt_alloc(g, st_pool_size); 
    ret = parse_wln(ptr, len, g); 
  }
  
  if (ret == ERR_NONE)
    ob_convert_wln_graph(mol, g);
  else 
    fprintf(stdout, "null\n");  

  gt_free(g); 
  free(g->idx_symbols); 
  return 1; 
}


int C_ReadWLNFile(FILE *fp, OpenBabel::OBMol* mol, OpenBabel::OBConversion *conv)
{   
  int ret = 0; 
  u16 st_pool_size = 128; 
  graph_t wln_graph; 
  graph_t *g = &wln_graph; 
  gt_alloc(g, st_pool_size); 
  
  u16 len; 
  char buffer[1024];
  memset(buffer,0,1024); // safety for len read  
  while (fgets(buffer, 1024, fp) != NULL) {
    // TODO: optimisation here, reuse the buffer
    len = strlen(buffer); // ignore nl
    g->idx_symbols = (symbol_t**)malloc(sizeof(symbol_t*) * len+1); 
    memset(g->idx_symbols, 0, sizeof(symbol_t*) * len+1); 
    
    ret = parse_wln(buffer, len, g); 
    while (ret == ERR_MEMORY && st_pool_size < 1024) {
      st_pool_size *= 2; 
      gt_free(g); 
      gt_alloc(g, st_pool_size); 
      ret = parse_wln(buffer, len, g); 
    }
    
    if (ret == ERR_NONE) {
      ob_convert_wln_graph(mol, g);
      conv->Write(mol, &std::cout); 

      gt_clear(g); 
      mol->Clear(); 
    }
    else 
      fprintf(stdout, "null\n");  

    free(g->idx_symbols); 
    memset(buffer,0,1024); // safety for len read  
  }
  
  gt_free(g);
  return 1; 
}

