/*********************************************************************
 
file: readwln.c 
author : Michael Blakey 2022, updated 2024 
description: WLN reader - write out SMILES etc from WLN

*********************************************************************/

#define OPENBABEL 1

#if OPENBABEL
/*********************************************************************
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
#endif

/* TODO: 
 *
 * - build into RDKit
 * - configure script for RDKit C++ build 
 * 
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if OPENBABEL
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/kekulize.h>

#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>
#endif 

#include "parser.h"

//#define DEBUG 1 // debug log - lvls: 0 - none, 1 - minimal, 2 - all

#define MAX_DEGREE 8
#define SYMBOL_MAX 256

// common bit fields
#define SPACE_READ  0x01 
#define DIGIT_READ  0x02 
#define DASH_READ   0x04

// standard parse fields
#define RING_READ   0x08
#define BIND_READ   0x10 // locant ring access
#define CHARGE_READ 0x20
#define DIOXO_READ  0x40 // this is an annoying one

// ring parse fields
#define SSSR_READ   0x08
#define MULTI_READ  0x10
#define PSEUDO_READ 0x20
#define AROM_READ   0x40

// Ring type "magic numbers"
#define LOCANT_BRIDGE 0x10  

// edge order as data packing 
#define HANGING_BOND 0x80 

// Error codes 
#define ERR_ABORT  0
#define ERR_NONE   1 // success

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

// 9 bytes
typedef struct {
  symbol_t* c; 
  u8 order; 
  u8 ring; 
} edge_t; 

// 4 + (8*9) = 76 bytes per symbol
// 80 = 16*5, multiple of 16 stack aligned
struct symbol_t {
  u8    atom_num;   
  char  arom; // -1 indicates ignore 
  char  charge;   
  u8    valence_pack;  // [  max u4    ][ curr     u4 ]  [0-8][0-8] 
  u8    n_bonds; 
  edge_t bonds[MAX_DEGREE]; // directional 
};

typedef struct locant {
  u8 hloc; // used for pathsolver 
  u8 r_pack; // [of 1b][ offL 1b ][ offR 1b ][ bridging 1b ][ dangling u4 ] 
  symbol_t *s; 
  locant   *off_path[2]; // 0 = -. 1 = '&' 
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
    char ref; // -1 for (ring_t*) else (symbol_t*) 
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
    if (g->stack[i].ref < 0)  
      free(g->stack[i].addr);
    g->stack[i].addr = 0; 
    g->stack[i].ref = 0;
  }
  g->stack_ptr = 0; 
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
    fprintf(stderr,"Warning: symbol limit reached (%d), this can be overidden in settings\n", g->s_max); 
    return (symbol_t*)0; 
  }
  else {
    symbol_t *s = 0; 
    
    if (!e->c)
      s = &g->symbols[g->s_num++];
    else 
      s = e->c; 

    s->atom_num     = id; 
    s->arom         = 0; 
    s->charge       = 0; 
    s->n_bonds      = 0;
    s->valence_pack = lim_valence; 
    s->valence_pack <<= 4;
    return s; 
  }
}

static symbol_t* new_symbol(graph_t *g, const u16 id, const u8 lim_valence)
{
  if (g->s_num == g->s_max) {
    fprintf(stderr,"Warning: symbol limit reached, reallocating...\n"); 
    return (symbol_t*)0; 
  }
  else {
    
    symbol_t *s = &g->symbols[g->s_num++];

    s->atom_num     = id; 
    s->arom         = 0; 
    s->charge       = 0; 
    s->n_bonds      = 0;
    s->valence_pack = lim_valence; 
    s->valence_pack <<= 4;
    return s; 
  }
}

/* used in the ring parse, leaves bonds alone */
static symbol_t* transmute_symbol(symbol_t *s, const u16 id, const u8 lim_valence)
{
  s->atom_num     = id; 
  s->charge       = 0; 
  s->valence_pack &= 0x0F; 
  s->valence_pack += lim_valence << 4; 
  return s; 
}

/* read off dash buffer, nearly all dash symbols 
 * are stack compatible atoms, 
 * transition metals are given a 6 valence limit - not allowed to expand
 * the octet. WLN does not have any rules for this, so this heuristic is
 * subject to interpretation 
*/

static __always_inline edge_t* next_virtual_edge(symbol_t *p)
{
  edge_t *e = &p->bonds[p->n_bonds++]; 
  e->order = (e->order | 0x01) & 0x03;  
  return e; 
}

static edge_t* check_bonding(edge_t *e, symbol_t *p, symbol_t *c)
{
  // TODO - nibble bit trick is definitely possible
  if ((p->valence_pack & 0x0F) > (p->valence_pack >> 4) || 
      (c->valence_pack & 0x0F) > (c->valence_pack >> 4)) 
  {
    fprintf(stderr,"Error: symbol reached WLN allowed valence - %d/%d & %d/%d\n",
            p->valence_pack & 0x0F, p->valence_pack >> 4,
            c->valence_pack & 0x0F, c->valence_pack >> 4); 
    return 0; 
  }
  else 
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
  return check_bonding(e, p, c); 
}


static void gt_load_stack_frame(symbol_t **p, edge_t **e, ring_t **r, graph_t *g)
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

static ring_t* rt_alloc(graph_t *g, const size_t size, symbol_t *inc, const u16 inc_pos) 
{
  ring_t *ring = 0; 
  ring = (ring_t*)malloc(sizeof(ring_t) + sizeof(locant)*(size-1)); 
  memset(ring, 0, sizeof(ring_t) + sizeof(locant)*(size-1)); 
  
  symbol_t *c = 0;  
  symbol_t *p = 0;
  edge_t *e   = 0; 

  for (u16 i=0; i<size; i++) {
    if (inc && i==inc_pos)
      c = inc; 
    else 
      c = new_symbol(g, CAR, 4); 
    
    if(!c)
      return (ring_t*)0; 
    else if (p) {
      e = next_virtual_edge(p); 
      e = set_virtual_edge(e, p, c); 
      e->ring = 1;  
    }

    ring->path[i].s = c; 
    ring->path[i].r_pack++; 
    ring->path[i].hloc = i+1; // point to the next locant 
    p = c; 
  }
  
  // end chain movement for pathfinderIII
  ring->path[0].r_pack++; 
  ring->path[size-1].r_pack++; 
  ring->path[size-1].hloc = size-1; 
  
  ring->size = size; 
  return ring; 
}


typedef struct  {
  u8 r_loc; // use to calculate size + add in off-branch positions
  u8 r_size;
  u8 arom; 
} r_assignment; 


static ring_t* pathsolverIII_fast(graph_t *g, ring_t *r, 
                                   r_assignment *SSSR, u8 SSSR_ptr) 
{
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
  */

  u8 steps, s_pos;
  edge_t *e; 
  locant *start, *end; 
  r_assignment *subcycle; 
  
  for (u16 i=0; i<SSSR_ptr; i++) {
    subcycle = &SSSR[i];   
    steps    = subcycle->r_size; 
    s_pos    = subcycle->r_loc; 
    start    = &r->path[s_pos]; 
    end = start; 
    
    // if used max times in ring, shift along path
    while ((start->r_pack & 0x0F) == 0 && s_pos < r->size) {
      start = &r->path[++s_pos];  
      steps--; 
    }
    
    for (u16 s=0; s<steps-1; s++) {
      end->s->arom |= subcycle->arom; 
      end = &r->path[end->hloc]; 
    }
    end->s->arom |= subcycle->arom; 

#if DEBUG 
    fprintf(stderr,"%d: %c --> %c (%d)\n",steps,start - &r->path[0] + 'A',end - &r->path[0] + 'A', subcycle->arom); 
#endif
    
    start->r_pack--; 
    end->r_pack--; 

    e = next_virtual_edge(start->s); 
    e = set_virtual_edge(e, start->s, end->s);
    if (!e)
      return 0; 
    else 
      e->ring = 1; 
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

static u8 add_oxy(graph_t *g, symbol_t *p, u8 ion)
{
  edge_t *e = next_virtual_edge(p); 
  symbol_t *c = next_symbol(g, e, OXY, 2); 
  
  if (!c)
    return ERR_ABORT; 
  else {
    if (!ion) {
      e->order++; 
      p->valence_pack++; 
    }
    else{
      p->valence_pack--; // let this be a free addition
      p->charge++; // WLN defines an immediate balance
      c->charge--; 
      c->valence_pack++; 
    }
    e = set_virtual_edge(e, p, c); 
  }
  if (!e)
    return ERR_ABORT; 
  else 
    return ERR_NONE; 
}

/* placeholder for charged alterations */ 
static u8 add_tauto_dioxy(graph_t *g, symbol_t *p)
{
  u8 ret = 0; 
  ret = add_oxy(g, p, 0); 
  if (ret != ERR_NONE)
    return ret; 
  else if (((p->valence_pack >> 4) - (p->valence_pack & 0x0F)) >= 2) {
    return add_oxy(g, p, 0); 
  }
  else 
    return add_oxy(g, p, 1); 
};

static u8 default_methyls(graph_t *g, symbol_t *c, const u8 n)
{
  edge_t *e; 
  symbol_t *m;
  u8 b_store = c->n_bonds; 
  for (u8 i=(c->valence_pack & 0x0F); i<n; i++) {
    e = next_virtual_edge(c); 
    m = next_symbol(g, e, CAR, 4); 
    if (!m)
      return ERR_ABORT; 
    else 
      e = set_virtual_edge(e, c, m); 
  }

  c->n_bonds = b_store;  
  return ERR_NONE; 
}

/* locants can be expanded or read as branches with '&' and '-'
 * chars. In order to move the state as the string is being read, 
 * these have to be handled per locant read. The logic this saves
 * outweighs the branches in the loop */
static u8 read_locant(const char *ptr, u16 *idx, u16 limit_idx)
{
  u8 locant = ptr[*idx] - 'A'; 
  for (u16 i=*idx+1; i<limit_idx; i++) {
    switch (ptr[i]) {
      case '&':
        locant += 23; 
        break; 

      case '-':
        fprintf(stderr,"off-branch reads need handling\n"); 
      
      default:
        return locant; 
    }
    (*idx)++; 
  }
  
  return locant; 
}


static ring_t* parse_cyclic(const char *ptr, const u16 start, u16 end, 
                            symbol_t *head, u16 head_loc, 
                            graph_t *g) 
{
  symbol_t  *c    = 0; 
  edge_t    *e    = 0; 
  ring_t    *ring = 0; 

  u8   locant_ch   = 0; 
  u8   arom_count = 0;
  u8   ring_size  = 0;  
  
  u8 SSSR_ptr  = 0; 
  r_assignment SSSR[32]; // this really is sensible for WLN
  
  u8 bridge_ptr  = 0; 
  u16 bridge_locants[16]; 

  u8 buff_ptr = 0; 
  unsigned char buffer[3]; 
    
  u8 *connection_table = 0; // stack allocates on pseudo read

  u8 state = 0; // bit field:                      
                 // [][][][pseudo][SSSR][dash][digit][space]

  unsigned char ch; 
  unsigned char ch_nxt; 

  state = SSSR_READ; 
  for (u16 sp=start; sp<end; sp++){
    ch = ptr[sp]; 
    ch_nxt = ptr[sp+1]; 

    if (state & PSEUDO_READ && (ch >= 'A' && ch <= 'Z')) {
      if (buff_ptr == 2) {
        fprintf(stderr,"Error: more than two pseudo bonds specified\n"); 
        goto ring_fail; 
      }
      else 
        buffer[buff_ptr++] = ch; 
    }
    else if (state & SPACE_READ && (ch >= 'A' && ch <= 'Z')) {
      locant_ch = read_locant(ptr, &sp, end);
      ch_nxt = ptr[sp+1]; // might get updated  
      
      // handle single letter bridges, requires check for ring build
      // bridges must come after multicyclic definition for correct syntax
      if (ch_nxt == 'J' || ch_nxt == 'T' || ch_nxt == ' ') {
        if (locant_ch > ring_size) {
          fprintf(stderr, "Error: out of bounds locant access"); 
          goto ring_fail; 
        }
        else 
          bridge_locants[bridge_ptr++] = locant_ch; 

        locant_ch = 0; 
      }
      state &= ~SPACE_READ; 
    }
    else if (state & MULTI_READ) {
      // must be the incoming ring size
      ring_size = read_locant(ptr, &sp, end) + 1; 
      ring = rt_alloc(g, ring_size, head, head_loc); 
      state &= ~MULTI_READ; 
    }
    else {
      // create the ring as necessary
      if (state & SSSR_READ && (ch >= 'A' && ch <= 'Z')){
        // end the SSSR read, start element reading
        if (ring || !SSSR_ptr) {
          fprintf(stderr, "Error: ring backbone failure - syntax ordering issue\n"); 
          goto ring_fail; 
        }
        else
          ring = rt_alloc(g, ring_size-bridge_ptr, head, head_loc); 

        state &= ~SSSR_READ; 
      } 

      switch (ch) {
        // special pi bonding
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
          if (state & SPACE_READ && state & SSSR_READ) {
            // multicyclic block start
            // move forward space to skip redundant block - should appear on space, 
            // if not, error. This is currently wrong

            for (u8 i=0; i < ch-'0'+1; i++) 
              read_locant(ptr, &++sp, end); 

            if (sp >= end || ptr[sp] != ' ') {
              fprintf(stderr, "Error: invalid format for multicyclic ring - %c\n", ptr[sp]);
              goto ring_fail; 
            }
            else {
              state |= MULTI_READ; 
              state &= ~SPACE_READ; 
              state &= ~SSSR_READ; 
            }
          }
          else if (state & SSSR_READ){
            ring_size += ch - '0' - (2*(ring_size>0)); 
            SSSR[SSSR_ptr].r_size   = ch - '0'; 
            SSSR[SSSR_ptr].arom     = 1; // default is aromatic 
            SSSR[SSSR_ptr++].r_loc  =  locant_ch; //read_locant(ptr, &sp, end); 
            locant_ch = 0; 
          }
          else {
            fprintf(stderr, "Error: digit used outside of known state\n"); 
            goto ring_fail; 
          }
          break; 
        
        case 'A':
        case 'C':
        case 'D':
        case 'J':
        case 'L':
        case 'R':
          fprintf(stderr, "Error: non-elemental symbols outside valid read states\n"); 
          goto ring_fail; 
          break; 
      
        case 'B':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {
            c = transmute_symbol(ring->path[locant_ch].s, BOR, 3); 
            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        case 'H':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {
            c = ring->path[locant_ch].s; 
            c->arom = -1; // will get ignored, bits overlap
            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        case 'K': 
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else { 
            c = transmute_symbol(ring->path[locant_ch].s, NIT, 4); 
            ring->path[locant_ch].r_pack++; 
            c->charge++; 
            
            c->bonds[0].order = HANGING_BOND; 
            c->bonds[1].order = HANGING_BOND;  
            c->bonds[2].order = HANGING_BOND; 

            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        case 'M':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else { 
            c = transmute_symbol(ring->path[locant_ch].s, NIT, 2); 
            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        case 'N':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {
            c = transmute_symbol(ring->path[locant_ch].s, NIT, 3); 
            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        case 'O':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {
            c = transmute_symbol(ring->path[locant_ch].s, OXY, 2); 
            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        case 'P':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else { 
            c = transmute_symbol(ring->path[locant_ch].s, PHO, 6); 
            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        case 'S':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {     
            c = transmute_symbol(ring->path[locant_ch].s, SUL, 6); 
            e = &c->bonds[0]; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 
  
        case 'T':
          state |= AROM_READ; 
          SSSR[arom_count++].arom = 0; 
          break; 

        // some nice ghost logic here. On a successful chain creation, the first bond will always to 
        // to the next symbol in the locant path. e.g A->B->C. therefore unsaturate bond[0]. 
        // when mixed with dash, this will be the next avaliable bond which is not in the chain, and MUST
        // be handled in the pathsolver algorithm, there place at bonds[1...n]. 
        case 'U':
          if (ch_nxt == '-') {
            

          }
          else if (e) {
            // means c must be live
            c->valence_pack++; 
            e->order++; 
            e->c->valence_pack++; 
          }
          else if (locant_ch == ring_size - 1) {
            // will wrap back onto the start, again, determined position likely at bonds[1], when logic shifts last bond, unsaturation on last position is technically undefined. 
            c = ring->path[0].s;  
            c->bonds[1].order += 2; // the edge is not yet live ~ 2
            c->valence_pack++; 
          }
          else if (locant_ch < ring_size) {
            c = ring->path[locant_ch].s;  
            c->bonds[0].order++; 
            c->valence_pack++; 
            ring->path[locant_ch+1].s->valence_pack++; 
            locant_ch++; 
          }
          break; 
          
        case 'V':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {
            c = transmute_symbol(ring->path[locant_ch].s, CAR, 4);
            e = &c->bonds[0]; 
            if (!add_oxy(g, c, 0))
              goto ring_fail; 
            else { 
              g->idx_symbols[sp+1] = c; 
              locant_ch++; 
            }
          }
          break; 

        case 'W':
          break; 

        case 'X':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {
            c = ring->path[locant_ch].s; 

            // purely virtual bonds that may or may not get tidyed up
            c->bonds[1].order = HANGING_BOND; 
            c->bonds[2].order = HANGING_BOND; 
            c->bonds[3].order = HANGING_BOND; 

            e = &c->bonds[0]; 
            ring->path[locant_ch].r_pack++; // the me
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

      
        case 'Y':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant access\n"); 
            goto ring_fail; 
          }
          else {
            c = ring->path[locant_ch].s; 
            c->valence_pack++; 
            e = &c->bonds[1]; // not involed in ring, if used more than once in ring graph, undefined behaviour
            e->order = 2; 
            g->idx_symbols[sp+1] = c; 
            locant_ch++; 
          }
          break; 

        
        // undefined terminator ring sequences
        case 'E':
        case 'F':
        case 'G':
        case 'I':
        case 'Q':
        case 'Z':
          fprintf(stderr, "Error: terminators as ring atoms are undefined symbols\n"); 
          goto ring_fail; 
          break; 

        case '&':
          if (state & MULTI_READ)
            ring->size += 23; 
          else if (state & PSEUDO_READ)
            buffer[buff_ptr-1] += 23; 
          else if (state & AROM_READ) 
            SSSR[arom_count++].arom = 1; 
          else if (locant_ch && locant_ch + 23 < ring_size)
            locant_ch += 23; 
          else {
            SSSR[arom_count++].arom = 1; 
            state |= AROM_READ; 
          }
          break; 

        case '/':
          if (state & SSSR_READ) {
            state &= ~(SSSR_READ); 
            buff_ptr = 0; 
            state |= PSEUDO_READ; 
          }
          else if (state & PSEUDO_READ) {
            if (buff_ptr != 2) {
              fprintf(stderr, "Error: pseudo locants must come in pairs\n"); 
              goto ring_fail; 
            }
            else 
              fprintf(stderr,"need impl\n"); 
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
              fprintf(stderr, "Error: pseudo locants must come in pairs\n"); 
              goto ring_fail; 
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
          e = 0; 
          break;

        default:
          fprintf(stderr,"Error: unknown symbol %c in ring parse\n",ch); 
          goto ring_fail; 
      }
    }
  }
  
  if (state & SSSR_READ) {
    // this must be a poly cyclic ring, 
    ring = rt_alloc(g, ring_size - bridge_ptr, head, head_loc); 
  }

  // resolve any bridges 
  for (u16 i=0; i<bridge_ptr; i++) {
    if (bridge_locants[i] >= ring->size) {
      fprintf(stderr,"Error: out of bounds locant access\n"); 
      goto ring_fail; 
    } 
    else 
      ring->path[bridge_locants[i]].r_pack--; 
  }


  // forward any single arom assignments
  for (u16 i=arom_count;i<SSSR_ptr;i++)
    SSSR[i].arom = SSSR[0].arom; 
  

  return pathsolverIII_fast(g, ring, SSSR, SSSR_ptr);     

// might be evil, tidies up the errors handling though
ring_fail:
  if (ring)
    free(ring);

  return (ring_t*)0; 
}

static u8 write_dash_symbol(symbol_t *c, u8 high, u8 low){
  switch (high){
    case 'A':
      switch (low) {
        case 'C':
          c->atom_num = 89; 
          break; 
        case 'G':
          c->atom_num = 47; 
          break; 
        case 'L':
          c->atom_num = 13; 
          break; 
        case 'M':
          c->atom_num = 95; 
          break; 
        case 'R':
          c->atom_num = 18; 
          break; 
        case 'S':
          c->atom_num = 33; 
          break; 
        case 'T':
          c->valence_pack = (1 << 4); 
          c->atom_num = 85; 
          break; 
        case 'U':
          c->atom_num = 79; 
          break; 
      }
      break; 

    case 'B':
      switch (low) {
        case 0:
          c->valence_pack = (4 << 4);
          c->atom_num = BOR; 
          break; 
        case 'A':
          c->valence_pack = (2 << 4);
          c->atom_num = 56; 
          break; 
        case 'E':
          c->valence_pack = (2 << 4);
          c->atom_num = 4; 
          break; 
        case 'H':
          c->atom_num = 107; 
          break; 
        case 'I':
          c->atom_num = 83; 
          break; 
        case 'K':
          c->atom_num = 97; 
          break; 
        case 'R':
          c->valence_pack = (1 << 4); 
          c->atom_num = BRO; 
          break; 
      }
      break; 

    case 'C':
      switch (low) {
        case 0:
          c->valence_pack = (5 << 4); 
          c->atom_num = CAR; // allow the texas carbon "easter egg"
          break; 
        case 'A':
          c->valence_pack = (2 << 4); 
          c->atom_num = 20;
          break; 
        case 'D':
          c->atom_num = 48;
          break; 
        case 'E':
          c->atom_num = 58;
          break; 
        case 'F':
          c->atom_num = 98;
          break; 
        case 'M':
          c->atom_num = 96;
          break; 
        case 'N':
          c->atom_num = 112;
          break; 
        case 'O':
          c->atom_num = 27;
          break; 
        case 'R':
          c->atom_num = 24;
          break; 
        case 'S':
          c->valence_pack = (1 << 4); 
          c->atom_num = 55;
          break; 
        case 'U':
          c->atom_num = 29;
          break; 
      }
      break; 

    case 'D':
      switch (low) {
        case 'B':
          c->atom_num = 105;
          break; 
        case 'S':
          c->atom_num = 110;
          break; 
        case 'Y':
          c->atom_num = 66;
          break; 
      }
      break;

    case 'E':
      switch (low) {
        case 0:
          c->valence_pack = (3 << 4); // hypervalent bromine
          c->atom_num = 35; 
          break; 
        case 'R':
          c->atom_num = 68; 
          break; 
        case 'S':
          c->atom_num = 99; 
          break; 
        case 'U':
          c->atom_num = 63; 
          break; 
      }

    case 'F':
      switch (low) {
        case 0:
          c->valence_pack = (3 << 4); // hypervalent flourine
          c->atom_num = FLU; 
          break; 
        case 'E':
          c->atom_num = 26; 
          break; 
        case 'L':
          c->atom_num = 114; 
          break; 
        case 'M':
          c->atom_num = 100; 
          break; 
        case 'R':
          c->atom_num = 87; 
          break; 
      }
      break;

    case 'G':
      switch (low) {
        case 0: 
          c->valence_pack = (3 << 4); // hypervalent chlorine
          c->atom_num = CHL; 
          break; 
        case 'A':
          c->valence_pack = (3 << 4); 
          c->atom_num = 31; 
          break; 
        case 'D':
          c->atom_num = 64; 
          break; 
        case 'E':
          c->valence_pack = (4 << 4); 
          c->atom_num = 32; 
          break; 
      }
      break;

    case 'H':
      switch (low) {
        case 'E':
          c->valence_pack = (1 << 4); 
          c->atom_num = 2; 
          break; 
        case 'F':
          c->atom_num = 72; 
          break; 
        case 'G':
          c->atom_num = 80; 
          break; 
        case 'O':
          c->atom_num = 67; 
          break; 
        case 'S':
          c->atom_num = 108; 
          break; 
      }
      break;

    case 'I':
      switch (low) {
        case 0: 
          c->valence_pack = (2 << 4); // hypervalent iodine
          c->atom_num = IOD;
          break; 
        case 'N':
          c->valence_pack = (3 << 4); 
          c->atom_num = 49; 
          break; 
        case 'R':
          c->atom_num = 77; 
          break; 
      }
      break;

    case 'K':
      switch (low) {
        case 0:  
          c->valence_pack = (4 << 4); // already hyper nitrogen 
          c->atom_num = NIT;
          c->charge++; 
          break; 
        case 'R':
          c->valence_pack = (1 << 4); 
          c->atom_num = 36; 
          break; 
        case 'A':
          c->valence_pack = (1 << 4); 
          c->atom_num = 19; 
          break; 
      }
      break;

    case 'L':
      switch (low) {
        case 'A':
          c->atom_num = 57; 
          break; 
        case 'I':
          c->valence_pack = (1 << 4); 
          c->atom_num = 3; 
          break; 
        case 'R':
          c->atom_num = 103; 
          break; 
        case 'U':
          c->atom_num = 71; 
          break; 
        case 'V':
          c->atom_num = 116; 
          break; 
      }
      break;

    case 'M':
      switch (low) {
        case 0:
          c->valence_pack = (2 << 4); // regular nitrogen
          c->atom_num = NIT; 
          break; 
        case 'C':
          c->atom_num = 115; 
          break; 
        case 'D':
          c->atom_num = 101; 
          break; 
        case 'G':
          c->valence_pack = (2 << 4); 
          c->atom_num = 12; 
          break; 
        case 'N':
          c->atom_num = 25; 
          break; 
        case 'O':
          c->atom_num = 42; 
          break; 
        case 'T':
          c->atom_num = 109; 
          break; 
      }
      break;

    case 'N':
      switch (low) {
        case 0:
          c->valence_pack = (2 << 4); // regular nitrogen
          c->atom_num = NIT; 
          break; 
        case 'A':
          c->valence_pack = (1 << 4); 
          c->atom_num = 11; 
          break; 
        case 'B':
          c->atom_num = 41; 
          break; 
        case 'D':
          c->atom_num = 60; 
          break; 
        case 'E':
          c->valence_pack = (1 << 4); 
          c->atom_num = 10; 
          break; 
        case 'H':
          c->atom_num = 113; 
          break; 
        case 'I':
          c->atom_num = 28; 
          break; 
        case 'O':
          c->atom_num = 102; 
          break; 
        case 'P':
          c->atom_num = 93; 
          break; 
      }
      break; 

    case 'O':
      switch (low) {
        case 0:
          c->valence_pack = (3 << 4); // hypervalent oxygen
          c->atom_num = OXY; 
          break; 
        case 'G':
          c->atom_num = 118; 
          break; 
        case 'S':
          c->atom_num = 76; 
          break; 
      }
      break;

    case 'P':
      switch (low) {
        case 0:
          c->valence_pack = (5 << 4); // phos
          c->atom_num = PHO; 
          break; 
        case 'A':
          c->atom_num = 91; 
          break; 
        case 'B':
          c->atom_num = 82; 
          break; 
        case 'D':
          c->atom_num = 46; 
          break; 
        case 'M':
          c->atom_num = 61; 
          break; 
        case 'O':
          c->atom_num = 84; 
          break; 
        case 'R':
          c->atom_num = 59; 
          break; 
        case 'T':
          c->atom_num = 78; 
          break; 
        case 'U':
          c->atom_num = 94; 
          break; 
      }
      break;
    
    case 'Q':
      c->valence_pack = (3 << 4); // hypervalent oxygen
      c->atom_num = OXY; 
      break; 

    case 'R':
      switch (low) {
        case 'A':
          c->atom_num = 88; 
          break; 
        case 'B':
          c->valence_pack = (1 << 4); 
          c->atom_num = 37; 
          break; 
        case 'E':
          c->atom_num = 75; 
          break; 
        case 'F':
          c->atom_num = 104; 
          break; 
        case 'G':
          c->atom_num = 111; 
          break; 
        case 'H':
          c->atom_num = 45; 
          break; 
        case 'N':
          c->atom_num = 86; 
          break; 
        case 'U':
          c->atom_num = 44; 
          break; 
      }
      break;

    case 'S':
      switch (low) {
        case 0: 
          c->valence_pack = (6 << 4); // regular sulphur
          c->atom_num = SUL; 
          break; 
        case 'B':
          c->valence_pack = (3 << 4); 
          c->atom_num = 51; 
          break; 
        case 'C':
          c->valence_pack = (3 << 4); 
          c->atom_num = 21; 
          break; 
        case 'E':
          c->valence_pack = (2 << 4); 
          c->atom_num = 34; 
          break; 
        case 'G':
          c->atom_num = 106; 
          break; 
        case 'I':
          c->valence_pack = (4 << 4); 
          c->atom_num = 14; 
          break; 
        case 'M':
          c->atom_num = 62; 
          break; 
        case 'N':
          c->valence_pack = (4 << 4); 
          c->atom_num = 50; 
          break; 
        case 'R':
          c->valence_pack = (2 << 4); 
          c->atom_num = 38; 
          break; 
      }
      break;

    case 'T':
      switch (low) {
        case 'A':
          c->atom_num = 73; 
          break; 
        case 'B':
          c->atom_num = 65; 
          break; 
        case 'C':
          c->atom_num = 43; 
          break; 
        case 'E':
          c->valence_pack = (2 << 4); 
          c->atom_num = 52; 
          break; 
        case 'H':
          c->atom_num = 90; 
          break; 
        case 'I':
          c->atom_num = 22; 
          break; 
        case 'L':
          c->valence_pack = (3 << 4); 
          c->atom_num = 81; 
          break; 
        case 'M':
          c->atom_num = 69; 
          break; 
        case 'S':
          c->atom_num = 117; 
          break; 
      }
      break;

    case 'U':
      if(low == 'R') {
        c->atom_num = 92; 
        break; 
      }
      break;

    case 'V':
      if(low == 'A') {
        c->atom_num = 23; 
        break; 
      }
      break;
  
    case 'W':
      if(low == 'T') {
        c->atom_num = 74; 
        break; 
      }
      break; 

    case 'X':
      if (low == 'E') {
        c->valence_pack = (1 << 4); 
        c->atom_num = 54; 
        break; 
      }
      break;

    case 'Y':
      if (low == 'T') {
        c->valence_pack = (3 << 4); 
        c->atom_num = 39; 
        break; 
      }
      else if (low == 'B') {
        c->atom_num = 70; 
        break; 
      }
      break;

    case 'Z':
      if (low == 'N') {
        c->atom_num = 30; 
        break; 
      }
      else if (low == 'R') {
        c->atom_num = 30; 
        break; 
      }
      break;
  }
  
  if (c->atom_num == DUM) {
    fprintf(stderr,"Error: invalid elemental code -%c%c-\n",high,low); 
    return ERR_ABORT; 
  }
  else 
    return ERR_NONE; 
}




/* returns the pending unsaturate level for the next symbol */
static edge_t* resolve_unsaturated_carbon(graph_t *g, edge_t *e, symbol_t *p, symbol_t *c, unsigned char ch_nxt)
{
  
  edge_t *ne = next_virtual_edge(c); // forward edge from 'C' 

  // first check can you double bond to the previous - 
  // the bond will be set to 1 so its >= 1
  if ((p->valence_pack>>4) - (p->valence_pack & 0x0F) >= 1) {
    
    // next predict what bonding is possible to whats coming next
    switch (ch_nxt) {
      // alkyl chains have a weird relationship where they can but shouldn't unsaturate here - murky waters with implied chemical information
      case '1':  
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
      case 'Z':
      case 'Q':
      case 'E':
      case 'F':
      case 'G':
      case 'H':
      case 'I':
        // force a triple bond to previous 
        e->order += 2; 
        p->valence_pack += 2;
        c->valence_pack += 2; 
        return check_bonding(ne, p, c); 

      // if a triple can be made, force a triple forward
      case 'N':
      case 'K':
      case 'X':
      case 'Y':
      case 'B':
      case 'C': // very murky 
      case 'P':
      case 'S':
      case '-':
        ne->order += 2; 
        c->valence_pack += 2; 
        return ne; 
      
      // split the a double bond between the two
      default:
        e->order++; 
        ne->order++; 
        p->valence_pack++;
        c->valence_pack += 2;
        return check_bonding(ne, p, c); 
    } 
  }
  else {
    // must be a triple to the next bond
    ne->order += 2; 
    c->valence_pack += 2; 
    return ne; 
  }
  
  fprintf(stderr,"Error: C symbol requires a mandatory double/triple bond - impossible bonding env\n"); 
  return (edge_t*)0; 
}

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

  u8 ring_chars = 0; 
  u16 locant_ch = 0; 
  u16 digit_n   = 0; 
  
  u8 state = 0; // bit field: 
                // [][dioxo][charge][ring locant][ring skip][dash][digit][space]
  
  u8 dash_ptr = 0; 
  unsigned char dash_chars[3] = {0}; // last byte is mainly for overflow 
 
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
      if (ptr[sp-1] >= 'A' && ch == 'J' && ch_nxt < '0') {
        // J can be used inside ring notation, requires lookahead 
        // condition 
        
        // note: The ptr passed in does not include the 
        //       starting L/T or ending J (<sp) 
        
        // ring must have at least one atom. 

        c = next_symbol(g, e, CAR, 4);
        if (!c)
          return ERR_ABORT; 
        else
          e = set_virtual_edge(e, p, c); 

        r = parse_cyclic(ptr, sp-ring_chars+1, sp, c, locant_ch, g);
        if (!r)
          return ERR_ABORT; 
        else {
          g->stack[g->stack_ptr].addr = r; 
          g->stack[g->stack_ptr++].ref = -1; 
          state &= ~(RING_READ);
          ring_chars = 0; 
        }

        state &= ~(BIND_READ); 
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

      if (!(state & BIND_READ)) {
        if (r && locant_ch < r->size) {
          c = r->path[locant_ch].s; 
          e = next_virtual_edge(c); 
          p = c; 
        } 
        else 
          return error("Error: out of bounds locant access"); 
      }
    }
    else {
      switch (ch) {

        case '0':
          digit_n *= 10; 

          if(state & DIOXO_READ) 
            return error("Error: dioxo attachment needs higher valence atom"); 
          else if (ch_nxt == '/') {
            // positive charge assignment
            if (digit_n > len || !g->idx_symbols[digit_n]) 
              return error("Error: charge assignment out of bounds"); 
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
              if (digit_n > len || !g->idx_symbols[digit_n]) 
                return error("Error: charge assignment out of bounds"); 
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
                  return ERR_ABORT; 
                else
                  e = set_virtual_edge(e, p, c); 

                if (!e)
                  return ERR_ABORT; 
                else {
                  p = c;
                  e = next_virtual_edge(p); 
                }
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

          if(state & DIOXO_READ) 
            return error("Error: dioxo attachment needs higher valence atom"); 
          else if (ch_nxt == '/') {
            // positive charge assignment
            if (digit_n > len || !g->idx_symbols[digit_n]) 
              return error("Error: charge assignment out of bounds"); 
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
              if (digit_n > len || !g->idx_symbols[digit_n])
                return error("Error: charge assignment out of bounds"); 
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
                  return ERR_ABORT; 
                else
                  e = set_virtual_edge(e, p, c); 

                if (!e)
                  return ERR_ABORT; 
                else {
                  p = c;
                  e = next_virtual_edge(p); 
                }
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
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {

            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, c); 
              state &= ~DIOXO_READ; 
            }
            else {
              g->stack[g->stack_ptr].addr = c; 
              g->stack[g->stack_ptr].ref  = 2;
              g->stack_ptr++; 
            }
            p = c;
            g->idx_symbols[sp+1] = p; 
            e = next_virtual_edge(c); 
          }
          break; 
        
        case 'C':
          c = next_symbol(g, e, CAR, 4);
          if (!c)
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = p; 
            e = resolve_unsaturated_carbon(g, e, p, c, ch_nxt); 
            if (!e)
              return ERR_ABORT; 
            else 
              p = c;
          }
          break; 

        // extrememly rare open chelate notation
        case 'D':
          fprintf(stderr,"chelate ring open needs handling\n"); 
          break; 

        case 'E':
          c = next_symbol(g, e, BRO, 1);
          if (!c)
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          if (!e)
            return ERR_ABORT; 
          else {
            // terminating symbol
            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                gt_load_stack_frame(&p, &e, &r, g); 
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
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          if (!e)
            return ERR_ABORT; 
          else {
            // terminating symbol
            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                gt_load_stack_frame(&p, &e, &r, g); 
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
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          if (!e)
            return ERR_ABORT; 
          else {
            // terminating symbol
            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                gt_load_stack_frame(&p, &e, &r, g); 
            }
            else {
              p = c; 
              e = next_virtual_edge(p); 
            }
          }
          break;

        case 'H':
          // treat hydrogen as a terminator
          c = next_symbol(g, e, 1, 1);
          if (!c)
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          if (!e)
            return ERR_ABORT; 
          else {

            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                gt_load_stack_frame(&p, &e, &r, g); 
            }
            else {
              // do not allow hydrogen to write back
              e = next_virtual_edge(p); 
            }
          }
          break;

        case 'I':
          c = next_symbol(g, e, IOD, 1);
          if (!c)
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          if (!e)
            return ERR_ABORT; 
          else {

            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                gt_load_stack_frame(&p, &e, &r, g); 
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
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {

            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, c); 
              state &= ~DIOXO_READ; 
            }
            else {
              default_methyls(g, c, 4);  

              g->stack[g->stack_ptr].addr = c; 
              g->stack[g->stack_ptr].ref  = 3;
              g->stack_ptr++; 
            }
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
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            
            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else {
              p = c;
              g->idx_symbols[sp+1] = p; 
              e = next_virtual_edge(c); 
            }

          }
          break;
       
        /* nitrogen symbols */
        case 'N':
          c = next_symbol(g, e, NIT, 3);
          if (!c)
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
  
            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, c); 
              state &= ~DIOXO_READ; 
            }
            else {
              g->stack[g->stack_ptr].addr = c; 
              g->stack[g->stack_ptr].ref  = 2;
              g->stack_ptr++; 
            }

            p = c;
            g->idx_symbols[sp+1] = p; 
            e = next_virtual_edge(c); 
          }
          break;
        
        case 'O':
          c = next_symbol(g, e, OXY, 2);
          if (!c)
            return ERR_ABORT; 
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
        
        case 'P':
          c = next_symbol(g, e, PHO, 5);
          if (!c)
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            // since these can expand valence, allow branch
            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, c); 
              state &= ~DIOXO_READ; 
            }

            g->stack[g->stack_ptr].addr = c; 
            g->stack[g->stack_ptr].ref  = 3;
            g->stack_ptr++; 

            p = c; 
            g->idx_symbols[sp+1] = c; 
            e = next_virtual_edge(c); 
          }
          break;


        case 'Q':
          c = next_symbol(g, e, OXY, 1);
          if (!c)
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          if (!e)
            return ERR_ABORT; 
          else {
            // terminating symbol
            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                gt_load_stack_frame(&p, &e, &r, g); 
            }
            else {
              p = c; 
              e = next_virtual_edge(p); 
            }
          }
          break;

        // shorthand benzene 
        case 'R': 
          // head symbol can come from virtual X|Y|K
          c = next_symbol(g, e, CAR, 4);
          if (!c)
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          r = rt_alloc(g, 6, c, 0); 
          
          if (!r)
            return ERR_ABORT; // only way this can fail 
          else {
            for (u8 i=0; i<6; i++)
              r->path[i].s->arom |= 1; 
            // set the ring. R objects are kekulised on stack pop
            p = r->path[0].s; 
            c = r->path[5].s; 
            e = next_virtual_edge(p); 
            e = set_virtual_edge(e, p, c);
            e->ring = 1; 
          }

          // add to ring stack
          g->stack[g->stack_ptr].addr = r; 
          g->stack[g->stack_ptr++].ref = -1; 

          g->idx_symbols[sp+1] = p; 
          e = next_virtual_edge(p); 
          break; 
        
        case 'S':
          c = next_symbol(g, e, SUL, 6);
          if (!c)
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, c); 
              state &= ~DIOXO_READ; 
            }
            // since these can expand valence, allow branch
            g->stack[g->stack_ptr].addr = c; 
            g->stack[g->stack_ptr].ref  = 3; // ?  
            g->stack_ptr++; 

            p = c; 
            g->idx_symbols[sp+1] = c; 
            e = next_virtual_edge(c); 
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
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {

            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (add_oxy(g, c, 0) == ERR_ABORT)
              return ERR_ABORT; 
            else {
              p = c; 
              g->idx_symbols[sp+1] = c; 
              e = next_virtual_edge(c); 
            }
          }
          break;
        
        // if not previous, create dummy carbon
        // really trying to avoid a bit state on this
        case 'W':
          if(state & DIOXO_READ) 
            return error("Error: double dioxo attachment needs is undefined"); 
          else if (p->atom_num != DUM) {
            if (add_tauto_dioxy(g, p) == ERR_ABORT)
              return ERR_ABORT; 
            else 
              e = next_virtual_edge(p);  
          }
          else 
            state |= DIOXO_READ; 
          break; 

        case 'X':
          c = next_symbol(g, e, CAR, 4);
          if (!c)
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            
            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, c); 
              state &= ~DIOXO_READ; 
            }
            else {
              default_methyls(g, c, 4);  

              g->stack[g->stack_ptr].addr = c; 
              g->stack[g->stack_ptr].ref  = 3;
              g->stack_ptr++; 
            }

            p = c; 
            g->idx_symbols[sp+1] = c; 
            e = next_virtual_edge(c); 
          }
          break;

        case 'Y':
          c = next_symbol(g, e, CAR, 4);
          if (!c)
            return ERR_ABORT; 
          else
            e = set_virtual_edge(e, p, c); 

          if (!e)
            return ERR_ABORT; 
          else {
            
            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, c); 
              state &= ~DIOXO_READ; 
            }
            else {
              default_methyls(g, c, 3);  

              g->stack[g->stack_ptr].addr = c; 
              g->stack[g->stack_ptr].ref  = 2;
              g->stack_ptr++; 
            }

            p = c; 
            g->idx_symbols[sp+1] = c; 
            e = next_virtual_edge(c); 
          }
          break;
        
        case 'Z':
          c = next_symbol(g, e, NIT, 1);
          if (!c)
            return ERR_ABORT; 
          else {
            g->idx_symbols[sp+1] = c; 
            e = set_virtual_edge(e, p, c); 
          }

          if (!e)
            return ERR_ABORT; 
          else {
            // terminating symbol
            if(state & DIOXO_READ) 
              return error("Error: dioxo attachment needs higher valence atom"); 
            else if (g->stack_ptr && g->stack[g->stack_ptr-1].ref != -1){
              g->stack[g->stack_ptr-1].ref--; 
              g->stack_ptr -= (g->stack[g->stack_ptr-1].ref==0); 

              // reseting block
              // note: terminators can empty the stack, '&' cannot
              if (!g->stack_ptr) {
                p = c; 
                e = next_virtual_edge(p); 
              }
              else 
                gt_load_stack_frame(&p, &e, &r, g); 
            }
            else {
              p = c; 
              e = next_virtual_edge(p); 
            }
          }
          break;
        
        case '-':
          if (state & DASH_READ) {
            c = next_symbol(g, e, DUM, 6); // create a blank symbol
            if (!c)
              return ERR_ABORT; 
            else {

              if (write_dash_symbol(c, dash_chars[0], dash_chars[1]) == ERR_ABORT)
                return ERR_ABORT; 
              else {
                e = set_virtual_edge(e, p, c); 

                if(state & DIOXO_READ) {
                  add_tauto_dioxy(g, c); 
                  state &= ~DIOXO_READ; 
                }
                
                memset(dash_chars,0,3); 
                state &= ~DASH_READ; 
                dash_ptr = 0; 

                // add to branch stack
                if ((c->valence_pack >> 4) >= 2) {
                  g->stack[g->stack_ptr].addr = c; 
                  g->stack[g->stack_ptr].ref  = c->valence_pack >> 4;
                  g->stack_ptr++; 
                }

                p = c; 
                g->idx_symbols[sp+1] = c; 
                e = next_virtual_edge(c); 
              }
            }
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
          if (ch_nxt == '&') {
            // all other states should make this the only place ions can be used. 
            gt_stack_flush(g); 
            sp++; 
            c = new_symbol(g, DUM, 1); 
            e = next_virtual_edge(c); 
            e->order = DUM; 
            p = c; 
          }
          else 
            state |= SPACE_READ; 
          break; 

        case '&':
          if(state & DIOXO_READ) 
            return error("Error: dioxo attachment requires an atomic symbol"); 
          else if (!g->stack_ptr) 
            return error("Error: empty stack - too many &?"); 
          else { 

            if (g->stack[g->stack_ptr-1].ref < 0) {
              // ring closures  
              g->stack_ptr--; 
              free(g->stack[g->stack_ptr].addr);
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
              gt_load_stack_frame(&p, &e, &r, g); 
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
  else     
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
  OpenBabel::OBBond *bond;   
  OpenBabel::OBAtom **amapping = (OpenBabel::OBAtom **)alloca(sizeof(OpenBabel::OBAtom*) * g->s_num); 
  for (u16 i=1; i<g->s_num; i++) {
    symbol_t *node = &g->symbols[i]; 
    if (node->atom_num != DUM) {
      u8 limit_val = node->valence_pack >> 4; 
      u8 actual_val = node->valence_pack & 0x0F; 
      
      if (node->arom == 1)
        actual_val += (limit_val != actual_val); 

      switch (node->atom_num) {
        case CAR: 
          atom = ob_add_atom(mol, node->atom_num, node->charge, 4 - actual_val); 
          break; 

        case NIT:
          atom = ob_add_atom(mol, node->atom_num, node->charge, 3 - actual_val); 
          break; 
        
        case OXY:
          node->charge |= node->arom; // ingvar input
          node->charge += -(!node->charge)*(limit_val - actual_val);
          node->valence_pack += abs(node->charge); 

          if (node->valence_pack >> 4 == 1)
            atom = ob_add_atom(mol, node->atom_num, node->charge, 1); 
          else 
            atom = ob_add_atom(mol, node->atom_num, node->charge, 0); 
          break; 

        case FLU:
        case CHL:
        case BRO:
        case IOD:
          // take a default -1 instead of H resolve
          node->charge += -(!node->charge)*(actual_val == 0); 
          atom = ob_add_atom(mol, node->atom_num, node->charge, 0); 
          break; 

        case SUL:
          atom = ob_add_atom(mol, node->atom_num, node->charge, (6 - actual_val) % 2); // sneaky H trick  
          break;

        case PHO: 
          atom = ob_add_atom(mol, node->atom_num, node->charge, (5 - actual_val) % 2); // sneaky H trick  
          break;

        default: // dont add hydrogens
          atom = ob_add_atom(mol, node->atom_num, node->charge, 0); 
      }
      
      amapping[i] = atom; 
      if (node->arom == 1)
        atom->SetAromatic(node->arom); 
    }
  }
  
  for (u16 i=1; i<g->s_num; i++) {
    symbol_t *node = &g->symbols[i];
    if (node->atom_num != DUM) {
      for (u16 j=0; j<MAX_DEGREE; j++) {
        u16 end; 
        u16 beg = i; 
        edge_t *e = &node->bonds[j];
        if (e->c) {
          end = e->c - g->symbols;
          bond = ob_add_bond(mol, amapping[beg], amapping[end], e->order);
          if (e->ring && node->arom == 1 && e->c->arom == 1)
            bond->SetAromatic(1); 
          else
            bond->SetAromatic(0);
        }
        // purely virtual additions which need tidying up
        else if (e->order == HANGING_BOND && (node->valence_pack & 0x0F) < (node->valence_pack>>4)) {
          atom = ob_add_atom(mol, CAR, 0, 3); 
          bond = ob_add_bond(mol, amapping[beg], atom, 1);
          amapping[beg]->SetImplicitHCount(amapping[beg]->GetImplicitHCount()-1); 
          node->valence_pack++;
        }
      }
    }
  }

    // WLN has no inherent stereochemistry, this can be a flag but should be off by default
  mol->SetChiralityPerceived(true);
  mol->SetAromaticPerceived(true); // let the toolkit do maximal matching
  mol->DeleteHydrogens();

  OpenBabel::OBKekulize(mol); 

  return 1; 
}

#if OPENBABEL

// openbabel format reader function
int C_ReadWLN(const char *ptr, OpenBabel::OBMol* mol)
{   
  int ret = 0; 
  graph_t wln_graph; 
  graph_t *g = &wln_graph; 
  const u16 len = strlen(ptr); 
  
  gt_alloc(g, SYMBOL_MAX); 
  g->idx_symbols = (symbol_t**)malloc(sizeof(symbol_t*) * len+1); 
  memset(g->idx_symbols, 0, sizeof(symbol_t*) * len+1); 

  ret = parse_wln(ptr, len, g); 
  if (ret == ERR_NONE)
    ob_convert_wln_graph(mol, g);
  else 
    fprintf(stdout, "null (ERR %d)\n", ret);  

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
  g->idx_symbols = (symbol_t**)malloc(sizeof(symbol_t*) * 128); 
  
  u16 len;
  u16 len_high = 128; 
  char buffer[1024];
  memset(buffer,0,1024); // safety for len read  
  while (fgets(buffer, 1024, fp) != NULL) {
    len = strlen(buffer); // ignore nl
    if (len > len_high) {
      g->idx_symbols = (symbol_t**)realloc(g->idx_symbols,sizeof(symbol_t*) * len_high); 
      len_high = len; 
    }
    memset(g->idx_symbols,0,sizeof(symbol_t*)*len_high); 
    
    ret = parse_wln(buffer, len, g); 
    if (ret == ERR_NONE) {
      ob_convert_wln_graph(mol, g);
      conv->Write(mol, &std::cout); 
      gt_clear(g); 
      mol->Clear(); 
    }
    else 
      fprintf(stdout, "null\n");  
     
    memset(buffer,0,1024); // safety for len read  
  }

  gt_free(g);
  free(g->idx_symbols); 
  return 1; 
}

#endif 
