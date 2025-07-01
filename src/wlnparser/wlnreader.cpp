/*********************************************************************
 
file: readwln.c 
author : Michael Blakey 2022, updated 2024 
description: WLN reader - write out SMILES etc from WLN

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
#include <stdbool.h>
#include <inttypes.h>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/kekulize.h>

#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>

#include "wlnparser.h"

//#define DEBUG 1 // debug log - lvls: 0 - none, 1 - minimal, 2 - all

#define MAX_DEGREE 8
#define SYMBOL_MAX 256

// common bit fields
#define SPACE_READ  0x01 
#define DIGIT_READ  0x02 
#define DASH_READ   0x04

// standard parse fields
#define RING_READ   0x08
#define BIND_READ   0x10 // locant_t ring access
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


static bool 
error(const char *message) 
{
  fprintf(stderr,"%s\n", message); 
  return false;  
}

#ifdef USING_OPENBABEL

#define graph_t  OpenBabel::OBMol
#define atom_t   OpenBabel::OBAtom
#define edge_t   OpenBabel::OBBond

#elif defined USING_RDKIT

#endif


typedef struct locant_t {
  uint8_t hloc; // used for pathsolver 
  uint8_t r_pack; // [of 1b][ offL 1b ][ offR 1b ][ bridging 1b ][ dangling u4 ] 
  atom_t *s; 
  locant_t *off_path[2]; // 0 = -. 1 = '&' 
} locant_t; 

/* can be called inline with other state tracks */
static locant_t* 
next_in_locant_tree(locant_t *locant, char nxt)
{
  if (!locant)
    return 0; 
  else
   return locant->off_path[(nxt=='-')]; 
}

typedef struct {
  uint8_t  size; 
  locant_t path[1]; // malloc sizeof(locant_t) * (size-1) + (1 byte for size)
} ring_t;  


static atom_t* 
new_symbol(graph_t *mol, uint64_t atomic_num)
{
  atom_t *atom; 
#if USING_OPENBABEL
  atom = mol->NewAtom();
  atom->SetAtomicNum(atomic_num); 
#elif defined USING_RDKIT

#endif
  return atom; 
}

static bool 
assign_symbol(atom_t *s, uint16_t atomic_num)
{
#if USING_OPENBABEL
  s->SetAtomicNum(atomic_num); 
#elif defined USING_RDKIT

#endif
  return true; 
}


static bool
assign_dash_symbol(atom_t *s, uint8_t fst_ch, uint8_t snd_ch) 
{
  switch (fst_ch){
    case 'A':
      switch (snd_ch) {
        case 'C': return assign_symbol(s,89); 
        case 'G': return assign_symbol(s,47); 
        case 'L': return assign_symbol(s,13); 
        case 'M': return assign_symbol(s,95); 
        case 'R': return assign_symbol(s,18); 
        case 'S': return assign_symbol(s,33); 
        case 'T': return assign_symbol(s,85); 
        case 'U': return assign_symbol(s,79); 
      }
      break; 

    case 'B':
      switch (snd_ch) {
        case 0: return assign_symbol(s,BOR); 
        case 'A': return assign_symbol(s,56); 
        case 'E': return assign_symbol(s,4); 
        case 'H': return assign_symbol(s,107); 
        case 'I': return assign_symbol(s,83); 
        case 'K': return assign_symbol(s,97); 
        case 'R': return assign_symbol(s,BRO); 
      }
      break; 

    case 'C':
      switch (snd_ch) {
        case 0: return assign_symbol(s,CAR); // allow the texas carbon "easter )egg"
        case 'A': return assign_symbol(s,20);
        case 'D': return assign_symbol(s,48);
        case 'E': return assign_symbol(s,58);
        case 'F': return assign_symbol(s,98);
        case 'M': return assign_symbol(s,96);
        case 'N': return assign_symbol(s,112);
        case 'O': return assign_symbol(s,27);
        case 'R': return assign_symbol(s,24);
        case 'S': return assign_symbol(s,55);
        case 'U': return assign_symbol(s,29);
      }
      break; 

    case 'D':
      switch (snd_ch) {
        case 'B': return assign_symbol(s,105);
        case 'S': return assign_symbol(s,110);
        case 'Y': return assign_symbol(s,66);
      }
      break;

    case 'E':
      switch (snd_ch) {
        case 0: return assign_symbol(s,35); 
        case 'R': return assign_symbol(s,68); 
        case 'S': return assign_symbol(s,99); 
        case 'U': return assign_symbol(s,63); 
      }

    case 'F':
      switch (snd_ch) {
        case 0: return assign_symbol(s,FLU); 
        case 'E': return assign_symbol(s,26); 
        case 'L': return assign_symbol(s,114); 
        case 'M': return assign_symbol(s,100); 
        case 'R': return assign_symbol(s,87); 
      }
      break;

    case 'G':
      switch (snd_ch) {
        case 0:  return assign_symbol(s,CHL); 
        case 'A': return assign_symbol(s,31); 
        case 'D': return assign_symbol(s,64); 
        case 'E': return assign_symbol(s,32); 
      }
      break;

    case 'H':
      switch (snd_ch) {
        case 'E': return assign_symbol(s,2); 
        case 'F': return assign_symbol(s,72); 
        case 'G': return assign_symbol(s,80); 
        case 'O': return assign_symbol(s,67); 
        case 'S': return assign_symbol(s,108); 
      }
      break;

    case 'I':
      switch (snd_ch) {
        case 0:  return assign_symbol(s,IOD);
        case 'N': return assign_symbol(s,49); 
        case 'R': return assign_symbol(s,77); 
      }
      break;

    case 'K':
      switch (snd_ch) {
        case 0:   return assign_symbol(s,NIT);
        case 'R': return assign_symbol(s,36); 
        case 'A': return assign_symbol(s,19); 
      }
      break;

    case 'L':
      switch (snd_ch) {
        case 'A': return assign_symbol(s,57); 
        case 'I': return assign_symbol(s,3); 
        case 'R': return assign_symbol(s,103); 
        case 'U': return assign_symbol(s,71); 
        case 'V': return assign_symbol(s,116); 
      }
      break;

    case 'M':
      switch (snd_ch) {
        case 0: return assign_symbol(s,NIT); 
        case 'C': return assign_symbol(s,115); 
        case 'D': return assign_symbol(s,101); 
        case 'G': return assign_symbol(s,12); 
        case 'N': return assign_symbol(s,25); 
        case 'O': return assign_symbol(s,42); 
        case 'T': return assign_symbol(s,109); 
      }
      break;

    case 'N':
      switch (snd_ch) {
        case 0: return assign_symbol(s,NIT); 
        case 'A': return assign_symbol(s,11); 
        case 'B': return assign_symbol(s,41); 
        case 'D': return assign_symbol(s,60); 
        case 'E': return assign_symbol(s,10); 
        case 'H': return assign_symbol(s,113); 
        case 'I': return assign_symbol(s,28); 
        case 'O': return assign_symbol(s,102); 
        case 'P': return assign_symbol(s,93); 
      }
      break; 

    case 'O':
      switch (snd_ch) {
        case 0: return assign_symbol(s,OXY); 
        case 'G': return assign_symbol(s,118); 
        case 'S': return assign_symbol(s,76); 
      }
      break;

    case 'P':
      switch (snd_ch) {
        case 0: return assign_symbol(s,PHO); 
        case 'A': return assign_symbol(s,91); 
        case 'B': return assign_symbol(s,82); 
        case 'D': return assign_symbol(s,46); 
        case 'M': return assign_symbol(s,61); 
        case 'O': return assign_symbol(s,84); 
        case 'R': return assign_symbol(s,59); 
        case 'T': return assign_symbol(s,78); 
        case 'U': return assign_symbol(s,94); 
      }
      break;
    
    case 'Q':
      switch (snd_ch) {
        case 0: return assign_symbol(s,OXY); 
      }
      break; 

    case 'R':
      switch (snd_ch) {
        case 'A': return assign_symbol(s,88); 
        case 'B': return assign_symbol(s,37); 
        case 'E': return assign_symbol(s,75); 
        case 'F': return assign_symbol(s,104); 
        case 'G': return assign_symbol(s,111); 
        case 'H': return assign_symbol(s,45); 
        case 'N': return assign_symbol(s,86); 
        case 'U': return assign_symbol(s,44); 
      }
      break;

    case 'S':
      switch (snd_ch) {
        case 0:  return assign_symbol(s,SUL); 
        case 'B': return assign_symbol(s,51); 
        case 'C': return assign_symbol(s,21); 
        case 'E': return assign_symbol(s,34); 
        case 'G': return assign_symbol(s,106); 
        case 'I': return assign_symbol(s,14); 
        case 'M': return assign_symbol(s,62); 
        case 'N': return assign_symbol(s,50); 
        case 'R': return assign_symbol(s,38); 
      }
      break;

    case 'T':
      switch (snd_ch) {
        case 'A': return assign_symbol(s,73); 
        case 'B': return assign_symbol(s,65); 
        case 'C': return assign_symbol(s,43); 
        case 'E': return assign_symbol(s,52); 
        case 'H': return assign_symbol(s,90); 
        case 'I': return assign_symbol(s,22); 
        case 'L': return assign_symbol(s,81); 
        case 'M': return assign_symbol(s,69); 
        case 'S': return assign_symbol(s,117); 
      }
      break;

    case 'U':
      switch (snd_ch) {
        case 'R': return assign_symbol(s,92); 
      }
      break;

    case 'V':
      switch (snd_ch) {
        case 'A': return assign_symbol(s,23); 
      }
      break;
  
    case 'W':
      switch (snd_ch) {
        case 'T': return assign_symbol(s,74); 
      }
      break; 

    case 'X':
      switch (snd_ch) {
        case 'E': return assign_symbol(s,54); 
      }
      break;

    case 'Y':
      switch (snd_ch) {
        case 'T': return assign_symbol(s,39); 
        case 'B': return assign_symbol(s,70); 
      }
      break;

    case 'Z':
      switch (snd_ch) {
        case 'N': return assign_symbol(s,30); 
        case 'R': return assign_symbol(s,30); 
      }
      break;
  }
  
  return false; 
}


// create a bond to a dummy atom type. allows bond modification without state hold 
static edge_t* 
add_bond(graph_t *g, atom_t *a, atom_t *b)
{
  edge_t *bptr = (edge_t*)0;
#ifdef USING_OPENBABEL
  if (!g->AddBond(a->GetIdx(), b->GetIdx(), 1)) {
    fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n", a->GetIdx(),b->GetIdx());
    return bptr;
  }
  else     
    bptr = g->GetBond(g->NumBonds() - 1);
#elif defined USING_RDKIT

#endif
  return bptr;
}


static bool
unsaturate_edge(edge_t *e)
{
  if (!e)
    return false; 
#ifdef USING_OPENBABEL
  if (e->GetBondOrder() < 3)
    e->SetBondOrder(e->GetBondOrder() + 1); 
  else return false; 
#elif USING_RDKIT

#endif
  return true; 
}

typedef struct  {
  uint8_t r_loc; // use to calculate size + add in off-branch positions
  uint8_t r_size;
  uint8_t arom; 
} r_assignment; 

#if 0
static ring_t* rt_alloc(graph_t *g, const size_t size, atom_t *inc, const uint16_t inc_pos) 
{
  ring_t *ring = 0; 
  ring = (ring_t*)malloc(sizeof(ring_t) + sizeof(locant_t)*(size-1)); 
  memset(ring, 0, sizeof(ring_t) + sizeof(locant_t)*(size-1)); 
  
  atom_t *c = 0;  
  atom_t *p = 0;
  edge_t *e   = 0; 

  for (uint16_t i=0; i<size; i++) {
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
    ring->path[i].hloc = i+1; // point to the next locant_t 
    p = c; 
  }
  
  // end chain movement for pathfinderIII
  ring->path[0].r_pack++; 
  ring->path[size-1].r_pack++; 
  ring->path[size-1].hloc = size-1; 
  
  ring->size = size; 
  return ring; 
}


static ring_t* pathsolverIII_fast(graph_t *g, ring_t *r, 
                                   r_assignment *SSSR, uint8_t SSSR_ptr) 
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

  uint8_t steps, s_pos;
  edge_t *e; 
  locant_t *start, *end; 
  r_assignment *subcycle; 
  
  for (uint16_t i=0; i<SSSR_ptr; i++) {
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
    
    for (uint16_t s=0; s<steps-1; s++) {
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

static ring_t* 
pathsolverIII(graph_t *g, ring_t *r, 
              r_assignment *SSSR, uint8_t SSSR_ptr, 
              uint8_t *connection_table)
{
  uint8_t steps, s_pos;
  edge_t *e; 
  locant_t *c, *p=0;
  locant_t *start, *end; 
  r_assignment *subcycle; 

  for (uint16_t i=0; i<r->size; i++) {
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
    c->hloc = i+1; // point to the next locant_t 
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
#endif

static bool 
add_oxy(graph_t *mol, atom_t *atom)
{
  atom_t *oxygen = new_symbol(mol, OXY); 
  edge_t *bptr = add_bond(mol, atom, oxygen);
  if (!bptr)  
    return false;
  else if (!unsaturate_edge(bptr))
    return false; 
  return true; 
}

#if 0
/* placeholder for charged alterations */ 
static bool 
add_tauto_dioxy(graph_t *g, atom_t *p)
{
  bool ret = add_oxy(g, p, 0); 
  if (ret)
    return ret; 
  if (((p->valence_pack >> 4) - (p->valence_pack & 0x0F)) >= 2) 
    return add_oxy(g, p, 0); 
  else 
    return add_oxy(g, p, 1); 
};
#endif

/* locants can be expanded or read as branches with '&' and '-'
 * chars. In order to move the state as the string is being read, 
 * these have to be handled per locant_t read. The logic this saves
 * outweighs the branches in the loop */
static uint8_t 
read_locant(const char *ptr, uint16_t *idx, uint16_t limit_idx)
{
  uint8_t locant_t = ptr[*idx] - 'A'; 
  for (uint16_t i=*idx+1; i<limit_idx; i++) {
    switch (ptr[i]) {
      case '&':
        locant_t += 23; 
        break; 
      case '-':
        fprintf(stderr,"off-branch reads need handling\n"); 
      default:
        return locant_t; 
    }
    (*idx)++; 
  }
  return locant_t; 
}

#if 0
static ring_t* 
parse_cyclic(const char *ptr, 
             uint16_t start, 
             uint16_t end, 
             atom_t *head, 
             uint16_t head_loc, 
             graph_t *g) 
{
  atom_t  *c    = 0; 
  edge_t    *e    = 0; 
  ring_t    *ring = 0; 

  uint8_t   locant_ch   = 0; 
  uint8_t   arom_count = 0;
  uint8_t   ring_size  = 0;  
  
  uint8_t SSSR_ptr  = 0; 
  r_assignment SSSR[32]; // this really is sensible for WLN
  
  uint8_t bridge_ptr  = 0; 
  uint16_t bridge_locants[16]; 

  uint8_t buff_ptr = 0; 
  unsigned char buffer[3]; 
    
  uint8_t *connection_table = 0; // stack allocates on pseudo read

  uint8_t state = 0; // bit field:                      
                 // [][][][pseudo][SSSR][dash][digit][space]

  unsigned char ch; 
  unsigned char ch_nxt; 

  state = SSSR_READ; 
  for (uint16_t sp=start; sp<end; sp++){
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
          fprintf(stderr, "Error: out of bounds locant_t access"); 
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
      if (state & SSSR_READ && (ch >= 'A' && ch <= 'Z')){
        // end the SSSR read, start element reading
        if (ring || !SSSR_ptr) {
          fprintf(stderr, "Error: ring backbone failure - syntax ordering issue\n"); 
          goto ring_fail; 
        }
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

            for (uint8_t i=0; i < ch-'0'+1; i++) 
              read_locant(ptr, &++sp, end); 

            if (sp >= end || ptr[sp] != ' ') {
              fprintf(stderr, "Error: invalid format for multicyclic ring - %c\n", ptr[sp]);
              goto ring_fail; 
            }

            state |= MULTI_READ; 
            state &= ~SPACE_READ; 
            state &= ~SSSR_READ; 
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
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, BOR, 3); 
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'H':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = ring->path[locant_ch].s; 
          c->arom = -1; // will get ignored, bits overlap
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'K': 
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, NIT, 4); 
          ring->path[locant_ch].r_pack++; 
          c->charge++; 
          
          c->bonds[0].order = HANGING_BOND; 
          c->bonds[1].order = HANGING_BOND;  
          c->bonds[2].order = HANGING_BOND; 

          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'M':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, NIT, 2); 
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'N':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, NIT, 3); 
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'O':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, OXY, 2); 
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'P':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, PHO, 6); 
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'S':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, SUL, 6); 
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 
  
        case 'T':
          state |= AROM_READ; 
          SSSR[arom_count++].arom = 0; 
          break; 

        // some nice ghost logic here. On a successful chain creation, the first bond will always to 
        // to the next symbol in the locant_t path. e.g A->B->C. therefore unsaturate bond[0]. 
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
            e = &c->bonds[0]; 
            e->order++; 
            c->valence_pack++;
            ring->path[locant_ch+1].s->valence_pack++; 
            locant_ch++; 
          }
          break; 
          
        case 'V':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = transmute_symbol(ring->path[locant_ch].s, CAR, 4);
          e = &c->bonds[0]; 
          if (!add_oxy(g, c, 0))
            goto ring_fail; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 

        case 'W':
          break; 

        case 'X':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          c = ring->path[locant_ch].s; 
          // purely virtual bonds that may or may not get tidyed up
          c->bonds[1].order = HANGING_BOND; 
          c->bonds[2].order = HANGING_BOND; 
          c->bonds[3].order = HANGING_BOND; 

          e = &c->bonds[0]; 
          ring->path[locant_ch].r_pack++; // the me
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
          break; 
      
        case 'Y':
          if (locant_ch > ring_size) {
            fprintf(stderr,"Error: out of bounds locant_t access\n"); 
            goto ring_fail; 
          }
          // treat as a carbon, WLN expects the U characters, this cannot be aromatic 
          c = ring->path[locant_ch].s; 
          c->arom = -1; 
          e = &c->bonds[0]; 
          g->idx_symbols[sp+1] = c; 
          locant_ch++; 
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
  for (uint16_t i=0; i<bridge_ptr; i++) {
    if (bridge_locants[i] >= ring->size) {
      fprintf(stderr,"Error: out of bounds locant_t access\n"); 
      goto ring_fail; 
    } 
    ring->path[bridge_locants[i]].r_pack--; 
  }

  // forward any single arom assignments
  for (uint16_t i=arom_count;i<SSSR_ptr;i++)
    SSSR[i].arom = SSSR[0].arom; 
  return pathsolverIII_fast(g, ring, SSSR, SSSR_ptr);     

ring_fail:
  if (ring)
    free(ring);
  return (ring_t*)0; 
}
#endif


#if 0
/* returns the pending unsaturate level for the next symbol */
static edge_t* 
resolve_unsaturated_carbon(graph_t *g, 
                           edge_t *e, 
                           atom_t *p, 
                           atom_t *c, 
                           unsigned char ch_nxt)
{
  edge_t *ne = next_virtual_edge(c); // forward edge from 'C' 

  // first check can you double bond to the previous - 
  // the bond will be set to 1 so its >= 1
  if ((p->valence_pack>>4) - (p->valence_pack & 0x0F) >= 1) {
    // next predict what bonding is possible to whats coming next
    switch (ch_nxt) {
      // alkyl chains have a weird relationship where they can but shouldn't unsaturate here - 
      // murky waters with implied chemical information
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
#endif


/*
 * -- Parse WLN Notation --
 */
bool 
ReadWLN(const char *wln, graph_t *molecule)
{
  edge_t *curr_edge=0; 
  atom_t *curr_atom=0;
  atom_t *dummy_atom=0;
  ring_t *curr_ring=0;

  uint8_t ring_chars = 0; 
  uint16_t locant_ch = 0; 
  uint16_t digit_n   = 0; 
  
  uint8_t state = 0; // bit field: 
                     // [][dioxo][charge][ring locant_t][ring skip][dash][digit][space]
  
  uint8_t dash_ptr = 0; 
  unsigned char dash_chars[3] = {0}; // last byte is mainly for overflow 
 
  // charges are assigned through string indexing - tf: need a strlen() 
    
  uint8_t stack_ptr = 0;  // WLN branch and ring dependency stack
  struct wlnrefaddr {
    void *addr; 
    char ref; // -1 for (ring_t*) else (atom_t*) 
  } dep_stack[64];  

#if 0
  // seperated in order to reuse the symbols memory on file read
  for (uint16_t i=0;i<stack_ptr;i++) {
    // if (stack[i].ref < 0)  
    //   free(stack[i].addr);
    stack[i].addr = 0; 
    stack[i].ref = 0;
  }
#endif

  // init conditions, make one dummy atom, and one bond - work of the virtual bond
  // idea entirely *--> grow...
  curr_atom = new_symbol(molecule, DUM); 
  
  unsigned char ch; 
  unsigned char ch_nxt; 
  const uint32_t len = strlen(wln);  
  for (uint16_t sp=0; sp<len; sp++) {
    ch     = wln[sp]; 
    ch_nxt = wln[sp+1]; // one lookahead is defined behaviour 
    
    if (state & RING_READ) {
#if 0
      if (wln[sp-1] >= 'A' && ch == 'J' && ch_nxt < '0') {
        // J can be used inside ring notation, requires lookahead 
        // condition 
        
        // note: The ptr passed in does not include the 
        //       starting L/T or ending J (<sp) 
        
        // ring must have at least one atom. 
        curr_atom = assign_symbol(g, curr_edge, CAR, 4);
        if (!curr_atom)
          return false; 
        curr_edge = set_virtual_edge(curr_edge, prev_atom, curr_atom); 
        curr_ring = parse_cyclic(wln, sp-ring_chars+1, sp, curr_atom, locant_ch, g);
        if (!curr_ring)
          return false; 
        g->stack[g->stack_ptr].addr = curr_ring; 
        g->stack[g->stack_ptr++].ref = -1; 
        state &= ~(RING_READ);
        ring_chars = 0; 
        state &= ~(BIND_READ); 
      }
      else
        ring_chars++; 
    }
    else if (state & DASH_READ && ch >= 'A' && ch <= 'Z') {
      if (dash_ptr == 3) 
        return error("Error: elemental code can only have 2 character symbols"); 
      dash_chars[dash_ptr++] = ch; 
    }
    else if (state & SPACE_READ) {
      // switch on packing - state popcnt() >= hit bits. 
      locant_ch = ch - 'A'; 
      state &= ~(SPACE_READ);

      if (!(state & BIND_READ)) {
        if (curr_ring && locant_ch < curr_ring->size) {
          curr_atom = curr_ring->path[locant_ch].s; 
          curr_edge = next_virtual_edge(curr_atom); 
          prev_atom = curr_atom; 
        } 
        else 
          return error("Error: out of bounds locant_t access"); 
      }
#endif
    }
    else {
      switch (ch) {
        case '0':
          digit_n *= 10; 

          if(state & DIOXO_READ) 
            return error("Error: dioxo attachment needs higher valence atom"); 
          else if (ch_nxt == '/') {
#if 0
            // positive charge assignment
            if (digit_n > len || !g->idx_symbols[digit_n]) 
              return error("Error: charge assignment out of bounds"); 
            else {
              g->idx_symbols[digit_n]->charge++; 
              state |= CHARGE_READ; 
              digit_n = 0; 
              sp++; // skip the slash, only ever used again in R notation
            }
#endif
          } 
          else if (state & (CHARGE_READ | DIGIT_READ) 
                   && (ch_nxt < '0' || ch_nxt > '9')) 
          {  
#if 0
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
              for (uint16_t i=0; i<digit_n; i++) {
                curr_atom = assign_symbol(g, curr_edge, CAR, 4);
                if (!curr_atom)
                  return false; 
                else
                  curr_edge = set_virtual_edge(curr_edge, prev_atom, curr_atom); 

                if (!curr_edge)
                  return false; 
                else {
                  prev_atom = curr_atom;
                  curr_edge = next_virtual_edge(prev_atom); 
                }
              }

              g->idx_symbols[sp+1] = prev_atom; 
              digit_n = 0; 
            }
#endif
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
#if 0
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
              for (uint16_t i=0; i<digit_n; i++) {
                curr_atom = assign_symbol(g, curr_edge, CAR, 4);
                if (!curr_atom)
                  return false; 
                else
                  curr_edge = set_virtual_edge(curr_edge, prev_atom, curr_atom); 

                if (!curr_edge)
                  return false; 
                else {
                  prev_atom = curr_atom;
                  curr_edge = next_virtual_edge(prev_atom); 
                }
              }
              g->idx_symbols[sp+1] = prev_atom; 
              digit_n = 0; 
            }
#endif
          }
          else
            state |= DIGIT_READ; 
          break;
        
        case 'A':
        case 'J':
          return error("Error: non-atomic symbol used in chain"); 
        
        case 'B':
          assign_symbol(curr_atom, BOR);
          
          dep_stack[stack_ptr].addr = curr_atom; 
          dep_stack[stack_ptr].ref  = 2;
          stack_ptr++; 

          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break; 
        
        case 'C':
          fprintf(stderr, "WLN symbol C currently unhandled\n"); 
          return false; 
          break; 

        // extrememly rare open chelate notation
        case 'D':
          fprintf(stderr,"chelate ring open needs handling\n"); 
          break; 
        
        // terminator symbol - no bond movement
        case 'E':
          assign_symbol(curr_atom, BRO);
          // note: terminators can act on an empty stack, '&' cannot
          if (!stack_ptr) {
            // the molecule has to end here. open charge notation
            if (ch_nxt != ' ' || ch_nxt != '\0') 
              return error("Error: terminating symbol closes molecule\n"); 
            else 
              state |= CHARGE_READ; 
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;

        case 'F':
          assign_symbol(curr_atom, BRO);
          // note: terminators can act on an empty stack, '&' cannot
          if (!stack_ptr) {
            // the molecule has to end here. open charge notation
            if (ch_nxt != ' ' || ch_nxt != '\0') 
              return error("Error: terminating symbol closes molecule\n"); 
            else 
              state |= CHARGE_READ; 
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;

        case 'G':
          assign_symbol(curr_atom, CHL);
          // note: terminators can act on an empty stack, '&' cannot
          if (!stack_ptr) {
            // the molecule has to end here. open charge notation
            if (ch_nxt != ' ' || ch_nxt != '\0') 
              return error("Error: terminating symbol closes molecule\n"); 
            else 
              state |= CHARGE_READ; 
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;

        case 'H':
          assign_symbol(curr_atom, 1);
          // note: terminators can act on an empty stack, '&' cannot
          if (!stack_ptr) {
            // the molecule has to end here. open charge notation
            if (ch_nxt != ' ' || ch_nxt != '\0') 
              return error("Error: terminating symbol closes molecule\n"); 
            else 
              state |= CHARGE_READ; 
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;

        case 'I':
          assign_symbol(curr_atom, IOD);
          // note: terminators can act on an empty stack, '&' cannot
          if (!stack_ptr) {
            // the molecule has to end here. open charge notation
            if (ch_nxt != ' ' || ch_nxt != '\0') 
              return error("Error: terminating symbol closes molecule\n"); 
            else 
              state |= CHARGE_READ; 
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;
       
        /* [N+](R)(R)(R)(R) */
        case 'K':
          assign_symbol(curr_atom, NIT);
          
          dep_stack[stack_ptr].addr = curr_atom; 
          dep_stack[stack_ptr].ref  = 3;
          stack_ptr++; 

          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;

        case 'L':
        case 'T':
          state |= RING_READ; 
          ring_chars++; 
          break; 
        
        /* NH(R)(R) */
        case 'M':
          assign_symbol(curr_atom, NIT);
          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;
       
        /* NR(R)(R) */
        case 'N':
          assign_symbol(curr_atom, NIT);
          
          dep_stack[stack_ptr].addr = curr_atom; 
          dep_stack[stack_ptr].ref  = 2;
          stack_ptr++; 

          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;
        
        /* OR(R) */
        case 'O':
          assign_symbol(curr_atom, OXY);
          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;
        
        case 'P':
          assign_symbol(curr_atom, PHO);
          
          dep_stack[stack_ptr].addr = curr_atom; 
          dep_stack[stack_ptr].ref  = 3;
          stack_ptr++; 
          
          dummy_atom = new_symbol(molecule, DUM);
          if (!dummy_atom) {
            fprintf(stderr, "wahat? %p - %p\n", curr_atom, dummy_atom); 
            return false; 
          }

          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;

        case 'Q':
          assign_symbol(curr_atom, OXY);
          // note: terminators can act on an empty stack, '&' cannot
          if (!stack_ptr) {
            // the molecule has to end here. open charge notation
            if (ch_nxt != ' ' || ch_nxt != '\0') 
              return error("Error: terminating symbol closes molecule\n"); 
            else 
              state |= CHARGE_READ; 
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;

        // shorthand benzene 
#if 0
        case 'R': 
          // head symbol can come from virtual X|Y|K
          curr_atom = assign_symbol(g, curr_edge, CAR, 4);
          if (!curr_atom)
            return false; 
          else {
            g->idx_symbols[sp+1] = curr_atom; 
            curr_edge = set_virtual_edge(curr_edge, prev_atom, curr_atom); 
          }

          curr_ring = rt_alloc(g, 6, curr_atom, 0); 
          
          if (!curr_ring)
            return false; // only way this can fail 
          else {
            for (uint8_t i=0; i<6; i++)
              curr_ring->path[i].s->arom |= 1; 
            // set the ring. R objects are kekulised on stack pop
            prev_atom = curr_ring->path[0].s; 
            curr_atom = curr_ring->path[5].s; 
            curr_edge = next_virtual_edge(prev_atom); 
            curr_edge = set_virtual_edge(curr_edge, prev_atom, curr_atom);
            curr_edge->ring = 1; 
          }

          // add to ring stack
          g->stack[g->stack_ptr].addr = curr_ring; 
          g->stack[g->stack_ptr++].ref = -1; 

          g->idx_symbols[sp+1] = prev_atom; 
          curr_edge = next_virtual_edge(prev_atom); 
          break; 
#endif
        
        case 'S':
          assign_symbol(curr_atom, SUL);
          
          dep_stack[stack_ptr].addr = curr_atom; 
          dep_stack[stack_ptr].ref  = 3;
          stack_ptr++; 

          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;
        
        case 'U':
          if (!unsaturate_edge(curr_edge))
            return error("Error: failed to unsaturate bond\n"); 
          break; 

        case 'V':
          assign_symbol(curr_atom, CAR);
          if (!add_oxy(molecule, curr_atom))
            return error("Error: failed to add =O group\n"); 
        
          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          if (!curr_edge)
            return error("Error: failed to add =O group\n"); 
          curr_atom = dummy_atom;
          break;
        
        // if not previous, create dummy carbon
        // really trying to avoid a bit state on this
        case 'W':
          return error("Error: W group needs supporting"); 
          break; 

        case 'X':
          assign_symbol(curr_atom, CAR);
          
          dep_stack[stack_ptr].addr = curr_atom; 
          dep_stack[stack_ptr].ref  = 3;
          stack_ptr++; 

          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;

        case 'Y':
          assign_symbol(curr_atom, CAR);
          
          dep_stack[stack_ptr].addr = curr_atom; 
          dep_stack[stack_ptr].ref  = 2;
          stack_ptr++; 

          dummy_atom = new_symbol(molecule, DUM); 
          curr_edge = add_bond(molecule, dummy_atom, curr_atom); 
          curr_atom = dummy_atom;
          break;
        
        case 'Z':
          assign_symbol(curr_atom, NIT);
          // note: terminators can act on an empty stack, '&' cannot
          if (!stack_ptr) {
            // the molecule has to end here. open charge notation
            if (ch_nxt != ' ' || ch_nxt != '\0') 
              return error("Error: terminating symbol closes molecule\n"); 
            else 
              state |= CHARGE_READ; 
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;
        
        case '-':
          return error("Error: dash code needs hooking in\n"); 
#if 0
          // dashes open up several possible states which needs a lot of information to reason 
          // 1. -XX- | -X- == symbol
          // 2. <locant> -& <space> <locant_t> <ring> == spiro ring defintion  
          // 3. <locant> (-|&)+ == high level access
          // ?  <locant> (-&) space <no ring> == methyl shorthand on the level access (solved with VIRTUAL edge packing)
          // seperating out 1 is trivial with a +1 lookahead. 
          // 2 vs 3: 
          if (state & DASH_READ) {
            curr_atom = assign_symbol(g, curr_edge, DUM, 6); // create a blank symbol
            if (!curr_atom)
              return false; 
            if (write_dash_symbol(curr_atom, dash_chars[0], dash_chars[1]) == false)
              return false; 

            curr_edge = set_virtual_edge(curr_edge, prev_atom, curr_atom); 
            if(state & DIOXO_READ) {
              add_tauto_dioxy(g, curr_atom); 
              state &= ~DIOXO_READ; 
            }
            
            memset(dash_chars,0,3); 
            state &= ~DASH_READ; 
            dash_ptr = 0; 

            // add to branch stack
            if ((curr_atom->valence_pack >> 4) >= 2) {
              g->stack[g->stack_ptr].addr = curr_atom; 
              g->stack[g->stack_ptr].ref  = curr_atom->valence_pack >> 4;
              g->stack_ptr++; 
            }

            prev_atom = curr_atom; 
            g->idx_symbols[sp+1] = curr_atom; 
            curr_edge = next_virtual_edge(curr_atom); 
          }
          else {
            if (ch_nxt == ' ')
              state |= BIND_READ;
            else {
              dash_ptr = 0; 
              state |= DASH_READ; 
            }
          }
#endif
          break; 

        case ' ':
          return error("Error: space code needs hooking in\n"); 
#if 0
          if (ch_nxt == '&') {
            // all other states should make this the only place ions can be used. 
            gt_stack_flush(g); 
            sp++; 
            curr_atom = new_symbol(g, DUM, 1); 
            curr_edge = next_virtual_edge(curr_atom); 
            curr_edge->order = DUM; 
            prev_atom = curr_atom; 
          }
          else 
            state |= SPACE_READ; 
#endif
          break; 

        case '&':
          if (state & DIOXO_READ) 
            return error("Error: dioxo attachment requires an atomic symbol"); 
          if (!stack_ptr) 
            return error("Error: empty dependency stack - too many &?"); 
                  
          if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            dummy_atom = new_symbol(molecule, DUM); 
            curr_edge = add_bond(molecule, (atom_t*)dep_stack[stack_ptr-1].addr, dummy_atom); 
            curr_atom = dummy_atom; 

            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
          break;

        case '\n':
          break; 
        case '/':
          return error("Error: slash seen outside of ring - multipliers currently unsupported");
        default:
          fprintf(stderr, "Error: invalid character read for WLN notation - %c\n", ch);
          return false; 
      }
    }
  }
  
  molecule->SetAromaticPerceived(false);
  molecule->AddHydrogens(); 
  return true; 
}

