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
#include <openbabel/obiter.h>

#include <openbabel/babelconfig.h>
#include <openbabel/obconversion.h>

#include "wlnparser.h"

//#define DEBUG 1 // debug log - lvls: 0 - none, 1 - minimal, 2 - all

#define MAX_DEGREE 8
#define SYMBOL_MAX 256

// common states 
#define INIT_READ   0x00
#define CLOSE_READ  0x01
#define SPACE_READ  0x02
#define ALKYL_READ  0x04 
#define DASH_READ   0x08

// ReadWLN states 
#define RING_READ   0x10
#define BIND_READ   0x20 // locant_t ring access
#define CHARGE_READ 0x40

// ring_parse states 
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
using namespace OpenBabel;
#define graph_t    OBMol
#define symbol_t   OBAtom
#define edge_t     OBBond

#define symbol_get_num(s)           s->GetAtomicNum()
#define symbol_change_num(s, n)     s->SetAtomicNum(n)
#define symbol_set_hydrogens(s, n)  s->SetImplicitHCount(n)
#define symbol_incr_hydrogens(s)    s->SetImplicitHCount(s->GetImplicitHCount()+1)

#define edge_unsaturate(e)        e->SetBondOrder(e->GetBondOrder()+1)

#define graph_delete_symbol(g,s)         g->DeleteAtom(s)

#elif defined USING_RDKIT

#endif

typedef struct locant_t {
  uint8_t hloc; // used for pathsolver 
  uint8_t r_pack; // [of 1b][ offL 1b ][ offR 1b ][ bridging 1b ][ dangling u4 ] 
  symbol_t *s; 
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


static symbol_t* 
symbol_create(graph_t *mol, uint64_t atomic_num)
{
  symbol_t *atom; 
#if USING_OPENBABEL
  atom = mol->NewAtom();
  atom->SetAtomicNum(atomic_num); 
#elif defined USING_RDKIT

#endif
  return atom; 
}


static bool
symbol_create_from_dash(symbol_t *s, uint8_t fst_ch, uint8_t snd_ch) 
{
  switch (fst_ch){
    case 'A':
      switch (snd_ch) {
        case 'C': symbol_change_num(s,89); return true;
        case 'G': symbol_change_num(s,47); return true;
        case 'L': symbol_change_num(s,13); return true;
        case 'M': symbol_change_num(s,95); return true;
        case 'R': symbol_change_num(s,18); return true;
        case 'S': symbol_change_num(s,33); return true;
        case 'T': symbol_change_num(s,85); return true;
        case 'U': symbol_change_num(s,79); return true;
      }
      break; 

    case 'B':
      switch (snd_ch) {
        case 0:   symbol_change_num(s,BOR); return true;
        case 'A': symbol_change_num(s,56);  return true;
        case 'E': symbol_change_num(s,4); return true;
        case 'H': symbol_change_num(s,107); return true;
        case 'I': symbol_change_num(s,83); return true;
        case 'K': symbol_change_num(s,97); return true;
        case 'R': symbol_change_num(s,BRO); return true;
      }
      break; 

    case 'C':
      switch (snd_ch) {
        case 0:   symbol_change_num(s,CAR); return true;
        case 'A': symbol_change_num(s,20);return true;
        case 'D': symbol_change_num(s,48);return true;
        case 'E': symbol_change_num(s,58);return true;
        case 'F': symbol_change_num(s,98);return true;
        case 'M': symbol_change_num(s,96);return true;
        case 'N': symbol_change_num(s,112);return true;
        case 'O': symbol_change_num(s,27);return true;
        case 'R': symbol_change_num(s,24);return true;
        case 'S': symbol_change_num(s,55);return true;
        case 'U': symbol_change_num(s,29);return true;
      }
      break; 

    case 'D':
      switch (snd_ch) {
        case 'B': symbol_change_num(s,105);return true;
        case 'S': symbol_change_num(s,110);return true;
        case 'Y': symbol_change_num(s,66);return true;
      }
      break;

    case 'E':
      switch (snd_ch) {
        case 0:   symbol_change_num(s,35); return true;
        case 'R': symbol_change_num(s,68); return true;
        case 'S': symbol_change_num(s,99); return true;
        case 'U': symbol_change_num(s,63); return true;
      }

    case 'F':
      switch (snd_ch) {
        case 0:   symbol_change_num(s,FLU); return true;
        case 'E': symbol_change_num(s,26); return true;
        case 'L': symbol_change_num(s,114); return true;
        case 'M': symbol_change_num(s,100); return true;
        case 'R': symbol_change_num(s,87); return true;
      }
      break;

    case 'G':
      switch (snd_ch) {
        case 0:   symbol_change_num(s,CHL); return true;
        case 'A': symbol_change_num(s,31); return true;
        case 'D': symbol_change_num(s,64); return true;
        case 'E': symbol_change_num(s,32); return true;
      }
      break;

    case 'H':
      switch (snd_ch) {
        case 'E': symbol_change_num(s,2); return true;
        case 'F': symbol_change_num(s,72); return true;
        case 'G': symbol_change_num(s,80); return true;
        case 'O': symbol_change_num(s,67); return true;
        case 'S': symbol_change_num(s,108); return true;
      }
      break;

    case 'I':
      switch (snd_ch) {
        case 0:  symbol_change_num(s,IOD);return true;
        case 'N': symbol_change_num(s,49); return true;
        case 'R': symbol_change_num(s,77); return true;
      }
      break;

    case 'K':
      switch (snd_ch) {
        case 0:   symbol_change_num(s,NIT);return true;
        case 'R': symbol_change_num(s,36); return true;
        case 'A': symbol_change_num(s,19); return true;
      }
      break;

    case 'L':
      switch (snd_ch) {
        case 'A': symbol_change_num(s,57); return true;
        case 'I': symbol_change_num(s,3); return true;
        case 'R': symbol_change_num(s,103); return true;
        case 'U': symbol_change_num(s,71); return true;
        case 'V': symbol_change_num(s,116); return true;
      }
      break;

    case 'M':
      switch (snd_ch) {
        case 0: symbol_change_num(s,NIT); 
        case 'C': symbol_change_num(s,115); return true;
        case 'D': symbol_change_num(s,101); return true;
        case 'G': symbol_change_num(s,12); return true;
        case 'N': symbol_change_num(s,25); return true;
        case 'O': symbol_change_num(s,42); return true;
        case 'T': symbol_change_num(s,109); return true;
      }
      break;

    case 'N':
      switch (snd_ch) {
        case 0: symbol_change_num(s,NIT); return true;
        case 'A': symbol_change_num(s,11); return true;
        case 'B': symbol_change_num(s,41); return true;
        case 'D': symbol_change_num(s,60); return true;
        case 'E': symbol_change_num(s,10); return true;
        case 'H': symbol_change_num(s,113); return true;
        case 'I': symbol_change_num(s,28); return true;
        case 'O': symbol_change_num(s,102); return true;
        case 'P': symbol_change_num(s,93); return true;
      }
      break; 

    case 'O':
      switch (snd_ch) {
        case 0: symbol_change_num(s,OXY); return true;
        case 'G': symbol_change_num(s,118); return true;
        case 'S': symbol_change_num(s,76); return true;
      }
      break;

    case 'P':
      switch (snd_ch) {
        case 0: symbol_change_num(s,PHO); return true;
        case 'A': symbol_change_num(s,91); return true;
        case 'B': symbol_change_num(s,82); return true;
        case 'D': symbol_change_num(s,46); return true;
        case 'M': symbol_change_num(s,61); return true;
        case 'O': symbol_change_num(s,84); return true;
        case 'R': symbol_change_num(s,59); return true;
        case 'T': symbol_change_num(s,78); return true;
        case 'U': symbol_change_num(s,94); return true;
      }
      break;
    
    case 'Q':
      switch (snd_ch) {
        case 0: symbol_change_num(s,OXY); return true;
      }
      break; 

    case 'R':
      switch (snd_ch) {
        case 'A': symbol_change_num(s,88); return true;
        case 'B': symbol_change_num(s,37); return true;
        case 'E': symbol_change_num(s,75); return true;
        case 'F': symbol_change_num(s,104); return true;
        case 'G': symbol_change_num(s,111); return true;
        case 'H': symbol_change_num(s,45); return true;
        case 'N': symbol_change_num(s,86); return true;
        case 'U': symbol_change_num(s,44); return true;
      }
      break;

    case 'S':
      switch (snd_ch) {
        case 0:  symbol_change_num(s,SUL); return true;
        case 'B': symbol_change_num(s,51); return true;
        case 'C': symbol_change_num(s,21); return true;
        case 'E': symbol_change_num(s,34); return true;
        case 'G': symbol_change_num(s,106); return true;
        case 'I': symbol_change_num(s,14); return true;
        case 'M': symbol_change_num(s,62); return true;
        case 'N': symbol_change_num(s,50); return true;
        case 'R': symbol_change_num(s,38); return true;
      }
      break;

    case 'T':
      switch (snd_ch) {
        case 'A': symbol_change_num(s,73); return true;
        case 'B': symbol_change_num(s,65); return true;
        case 'C': symbol_change_num(s,43); return true;
        case 'E': symbol_change_num(s,52); return true;
        case 'H': symbol_change_num(s,90); return true;
        case 'I': symbol_change_num(s,22); return true;
        case 'L': symbol_change_num(s,81); return true;
        case 'M': symbol_change_num(s,69); return true;
        case 'S': symbol_change_num(s,117); return true;
      }
      break;

    case 'U':
      switch (snd_ch) {
        case 'R': symbol_change_num(s,92); return true;
      }
      break;

    case 'V':
      switch (snd_ch) {
        case 'A': symbol_change_num(s,23); return true;
      }
      break;
  
    case 'W':
      switch (snd_ch) {
        case 'T': symbol_change_num(s,74); return true;
      }
      break; 

    case 'X':
      switch (snd_ch) {
        case 'E': symbol_change_num(s,54); return true;
      }
      break;

    case 'Y':
      switch (snd_ch) {
        case 'T': symbol_change_num(s,39); return true;
        case 'B': symbol_change_num(s,70); return true;
      }
      break;

    case 'Z':
      switch (snd_ch) {
        case 'N': symbol_change_num(s,30); return true;
        case 'R': symbol_change_num(s,30); return true;
      }
      break;
  }
  return false; 
}


static void 
symbol_safeset_hydrogens(symbol_t *s, char n)
{
  if (n < 0)
    n = 0; 
  symbol_set_hydrogens(s, n);
}


static void 
graph_sanitize(graph_t *g)
{
#ifdef USING_OPENBABEL
  FOR_ATOMS_OF_MOL(a,g) {
    unsigned char arom = a->IsAromatic(); 
    unsigned char valence = a->GetExplicitValence();
    unsigned char modifier = valence + arom; 
    if (a->GetFormalCharge() == 0 && 
        a->GetImplicitHCount() == 0) {
      switch (a->GetAtomicNum()) {
        case 6: symbol_safeset_hydrogens(a, 4 - modifier); break;
        case 7: symbol_safeset_hydrogens(a, 3 - modifier); break;
        case 8: symbol_safeset_hydrogens(a, 2 - modifier); break;
      }
    }
  }; 
#elif defined USING_RDKIT

#endif
}


// create a bond to a dummy atom type. allows bond modification without state hold 
static edge_t* 
edge_create(graph_t *g, symbol_t *curr, symbol_t *prev)
{
  edge_t *bptr = (edge_t*)0;
#ifdef USING_OPENBABEL
  if (!g->AddBond(curr->GetIdx(), prev->GetIdx(), 1)) {
    fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n", curr->GetIdx(),prev->GetIdx());
    return bptr;
  }
  else 
    bptr = g->GetBond(g->NumBonds() - 1);
#elif defined USING_RDKIT

#endif
  return bptr;
}


static bool 
add_oxy(graph_t *mol, symbol_t *atom)
{
  symbol_t *oxygen = symbol_create(mol, OXY); 
  edge_t *bptr = edge_create(mol, atom, oxygen);
  if (!bptr)  
    return false;
  return true; 
}


typedef struct  {
  uint8_t r_loc; // use to calculate size + add in off-branch positions
  uint8_t r_size;
  uint8_t arom; 
} r_assignment; 

#if 0
static ring_t* rt_alloc(graph_t *g, const size_t size, symbol_t *inc, const uint16_t inc_pos) 
{
  ring_t *ring = 0; 
  ring = (ring_t*)malloc(sizeof(ring_t) + sizeof(locant_t)*(size-1)); 
  memset(ring, 0, sizeof(ring_t) + sizeof(locant_t)*(size-1)); 
  
  symbol_t *c = 0;  
  symbol_t *p = 0;
  edge_t *e   = 0; 

  for (uint16_t i=0; i<size; i++) {
    if (inc && i==inc_pos)
      c = inc; 
    else 
      c = symbol_create(g, CAR, 4); 
    
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
      c->s = symbol_create(g, CAR, 4); 

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

#if 0
/* placeholder for charged alterations */ 
static bool 
add_tauto_dioxy(graph_t *g, symbol_t *p)
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
parse_cyclic(const char *&ptr, 
             uint16_t start, 
             uint16_t end, 
             symbol_t *head, 
             uint16_t head_loc, 
             graph_t *g) 
{
  symbol_t  *c    = 0; 
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


/*
 * -- Parse WLN Notation --
 */
bool 
ReadWLN(const char *wln, graph_t *molecule)
{
  edge_t *curr_edge=0; 
  symbol_t *curr_symbol=0;
  symbol_t *prev_symbol=0;
  ring_t *curr_ring=0;

  uint8_t ring_chars = 0; 
  uint16_t locant_ch = 0; 
  uint16_t digit_n   = 0; 
  int8_t   h_modifier=0; // can be negative to add hydrogens
  
  uint8_t state = 0; // bit field: 
                     // [][dioxo][charge][ring locant_t][ring skip][dash][digit][space]
  
  uint8_t dash_ptr = 0; 
  unsigned char dash_chars[3] = {0}; // last byte is mainly for overflow 
 
  uint8_t stack_ptr = 0;  // WLN branch and ring dependency stack
  struct wlnrefaddr {
    void *addr; 
    char ref; // -1 for (ring_t*) else (symbol_t*) 
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
  curr_symbol = symbol_create(molecule, DUM); 
  
  unsigned char ch; 
  unsigned char ch_nxt; 
  while (*wln) {
    ch     = *wln++;
    ch_nxt = *wln; 
    
    if (state == INIT_READ) switch (ch) {
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
        digit_n *= 10; 
        digit_n += ch - '0'; 

        if (ch_nxt < '0' || ch_nxt > '9') {
          symbol_change_num(curr_symbol, CAR);
          h_modifier = 0; 
          for (uint16_t i=0; i<digit_n-1; i++) {
            prev_symbol = curr_symbol; 
            curr_symbol = symbol_create(molecule, CAR);
            edge_create(molecule, curr_symbol, prev_symbol);  
          }
          digit_n = 0; 
          state &= ~ALKYL_READ; 

          prev_symbol = curr_symbol; 
          curr_symbol = symbol_create(molecule, DUM); 
          curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        }
        else
          state |= ALKYL_READ;
        break;
      
      case 'A':
      case 'J':
        return error("Error: non-atomic symbol used in chain"); 
      
      case 'B':
        symbol_change_num(curr_symbol, BOR);
        
        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 2;
        stack_ptr++; 
        
        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
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
        symbol_change_num(curr_symbol, BRO);
        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          // the molecule has to end here. open charge notation
          
          if (!prev_symbol) {
            prev_symbol = curr_symbol;
            curr_symbol = symbol_create(molecule, DUM); 
            curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          }
          else if (ch_nxt != ' ' || ch_nxt != '\0') 
            return error("Error: terminating symbol closes molecule"); 
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
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          curr_symbol = symbol_create(molecule, DUM); 
          curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 

          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;

      case 'F':
        symbol_change_num(curr_symbol, BRO);
        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          // the molecule has to end here. open charge notation
          if (!prev_symbol) {
            prev_symbol = curr_symbol;
            curr_symbol = symbol_create(molecule, DUM); 
            curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          }
          else if (ch_nxt != ' ' || ch_nxt != '\0') 
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
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          curr_symbol = symbol_create(molecule, DUM); 
          curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 

          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;

      case 'G':
        symbol_change_num(curr_symbol, CHL);
        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          // the molecule has to end here. open charge notation
          if (!prev_symbol) {
            prev_symbol = curr_symbol;
            curr_symbol = symbol_create(molecule, DUM); 
            curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          }
          else if (ch_nxt != ' ' || ch_nxt != '\0') 
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
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          curr_symbol = symbol_create(molecule, DUM); 
          curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 

          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;

      case 'H':
        symbol_incr_hydrogens(prev_symbol); 
        break;

      case 'I':
        symbol_change_num(curr_symbol, IOD);
        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          // the molecule has to end here. open charge notation
          if (!prev_symbol) {
            prev_symbol = curr_symbol;
            curr_symbol = symbol_create(molecule, DUM); 
            curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          }
          else if (ch_nxt != ' ' || ch_nxt != '\0') 
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
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          curr_symbol = symbol_create(molecule, DUM); 
          curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          
          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;
     
      /* [N+](R)(R)(R)(R) */
      case 'K':
        symbol_change_num(curr_symbol, NIT);
        
        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 3;
        stack_ptr++; 
        
        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        break;

      case 'L':
      case 'T':
        state |= RING_READ; 
        ring_chars++; 
        break; 
      
      /* NH(R)(R) */
      case 'M':
        symbol_change_num(curr_symbol, NIT);
        prev_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, prev_symbol, curr_symbol); 
        curr_symbol = prev_symbol;
        break;
     
      /* NR(R)(R) */
      case 'N':
        symbol_change_num(curr_symbol, NIT);
        
        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 2;
        stack_ptr++; 

        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        break;
      
      /* OR(R) */
      case 'O':
        symbol_change_num(curr_symbol, OXY);
        h_modifier=0; 

        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
        break;
      
      case 'P':
        symbol_change_num(curr_symbol, PHO);
        h_modifier=0;
        
        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 3;
        stack_ptr++; 
        
        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        break;

      case 'Q':
        symbol_change_num(curr_symbol, OXY);
        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          // the molecule has to end here. open charge notation
          if (!prev_symbol) {
            prev_symbol = curr_symbol;
            curr_symbol = symbol_create(molecule, DUM); 
            curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          }
          else if (ch_nxt != ' ' && ch_nxt == '\0') 
            return error("Error: terminater used to close molecule"); 
        }
        else if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          curr_symbol = symbol_create(molecule, DUM); 
          curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          
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
        symbol_change_num(curr_symbol, SUL);
        
        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 3;
        stack_ptr++; 

        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        break;
      
      case 'U':
        if (!curr_edge)
          return error("Error: failed to unsaturate bond\n"); 
        edge_unsaturate(curr_edge); 
        break; 

      case 'V':
        symbol_change_num(curr_symbol, CAR);
        if (!add_oxy(molecule, curr_symbol))
          return error("Error: failed to add =O group\n"); 
      
        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        if (!curr_edge)
          return error("Error: failed to add =O group\n"); 
        break;
      
      // if not previous, create dummy carbon
      // really trying to avoid a bit state on this
      case 'W':
        return error("Error: W group needs supporting"); 
        break; 

      case 'X':
        symbol_change_num(curr_symbol, CAR);
        
        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 3;
        stack_ptr++; 

        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        break;

      case 'Y':
        symbol_change_num(curr_symbol, CAR);
        
        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 2;
        stack_ptr++; 

        prev_symbol = curr_symbol; 
        curr_symbol = symbol_create(molecule, DUM); 
        curr_edge = edge_create(molecule, curr_symbol, prev_symbol); 
        break;
      
      case 'Z':
        symbol_change_num(curr_symbol, NIT);
        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          // the molecule has to end here. open charge notation
          if (!prev_symbol) {
            prev_symbol = curr_symbol;
            curr_symbol = symbol_create(molecule, DUM); 
            curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 
          }
          else if (ch_nxt != ' ' || ch_nxt != '\0') 
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
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          curr_symbol = symbol_create(molecule, DUM); 
          curr_edge   = edge_create(molecule, curr_symbol, prev_symbol); 

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
        if (ch_nxt == '&') {
          prev_symbol = NULL;
          stack_ptr   = 0;
          curr_symbol = symbol_create(molecule, DUM);
          state = CLOSE_READ;
          wln++;
          fprintf(stderr, "going to close\n"); 
        }
        else if (curr_ring) {

          state |= SPACE_READ; 
        }
        else return error("Error: space used outside locant|ionic|mixture syntax"); 
        break; 

      case '&':
        if (!stack_ptr) 
          return error("Error: empty dependency stack - too many &?"); 
                
        if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = symbol_create(molecule, DUM); 
          curr_edge = edge_create(molecule, (symbol_t*)dep_stack[stack_ptr-1].addr, prev_symbol); 
          curr_symbol = prev_symbol; 

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
    else if (state == CLOSE_READ) switch (*wln) {
    }
  }
  
  // clean up the hanging bond, one branch in exchange for many
  if (symbol_get_num(curr_symbol) == DUM) 
    graph_delete_symbol(molecule, curr_symbol); 
  
  graph_sanitize(molecule);
  return true; 
}
