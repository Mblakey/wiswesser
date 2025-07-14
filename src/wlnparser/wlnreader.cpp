/*********************************************************************
file: readwln.c 
author : Michael Blakey 2022, updated 2025
description: WLN reader - write out SMILES etc from WLN
***********************************************************************/

#ifdef USING_OPENBABEL
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
  fprintf(stderr,"%s", message); 
  return false;  
}


static void* 
ring_error(const char *message) 
{
  fprintf(stderr,"%s", message); 
  return NULL;  
}

#ifdef USING_OPENBABEL
using namespace OpenBabel;
#define graph_t    OBMol
#define symbol_t   OBAtom
#define edge_t     OBBond

#define symbol_get_num(s)           s->GetAtomicNum()
#define symbol_set_num(s, n)        s->SetAtomicNum(n)
#define symbol_get_charge(s)        s->GetFormalCharge()
#define symbol_set_charge(s, n)     s->SetFormalCharge(n)
#define symbol_set_hydrogens(s, n)  s->SetImplicitHCount(n)
#define symbol_get_hydrogens(s)     s->GetImplicitHCount()
#define symbol_set_aromatic(s, b)   s->SetAromatic(b)
#define symbol_get_aromatic(s)      s->IsAromatic()
#define symbol_incr_hydrogens(s)    s->SetImplicitHCount(s->GetImplicitHCount()+1)
#define symbol_get_valence(s)       s->GetExplicitValence()

#define edge_unsaturate(e, n)       e->SetBondOrder(1 + n)
#define edge_set_order(e, o)        e->SetBondOrder(o);
#define edge_set_aromatic(e, b)     e->SetAromatic(b)

#define graph_new_symbol(g)         g->NewAtom()
#define graph_get_bond(x,y)         g->GetBond(x,y)
#define graph_delete_symbol(g,s)    g->DeleteAtom(s)

#define graph_symbol_iter(s, g)     FOR_ATOMS_OF_MOL(s,g)
#define graph_edge_iter(e, g)       FOR_BONDS_OF_MOL(e,g)

#elif defined USING_RDKIT

#endif


static symbol_t* 
symbol_create(graph_t *mol, uint64_t atomic_num)
{
  symbol_t *atom = graph_new_symbol(mol); 
  symbol_set_num(atom, atomic_num);
  return atom; 
}


static symbol_t*
symbol_create_from_dash(graph_t *mol, unsigned char fst_ch, unsigned char snd_ch) 
{
  switch (fst_ch){
    case 'A':
      switch (snd_ch) {
        case 'C': return symbol_create(mol,89);
        case 'G': return symbol_create(mol,47);
        case 'L': return symbol_create(mol,13);
        case 'M': return symbol_create(mol,95);
        case 'R': return symbol_create(mol,18);
        case 'S': return symbol_create(mol,33);
        case 'T': return symbol_create(mol,85);
        case 'U': return symbol_create(mol,79);
      }
      break; 

    case 'B':
      switch (snd_ch) {
        case 0:   return symbol_create(mol,BOR);
        case 'A': return symbol_create(mol,56);
        case 'E': return symbol_create(mol,4);
        case 'H': return symbol_create(mol,107);
        case 'I': return symbol_create(mol,83);
        case 'K': return symbol_create(mol,97);
        case 'R': return symbol_create(mol,BRO);
      }
      break; 

    case 'C':
      switch (snd_ch) {
        case 0:   return symbol_create(mol,CAR);
        case 'A': return symbol_create(mol,20);
        case 'D': return symbol_create(mol,48);
        case 'E': return symbol_create(mol,58);
        case 'F': return symbol_create(mol,98);
        case 'M': return symbol_create(mol,96);
        case 'N': return symbol_create(mol,112);
        case 'O': return symbol_create(mol,27);
        case 'R': return symbol_create(mol,24);
        case 'S': return symbol_create(mol,55);
        case 'U': return symbol_create(mol,29);
      }
      break; 

    case 'D':
      switch (snd_ch) {
        case 'B': return symbol_create(mol,105);
        case 'S': return symbol_create(mol,110);
        case 'Y': return symbol_create(mol,66);
      }
      break;

    case 'E':
      switch (snd_ch) {
        case 0:   return symbol_create(mol,35);
        case 'R': return symbol_create(mol,68);
        case 'S': return symbol_create(mol,99);
        case 'U': return symbol_create(mol,63);
      }

    case 'F':
      switch (snd_ch) {
        case 0:   return symbol_create(mol,FLU);
        case 'E': return symbol_create(mol,26);
        case 'L': return symbol_create(mol,114);
        case 'M': return symbol_create(mol,100);
        case 'R': return symbol_create(mol,87);
      }
      break;

    case 'G':
      switch (snd_ch) {
        case 0:   return symbol_create(mol,CHL);
        case 'A': return symbol_create(mol,31);
        case 'D': return symbol_create(mol,64);
        case 'E': return symbol_create(mol,32);
      }
      break;

    case 'H':
      switch (snd_ch) {
        case 'E': return symbol_create(mol,2);
        case 'F': return symbol_create(mol,72);
        case 'G': return symbol_create(mol,80);
        case 'O': return symbol_create(mol,67);
        case 'S': return symbol_create(mol,108);
      }
      break;

    case 'I':
      switch (snd_ch) {
        case 0:  return symbol_create(mol,IOD);
        case 'N': return symbol_create(mol,49);
        case 'R': return symbol_create(mol,77);
      }
      break;

    case 'K':
      switch (snd_ch) {
        case 0:   return symbol_create(mol,NIT);
        case 'R': return symbol_create(mol,36);
        case 'A': return symbol_create(mol,19);
      }
      break;

    case 'L':
      switch (snd_ch) {
        case 'A': return symbol_create(mol,57);
        case 'I': return symbol_create(mol,3);
        case 'R': return symbol_create(mol,103);
        case 'U': return symbol_create(mol,71);
        case 'V': return symbol_create(mol,116);
      }
      break;

    case 'M':
      switch (snd_ch) {
        case 0: return symbol_create(mol,NIT);
        case 'C': return symbol_create(mol,115);
        case 'D': return symbol_create(mol,101);
        case 'G': return symbol_create(mol,12);
        case 'N': return symbol_create(mol,25);
        case 'O': return symbol_create(mol,42);
        case 'T': return symbol_create(mol,109);
      }
      break;

    case 'N':
      switch (snd_ch) {
        case 0: return symbol_create(mol,NIT);
        case 'A': return symbol_create(mol,11);
        case 'B': return symbol_create(mol,41);
        case 'D': return symbol_create(mol,60);
        case 'E': return symbol_create(mol,10);
        case 'H': return symbol_create(mol,113);
        case 'I': return symbol_create(mol,28);
        case 'O': return symbol_create(mol,102);
        case 'P': return symbol_create(mol,93);
      }
      break; 

    case 'O':
      switch (snd_ch) {
        case 0: return symbol_create(mol,OXY);
        case 'G': return symbol_create(mol,118);
        case 'S': return symbol_create(mol,76);
      }
      break;

    case 'P':
      switch (snd_ch) {
        case 0: return symbol_create(mol,PHO);
        case 'A': return symbol_create(mol,91);
        case 'B': return symbol_create(mol,82);
        case 'D': return symbol_create(mol,46);
        case 'M': return symbol_create(mol,61);
        case 'O': return symbol_create(mol,84);
        case 'R': return symbol_create(mol,59);
        case 'T': return symbol_create(mol,78);
        case 'U': return symbol_create(mol,94);
      }
      break;
    
    case 'Q':
      switch (snd_ch) {
        case 0: return symbol_create(mol,OXY);
      }
      break; 

    case 'R':
      switch (snd_ch) {
        case 'A': return symbol_create(mol,88);
        case 'B': return symbol_create(mol,37);
        case 'E': return symbol_create(mol,75);
        case 'F': return symbol_create(mol,104);
        case 'G': return symbol_create(mol,111);
        case 'H': return symbol_create(mol,45);
        case 'N': return symbol_create(mol,86);
        case 'U': return symbol_create(mol,44);
      }
      break;

    case 'S':
      switch (snd_ch) {
        case 0:  return symbol_create(mol,SUL);
        case 'B': return symbol_create(mol,51);
        case 'C': return symbol_create(mol,21);
        case 'E': return symbol_create(mol,34);
        case 'G': return symbol_create(mol,106);
        case 'I': return symbol_create(mol,14);
        case 'M': return symbol_create(mol,62);
        case 'N': return symbol_create(mol,50);
        case 'R': return symbol_create(mol,38);
      }
      break;

    case 'T':
      switch (snd_ch) {
        case 'A': return symbol_create(mol,73);
        case 'B': return symbol_create(mol,65);
        case 'C': return symbol_create(mol,43);
        case 'E': return symbol_create(mol,52);
        case 'H': return symbol_create(mol,90);
        case 'I': return symbol_create(mol,22);
        case 'L': return symbol_create(mol,81);
        case 'M': return symbol_create(mol,69);
        case 'S': return symbol_create(mol,117);
      }
      break;

    case 'U':
      switch (snd_ch) {
        case 'R': return symbol_create(mol,92);
      }
      break;

    case 'V':
      switch (snd_ch) {
        case 'A': return symbol_create(mol,23);
      }
      break;
  
    case 'W':
      switch (snd_ch) {
        case 'T': return symbol_create(mol,74);
      }
      break; 

    case 'X':
      switch (snd_ch) {
        case 'E': return symbol_create(mol,54);
      }
      break;

    case 'Y':
      switch (snd_ch) {
        case 'T': return symbol_create(mol,39);
        case 'B': return symbol_create(mol,70);
      }
      break;

    case 'Z':
      switch (snd_ch) {
        case 'N': return symbol_create(mol,30);
        case 'R': return symbol_create(mol,30);
      }
      break;
  }
  return NULL; 
}


static void 
symbol_safeset_hydrogens(symbol_t *s, char n)
{
  if (n < 0)
    n = 0; 
  symbol_set_hydrogens(s, n);
}


static void 
graph_cleanup_hydrogens(graph_t *g)
{
  graph_symbol_iter(a, g) {
    unsigned char arom     = symbol_get_aromatic(a); 
    unsigned char valence  = symbol_get_valence(a); 
    unsigned char modifier = valence + arom;
    if (symbol_get_charge(a) == 0 && symbol_get_hydrogens(a) == 0) {
      switch (a->GetAtomicNum()) {
        case 6: symbol_safeset_hydrogens(a, 4 - modifier); break;
        case 7: symbol_safeset_hydrogens(a, 3 - modifier); break;
        case 8: symbol_safeset_hydrogens(a, 2 - modifier); break;
      }
    }
  }; 
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
  edge_unsaturate(bptr, 1);
  return true; 
}


static void 
add_methyl(graph_t *mol, symbol_t *atom) {
  symbol_t *methyl = symbol_create(mol, CAR);
  edge_t   *e = edge_create(mol, methyl, atom);
}


struct wlnlocant {
  symbol_t *s; 
  struct wlnlocant *lft; 
  struct wlnlocant *rgt; 
  bool    is_bridge;
};

struct wlnsubcycle {
  unsigned int size;
  unsigned int locant;
  bool aromatic;
}; 

struct wlnpath {
  uint8_t  size; 
  wlnlocant path[1]; // malloc sizeof(locant_t) * (size-1) + (1 byte for size)
}; 


static struct wlnpath* 
wlnpath_alloc(graph_t *g, unsigned int size) 
{
  struct wlnpath *ring = 0; 
  ring = (struct wlnpath*)malloc(sizeof(struct wlnpath) + sizeof(struct wlnlocant)*(size-1)); 
  memset(ring, 0, sizeof(struct wlnpath) + sizeof(struct wlnlocant)*(size-1)); 
   
  ring->path[0].s = symbol_create(g, CAR);
  for (uint16_t i = 1; i < size; i++) {
    struct wlnlocant *curr = &ring->path[i];
    struct wlnlocant *prev = &ring->path[i-1];

    curr->s = symbol_create(g, CAR);
    if (!edge_create(g, curr->s, prev->s)) {
      free(ring);
      return NULL;
    }
  }
  ring->size = size; 
  return ring; 
}


static bool
pathsolverIII_fast(graph_t *g, 
                   const struct wlnpath *r, 
                   struct wlnsubcycle SSSR[], 
                   uint8_t SSSR_ptr) 
{
  uint8_t last_idx = r->size-1; 
  struct pathmapping {
    uint8_t nxt_locant; // used for pathsolver 
    uint8_t nlocants;
  } *path = (struct pathmapping*)alloca(r->size*sizeof(struct pathmapping));

  for (uint16_t i = 1; i < r->size; i++) {
    path[i].nlocants = 2;
    path[i-1].nxt_locant = i;
  }
  path[0].nlocants = 3; 
  path[last_idx].nlocants = 3; 
  path[last_idx].nxt_locant = last_idx;

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

  uint8_t steps, start, end;
  bool arom;
  struct wlnsubcycle *subcycle; 
  
  for (uint16_t i=0; i<SSSR_ptr; i++) {
    subcycle = &SSSR[i];   
    steps    = subcycle->size; 
    start    = subcycle->locant; 
    arom     = subcycle->aromatic;
    end      = start; 
    
    // if used max times in ring, shift along path
    while (path[start].nlocants == 0 && start < r->size) {
      start++;
      steps--; 
    }

    for (uint16_t s = 0; s < steps-1; s++) {
      unsigned char nxt = path[end].nxt_locant; 
      symbol_set_aromatic(r->path[end].s, arom);
      edge_t *e = graph_get_bond(r->path[end].s, r->path[nxt].s);
      edge_set_aromatic(e, arom);
      end = nxt; 
    }
    symbol_set_aromatic(r->path[end].s, arom);

    fprintf(stderr,"%d: %c --> %c (%d)\n",steps, start + 'A', end + 'A', subcycle->aromatic); 
    
    path[start].nlocants--;
    path[end].nlocants--;
    path[start].nxt_locant = end;

    edge_t *cross = edge_create(g, r->path[start].s, r->path[end].s);
    edge_set_aromatic(cross, arom);
  }
  return r;   
}

#if 0
static struct ring_t* 
pathsolverIII(graph_t *g, struct ring_t *r, 
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


/* keep the struct scoped, no leaking */
static void depstack_cleanup(graph_t *mol, const void *_depstack, unsigned int stack_ptr)
{
  struct wlnrefaddr {
    void *addr; 
    char ref; 
  } *dep_stack = (struct wlnrefaddr*)_depstack;

  // clean up the dep stack
  for (unsigned int i = 0; i < stack_ptr; i++) {
    symbol_t *potential_symbol = (symbol_t*)dep_stack[i].addr;
    if (dep_stack[i].ref == -1)
      free((struct wlnpath*)dep_stack[i].addr);
#if 0
    else if (symbol_get_num(potential_symbol) == 6 || 
             (symbol_get_num(potential_symbol) == 7 && 
              symbol_get_charge(potential_symbol) == +1)) 
    {
      for (unsigned int j=0;j<dep_stack[i].ref;j++) 
        add_methyl(mol, potential_symbol);  
    }
#endif
  }
}


static symbol_t*
parse_opening_terminator(graph_t *mol, unsigned char ch) {
  switch (ch) {
    case 'E': return symbol_create(mol, BRO);
    case 'F': return symbol_create(mol, FLU);
    case 'G': return symbol_create(mol, CHL);
    case 'I': return symbol_create(mol, IOD);
    case 'Q': return symbol_create(mol, OXY);
    case 'Z': return symbol_create(mol, NIT);
  } 
  return NULL;
}

// ring_parse states 
#define SSSR_READ   0x01
#define MULTI_READ  0x02
#define PSEUDO_READ 0x04


static bool 
parse_aromaticity(const char **wln, 
                  wlnsubcycle SSSR[], 
                  unsigned short SSSR_ptr)
{
  const char *ptr = *wln;
  unsigned char ch = *ptr; 

  unsigned short assigments = 0;
  while (ch) {
    if (*ptr == 'J') {
      if (assigments == 1 && !SSSR[0].aromatic) {
        for (unsigned int i=1; i < SSSR_ptr; i++)
          SSSR[i].aromatic = false;
      }
      else if (assigments != SSSR_ptr) {
        fprintf(stderr, "Error: not enough aromaticity assignments for wln ring\n");
        return false;
      }
      *wln = ptr;
      return true;
    }

    ch = *ptr++; 
    if (ch == 'T' || ch == '&') {
      if (assigments == SSSR_ptr) {
        fprintf(stderr, "Error: too many aromaticity assignments for wln ring\n");
        return false; 
      }
      SSSR[assigments++].aromatic = (ch == '&');
    }
    else {
      fprintf(stderr, "Error: invalid character in aromaticity parse\n - %c", ch);
      return false;
    }
  }
  return false;
}


static struct wlnpath* 
parse_cyclic(const char **wln, graph_t *g) 
{
  const char *ptr = *wln; // set at the end
  struct wlnpath   *ring = 0; 
  
  bool      expecting_locant = false;
  uint8_t   locant_ch   = 0; 
  uint8_t   max_path_size  = 0;  
  
  uint8_t SSSR_ptr  = 0; 
  struct wlnsubcycle SSSR[32]; // this really is sensible for WLN

  unsigned char ch; 
  while (*ptr) {
    ch = *ptr++; 
    if (ch >= '0' && ch <= '9') {

      if (expecting_locant) 
        ; // MOVE TO MULTI BLOCK 
          //
      else {
        max_path_size += ch - '0' - (2*(max_path_size>0)); 
        SSSR[SSSR_ptr].size = ch - '0'; 
        SSSR[SSSR_ptr].aromatic = true; // default is aromatic 
        SSSR[SSSR_ptr++].locant  =  locant_ch; 
        locant_ch = 0; 
        while ((ch = *ptr) && ch >= '0' && ch <= '9') {
          max_path_size += ch - '0' - (2*(max_path_size>0)); 
          SSSR[SSSR_ptr].size = ch - '0'; 
          SSSR[SSSR_ptr].aromatic = true; // default is aromatic 
          SSSR[SSSR_ptr++].locant  =  locant_ch; 
          locant_ch = 0; 
          ptr++;
        }
      }
    }
    else if (expecting_locant) {
      if (ch >= 'A' && ch <= 'Z') {
        locant_ch = ch - 'A';
        expecting_locant = false;
      }
      else return (struct wlnpath*)ring_error("Error: ring locant requires a space before the character\n"); 
    }
    else switch (ch) {
      case 'J': 
        *wln = ptr; 
        ring = wlnpath_alloc(g, max_path_size); 
        if (!pathsolverIII_fast(g, ring, SSSR, SSSR_ptr)) {
          free(ring);
          return (struct wlnpath*)ring_error("Error: failed on path solver algorithm");  
        }
        else return ring; 

      case 'T':
      case '&':
        if (!SSSR_ptr) 
          return (struct wlnpath*)ring_error("Error: no rings designated before ring close\n"); 
        if (!parse_aromaticity(&(--ptr), SSSR, SSSR_ptr)) {
          free(ring);
          return NULL;
        }
        break;

      case ' ':
        if (expecting_locant)
          return (struct wlnpath*)ring_error("Error: double space seen in SSSR ring block\n"); 
        expecting_locant = true; 
        break;
      default: return (struct wlnpath*)ring_error("Error: invalid character in SSSR ring block\n"); 
    }
  }

  *wln = ptr;
  return ring; 
}


/*
 * -- Parse WLN Notation --
 */
bool 
ReadWLN(const char *wln, graph_t *mol)
{
#ifdef USING_OPENBABEL
  mol->BeginModify(); 
  mol->SetAromaticPerceived(true);
  mol->SetChiralityPerceived(true); // no stereo for WLN
#endif

  edge_t *curr_edge=0; 
  symbol_t *init_symbol=0;
  symbol_t *curr_symbol=0;
  symbol_t *prev_symbol=0;
  struct wlnpath *curr_ring=0;

  uint8_t ring_chars = 0; 
  uint16_t locant_ch = 0; 
  uint16_t counter   = 0; 
  uint8_t  unsaturation=0; 
  
  uint8_t state = 0; // bit field: 
                     // [][dioxo][charge][ring locant_t][ring skip][dash][digit][space]
  
  uint8_t dash_ptr = 0; 
  unsigned char dash_chars[3] = {0}; // last byte is mainly for overflow 
 
  uint8_t stack_ptr = 0;  // WLN branch and ring dependency stack
  struct wlnrefaddr {
    void *addr; 
    char ref; // -1 for (struct ring_t*) else (symbol_t*) how many extra branches >2 can the symbol have
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
  // idea entirely *--> grow..., delete at the end to save branches
  init_symbol = prev_symbol = symbol_create(mol, DUM); 
  curr_symbol = parse_opening_terminator(mol, *wln);
  if (curr_symbol) {
    curr_edge   = edge_create(mol, curr_symbol, prev_symbol);
    prev_symbol = curr_symbol; 
    wln++;
  }
  
  unsigned char ch = *wln; 
  while (*wln) {
    ch  = *wln++;
    switch (ch) {
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
        counter = ch - '0';
        while ((ch = *wln) && ch >= '0' && ch <= '9') {
          counter *= 10; 
          counter += ch - '0'; 
          wln++;
        }
        for (uint16_t i=0; i<counter; i++) {
          curr_symbol = symbol_create(mol, CAR);
          curr_edge = edge_create(mol, curr_symbol, prev_symbol);  
          edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
          prev_symbol = curr_symbol; 
        }
        break;
      
      case 'A':
      case 'J':
        return error("Error: non-atomic symbol used in chain"); 
      
      case 'B':
        curr_symbol = symbol_create(mol, BOR);

        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 1;
        stack_ptr++; 
        
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol);
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        prev_symbol = curr_symbol; 
        break; 
      
      case 'C':
        return error("Error: WLN symbol C currently unhandled\n"); 

      // extrememly rare open chelate notation
      case 'D':
        return error("Error: WLN symbol D (chelate) currently unhandled\n"); 
      
      // terminator symbol - no bond movement
      case 'E':
        curr_symbol = symbol_create(mol, BRO);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        if (unsaturation) return error("Error: unsaturation on a terminator is not allowed\n");

        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          if (*wln != ' ' && *wln != '\0')
            return error("Error: terminator character closes molecule\n");
        }
        else if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;

      case 'F':
        curr_symbol = symbol_create(mol, FLU);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        if (unsaturation) return error("Error: unsaturation on a terminator is not allowed\n");

        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          if (*wln != ' ' && *wln != '\0')
            return error("Error: terminator character closes molecule\n");
        }
        else if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;

      case 'G':
        curr_symbol = symbol_create(mol, CHL);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        if (unsaturation) return error("Error: unsaturation on a terminator is not allowed\n");

        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          if (*wln != ' ' && *wln != '\0')
            return error("Error: terminator character closes molecule\n");
        }
        else if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;

      case 'H':
        symbol_incr_hydrogens(prev_symbol); 
        break;

      case 'I':
        curr_symbol = symbol_create(mol, IOD);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        if (unsaturation) return error("Error: unsaturation on a terminator is not allowed\n");

        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          if (*wln != ' ' && *wln != '\0')
            return error("Error: terminator character closes molecule\n");
        }
        else if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;
     
      /* [N+](R)(R)(R)(R) */
      case 'K':
        curr_symbol = symbol_create(mol, NIT);
        symbol_set_charge(curr_symbol, +1); 
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        
        counter = 0;
        if (*wln == '&') while (*wln == '&' && counter < 2) {
          add_methyl(mol, curr_symbol); 
          counter++; 
          wln++;
        }
        if (counter < 2) {
          dep_stack[stack_ptr].addr = curr_symbol; 
          dep_stack[stack_ptr].ref  = 2 - counter;
          stack_ptr++; 
        }
       
        // final contraction and stack movement if possible
        if (*wln == '&') {
          add_methyl(mol, curr_symbol);
          wln++;
          if (!stack_ptr) {
            if (*wln != ' ' && *wln != '\0')
              return error("Error: terminator character closes molecule\n");
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
        }
        else 
          prev_symbol = curr_symbol; 
        break;

      case 'L':
      case 'T':
        curr_ring = parse_cyclic(&wln, mol);
        if (!curr_ring) return false;
        break; 
      
      /* NH(R)(R) */
      case 'M':
        curr_symbol = symbol_create(mol, NIT); 
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        prev_symbol = curr_symbol; 
        break;
     
      /* NR(R)(R) */
      case 'N':
        curr_symbol = symbol_create(mol, NIT);

        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 1;
        stack_ptr++; 
        
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        prev_symbol = curr_symbol; 
        break;
      
      /* OR(R) */
      case 'O':
        curr_symbol = symbol_create(mol, OXY); 
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        prev_symbol = curr_symbol; 
        break;
      
      case 'P':
        curr_symbol = symbol_create(mol, PHO);

        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 1;
        stack_ptr++; 
        
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        prev_symbol = curr_symbol; 
        break;

      case 'Q':
        curr_symbol = symbol_create(mol, OXY);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        if (unsaturation) return error("Error: unsaturation on a terminator is not allowed\n");

        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          if (*wln != ' ' && *wln != '\0')
            return error("Error: terminator character closes molecule\n");
        }
        else if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
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
        curr_symbol = symbol_create(mol, SUL);

        dep_stack[stack_ptr].addr = curr_symbol; 
        dep_stack[stack_ptr].ref  = 3;
        stack_ptr++; 
        
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        prev_symbol = curr_symbol; 
        break;
      
      case 'U':
        unsaturation++; 
        break; 

      case 'V':
        curr_symbol = symbol_create(mol, CAR);
        if (!add_oxy(mol, curr_symbol))
          return error("Error: failed to add =O group\n"); 
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        prev_symbol = curr_symbol; 
        break;
      
      // if not previous, create dummy carbon
      // really trying to avoid a bit state on this
      case 'W':
        return error("Error: W group needs supporting"); 
        break; 

      case 'X':
        curr_symbol = symbol_create(mol, CAR);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        
        // default methyl contraction for X,Y,K
        counter = 0;
        if (*wln == '&') while (*wln == '&' && counter < 2) {
          add_methyl(mol, curr_symbol); 
          counter++; 
          wln++;
        }

        if (counter < 2) {
          dep_stack[stack_ptr].addr = curr_symbol; 
          dep_stack[stack_ptr].ref  = 2 - counter;
          stack_ptr++; 
        }

        // final contraction and stack movement if possible
        if (*wln == '&') {
          add_methyl(mol, curr_symbol);
          wln++;
          if (!stack_ptr) {
            if (*wln != ' ' && *wln != '\0')
              return error("Error: terminator character closes molecule\n");
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
        }
        else 
          prev_symbol = curr_symbol; 
        break;

      case 'Y':
        curr_symbol = symbol_create(mol, CAR);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        edge_unsaturate(curr_edge, unsaturation); unsaturation = 0;
        
        // default methyl contraction for X,Y,K
        if (*wln == '&') { 
          add_methyl(mol, curr_symbol); 
          counter++; 
          wln++;
        }
        else {
          dep_stack[stack_ptr].addr = curr_symbol; 
          dep_stack[stack_ptr].ref  = 1;
          stack_ptr++; 
        }

        // final contraction and stack movement if possible
        if (*wln == '&') {
          add_methyl(mol, curr_symbol);
          wln++;
          if (!stack_ptr) {
            if (*wln != ' ' && *wln != '\0')
              return error("Error: terminator character closes molecule\n");
          }
          else if (dep_stack[stack_ptr-1].ref == -1) {
            // ring is the next object to resolve
            fprintf(stderr, "hook in ring code\n"); 
            return false; 
          }
          else {
            // a branch must be open.
            prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
            if (--dep_stack[stack_ptr-1].ref == 0)
              stack_ptr--; 
          }
        }
        else 
          prev_symbol = curr_symbol; 
        break;
      
      case 'Z':
        curr_symbol = symbol_create(mol, NIT);
        curr_edge   = edge_create(mol, curr_symbol, prev_symbol); 
        if (unsaturation) return error("Error: unsaturation on a terminator is not allowed\n");

        // note: terminators can act on an empty stack, '&' cannot
        if (!stack_ptr) {
          if (*wln != ' ' && *wln != '\0')
            return error("Error: terminator character closes molecule\n");
        }
        else if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          // a branch must be open.
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr; 
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
#if 0
        if (ch_nxt == '&') {
          prev_symbol = NULL;
          stack_ptr   = 0;
          curr_symbol = symbol_create(molecule, DUM);
          state = BRANCH_CLOSE;
          wln++;
          fprintf(stderr, "going to close\n"); 
        }
#endif
        if (0);
        else if (curr_ring) {

        }
        else return error("Error: space used outside locant|ionic|mixture syntax\n"); 
        break; 

      case '&':
        if (!stack_ptr) 
          return error("Error: empty dependency stack - too many &?\n"); 
                
        if (dep_stack[stack_ptr-1].ref == -1) {
          // ring is the next object to resolve
          fprintf(stderr, "hook in ring code\n"); 
          return false; 
        }
        else {
          prev_symbol = (symbol_t*)dep_stack[stack_ptr-1].addr;  
          if (--dep_stack[stack_ptr-1].ref == 0)
            stack_ptr--; 
        }
        break;

      case '\n':
        break; 
      case '/':
        return error("Error: slash seen outside of ring - multipliers currently unsupported\n");
      default:
        fprintf(stderr, "Error: invalid character read for WLN notation - %c(%u)\n", ch, ch);
        return false; 
    }
  }

  // clean up the hanging bond, one branch in exchange for many
  graph_delete_symbol(mol, init_symbol); 
  depstack_cleanup(mol, dep_stack, stack_ptr);  
  graph_cleanup_hydrogens(mol);

#ifdef USING_OPENBABEL
  OBKekulize(mol); 
  //mol->EndModify(); // do not retrigger chirality 
#endif
  return true; 
}
