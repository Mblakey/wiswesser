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
#include <stdbool.h>
#include <stdint.h>
#include <stdarg.h>

#include <errno.h>
#include <string.h>

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

#define WLN_ERROR -1
#define WLN_OK 0

char *wln_ptr; 

static int wln_error(const char *format, ...) 
{
  va_list args;
  va_start(args, format);
  fprintf(stderr, "Error: ");  
  vfprintf(stderr, format, args); 
  va_end(args); 
  return WLN_ERROR;  
}



#ifdef USING_OPENBABEL
using namespace OpenBabel;
#define graph_t    OBMol
#define symbol_t   OBAtom
#define edge_t     OBBond
#define ring_t struct wlnpath

#define symbol_get_id(s)            s->GetId()
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
#define symbol_nbor_iter(a, s)      FOR_NBORS_OF_ATOM(a,s)

#define edge_unsaturate(e)          e->SetBondOrder(1+e->GetBondOrder())
#define edge_set_order(e, o)        e->SetBondOrder(o);
#define edge_set_aromatic(e, b)     e->SetAromatic(b)

#define edge_set_begin(e, s)       e->SetBegin(s)
#define edge_set_end(e, s)         e->SetEnd(s)
#define edge_get_end(e)            e->GetEndAtom()

#define graph_new_symbol(g)         g->NewAtom()
#define graph_new_edge(g)           g->NewBond()
#define graph_set_edge(g,e)         g->AddBond(*e)
#define graph_get_edge(g,x,y)       g->GetBond(x,y)
#define graph_delete_symbol(g,s)    g->DeleteAtom(s)
#define graph_delete_edge(g,e)      g->DeleteBond(e)
#define graph_num_atoms(g)          g->NumAtoms()

#define graph_symbol_iter(s, g)     FOR_ATOMS_OF_MOL(s,g)
#define graph_edge_iter(e, g)       FOR_BONDS_OF_MOL(e,g)

#elif defined USING_RDKIT

#endif

static int start_wln_parse(graph_t *mol); 
static int branch_recursive_parse(graph_t *mol, edge_t *edge, ring_t *ring); 
static int parse_ring_locants(graph_t *mol, ring_t *ring); 

static symbol_t* 
symbol_create(graph_t *mol, uint16_t atomic_num)
{
  symbol_t *atom = graph_new_symbol(mol); 
  symbol_set_num(atom, atomic_num);
  return atom; 
}


static edge_t* edge_create(graph_t *mol, symbol_t *parent) {
  edge_t *bond = graph_new_edge(mol); 
  edge_set_begin(bond, parent); 
  edge_set_order(bond, 1);
  return bond; 
}


static bool edge_bond(graph_t *mol, edge_t *bond, symbol_t *child)
{
  edge_set_end(bond, child); 
  graph_set_edge(mol, bond);
  return true; 
}


static symbol_t*
symbol_change(symbol_t *s, uint16_t atomic_num)
{
  symbol_set_num(s, atomic_num);
  return s; 
}


static symbol_t*
dash_symbol_create(graph_t *mol, char **wln) 
{
  char *ptr = *wln;
  unsigned char fst_ch = *ptr++; // this is known to not be null
  unsigned char snd_ch = *ptr++; // this could be null (C-string safe)

  if (snd_ch && *ptr != '-')
    return NULL;
  else 
    *wln = ++ptr;

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
graph_cleanup_hydrogens(graph_t *mol)
{
  graph_symbol_iter(a, mol) {
    symbol_t *s = &(*a); 
    unsigned char arom     = symbol_get_aromatic(s); 
    unsigned char valence  = symbol_get_valence(s); 
    unsigned char modifier = valence + arom;
    if (symbol_get_charge(s) == 0 && symbol_get_hydrogens(s) == 0) {
      switch (a->GetAtomicNum()) {
        case 6: symbol_safeset_hydrogens(s, 4 - modifier); break;
        case 7: symbol_safeset_hydrogens(s, 3 - modifier); break;
        case 8: symbol_safeset_hydrogens(s, 2 - modifier); break;
      }
    }
  }; 
}


// create a bond to a dummy atom type. allows bond modification without state hold 
static edge_t* 
edge_create(graph_t *mol, symbol_t *curr, symbol_t *prev)
{
  edge_t *bptr = (edge_t*)0;
#ifdef USING_OPENBABEL
  if (!mol->AddBond(curr->GetIdx(), prev->GetIdx(), 1)) {
    fprintf(stderr, "Error: failed to make bond betweens atoms %d --> %d\n", curr->GetIdx(),prev->GetIdx());
    return bptr;
  }
  else 
    bptr = mol->GetBond(mol->NumBonds() - 1);
#elif defined USING_RDKIT

#endif
  return bptr;
}


static int add_oxy(graph_t *mol, symbol_t *atom)
{
  symbol_t *oxygen = symbol_create(mol, OXY); 
  edge_t *bptr = edge_create(mol, atom, oxygen);
  if (!bptr) 
    return WLN_ERROR;

  edge_set_order(bptr, 2);
  return WLN_OK; 
}


static bool 
add_dioxo(graph_t *mol, symbol_t *atom) {
  if (!add_oxy(mol, atom))
    return false;

  symbol_t *oxygen = symbol_create(mol, OXY); 
  switch (symbol_get_num(atom)) {
    case NIT:
      symbol_set_charge(oxygen, -1);
      symbol_set_charge(atom, +1);
      break;
  }

  edge_t *bptr = edge_create(mol, atom, oxygen);
  if (!bptr)  
    return false;
  return true;
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


static ring_t* ring_create(graph_t *mol, unsigned int size) 
{
  struct wlnpath *ring = 0; 
  ring = (struct wlnpath*)malloc(sizeof(struct wlnpath) + sizeof(struct wlnlocant)*(size-1)); 
  memset(ring, 0, sizeof(struct wlnpath) + sizeof(struct wlnlocant)*(size-1)); 
   
  ring->path[0].s = symbol_create(mol, CAR);
  for (uint16_t i = 1; i < size; i++) {
    struct wlnlocant *curr = &ring->path[i];
    struct wlnlocant *prev = &ring->path[i-1];

    curr->s = symbol_create(mol, CAR);
    if (!edge_create(mol, curr->s, prev->s)) {
      free(ring);
      return NULL;
    }
  }
  ring->size = size; 
  return ring; 
}


static ring_t* ring_create_benzene(graph_t *mol) 
{
  wlnpath *benzene = ring_create(mol, 6);
  edge_t *bond = edge_create(mol, benzene->path[0].s, benzene->path[5].s);
  edge_set_aromatic(bond, true);

  for (unsigned int i=1; i<6; i++) {
    symbol_t *f = benzene->path[i].s;
    symbol_t *p = benzene->path[i-1].s;

    bond = graph_get_edge(mol, f, p);

    symbol_set_aromatic(f, true);
    symbol_set_aromatic(p, true);
    edge_set_aromatic(bond, true);
  }

  return benzene;
}


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
static void pathsolverIII_fast(graph_t *mol, 
                               ring_t *r, 
                               struct wlnsubcycle *SSSR, 
                               uint8_t nSSSR)
{
  uint8_t last_idx = r->size-1; 
  struct pathmapping {
    uint8_t nxt_locant; // used for pathsolver 
    uint8_t nlocants;
  } *path = (struct pathmapping*)alloca(r->size*sizeof(struct pathmapping));

  for (uint16_t i = 1; i < r->size; i++) {
    path[i].nlocants = 1;
    path[i-1].nxt_locant = i;
  }
  path[0].nlocants = 2; 
  path[last_idx].nlocants = 2; 
  path[last_idx].nxt_locant = last_idx;

  uint8_t steps, start, end;
  struct wlnsubcycle *subcycle; 
  bool arom;
  
  for (uint16_t i=0; i<nSSSR; i++) {
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
      edge_t *e = graph_get_edge(mol, r->path[end].s, r->path[nxt].s);
      edge_set_aromatic(e, arom);
      end = nxt; 
    }
    symbol_set_aromatic(r->path[end].s, arom);

    fprintf(stderr,"%d: %c --> %c (%d)\n",steps, start + 'A', end + 'A', subcycle->aromatic); 
    
    path[start].nlocants--;
    path[start].nxt_locant = end;

    edge_t *cross = edge_create(mol, r->path[start].s, r->path[end].s);
    edge_set_aromatic(cross, arom);
  }
}


static int pathsolver_recursive_floodfill(struct wlnpath *r, 
                                          symbol_t *s, 
                                          bool *seen, 
                                          int *path, 
                                          int *best_path,
                                          int n)
{
  if (n==0) { 
    int end; 
    for (end = 0; end < r->size; end++) {
      unsigned int id = symbol_get_id(s); 
      if (r->path[end].s == s) {
        seen[id] = false;
        break; 
      }
    }
    path[n] = end; 
    if (end > best_path[0])
      memcpy(best_path, path, r->size*sizeof(int));
    return end; 
  }
  
  int max = 0; 
  symbol_nbor_iter(a, s) {
    symbol_t *nbr = &(*a); 
    unsigned int id = symbol_get_id(nbr); 
    if (!seen[id]) {
      for (int i = 0; i < r->size; i++) {
        if (r->path[i].s == s) {
          path[n] = i; 
          break; 
        }
      }
      seen[id] = true; 
      int loc  = pathsolver_recursive_floodfill(r, nbr, seen, path, best_path, n-1);
      if (loc > max)
        max = loc; 
      seen[id] = false; 
    }
  }
  return max;  
}

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
static void pathsolverIII(graph_t *mol, 
                          ring_t *r, 
                          struct wlnsubcycle *SSSR, 
                          uint8_t nSSSR)
{
  const unsigned int natoms = graph_num_atoms(mol); 
  const uint8_t last_idx = r->size-1; 
  
  uint8_t *nlocants = (uint8_t*)alloca(r->size); 

  int *path = (int*)alloca(r->size*sizeof(int)); 
  int *best_path = (int*)malloc(r->size*sizeof(int)); 

  bool *seen = (bool*)malloc(natoms); 
  
  for (uint16_t i = 1; i < r->size; i++) 
    nlocants[i] = 1;
  nlocants[0]        = 2; 
  nlocants[last_idx] = 2; 
  
  uint8_t steps, start, end;
  bool arom;
  struct wlnsubcycle *subcycle; 
  
  for (uint16_t i=0; i<nSSSR; i++) {
    subcycle = &SSSR[i];   
    steps    = subcycle->size; 
    start    = subcycle->locant; 
    arom     = subcycle->aromatic;
    
    memset(seen, 0, natoms); 
    memset(best_path, 0, r->size*sizeof(int)); 
    
    // if used max times in ring, shift along path
    while (nlocants[start] == 0 && start < r->size) {
      start++;
      steps--; 
    }
    
    symbol_t *s_sym = r->path[start].s; 
    seen[symbol_get_id(s_sym)] = true; 
    path[last_idx] = start; 
    end = pathsolver_recursive_floodfill(r, s_sym, seen, path, best_path, steps-1); 
    seen[symbol_get_id(s_sym)] = false; 

    edge_t *e = edge_create(mol, r->path[end].s, r->path[start].s);
    edge_set_aromatic(e, arom);

    for (int i=1; i < subcycle->size; i++) {
      symbol_set_aromatic(r->path[best_path[i]].s, arom);
      symbol_set_aromatic(r->path[best_path[i-1]].s, arom);
      edge_t *e = graph_get_edge(mol, r->path[best_path[i]].s, r->path[best_path[i-1]].s);
      edge_set_aromatic(e, arom);
    }

    fprintf(stderr,"%d: %c --> %c (%d)\n",steps, start + 'A', end + 'A', subcycle->aromatic); 
    
    nlocants[start]--;
    nlocants[end]--;
  }
  
  free(best_path); 
  free(seen); 
}


/* ################# Dispatch Functions ################# */

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


static int parse_aromaticity(wlnsubcycle SSSR[], 
                             unsigned short nSSSR)
{
  unsigned char ch = *wln_ptr; 
  unsigned short assigments = 0;

  while (ch) {
    if (*wln_ptr == 'J') {
      if (assigments == 1 && !SSSR[0].aromatic) {
        for (unsigned int i=1; i < nSSSR; i++)
          SSSR[i].aromatic = false;
      }
      else if (assigments != nSSSR) 
        return wln_error("not enough aromaticity assignments for wln ring\n");
    }

    ch = *wln_ptr++; 
    if (ch == 'T' || ch == '&') {
      if (assigments == nSSSR) {
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


// -1 for error
static int parse_locant() 
{
  int locant_ch   = -1;
  unsigned char ch = *wln_ptr; 

  while (*wln_ptr) {
    if ((ch >= 'A' && ch <= 'Z')) {
      if (locant_ch == -1)
        locant_ch = ch - 'A';
      else break;
    }
    else if ((ch >= '0' && ch <= '9')) {
      if (locant_ch == -1)
        return WLN_ERROR;
      else break;
    }
    else if (ch == '-')
      ;
    else if (ch == '&')
      ;
    else 
      return WLN_ERROR;
    ch = *(++wln_ptr);
  }
  
  return locant_ch;
}



/* this expects the symbol after the opening L|T */
static int cyclic_recursive_parse(graph_t *mol) 
{
  bool seen_pseudo = false; 
  struct wlnpath   *ring = NULL; 
  int locant_ch     = 0; 
  int max_path_size = 0;  
  
  uint8_t nSSSR  = 0; 
  struct wlnsubcycle SSSR[32]; // this really is sensible for WLN
  unsigned char ch;   
ring_parse_sssr:
  while (*wln_ptr) {
    ch = *wln_ptr++; 
    switch (ch) {
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        max_path_size += (ch - '0') - (2 * (max_path_size > 0)); 
        SSSR[nSSSR].size      = ch - '0'; 
        SSSR[nSSSR].aromatic  = true; // default is aromatic 
        SSSR[nSSSR++].locant  =  locant_ch; 
        locant_ch = 0; 
        break; 

      case ' ':
        if (*wln_ptr >= 'A' && *wln_ptr <= 'Z') {
          if ((locant_ch = parse_locant()) == WLN_ERROR)
            return WLN_ERROR; 
        }
        else if (*wln_ptr >= '1' && *wln_ptr <= '9')
          goto ring_parse_multi;
        else 
          return WLN_ERROR; 
        break; 

      case '&':
        ring = ring_create(mol, max_path_size); 
        goto ring_parse_arom;
      case 'T': 
        ring = ring_create(mol, max_path_size); 
        goto ring_parse_arom;

      case 'J': 
        ring = ring_create(mol, max_path_size); 
        goto ring_parse_end;

      case 'B':
      case 'O':
      case 'S':
      case 'N':
      case 'K':
      case 'M':
        ring = ring_create(mol, max_path_size); 
        goto ring_parse_hetero; 
    }
  }

ring_parse_multi:
  ch = *wln_ptr; 
  wln_ptr += ch - '0' + 1; // skip the multi-block
  if (*wln_ptr != ' ')
    return wln_error("invalid notation for ring multi block\n"); 
  wln_ptr++; 
  max_path_size = parse_locant(); 
  if (max_path_size == WLN_ERROR)
    return WLN_ERROR; 
  
  ring = ring_create(mol, max_path_size+1); 
  ch = *wln_ptr; 
  switch (ch) {
      case ' ': goto ring_parse_hetero; 
      case 'J': goto ring_parse_end;

      case '&':
        goto ring_parse_arom;
      case 'T': 
        goto ring_parse_arom;
    
      case 'B':
      case 'O':
      case 'S':
      case 'N':
      case 'K':
      case 'M':
        goto ring_parse_hetero; 
  }

  while (*wln_ptr) {
    ch = *wln_ptr++; 
ring_parse_hetero:
    switch (ch) {
      case ' ':
        if (*wln_ptr >= 'A' && *wln_ptr <= 'Z') {
          if ((locant_ch = parse_locant()) == WLN_ERROR)
            return WLN_ERROR; 
          if (locant_ch >= ring->size)
            return wln_error("hetero-atom assignment out of bounds\n");
        }
        else if (*wln_ptr >= '1' && *wln_ptr <= '9')
          goto ring_parse_multi;
        else 
          return WLN_ERROR; 

      case 'H': break;
      case 'B': symbol_change(ring->path[locant_ch++].s, BOR); break;

      case 'K':
        symbol_change(ring->path[locant_ch].s, NIT);
        symbol_set_charge(ring->path[locant_ch++].s, +1);
        break;

      case 'N':
      case 'M': symbol_change(ring->path[locant_ch++].s, NIT); break;
      case 'O': symbol_change(ring->path[locant_ch++].s, OXY); break;
      case 'P': symbol_change(ring->path[locant_ch++].s, PHO); break;
      case 'S': symbol_change(ring->path[locant_ch++].s, SUL); break;

      case 'U':
#if 0
        if (locant_ch > ring->size) 
          return wln_error("could not unsaturate bond due - locant out of bounds");
        if (locant_ch == ring->size)
          bond = graph_get_edge(mol, ring->path[locant_ch].s, 
                                ring->path[0].s);
        else
          bond = graph_get_edge(mol, ring->path[locant_ch].s, 
                                ring->path[locant_ch+1].s);
        edge_unsaturate(bond);
#endif
        break;

      case '-':
        break;
      case 'X':
      case 'Y':
        break;

      case 'V': if (add_oxy(mol, ring->path[locant_ch++].s) != WLN_OK) 
                  return WLN_ERROR; 
                else break; 
      
      case '&':
        goto ring_parse_arom;
      case 'T': 
        goto ring_parse_arom;

      case 'J': goto ring_parse_end;

      default:
        return wln_error("unhandled symbol - %c\n", ch); 
    }
  }


ring_parse_arom:


ring_parse_end:
  if (seen_pseudo)
    pathsolverIII(mol, ring, SSSR, nSSSR);  
  else
    pathsolverIII_fast(mol, ring, SSSR, nSSSR);  

  wln_ptr++; 
  if (*wln_ptr == ' ') {
    wln_ptr++; 
    if (*wln_ptr == '&') {
      free(ring); 
      return start_wln_parse(mol); 
    }
    else if (parse_ring_locants(mol, ring) != WLN_OK)
      return WLN_ERROR; 
  }
  free(ring); 
  return WLN_OK; 
}


/* this expects the entry character to be the character after the locant space */
static int parse_ring_locants(graph_t *mol, ring_t *ring)
{
  unsigned char ch = *wln_ptr; 
  if (!*wln_ptr)
    return WLN_ERROR;

  while (*wln_ptr) {
    if (*wln_ptr == '&') {
      wln_ptr++; 
      return start_wln_parse(mol); 
    }

    int locant = parse_locant(); 

    if (locant == WLN_ERROR) 
      return wln_error("failed to parse locant\n");  
    if (locant >= ring->size) 
      return wln_error("locant out of range of ring\n"); 
  
    symbol_t *curr_symbol = ring->path[locant].s; 
    edge_t *curr_edge = edge_create(mol, curr_symbol); 

    switch (ch) {
      case '-':
        /* inline ring definition */
        break; 
      
      default:
        if (branch_recursive_parse(mol, curr_edge, ring) == WLN_ERROR)
          return WLN_ERROR; 
        if (!edge_get_end(curr_edge)) {
          symbol_t *methyl = symbol_create(mol, CAR);
          edge_bond(mol, curr_edge, methyl);  
        }

        /* move to the next locant */
        if (*wln_ptr == ' ')
          wln_ptr++; 
        if (*wln_ptr == '&') {
          /* move out of the current ring */
          wln_ptr++; 
          return WLN_OK; 
        }
        break; 
    }
  }

  return WLN_OK; 
}


static int branch_recursive_parse(graph_t *mol, edge_t *edge, ring_t *ring)
{
  int ret; 
  unsigned char ch = *wln_ptr; 
  edge_t *curr_edge = edge;
  symbol_t *curr_symbol;
  ring_t *benzene; 
  uint16_t counter   = 0; 

  if (!*wln_ptr)
    return WLN_OK;

  while (*wln_ptr) {
    ch  = *wln_ptr++;
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
        while ((ch = *wln_ptr) && ch >= '0' && ch <= '9') {
          counter *= 10; 
          counter += ch - '0'; 
          wln_ptr++;
        }
        for (uint16_t i=0; i<counter; i++) {
          symbol_t *carbon = symbol_create(mol, CAR);
          edge_bond(mol, curr_edge, carbon); 
          curr_edge = edge_create(mol, carbon); 
          curr_symbol = carbon; 
        }
        break;
      
      case 'A':
      case 'J':
        return wln_error("non-atomic symbol used in chain"); 
      
      case 'B':
        curr_symbol = symbol_create(mol, BOR);
        edge_bond(mol, curr_edge, curr_symbol); 

        for (unsigned int i=0; i<2; i++) {
          curr_edge = edge_create(mol, curr_symbol); 
          ret = branch_recursive_parse(mol, curr_edge, ring); 
          if (ret == WLN_ERROR)
            return ret; 
          if (!*wln_ptr)
            break; 
        }
        return WLN_OK; 
      
      case 'C':
        return wln_error("WLN symbol C currently unhandled\n"); 

      // extrememly rare open chelate notation
      case 'D':
        return wln_error("WLN symbol D (chelate) currently unhandled\n"); 
      
      // terminator symbol 
      case 'E':
        curr_symbol = symbol_create(mol, BRO);
        edge_bond(mol, curr_edge, curr_symbol); 
        return WLN_OK; 

      case 'F':
        curr_symbol = symbol_create(mol, FLU);
        edge_bond(mol, curr_edge, curr_symbol); 
        return WLN_OK; 

      case 'G':
        curr_symbol = symbol_create(mol, CHL);
        edge_bond(mol, curr_edge, curr_symbol); 
        return WLN_OK; 

      case 'H':
        break;

      case 'I':
        curr_symbol = symbol_create(mol, IOD);
        edge_bond(mol, curr_edge, curr_symbol); 
        return WLN_OK; 
     
      /* [N+](R)(R)(R)(R) */
      case 'K':
        curr_symbol = symbol_create(mol, NIT);
        edge_bond(mol, curr_edge, curr_symbol); 
        symbol_set_charge(curr_symbol, +1); 

        for (unsigned int i=0; i<3; i++) {
          curr_edge = edge_create(mol, curr_symbol); 
          ret = branch_recursive_parse(mol, curr_edge, ring); 
          if (ret == WLN_ERROR)
            return ret; 
          if (!edge_get_end(curr_edge)) {
            symbol_t *methyl = symbol_create(mol, CAR);
            edge_bond(mol, curr_edge, methyl);  
          }
        }
        return WLN_OK; 

      case 'L':
      case 'T':
        /* impossible from here */
        return wln_error("ring notation must start the molecule to be used\n");
      
      /* NH(R)(R) */
      case 'M':
        curr_symbol = symbol_create(mol, NIT); 
        edge_bond(mol, curr_edge, curr_symbol); 
        curr_edge = edge_create(mol, curr_symbol); 
        break;
     
      /* NR(R)(R) */
      case 'N':
        curr_symbol = symbol_create(mol, NIT);
        edge_bond(mol, curr_edge, curr_symbol); 
        for (unsigned int i=0; i<2; i++) {
          curr_edge = edge_create(mol, curr_symbol); 
          ret = branch_recursive_parse(mol, curr_edge, ring); 
          if (ret == WLN_ERROR)
            return ret; 
          if (!*wln_ptr)
            break; 
        }
        return WLN_OK; 
      
      /* OR(R) */
      case 'O':
        curr_symbol = symbol_create(mol, OXY); 
        edge_bond(mol, curr_edge, curr_symbol); 
        curr_edge = edge_create(mol, curr_symbol); 
        break;
      
      case 'P':
        curr_symbol = symbol_create(mol, PHO);
        edge_bond(mol, curr_edge, curr_symbol); 

        for (unsigned int i=0; i<3; i++) {
          curr_edge = edge_create(mol, curr_symbol); 
          ret = branch_recursive_parse(mol, curr_edge, ring); 
          if (ret == WLN_ERROR)
            return ret; 
          if (!*wln_ptr)
            break; 
        }
        return WLN_OK; 

      case 'Q':
        curr_symbol = symbol_create(mol, OXY);
        edge_bond(mol, curr_edge, curr_symbol); 
        return WLN_OK; 

      // shorthand benzene 
      case 'R': 
        benzene = ring_create_benzene(mol);
        curr_symbol = benzene->path[0].s; 
        edge_bond(mol, curr_edge, curr_symbol); 
        if (*wln_ptr == ' ') {
          wln_ptr++; 
          return parse_ring_locants(mol, benzene); 
        }
        else
          curr_edge = edge_create(mol, curr_symbol); 
        break;

      case 'S':
        curr_symbol = symbol_create(mol, SUL);
        edge_bond(mol, curr_edge, curr_symbol); 

        for (unsigned int i=0; i<3; i++) {
          curr_edge = edge_create(mol, curr_symbol); 
          ret = branch_recursive_parse(mol, curr_edge, ring); 
          if (ret == WLN_ERROR)
            return ret; 
          if (!*wln_ptr)
            break; 
        }
        return WLN_OK; 
      
      case 'U':
        edge_unsaturate(curr_edge); 
        break; 

      case 'V':
        curr_symbol = symbol_create(mol, CAR);
        edge_bond(mol, curr_edge, curr_symbol); 
        if (add_oxy(mol, curr_symbol) != WLN_OK)
          return wln_error("failed to add =O group\n"); 
        curr_edge = edge_create(mol, curr_symbol); 
        break;
      
      // if not previous, create dummy carbon
      // really trying to avoid a bit state on this
      case 'W':
#if 0
        if (!add_dioxo(mol, prev_symbol))
          return error("failed to add dioxo group with W\n"); 
#endif
        break; 

      case 'X':
        curr_symbol = symbol_create(mol, CAR);
        edge_bond(mol, curr_edge, curr_symbol); 

        for (unsigned int i=0; i<3; i++) {
          curr_edge = edge_create(mol, curr_symbol); 
          ret = branch_recursive_parse(mol, curr_edge, ring); 
          if (ret == WLN_ERROR)
            return ret; 
          if (!edge_get_end(curr_edge)) {
            symbol_t *methyl = symbol_create(mol, CAR);
            edge_bond(mol, curr_edge, methyl);  
          }
        }

        return WLN_OK; 

      case 'Y':
        curr_symbol = symbol_create(mol, CAR);
        edge_bond(mol, curr_edge, curr_symbol); 
        
        for (unsigned int i=0; i<2; i++) {
          curr_edge = edge_create(mol, curr_symbol); 
          int ret = branch_recursive_parse(mol, curr_edge, ring); 
          if (ret == WLN_ERROR)
            return ret; 
          if (!edge_get_end(curr_edge)) {
            symbol_t *methyl = symbol_create(mol, CAR);
            edge_bond(mol, curr_edge, methyl);  
          }
        }
        return WLN_OK; 
      
      case 'Z':
        curr_symbol = symbol_create(mol, NIT);
        edge_bond(mol, curr_edge, curr_symbol); 
        return WLN_OK; 
      
      case '-':
#if 1
        if (*wln_ptr == ' ') {
          /* must be an inline ring */
        }
#if 0
          if (!curr_ring) return false; 
          stack_ptr = depstack_push_ring(dep_stack, stack_ptr, curr_ring);  
#endif
        else {
          curr_symbol = dash_symbol_create(mol, &wln_ptr); 
          edge_bond(mol, curr_edge, curr_symbol); 
          if (!curr_symbol)
            return wln_error("invalid elemental code - %s\n", wln_ptr);

          for (unsigned int i=0; i<3; i++) {
            curr_edge = edge_create(mol, curr_symbol); 
            ret = branch_recursive_parse(mol, curr_edge, ring); 
            if (ret == WLN_ERROR)
              return ret; 
            if (!*wln_ptr)
              break; 
          }
          return WLN_OK; 
        }
#endif
        break;

      case ' ':
        if (*wln_ptr == '&') {
          wln_ptr++; 
          return start_wln_parse(mol); 
        }
        else {
          if (!ring) 
            return wln_error("opening ring notation without a prior ring\n"); 
          return parse_ring_locants(mol, ring); 
        }
        break; 

      case '&':
        /* branch must of closed */
        return WLN_OK;

      case '\n':
        break; 
      case '/':
        return wln_error("slash seen outside of ring - multipliers currently unsupported\n");
      default:
        fprintf(stderr, "invalid character read for WLN notation - %c(%u)\n", ch, ch);
        return false; 
    }
  }
  
  return WLN_OK; 
}


// init conditions, make one dummy atom, and one bond - work of the virtual bond
// idea entirely *--> grow..., delete at the end to save branches
static int start_wln_parse(graph_t *mol)
{
  symbol_t *init_symbol = symbol_create(mol, DUM); 
  edge_t *init_edge = graph_new_edge(mol); 
  edge_t *null_edge = init_edge;
  edge_set_begin(init_edge, init_symbol); 

  symbol_t *open_term = parse_opening_terminator(mol, *wln_ptr);
  if (open_term) {
    edge_bond(mol, init_edge, open_term); 
    init_edge = edge_create(mol, open_term); 
    wln_ptr++; 
  }

  switch (*wln_ptr) {
    case 'L':
    case 'T':
      if (cyclic_recursive_parse(mol) == WLN_ERROR)
        return WLN_ERROR; 
      break; 
    default:
      if (branch_recursive_parse(mol, init_edge, NULL) == WLN_ERROR)
        return WLN_ERROR; 
  }

  graph_delete_symbol(mol, init_symbol); 
  if (graph_num_atoms(mol) == 0)
    return wln_error("empty molecule from wln parse\n");
  return WLN_OK;
}


/*
 * -- Parse WLN Notation --
 */
bool ReadWLN(const char *wln, graph_t *mol)
{
#ifdef USING_OPENBABEL
  mol->BeginModify(); 
  mol->SetAromaticPerceived(true);
  mol->SetChiralityPerceived(true); // no stereo for WLN
#endif
  
  wln_ptr = (char*)wln; 
  if (start_wln_parse(mol) == WLN_ERROR)
    return false; 

  if (*wln_ptr) {
    wln_error("parse ended before end of notation - %s\n", wln_ptr); 
    return false; 
  }

  graph_cleanup_hydrogens(mol);
  
#ifdef USING_OPENBABEL
  OBKekulize(mol); 
#endif
  return true; 
}


