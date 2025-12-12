
/*********************************************************************
 
Author : Michael Blakey

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

#include <openbabel/mol.h>
#include <openbabel/plugin.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/kekulize.h>
#include <openbabel/ring.h>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include "openbabel/stereo/stereo.h"
#include <openbabel/stereo/tetrahedral.h>

#include "wlnparser.h"

#define DEBUG
#define WLN_OK 0
#define WLN_ERROR -1

#ifdef USING_OPENBABEL
using namespace OpenBabel;
#define graph_t    OBMol
#define symbol_t   OBAtom
#define edge_t     OBBond
#define ring_t     OBRing 

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
#define symbol_get_degree(s)        s->GetExplicitDegree()
#define symbol_is_cyclic(s)         s->IsInRing()

#define edge_unsaturate(e)          e->SetBondOrder(1+e->GetBondOrder())
#define edge_set_order(e, o)        e->SetBondOrder(o)
#define edge_set_aromatic(e, b)     e->SetAromatic(b)
#define edge_get_order(e)           e->GetBondOrder()

#define edge_set_begin(e, s)       e->SetBegin(s)
#define edge_set_end(e, s)         e->SetEnd(s)
#define edge_get_begin(e)          e->GetBeginAtom()
#define edge_get_end(e)            e->GetEndAtom()
#define edge_get_nbor(e, s)        e->GetNbrAtom(s)

#define ring_get_size(r)            r->Size()
#define ring_get_first(r, g)        g->GetAtom(r->_path[0])
#define ring_get_id(r)              r->ring_id
#define ring_symbol_member(r,s)     r->IsMember(s)
#define ring_edge_member(r,e)       r->IsMember(e)

#define graph_new_symbol(g)         g->NewAtom()
#define graph_new_edge(g)           g->NewBond()
#define graph_set_edge(g,e)         g->AddBond(*e)
#define graph_get_edge(g,x,y)       g->GetBond(x,y)
#define graph_delete_symbol(g,s)    g->DeleteAtom(s)
#define graph_delete_edge(g,e)      g->DeleteBond(e)
#define graph_num_atoms(g)          g->NumAtoms()
#define graph_num_cycles(g)         g->GetSSSR().size()

#define symbol_nbor_iter(a, s)      FOR_NBORS_OF_ATOM(a,s)
#define symbol_bond_iter(b, s)      FOR_BONDS_OF_ATOM(b,s)

#define graph_symbol_iter(s, g)     FOR_ATOMS_OF_MOL(s,g)
#define graph_edge_iter(e, g)       FOR_BONDS_OF_MOL(e,g)
#define graph_ring_iter(r, g)       FOR_RINGS_OF_MOL(r,g)

#elif defined USING_RDKIT

#endif

#define int_to_locant(X) (X+64)
#define locant_to_int(X) (X-64)

char *wln_out; 
unsigned int wln_len; 
bool *seen; 

#define wln_push(ch)        wln_out[wln_len++] = ch
#define wln_pop()           wln_out[--wln_len] = '\0'
#define wln_back()          wln_out[wln_len-1]
#define wln_terminate()     wln_out[wln_len] = '\0'
#define wln_length()        wln_len
#define wln_at(n)           wln_out[n]
#define wln_set(ch, n)      wln_out[n] = ch
#define wln_ptr(n)          &wln_out[n]
#define wln_end_ptr()       &wln_out[wln_len-1]


static int locant_recursive_write(struct wln_ring *wln_ring, graph_t *mol); 
static int branch_recursive_write(graph_t *mol, symbol_t *atom); 


static bool is_wln_V(symbol_t *atom)
{
  if (symbol_get_valence(atom) == 4 &&
      symbol_get_degree(atom) == 3  &&
      symbol_get_charge(atom) == 0) 
  {
    
    symbol_bond_iter(b, atom) {
      edge_t *edge = &(*b); 
      symbol_t *child = edge_get_nbor(edge, atom); 
      if (edge_get_order(edge) == 2 &&
          symbol_get_num(child) == 8) 
      {
        seen[symbol_get_id(child)] = true; 
        return true; 
      }
    }      
  } 

  return false; 
}

// if modern, charges are completely independent apart from assumed K
static void wln_write_element_symbol(OBAtom* atom) {
  const unsigned int neighbours = atom->GetExplicitDegree(); 
  const unsigned int orders     = atom->GetExplicitValence(); 
  const unsigned int hcount     = atom->GetImplicitHCount(); 
  const int charge              = atom->GetFormalCharge(); 

  switch(atom->GetAtomicNum()){
    case 1: wln_push('H'); break; 

    case 5:
      if (orders > 3) {
        wln_push('-');  
        wln_push('B');  
        wln_push('-');  
      }
      else 
        wln_push('H'); 
      break; 

    case 6:
      if (is_wln_V(atom)) 
        wln_push('V'); 
      else if (neighbours <= 2)
        wln_push('1');
      else if(neighbours == 3)
        wln_push('Y');
      else if(neighbours == 4)
        wln_push('X');
      break; 
    
    case 7:
        if (orders <= 1 && hcount == 2)
          wln_push('Z');
        else if (orders == 2 && hcount == 1)
          wln_push('M');
        else if (charge == +1 && orders == 4)
          wln_push('K');
        else if (orders >= 4) {
          wln_push('-');  
          wln_push('N');  
          wln_push('-');  
        }
        else 
          wln_push('N');  
        break; 
    
    case 8:
      if (neighbours == 1 && orders==1 && charge == 0)
        wln_push('Q');
      else if (neighbours == 0 && charge != -2)
        wln_push('Q');
      else if (orders > 2) {
        wln_push('-');  
        wln_push('O');  
        wln_push('-');  
      }
      else
        wln_push('O');  
      break; 
        
    case 9:
      if (neighbours > 1) {
        wln_push('-');  
        wln_push('F');  
        wln_push('-');  
      }
      else 
        wln_push('F');  
      break; 

    case 15:
      if (neighbours > 5) {
        wln_push('-');  
        wln_push('P');  
        wln_push('-');  
      }
      else 
        wln_push('P');  
      break; 

    case 16: wln_push('S'); break; 

    case 17:
      if (neighbours > 1) {
        wln_push('-');  
        wln_push('G');  
        wln_push('-');  
      }
      else 
        wln_push('G');  
      break; 

    case 35:
      if (neighbours > 1) {
        wln_push('-');  
        wln_push('E');  
        wln_push('-');  
      }
      else 
        wln_push('E');  
      break; 

    case 53:
      if (neighbours > 1) {
        wln_push('-');  
        wln_push('I');  
        wln_push('-');  
      }
      else 
        wln_push('I');  
      break; 

      case 89: wln_push('-'); wln_push('A'); wln_push('C'); wln_push('-'); break;
      case 47: wln_push('-'); wln_push('A'); wln_push('G'); wln_push('-'); break;
      case 13: wln_push('-'); wln_push('A'); wln_push('L'); wln_push('-'); break;
      case 95: wln_push('-'); wln_push('A'); wln_push('M'); wln_push('-'); break;
      case 18: wln_push('-'); wln_push('A'); wln_push('R'); wln_push('-'); break;
      case 33: wln_push('-'); wln_push('A'); wln_push('S'); wln_push('-'); break;
      case 85: wln_push('-'); wln_push('A'); wln_push('T'); wln_push('-'); break;
      case 79: wln_push('-'); wln_push('A'); wln_push('U'); wln_push('-'); break;
      case 56: wln_push('-'); wln_push('B'); wln_push('A'); wln_push('-'); break;
      case 4: wln_push('-'); wln_push('B'); wln_push('E'); wln_push('-'); break;
      case 107: wln_push('-'); wln_push('B'); wln_push('H'); wln_push('-'); break;
      case 83: wln_push('-'); wln_push('B'); wln_push('I'); wln_push('-'); break;
      case 97: wln_push('-'); wln_push('B'); wln_push('K'); wln_push('-'); break;
      case 20: wln_push('-'); wln_push('C'); wln_push('A'); wln_push('-'); break;
      case 48: wln_push('-'); wln_push('C'); wln_push('D'); wln_push('-'); break;
      case 58: wln_push('-'); wln_push('C'); wln_push('E'); wln_push('-'); break;
      case 98: wln_push('-'); wln_push('C'); wln_push('F'); wln_push('-'); break;
      case 96: wln_push('-'); wln_push('C'); wln_push('N'); wln_push('-'); break;
      case 112: wln_push('-'); wln_push('C'); wln_push('N'); wln_push('-'); break;
      case 27: wln_push('-'); wln_push('C'); wln_push('O'); wln_push('-'); break;
      case 24: wln_push('-'); wln_push('C'); wln_push('R'); wln_push('-'); break;
      case 55: wln_push('-'); wln_push('C'); wln_push('S'); wln_push('-'); break;
      case 29: wln_push('-'); wln_push('C'); wln_push('U'); wln_push('-'); break;
      case 105: wln_push('-'); wln_push('D'); wln_push('B'); wln_push('-'); break;
      case 110: wln_push('-'); wln_push('D'); wln_push('S'); wln_push('-'); break;
      case 66: wln_push('-'); wln_push('D'); wln_push('Y'); wln_push('-'); break;
      case 68: wln_push('-'); wln_push('E'); wln_push('R'); wln_push('-'); break;
      case 99: wln_push('-'); wln_push('E'); wln_push('S'); wln_push('-'); break;
      case 63: wln_push('-'); wln_push('E'); wln_push('U'); wln_push('-'); break;
      case 26: wln_push('-'); wln_push('F'); wln_push('E'); wln_push('-'); break;
      case 114: wln_push('-'); wln_push('F'); wln_push('L'); wln_push('-'); break;
      case 100: wln_push('-'); wln_push('F'); wln_push('M'); wln_push('-'); break;
      case 87: wln_push('-'); wln_push('F'); wln_push('R'); wln_push('-'); break;
      case 31: wln_push('-'); wln_push('G'); wln_push('A'); wln_push('-'); break;
      case 64: wln_push('-'); wln_push('G'); wln_push('D'); wln_push('-'); break;
      case 32: wln_push('-'); wln_push('G'); wln_push('E'); wln_push('-'); break;
      case 2: wln_push('-'); wln_push('H'); wln_push('E'); wln_push('-'); break;
      case 72: wln_push('-'); wln_push('H'); wln_push('F'); wln_push('-'); break;
      case 80: wln_push('-'); wln_push('H'); wln_push('G'); wln_push('-'); break;
      case 67: wln_push('-'); wln_push('H'); wln_push('O'); wln_push('-'); break;
      case 108: wln_push('-'); wln_push('H'); wln_push('S'); wln_push('-'); break;
      case 49: wln_push('-'); wln_push('I'); wln_push('N'); wln_push('-'); break;
      case 77: wln_push('-'); wln_push('I'); wln_push('R'); wln_push('-'); break;
      case 36: wln_push('-'); wln_push('K'); wln_push('R'); wln_push('-'); break;
      case 19: wln_push('-'); wln_push('K'); wln_push('A'); wln_push('-'); break;
      case 57: wln_push('-'); wln_push('L'); wln_push('A'); wln_push('-'); break;
      case 3: wln_push('-'); wln_push('L'); wln_push('I'); wln_push('-'); break;
      case 103: wln_push('-'); wln_push('L'); wln_push('R'); wln_push('-'); break;
      case 71: wln_push('-'); wln_push('L'); wln_push('U'); wln_push('-'); break;
      case 116: wln_push('-'); wln_push('L'); wln_push('V'); wln_push('-'); break;
      case 115: wln_push('-'); wln_push('M'); wln_push('C'); wln_push('-'); break;
      case 101: wln_push('-'); wln_push('M'); wln_push('D'); wln_push('-'); break;
      case 12: wln_push('-'); wln_push('M'); wln_push('G'); wln_push('-'); break;
      case 25: wln_push('-'); wln_push('M'); wln_push('N'); wln_push('-'); break;
      case 42: wln_push('-'); wln_push('M'); wln_push('O'); wln_push('-'); break;
      case 109: wln_push('-'); wln_push('M'); wln_push('T'); wln_push('-'); break;
      case 11: wln_push('-'); wln_push('N'); wln_push('A'); wln_push('-'); break;
      case 41: wln_push('-'); wln_push('N'); wln_push('B'); wln_push('-'); break;
      case 60: wln_push('-'); wln_push('N'); wln_push('D'); wln_push('-'); break;
      case 10: wln_push('-'); wln_push('N'); wln_push('E'); wln_push('-'); break;
      case 113: wln_push('-'); wln_push('N'); wln_push('H'); wln_push('-'); break;
      case 28: wln_push('-'); wln_push('N'); wln_push('I'); wln_push('-'); break;
      case 102: wln_push('-'); wln_push('N'); wln_push('O'); wln_push('-'); break;
      case 93: wln_push('-'); wln_push('N'); wln_push('P'); wln_push('-'); break;
      case 118: wln_push('-'); wln_push('O'); wln_push('G'); wln_push('-'); break;
      case 76: wln_push('-'); wln_push('O'); wln_push('S'); wln_push('-'); break;
      case 91: wln_push('-'); wln_push('P'); wln_push('A'); wln_push('-'); break;
      case 82: wln_push('-'); wln_push('P'); wln_push('B'); wln_push('-'); break;
      case 46: wln_push('-'); wln_push('P'); wln_push('D'); wln_push('-'); break;
      case 61: wln_push('-'); wln_push('P'); wln_push('M'); wln_push('-'); break;
      case 84: wln_push('-'); wln_push('P'); wln_push('O'); wln_push('-'); break;
      case 59: wln_push('-'); wln_push('P'); wln_push('R'); wln_push('-'); break;
      case 78: wln_push('-'); wln_push('P'); wln_push('T'); wln_push('-'); break;
      case 94: wln_push('-'); wln_push('P'); wln_push('U'); wln_push('-'); break;
      case 88: wln_push('-'); wln_push('R'); wln_push('A'); wln_push('-'); break;
      case 37: wln_push('-'); wln_push('R'); wln_push('B'); wln_push('-'); break;
      case 75: wln_push('-'); wln_push('R'); wln_push('E'); wln_push('-'); break;
      case 104: wln_push('-'); wln_push('R'); wln_push('F'); wln_push('-'); break;
      case 111: wln_push('-'); wln_push('R'); wln_push('G'); wln_push('-'); break;
      case 45: wln_push('-'); wln_push('R'); wln_push('H'); wln_push('-'); break;
      case 86: wln_push('-'); wln_push('R'); wln_push('N'); wln_push('-'); break;
      case 44: wln_push('-'); wln_push('R'); wln_push('U'); wln_push('-'); break;
      case 51: wln_push('-'); wln_push('S'); wln_push('B'); wln_push('-'); break;
      case 21: wln_push('-'); wln_push('S'); wln_push('C'); wln_push('-'); break;
      case 34: wln_push('-'); wln_push('S'); wln_push('E'); wln_push('-'); break;
      case 106: wln_push('-'); wln_push('S'); wln_push('G'); wln_push('-'); break;
      case 14: wln_push('-'); wln_push('S'); wln_push('I'); wln_push('-'); break;
      case 62: wln_push('-'); wln_push('S'); wln_push('M'); wln_push('-'); break;
      case 50: wln_push('-'); wln_push('S'); wln_push('N'); wln_push('-'); break;
      case 38: wln_push('-'); wln_push('S'); wln_push('R'); wln_push('-'); break;
      case 73: wln_push('-'); wln_push('T'); wln_push('A'); wln_push('-'); break;
      case 65: wln_push('-'); wln_push('T'); wln_push('B'); wln_push('-'); break;
      case 43: wln_push('-'); wln_push('T'); wln_push('C'); wln_push('-'); break;
      case 52: wln_push('-'); wln_push('T'); wln_push('E'); wln_push('-'); break;
      case 90: wln_push('-'); wln_push('T'); wln_push('H'); wln_push('-'); break;
      case 22: wln_push('-'); wln_push('T'); wln_push('I'); wln_push('-'); break;
      case 81: wln_push('-'); wln_push('T'); wln_push('L'); wln_push('-'); break;
      case 69: wln_push('-'); wln_push('T'); wln_push('M'); wln_push('-'); break;
      case 117: wln_push('-'); wln_push('T'); wln_push('S'); wln_push('-'); break;
      case 92: wln_push('-'); wln_push('U'); wln_push('R'); wln_push('-'); break;
      case 23: wln_push('-'); wln_push('V'); wln_push('A'); wln_push('-'); break;
      case 74: wln_push('-'); wln_push('W'); wln_push('T'); wln_push('-'); break;
      case 54: wln_push('-'); wln_push('X'); wln_push('E'); wln_push('-'); break;
      case 39: wln_push('-'); wln_push('Y'); wln_push('T'); wln_push('-'); break;
      case 70: wln_push('-'); wln_push('Y'); wln_push('B'); wln_push('-'); break;
      case 30: wln_push('-'); wln_push('Z'); wln_push('N'); wln_push('-'); break;
      case 40: wln_push('-'); wln_push('Z'); wln_push('R'); wln_push('-'); break;
  }
}


struct wln_ring {
  unsigned int size;
  unsigned int nsssr;
  unsigned int fsum;
  unsigned long long multi; // bitset 
  bool hetero;
  ring_t *sssr[32]; 
  symbol_t *locants[1];
}; 


static struct wln_ring* wln_ring_alloc(size_t size)
{
  const size_t bytes = sizeof(struct wln_ring) + (size-1)*sizeof(symbol_t*); 
  struct wln_ring *r = (struct wln_ring*)malloc(bytes); 

  r->size = size; 
  r->nsssr  = 0; 
  r->fsum   = 0xFFFFFFFF; 
  r->multi  = 0; 
  r->hetero = false; 
  return r; 
}


static void wln_ring_free(struct wln_ring *r)
{
  free(r); 
}


static void wln_ring_add_subcycle(struct wln_ring *r, ring_t *ring)
{
  r->sssr[r->nsssr++] = ring; 
}


static unsigned int symbol_ring_share_count(struct wln_ring *r, symbol_t *s) 
{
  unsigned int memb = 0; 
  for (unsigned int i=0; i<r->nsssr;i++) {
    ring_t *ring = r->sssr[i]; 
    memb += ring_symbol_member(ring, s); 
  }
  return memb; 
}


/* will return the size of the local SSSR sum */
static int walk_ring_recursive(struct wln_ring *wln_ring, 
                               graph_t *mol, 
                               symbol_t *parent, 
                               bool *ring_set, 
                               int ratoms)
{
  if (symbol_get_num(parent) != 6)
    wln_ring->hetero = true; 

  seen[symbol_get_id(parent)] = true; 
  symbol_nbor_iter(a, parent) {
    symbol_t *nbor = &(*a); 
    const unsigned int nid = symbol_get_id(nbor); 
    if (symbol_is_cyclic(nbor) && !seen[nid]) {
      /* get the ring membership */
      graph_ring_iter(r, mol) {
        ring_t *sssr_ring = &(*r);
        if (!ring_symbol_member(sssr_ring, nbor))
          continue; 

        if (!ring_set[ring_get_id(sssr_ring)]) {
          wln_ring_add_subcycle(wln_ring, sssr_ring); 
          ring_set[ring_get_id(sssr_ring)] = true; 
          break;
        }
      }
      wln_ring->locants[ratoms] = nbor;
      ratoms = walk_ring_recursive(wln_ring, mol, nbor, ring_set, ratoms+1);
    }
  }
  
  return ratoms;
}


static void wln_ring_fill_sssr(struct wln_ring *wln_ring,
                               graph_t *mol,   
                               symbol_t *init_atom)
{
  const unsigned int ncycles = graph_num_cycles(mol); 
  
  bool *added = (bool*)alloca(ncycles);
  memset(added, false, ncycles); 
  
  wln_ring->locants[0] = init_atom;
  unsigned int size = walk_ring_recursive(wln_ring, mol, init_atom, added, 1); 
  wln_ring->size = size; 

  for (unsigned int i=0; i<size; i++)
    seen[symbol_get_id(wln_ring->locants[i])] = true; 

#ifdef DEBUG
  fprintf(stderr, "WLN cycle: %lu atoms, %lu SSSR\n", size, wln_ring->nsssr); 
#endif
}


/*
Fusion sums are a way of comparing two unique locant paths in order to provide a unique solution. 
The sum is calculated by taking each individual sub-cycle in the local SSSR, and summing the lowest 
locant value contained in that cycle from the calculated path. For example if a sub-cycle contained the 
locant 'A', its locant sum value would be zero.
*/
static int fusion_sum_score_path(struct wln_ring *r, symbol_t **path)
{
  unsigned int sum = 0; 
  for (unsigned int i=0; i<r->nsssr; i++) {
    ring_t *subcycle = r->sssr[i]; 
    for (unsigned int j=0; j<r->size; j++) {
      symbol_t *atom = path[j]; 
      if (ring_symbol_member(subcycle, atom)) {
        sum += j;
        break;
      } 
    }
  }
  return sum; 
}


static int ring_share_score_path(struct wln_ring *r, symbol_t **path, unsigned int *shares)
{
  unsigned int sum = 0; 
  for (unsigned int j=0; j<r->size; j++) 
    sum += j * shares[j]; 
  return sum; 
}


/* flood fill style, take every path, and minimise the fusion sum */
static bool fusion_sum_traverse_recursive(struct wln_ring *r, 
                                          graph_t *mol, 
                                          symbol_t *parent, 
                                          symbol_t **path,
                                          unsigned int *shares,
                                          bool *local_seen, 
                                          int id)
{
  if (id == r->size) {
    /* score and copy path */
    unsigned int fsum = fusion_sum_score_path(r, path); 
    if (fsum < r->fsum) {
      memcpy(r->locants, path, r->size*sizeof(symbol_t*)); 
      r->fsum = fsum; 
    }
    
    /* break ties with share sum */
    if (fsum == r->fsum) {
      unsigned int orig_share = ring_share_score_path(r, r->locants, shares); 
      if (ring_share_score_path(r, path, shares) < orig_share) 
        memcpy(r->locants, path, r->size*sizeof(symbol_t*)); 
    }
    return true; 
  }

  bool ret = false;
  local_seen[symbol_get_id(parent)] = true; 

  symbol_nbor_iter(a, parent) {
    symbol_t *nbor = &(*a); 
    const unsigned int nid = symbol_get_id(nbor); 
    if (symbol_is_cyclic(nbor) && !local_seen[nid]) {
      local_seen[nid] = true;
      
      path[id++] = nbor;
      ret |= fusion_sum_traverse_recursive(r, mol, nbor, path, shares, local_seen, id); 
      id--;  
      local_seen[nid] = false;
    }
  }

  return ret; 
}


static bool wln_ring_fill_locant_path(struct wln_ring *r, 
                                      graph_t *mol)
{
  const unsigned int size = r->size;
  unsigned int *shares = (unsigned int*)malloc(sizeof(unsigned int)*size); 
  memset(shares, 0, sizeof(unsigned int)*size);  
  
  const unsigned int natoms = graph_num_atoms(mol); 

  bool *local_seen = (bool*)malloc(natoms); 
  memset(local_seen, false, natoms); 

  symbol_t **ordered_path = (symbol_t**)malloc(sizeof(symbol_t*)*size); 
  
  symbol_t *start_symbol; 
  unsigned int max_share = 0; 
  for (unsigned int i=0; i<size; i++) {
    shares[i] = symbol_ring_share_count(r, r->locants[i]); 
    if (max_share < shares[i]) {
      max_share = shares[i];
      start_symbol = r->locants[i]; 
    }
  }

  local_seen[symbol_get_id(start_symbol)] = true; 
  ordered_path[0] = start_symbol; 
  if (!fusion_sum_traverse_recursive(r, mol, start_symbol, 
                                     ordered_path, 
                                     shares, 
                                     local_seen, 1)) 
  {
    fprintf(stderr, "Error: no locant path possible for ring of size: %u\n", size); 
    return false; 
  }
  
  /* set up multi bit set */
  for (unsigned int i=0; i<size; i++) {
    if (symbol_ring_share_count(r, r->locants[i]) > 2) 
      r->multi |= (1 << i); 
  }
  
  free(local_seen); 
  free(ordered_path); 
  free(shares); 
  return true; 
}


static int lowest_fusion_locant(struct wln_ring *r, ring_t *subcycle)
{
  for (unsigned int j=0; j<r->size; j++) {
    symbol_t *atom = r->locants[j]; 
    if (ring_symbol_member(subcycle, atom)) 
      return j; 
  }
  return 0; 
}


static bool wln_write_cycle(struct wln_ring *r, 
                            graph_t *mol)
{
  if (r->hetero) 
    wln_push('T');
  else 
    wln_push('L');
  
  const unsigned int sssr_size = r->nsssr; 
  
  int *atoms_seen = (int*)malloc(sizeof(int)*sssr_size);
  for (unsigned int i=0; i<sssr_size; i++) 
    atoms_seen[i] = ring_get_size(r->sssr[i]);  
  
  for (unsigned int i=0; i<r->size; i++) {
    symbol_t *atom = r->locants[i]; 
    for (unsigned int j=0; j<sssr_size; j++) {
      ring_t *subcycle = r->sssr[j]; 
      if (ring_symbol_member(subcycle, atom)) {
        if (--atoms_seen[j] == 0) {
          int locant = lowest_fusion_locant(r, subcycle);
          if (locant != 0) {
            wln_push(' '); 
            wln_push((locant + 'A')); 
          }
          wln_push((ring_get_size(subcycle) + '0')); 
        }
      }
    }
  }  
  
  /* write the multi block and size */
  if (r->multi) {
    wln_push(' ');
    wln_push((__builtin_popcount(r->multi) + '0')); 

    unsigned int pos = 0;
    unsigned long long n = r->multi; 
    while (n) {
    if (n & 1)
      wln_push((pos + 'A'));
    n >>= 1; 
    pos++; 
    }
    wln_push(' ');
    wln_push((r->size-1 + 'A')); 
  }

  /* write any heteroatom assignments */ 
  if (r->hetero) {
    unsigned int last_pos = -1; 
    for (unsigned int i=0; i<r->size; i++) {
      symbol_t *atom = r->locants[i]; 
      if (symbol_get_num(atom) != 6) {
        if (i != last_pos + 1) {
          wln_push(' ');
          wln_push((i + 'A'));
        }
        wln_write_element_symbol(atom); 
        last_pos = i; 
      }
    } 
  }
  
  wln_push('J'); 
  free(atoms_seen); 
  return true; 
}


static int branch_recursive_write(graph_t *mol, symbol_t *atom)
{
  seen[symbol_get_id(atom)] = true; 
  wln_write_element_symbol(atom);
  
  unsigned int nbranch = 1; /* already came from one branch */ 
  const unsigned int ndegree = symbol_get_degree(atom); 

  symbol_nbor_iter(s, atom) {
    symbol_t *nbor = &(*s); 
    edge_t *edge = graph_get_edge(mol, nbor, atom);

    if (!seen[symbol_get_id(nbor)]) {
      nbranch++; 
      for (unsigned int o=1; o < edge_get_order(edge); o++)
        wln_push('U'); 

      if (symbol_is_cyclic(nbor)) {
        wln_push('-'); 
        wln_push(' '); 

        struct wln_ring *wln_ring = wln_ring_alloc(graph_num_atoms(mol));         
        wln_ring_fill_sssr(wln_ring, mol, nbor);  
        if (!wln_ring_fill_locant_path(wln_ring, mol))
          return WLN_ERROR; 
        
        for (unsigned l=0; l<wln_ring->size; l++) {
          if (wln_ring->locants[l] == nbor) {
            wln_push((l + 'A')); 
            break; 
          }
        }

        wln_write_cycle(wln_ring, mol); 
        if (locant_recursive_write(wln_ring, mol) != WLN_OK)
          return WLN_ERROR; 
        wln_ring_free(wln_ring); 
        wln_push('&'); 
      }
      else {
        if (branch_recursive_write(mol, nbor) != WLN_OK)
          return WLN_ERROR; 
    
        if (ndegree > 2) switch (wln_back()) {
          case 'Q':
          case 'E':
          case 'F':
          case 'G':
          case 'I':
          case 'Z':
            break;
          default:
            if (nbranch != ndegree)
              wln_push('&'); 
        }
      }
    }
  }  
  return WLN_OK;
}


static int locant_recursive_write(struct wln_ring *wln_ring, 
                                  graph_t *mol)
{
  
  for (unsigned int i=0; i<wln_ring->size; i++) {
    symbol_t *locant = wln_ring->locants[i]; 
    symbol_bond_iter(b, locant) {
      edge_t *edge = &(*b); 
      symbol_t *nbor = edge_get_nbor(edge, locant); 
      if (!symbol_is_cyclic(nbor) && !seen[symbol_get_id(nbor)]) {
        wln_push(' ');
        wln_push((i + 'A'));

        for (unsigned int o = 1; o < edge_get_order(edge); o++)
          wln_push('U'); 

        if (branch_recursive_write(mol, nbor) != WLN_OK)
          return WLN_ERROR; 
      }
    }
  }
  return WLN_OK; 
}


/* all sequences of 1's can be folded into their singular numbers */
static void fold_carbon_chains()
{
  char aux[4096]; 
  unsigned int naux = 0; 

  int start = -1; 
  unsigned int chain_len = 0;  
  for (unsigned int i=0; i<wln_length(); i++) {
    unsigned char ch = wln_at(i); 
    if (ch == '1') {
      if (start == -1) {
        start = i; 
        chain_len = 1; 
      }
      else 
        chain_len++; 
    }
    else {
      if (start != -1) {
        while (chain_len > 0) {
          unsigned char digit = (chain_len % 10) + '0';  // get rightmost digit
          aux[naux++] = digit; 
          chain_len = chain_len / 10;    // remove rightmost digit
        }
      }
      aux[naux++] = ch;  
      start = -1;
    }
  }

  while (chain_len > 0) {
    unsigned char digit = (chain_len % 10) + '0';  // get rightmost digit
    aux[naux++] = digit; 
    chain_len = chain_len / 10;    // remove rightmost digit
  }
  
  memcpy(wln_out, aux, naux);
  wln_length() = naux; 
  wln_terminate(); 
}


bool WriteWLN(char *buffer, unsigned int nbuffer, OBMol* mol)
{   
  wln_out = buffer; 
  wln_len = 0; 

  seen = (bool*)malloc(graph_num_atoms(mol)); 
  memset(seen, 0, graph_num_atoms(mol)); 

  char new_mol = 0; 
  bool is_cyclic = graph_num_cycles(mol) > 0; 

  if (!is_cyclic) {
    /* find an atom that has only one bond, if not cyclic this must be true */
    symbol_t *seed; 
    graph_symbol_iter(s, mol) {
      seed = &(*s); 
      /* from the seed atom, recursively build the branch */
      if (!seen[symbol_get_id(s)] && symbol_get_degree(seed) == 1) {
        if (new_mol) {
          wln_push(' ');
          wln_push('&');
        }
        
        if (branch_recursive_write(mol, seed) != WLN_OK)
          return false; 
        new_mol = 1; 
      }
    }  
  }
  else {
#ifdef USING_OPENBABEL
    /* 
     * openbabel rings have an internal id assignement, 
     * this is naughty, but stops the need for a address hash 
     */
    unsigned int nrings = 0; 
    graph_ring_iter(r, mol) {
      ring_t *ring = &(*r); 
      ring->ring_id = nrings++; 
    }
#endif

    /* get the first cycle, and create ring */
    graph_ring_iter(r, mol) {
      ring_t *ring = &(*r); 
      symbol_t *root = ring_get_first(ring, mol); 
      if (!seen[symbol_get_id(root)]) {
        seen[symbol_get_id(root)] = true; 
        struct wln_ring *wln_ring = wln_ring_alloc(graph_num_atoms(mol));         
        wln_ring_fill_sssr(wln_ring, mol, root);  
        if (!wln_ring_fill_locant_path(wln_ring, mol))
          return false; 
        wln_write_cycle(wln_ring, mol);
        if (locant_recursive_write(wln_ring, mol) != WLN_OK)
          return false; 
        wln_ring_free(wln_ring); 
      }
    }
  }

  if (wln_length() == 0)
    return false; 

  fold_carbon_chains(); 
  while (wln_back() == '&')
    wln_pop(); 
  
  fprintf(stderr, "%s\n", wln_out); 

  free(seen); 
  return true; 
}


