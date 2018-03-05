#ifndef _PAIR_MAT_HH
#define _PAIR_MAT_HH
/**
 * From pair_mat.h in Vienna RNA package 2.0.7.
 */

#define NBASES (8)
#define NBPAIRS (7)

// static int BP_pair[5][5]=
// /* _  A  C  G  U*/
// {{ 0, 0, 0, 0, 0},
//  { 0, 0, 0, 0, 5},
//  { 0, 0, 0, 1, 0},
//  { 0, 0, 2, 0, 3},
//  { 0, 6, 0, 4, 0}};

static int BP_pair[NBASES][NBASES]=
/* _  A  C  G  U  X  K  I */
{{ 0, 0, 0, 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5, 0, 0, 5},
 { 0, 0, 0, 1, 0, 0, 0, 0},
 { 0, 0, 2, 0, 3, 0, 0, 0},
 { 0, 6, 0, 4, 0, 0, 0, 6},
 { 0, 0, 0, 0, 0, 0, 2, 0},
 { 0, 0, 0, 0, 0, 1, 0, 0},
 { 0, 6, 0, 0, 5, 0, 0, 0}};

static int rtype[7] = {0, 2, 1, 4, 3, 6, 5};

#endif
