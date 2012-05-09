/*
  This file is part of INDDGO.

  Copyright (C) 2012, Oak Ridge National Laboratory 

  This product includes software produced by UT-Battelle, LLC under Contract No. 
  DE-AC05-00OR22725 with the Department of Energy. 

  This program is free software; you can redistribute it and/or modify
  it under the terms of the New BSD 3-clause software license (LICENSE). 
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  LICENSE for more details.

  For more information please contact the INDDGO developers at: 
  inddgo-info@googlegroups.com

*/

#ifndef __PATHS_H
#define __PATHS_H

#define CC_SAFE_MALLOC(nnum,type) (type *) CCutil_allocrus (((size_t) (nnum)) * sizeof (type))

#define CC_FREE(object,type) { CCutil_freerus ((void *) (object)); object = (type *) NULL; }

#define CC_IFFREE(object,type) { if ((object)) CC_FREE ((object),type);}


/*
#define CC_SAFE_MALLOC(nnum,type)					
    (type *) malloc (((size_t) (nnum)) * sizeof (type))

#define CC_FREE(object,type) {						
	free ((void *) (object));					
	object = (type *) NULL;						
}

#define CC_IFFREE(object,type) {                                           
    if ((object)) CC_FREE ((object),type);                                 
}
*/
typedef struct BW_pedgeptr {
    struct BW_pedge    *this_one;
    struct BW_pnode    *otherend;
    struct BW_pedgeptr *next;
} BW_pedgeptr;
 
typedef struct BW_pnode {
    struct BW_pedgeptr *adj;
    struct BW_pnode    *pred;
    struct BW_pnode    *pred2;
    struct BW_pnode    *hitpred;
    int                 hit;
    int                 mark;
    int                 mark2;
} BW_pnode;
 
typedef struct BW_pedge {
    BW_pedgeptr         pends[2];
    struct BW_pnode    *ends[2];
    struct BW_pedge    *next;
    int                 mark;
} BW_pedge;
 
typedef struct BW_pgraph {
    struct BW_pnode    *nodelist;
    struct BW_pedge    *edgelist;
    struct BW_pnode   **nodestack;
    struct BW_pnode   **top;
    struct BW_pnode    *source;
    struct BW_pnode    *sink;
    int                 ncount;
    int                 ecount;
    int                 HIT;
    int                 MARK;
} BW_pgraph;


void
    BW_init_pgraph (BW_pgraph *G),
    BW_free_pgraph (BW_pgraph *G),
    BW_internally_disjoint_paths_cut (BW_pgraph *G, BW_pnode *s,
        BW_pnode *t, int *ccount, int *cut, int requirement);
int
    BW_build_pgraph (BW_pgraph *G, int ncount, int ecount, int *elist),
    BW_disjoint_flow_cut (BW_pgraph *G, int scount, int *slist, int tcount,
        int *tlist, int *ccount, int *cut),
    BW_disjoint_paths_cut (BW_pgraph *G, int scount, int *slist, int tcount,
        int *tlist, int *ccount, int *cut, int requirement);


int 
    BW_connected_test (int ncount, int ecount, int *elist, int *yes_no),
    BW_biconnected_test (int ncount, int ecount, int *elist, int *yes_no);


#endif  /* __PATHS_H */
