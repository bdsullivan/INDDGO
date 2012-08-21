#ifndef SUPERWRAP_H_
#define SUPERWRAP_H_

#include <cholmod.h>


/* On normal machines, SuiteSparse_long is just a "long".  On Windows, though,
 * it's something else.  See SuiteSparse_config for details. */
#define Long SuiteSparse_long

/* SuperTree: a supernodal tree decomposition */
typedef struct SuperTree_struct
{
    Long n ;        /* number of nodes in the graph */
    Long ns ;       /* number of supernodes in the tree */
    Long maxsuper ; /* size of largest supernode (# of nonzero row indices) */
    Long *Sparent ; /* size ns.  Sparent[i] = j if supernode j is the
                        parent of supernode i */
    Long *Perm ;    /* size n, permutation used to create the tree */
    Long *Sp ;      /* size ns+1 */
    Long *Si ;      /* size Sp [ns].  The nodes in the kth supernode are
                        given by the list Si [Sp [k] ... Sp [k+1]-1]. */
} SuperTree ;

SuperTree *SuperWrap        /* returns a tree decomposition or NULL if error */
(
    /* inputs, not modified on output */
    Long n,         /* number of nodes in the input graph */
    Long *Ap,       /* column pointers for input graph, size n+1 */
    Long *Ai,       /* row indices for input graph, size nz = Ap [n]. */
    Long *UserPerm, /* user-provided permutation; ignored if NULL */

    /* input/output: esoteric parameters, status, and workspace */
    cholmod_common *Common
) ;

void SuperFree              /* frees a tree decomposition */
(
    SuperTree **SHandle,    /* SuperTree to free */
    cholmod_common *Common
) ;

#endif
