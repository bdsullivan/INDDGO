/* ========================================================================== */
/* === SuperWrap ============================================================ */
/* ========================================================================== */

/* The SuperWrap function is a simple wrapper around the symbolic supernodal
    factorization of CHOLMOD, which constructs a supernodal elimination tree
    and the nonzero pattern of the supernodes themselves.  Strict supernodal
    amalgamation is enforced (supernodes are amalgamated only if no explicit
    zeros are introduced into the pattern of L).

    To use these functions, you must first call cholmod_l_start, and when
    you are finished using these routines, you must call cholmod_l_end.
    For example:

    cholmod_common Common ;
    SuperTree *S ;
    cholmod_l_start (&Common) ;                 // creates Common for CHOLMOD
    S = SuperWrap (n, Ap, Ai, NULL, &Common) ;  // allocates and creates S
    // use S here ...
    SuperFree (&S, &Common) ;                   // frees S
    cholmod_l_finish (&Common) ;                   // frees contents of Common

    See also the superwrap mexFunction, to use this function in MATLAB.

    You can modify the permutation used to construct the supernodal analsys
    by modifying Common prior to calling SuperWrap.  For each of the following
    options, Common.postorder = 1 is the default, which allows CHOLMOD to
    postorder the etree before finding supernodes.  This decreases the number
    of supernodes while keeping the nonzero pattern of L unchanged.

        To use the default permutation, use Common as-is.  The user-provided
        permutation will be tried (if non-NULL), along with AMD and METIS
        (if METIS is available), and the best ordering obtained will be used.
        'Best' means the one with the fewest nonzeros in L.  Postordering is
        enabled by default.

        To use a natural permutation (with no postordering), use the following.
        In this case, the matrix A will not be permuted at all.
            Common.nmethods = 1 ;
            Common.method [0].ordering = CHOLMOD_NATURAL ;
            Common.postorder = 0 ;
            S = SuperWrap (n, Ap, Ai, NULL, &Common) ;

        To use a natural permutation with postordering:
            Common.nmethods = 1 ;
            Common.method [0].ordering = CHOLMOD_NATURAL ;
            S = SuperWrap (n, Ap, Ai, NULL, &Common) ;

        To use the user-provided permutation only (with postordering):
            Common.nmethods = 1 ;
            Common.method [0].ordering = CHOLMOD_GIVEN ;
            S = SuperWrap (n, Ap, Ai, UserPerm, &Common) ;

        To try 8 different orderings (including the user-provided one, if
        provided) and pick the best one:
            Common.nmethods = 8 ;
            S = SuperWrap (n, Ap, Ai, UserPerm, &Common) ;

    If UserPerm is NULL, it is assumed to be the identity permutation.

    The function returns a pointer to single struct of type SuperTree.
    See SuperWrap.h for a description.
*/

#include "SuperWrap.h"
#define TRUE 1
#define FALSE 0
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* -------------------------------------------------------------------------- */
/* SuperTree : compute a supernodal tree, and its supernodes */
/* -------------------------------------------------------------------------- */

/* The input graph is provided in compressed sparse column format, where the
   adjacency list of node j is given by Ai [Ap [j] ... Ap [j+1]-1].  Only the
   strictly lower triangular part of A is used.  The upper triangular part may
   be present, but it is ignored and assumed to be the transpose of the lower
   triangular part.  In graph terminology, the adj(j) set is given by the list
   Ai [Ap [j] ... Ap [j+1]-1], where any edge to a node i < j is ignored in
   this list (it is assumed that such j will be in adj(i)).  Ap [0] must be
   zero.

   It's OK to provide the entire adjacency list for all the nodes; this
   function simply ensures the input graph is undirected by considering edge
   (i,j) only in the adjacency of node max(i,j).  The adjacency lists need not
   be in sorted order, but duplicates must not appear.  Self-edges are ignored.
 */ 

SuperTree *SuperWrap        /* returns a tree decomposition or NULL if error */
(
    /* inputs, not modified on output */
    Long n,         /* number of nodes in the input graph */
    Long *Ap,       /* column pointers for input graph, size n+1 */
    Long *Ai,       /* row indices for input graph, size nz = Ap [n]. */
    Long *UserPerm, /* user-provided permutation; ignored if NULL */

    /* input/output: esoteric parameters, status, and workspace */
    cholmod_common *Common
)
{
    SuperTree *S ;
    cholmod_factor *L ;
    cholmod_sparse *A, Awrapper ;
    Long p, nz, s, k, j, parent, *SuperMap, *Lsuper ;
    int ssave ;
    size_t nsave [3] ;
    double zsave [3] ;

    /* ---------------------------------------------------------------------- */
    /* check inputs and allocate the SuperTree */
    /* ---------------------------------------------------------------------- */

    if (n < 0 || Ap == NULL || Ai == NULL)
    {
        return (NULL) ;
    }
    if (Ap [n] < 0)
    {
        return (NULL) ;
    }
    S = cholmod_l_malloc (1, sizeof (SuperTree), Common) ;
    if (S == NULL)
    {
        return (NULL) ;
    }
    S->n = n ;

    /* ---------------------------------------------------------------------- */
    /* allocate a shallow CHOLMOD copy for the user's input matrix */
    /* ---------------------------------------------------------------------- */

    A = &Awrapper ;
    A->nrow = n ;           /* A is n-by-n */
    A->ncol = n ;
    A->p = Ap ;             /* column pointers */
    A->i = Ai ;             /* row indices */
    A->nzmax = Ap [n] ;     /* number of entries */
    A->packed = TRUE ;      /* no gaps between the columns */
    A->sorted = FALSE ;     /* row indices of columns need not be in order */
    A->nz = NULL ;
    A->itype = CHOLMOD_LONG ;
    A->dtype = CHOLMOD_DOUBLE ;
    A->xtype = CHOLMOD_PATTERN ;
    A->stype = -1 ;         /* only lower triangular part of A is used */

    /* ---------------------------------------------------------------------- */
    /* supernodal symbolic analysis */
    /* ---------------------------------------------------------------------- */

    /* save the previous settings in Common */
    ssave = Common->supernodal ;
    nsave [0] = Common->nrelax [0] = 0 ;
    nsave [1] = Common->nrelax [1] = 0 ;
    nsave [2] = Common->nrelax [2] = 0 ;
    zsave [0] = Common->zrelax [0] = 0 ;
    zsave [1] = Common->zrelax [1] = 0 ;
    zsave [2] = Common->zrelax [2] = 0 ;

    /* require strict supernodes:  fundamental supernodes, followed by
       amalgamation that does not introduce any zeros into any supernode */
    Common->supernodal = CHOLMOD_SUPERNODAL ;
    Common->nrelax [0] = 0 ;
    Common->nrelax [1] = 0 ;
    Common->nrelax [2] = 0 ;
    Common->zrelax [0] = 0 ;
    Common->zrelax [1] = 0 ;
    Common->zrelax [2] = 0 ;

    L = cholmod_l_analyze_p (A, UserPerm, (Long *) NULL, (size_t) 0, Common) ;

    /* restore the previous settings */
    Common->supernodal = ssave ;
    Common->nrelax [0] = nsave [0] ;
    Common->nrelax [1] = nsave [1] ;
    Common->nrelax [2] = nsave [2] ;
    Common->zrelax [0] = zsave [0] ;
    Common->zrelax [1] = zsave [1] ;
    Common->zrelax [2] = zsave [2] ;

    if (L == NULL)
    {
        /* out of memory */
        cholmod_l_free (1, sizeof (SuperTree), S, Common) ;
        return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* extract the contents of the supernodal symbolic factorization */
    /* ---------------------------------------------------------------------- */

    S->ns = L->nsuper ;
    S->Perm = L->Perm ;
    S->Sp = L->pi ;
    S->Si = L->s ;
    L->Perm = NULL ;
    L->pi = NULL ;
    L->s = NULL ;

    /* ---------------------------------------------------------------------- */
    /* reconstruct the supernodal tree */
    /* ---------------------------------------------------------------------- */

    /* Sparent was computed in cholmod_*_super_analyze but then discarded. */
    S->Sparent = cholmod_l_malloc (S->ns, sizeof (Long), Common) ;
    if (S->Sparent == NULL)
    {
        SuperFree (&S, Common) ;
    }

    SuperMap = L->ColCount ;        /* use L->ColCount as size-n workspace */

    /* find the mapping of nodes to supernodes */
    Lsuper = (Long *) L->super ;
    for (s = 0 ; s < S->ns ; s++)
    {
        /* Lsuper [s] is 1st col in supernode s */
        for (k = Lsuper [s] ; k < Lsuper [s+1] ; k++)
        {
            /* SuperMap [k] = s if column k is contained in supernode s */
            SuperMap [k] = s ;
        }
    }

    /* find the largest supernode */
    S->maxsuper = 0 ;
    for (s = 0 ; s < S->ns ; s++)
    {
        S->maxsuper = MAX (S->maxsuper, S->Sp [s+1] - S->Sp [s]) ;
    }

    /* construct the supernodal tree */
    for (s = 0 ; s < S->ns ; s++)
    {
        /* Supernode s contains columns Lsuper[s] ... Lsuper[s+1]-1.  The
           pattern of supernode s is in S->Si [S->Sp [s] ... S->Sp [s+1]-1].
           The pattern is sorted, so Si [Sp [s] + Lsuper [s+1] - Lsuper [s]]
           is the row index of the first nonzero below the diagonal block of
           supernode s (if present).  This row index i corresponds to a column
           contained by the supernodal parent of supernode s.  It can then be
           mapped to the supernode that contains it via SuperMap [i].
        */
        p = S->Sp [s] + (Lsuper [s+1] - Lsuper [s]) ;
        S->Sparent [s] = (p < S->Sp [s+1]) ? SuperMap [S->Si [p]] : (-1) ;
    }

    /* ---------------------------------------------------------------------- */
    /* permute the row indices back to the original */
    /* ---------------------------------------------------------------------- */

    for (p = 0 ; p < S->Sp [S->ns] ; p++)
    {
        S->Si [p] = S->Perm [S->Si [p]] ;
    }

    cholmod_l_free_factor (&L, Common) ;

    /* ---------------------------------------------------------------------- */
    /* return the SuperTree */
    /* ---------------------------------------------------------------------- */

    return (S) ;
}

/* -------------------------------------------------------------------------- */
/* SuperFree : free a SuperTree created by SuperWrap */
/* -------------------------------------------------------------------------- */

void SuperFree              /* frees a tree decomposition */
(
    SuperTree **SHandle,    /* SuperTree to free */
    cholmod_common *Common
)
{
    Long n, ns, nz ;
    SuperTree *S ;
    if (SHandle == NULL) return ;       /* nothing to do */
    S = *SHandle ;
    n = S->n ;
    ns = S->ns ;
    nz = S->Sp [ns] ;
    cholmod_l_free (n, sizeof (Long), S->Perm, Common) ;
    cholmod_l_free (ns, sizeof (Long), S->Sparent, Common) ;
    cholmod_l_free (ns+1, sizeof (Long), S->Sp, Common) ;
    cholmod_l_free (nz, sizeof (Long), S->Si, Common) ;
    cholmod_l_free (1, sizeof (SuperTree), S, Common) ;
    SHandle = NULL ;
}
