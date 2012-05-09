/***************************************************************************/
/*                                                                         */
/*                    Code for Minimum Node Cuts                           */
/*                                                                         */
/*  Written by:  William Cook                                              */
/*  Date:  Fall 1998                                                       */
/*                                                                         */
/*  EXPORTED FUNCTIONS                                                     */
/*                                                                         */
/*    void BW_init_pgraph (BW_pgraph *G)                                   */
/*    void BW_free_pgraph (BW_pgraph *G)                                   */
/*    int BW_build_pgraph (BW_pgraph *G, int ncount, int ecount,           */ 
/*        int *elist)                                                      */
/*    int BW_disjoint_flow_cut (BW_pgraph *G, int scount, int *slist,      */
/*        int tcount, int *tlist, int *ccount, int *cut)                   */
/*      FINDS a min cardinality node cut from S to T                       */
/*       - scount and slist specify S and tcount and tlist specify T       */
/*       - ccount returns the number of nodes in the cut, and cut returns  */
/*         an array of the indices of the nodes in the cut.                */
/*    int BW_disjoint_paths_cut (BW_pgraph *G, int scount, int *slist,     */
/*        int tcount, int *tlist, int *ccount, int *cut, int requirement)  */
/*      FINDS a min cardinality cut from S to T or determines that there   */
/*        are at least requirement number of disjoint paths from S to T    */
/*        (in the later case, ccount is set to 9).                         */
/*    void BW_internally_disjoint_paths_cut (BW_pgraph *G, BW_pnode *s,    */
/*        BW_pnode *t, int *ccount, int *cut, int requirement)             */
/*      SPECIAL case of disjoint_paths_cut where S and T are single nodes. */ 
/*                                                                         */
/*  NOTES                                                                  */
/*                                                                         */
/***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "paths.h"

static void
    insert_edge (BW_pedge *e),
    delete_edge (BW_pedge *e),
    internally_disjoint_paths (BW_pgraph *G, BW_pnode *s, BW_pnode *t,
    int requirement, int *count),
    nextpath (BW_pgraph *G, int *hit),
    push (BW_pgraph *G, BW_pnode *v, int *hit),
    retrace (BW_pgraph *G);
static int
    disjoint_flow (BW_pgraph *G, int scount, int *slist, int tcount,
    int *tlist, int *count),
    disjoint_paths (BW_pgraph *G, int scount, int *slist, int tcount,
    int *tlist, int requirement, int *count),
    internally_disjoint_flow (BW_pgraph *G, BW_pnode *s, BW_pnode *t,
    int *count),
    create_super_source_sink (BW_pgraph *G, BW_pnode *s, BW_pnode *t,
    BW_pedge **newedges, int scount, int *slist, int tcount, int *tlist);

static BW_pedge
    *checkedge (BW_pnode *n1, BW_pnode *n2);

void *CCutil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf (stderr, "Warning: 0 bytes allocated\n");
    }

    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void CCutil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }

    free (p);
}


void BW_init_pgraph (BW_pgraph *G)
{
    if (G) {
        G->nodelist  = (BW_pnode *) NULL;
        G->edgelist  = (BW_pedge *) NULL;
        G->nodestack = (BW_pnode **) NULL;
        G->ncount    = 0;
        G->ecount    = 0;
        G->HIT       = 0;
        G->MARK      = 0;
    }
}

void BW_free_pgraph (BW_pgraph *G)
{
    if (G) 
    {
#if 0
        if(G->nodelist)
        {
            free(G->nodelist);
            G->nodelist=NULL;
        }
        if(G->edgelist)
        {
            free(G->edgelist);
            G->edgelist=NULL;
        }
        if(G->nodestack)
        {
            free(G->nodestack);
            G->nodestack=NULL;
        }
#else
        CC_IFFREE (G->nodelist, BW_pnode);
        CC_IFFREE (G->edgelist, BW_pedge);
        CC_IFFREE (G->nodestack, BW_pnode *);
#endif
    }
}

int BW_build_pgraph (BW_pgraph *G, int ncount, int ecount, int *elist)
{
    int i;
    int rval = 0;

    BW_init_pgraph (G);

    G->ncount = ncount;
    G->ecount = ecount;
    G->nodelist = CC_SAFE_MALLOC (ncount, BW_pnode);
    G->edgelist = CC_SAFE_MALLOC (ecount, BW_pedge);
    if (!G->nodelist || !G->edgelist) {
        fprintf (stderr, "out of memory in BW_build_pgraph\n");
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < ecount; i++) {
        G->edgelist[i].ends[0] = G->nodelist + elist[2*i];
        G->edgelist[i].ends[1] = G->nodelist + elist[2*i+1];
    }

    for (i = 0; i <ncount; i++) {
        G->nodelist[i].adj   = (BW_pedgeptr *) NULL;
        G->nodelist[i].hit   = 0;
        G->nodelist[i].mark  = 0;
        G->nodelist[i].mark2 = 0;
    }

    for (i = 0; i < ecount; i++) {
        insert_edge (&(G->edgelist[i]));
    }

    G->HIT = G->MARK = 1;
    G->nodestack = CC_SAFE_MALLOC (G->ncount + 3, BW_pnode *);
    if (!G->nodestack) {
        fprintf (stderr, "out of memory in BW_build_pgraph\n");
        rval = 1;  goto CLEANUP;
    }

CLEANUP:

    if (rval) BW_free_pgraph (G);
    return rval;
}


static void insert_edge (BW_pedge *e)
{
    BW_pedgeptr *ep;
    BW_pnode *n0 = e->ends[0];
    BW_pnode *n1 = e->ends[1];

    ep = &(e->pends[0]);
    ep->this_one = e;
    ep->otherend = n1;
    ep->next = n0->adj;
    n0->adj = ep;

    ep = &(e->pends[1]);
    ep->this_one = e;
    ep->otherend = n0;
    ep->next = n1->adj;
    n1->adj = ep;
}

static BW_pedge *checkedge (BW_pnode *n1, BW_pnode *n2)
{
    BW_pedgeptr *ep;

    for (ep = n1->adj; ep; ep = ep->next) {
        if (ep->otherend == n2) {
            return ep->this_one;
        }
    }
    return (BW_pedge *) NULL;
}

static void delete_edge (BW_pedge *e)
{
    BW_pedgeptr *ep;
    BW_pedgeptr *prev, *next;
    int i;

    /* edge e should be in the adj lists at this point */

    for (i = 0; i < 2; i++) {
        ep = e->ends[i]->adj;
        if (ep->this_one == e) {
            e->ends[i]->adj = ep->next;
        } else {
            prev = ep;
            for (ep = ep->next; ep; ep = next) {
                next = ep->next;
                if (ep->this_one == e) {
                    prev->next = ep->next;
                    break;
                }
                prev = ep;
            }
        }
    }
}

static int create_super_source_sink (BW_pgraph *G, BW_pnode *s, BW_pnode *t,
    BW_pedge **newedges, int scount, int *slist, int tcount, int *tlist)
{
    int rval = 0;
    int i;
    BW_pnode *n;
    BW_pedge *e;

    *newedges = (BW_pedge *) NULL;

    s->hit = s->mark = s->mark2 = 0;
    t->hit = t->mark = t->mark2 = 0;
    s->adj = (BW_pedgeptr *) NULL;
    t->adj = (BW_pedgeptr *) NULL;

    for  (i = 0; i < scount; i++) {
        n = G->nodelist + slist[i];
        e = CC_SAFE_MALLOC (1, BW_pedge);
        if (!e) {
            rval = 1; goto CLEANUP;
        }
        e->ends[0] = s;
        e->ends[1] = n;
        // CSG
        fprintf(stderr,"Inserting source edge (%d,%d)\n",s,n);
        insert_edge (e);
        e->next = *newedges;
        *newedges = e;
    }
    for (i = 0; i < tcount; i++) {
        n = G->nodelist + tlist[i];
        e = CC_SAFE_MALLOC (1, BW_pedge);
        if (!e) {
            rval = 1; goto CLEANUP;
        }
        e->ends[0] = t;
        e->ends[1] = n;
        // CSG
        fprintf(stderr,"Inserting sink edge (%d,%d)\n",t,n);
        insert_edge (e);
        e->next = *newedges;
        *newedges = e;
    }

CLEANUP:

    return rval;
}

int BW_disjoint_flow_cut (BW_pgraph *G, int scount, int *slist, int tcount,
    int *tlist, int *ccount, int *cut)
{
    int i, k;
    int rval = 0;

    *ccount = 0;

    rval = disjoint_flow (G, scount, slist, tcount, tlist, &k);
    // CSG
    fprintf(stderr,"disjoint_flow set k=%d\n",k);
    if (rval) {
        fprintf (stderr, "disjoint_flow failed\n"); goto CLEANUP;
    }

    if (k == 1 && scount < 2) 
    {
        cut[(*ccount)++] = slist[0];
    } 
    else 
    {
        // CSG
        fprintf(stderr,"Computing cut\nG->MARK=%d\nG->HIT=%d\n",G->MARK,G->HIT);
        fprintf(stderr,"source->mark=%d\nsink->mark=%d\n",G->source->mark,G->sink->mark);
        for (i = 0; i < G->ncount; i++)
        {
            fprintf(stderr,"i=%d, mark=%d, hit=%d\n",i,G->nodelist[i].mark,
                G->nodelist[i].hit);
            if (G->nodelist[i].mark == G->MARK) 
            {
                cut[(*ccount)++] = i;
            }
        }
    }

CLEANUP:

    return rval;
}

static int disjoint_flow (BW_pgraph *G, int scount, int *slist, int tcount,
    int *tlist, int *count)
{
    BW_pnode s, t;
    BW_pedge *e, *enext;
    BW_pedge *newedges = (BW_pedge *) NULL;
    int rval = 0;

    *count = 0;

    rval = create_super_source_sink (G, &s, &t, &newedges, scount, slist,
        tcount, tlist);
    if (rval) {
        fprintf (stderr, "create_super_source_sink failed\n");
        goto CLEANUP;
    }

    rval = internally_disjoint_flow (G, &s, &t, count);
    if (rval) {
        fprintf (stderr, "internally_disjoint_flow failed\n"); goto CLEANUP;
    }
    fprintf(stderr,"internally_disjoint_flow setting count=%d\n",*count);

CLEANUP:

    for (e = newedges; e; e = enext) {
        enext = e->next;
        delete_edge (e);
        CC_FREE (e, BW_pedge);
    }
    return rval;
}

int BW_disjoint_paths_cut (BW_pgraph *G, int scount, int *slist, int tcount,
    int *tlist, int *ccount, int *cut, int requirement)
{
    int i, k;
    int rval = 0;

    *ccount = 0;

    rval = disjoint_paths (G, scount, slist, tcount, tlist, requirement, &k);
    if (rval) {
        fprintf (stderr, "disjoint_paths failed\n"); goto CLEANUP;
    }

    if (k < requirement) 
    {
        if (scount < 2 && k == 1) 
        {
            cut[(*ccount)++] = slist[0];
        } 
        else  
        {
            for (i = 0; i < G->ncount; i++) {
                if (G->nodelist[i].mark == G->MARK) {
                    cut[(*ccount)++] = i;
                }
            }
        }
    }

CLEANUP:

    return rval;
}

static int disjoint_paths (BW_pgraph *G, int scount, int *slist, int tcount,
    int *tlist, int requirement, int *count)
{
    BW_pnode s, t;
    BW_pedge *e, *enext;
    BW_pedge *newedges = (BW_pedge *) NULL;
    int rval = 0;

    rval = create_super_source_sink (G, &s, &t, &newedges, scount, slist,
        tcount, tlist);
    if (rval) {
        fprintf (stderr, "create_super_source_sink failed\n");
        goto CLEANUP;
    }

    internally_disjoint_paths (G, &s, &t, requirement, count);

CLEANUP:

    for (e = newedges; e; e = enext) {
        enext = e->next;
        delete_edge (e);
        CC_FREE (e, BW_pedge);
    }

    return rval;
}

static int internally_disjoint_flow (BW_pgraph *G, BW_pnode *s, BW_pnode *t,
    int *count)
{
    int k = -1, do_delete = 0;
    int hit;
    BW_pedge *e;

    G->source = s;
    G->sink = t;
    G->HIT++;

    if ((e = checkedge (G->source, G->sink)) != (BW_pedge *) NULL) 
    {
        //CSG
        fprintf(stderr,"There is an edge from source to sink?\n");
        k++;
        delete_edge (e);
        do_delete = 1;
    }
    do 
    {
        k++;
        G->top    = G->nodestack;
        *(G->top) = (BW_pnode *) NULL;
        nextpath (G, &hit);
    } while (hit);

    if (do_delete) {
        insert_edge (e);
    }
    *count = k;
    fprintf(stderr,"Set count=%d in internally_disjoint_flow\n",k);
    return 0;
}

static void internally_disjoint_paths (BW_pgraph *G, BW_pnode *s, BW_pnode *t,
    int requirement, int *count)
{
    int i = 0, do_delete = 0, hit = 0;
    BW_pedge *e;

    G->source = s;
    G->sink = t;
    G->HIT++;

    if ((e = checkedge (G->source, G->sink)) != (BW_pedge *) NULL) {
        i++;
        delete_edge (e);
        do_delete = 1;
    }

    while (i < requirement) {
        G->top    = G->nodestack;
        *(G->top) = (BW_pnode *) NULL;
        nextpath (G, &hit);
        if (!hit) {
            if (do_delete) insert_edge (e);
            *count = i;
            return;
        } else {
            i++;
        }
    }
    if (do_delete) insert_edge (e);
    *count = i;
}

void BW_internally_disjoint_paths_cut (BW_pgraph *G, BW_pnode *s,
    BW_pnode *t, int *ccount, int *cut, int requirement)
{
    int i, k;

    *ccount = 0;
    internally_disjoint_paths (G, s, t, requirement, &k);

    if (k < requirement) {
        for (i = 0; i < G->ncount; i++) {
            if (G->nodelist[i].mark == G->MARK) {
                cut[(*ccount)++] = i;
            }
        }
    }
}

static void nextpath (BW_pgraph *G, int *hit)
{
    BW_pnode *v;
    int pushhit;

    G->MARK++;
    G->source->mark = G->MARK;
    *(++(G->top)) = G->source;
    while (*(G->top) != (BW_pnode *) NULL) {
        v = *((G->top)--);
        push (G, v, &pushhit);
        if (pushhit) {
            retrace (G);
            *hit = 1;
            return;
        }
    }
    *hit = 0;
}

static void push (BW_pgraph *G, BW_pnode *v, int *hit)
{
    BW_pedgeptr *ep;
    BW_pnode *w, *u;

    for (ep = v->adj; ep; ep = ep->next) {
        w = ep->otherend;
        if (w == G->sink) {
            w->pred = v;
            *hit = 1;
            return;
        }
        if (w->mark < G->MARK) 
        {
            if (w->hit < G->HIT) 
            {
                w->mark = G->MARK;
                w->pred = v;
                *(++(G->top)) = w;
            } else if (w->mark2 < G->MARK && (v->hit < G->HIT ||
                (v->hitpred != w && w->hitpred != v))) {
                    w->pred2 = v;
                    w->mark2 = G->MARK;
                    u = w->hitpred;
                    if (u->mark < G->MARK) {
                        u->mark = G->MARK;
                        u->pred = w;
                        *(++(G->top)) = u;
                    }
            } else if (v->hit == G->HIT && v->hitpred == w) {
                w->mark = G->MARK;
                w->pred = v;
                *(++(G->top)) = w;
            }
        }
    }
    *hit = 0;
}

static void retrace (BW_pgraph *G)
{
    BW_pnode *p, *q;

    p = G->sink;
    q = p->pred;
    while (q != G->source) {
        if (q->hit < G->HIT) {
            q->hit = G->HIT;
            p->hitpred = q;
            p = q;
            q = q->pred;
        } else {
            p->hitpred = q;
            p = q;
            q = q->pred;
            while (q->mark2 != G->MARK) {
                q->hit = 0;
                p = q;
                q = q->pred;
            }
            p = q;
            q = q->pred2;
        }
    }
    p->hitpred = q;
}
