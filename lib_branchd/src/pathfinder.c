#include <stdio.h>
#include "paths.h"

#define NUM_EDGES 8
#define NUM_NODES 8

int main(int argc, char **argv)
{

    // Declare a BW_pgraph
    BW_pgraph G;
    /*int elist[2*NUM_EDGES]={
        1, 2, 
        1, 6,
        1, 9,
        1, 14,
        2, 3,
        2, 5,
        2, 7,
        2, 11,
        2, 13,
        3, 5,
        3, 10,
        3, 15,
        4, 8,
        4, 12,
        5, 6,
        5, 11,
        5, 15,
        6, 7,
        6, 8,
        6, 11
    };*/
    int elist[2*NUM_EDGES]={
        0,1,
        2,3,
        1,4,
        3,5,
        6,0,
        6,2,
        4,7,
        5,7
    };

  
       /* int elist[2*NUM_EDGES]={ 
            1, 2, 
            1,3,
            1,4,
            1,5,
            2,3,
            2,4,
            2,5,
            3,4,
            3,6,
            4,6,
            4,7,
            4,8,
            5,7,
            5,8,
            6,7,
            6,8,
            7,8};*/


    int i, cutsize, source[4], sink[4], cut[2*NUM_EDGES];

    // Make the edge list zero based
    //for(i=0;i<2*NUM_EDGES;i++)
    //    elist[i]--;
    source[0]=1; source[1]=2; source[2]=3; source[3]=4;
    sink[0]=3; sink[1]=4; sink[2]=5; sink[3]=6;

    BW_build_pgraph (&G, NUM_NODES, NUM_EDGES, elist);

    fprintf(stderr,"Graph built\n");
    //BW_disjoint_paths_cut(&G,4,source,4,sink,&cutsize,cut,0);
    BW_internally_disjoint_paths_cut (&G, G.nodelist+6, G.nodelist+7,&cutsize,cut,10);
        //BW_disjoint_flow_cut(&G,2,source,2,sink,&cutsize,cut);
    /*      FINDS a min cardinality node cut from S to T                       */
    /*       - scount and slist specify S and tcount and tlist specify T       */
    /*       - ccount returns the number of nodes in the cut, and cut returns  */
    /*         an array of the indices of the nodes in the cut.                */
    fprintf(stderr,"cutsize=%d\n",cutsize);
    for(i=0;i<cutsize;i++)
        fprintf(stderr,"Cut node %d=%d\n",i,cut[i]);

    BW_free_pgraph(&G);

    return 1;

}