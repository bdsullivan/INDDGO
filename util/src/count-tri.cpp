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

#include "GraphDecomposition.h"
#include "GraphProperties.h"
#include "Log.h"
#include "VertexWeightedGraph.h"
#include "GraphException.h"
#include "NewGraphReader.h"
#include "NewGraphWriter.h"

using namespace std;

void usage(const char *s){
    fprintf(stderr,"Usage: %s [options]\n",s);
    fprintf(stderr,
            "\t -t [#] sets number of graphs to generate (t)\n"
            "\t -n [#] Sets the number of vertices (n)\n"
            "\t -k [#] sets value for tree-width (k)\n"
            "\t -p [#] sets percent of edges to keep from complete k-tree (p)\n"
            "\t -s seed : uses the given RNG seed\n"
            "\t -fn [name] sets the prefix for the filenames written out\n"
            "\t -r : randomizes the vertex labels before writing to file.\n"
            "\t -scotch : writes a .scotch file in Scotch graph format as well\n"
            );
}

int main(int argc, char **argv){
    // Check for a cry for help
    if((argc == 1) || ((argc == 2) && (strcmp(argv[1],"-h") == 0)) || ((argc == 2) && (strcmp(argv[1],"--help") == 0))
       || ((argc == 2) && (strcmp(argv[1],"--h") == 0) ) ){
        usage(argv[0]);
        exit(-1);
    }

    Graph::Graph *g;
    int seed = 0;


    Graph::GraphProperties prop;
    Graph::NewGraphReader ngr;

    g = new Graph::Graph();

    //fprintf(stderr,"g before: 0x%x\n", g);
    ngr.read_graph(g, argv[1], "Edge", false);
    //fprintf(stderr,"g after: 0x%x\n", g);
    //
    fprintf(stderr, "Read %d vertices and %d edges\n", g->get_num_nodes(), g->get_num_edges());

    vector<long int> triangles(g->get_num_nodes(), 0);

    //prop.all_triangles(g, triangles, 5);
    prop.all_triangles_edge_listing(g, triangles);

    for(int i=0; i<g->get_num_nodes(); i++){
        fprintf(stderr, "vertex: %d: %ld\n", i, triangles[i]);
    }


    return 0;
} // main

