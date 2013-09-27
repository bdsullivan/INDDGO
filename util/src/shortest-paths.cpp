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

#include "GraphReader.h"
#include "GraphProperties.h"
#include <ctime>

using namespace std;

void usage(const char *s){
    fprintf(stderr, "Usage: %s filename\n", s);
}

int main(int argc, char **argv){
    //usage check
    if((argc == 1) ||
       ((argc == 2) && (strcmp(argv[1],"-h") == 0)) ||
       ((argc == 2) && (strcmp(argv[1], "--help") == 0))){
        usage(argv[0]);
        exit(-1);
    }

    //check if too many args?

    Graph::Graph *g;

    clock_t begin, end;

    Graph::GraphProperties prop;
    Graph::GraphReader ngr;

    //create the graph object
    g = new Graph::Graph();

    //read the graph from the filename, assume it is an edgelist
    ngr.read_graph(g, argv[1], "Edge", false);
    //ngr.read_graph(g, argv[1], "ADJLIST", false);

    printf("Read %d vertices and %d edges\n", g->get_num_nodes(), g->get_num_edges());

    printf("Simplifying graph\n");
    begin = clock();
    prop.make_simple(g);    //remove self loops and duplicate edges
    end = clock();
    printf("Time: %f\n", double(end - begin) / CLOCKS_PER_SEC);

    //compute shortest paths from source node
    begin = clock();
    //for each node run it
    int source = 8;
    vector<int> onePaths;

    printf("\nRunning Dijsktra for source %d\n",source);
    prop.paths_dijkstra_single(g, onePaths, source);
    end = clock();

    int i;
    for(i = 0; i < onePaths.size(); i++){
        printf("%d: %d\n", i, onePaths[i]);
    }
    printf("Alg Time (single): %f\n", double(end - begin) / CLOCKS_PER_SEC);

    //compute all pairs shortest paths
    begin = clock();
//    vector< vector<int> > allPaths;

    //   printf("\nRunning Dijsktra for all pairs \n");
    //  prop.paths_dijkstra_all(g,allPaths); //alternate: vector< vector<int> > pt2 =  g->get_shortest_path_dist_ref();
    // end = clock();
    //printf("Alg Time (all): %f\n", double(end - begin) / CLOCKS_PER_SEC);
} // main

