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

    //compute all pairs shortest paths
    begin = clock();
    vector<int> ecc;
    prop.eccentricity(g,ecc);
    end = clock();
    printf("Alg Time (eccentricity): %f\n", double(end - begin) / CLOCKS_PER_SEC);

    begin = clock();
    vector<float> ecc_dist;
    prop.eccentricity_dist(g,ecc, ecc_dist);
    end = clock();
    printf("Alg Time (freq dist eccentricity): %f\n", double(end - begin) / CLOCKS_PER_SEC);
} // main

