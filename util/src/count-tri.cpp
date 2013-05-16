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
#include "GraphReader.h"
#include "GraphWriter.h"
#include <numeric>
#include <ctime>

using namespace std;

void usage(const char *s){
    fprintf(stderr,"Usage: %s filename\n",s);
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
    int sum;

    clock_t begin, end;

    Graph::GraphProperties prop;
    Graph::GraphReader ngr;

    // Create the graph object
    g = new Graph::Graph();

    // read the graph from the filename, assume it is an edgelist
    ngr.read_graph(g, argv[1], "Edge", false);
    printf("Read %d vertices and %d edges\n", g->get_num_nodes(), g->get_num_edges());

    printf("Simplifying graph\n");
    begin = clock();
    prop.make_simple(g);
    end = clock();
    printf("Time: %lf\n", double(end - begin) / CLOCKS_PER_SEC);
    printf("After simplification: %d vertices and %d edges\n", g->get_num_nodes(), g->get_num_edges());

    //Now, do our calculations
    vector<long int> triangles(g->get_num_nodes(), 0);

    printf("Calculating triangles using compact-forward method\n");
    begin = clock();
    prop.all_triangles_compact_forward(g, triangles);
    end = clock();
    sum = std::accumulate(triangles.begin(), triangles.end(), 0);
    printf("Total triangles (compact-forward): %d\n", sum);
    printf("Time: %lf\n", double(end - begin) / CLOCKS_PER_SEC);

    triangles.assign(g->get_num_nodes(), 0);
    printf("Calculating triangles using edge-listing method\n");
    begin = clock();
    prop.all_triangles_edge_listing(g, triangles);
    end = clock();
    sum = std::accumulate(triangles.begin(), triangles.end(), 0);
    printf("Total triangles (edge-listing): %d (%d)\n", sum / 3, sum);
    printf("Time: %lf\n", double(end - begin) / CLOCKS_PER_SEC);

    double g_cc, a_cc;
    vector<double> l_cc;
    prop.clustering_coefficients(g, g_cc, a_cc, l_cc);

    printf("Local CCs:");
    int i;
    for(i = 0; i < g->get_num_nodes(); i++){
        printf(" %d:%lf",i,l_cc[i]);
    }
    printf("\n");

    printf("Global cc: %lf\nAvg cc: %lf", g_cc, a_cc);

    delete g;

    return 0;
} // main

