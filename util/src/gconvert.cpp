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

#if !WIN32 && !CYGWIN
  #include "orbconfig.h"
  #include "orbtimer.h"
using namespace std;

void print_time(string prefix, ORB_t start, ORB_t end){
    cout << prefix + ": " << ORB_seconds(end, start) << "\n";
}

#endif

void usage(const char *s){
    fprintf(stderr,"Usage: %s intype outtype infile outfile\n",s);
}

int main(int argc, char **argv){
    #if !WIN32 & !CYGWIN
    ORB_t t1, t2;

    char *intype, *outtype, *infile, *outfile;

    // Check for a cry for help
    if((argc == 1) || ((argc == 2) && (strcmp(argv[1],"-h") == 0)) || ((argc == 2) && (strcmp(argv[1],"--help") == 0))
       || ((argc == 2) && (strcmp(argv[1],"--h") == 0) ) ){
        usage(argv[0]);
        exit(0);
    }

    if(argc != 5){
        usage(argv[0]);
        exit(1);
    }

    intype = argv[1];
    outtype = argv[2];
    infile = argv[3];
    outfile = argv[4];

    Graph::Graph *g;
    int seed = 0;

    cout << "calibrating timers\n";

    ORB_calibrate();
    if(!seed){
        // Set the seed to a rand int in 0,2^24
        seed = Graph::rand_int(0,0xffffff);
    }
    // Spin the RNG seed times
    for(int ss = 0; ss < seed; ss++){
        Graph::lcgrand(0);
    }

    Graph::GraphCreatorFile *gcf;

    Graph::GraphProperties prop;
    Graph::GraphReader ngr;
    Graph::GraphWriter writer;

    g = new Graph::Graph();

    cout << "Input type : " << intype << endl;
    cout << "Output type: " << outtype << endl;
    cout << "Reading graph" << endl;
    ORB_read(t1);
    ngr.read_graph(g, infile, intype, false);
    ORB_read(t2);
    print_time("Time(read_graph)", t1, t2);

    // if we don't get rid of duplicate edges, bad things happen
    // when trying to output the graph
    //prop.make_simple(g);

    fprintf(stderr, "edges read in: %d nodes read in: %d\n", g->get_num_edges(), g->get_num_nodes());

    cout << "Writing graph\n";
    ORB_read(t1);
    writer.write_graph(g, outfile, outtype);
    ORB_read(t2);
    print_time("Time(write_graph)", t1, t2);

    return 0;
    #else // if !WIN32 & !CYGWIN
    fprintf(stderr,"Can't build under Cygwin or Windows\n");
    return 1;
    #endif // if !WIN32 & !CYGWIN
} // main

