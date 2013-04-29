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
#include "WeightedMutableGraph.h"
#include "GraphException.h" 
#include "NewGraphReader.h"

using namespace std;

void usage(const char *s)
{
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

int main(int argc, char **argv)
{
    // Check for a cry for help
    if(argc==1 || (argc==2 && strcmp(argv[1],"-h")==0) || (argc==2 && strcmp(argv[1],"--help")==0) 
        || (argc==2 && strcmp(argv[1],"--h")==0) )
    {
        usage(argv[0]);
        exit(-1);
    }

    Graph::MutableGraph *g=NULL;
    int seed;
    
    seed=Graph::rand_int(0,0xffffff);
    // Spin the RNG seed times
    for(int ss=0;ss<seed;ss++)
        Graph::lcgrand(0);

    
    Graph::GraphCreatorFile *gcf;
    Graph::GraphWriter *writer;
    Graph::GraphReaderWriterFactory rwf;
    Graph::GraphProperties prop;
    Graph::NewGraphReader ngr;

    g = ngr.read_graph<Graph::WeightedMutableGraph>(argv[1], "MeTiS");

    // mucking with stuff
    //
    Graph::MutableGraph *mg;
    mg = new Graph::MutableGraph();
    //fprintf(stderr,"mg number of nodes (before): %d\n", mg->get_num_nodes());
    //fprintf(stderr,"mg before: 0x%x\n", mg);
    ngr.new_readgraph(mg, argv[1], "MeTiS");
    fprintf(stderr,"mg after: 0x%x\n", mg);
    fprintf(stderr,"mg number of nodes (after) : %d\n", mg->get_num_nodes());


    // if we don't get rid of duplicate edges, bad things happen
    // when trying to output the graph
    //prop.make_simple(g);

    fprintf(stderr, "edges read in: %d nodes read in: %d\n", g->get_num_edges(), g->get_num_nodes());

    //reader->read_graph(g, argv[1]);
    gcf = new Graph::GraphCreatorFile(argv[1], "METIS");
   
    writer = rwf.create_writer("DIMACS", "t");

    writer->set_out_file_name(argv[2]);
    writer->write_graph(g);

    return 0;
}



