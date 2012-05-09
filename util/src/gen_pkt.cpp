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
#include "Log.h"
#include "WeightedMutableGraph.h"
#include "GraphException.h" 

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
    // This is a utility to generate partial k-trees with specified parameters for use in other applications.

    // Check for a cry for help
    if(argc==1 || (argc==2 && strcmp(argv[1],"-h")==0) || (argc==2 && strcmp(argv[1],"--help")==0) 
        || (argc==2 && strcmp(argv[1],"--h")==0) )
    {
        usage(argv[0]);
        exit(-1);
    }

    int i;
    int t = 1;
    int ktree_n=-1, ktree_k=-1, ktree_p=-1;
    int seed=0;
	Graph::WeightedMutableGraph *G=NULL, *H=NULL; 
    time_t start, stop;
    char filename[100];
    char *prefix= (char *)"pkt";
    
    bool random = false, write_scotch=false;
    char lfname[100];
    char efname[100];
    sprintf(lfname, "genpkt.log");
    sprintf(efname, "genpkt.log");
    //0: debug, 5: critical
    LOG_INIT(lfname, efname, 0);

    // Parse arguments
    for(i=0;i<argc;i++)
    {
        if(strcmp(argv[i],"-t")==0)
        {
            t=atoi(argv[i+1]);
        }
        if(strcmp(argv[i],"-k")==0)
        {
            ktree_k=atoi(argv[i+1]);
        }
        if(strcmp(argv[i],"-n")==0)
        {
            ktree_n=atoi(argv[i+1]);
        }
        if(strcmp(argv[i],"-p")==0)
        {
            ktree_p=atoi(argv[i+1]);
        }
        if(strcmp(argv[i],"-s")==0)
        {
            seed=atoi(argv[i+1]);
        }
        if(strcmp(argv[i],"-fn")==0)
            prefix = argv[i+1];
        if(strcmp(argv[i],"-r")==0)
            random = true;
        if(strcmp(argv[i],"-scotch")==0)
            write_scotch = true;
    }

    if(!seed)
        // Set the seed to a rand int in 0,2^24
        seed=Graph::rand_int(0,0xffffff);
    // Spin the RNG seed times
    for(int ss=0;ss<seed;ss++)
        Graph::lcgrand(0);

    //Make sure we set all the parameters to viable values. 
    if(ktree_k < 0 || ktree_n < 0 || ktree_p < 0)
	{
	    fatal_error("Failed to input viable parameters for n, k, p (%d,%d,%d)\n",
            ktree_k, ktree_n, ktree_p);
	}
    
    Graph::GraphCreatorFile creator;
    Graph::GraphWriter *writer;
    Graph::GraphReaderWriterFactory rwf;
    Graph::GraphProperties prop;

    writer = rwf.create_writer("DIMACS", "t");

    DEBUG("Graph generation loop\n");
    DEBUG("n : %d k: %d\n", ktree_n, ktree_p);
    for(i = 0; i < t; i++)
    {
	    // Create the Ktree
	    //H = new Graph::Graph(ktree_n, ktree_k);
        H = creator.initialize_ktree(ktree_n, ktree_k);
	    // Generate an n-node Graph G that is a partial k-tree derived from the k-tree
	    start=clock();
	    G= creator.create_random_edge_subgraph(H, ktree_p);
	    stop=clock();
	    fprintf(stderr,"Generated random partial k-tree (%d,%d) in %f sec.\n",
		   ktree_n, ktree_k, ((double)(stop-start))/CLOCKS_PER_SEC);
	    //write it to a file
	    sprintf(filename, "%s.%d.%d.%d_%d.dimacs", prefix, ktree_n,ktree_k,ktree_p, i);
        print_message(0,"Writing file %s\n",filename);

        writer->set_out_file_name(filename);

        if (random)
            writer->shuffle(G, seed);

        writer->write_graph(G);

 	    delete G; 
	    delete H; 
	}

    LOG_CLOSE();
    return 0;
}



