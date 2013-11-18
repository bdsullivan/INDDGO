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
#include "VertexWeightedGraph.h"
#include "GraphException.h"
#include "GraphWriter.h"

#if !WIN32 && !CYGWIN
  #include "orbconfig.h"
  #include "orbtimer.h"
#else
  #define ORB_t clock_t
  #define ORB_seconds(x,y) "lame"
  #define ORB_calibrate()
  #define ORB_read(x)
#endif

using namespace std;

void usage(const char *s){
    fprintf(stderr,"Usage: %s [options]\n",s);
    fprintf(stderr,
            "\t -t [#] sets number of graphs to generate (t)\n"
            "\t -n [#] Sets the number of vertices (n)\n"
            "\t -k [#] sets value for tree-width (k)\n"
            "\t -p [#] sets percent of edges to keep from complete k-tree (p)\n"
            "\t -s seed : uses the given RNG seed\n"
            "\t -m : print more precise timings to stdout\n"
            "\t -fn [name] sets the prefix for the filenames written out\n"
            "\t -fe [name] sets the exact filename written out (can only be used with one output graph)\n"
            "\t -o [type] sets output file type - valid values are: dimacs (default), adjlist, adjlist\n"
            "\t -r : randomizes the vertex labels before writing to file.\n"
            "\t -scotch : writes a .scotch file in Scotch graph format as well\n"
            );
}

void print_time(string prefix, ORB_t start, ORB_t end){
    cout << prefix + ": " << ORB_seconds(end, start) << "\n";
}

int main(int argc, char **argv){
    // This is a utility to generate partial k-trees with specified parameters for use in other applications.

    // Check for a cry for help
    if((argc == 1) || ((argc == 2) && (strcmp(argv[1],"-h") == 0)) || ((argc == 2) && (strcmp(argv[1],"--help") == 0))
       || ((argc == 2) && (strcmp(argv[1],"--h") == 0) ) ){
        usage(argv[0]);
        exit(-1);
    }

    ORB_t t1 = 0, t2 = 0, t3 = 0;
    int i;
    int t = 1;
    int ktree_n = -1, ktree_k = -1, ktree_p = -1;
    int seed = 0;
    Graph::VertexWeightedGraph *G = NULL, *H = NULL;
    time_t start, stop;
    char filename[100];
    char *prefix = (char *)"pkt";
    string out_format = "DIMACS";

    bool random = false;
    bool timings = false;
    bool exact_filename = false;
    char lfname[100];
    char efname[100];
    sprintf(lfname, "genpkt.log");
    sprintf(efname, "genpkt.log");
    //0: debug, 5: critical
    LOG_INIT(lfname, efname, 0);

    // Parse arguments
    for(i = 0; i < argc; i++){
        if(strcmp(argv[i],"-t") == 0){
            t = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i],"-k") == 0){
            ktree_k = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i],"-n") == 0){
            ktree_n = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i],"-p") == 0){
            ktree_p = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i],"-o") == 0){
            out_format = string(argv[i + 1]);
        }
        if(strcmp(argv[i],"-s") == 0){
            seed = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i],"-fn") == 0){
            prefix = argv[i + 1];
        }
        if(strcmp(argv[i],"-fe") == 0){
            prefix = argv[i + 1];
            exact_filename = true;
        }
        if(strcmp(argv[i],"-m") == 0){
            timings = true;
        }
        if(strcmp(argv[i],"-r") == 0){
            random = true;
        }
    }

    if((true == exact_filename) && (t != 1)){
        cerr << "Error: cannot specify both -t and -fe flags\n";
        exit(1);
    }

    if(("adjlist" != out_format) && ("dimacs" != out_format)){
        cerr << "Error: only allowed output formats are: adjlist, dimacs\n";
        exit(1);
    }

    if(!seed){
        // Set the seed to a rand int in 0,2^24
        seed = Graph::rand_int(0,0xffffff);
    }
    // Spin the RNG seed times
    for(int ss = 0; ss < seed; ss++){
        Graph::lcgrand(0);
    }

    //Make sure we set all the parameters to viable values.
    if((ktree_k < 0) || (ktree_n < 0) || (ktree_p < 0)){
        fatal_error("Failed to input viable parameters for n, k, p (%d,%d,%d)\n",
                    ktree_k, ktree_n, ktree_p);
    }

    Graph::GraphCreatorFile creator;
    Graph::GraphProperties prop;
    Graph::GraphWriter writer;

    ORB_calibrate();

    DEBUG("Graph generation loop\n");
    DEBUG("n : %d k: %d\n", ktree_n, ktree_p);
    for(i = 0; i < t; i++){
        // Create the Ktree
        //H = new Graph::Graph(ktree_n, ktree_k);
        cout << "Generating k-tree\n";
        H = creator.initialize_ktree(ktree_n, ktree_k);
        // Generate an n-node Graph G that is a partial k-tree derived from the k-tree
        start = clock();
        cout << "Deriving partial k-tree\n";
        ORB_read(t1);
        G = creator.create_random_edge_subgraph(H, ktree_p);
        ORB_read(t2);
        stop = clock();
        if(timings){
            print_time("Generation time", t1, t2);
        }
        fprintf(stderr,"Generated random partial k-tree (%d,%d) in %f sec.\n",
                ktree_n, ktree_k, ((double)(stop - start)) / CLOCKS_PER_SEC);
        //write it to a file
        if(!exact_filename){
            sprintf(filename, "%s.%d.%d.%d_%d.dimacs", prefix, ktree_n,ktree_k,ktree_p, i);
        }
        else {
            strncpy(filename,prefix, 99);
        }

        print_message(0,"Writing file %s\n",filename);

        //writer->set_out_file_name(filename);

        if(random){
            writer.set_shuffle(true);
            writer.set_shuffle_seed(seed);
        }

        ORB_read(t2);
        writer.write_graph(G,filename, out_format);
        ORB_read(t3);
        if(timings){
            print_time("Output time", t2, t3);
            print_time("Total time", t1, t3);
            cout << "Output edges: " << G->get_num_edges() << "\nOutput vertices: " << G->get_num_nodes();
            cout << "\nOutput format: " << out_format << "\nOutput filename: " << filename << "\n";
        }

        delete G;
        delete H;
    }

    LOG_CLOSE();
    exit(0);
} // main

