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

#include "MetisGraphWriter.h"
#include "Debug.h"
#include <stdio.h>

namespace Graph {
    MetisGraphWriter::MetisGraphWriter(){
    }

    MetisGraphWriter::MetisGraphWriter(string & filename) :
        GraphWriter(filename){
    }

    MetisGraphWriter::~MetisGraphWriter(){
    }

    /**
     * Writes the graph to a 1-based Metis file.
     */
    void MetisGraphWriter::write_graph(Graph *g){
        int capacity = g->get_capacity();
        int num_nodes = g->get_num_nodes();
        int num_edges = g->get_num_edges();
        vector<Node> nodes = g->get_nodes();

        FILE *out;
        if((out = fopen(out_file_name.c_str(), "w")) == NULL){
            fatal_error("%s:  Couldn't open %s for writing\n", __FUNCTION__,
                        out_file_name.c_str());
        }

        if(capacity != num_nodes){
            fprintf(stderr, "%s:  Warning! Metis file may have holes.\n",
                    __FUNCTION__);
            // Print a comment warning in the file
            fprintf(
                out,
                "c Warning: capacity=%d, num_nodes=%d. Resulting graph may have holes\n",
                capacity, num_nodes);
        }

        fprintf(out, "%d %d\n", num_nodes, num_edges);
        list<int>::iterator ii;
        int edges_written = 0;
        list<int> *nbrs;
        for(int i = 0; i < num_nodes; i++){
            nbrs = nodes[i].get_nbrs_ptr();
            for(ii = nbrs->begin(); ii != nbrs->end(); ++ii){
                fprintf(out, "%d ", (*ii) + 1);
            }
            fprintf(out, "\n");
        }

        fflush(out);
        fclose(out);
        return;
    } // write_graph
}
