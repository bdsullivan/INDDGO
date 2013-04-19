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

#include "DIMACSGraphWriter.h"
#include "Debug.h"
#include <stdio.h>
//#include <random>
#include <algorithm>

namespace Graph {
    DIMACSGraphWriter::DIMACSGraphWriter(){
    }

    DIMACSGraphWriter::DIMACSGraphWriter(string & filename) :
        GraphWriter(filename){
    }

    DIMACSGraphWriter::~DIMACSGraphWriter(){
    }

    /**
     * Writes the graph to a 1-based DIMACS file. If randomize == true,
     * we permute the vertex labels (at random) before writing.
     */
    void DIMACSGraphWriter::write_graph(Graph *g){
        int capacity = g->get_capacity();
        int num_nodes = g->get_num_nodes();
        int num_edges = g->get_num_edges();
        vector<Node> nodes = g->get_nodes();

        if(perm.size() != capacity){
            perm.clear();
            for(int i = 0; i < capacity; i++){
                perm.push_back(i);
            }
        }

        FILE *out;
        if((out = fopen(out_file_name.c_str(), "w")) == NULL){
            fatal_error("%s:  Couldn't open %s for writing\n", __FUNCTION__,
                        out_file_name.c_str());
        }

        // Print out a comment line
        //fprintf(out,"c Graph has %d connected components\n",g->num_connected_components);

        // Will write out non-canonical - note that we may have "holes" in the numbering
        // Print a warning to stderr and the file in this case
        if(capacity != num_nodes){
            fprintf(stderr, "%s:  Warning! DIMACS file may have holes.\n",
                    __FUNCTION__);
            // Print a comment warning in the file
            fprintf(
                out,
                "c Warning: capacity=%d, num_nodes=%d. Resulting graph may have holes\n",
                capacity, num_nodes);
        }

        //EDITED TO PRINT p edge instead of p e for compatibility with Koster.
        fprintf(out, "p edge %d %d\n", num_nodes, num_edges);
        // Write edges out using the label field from the nodes
        // Only write out edges u-v with u<=v - so file is at least simple
        list<int>::iterator ii;
        int edges_written = 0;
        list<int> nbrs;
        for(int i = 0; i < num_nodes; i++){
            // CSG- could add sort here although that won't work for the randomized labels
            // g->nodes[perm[i]].nbrs.sort();
            nbrs = nodes[i].get_nbrs();
            for(ii = nbrs.begin(); ii != nbrs.end(); ++ii){
                //if(i<=*ii)
                // CSG - we want the edges in the file to be u-v with u<=v
                if(nodes[perm[i]].get_label() <= nodes[perm[*ii]].get_label()){
                    // Changing june 21 for testing - not using labels
                    //if(i<=*ii)
                    edges_written++;
                    fprintf(out, "e %d %d\n", nodes[perm[i]].get_label(),
                            nodes[perm[*ii]].get_label());
                    //,i+1,*ii+1);//
                }
            }
        }
        if(edges_written != num_edges){
            fatal_error("%s:  Didn't write the correct number of edges! (%d!=%d)\n",
                        __FUNCTION__, edges_written, num_edges);
        }

        fflush(out);
        fclose(out);

        return;
    } // write_graph
}
