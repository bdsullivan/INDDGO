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

#include "AdjMatrixGraphReader.h"
#include "Debug.h"
#include <stdio.h>

namespace Graph {
    AdjMatrixGraphReader::AdjMatrixGraphReader(){
    }

    AdjMatrixGraphReader::~AdjMatrixGraphReader(){
    }

    /**
     * Reads the file and populates the Graph data fields.
     * Ignores any self-loop edges.
     * If an error is found, some hopefully relevant information is printed to
     * stderr.
     * The symmetric flag should be true if you want symmetric adjacency lists and false for anti_symmetric ones.
     */
    void AdjMatrixGraphReader::read_graph(const char *mat_file){
        int i, j, k, m, n, retval;
        FILE *in;

        if((in = fopen(mat_file, "r")) == NULL){
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        mat_file);
        }

        fscanf(in, "%d\n", &n);
        //printf("%s:  Graph has %d nodes\n", __FUNCTION__, n);
        // Each row will have n binary entries, and there will be n rows.
        this->degree.resize(n, 0);
        this->nodes.resize(n);
        this->capacity = n;
        i = 0;         // Use i to keep track of row
        m = 0;         // To count edges
        while(!feof(in)){
            if(i >= n){
                break;
            }
            for(k = 0; k < n; k++){
                retval = fscanf(in, "%d ", &j);
                //printf("%s:  Read in %d\n", __FUNCTION__, j);
                if(retval != 1){
                    fatal_error("Didn't read digit properly from matrix\n");
                }
                if((j != 0) && (j != 1) ){
                    fatal_error("%s:  Encountered entry %d!=0 or 1!\n",
                                __FUNCTION__, j);
                }
                if(j == 1){
                    //ignore self loops!
                    if(i != k){
                        nodes[i].add_nbr(k);
                        degree[i]++;
                        m++;
                    }
                }
            }
            // Set 1-based labels?
            nodes[i].set_label(i + 1);
            // Advance to next row
            i++;
        }
        fclose(in);

        this->num_nodes = n;
        if(2 * (m / 2) != m){
            fatal_error(
                "%s: Found an odd number of ones (%d). Non-symmetric matrix. Please correct.\n",
                __FUNCTION__, m);
        }
        this->num_edges = m / 2;         // Assume we were given a symmetric matrix!!
        this->capacity = n;

        for(i = 0; i < this->capacity; i++){
            if(nodes[i].get_label() == -1){
                nodes[i].set_label(i + 1);
            }
        }

        return;
    } // read_graph

    vector<int> AdjMatrixGraphReader::get_degree(){
        return degree;
    }

    vector<Node> AdjMatrixGraphReader::get_nodes(){
        return nodes;
    }

    int AdjMatrixGraphReader::get_num_edges() const {
        return num_edges;
    }

    vector<int> AdjMatrixGraphReader::get_weights(){
        return weights;
    }

    int AdjMatrixGraphReader::get_capacity() const {
        return capacity;
    }

    void AdjMatrixGraphReader::set_capacity(int capacity){
    }
}
