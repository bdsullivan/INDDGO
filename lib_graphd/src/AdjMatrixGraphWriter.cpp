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

#include "AdjMatrixGraphWriter.h"
#include "Debug.h"
#include <iostream>
#include <stdio.h>

using namespace std;

namespace Graph {
    AdjMatrixGraphWriter::AdjMatrixGraphWriter(){
    }

    AdjMatrixGraphWriter::~AdjMatrixGraphWriter(){
    }

    AdjMatrixGraphWriter::AdjMatrixGraphWriter(string &filename) :
        GraphWriter(filename){
    }

    void AdjMatrixGraphWriter::write_graph(Graph *g){
        FILE *fout = fopen(out_file_name.c_str(), "w");
        if(!fout){
            fatal_error("can not open %s for write", out_file_name.c_str());
        }

        int n = g->get_num_nodes();
        int i = 0;
        int j = 0;
        int val = 0;

        fprintf(fout, "%d \n", n);

        for(i = 0; i < n; i++){
            for(j = 0; j <  n; j++){
                val = 0;
                if(g->is_edge(i, j)){
                    val = 1;
                }
                fprintf(fout, "%d ", val);
            }
            fprintf(fout, "\n");
        }

        fclose(fout);
    } // write_graph
}
