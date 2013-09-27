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
#include "TreeDecomposition.h"
#include "Log.h"
#include "GraphException.h"
#include <fstream>
#include <string>

#include "VertexWeightedGraph.h"

void usage(char *str){
    fprintf(stderr,"%s <DIMACS_file>\n\tWrites the complement graph of the provided\n"
                   "\tDIMACS file to <DIMACS_file>.comp.norm\n"
                   "\tin normalized DIMACS form\n", str);
    exit(-1);
};

int main(int argc, char **argv){
    if((argc != 2) || (strcmp(argv[1],"-h") == 0) || (strcmp(argv[1],"--help") == 0) ){
        usage(argv[0]);
    }

    Graph::VertexWeightedGraph *G;

    Graph::GraphCreatorFile creator;
    creator.set_file_name(argv[1]);
    creator.set_graph_type("DIMACS");
    G = creator.create_weighted_mutable_graph();
    G->complement();
    Graph::GraphWriter writer;

    char comp_DIMACS_file[200];
    char normalized_comp_DIMACS_file[200];
    sprintf(comp_DIMACS_file,"%s.comp",argv[1]);
    writer.write_graph(G,comp_DIMACS_file,"DIMACS",true,false);
    // Now normalize it
    sprintf(normalized_comp_DIMACS_file,"%s.norm",comp_DIMACS_file);
    normalize_DIMACS_file(comp_DIMACS_file,normalized_comp_DIMACS_file);
    fprintf(stderr,"Complement graph of %s written to %s in normalized DIMACS form\n",
            argv[1],normalized_comp_DIMACS_file);

    delete G;
    return 1;
} // main

