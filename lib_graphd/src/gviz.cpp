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

/**
   Main file which takes DIMACS file as input and output GraphViz
   file.
   Date: 03/02/2010
 **/

#include <iostream>
#include "NewGraphWriter.h"
#include <Log.h>
#include <GraphCreator.h>
#include <GraphCreatorFile.h>

using namespace Graph;

int main(int argc, char **argv){
    LOG_INIT("gviz.log", "gviz.log", 0);

    if(argc < 3){
        fprintf(stderr, "usage:%s dimacs_file output.gviz\n", argv[0]);
        FERROR("usage:%s dimacs_file output.gviz\n", argv[0]);
        return 0;
    }
    GraphCreator *gc;
    NewGraphWriter writer;
    gc = new GraphCreatorFile(argv[1], "DIMACS");
    VertexWeightedGraph *vwg = gc->create_vertex_weighted_graph();
    writer.write_graph(vwg, argv[2], "graphviz", false);

    LOG_CLOSE();
    delete gc;
    return 0;
} // main

