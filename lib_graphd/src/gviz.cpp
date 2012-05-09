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
#include <GraphReaderWriterFactory.h>
#include <GraphWriter.h>
#include <Log.h>
#include <GraphCreator.h>
#include <GraphCreatorFile.h>

using namespace Graph;

int main(int argc, char **argv)
{
    LOG_INIT("gviz.log", "gviz.log", 0);

    if (argc < 3)
    {
        fprintf(stderr, "usage:%s dimacs_file output.gviz\n", argv[0]);
        FERROR("usage:%s dimacs_file output.gviz\n", argv[0]);
        return 0;
    }
    GraphCreator *gc;
    GraphWriter *writer;
    GraphReaderWriterFactory rwf;
    gc = new GraphCreatorFile(argv[1], "DIMACS");
    WeightedMutableGraph *wmg = gc->create_weighted_mutable_graph();
    writer = rwf.create_writer("graphviz");
    writer->set_out_file_name(argv[2]);
    writer->write_graph(wmg);

    LOG_CLOSE();
    delete  gc;
    delete writer;
    return 0;
}
