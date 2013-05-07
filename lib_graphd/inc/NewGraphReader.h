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

#ifndef NEWGRAPHREADER_H_
#define NEWGRAPHREADER_H_

#include "Node.h"
#include "Graph.h"
#include "VertexWeightedGraph.h"
#include "Log.h"
#include "Debug.h"
#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

namespace Graph {
    class NewGraphReader {
public:

    NewGraphReader();
    virtual ~NewGraphReader();
    /**
     * \brief Reads in a graph of the specified type
     *
     * some more comments
     * on multiple lines
     * \param g Pointer to a graph
     * \param filename file to read graph from
     * \param type string specifying the type of the graph
     */
    int read_graph(Graph *g, const char *filename, char *type, bool read_vertex_weights);
private:
    int read_edgelist(Graph *g, const char *filename);
    int read_dimacs(Graph *g, const char *filename, bool read_vertex_weights);
    int read_adjmatrix(Graph *g, const char *filename);
    int read_metis(Graph *g, const char *filename);
    void split(const std::string& s, char sep, vector<int>& v);      // used in metisgraph reader
    };
}

#endif /* GRAPHREADER_H_ */
