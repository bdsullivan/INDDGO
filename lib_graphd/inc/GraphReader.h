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

#ifndef GRAPHREADER_H_
#define GRAPHREADER_H_

#include "Node.h"
#include "Graph.h"
#include "VertexWeightedGraph.h"
#include "Log.h"
#include "Debug.h"
#include <stdio.h>
#if !WIN32
  #include <strings.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

namespace Graph {
    class GraphReader {
public:

    GraphReader();
    virtual ~GraphReader();
    /**
     * \brief Reads in a graph of the specified type
     *
     * some more comments
     * on multiple lines
     * \param g Pointer to a graph
     * \param filename file to read graph from
     * \param type string specifying the type of the graph
     */
    int read_graph(Graph *g, const string filename, string type, bool read_vertex_weights);
private:
    /** \brief read an edgelist */
    int read_edgelist(Graph *g, const string filename);
    /** \brief read a DIMACS file */
    int read_dimacs(Graph *g, const string filename, bool read_vertex_weights);
    /** \brief read an adjacency list */
    int read_adjlist(Graph *g, const string filename);
    /** \brief read an adjacency matrix */
    int read_adjmatrix(Graph *g, const string filename);
    /** \brief read a METIS file */
    int read_metis(Graph *g, const string filename);
    };
}

#endif /* GRAPHREADER_H_ */
