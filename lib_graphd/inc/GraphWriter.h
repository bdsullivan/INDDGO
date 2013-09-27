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

#ifndef GRAPHWRITER_H_
#define GRAPHWRITER_H_

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
    class GraphWriter {
public:

    GraphWriter();
    virtual ~GraphWriter();
    /**
     * \brief Writes a graph of the specified type
     */
    int write_graph(Graph *g, const string filename, const string type, bool write_vertex_weights = false, bool shuffle = false);

    /**
     * \brief Sets whether to shuffle output graphs (for formats where we support it)
     */
    void set_shuffle(bool shuf);
    /**
     * \brief Sets the seed for shuffling
     */
    void set_shuffle_seed(int seed);
private:
    bool shuffle;
    int shuffle_seed;
    /**
     * \brief Writes a given graph to file in the METIS format
     */
    int write_metis(Graph *g, const string filename);
    /**
     * \brief Writes a given graph to file in the DIMACS format
     */
    int write_dimacs(Graph *g, const string filename, bool write_vertex_weights);
    /**
     * \brief Writes a given graph to file in our adjacency matrix format
     */
    int write_adjmatrix(Graph *g, const string filename);
    /**
     * \brief Writes a given graph to file in GraphViz format
     */
    int write_graphviz(Graph *g, const string filename);
    /**
     * \brief Writes a given graph to file in Adjacency List format
     */
    int write_adjlist(Graph *g, const string filename);
    };
}

#endif /* GRAPHWRITER_H_ */
