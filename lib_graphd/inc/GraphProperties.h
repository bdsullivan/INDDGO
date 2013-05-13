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

#ifndef GRAPHPROPERTIES_H_
#define GRAPHPROPERTIES_H_
#include "GraphDecomposition.h"

using namespace std;

namespace Graph {
    class GraphProperties
    {
public:
    GraphProperties();
    virtual ~GraphProperties();
    void make_simple(Graph *g);             //Removes all loops and duplicate edges. Sets the simple flag to true.
    void make_canonical(Graph *g);
    int make_clique(Graph *g, list<int> *vertices);
    int fill_adj_vec(Graph *g, int v);
    bool check_simple(Graph *g);
    // Connectivity functions

    bool is_connected(Graph *g);

    bool is_clique(Graph *g, list<int> *vertices);
    bool is_independent_set(Graph *g, list<int> *vertices);
    bool is_independent_set(VertexWeightedGraph *mg, list<int> *vertices,
                            int *val);
    bool is_path(Graph *g, int start, int end, bool *t);             // Determines existence of a path between start
    // and end involving only vertices v with t[v]=true
    bool is_path(Graph *g, int start, int end);
    };
}

#endif /* GRAPHPROPERTIES_H_ */
