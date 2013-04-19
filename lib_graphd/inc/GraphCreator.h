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

#ifndef GRAPHCREATOR_H_
#define GRAPHCREATOR_H_

#include "Graph.h"
#include "WeightedGraph.h"
#include "MutableGraph.h"
#include "WeightedMutableGraph.h"

using namespace std;

namespace Graph {
    class GraphCreator
    {
public:
    virtual WeightedMutableGraph *create_random_edge_subgraph(WeightedMutableGraph *wmg, int percent_edges);
    virtual WeightedMutableGraph *initialize_ktree(int n, int tw);
    // Creates the induced subgraph on the given list of vertices (by position in adj_list).
    // Does not check connectedness. Caller must allocate space for the graph G.
    virtual WeightedMutableGraph *create_induced_subgraph(WeightedMutableGraph *wmg, list<int> *members, bool make_simple);

    // Specialized induced subgraph functionality for components. Sets num_connected_components to 1.
    virtual WeightedMutableGraph *create_component(WeightedMutableGraph *g,
                                                   list<int> *members, bool make_simple);

    // Recursive and non-recursive functions that create all connected components of the graph,
    // which it places in the given vector
    // (each Graph is allocated with new inside the function). Labels transfer to adjlists.
    // Returns the number of components.
    virtual list<WeightedMutableGraph *> create_rec_all_components(WeightedMutableGraph *g, bool make_simple);
    virtual list<WeightedMutableGraph *> create_all_components(WeightedMutableGraph *g, bool make_simple);

    virtual Graph *create_graph() = 0;
    virtual WeightedGraph *create_weighted_graph() = 0;
    virtual MutableGraph *create_mutable_graph() = 0;
    virtual WeightedMutableGraph *create_weighted_mutable_graph() = 0;

    GraphCreator();
    virtual ~GraphCreator();
    };
}

#endif /* GRAPHCREATOR_H_ */
