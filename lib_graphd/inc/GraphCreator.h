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
#include "VertexWeightedGraph.h"

using namespace std;

namespace Graph {
    class GraphCreator
    {
public:
    virtual VertexWeightedGraph *create_random_edge_subgraph(VertexWeightedGraph *wmg, int percent_edges);
    virtual VertexWeightedGraph *initialize_ktree(int n, int tw);
    virtual Graph *initialize_rig(int n, int seed, list<double> *probs);
    virtual Graph *initialize_rmat(int l, int m, double *probs, int seed, bool self_loop);

    // Creates the induced subgraph on the given list of vertices (by position in adj_list).
    // Does not check connectedness. Caller must allocate space for the graph G.
    virtual VertexWeightedGraph *create_induced_subgraph(VertexWeightedGraph *wmg, list<int> *members, bool make_simple);

    // Specialized induced subgraph functionality for components. Sets num_connected_components to 1.
    virtual VertexWeightedGraph *create_component(VertexWeightedGraph *g,
                                                  list<int> *members, bool make_simple);
    /**
     * \brief Create an induced subgraph on g
     */
    virtual Graph *create_induced_subgraph(Graph *g, list<int> *members, bool make_simple);
    /**
     * \brief Create a component graph of g
     */
    virtual Graph *create_component(Graph *g, list<int> *members, bool make_simple);

    // Recursive and non-recursive functions that create all connected components of the graph,
    // which it places in the given vector
    // (each Graph is allocated with new inside the function). Labels transfer to adjlists.
    // Returns the number of components.
    virtual list<VertexWeightedGraph *> create_rec_all_components(VertexWeightedGraph *g, bool make_simple);
    virtual list<VertexWeightedGraph *> create_all_components(VertexWeightedGraph *g, bool make_simple);

    virtual Graph *create_graph() = 0;
    virtual Graph *create_mutable_graph() = 0;
    virtual VertexWeightedGraph *create_vertex_weighted_graph() = 0;

    GraphCreator();
    virtual ~GraphCreator();
    };
}

#endif /* GRAPHCREATOR_H_ */
