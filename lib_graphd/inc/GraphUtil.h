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

#ifndef GRAPHUTIL_H_
#define GRAPHUTIL_H_

#include "GraphDecomposition.h"

namespace Graph {
    class GraphUtil
    {
public:
    GraphUtil();
    virtual ~GraphUtil();

    void recompute_degrees(Graph *g);
    // Returns the index of one of the vertices with maximal degree
    int get_random_high_degree_vertex(Graph *g) const;
    int get_random_low_degree_vertex(Graph *g) const;
    int label_component(Graph *g, int v, int label,
                        vector<int> *components);
    int rec_label_component(Graph *g, int v, int label,
                            vector<int> *components);

    // Recursive and non-recursive functions that label all connected
    // components of the graph by populating the given vector (of size capacity)
    // with component labels (in [0, num_connected_components-1]). Sets to GD_UNDEFINED if adj_list[i]
    // is not an active vertex. Returns the number of components.
    int label_all_components(Graph *g, vector<int> *components);
    int rec_label_all_components(Graph *g,
                                 vector<int> *components);

    // Recursive and non-recursive functions that find the vertices in the
    // connected component of v (by position in nodes[])
    // and fills the members list with their positions in the nodes[] array.
    // Returns the number of vertices in the component (including v).
    int rec_find_component(Graph *g, int v, list<int> *members);
    int find_component(Graph *g, int v, list<int> *members);

    // Non-recursive function that fills the vector with pointers to lists of
    // component members
    int find_all_components(Graph *g, vector<list<int> *> *members);

    int find_isolated_nodes(Graph *g, list<int> *isolated_nodes);

    void populate_CRS(Graph *g);
    void free_CRS(Graph *g);

    #ifdef HAS_BOOST
    /**
     * \brief Populate the boost version of the graph
     */
    void populate_boost(Graph *g);
    #endif

    //uses V as a vertex separator of G. Returns the number of components in G\V, and fills in members with their vertex lists..
    int vertex_separator(Graph *g, list<int> *V,
                         vector<list<int> *> *members);

    //runs a BFS from start using
    //only vertices with allowed[v] = true.
    //you need to delete the memory in the returned bool array.
    //num_reached is the number of vertices which are reachable. This does not include
    //disallowed vertices, and does include the start.
    bool *bfs(Graph *g, int start, bool *allowed, int *num_reached);
    bool *bfs_old(Graph *g, int start, bool *allowed, int *num_reached);
    int *bfs_dist(Graph *g, int start, bool *allowed, int *num_reached);
    int *bfs_dist(Graph *g, int start, bool *allowed, int *num_reached, int *ecc);

    //Find the eccentricity of each vertex and store it in ecc.
    void find_ecc(Graph *g, vector<int> *ecc);

    /**
     * \brief Find the k-core number of each vertex and store it in kcore.
     */
    int find_kcore(Graph *g, vector<int> *kcore);

    //Calculate the maximum distance between nodes within a subset of vertices
    //given as a list
    int subset_max_dist(Graph *g, vector<int> subset);
    /**
     * \brief Populate adjency matrix as a Petsc Matrix
     */
    void populate_PetscMat(Graph *g);
    /**
     * \brief Releases PetscMat.
     */
    void free_PetscMat(Graph *g);

    /**
     * \brief returns a new graph that is the largest component of g
     */
    Graph *get_largest_component_graph(Graph *g);
    };

    void create_largestcomponent_graph(char *graph_file, VertexWeightedGraph *&G);
}

#endif /* GRAPHUTIL_H_ */
