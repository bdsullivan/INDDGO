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

    /**
     * \brief Sorts a nodes neighor list by degree
     */
    void sort_nbrs_by_map(vector<int> map, Node *n);

    /**
     * \brief Returns an array of all triangles
     */
    void all_triangles_compact_forward(Graph *g, vector<long int> &t);
    /**
     * \brief Returns an array of all triangles
     */
    void all_triangles_edge_listing(Graph *g, vector<long int> &t);

    /**
     * \brief returns clustering coeffecients
     */
    void clustering_coefficients(Graph *g,double &global_cc, double &avg_cc, vector<double> &local_ccs);

    /**
     * \brief returns shortest paths from source to all other nodes
     */
    void paths_dijkstra_single(Graph *g, vector<int> &p, int source);

    /**
     * \brief returns shortest paths from all nodes to all nodes
     */
    void paths_dijkstra_all(Graph *g, vector< vector<int> > &p);

    /**
     * \brief returns eccentricity for all nodes
     */
    void eccentricity(Graph *g, vector<int> &ecc);

    /**
     * \brief returns frequency distributions of eccentricities
     */
    void eccentricity_dist(Graph *g, vector<int> &ecc, vector<double> &freq_ecc);

    /**
     * \brief returns normalized expansion (distance distribution)
     */
    void expansion(Graph *g, vector<double> &norm_hops);

    /**
     * \brief Calculates the diameter of the specified graph
     */
    void diameter(Graph *g, int &diam);

    /**
     * \brief Calculates the effective diameter of the specified graph
     */
    void effective_diameter(Graph *g, float &ediam, float perc = 0.9);

    /**
     * \brief Calculates the edge density of the specified graph
     */
    void edge_density(Graph *g, float &ed);

    /**
     * \brief Calculates the average degree of the specified graph
     */
    void avg_degree(Graph *g, float &ad);

    /**
     * \brief Calculates the average shortest path length of the specified graph
     */
    void avg_path_length(Graph *g, double &pl);

    /**
     * \brief Calculates the degree distribution of the specified graph
     */
    void deg_dist(Graph *g, vector<int> &dist);

    /**
     * \brief Returns the degree assortativity coefficient of a graph
     */
    void deg_assortativity(Graph *g, double &coeff);
    };
}

#endif /* GRAPHPROPERTIES_H_ */
