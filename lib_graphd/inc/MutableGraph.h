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

#ifndef MUTABLEGRAPH_H_
#define MUTABLEGRAPH_H_

#include "Graph.h"

namespace Graph {
    class MutableGraph : virtual public Graph
    {
protected:
    bool simple;
    bool canonical;
    int key;
    // For CRS format, a la Metis
    vector<int> xadj;
    vector<int> adjncy;
    vector<int> adj_vec;

public:
    MutableGraph();
    MutableGraph(int n);
    virtual ~MutableGraph();
    void add_edge(int u, int v);
    void add_edge_advance(int u, int v);
    bool remove_edge(int u, int v);
    void remove_vertex(int u);
    int contract_edge(int u, int v);
    bool is_canonical() const;
    bool is_simple() const;
    void set_canonical(bool canonical);
    void set_simple(bool simple);
    void resize_adj_vec(int n);
    void eliminate_vertex(int v, list<int> *forward_neighbors, bool remove);
    vector<int> get_adjncy() const;
    vector<int> get_xadj() const;
    void set_adjncy(vector<int> adjncy);
    void set_xadj(vector<int> xadj);
    vector<int> get_adj_vec() const;
    vector<int> *get_adj_vec_ptr();
    void initialize_params();
    void initialize_params(bool simple, bool canonical, int num_comp);

    friend class GraphUtil;
    friend class GraphProperties;
    friend class GraphDisplay;
    friend class GraphEOUtil;
    friend class GraphCreator;
    };
}

#endif /* MUTABLEGRAPH_H_ */
