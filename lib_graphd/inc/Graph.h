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

#ifndef GRAPH_H_
#define GRAPH_H_

#include "GraphInterface.h"
#include "Node.h"
#include <string>

using namespace std;

namespace Graph {
    class Graph : public GraphInterface
    {
protected:
    int num_nodes;
    int num_edges;
    int capacity;
    int next_label;
    int num_connected_components;
    string graph_type;
    vector<int> degree;
    vector<Node> nodes;
    string input_file;
    bool simple;

    bool canonical;
    int key;
    // For CRS format, a la Metis
    vector<int> xadj;
    vector<int> adjncy;
    vector<int> adj_vec;

public:
    Graph();
    Graph(int n);
    virtual ~Graph();

    /* accessor functions */
    /** \brief set the nodes of the graph to the specified nodes */
    void set_canonical(bool canonical);
    /** \brief set the degree of the graph */
    void set_degree(vector<int> degree);
    void set_capacity(int capacity);
    void set_graph_type(string graphType);
    void set_input_file(string input_file);
    void set_nodes(vector<Node> nodes);
    void set_num_edges(int numEdges);
    void set_num_nodes(int numNodes);
    void set_simple(bool simple);
    void set_next_label(int nextLabel);
    void set_num_components(int num_components);
    void set_adjncy(vector<int> adjncy);
    void set_xadj(vector<int> xadj);

    vector<int> get_adj_vec() const;
    vector<int> *get_adj_vec_ptr();
    vector<int> get_adjncy() const;
    virtual int get_capacity() const;
    vector<int> get_degree() const;
    virtual int get_degree(int v) const;
    string get_graph_type() const;
    string get_input_file();
    int get_next_label() const;
    Node *get_node(int i);
    vector<Node> get_nodes() const;
    int get_num_components() const;
    int get_num_edges() const;
    int get_num_edges_in_subgraph(list<int> *vertices);
    int get_num_nodes() const;
    vector<int> get_xadj() const;

    virtual bool is_edge(int i, int j) const;
    bool is_canonical() const;
    bool is_simple() const;
    virtual int *serialize();
    virtual void deserialize(int *buffer);
    void complement();

    /* other methods from mutable */
    void add_edge(int u, int v);
    void add_edge_advance(int u, int v);
    bool remove_edge(int u, int v);
    void remove_vertex(int u);
    int contract_edge(int u, int v);
    void resize_adj_vec(int n);
    void eliminate_vertex(int v, list<int> *forward_neighbors, bool remove);
    void initialize_params();
    void initialize_params(bool simple, bool canonical, int num_comp);

    friend class GraphUtil;
    friend class GraphProperties;
    friend class GraphDisplay;
    friend class GraphEOUtil;
    friend class GraphCreator;
    };
}

#endif /* GRAPH_H_ */

