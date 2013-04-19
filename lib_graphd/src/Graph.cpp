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

#include "Graph.h"
#include "Log.h"
#include "GraphException.h"
#include <list>
#include <algorithm>

/* GRAPH_H_ */
namespace Graph {
    Graph::Graph(){
        this->num_nodes = 0;
        this->num_edges = 0;
        this->graph_type = "DIMACS";
        this->next_label = 1;
        this->num_connected_components = 1;
    }

    vector<int> Graph::get_degree() const {
        return degree;
    }

    vector<Node> Graph::get_nodes() const {
        return nodes;
    }

    void Graph::set_degree(vector<int> degree){
        this->degree = degree;
    }

    string Graph::get_graph_type() const {
        return graph_type;
    }

    int Graph::get_num_edges() const {
        return num_edges;
    }

    int Graph::get_num_nodes() const {
        return num_nodes;
    }

    void Graph::set_graph_type(string graphType){
        this->graph_type = graphType;
    }

    void Graph::set_num_edges(int numEdges){
        this->num_edges = numEdges;
    }

    void Graph::set_num_nodes(int numNodes){
        this->num_nodes = numNodes;
    }

    void Graph::set_nodes(vector<Node> nodes){
        this->nodes = nodes;
    }

    Graph::~Graph(){
    }

    int Graph::get_degree(int v) const {
        return this->degree[v];
    }

    bool Graph::is_edge(int u, int v) const {
        //      Assumes graph is symmetric
        list<int> nbrs = nodes[u].get_nbrs();
        list<int>::iterator it;
        for(it = nbrs.begin(); it != nbrs.end(); ++it){
            if(*it == v){
                return true;
            }
        }
        return false;
    }

    int *Graph::serialize(){
        return 0;
    }

    void Graph::deserialize(int *buffer){
    }

    Node *Graph::get_node(int i){
        if(i > capacity){
            FERROR("%s: element is out of bounds", __FUNCTION__);
            throw GraphException("element is out of bounds\n");
        }
        return &nodes[i];
    }

    int Graph::get_next_label() const {
        return next_label;
    }

    int Graph::get_num_components() const {
        return num_connected_components;
    }

    Graph::Graph(int n){
        nodes.resize(n);
        degree.resize(n, 0);
        num_nodes = n;
        capacity = n;
        num_connected_components = 1;
        num_edges = 0;
        for(int i = 0; i < n; i++){
            nodes[i].set_label(i + 1);
        }
    }

    void Graph::set_num_components(int num_components){
        this->num_connected_components = num_components;
    }

    void Graph::set_next_label(int nextLabel){
        this->next_label = nextLabel;
    }

    void Graph::set_capacity(int capacity){
        this->capacity = capacity;
    }

    int Graph::get_capacity() const {
        return this->capacity;
    }

    void Graph::set_input_file(string input_file){
        this->input_file = input_file;
    }

    string Graph::get_input_file(){
        return this->input_file;
    }

    int Graph::get_num_edges_in_subgraph(list<int> *vertices){
        vertices->sort();
        int highest_index = vertices->back();
        vector<bool> v(this->num_nodes,false);
        list<int>::iterator ii;
        for(ii = vertices->begin(); ii != vertices->end(); ++ii){
            v[*ii] = true;
        }
        list<int> nbrs;
        int cnt = 0;
        for(int i = 0; i < this->num_nodes; i++){
            if(v[i]){
                // This is a totally unnecessary copy
                nbrs.clear();
                nbrs = this->nodes[i].get_nbrs();
                for(ii = nbrs.begin(); ii != nbrs.end(); ++ii){
                    if(v[*ii]){
                        // This is an edge where both endpoints are in the input list
                        cnt++;
                    }
                }
            }
        }
        return (cnt >> 1);
    } // get_num_edges_in_subgraph

    void Graph::complement(){
        Node *n;

        list<int> elements;
        list<int> newnbrs;
        list<int> nbrs;
        list<int>::iterator it;
        int edge = 0;

        for(int i = 0; i < this->capacity; i++){
            elements.push_back(i);
        }

        for(int i = 0; i < this->capacity; i++){
            n = this->get_node(i);
            nbrs = n->get_nbrs();
            newnbrs = elements;
            nbrs.sort();

            newnbrs.remove(i);
            for(it = nbrs.begin(); it != nbrs.end(); ++it){
                newnbrs.remove(*it);
            }

            this->degree[i] = newnbrs.size();
            edge += this->degree[i];
            n->set_nbr(newnbrs);
        }

        this->set_num_edges(edge / 2);
    } // complement
}
