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
#include "Debug.h"
#include "Log.h"
#include "GraphException.h"
#include "GraphProperties.h"
#include <list>
#include <algorithm>
#include <string.h>


#ifndef _OPENMP
    #ifdef HAS_METIS
        void omp_set_num_threads(int num_threads) { return; }
        int omp_get_num_threads() { return 1; }
        int omp_get_max_threads(void) { return 1; }
        int omp_get_thread_num(void) { return 0; }
        int omp_get_num_procs(void) { return 1; }
        int omp_in_parallel(void) { return 0; }
        void omp_set_dynamic(int num_threads) { return; }
        int omp_get_dynamic(void) { return 0; }
        void omp_set_nested(int nested) { return; }
        int omp_get_nested(void) { return 0; }
    #endif
#endif

/* GRAPH_H_ */
namespace Graph {
    Graph::Graph(){
        this->num_nodes = 0;
        this->num_edges = 0;
        this->graph_type = "Graph";
        this->next_label = 1;
        this->num_connected_components = 1;
        this->simple = true;   //FIXME: check if this is correct behavior
        this->canonical = true;
        this->key = 0;
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
        this->next_label = n + 1;
    }

    Graph::~Graph(){
    }

    void Graph::set_canonical(bool c){
        canonical = c;
    }

    void Graph::set_degree(vector<int> degree){
        this->degree = degree;
    }

    vector<int> Graph::get_degree() const {
        return degree;
    }

    const vector<int> &Graph::get_degree_ref() const {
        return degree;
    }

    vector<Node> Graph::get_nodes() const {
        return nodes;
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

    int Graph::get_degree(int v) const {
        return this->degree[v];
    }

    /**
     * \param[in] n number of elements
     * */
    void Graph::resize(int n){
        nodes.reserve(n);
        nodes.resize(n);
        degree.reserve(n);
        degree.resize(n);
        num_nodes = n;
        capacity = n;
    }

    bool Graph::is_edge(int u, int v) const {
        //      Assumes graph is symmetric
        const list<int> &nbrs = nodes[u].get_nbrs_ref();  // passing a ref is much faster than copying a list
        list<int>::const_iterator it;
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

    bool Graph::is_canonical() const {
        return canonical;
    }

    bool Graph::is_simple() const {
        return simple;
    }

    void Graph::add_edge_advance(int u, int v){
        if((v < 0) || (u < 0) || (v >= capacity) || (u >= capacity) ){
            fatal_error(
                "%s:  add_edge_advance() called with vertices %d and %d but there are %d connected nodes\n",
                __FUNCTION__, u, v, capacity);
        }

        if(nodes[u].label == -1){
            fatal_error("%s:  Tried to use  GD_UNDEFINED vertex %d\n",
                        __FUNCTION__, u);
        }

        if(nodes[v].label == -1){
            fatal_error("%s:  Tried to use GD_UNDEFINED vertex %d\n", __FUNCTION__,
                        v);
        }

        if(simple == true){
            list<int>::iterator it;
            if(u == v){
                simple = false;
                canonical = false;
            }
            else {
                it = nodes[u].nbrs.begin();
                while(it != nodes[u].nbrs.end()){
                    if(*it == v){
                        simple = false;
                        canonical = false;
                        break;
                    }
                    ++it;
                }

                if(simple == true){
                    it = nodes[v].nbrs.begin();
                    while(it != nodes[v].nbrs.end()){
                        if(*it == u){
                            simple = true;
                            canonical = false;
                            break;
                        }
                        ++it;
                    }
                }
            }
        }

        nodes[u].nbrs.push_back(v);
        if(u != v){
            nodes[v].nbrs.push_back(u);
        }

        degree[u]++;
        if(u != v){
            degree[v]++;
        }

        num_edges++;
        if(num_connected_components > 1){
            num_connected_components = -1;
        }
    } // add_edge_advance

    bool Graph::remove_edge(int u, int v){
        if((v < 0) || (u < 0) || (v >= capacity) || (u >= capacity) ){
            fatal_error(
                "%s:  remove_edge() called with vertices %d and %d but there are %d connected nodes\n",
                __FUNCTION__, u, v, capacity);
        }
        if(nodes[u].label == -1){
            fatal_error("%s:  Tried to remove GD_UNDEFINED vertex %d\n",
                        __FUNCTION__, u);
        }

        if(nodes[v].label == -1){
            fatal_error("%s:  Tried to remove GD_UNDEFINED vertex %d\n",
                        __FUNCTION__, v);
        }

        list<int>::iterator it;
        bool found = false;
        it = nodes[u].nbrs.begin();
        while(it != nodes[u].nbrs.end()){
            if(*it == v){
                nodes[u].nbrs.erase(it);
                found = true;
                break;
            }
            ++it;
        }

        it = nodes[v].nbrs.begin();
        while(it != nodes[v].nbrs.end()){
            if(*it == u){
                nodes[v].nbrs.erase(it);
                found = true;
                break;
            }
            ++it;
        }

        num_connected_components = -1;
        if(found){
            degree[u]--;
            degree[v]--;
            num_edges--;
            return true;
        }
        return false;
    } // remove_edge

    void Graph::remove_vertex(int v){
        if(!canonical){
            print_message(0, "%s:  Graph must be in canonical form\n", __FUNCTION__);
        }

        if((v < 0) || (v >= capacity) ){
            fatal_error(
                "%s:  remove_vertex() called with vertex %d and there are %d connected nodes\n",
                __FUNCTION__, v, capacity);
        }

        if(nodes[v].label == -1){
            fatal_error("%s:  Tried to remove GD_UNDEFINED vertex %d\n",
                        __FUNCTION__, v);
        }

        if(nodes[v].nbrs.size() != 0){
            list<int>::iterator i, j;
            int u, len;
            num_edges = num_edges - nodes[v].nbrs.size();
            for(i = nodes[v].nbrs.begin(); i != nodes[v].nbrs.end(); ++i){
                u = *i;
                len = nodes[u].nbrs.size();
                nodes[u].nbrs.remove(v);
                if((int) ((((nodes[u].nbrs.size())))) == len){
                    fatal_error(
                        "%s: Didn't find node %d (mask=%d) in adjacency list of %d\n",
                        __FUNCTION__, v, u);
                }

                degree[u] = degree[u] - (len - (int) ((((nodes[u].nbrs.size())))));
            }
            nodes[v].nbrs.clear();
        }

        degree[v] = 0;
        nodes[v].label = -1;
        num_nodes--;
    } // remove_vertex

    int Graph::contract_edge(int u, int v){
        int i;
        if(simple != true){
            fatal_error("%s: called on a non-simple graph!\n", __FUNCTION__);
        }
        if((nodes[u].label == -1) || (nodes[v].label == -1) ){
            fatal_error(
                "%s: Cannot remove edge (%d, %d) as one of its vertices is undefined!\n",
                __FUNCTION__, u, v);
        }
        vector<bool> neighbors(capacity);
        fill(neighbors.begin(), neighbors.end(), false);
        list<int>::iterator it;
        bool foundv = false;
        it = nodes[u].nbrs.begin();
        while(it != nodes[u].nbrs.end()){
            neighbors[*it] = true;
            if(*it == v){
                foundv = true;
            }

            ++it;
        }
        if(foundv == false){
            return false;
        }

        it = nodes[v].nbrs.begin();
        while(it != nodes[v].nbrs.end()){
            neighbors[*it] = true;
            ++it;
        }
        remove_vertex(u);
        remove_vertex(v);
        nodes[u].label = next_label;
        next_label++;
        num_nodes++;
        for(i = 0; i < capacity; i++){
            if((i != u) && (i != v) ){
                if(neighbors[i] == true){
                    nodes[i].nbrs.push_back(u);
                    nodes[u].nbrs.push_back(i);
                    degree[u]++;
                    degree[i]++;
                    num_edges++;
                }
            }
        }

        return u;
    } // contract_edge

    void Graph::set_simple(bool si){
        simple = si;
    }

    void Graph::add_edge(int u, int v){
        if((v < 0) || (u < 0) || (v >= capacity) || (u >= capacity) ){
            fatal_error( "%s: called with vertices %d and %d but there are %d connected nodes\n",
                         __FUNCTION__, u, v, capacity);
        }

        if(nodes[u].label == -1){
            fatal_error("%s: Tried to use  GD_UNDEFINED vertex %d\n",
                        __FUNCTION__, u);
        }

        if(nodes[v].label == -1){
            fatal_error("%s: Tried to use GD_UNDEFINED vertex %d\n", __FUNCTION__,
                        v);
        }

        nodes[u].nbrs.push_back(v);
        if(u != v){
            nodes[v].nbrs.push_back(u);
        }

        degree[u]++;
        if(u != v){
            degree[v]++;
        }

        num_edges++;
    } // add_edge

    /**
     * does some stuff
     * and other stuff
     *
     */
    int Graph::add_vertices(int n){
        if(n < 1){
            return -1;
        }
        int old_size = this->num_nodes;
        this->num_nodes = this->num_nodes + n;
        this->nodes.resize(this->num_nodes);
        this->degree.resize(this->num_nodes);
        this->capacity = this->num_nodes;
        for(int i = old_size; i < this->num_nodes; i++){
            this->nodes[i].set_label(this->next_label);
            this->next_label++;
        }
        return this->num_nodes - 1;
    }

    /**
     * Convenience function to add a single vertex
     */
    int Graph::add_vertex(){
        return this->add_vertices(1);
    }

    void Graph::resize_adj_vec(int n){
        adj_vec.resize(n, 0);
    }

    vector<int> Graph::get_adj_vec() const {
        return adj_vec;
    }

    vector<int> *Graph::get_adj_vec_ptr(){
        return (&(this->adj_vec));
    }

    vector<int> Graph::get_adjncy() const {
        return adjncy;
    }

    vector<int> Graph::get_xadj() const {
        return xadj;
    }

    void Graph::set_adjncy(vector<int> adjncy){
        this->adjncy = adjncy;
    }

    void Graph::set_xadj(vector<int> xadj){
        this->xadj = xadj;
    }

    void Graph::eliminate_vertex(int v, list<int> *forward_neighbors,
                                 bool remove){
        if(!this->canonical){
            print_message(0, "%s:  Graph must be in canonical form\n", __FUNCTION__);
        }

        if((v < 0) || (v >= this->capacity) ){
            fatal_error(
                "%s:  eliminate_vertex() called with vertex %d and there are %d possible nodes\n",
                __FUNCTION__, v, this->capacity);
        }

        if(this->degree[v] == 0){
            return;
        }

        int i, j;
        GraphProperties properties;
        if(this->degree[v] == 1){
            print_message(1, "Degree[%d]=1\n", v);
            j = this->nodes[v].nbrs.front();
            this->nodes[j].nbrs.remove(v);
            if(remove){
                this->nodes[v].nbrs.pop_front();
            }

            this->num_edges--;
            this->degree[v]--;
            this->degree[j]--;
            return;
        }
        int *v_neighbors;
        char *forward_neighbor_vec;
        list<int>::iterator k;
        list<int>::iterator w;
        if(forward_neighbors != NULL){
            forward_neighbor_vec = new char[this->capacity];
            memset(forward_neighbor_vec, 0, (this->capacity) * sizeof(char));

            // Create the incidence vector
            for(w = forward_neighbors->begin(); w != forward_neighbors->end(); ++w){
                forward_neighbor_vec[*w] = 1;
            }
        }
        else {
            forward_neighbor_vec = NULL;
        }
        v_neighbors = new int[this->degree[v]];
        memset(v_neighbors, 0, this->degree[v] * sizeof(int));
        j = 0;
        for(k = this->nodes[v].nbrs.begin(); k != this->nodes[v].nbrs.end(); ++k){
            v_neighbors[j] = *k;
            j++;
        }
        print_message(100, "Found %d neighbors of %d\nDeg is %d\n", j, v,
                      this->degree[v]);
        for(i = 0; i < j; i++){
            print_message(100, "v_neighbors[%d]=%d\n", i, v_neighbors[i]);
        }

        if(forward_neighbor_vec){
            for(i = 0; i < this->degree[v]; i++){
                if(forward_neighbor_vec[v_neighbors[i]] == 1){
                    int current_key = properties.fill_adj_vec(this, v_neighbors[i]);
                    for(j = i + 1; j < this->degree[v]; j++){
                        if(forward_neighbor_vec[v_neighbors[j]] == 1){
                            if(this->adj_vec[v_neighbors[j]] != current_key){
                                this->nodes[v_neighbors[i]].nbrs.push_back(
                                    v_neighbors[j]);
                                this->degree[v_neighbors[i]]++;
                                this->nodes[v_neighbors[j]].nbrs.push_back(
                                    v_neighbors[i]);
                                this->degree[v_neighbors[j]]++;
                                this->num_edges++;
                            }
                        }
                    }
                }

                if(remove){
                    this->nodes[v_neighbors[i]].nbrs.remove(v);
                    this->degree[v_neighbors[i]]--;
                    this->num_edges--;
                }
            }
        }
        else {
            for(i = 0; i < this->degree[v]; i++){
                int current_key = properties.fill_adj_vec(this, v_neighbors[i]);
                for(j = i + 1; j < this->degree[v]; j++){
                    if(this->adj_vec[v_neighbors[j]] != current_key){
                        this->nodes[v_neighbors[i]].nbrs.push_back(v_neighbors[j]);
                        this->degree[v_neighbors[i]]++;
                        this->nodes[v_neighbors[j]].nbrs.push_back(v_neighbors[i]);
                        this->degree[v_neighbors[j]]++;
                        this->num_edges++;
                    }
                }

                if(remove){
                    this->nodes[v_neighbors[i]].nbrs.remove(v);
                    this->degree[v_neighbors[i]]--;
                    this->num_edges--;
                }
            }
        }

        if(remove){
            this->nodes[v].nbrs.clear();
            this->degree[v] = 0;
        }
        delete[] v_neighbors;
        if(forward_neighbor_vec){
            delete[] forward_neighbor_vec;
        }

        return;
    } // eliminate_vertex

    /*
     * Calculates number of connected components, checks whether the graph is simple/canonical.
     * Also checks and sets num_components.
     *
     */
    void Graph::initialize_params(){
        fatal_error("Initialize params not yet implemented\n");

        //Graph::GraphProperties properties;
        //properties.make_canonical(G);
        //if (!properties.is_connected(G))
        //	list<Graph::VertexWeightedGraph *> C;
        //	C = creator.create_all_components(G, true);
        //  print_message(0, "Found %d components\n", C.size());

        return;
    }

    /*
     * For advanced users - allows you to set number of connected components, and whether the graph is simple/canonical
     * Will not verify that the values passed in are accurate.
     */
    void Graph::initialize_params(bool simple, bool canonical, int num_comp){
        this->simple = simple;
        this->canonical = canonical;
        this->num_connected_components = num_comp;
        this->next_label = this->num_nodes + 1;

        return;
    }

    /**
     * \param[in] apsp_dist all pairs shortest paths distances
     **/
    void Graph::set_shortest_path_dist(vector< vector<int> > apsp_dist){
        this->apsp_dist = apsp_dist;
    }

    /**
     * Checks if distance matrix has been computed (or needs to be recomputed)
     * if yes, returns the all pairs shortest paths distances
     * otherwise computes then returns
     **/
    const vector< vector<int> > &Graph::get_shortest_path_dist_ref(){
        if(this->apsp_dist.empty()){
            cout << "Empty -- calling function to compute shortest paths" << endl;
            GraphProperties properties;
            properties.paths_dijkstra_all(this,this->apsp_dist);   //sets this>apsp_dist with values
            return this->apsp_dist;
        }
        return this->apsp_dist;
    }

    const vector<int> &Graph::get_u_shortest_path_dist(int u){
        if(this->apsp_dist[u].empty()){
            cout << "u Empty -- calling function to compute shortest paths" << endl;
            GraphProperties properties;
            properties.paths_dijkstra_single(this,this->apsp_dist[u], u);   //sets this>apsp_dist[u] with values

            return this->apsp_dist[u];
        }

        return this->apsp_dist[u];
    }
}
