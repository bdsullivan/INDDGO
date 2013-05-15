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

#include "GraphDecomposition.h"

namespace Graph {
    GraphProperties::GraphProperties(){
    }

    GraphProperties::~GraphProperties(){
        // TODO Auto-generated destructor stub
    }

    /**
     * Removes all loops and duplicate edges. Sets the simple flag to true.
     */
    void GraphProperties::make_simple(Graph *g){
        int i;
        list<int>::iterator it;
        GraphUtil graph_util;

        for(i = 0; i < g->capacity; i++){
            //loop through the active nodes
            if(g->nodes[i].label != -1){
                //necessary prereq to using unique()
                g->nodes[i].nbrs.sort();
                //remove duplicate edges
                g->nodes[i].nbrs.unique();
                //remove loops
                g->nodes[i].nbrs.remove(i);
            }
        }

        g->simple = true;
        //update degrees and number of edges
        graph_util.recompute_degrees(g);
        return;
    } // make_simple

    /**
     * g forces the graph to have symmetric adjacency lists and be simple (no loops/duplicate edges) and updates flags accordingly.
     */
    void GraphProperties::make_canonical(Graph *g){
        // CSG -Should we also have g->capacity==g->num_nodes?
        make_simple(g);
        g->canonical = true;

        return;
    }

    /**
     * Increments the value of Graph.key and sets adj_vec[u]=key for all neighbors u of v.
     * returns the value of the key used.
     */
    int GraphProperties::fill_adj_vec(Graph *g, int v){
        if(!g->canonical){
            fatal_error("%s:  Graph must be in canonical form\n", __FUNCTION__);
        }

        //FIXME: this will surely lead to 64 bit issues
        if(g->key == 1 << 31){
            // The key is big. Reset the key and zero out the adj_vec
            g->key = 0;
            for(int i = 0; i < g->capacity; i++){
                g->adj_vec[i] = 0;
            }
        }

        g->key++;
        for(list<int>::iterator L = g->nodes[v].nbrs.begin(); L
            != g->nodes[v].nbrs.end(); ++L){
            g->adj_vec[*L] = g->key;
        }

        return g->key;
    } // fill_adj_vec

    /**
     * Adds edges as necessary so that the vertices form a clique.  Returns the number
     * of edges added to the graph.
     */
    int GraphProperties::make_clique(Graph *g, list<int> *vertices){
        int m = 0;
        list<int>::iterator ii, jj;

        for(ii = vertices->begin(); ii != vertices->end(); ++ii){
            int current_key = fill_adj_vec(g, *ii);
            for(jj = ii; ++jj != vertices->end(); ){
                if(g->adj_vec[*jj] != current_key){
                    // Add the edge *ii-*jj
                    m++;
                    g->add_edge(*ii, *jj);
                }
            }
        }

        // Return the # of edges added
        return m;
    } // make_clique

    /**
     * Returns true of the list of vertices forms a clique in the
     * graph, false otherwise.
     */
    bool GraphProperties::is_clique(Graph *g, list<int> *vertices){
        list<int>::iterator ii, jj;

        for(ii = vertices->begin(); ii != vertices->end(); ++ii){
            jj = vertices->begin();
            int current_key = fill_adj_vec(g, *ii);
            while(jj != vertices->end()){
                // CSG - ignoring self loops here!
                if((g->adj_vec[*jj] != current_key) && (*ii != *jj) ){
                    //				print_message(0, "Did not find required clique edge %d-%d\n",
                    //						*ii, *jj);
                    return false;
                }
                ++jj;
            }
        }

        return true;
    } // is_clique

    /**
     * Manually verify if current graph is simple - time intensive with data copies (checks adj lists for duplicates & loops).
     * Returns true or false.
     */
    bool GraphProperties::check_simple(Graph *g){
        int i;
        list<int> *tmp;
        size_t newdeg;

        for(i = 0; i < g->capacity; i++){
            //loop through the active nodes
            if(g->nodes[i].label != -1){
                //copy the list of neighbors
                tmp = new list<int> (g->nodes[i].nbrs);
                //sort it and remove duplicates and loops
                tmp->sort();
                tmp->unique();
                tmp->remove(i);
                newdeg = tmp->size();
                delete tmp;
                //check to see if anything changed & return false if it did (graph was not simple)
                if(newdeg != g->nodes[i].nbrs.size()){
                    return false;
                }
            }
        }

        return true;
    } // check_simple

    bool GraphProperties::is_connected(Graph *g){
        list<int> m;
        GraphUtil util;
        if(util.find_component(g, 0, &m) == g->num_nodes){
            g->num_connected_components = 1;
            return true;
        }
        else {
            return false;
        }
    }

    bool GraphProperties::is_independent_set(Graph *g, list<int> *vertices){
        /**
         * Returns true of the list of vertices forms an independent set in the
         * graph, false otherwise.
         */

        if(vertices->size() == 1){
            return true;
        }

        list<int>::iterator ii, jj;

        for(ii = vertices->begin(); ii != vertices->end(); ++ii){
            int current_key = fill_adj_vec(g, *ii);
            for(jj = ii; ++jj != vertices->end(); ){
                if(g->adj_vec[*jj] == current_key){
                    return false;
                }
            }
        }

        // Found no adjacent vertices
        return true;
    } // is_independent_set

    /**
     * Returns true of the list of vertices forms an independent set in the
     * graph, false otherwise. Sets val to be the weight of the set.
     */
    bool GraphProperties::is_independent_set(VertexWeightedGraph *wg,
                                             list<int> *vertices, int *val){
        if(vertices->size() == 1){
            return true;
        }

        list<int>::iterator ii, jj;

        for(ii = vertices->begin(); ii != vertices->end(); ++ii){
            int current_key = fill_adj_vec(wg, *ii);
            for(jj = ii; ++jj != vertices->end(); ){
                if(wg->adj_vec[*jj] == current_key){
                    return false;
                }
            }

            *val += wg->weight[*ii];
        }

        // Found no adjacent vertices
        return true;
    } // is_independent_set

    /**
     * Uses BFS to determine if there is a path between start and end in the
     * graph that passes through only those vertices k such that t[k]=true.
     */

    bool GraphProperties::is_path(Graph *g, int start, int end, bool *t){
        // Returns true if we find a path b/w start and end that uses only vertices v
        // such that t[v]=true, false o/w
        // We need the graph to be symmetric, as we're going to walk through a neighbor list
        if(!g->canonical){
            fatal_error("%s:  must be in canonical format\n", __FUNCTION__);
        }

        int j;
        bool *visited = new bool[g->capacity];
        for(j = 0; j < g->capacity; j++){
            visited[j] = false;
        }

        // S is a stack that will contain the nodes we have visited -
        // take the first one off, and then add its unvisited neighbors to the end
        // If the stack ever gets empty and we didn't see end, then return false
        // since there is no such path!
        list<int> S;
        list<int>::iterator ii;

        // Put start on the stack
        S.clear();
        S.push_front(start);
        while(!(S.empty())){
            // Remove the oldest element from the stack
            j = S.front();
            // See if j is the end
            if(j == end){
                delete[] visited;
                return true;
            }
            S.pop_front();
            if(visited[j] == false){
                // We have now visited j
                visited[j] = true;

                // Check j's neighbors
                for(ii = g->nodes[j].nbrs.begin(); ii != g->nodes[j].nbrs.end(); ++ii){
                    if(*ii == end){
                        delete[] visited;
                        return true;
                    }
                    if((visited[*ii] == false) && (t[*ii] == true) ){
                        // Note - g is the only place that we refer to the t[] vector
                        // We haven't seen *ii before and it is an "acceptable" vertex,
                        // so it is a candidate to be in the path - add it to the Stack
                        S.push_back(*ii);                         // used to be push_front - does it matter?
                    }
                }
            }
        }

        // We emptied the stack and never found the end --> no path!
        delete[] visited;
        return false;
    } // is_path

    /**
     * Uses BFS to determine if there is a path between start and end in the
     * graph.
     */
    bool GraphProperties::is_path(Graph *g, int start, int end){
        // Returns true if we find a path b/w start and end that uses only vertices v
        // such that t[v]=true, false o/w

        // We need the graph to be symmetric, as we're going to walk through a neighbor list
        if(!g->canonical){
            fatal_error("%s:  must be in canonical format\n", __FUNCTION__);
        }

        int j;
        bool *visited = new bool[g->capacity];
        for(j = 0; j < g->capacity; j++){
            visited[j] = false;
        }

        // S is a stack that will contain the nodes we have visited -
        // take the first one off, and then add its unvisited neighbors to the end
        // If the stack ever gets empty and we didn't see end, then return false
        // since there is no such path!
        list<int> S;
        list<int>::iterator ii;

        // Put v on the stack
        S.push_front(start);
        while(!(S.empty())){
            // Remove the oldest element from the stack
            j = S.front();
            S.pop_front();
            if(visited[j] == false){
                // We have now visited j
                visited[j] = true;

                for(ii = g->nodes[j].nbrs.begin(); ii != g->nodes[j].nbrs.end(); ++ii){
                    if(*ii == end){
                        delete[] visited;
                        return true;
                    }
                    if(visited[*ii] == false){
                        // We haven't seen *ii before and it is an "acceptable" vertex,
                        // so it is a candidate to be in the path - add it to the Stack
                        S.push_back(*ii);                         // used to be push_front - does it matter?
                    }
                }
            }
        }

        // We emptied the stack and never found the end --> no path!
        delete[] visited;
        return false;
    } // is_path

    /**
     * Equivalent to Algorithm 8 in Latapy 'new-vertex-listing'.  Increments elements in the triangle array.
     * \param[in] g the graph
     * \param[in] v the vertex we are listing triangles for
     * \param[in,out] tr the counts of our triangles
     */
    void GraphProperties::vertex_listing(Graph *g, int v, vector<long int> &t){
        vector<bool> A(g->get_num_nodes(), false);
        int u;
        int i;
        Node *n;
        list<int> *nbrs;
        list<int> *nbr_nbrs;
        list<int>::iterator inner, outer;
        list<int>::const_reverse_iterator crit, stuff;

        fprintf(stderr,"vertex_listing(%d)\n", v);

        nbrs = g->get_node(v)->get_nbrs_ptr();

        // set A[u] true for all neighbors of v
        for(outer = nbrs->begin(); outer != nbrs->end(); ++outer){
            A[*outer] = true;
        }

        fprintf(stderr,"  vertex_listing v->nbrs:");
        for(inner = nbrs->begin(); inner != nbrs->end(); ++inner){
            fprintf(stderr, " %d", *inner);
        }
        fprintf(stderr, "\n");

        for(outer = nbrs->begin(); outer != nbrs->end(); ++outer){
            n = g->get_node(*outer);
            nbr_nbrs = n->get_nbrs_ptr();

            fprintf(stderr,"    vertex_listing %d->nbrs:", *outer);
            for(inner = nbr_nbrs->begin(); inner != nbr_nbrs->end(); ++inner){
                fprintf(stderr, " %d", *inner);
            }
            fprintf(stderr, "\n");

            //    for(inner = nbr_nbrs->begin(); inner != nbr_nbrs->end(); ++inner){
            stuff = n->get_largest_neighbor_below(*outer);
            fprintf(stderr, "Found closest to %d as %d\n", *outer, *stuff);
            for(crit = stuff; crit != nbr_nbrs->rend(); crit--){
                if(*crit != v){
                    fprintf(stderr, "      vertex_listing: checking triangle (%d,%d,%d)\n", v,*crit,*outer);
                    if(A[*crit] == true){
                        t[*crit]++;
                        t[*outer]++;
                        t[v]++;
                        fprintf(stderr,"      vertex_listing: ---------Found triangle (%d,%d,%d)\n",v,*crit,*outer);
                    }
                }
            }
        }
    } //vertex_listing

    /**
     * Equivalent to Algorithm 3 in Latapy 'edge-iterator'.  Increments elements in the triangle array.
     * \param[in] g the graph
     * \param[in] u first end of the we are listing triangles for
     * \param[in] v second end of the we are listing triangles for
     * \param[in,out] tr the counts of our triangles
     */
    void GraphProperties::edge_listing(Graph *g, int u, int v, vector<long int> &t, int number_high){
        vector<int>::iterator it;
        vector<int> intersection;

        fprintf(stderr,"edge_listing(%d-%d)\n", u, v);

        list<int> *un = g->get_node(u)->get_nbrs_ptr();
        list<int> *vn = g->get_node(v)->get_nbrs_ptr();

        // take adavantage of STL stuff
        // back_inserter allows us to just declare a vector and not pre-determine size, saves space
        // can remove and make vector the size of the smaller neighbor list if we find it's too slow
        std::set_intersection(un->begin(), un->end(), vn->begin(), vn->end(), std::back_inserter(intersection));

        for(it = intersection.begin(); it < intersection.end(); ++it){
            if(*it >= number_high){
                t[*it]++;
            }
        }
    } // edge_listing

    //FIXME: figure out where this should really live
    struct sort_pair {
        bool operator()(const pair<int, int> &left, const pair<int, int> &right){
            return left.second > right.second;
        }
    };

    void GraphProperties::all_triangles(Graph *g, vector<long int> &t, int number_high){
        int i, j;
        int u, v;
        int retcode;
        Node *vn;

        std::list<int>::const_reverse_iterator rit;

        if(number_high > g->get_num_nodes()){
            number_high = g->get_num_nodes();
        }

        vector<pair <int, int> > sorted_indices(g->get_num_nodes());

        // we want our list of vertices sorted by degree, with higest degree in element 0
        // this is a goofy way to handle it, but that's life
        for(i = 0; i < g->get_num_nodes(); i++){
            sorted_indices[i].first = i;
            sorted_indices[i].second = g->get_degree(i);
        }

        fprintf(stderr, "Before:\n");
        for(i = 0; i < g->get_num_nodes(); i++){
            fprintf(stderr, " vertex: %d degree %d\n",i,g->get_degree(i));
        }

        std::sort(sorted_indices.begin(), sorted_indices.end(), sort_pair());
        fprintf(stderr, "After:\n");
        for(i = 0; i < g->get_num_nodes(); i++){
            fprintf(stderr, " vertex: %d(%d) degree %d\n",i,sorted_indices[i].first, g->get_degree(sorted_indices[i].first));
        }

        // we need sorted neighbor lists for edge_listing
        for(i = 0; i < g->get_num_nodes(); i++){
            g->get_node(i)->sort_nbr();  //FIXME: make sure this actually does what i think it does -jkl
            g->get_node(i)->reverse_nbr();
        }

        // count triangles using vertex_listing for 'high degree' vertices (1a)
        for(i = 0; i < number_high; i++){
            vertex_listing(g, sorted_indices[i].first, t);
        }

        for(v = g->get_num_nodes() - 1; v >= number_high; v--){
            vn = g->get_node(sorted_indices[v].first);
            for(rit = vn->get_largest_neighbor_below(v); rit != vn->get_nbrs_ptr()->rbegin(); rit--){
                if((g->get_degree(*rit) >= number_high)){
                    edge_listing(g, *rit, sorted_indices[v].first, t, number_high);
                }
            }
        }
        // count triangles using edge_listing for 'low degree' vertices (2)
    } // all_triangles

    void GraphProperties::all_triangles_compact_forward(Graph *g, vector<long int> &t){
        int i, j;
        int u, v;
        int retcode;
        Node *vn;

        std::list<int>::const_reverse_iterator rit;

        vector<pair <int, int> > sorted_indices(g->get_num_nodes());

        // we want our list of vertices sorted by degree, with higest degree in element 0
        // this is a goofy way to handle it, but that's life
        for(i = 0; i < g->get_num_nodes(); i++){
            sorted_indices[i].first = i;
            sorted_indices[i].second = g->get_degree(i);
        }

        // here we've basically renumbered the array using the
        // injective function as required by Algorithm 7 in Latapy
        fprintf(stderr, "Before:\n");
        for(i = 0; i < g->get_num_nodes(); i++){
            fprintf(stderr, " vertex: %d degree %d\n",i,g->get_degree(i));
        }

        std::sort(sorted_indices.begin(), sorted_indices.end(), sort_pair());
        fprintf(stderr, "After:\n");
        for(i = 0; i < g->get_num_nodes(); i++){
            fprintf(stderr, " vertex: %d(%d) degree %d\n",i,sorted_indices[i].first, g->get_degree(sorted_indices[i].first));
        }

        // we need sorted neighbor lists for edge_listing
        for(i = 0; i < g->get_num_nodes(); i++){
            g->get_node(i)->sort_nbr();  //FIXME: make sure this actually does what i think it does -jkl
        }

        list<int> *nbrs;
        Node *nv, *nu;
        int vp, up;
        for(v = 0; v < g->get_num_nodes(); v++){ //3
            nv = g->get_node(v);
            for(std::list<int>::iterator it = nbrs->begin(); it != nbrs->end(); it++){ //3a
                nu = g->get_node(*it);
                if(nu->get_degree() > nv->get_degree()){ //3a
                }
            }
        }
    } // all_triangles_compact_forward

    /*
     * Counts triangles using the edge-listing algorithm in Latapy
     * \param[in] g input graph
     * \param[out] t vector of long ints, length |V|, returns 3x number of triangles for each vertex
     */
    void GraphProperties::all_triangles_edge_listing(Graph *g, vector<long int> &t){
        int i, j, u, v;
        std::list<int> *u_n, *v_n;
        vector<int>::iterator it;
        list<int>::const_iterator cit;
        list<int>::iterator lt;

        // all the edgelists must be sorted
        for(i = 0; i < g->get_num_nodes(); i++){
            g->get_node(i)->sort_nbr();
        }

        for(u = 0; u < g->get_num_nodes(); u++){
            const list<int> &c_nbrs = g->get_node(u)->get_nbrs_ref();
            for(cit = c_nbrs.begin(); cit != c_nbrs.end(); ++cit){
                v = *cit;
                if(v > u){
                    vector<int> intersection;
                    u_n = g->get_node(u)->get_nbrs_ptr();
                    v_n = g->get_node(v)->get_nbrs_ptr();
                    //printf("looking at edge: %d-%d\n", u, v);
                    std::set_intersection(u_n->begin(), u_n->end(), v_n->begin(), v_n->end(), std::back_inserter(intersection));
                    /*
                       printf("  nbrs_u(%d): ", u);
                       for(lt = u_n->begin(); lt != u_n->end(); ++lt){
                        printf(" %d", *lt);
                       }
                       //printf("\n  nbrs_v(%d): ", v);
                       for(lt = v_n->begin(); lt != v_n->end(); ++lt){
                        printf(" %d", *lt);
                       }
                       //printf("\n    Intersection: ");
                       for(it = intersection.begin(); it != intersection.end(); it++){
                        printf(" %d", *it);
                       }
                       printf("\n");
                     */

                    for(it = intersection.begin(); it != intersection.end(); it++){
                        //printf("found triangle: (%d,%d,%d)\n",*it, u, v);
                        t[*it]++;
                        t[u]++;
                        t[v]++;
                    }
                }
            }
        }
    } // all_triangles_edge_listing
}
