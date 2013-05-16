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


    void GraphProperties::edgeDensity(Graph *g, float &ed){
        int V = g->num_nodes;
        int E = g->num_edges;
        
        ed = (2.0*E)/(V*(V-1.0));
    }


    void GraphProperties::graphSize(Graph *g, int &size){
        size = g->num_nodes;
    }


    void GraphProperties::degDist(Graph *g, vector<int> &dist){
        int i, size=g->degree.size();

        dist.resize(size);
        for(i=0; i<size; i++){
            dist[i] = 0;
        }

        for(i=0; i<size; i++){
            dist[g->degree[i]] += 1;
        }
    }
}
