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
#include <numeric>

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

    //FIXME: figure out where this should really live
    struct sort_pair {
        bool operator()(const pair<int, int> &left, const pair<int, int> &right){
            return left.second > right.second;
        }
    };

    /**
     * Counts all triangles in g using compact-forward algorithm from Latapy
     * \param[in] g input graph
     * \param[out] t output vector of per-vertex triangle counts
     */
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
            /*
               fprintf(stderr, "  nbrs(%d): ", i);
               for(list<int>::iterator myit=g->get_node(i)->get_nbrs_ptr()->begin(); myit != g->get_node(i)->get_nbrs_ptr()->end(); ++myit){
                fprintf(stderr, " %d(%d)", *myit, g->get_degree(*myit));
               }
               fprintf(stderr,"\n");
             */
        }

        // here we've basically renumbered the array using the
        // injective function as required by Algorithm 7 in Latapy
        /*
           fprintf(stderr, "Before:\n");
           for(i = 0; i < g->get_num_nodes(); i++){
            fprintf(stderr, " vertex: %d degree %d\n",i,g->get_degree(i));
           }
         */
        std::sort(sorted_indices.begin(), sorted_indices.end(), sort_pair());
        /*
           fprintf(stderr, "After:\n");
           for(i = 0; i < g->get_num_nodes(); i++){
            fprintf(stderr, " vertex: %d(%d) degree %d\n",i,sorted_indices[i].first, g->get_degree(sorted_indices[i].first));
           }
         */

        Node *nv, *nu;
        int vp, up;
        int iu, iv;
        int fakev;
        int uprime, vprime;
        const int n = g->get_num_nodes();
        // FIXME: this is a terrible hack
        //
        vector<int> revmap(n, -1);
        for(i = 0; i < n; i++){
            revmap[sorted_indices[i].first] = n - (i + 1);
        }
        for(i = 0; i < n; i++){
            sort_nbrs_by_map(revmap, g->get_node(i));
        }

        for(fakev = 1; fakev < n; fakev++){   //3
            v = sorted_indices[fakev].first;
            //fprintf(stderr,"fakev: %d v: %d\n",fakev, v);
            nv = g->get_node(v);
            const list<int> &nbrs_v = nv->get_nbrs_ref();
            for(std::list<int>::const_iterator cit = nbrs_v.begin(); cit != nbrs_v.end(); ++cit){ //3a
                u = *cit;
                //fprintf(stderr,"    looking at u=%d v=%d\n", u, v);
                //fprintf(stderr,"      n(u):%d n(v):%d\n",revmap[u],revmap[v]);
                if(revmap[u] > revmap[v]){ //3a
                    nu = g->get_node(u);
                    const list<int> &nbrs_u = nu->get_nbrs_ref();
                    list<int>::const_iterator upit = nbrs_u.begin(); //3aa
                    list<int>::const_iterator vpit = nbrs_v.begin(); //3aa
                    uprime = *upit;
                    vprime = *vpit;
                    //fprintf(stderr, "       u':%d v':%d ",uprime, vprime);
                    //fprintf(stderr, "n(v):%d n(u'):%d n(v'):%d\n", revmap[v], revmap[uprime], revmap[vprime]);
                    while((upit != nbrs_u.end()) && (vpit != nbrs_v.end() &&
                                                     (revmap[uprime] < revmap[v]) && (revmap[vprime] < revmap[v]))){ //3ab
                        uprime = *upit;
                        vprime = *vpit;
                        //fprintf(stderr, "       u':%d v':%d ",uprime, vprime);
                        //fprintf(stderr, "n(v):%d n(u'):%d n(v'):%d\n", revmap[v], revmap[uprime], revmap[vprime]);
                        if(revmap[uprime] < revmap[vprime]){ //3aba
                            upit++;
                            //fprintf(stderr, "        Failed triangle: (%d,%d,%d) incrementing u'\n", v, u, uprime);
                        }
                        else if(revmap[uprime] > revmap[vprime]){ //3abb
                            vpit++;
                            //fprintf(stderr, "        Failed triangle: (%d,%d,%d) incrementing v'\n", v, u, uprime);
                        }
                        else {   //3abc
                            //fprintf(stderr, "Found triangle: (%d,%d,%d)\n", v, u, uprime);
                            t[v]++;
                            t[u]++;
                            t[uprime]++;
                            upit++;
                            vpit++;
                        }
                    }
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
                    const std::list<int> &u_n = g->get_node(u)->get_nbrs_ref();
                    const std::list<int> &v_n = g->get_node(v)->get_nbrs_ref();
                    vector<int> intersection (max(u_n.size(), v_n.size()));
                    //printf("looking at edge: %d-%d\n", u, v);
                    it = std::set_intersection(u_n.begin(), u_n.end(), v_n.begin(), v_n.end(), intersection.begin());
                    intersection.resize(it - intersection.begin());
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
                        t[*it]++;
                        t[u]++;
                        t[v]++;
                    }
                }
            }
        }
    } // all_triangles_edge_listing

    /**
     * Functor to compare degrees and sort by them
     */
    struct compare_inject : std::binary_function<size_t, size_t, bool>{
        compare_inject(const std::vector<int> &inject) : m_inject(inject){
        }

        bool operator()(size_t left, size_t right) const {
            return m_inject[left] < m_inject[right];
        }

        const std::vector<int>& m_inject;
    };

    void GraphProperties::sort_nbrs_by_map(vector<int> map, Node *n){
        list<int> nbrs(n->get_nbrs());

        //std::sort(nbrs.first(), nbrs.last(), compare_degrees(g->get_degree_ref()));
        nbrs.sort(compare_inject(map));
        n->set_nbr(nbrs);
    }

    /**
     * Returns the clustering coefficients of g.  Uses compact-forward under the hood for now
     * \param[in] g input graph
     * \param[out] global_cc global clustering coefficient, defined as (3*num_triangles)/num_possible_triangles
     * \param[out] avg_cc average clustering coefficient, defined as the average of the local clustering coefficients
     * \param[out] local_ccs local clustering coeffients
     */

    void GraphProperties::clustering_coefficients(Graph *g, double &global_cc, double &avg_cc, vector<double> &local_ccs){
        int i, d_i;
        const int n = g->get_num_nodes();
        vector<long int> triangles(n, 0);
        int total_possible_triangles = 0;

        all_triangles_compact_forward(g, triangles);

        global_cc = 0;
        avg_cc = 0;
        local_ccs.resize(n, 0);

        const int num_triangles = std::accumulate(triangles.begin(), triangles.end(), 0);

        // Calculate local CCs and running sum for global cc calc later
        for(i = 0; i < n; i++){
            d_i = g->get_degree(i);
            total_possible_triangles += (d_i * (d_i - 1) / 2);
            if(triangles[i] > 0){
                local_ccs[i] = 2.0 * triangles[i] / ((double)d_i * (double)(d_i - 1));
            }
        }

        avg_cc = std::accumulate(local_ccs.begin(), local_ccs.end(), 0.0) / (double)n;
        global_cc = (double)num_triangles / (double)total_possible_triangles;
    } // clustering_coefficients

    //==================================================================
    void GraphProperties::paths_dijkstra_single(Graph *g, vector<int> &p, int source){
        int inf = 100000;
        int nVisited = 0;
        int minD;
	int nvisiting;
        const int n = g->get_num_nodes(); //number of nodes in graph
        
        //initialize empty set of vertices
        vector<int> dist(n, inf);
        vector<int> visited(n, 0);
	Node *nv;

        /*
        printf("Initial p is:\n");
        for(int i = 0; i < n; i++){
	    printf("p[i] = %d ;",p[i]);

            nv = g->get_node(i); //this is the node, not the ID!
            printf(" Node %d with %d neighbors \n",nv->label, nv->nbrs.size());
            const list<int> &mynbrs = nv->get_nbrs_ref();

            printf("\tneighbors are ");
            for(list<int>::const_iterator cit = mynbrs.begin(); cit != mynbrs.end(); ++cit){
                printf("%d",*cit);
            }
            printf("\n");
        }
	
        */
        
        //initialize
        dist[source] = 0;
	nvisiting = 0;

        //for(int i = 0; i < (int) dist.size(); i++){
        //    cout << dist[i] << " ";
        //}

        while(nVisited < n){
            //get list of neighbors
            nv = g->get_node(nvisiting);  //"true" name is tmp->label
            const list<int> &mynbrs = nv->get_nbrs_ref();

            //for each neighbor of the current vertex considered
            for(list<int>::const_iterator nb = mynbrs.begin(); nb != mynbrs.end(); ++nb){
                //consider each neighbor not already visited and "relax" dist
                if(!visited[*nb]){
                    dist[*nb] = min(dist[nvisiting] + 1, dist[*nb]); //weight=1 b/c unweighted graph
                }
            }

            //mark node as visited
            visited[nvisiting] = 1;

            //find next move
            minD = inf;
            for(int i = 0; i < n; i++){
                if(!visited[i] && (dist[i] < minD)){
                    minD = dist[i];
                    nvisiting = i;
                }
            }

	    /*
            printf("After pass d\n\t");
            for(int i = 0; i < (int) dist.size(); i++){
                cout << dist[i] << " ";
            }
            printf("\n\t");
            for(int i = 0; i < (int) visited.size(); i++){
                cout << visited[i] << " ";
            }
            printf("\n");
            printf("*****Next to visit is %d******\n",source);
	    */

            nVisited++;
        }

        for(int i=0;i<n;i++) {
	    printf("%d,  ",dist[i]);
	}
	printf("\n");
    } //paths_dijkstra_single
  
 void GraphProperties::paths_dijkstra_all(Graph *g, vector< vector<int> > &p){
        int inf = 100000; //assume this is larger than any weights on the graph
        int nVisited = 0;
        int minD = inf;
        const int n = g->get_num_nodes(); 
        int nvisiting;
	
        Node *nv;
        vector<int> dist(n, inf);
        vector<int> visited(n, 0);

	//loop over all vertices
        for (int v = 0; v < n; v++) {
	    //printf("\n***Source is now v %d \n ",v);
	    nvisiting = v;        
	    nVisited = 0;   

            //reset all distances to INF and mark all vertices as unvisited
	    fill(dist.begin(),dist.end(),inf);
	    fill(visited.begin(),visited.end(),0);

            //initialize
            dist[nvisiting] = 0;

            //loop until all have been marked visited
            while(nVisited < n){
                //get list of neighbors
                nv = g->get_node(nvisiting);  //"true" name is tmp->label
                const list<int> &mynbrs = nv->get_nbrs_ref();

                //for each neighbor of the current vertex considered
                for(list<int>::const_iterator nb = mynbrs.begin(); nb != mynbrs.end(); ++nb){
                    //consider each neighbor not already visited and "relax" dist
                    if(!visited[*nb]){
                        dist[*nb] = min(dist[nvisiting] + 1, dist[*nb]); //weight=1 b/c unweighted graph
                    }
                }

                //mark node as visited
                visited[nvisiting] = 1;

                //find next vertex to move to
                minD = inf;
                for(int i = 0; i < n; i++){
                    if(!visited[i] && (dist[i] < minD)){
                        minD = dist[i];
                        nvisiting = i;
                    }
                }
                
                nVisited++;
	    }//end while

    	    //store shortest paths from this vertex to all
	    p.push_back(dist);
	    //printf(" with distances: ");
	    for(int i=0;i<n;i++) {
	      printf("%d,  ",dist[i]);
	    }
	    printf("\n");
	    
        }//end loop over vertices
	
    } //paths_dijkstra_all
}
