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
#include "Log.h"
#include <sys/types.h>
#if !WIN32
  #include <sched.h>
#endif

#include <algorithm>
#include <vector>
#include <cmath>

#include <numeric>
#ifdef HAS_SLEPC
  #include <slepceps.h>
#endif

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

        #pragma omp parallel for schedule(dynamic, 8) default(none) shared(g)
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
        int i;
        int u, v;

        std::list<int>::const_reverse_iterator rit;

        vector<pair <int, int> > sorted_indices(g->get_num_nodes());

        // we want our list of vertices sorted by degree, with higest degree in element 0
        // this is a goofy way to handle it, but that's life
        #pragma omp parallel for default(none) shared(g, sorted_indices)
        for(i = 0; i < g->get_num_nodes(); i++){
            sorted_indices[i].first = i;
            sorted_indices[i].second = g->get_degree(i);
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
        int fakev;
        int uprime, vprime;
        const int n = g->get_num_nodes();
        // FIXME: this is a terrible hack
        //
        vector<int> revmap(n, -1);

        vector< vector<long int> > local_t;

        #pragma omp parallel default(none) shared(sorted_indices, revmap, g, t, local_t) private(fakev,nv,v,u,uprime,vprime,nu)
        {
            #pragma omp single
            {
                local_t.resize(omp_get_num_threads());
            }
            local_t[omp_get_thread_num()].resize(t.size());
            //pid_t tid = syscall(SYS_gettid);
            //cpu_set_t cpuset;
            //printf("OMP thread %d has Linux TID %ld running on CPU: %d\n", omp_get_thread_num(), tid, sched_getcpu());
            //sched_getaffinity(tid, sizeof(cpu_set_t), &cpuset);
            //for(int k=0; k<47; k++){
            //    if(CPU_ISSET(k, &cpuset)){
            //        printf(" %d",k);
            //    }
            // }
            // printf("\n");

            #pragma omp for schedule(dynamic, 8)
            for(i = 0; i < n; i++){
                revmap[sorted_indices[i].first] = n - (i + 1);
            }

            #pragma omp for  schedule(dynamic, 64)
            for(i = 0; i < n; i++){
                sort_nbrs_by_map(revmap, g->get_node(i));
            }

            int omp_tid = omp_get_thread_num();
            #pragma omp for schedule(dynamic, 4)
//#pragma omp parallel for  default(none) schedule(dynamic, 16) shared(t, revmap, sorted_indices, g) private(fakev,nv,v,u,uprime,vprime,nu)
            for(fakev = 1; fakev < n; fakev++){ //3
                v = sorted_indices[fakev].first;
                //fprintf(stderr,"fakev: %d v: % schedule(dynamic, 8)d\n",fakev, v);
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
                            else { //3abc
                                   //fprintf(stderr, "Found triangle: (%d,%d,%d)\n", v, u, uprime);
                                local_t[omp_tid][v]++;
                                local_t[omp_tid][u]++;
                                local_t[omp_tid][uprime]++;
                                upit++;
                                vpit++;
                            }
                        }
                    } //if revmap
                } //for vtxs
            } //for fakev
            int tsize = (int)t.size();
            #pragma omp for
            for(i = 0; i < tsize; i++){
                for(int j = 0; j < omp_get_num_threads(); j++){
                    t[i] += local_t[j][i];
                }
            }
        } //parallel
    } // all_triangles_compact_forward

    /*
     * Counts triangles using the edge-listing algorithm in Latapy
     * \param[in] g input graph
     * \param[out] t vector of long ints, length |V|, returns 3x number of triangles for each vertex
     */
    void GraphProperties::all_triangles_edge_listing(Graph *g, vector<long int> &t){
        int i, u, v;
        vector<int>::iterator it;
        list<int>::const_iterator cit;
        list<int>::iterator lt;

        // all the edgelists must be sorted
        #pragma omp parallel for schedule(dynamic, 8) default(none) shared(g)
        for(i = 0; i < g->get_num_nodes(); i++){
            g->get_node(i)->sort_nbr();
        }

        #pragma omp parallel for schedule(dynamic, 8) default(none) shared(g,t) private(u, v, it, cit)
        for(u = 0; u < g->get_num_nodes(); u++){
            const list<int> &c_nbrs = g->get_node(u)->get_nbrs_ref();
//#pragma omp parallel for schedule(dynamic, 8) default(none) shared(g,t) private(u, v, it, cit)
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
                        #pragma omp atomic
                        t[*it]++;
                        #pragma omp atomic
                        t[u]++;
                        #pragma omp atomic
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
        list<int> *nbrs(n->get_nbrs_ptr());

        //std::sort(nbrs.first(), nbrs.last(), compare_degrees(g->get_degree_ref()));
        nbrs->sort(compare_inject(map));
        //n->set_nbr(*nbrs);
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

    #ifdef HAS_BOOST
    /**
     * We assume that edge weights are all 1
     *
     * \param[in] g input graph
     * \param[in] source node id
     * \param[out] p path distances from source to all other vertices
     */
    void GraphProperties::paths_dijkstra_boost_single(Graph *g, vector<int> &dists, vector<vertex_descriptor> &preds, int source){
        BoostUndirected *bg = g->boost_graph;
        //std::vector<vertex_descriptor> p(boost::num_vertices(*bg));
        dists.resize(g->get_num_nodes());
        //std::vector<int> d(boost::num_vertices(*bg));
        vertex_descriptor s = vertex(source, *bg);
        boost::dijkstra_shortest_paths(*bg, s, boost::predecessor_map(&preds[0]).distance_map(&dists[0]).distance_inf(INDDGO_INFINITY));
        int i;
        //for(i=0;i<preds.size();i++){
        //    cout << "node " << i << " pred " << preds[i] << "\n";
        //}
    }

    /**
     *  All pairs shortest paths
     * \param[in] g input graph
     * \param[out] p multidimentional list of all pairs shortest paths
     */

    void GraphProperties::paths_dijkstra_boost_all(Graph *g, vector< vector<int> > &pAll){
        int inf = INDDGO_INFINITY;
        int minD = inf;

        const int n = g->get_num_nodes();

        std::vector<uint64_t> betweenness_counts(n, 0);
        std::vector<double> betweenness(n, 0.0);

        pAll.resize(n);

        //#pragma omp parallel for default(none) shared(g, inf, pAll) private(nvisiting, nVisited, nv) firstprivate(dist, minD, visited)
        #pragma omp parallel default(none) shared(g, inf, betweenness_counts, pAll)
        {
            std::vector<vertex_descriptor> p(boost::num_vertices(*(g->boost_graph)));
            int i;
            #pragma omp for schedule(dynamic, 8)
            for(int v = 0; v < n; v++){
                //reset all distances to INF and mark all vertices as unvisited
                fill(pAll[v].begin(),pAll[v].end(),inf);
                paths_dijkstra_boost_single(g, pAll[v], p, v); //stores shortest paths from this vertex to all in pAll[v]
            }
        }
        //store the results
        g->set_shortest_path_dist(&pAll);
    } // paths_dijkstra_boost_all

    /**
     * Calculate the relative betweenness centrality for all nodes
     * \param[in] g input graph
     * \param[out] bc vector<double> with one BC value for each vertex
     */

    void GraphProperties::betweenness_centrality(Graph *g, vector<double> &bc){
        BoostUndirected *bg = g->boost_graph;
        bc.resize(g->get_num_nodes());
        //boost::brandes_betweenness_centrality(*bg, get(boost::vertex_centrality, *bg));
        boost::brandes_betweenness_centrality(*bg,
                                              boost::centrality_map(
                                                  boost::make_iterator_property_map(bc.begin(), get(boost::vertex_index, *bg), double())).vertex_index_map(get(boost::vertex_index, *bg)
                                                                                                                                                           )
                                              );

        boost::relative_betweenness_centrality(*bg,
                                               boost::make_iterator_property_map(bc.begin(), get(boost::vertex_index, *bg), double()));
        g->set_betweenness(bc);
    }

    #endif // ifdef HAS_BOOST

    /**
     * \param[in] g input graph
     * \param[in] source node id
     * \param[out] p path distances from source to all other vertices
     */

    void GraphProperties::paths_dijkstra_single(Graph *g, vector<int> &p, int source){
        int inf = INDDGO_INFINITY;
        int nVisited = 0;
        int minD;
        int nvisiting;
        const int n = g->get_num_nodes(); //number of nodes in graph

        //initialize empty set of vertices
        vector<int> dist(n, inf);
        vector<int> visited(n, 0);
        Node *nv;

        //initialize
        dist[source] = 0;
        nvisiting = source;

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

            nVisited++;
        }

        //make sure to pass back correct values
        p = dist;

        /*
           printf("In single: \n");
           for(int i = 0; i < n; i++){
            printf("%d,  ",dist[i]);
           }
           printf("\n");
         */
    } //paths_dijkstra_single

    /**
     *  All pairs shortest paths
     * \param[in] g input graph
     * \param[out] p multidimentional list of all pairs shortest paths
     */

    void GraphProperties::paths_dijkstra_all(Graph *g, vector< vector<int> > &pAll){
        int inf = INDDGO_INFINITY;

        int minD = inf;
        const int n = g->get_num_nodes();

        pAll.resize(n);

        //#pragma omp parallel for default(none) shared(g, inf, pAll) private(nvisiting, nVisited, nv) firstprivate(dist, minD, visited)
        #pragma omp parallel for schedule(dynamic, 8)  default(none) shared(g, inf, pAll)
        //loop over all vertices
        for(int v = 0; v < n; v++){ //0; v < n; v++){
            //reset all distances to INF and mark all vertices as unvisited
            fill(pAll[v].begin(),pAll[v].end(),inf);
            paths_dijkstra_single(g, pAll[v], v); //stores shortest paths from this vertex to all in pAll[v]
        } //end loop over vertices

        //store the results
        g->set_shortest_path_dist(&pAll);

        //print out results
        //for(int i = 0; i < n; i++){
        //   for (int j = 0; j < n; j++) {
        //        printf("%d,  ",pAll[i][j]);
        //   }
        //    printf("\n");
        //}
        //printf("\n");
    } //paths_dijkstra_all

    /**
     * Calcuates the eccentricity for each vertex (max dist to any other vertex)
     * \param[in] g input graph
     * \param[out] ecc vector of eccentricies
     */
    void GraphProperties::eccentricity(Graph *g, vector<int> &ecc){
        const int n = g->get_num_nodes();
        const vector< vector<int> > &short_paths = g->get_shortest_path_dist_ref();   //computes if not already
        int bestMax = 0;
        ecc.resize(n);

        #pragma omp parallel for default(none) shared(short_paths, ecc) private(bestMax)
        //compute diameter of each vertex
        for(int i = 0; i < n; i++){
            bestMax = 0;
            for(int j = 0; j < n; j++){
                if(short_paths[i][j] > bestMax){
                    bestMax = short_paths[i][j];
                }
            }
            ecc[i] = bestMax;       //eccentricity of vertex i
            //printf("Eccentricity of node %d is %d\n", i, ecc[i]);
        }
    } //eccentricity

    /**
     * Calcuates the frequency of each eccentricity value for all vertices [ref. Takes, Kostes 2013]
     * \param[in] g input graph
     * \param[in] ecc vector of eccentricities (empty or pre-computed)
     * \param[out] freq_ecc vector of eccentricity frequencies
     */
    void GraphProperties::eccentricity_dist(Graph *g, vector<int> &ecc, vector<double> &freq_ecc){
        const int n = g->get_num_nodes();
        int bestAll = 0;

        if(ecc.empty()){
            cout << "Empty -- calling function to compute eccentricities" << endl;
            eccentricity(g,ecc);
        }

        freq_ecc.resize(*(std::max_element(ecc.begin(), ecc.end())) + 1);

        //compute diameter of each vertex
        int eccsize = ecc.size();
        for(int i = 0; i < eccsize; i++){
            freq_ecc[ecc[i]]++; //add to tally for this diameter size
        }
        //printf("Graph diameter is %d\n", freq_ecc.size()-1);

        int freq_ecc_size = freq_ecc.size();
        #pragma omp parallel for default(none) shared(freq_ecc, freq_ecc_size)
        for(int i = 0; i <= freq_ecc_size - 1; i++){
            freq_ecc[i] = freq_ecc[i] / n;
            //printf("i=%d and n=%d with freq eccentricity %f\n",i,n,freq_ecc[i]);
        }
    } //eccentricity_dist

    /**
     * Calcuates the expansion (hop-distance) from each vertex: [see Tangmunarunkit (2002)]
     * compute #nodes reachable in h hops starting at each vertex, average, normalize
     * \param[in] g input graph
     * \param[out] norm_hops normalized distribution
     */
    void GraphProperties::expansion(Graph *g, vector<double> &norm_hops){
        const int n = g->get_num_nodes();
        const vector< vector<int> > &short_paths = g->get_shortest_path_dist_ref();   //computes if not already
        vector <int> hops(n,0);  //largest possible dist is n-1
        norm_hops.resize(n);
        int k;

        //for each vertex, find the number of vertices reachable for all hops; add to tally
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                k = short_paths[i][j];
                if(k != INDDGO_INFINITY){
                    hops[k]++;
                }
            }
        }

        //average (divide by n) and normalize (divide by n-1)
        norm_hops[0] = 0.0; //no one is reachable in 0 hops (not counting self)

        #pragma omp parallel for default(none) shared(norm_hops, hops)
        for(int h = 1; h < n; h++){
            norm_hops[h] = (double)hops[h] / ((double)n * (n - 1));
            //printf("h = %d and number is %d; norm value is %f\n",h,hops[h],norm_hops[h]);
        }
    } //expansion

    /**
     * Calculates the diameter of graph g
     * \param[in] g Pointer to a graph
     * \param[out] diam Floating point value holding the calculated diamter
     */
    void GraphProperties::diameter(Graph *g, int &diam){
        int i, j, size, temp;
        const vector< vector<int> > &dist_mat = g->get_shortest_path_dist_ref();   //computes if not already

        size = dist_mat.size();
        diam = dist_mat[0][0];
        for(i = 0; i < size; i++){
            for(j = 0; j < i; j++){
                temp = dist_mat[i][j];

                if((temp > diam) && (temp < INDDGO_INFINITY)){
                    diam = temp;
                }
            }
        }
    } // diameter

    /**
     * Calculates the effective diameter of graph g
     * \param[in] g Pointer to a graph
     * \param[out] ediam Floating point value holding the calculated effective diameter
     * \param[in] perc Percentage for distances over total connected pairs, defaults to 0.9
     */
    void GraphProperties::effective_diameter(Graph *g, float &ediam, float perc){
        int i, j, d, size, temp, diam = 0;
        int n0, numer, tot_con_pairs = 0;
        float gd;
        const vector< vector<int> > &dist_mat = g->get_shortest_path_dist_ref();   //computes if not already
        vector<int> bins (g->num_nodes, 0);

        size = dist_mat.size();
        for(i = 0; i < size; i++){
            for(j = 0; j < i; j++){
                temp = dist_mat[i][j];
                if(temp < INDDGO_INFINITY){
                    tot_con_pairs++;
                    if(temp > diam){
                        diam = temp;
                    }
                    bins[temp]++;
                }
            }
        }

        bins.resize(diam + 1);
        numer = 0;
        gd = 0.0;
        d = 0;
        while(d <= diam && gd < perc){
            d++;
            n0 = numer;
            numer += bins[d];
            gd = numer / float(tot_con_pairs);
        }

        ediam = (tot_con_pairs * 0.9 - n0) / float(numer - n0) + d - 1;
    } // effective_diameter

    /**
     * Calculates the edge density of graph g
     * \param[in] g Pointer to a graph
     * \param[out] ed Floating point value holding the calculated edge density
     */
    void GraphProperties::edge_density(Graph *g, float &ed){
        int V = g->num_nodes;
        int E = g->num_edges;

        ed = (2.0 * E) / (V * (V - 1.0));
    }

    /**
     * Calculates the average degree of graph g
     * \param[in] g Pointer to a graph
     * \param[out] ad Floating point value holding the average degree of the graph
     */
    void GraphProperties::avg_degree(Graph *g, float &ad){
        int V = g->num_nodes;
        int E = g->num_edges;

        ad = (2.0 * E) / V;
    }

    void GraphProperties::avg_path_length(Graph *g, double &pl){
        const uint64_t n = g->get_num_nodes();
        const vector< vector<int> > &short_paths = g->get_shortest_path_dist_ref();
        double sum = 0;
        uint64_t intermediate = 0;
        int i, j;
        int inf_path = 0;

        #pragma omp parallel for default(none) reduction(+:inf_path) private(j) shared(short_paths)
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                if(INDDGO_INFINITY == short_paths[i][j]){
                    inf_path++;
                }
            }
        }

        #pragma omp parallel for default(none) private(j, intermediate) reduction(+:sum) shared(short_paths, inf_path, std::cout)
        for(i = 0; i < n; i++){
            intermediate = 0;
            for(j = 0; j < n; j++){
                if(INDDGO_INFINITY != short_paths[i][j]){
                    intermediate += short_paths[i][j];
                }
            }
            sum += (double)(intermediate / (double)((n * (n - 1)) - inf_path));
            //cout << "For node " << i << " got sum " << sum << endl;
            //sum = sum /  (double)((n * (n - 1)) - inf_path);
        }
        cout << "Got a SUM: " << sum << endl;
        pl = sum;
    } // avg_path_length

    /**
     * Calculates the degree distribution of graph g
     * \param[in] g Pointer to a graph
     * \param[out] dist Integer vector holding the degree distribution of the graph
     */
    void GraphProperties::deg_dist(Graph *g, vector<int> &dist){
        int i, size = g->degree.size();

        dist.resize(size);
        for(i = 0; i < size; i++){
            dist[i] = 0;
        }

        for(i = 0; i < size; i++){
            dist[g->degree[i]] += 1;
        }
    }

    #ifdef HAS_BOOST
    /**
     * Fits the degree distribution of g to a power law distribution
     * \param[in] g Pointer to a graph
     * \param[out] xmin Integer value holding the degree at which the data begins behaving like a power law
     * \param[out] alpha Double holding the most likely value for the exponent
     * \param[out] KS Double holding the Kolmogorov-Smirnov statistic
     * \param[in] start Lower bound for exponent, defaults to 1.5
     * \param[in] inc Determines granularity of the exponent search range, defaults to 0.01
     * \param[in] end Upper bound for exponent, defaults to 3.5
     */
    void GraphProperties::powerlaw(Graph *g, int &xmin, double &alpha, double &KS, double start, double inc, double end){
        int xm, i, k, n, I;
        int prev, pprev;
        int size = 1 + (int)(end - start) / inc;
        double slogz, L, MLE;
        double znorm, denom;
        double f, fn, D, tD;
        double exp;

        vector<int> x(g->degree);

        //Holds the possible \alpha values from start to end with a delta of inc
        vector< double > vec( size );
        //Holds the zeta of the \alpha values, since it is used a lot, no need to keep calculating
        vector< double > zvec( size );
        f = start;
        for(i = 0; i < size; i++ ){
            vec[i] = f;
            zvec[i] = boost::math::zeta<double>(f);
            f += inc;
        }

        sort(x.begin(), x.end());

        KS = -1;
        prev = 0;
        for( xm = 0; xm < x.size(); xm++ ){
            if(x[xm] == prev){
                continue;
            }

            n = x.size() - xm;

            slogz = 0;
            for(i = xm; i < x.size(); i++){
                slogz += log(x[i]);
            }

            //MLE loop: for each possible \alpha value in the specified grid with the current test xmin
            //          calculate the log-likelihood, normalizing the Hurwitz zeta to the Riemann zeta
            MLE = -1;
            I = 0;
            for( k = 0; k < vec.size(); k++ ){
                exp = vec[k];

                znorm = 0;
                for(i = 1; i < x[xm]; i++){
                    znorm += pow(i, -exp);
                }

                L = -exp * slogz - n * log( zvec[k] - znorm );

                if(MLE != -1){
                    if(L > MLE){
                        MLE = L;
                        I = k;
                    }
                }
                else {
                    MLE = L;
                    I = k;
                }
            }

            //compute KS statistic
            exp = vec[I];

            znorm = 0;
            for(i = 1; i < x[xm]; i++){
                znorm += pow(i, -exp);
            }
            denom = zvec[I] - znorm;

            pprev = 0;
            fn = f = 0;
            D = 0;
            for(i = xm; i < x.size(); i++){
                if(x[i] == pprev){
                    continue;
                }

                //CDF (f) and EDF (fn)
                for(k = i; k < x.size() && x[k] <= x[i]; k++){
                    fn += (1.0) / n;
                }
                f += pow(x[i], -exp) / denom;

                tD = abs(fn - f);
                if(tD > D){
                    D = tD;
                }

                pprev = x[i];
            }

            if(KS > -1){
                if(D < KS){
                    KS = D;
                    xmin = x[xm];
                    alpha = exp;
                }
            }
            else {
                KS = D;
                xmin = x[xm];
                alpha = exp;
            }

            prev = x[xm];
        }
    } // powerlaw

    #endif // ifdef HAS_BOOST

    /**
     * Calculates the degree assortativity of a graph g
     * \param[in] g input graph
     * \param[out] coeff the degree assortativity coefficient of g
     */
    void GraphProperties::deg_assortativity(Graph *g, double &coeff){
        double n1 = 0.0, n2 = 0.0, de = 0.0;
        int m = g->get_num_edges();
        const vector<int> &degrees = g->get_degree_ref();
        int i;
        double di, dc;
        std::list<int>::const_iterator cit;
        Node *node;

        #pragma omp parallel for schedule(dynamic, 8) default(none) private(i,di,dc,node,cit) reduction(+:n1,n2,de) shared(g,degrees)
        for(i = 0; i < g->get_num_nodes(); i++){
            node = g->get_node(i);
            const list<int> &nbrs = node->get_nbrs_ref();
            for(cit = nbrs.begin(); cit != nbrs.end(); ++cit){
                if(*cit > i){
                    di = (double)degrees[i] - 1;
                    dc = (double)degrees[*cit] - 1;

                    n1 += di * dc;
                    n2 += di + dc;
                    de += pow(di,2) + pow(dc,2);
                }
            }
        }

        n1 = n1 / (double)m;
        de = de / (2.0 * m);
        n2 = pow(n2 / (2.0 * m), 2);

        coeff = (n1 - n2) / (de - n2);
    } // deg_assortativity

    /**
     * Calculates the delta hyperbolicity of a graph g and returns the hyperbolicity distribution
     * \param[in] g Pointer to a graph
     * \param[out] mex_delta Double value holding the maximum delta value encountered
     * \param[out] delta Vector of vectors of doubles holding the delta hyperbolicity distribution
     */
    void GraphProperties::delta_hyperbolicity(Graph *g, double &max_delta, vector<vector<double> > &delta){
        const vector< vector<int> > &dist_mat = g->get_shortest_path_dist_ref();
        int size = dist_mat.size();
        int mat_size = dist_mat.size();
        int counter = 0;
        int row = 0;
        int diam;

        if(g->num_nodes < 4){
            cerr << "Graph is too small for hyperbolicity, returning zero\n";
            vector<double> delta_vec(1, 0);
            delta.push_back(delta_vec);
            return;
        }

        if(!is_connected(g)){
            cerr << "Graph passed to calc_gromov has multiple components\n";

            //Find the largest connected component
            int temp = 0;
            size = 0;
            for(int i = 0; i < mat_size; i++){
                temp = 0;
                for(int j = 0; j < mat_size; j++){
                    if(dist_mat[i][j] != INDDGO_INFINITY){
                        temp++;
                    }
                }
                //temp++; //count the vertex currently being tested
                if(temp > size){
                    size = temp;
                    row = i;
                }
            }
        }

        //build vertex set of largest connected component and find its diameter
        int temp_diam = 0;
        vector<int> vert_vec(size);
        for(int j = 0; j < mat_size; j++){
            if(dist_mat[row][j] != INDDGO_INFINITY){
                vert_vec[counter] = j;
                counter++;
                if(dist_mat[row][j] > temp_diam){
                    temp_diam = dist_mat[row][j];
                }
            }
            else if(j == row){
                vert_vec[counter] = j;
                counter++;
            }
        }
        diam = temp_diam;

        max_delta = 0;

        for(int i = 0; i < diam; i++){
            vector<double> zero_dist(2 * diam + 1,0);
            delta.push_back(zero_dist);
        }

        int num_threads = omp_get_max_threads();

        ////////////////////////////////////////////////////////////
        //These are the shared variables between threads.  Each thread has access to a certain pointer.
        //The pointer should then point to a copy of the variable for the thread's private (local) use.
        //Shared max_delta and shared_delta are updated by each thread, must be careful to avoid contention.
        ////////////////////////////////////////////////////////////

        int size_of_array = 64 / sizeof(double *);  //64 is cacheline size, needs to be updated on different machines for performance
        int type_size = sizeof(double *);
        if(type_size < sizeof(vector< vector<double> > *)){
            type_size = sizeof(vector< vector<double> > *);
        }

        counter = 1;
        while(num_threads * type_size > size_of_array){
            ++counter;
            size_of_array = counter * 64 / sizeof(double *);
        }

        ////////////////////////////////////////////////////////////
        //Pointers to shared data arrays that have copy of the delta information for each thread,
        //ie, if there are 8 threads, there will be an 8 (pointers to) arrays each accessed by only one thread,
        //At the end of the parallel section of code, the 8 arrays will be combined into a single array
        ///////////////////////////////////////////////////////////

        double **shared_max_delta;
        shared_max_delta = new double * [size_of_array];

        vector< vector<double> > **shared_delta;
        shared_delta = new vector< vector<double> > * [size_of_array];

        #if defined (TASK_PROFILE) || defined(TASK_PROFILE_OVERHEAD)
        int num_tasks1[1024][1024];
        int num_tasks2[1024][1024];
        int num_tasks1_switch[1024][1024];
        int num_tasks2_switch[1024][1024];
        int total_threads = 0;
        double task_time[1024][512];
        double time_starts = 0.0;
        double time_stops = 0.0;
        #endif

        //OMP parallel region, this current implementation generates tasks at two of the four levels of the nested for loop,
        //    these tasks are then passed to threads which execute them and update their private copy of the results array
        //After the tasks have all been generated and executed the results are collated
        #pragma omp parallel shared(delta, max_delta, shared_max_delta, diam, shared_delta, dist_mat, vert_vec, size)
        {
            double max_delta_loc = 0;

            vector<double> zero_dist(2 * diam + 1,0);
            zero_dist.reserve(256 / sizeof(double));

            vector< vector<double> > delta_loc(diam,zero_dist);

            int thread_id = omp_get_thread_num();

            #ifdef TASK_PROFILE
            num_tasks1[thread_id][0] = 0;
            num_tasks2[thread_id][0] = 0;
            num_tasks1_switch[thread_id][0] = 0;
            num_tasks2_switch[thread_id][0] = 0;
            int total_threads = 0;
            #endif
            #ifdef TASK_PROFILE_OVERHEAD
            task_time[thread_id][0] = 0.0;
            #endif

            shared_max_delta[thread_id] = &max_delta_loc;

            shared_delta[thread_id] = &delta_loc;

            //size_t size = vert_vec.size();

            #pragma omp single
            {
                #ifdef TASK_PROFILE_OVERHEAD
                time_starts = omp_get_wtime();
                #endif

                for(int i = 0; i < size; ++i){
                    #if defined (TASK_PROFILE) || defined (TASK_PROFILE_OVERHEAD)
                    #pragma omp task shared(shared_delta, shared_max_delta, dist_mat, vert_vec, size, num_tasks1, num_tasks2, task_time)
                    #else
                    #pragma omp task shared(shared_delta, shared_max_delta, dist_mat, vert_vec, size)
                    #endif
                    {
                        #ifdef TASK_PROFILE
                        int thread_id = omp_get_thread_num();
                        num_tasks1[thread_id][0]++;
                        #endif
                        #ifdef TASK_PROFILE_OVERHEAD
                        double tt1 = omp_get_wtime();
                        #endif

                        for(int j = i + 1; j < size; ++j){
                            #if defined (TASK_PROFILE) || defined (TASK_PROFILE_OVERHEAD)
                            #pragma omp task shared(shared_delta, shared_max_delta, dist_mat, vert_vec, size, num_tasks1, num_tasks2, task_time)
                            #else
                            #pragma omp task shared(shared_delta, shared_max_delta, dist_mat, vert_vec, size)
                            #endif
                            {
                                int thread_num = omp_get_thread_num();

                                #ifdef TASK_PROFILE
                                num_tasks2[thread_num][0]++;
                                #endif
                                #ifdef TASK_PROFILE_OVERHEAD
                                double ttt1 = omp_get_wtime();
                                #endif

                                vector<vector<double> > &loc_delta = *shared_delta[thread_num];

                                double &max_delta_loc = *shared_max_delta[thread_num];

                                for(int k = j + 1; k < size; ++k){
                                    for(int l = k + 1; l < size; ++l){
                                        const int ui = vert_vec[i];

                                        const int vi = vert_vec[j];
                                        const int uv = (dist_mat)[ui][vi];

                                        const int xi = vert_vec[k];
                                        const int ux = (dist_mat)[ui][xi];
                                        const int vx = (dist_mat)[vi][xi];

                                        const int yi = vert_vec[l];
                                        const int uy = (dist_mat)[ui][yi];
                                        const int vy = (dist_mat)[vi][yi];
                                        const int xy = (dist_mat)[xi][yi];

                                        const int s1 = uv + xy;
                                        const int s2 = ux + vy;
                                        const int s3 = uy + vx;

                                        double currDelta = 0;
                                        unsigned long maxPair = 1;

                                        //find max pair distance
                                        maxPair = uv;

                                        if(xy > maxPair){
                                            maxPair = xy;
                                        }

                                        if(ux > maxPair){
                                            maxPair = ux;
                                        }

                                        if(vy > maxPair){
                                            maxPair = vy;
                                        }

                                        if(uy > maxPair){
                                            maxPair = uy;
                                        }

                                        if(vx > maxPair){
                                            maxPair = vx;
                                        }

                                        //Delta calculations

                                        if((s1 >= s2) && (s3 >= s2) ){
                                            currDelta = abs(s1 - s3) / 2.0; //works since s1 and s2 are positive
                                        }
                                        else if((s2 >= s1) && (s3 >= s1) ){
                                            currDelta = abs(s2 - s3) / 2.0;
                                        }
                                        else {
                                            currDelta = abs(s1 - s2) / 2.0;
                                        }

                                        ++(loc_delta)[maxPair - 1][2 * currDelta + 1];

                                        if(currDelta > max_delta_loc){
                                            max_delta_loc = currDelta;
                                        }

                                        ++(loc_delta)[maxPair - 1][0];
                                    }
                                }

                                #ifdef TASK_PROFILE_OVERHEAD
                                double ttt2 = omp_get_wtime();
                                task_time[thread_num][0] += (ttt2 - ttt1);
                                #endif
                                #ifdef TASK_PROFILE
                                int thread_id_end = omp_get_thread_num();
                                if(thread_id != thread_id_end){
                                    num_tasks2_switch[thread_id_end][0]++;
                                }
                                #endif
                            } // end task
                        }

                        #ifdef TASK_PROFILE_OVERHEAD
                        //double tt2 = omp_get_wtime();
                        //task_time[thread_id][0]+= (tt2 - tt1);
                        #endif
                        #ifdef TASK_PROFILE
                        int thread_id_end = omp_get_thread_num();
                        if(thread_id != thread_id_end){
                            num_tasks1_switch[thread_id_end][0]++;
                        }
                        #endif
                    } // end task
                }
            } // end single

            #ifdef TASK_PROFILE_OVERHEAD
            #pragma omp single nowait
            {
                time_stops = omp_get_wtime();
            }
            #endif

            double t1 = time(NULL);

            #pragma omp critical (collate_delta)
            {
                int thread_num = omp_get_thread_num();
                vector<vector<double> > *delta_point = shared_delta[thread_num];
                vector<vector<double> > &loc_delta = *delta_point;

                int i_size = loc_delta.size();
                for(int i = 0; i < i_size; ++i){
                    int j_size = loc_delta[i].size();
                    for(int j = 0; j < j_size; ++j){
                        delta[i][j] += loc_delta[i][j];
                    }
                }

                if(max_delta_loc > max_delta){
                    max_delta = max_delta_loc;
                }
            }

            double t2 = time(NULL);

            #pragma omp single
            {
                cout << "Critical Region time: " << t2 - t1 << "\n";
            }
        } //END OF FIRST PRAGMA

        //take care of pointer!
        delete [] shared_delta;
        delete [] shared_max_delta;

        //Needs removed before parallel code is implemented
        for(int i = 0; i < diam; ++i){
            delta[i].erase(delta[i].begin() + 2 * max_delta + 2, delta[i].end());
        }

        #if defined (TASK_PROFILE) || defined (TASK_PROFILE_OVERHEAD)

        #pragma omp parallel shared(total_threads)
        {
            #pragma omp master
            {
                total_threads = omp_get_num_threads();
            }
        }

        cout << endl << endl << "******Tasks Statistics: (num_task, num_tasks, migrated_tasks, switched_task, task_time)" << endl;

        for(int i = 0; i < total_threads; i++){
            cout << i
        #endif
        #ifdef TASK_PROFILE
            << "," << num_tasks1[i][0] << "," << num_tasks2[i][0] << "," << num_tasks1_switch[i][0] << "," << num_tasks2_switch[i][0]
            #endif
            #ifdef TASK_PROFILE_OVERHEAD
            << "," << task_time[i][0]
            #endif
            #if defined (TASK_PROFILE) || defined (TASK_PROFILE_OVERHEAD)
            << endl;
    } // delta_hyperbolicity

    cout << endl;

        #endif
        #ifdef TASK_PROFILE_OVERHEAD
        cout << "Time in task region = " << time_stops - time_starts << endl;
        #endif
    } // delta_hyperbolicity

    #ifdef HAS_PETSC
    /**
     * Built in self tester for the eigen solver. Looks at Ax - yx and looks for the
     * largest magnitude of error.
     */
    void check_eigen_values(Mat A, Vec eigen_vec, PetscScalar eigen_val){
        Vec xr, xi, xr2, result1, result2;
        MatGetVecs(A, PETSC_NULL, &xr2);
        MatGetVecs(A, PETSC_NULL, &result1);
        MatGetVecs(A, PETSC_NULL, &result2);

        cout << "Eigen vector is:\n";
        VecView(eigen_vec,PETSC_VIEWER_STDOUT_SELF);
        cout << "Eigen value is: " << eigen_val << endl;

        MatMult(A,eigen_vec,result1);
        VecSet(xr2,eigen_val);
        VecPointwiseMult(result2,xr2,eigen_vec);
        cout << "Ax = " << endl;
        VecView(result1,PETSC_VIEWER_STDOUT_SELF);
        cout << "yx = " << endl;
        VecView(result1,PETSC_VIEWER_STDOUT_SELF);
        PetscScalar *a;
        PetscScalar *b;
        VecGetArray(result1, &a);
        VecGetArray(result2, &b);
        PetscInt size;
        VecGetLocalSize(result1, &size);
        double max = 0.0,current;
        for(int idx = 0; idx < size; idx++){
            current = fabs(a[idx] - b[idx]);
            max = current > max ? current : max;
        }
        cout << "Magnitude of greatest error is: " << max << endl;
        VecRestoreArray(result1, &a);
        VecRestoreArray(result2, &b);
    } // check_eigen_values

    /**
     * looks for eigen values in the adjacency matrix
     * \param[in] g input graph
     * \param[out] eigen_values vector of eigen values found
     * \param[in] spread the number of values wanted, spread high values and spred low values, so eigen_values can be upto 2*spread in size
     */
    void GraphProperties::eigen_spectrum(Graph *g, vector<double> &eigen_values, int spread){
        #ifndef HAS_SLEPC
        fatal_error("Called SLEPC eigen solvers without HAS_SLEPC.\n");
        #else
        GraphUtil graph_util;
        graph_util.populate_PetscMat(g);
        EPS eps;
        PetscInt nconv;
        PetscScalar kr,ki;
        EPSCreate(PETSC_COMM_WORLD, &eps);
        EPSSetType(eps, EPSPOWER);
        EPSSetOperators(eps,g->PetscMat,PETSC_NULL);
        EPSSetLeftVectorsWanted(eps,PETSC_TRUE);
        EPSSetFromOptions(eps);
        EPSSetDimensions(eps,spread * 2,spread * 8,spread * 8);
        EPSSolve(eps);
        EPSGetConverged(eps,&nconv);
        EPSGetConverged(eps,&nconv);
        eigen_values.resize(nconv);
        for(int idx = 0; idx < nconv; idx++){
            EPSGetEigenvalue(eps,idx,&kr,&ki);
            eigen_values[idx] = kr;
            #ifdef EIGENSOLVER_SELFTEST
            //built in self tester. Don't use in production runs.
            Vec xr, xi;
            MatGetVecs(g->PetscMat, PETSC_NULL, &xr);
            MatGetVecs(g->PetscMat, PETSC_NULL, &xi);
            EPSGetEigenpair(eps,idx,&kr,&ki, xr, xi);
            check_eigen_values(g->PetscMat, xr, kr);
            #endif
        }
        EPSDestroy(&eps);

        #endif // ifndef HAS_SLEPC
    } // eigen_spectrum

    #endif // ifdef HAS_PETSC
}
