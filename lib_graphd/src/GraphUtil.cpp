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

#include "GraphUtil.h"
#include "GraphDecomposition.h"

#ifdef HAS_BOOST
  #include <iostream>
  #include <deque>
  #include <iterator>

  #include "boost/graph/adjacency_list.hpp"
  #include "boost/graph/topological_sort.hpp"
#endif

namespace Graph {
    GraphUtil::GraphUtil(){
    }

    GraphUtil::~GraphUtil(){
    }

    /**
     * Recomputes the entries in the degree[] array using the size of
     * the adjacency lists. Updates the number of edges in the graph.
     */
    void GraphUtil::recompute_degrees(Graph *g){
        g->num_edges = 0;
        int i;
        //set the degrees to 0 - should g be -1 if the nodes label is -1?
        // CSG - all these loops used to start at i=1?
        fill(g->degree.begin(), g->degree.end(), 0);

        //not strictly necessary to symmetrize, but much easier - lists are sorted now.
        if(g->simple == true){
            for(i = 0; i < g->capacity; i++){
                g->degree[i] = g->nodes[i].nbrs.size();
                g->num_edges += g->degree[i];
            }
            // We double counted all the edges;
            g->num_edges = g->num_edges / 2;
        }
        else {
            print_message(0,"Not simple??\n");
            // we could have loops or multiedges - the former are the only thing that screw up our accounting.
            size_t origsize, newsize;
            int loops = 0;
            int j;
            for(i = 0; i < g->capacity; i++){
                origsize = g->nodes[i].nbrs.size();
                g->nodes[i].nbrs.remove(i);
                newsize = g->nodes[i].nbrs.size();
                if(g->nodes[i].nbrs.size() < origsize){
                    loops += (int) (origsize - newsize);
                    g->num_edges += (int) (newsize);                         //all but the loop are normal edges
                    for(j = 0; j < (int) (origsize - newsize); j++){
                        g->nodes[i].nbrs.push_back(i);                         //add it back
                    }
                    g->degree[i] = g->nodes[i].nbrs.size();
                }
                else {
                    g->degree[i] = origsize;
                    g->num_edges += g->degree[i];                     //all are normal edges
                }
            }

            // We double counted all the normal edges & need to add the loops
            g->num_edges = g->num_edges / 2 + loops;
        }

        return;
    } // recompute_degrees

    /**
     * Returns the index of a vertex with the lowest degree. If there are multiple minimum degree vertices,
     * then a random representative is selected.
     */
    int GraphUtil::get_random_high_degree_vertex(Graph *g) const {
        vector<int> high_degree_vs(g->capacity, -1);
        int max_deg = -INT_MAX;
        int i, j;
        j = 0;
        for(i = 0; i < g->capacity; i++){
            if((g->nodes[i].label != -1)
               && ( (int) (g->nodes[i].nbrs.size())
                    >= max_deg) ){
                max_deg = g->nodes[i].nbrs.size();
                high_degree_vs[j] = i;
                //print_message(0, "High[%d]=%d(deg=%d)\n", j, i, max_deg);
                j++;
            }
        }

        // We now have j high degree vertices in positions 0,1,...,j-1
        i = j - 1;
        while(high_degree_vs[i] == high_degree_vs[i - 1]){
            i--;
        }

        // Pick a random integer k in {i,i+1,...,j-1} and return
        // high_degree_vs[k]
        int k = rand_int(i, j - 1);
        k = high_degree_vs[k];
        // sanity check
        if(k == -1){
            fatal_error("%s: Didn't find valid low degree vertex!\n", __FUNCTION__);
        }

        return k;
    } // get_random_high_degree_vertex

    /**
     * Returns the index of a vertex with the lowest degree. If there are multiple minimum degree vertices,
     * then a random representative is selected.
     */
    int GraphUtil::get_random_low_degree_vertex(Graph *g) const {
        vector<int> low_degree_vs(g->capacity, -1);
        int min_deg = INT_MAX;
        int i, j;
        j = 0;
        for(i = 0; i < g->capacity; i++){
            if((g->nodes[i].label != -1) && (( (int) (g->nodes[i].nbrs.size())) <= min_deg) ){
                min_deg = g->nodes[i].nbrs.size();
                low_degree_vs[j] = i;
                print_message(1, "Low[%d]=%d(deg=%d)\n", j, i, min_deg);
                j++;
            }
        }

        //BDS - fixed Aug 5.
        if(j == 1){
            return low_degree_vs[0];
        }
        // We now have j low degree vertices in positions 0,1,...,j-1
        i = j - 1;
        while(low_degree_vs[i] == low_degree_vs[i - 1]){
            i--;
        }

        // Pick a random integer k in {i,i+1,...,j-1} and return
        // low_degree_vs[k]
        int k = rand_int(i, j - 1);
        k = low_degree_vs[k];
        // sanity check
        if(k == -1){
            fatal_error("%s: Didn't find valid low degree vertex!\n", __FUNCTION__);
        }

        return k;
    } // get_random_low_degree_vertex

    /**
     * Given the components vector of size capacity, g uses a non-recursive
     * procedure to find the positions in the nodes[] array of those nodes that
     * belong to the same component C as v and sets component[i]=label for
     * all i in C.  The value of label should be non-negative and the components
     * vector should be initially set with all entries -1.  Returns the
     * number of nodes in g component.
     */
    int GraphUtil::label_component(Graph *g, int v, int label,
                                   vector<int> *components){
        // g assumes that components vector is initially set to all -1's (-1's)
        if(g->nodes[v].label == -1){
            fatal_error("%s:  Tried to find component for -1 position??\n",
                        __FUNCTION__);
        }

        if(label < 0){
            fatal_error("%s:  The value of label should be non-negative\n",
                        __FUNCTION__);
        }

        print_message(1, "Labeling component of vertex %d with %d\n", v, label);
        if((int) (components->size()) < g->capacity){
            fatal_error("%s:  Component vector does not have enough space.\n",
                        __FUNCTION__);
        }

        //components->resize(g->capacity,-1);
        // Just return 0 if we have already labeled v in the components array
        if(components->at(v) != -1){
            return 0;
        }

        int j, cnt = 0;
        list<int> S;
        list<int>::iterator ii;
        // Put v on the stack
        S.push_front(v);
        while(!(S.empty())){
            j = S.front();
            S.pop_front();
            if(components->at(j) != label){
                components->at(j) = label;
                cnt++;
                for(ii = g->nodes[j].nbrs.begin(); ii != g->nodes[j].nbrs.end();
                    ++ii){
                    if(components->at(*ii) != label){
                        print_message(100, "Pushing %d\n", *ii);
                        S.push_front(*ii);
                    }
                }
            }
        }

        return cnt;
    } // label_component

    /**
     * Given the components vector of size capacity, g uses a recursive
     * procedure to find the positions in the nodes[] array of those nodes that
     * belong to the same component C as v and sets component[i]=label for
     * all i in C.  The value of label should be non-negative and the components
     * vector should be initially set with all entries -1.  Returns the
     * number of nodes in g component.
     */
    int GraphUtil::rec_label_component(Graph *g, int v, int label,
                                       vector<int> *components){
        // g assumes that components vector is initially set to all -1's (-1's)
        if(g->nodes[v].label == -1){
            fatal_error("%s:  Tried to find component for -1 position??\n",
                        __FUNCTION__);
        }

        if(label < 0){
            fatal_error("%s:  The value of label should be non-negative\n",
                        __FUNCTION__);
        }

        print_message(1, "Labeling component of vertex %d with %d\n", v, label);
        if((int) (components->size()) < g->capacity){
            fatal_error("%s:  Component vector does not have enough space.\n",
                        __FUNCTION__);
        }

        //components->resize(g->capacity,-1);
        // Just return 0 if we have already labeled v in the components array
        if(components->at(v) != -1){
            return 0;
        }

        // Mark v
        components->at(v) = label;
        list<int>::iterator ii;
        // Recursive DFS
        for(ii = g->nodes[v].nbrs.begin(); ii != g->nodes[v].nbrs.end(); ++ii){
            print_message(1, "Neighbor %d of %d\n", *ii, v);
            // Put if check here to eliminate unnecessary function calls?
            if(components->at(*ii) == -1){
                rec_label_component(g, *ii, label, components);
            }
        }
        // Count the # of nodes in g component
        int cnt = 0;
        for(int i = 0; i < g->capacity; i++){
            if(components->at(i) == label){
                cnt++;
            }
        }

        return cnt;
    } // rec_label_component

    /**
     * Non-recursive function that labels all connected components in the graph
     * by filling the components vector with the component label for each node.
     * Every entry in the components vector is initially set to
     * -1 and is eventually populated with integers in the range
     * 0,1,...,num_components-1.  If components[i]=k, then g means
     * that the node in position i in the nodes[] array is in component k.
     * Returns the number of components found.
     */
    int GraphUtil::label_all_components(Graph *g,
                                        vector<int> *components){
        // Make sure the components vector is the right size
        if((int) (components->size()) < g->capacity){
            components->resize(g->capacity, -1);
        }

        else {
            // Initialize components vector to all -1's
            fill(components->begin(), components->end(), -1);
        }

        int i = 0, j = 0;
        // Find the first "active" node
        while(g->nodes[j].label == -1){
            j++;
        }

        list<int> S;
        list<int>::iterator ii;
        // Use S as a stack - use only push_front, pop_front
        bool *visited = new bool[g->capacity];
        int *origin = new int[g->capacity];
        int *component_labels = new int[g->capacity];
        for(i = 0; i < g->capacity; i++){
            visited[i] = false;
            origin[i] = component_labels[i] = -1;
        }
        int numc = 0;
        for(i = 0; i < g->capacity; i++){
            if(!(visited[i]) && (g->nodes[i].label != -1) ){
                S.push_front(i);
                while(!(S.empty())){
                    j = S.front();
                    S.pop_front();
                    if(!visited[j]){
                        visited[j] = true;
                        origin[j] = i;
                        if(origin[j] == j){
                            // g indicates a new component
                            // Any node u in g component will have origin[u]=j
                            // Keep track of these "special" values by marking position
                            // j in the component_labels[] array and incrementing numc
                            component_labels[j] = numc;
                            numc++;
                        }
                        // Check out j's neighbors - assume these are not -1 nodes!
                        for(ii = g->nodes[j].nbrs.begin();
                            ii != g->nodes[j].nbrs.end(); ++ii){
                            if(!visited[*ii]){
                                S.push_front(*ii);
                            }
                        }
                    }
                }
            }

            if(g->nodes[i].label != -1){
                // i is visited once we get here, so we know the value of components[i]
                components->at(i) = component_labels[origin[i]];
            }
        }
        // Clean up
        delete[] visited;
        delete[] origin;
        delete[] component_labels;
        // Set the value of and return the # of connected components
        g->num_connected_components = numc;
        return g->num_connected_components;
    } // label_all_components

    /**
     * Recursive function that labels all connected components in the graph
     * by filling the components vector with the component label for each node.
     * Every entry in the components vector is initially set to
     * -1 and is eventually populated with integers in the range
     * 0,1,...,num_components-1.  If components[i]=k, then g means
     * that the node in position i in the nodes[] array is in component k.
     * Returns the number of components found.
     */
    int GraphUtil::rec_label_all_components(Graph *g,
                                            vector<int> *components){
        // Make sure the components vector is the right size
        if((int) (components->size()) < g->capacity){
            components->resize(g->capacity, -1);
        }

        else {
            // Initialize components vector to all -1's
            fill(components->begin(), components->end(), -1);
        }

        int i = 0, label = 0, j = 0;
        // Find the first "active" node
        while(g->nodes[j].label == -1){
            j++;
        }

        // Label the component containing node i
        rec_label_component(g, j, label, components);
        for(i = j + 1; i < g->capacity; i++){
            if((components->at(i) == -1) && (g->nodes[i].label != -1) ){
                // We found an active node that we haven't seen before
                label++;
                rec_label_component(g, i, label, components);
            }
        }

        // We now know number of components
        g->num_connected_components = label + 1;
        // components[] will now contain integers in 0,1,...,label-1
        // where components[v]=c means that node v is in component c
        // If components[v]=-1, then g node must have an -1 label
        // and is therefore not in a component of the graph
        // Return the # of components we discovered
        return g->num_connected_components;
    } // rec_label_all_components

    /**
     * Uses a recursive function to populate the members list with the
     * positions in the nodes[] array of those
     * nodes that belong to the saem component as v.
     * Returns the number of nodes in g component.
     */
    int GraphUtil::rec_find_component(Graph *g, int v,
                                      list<int> *members){
        // Create a vector of -1's and then label v's component w/ 0
        vector<int> components(g->capacity, -1);
        rec_label_component(g, v, 0, &components);
        members->clear();
        for(int i = 0; i < g->capacity; i++){
            if(components.at(i) == 0){
                members->push_back(i);
            }
        }
        // Sort and return the component size
        members->sort();
        return members->size();
    }

    /**
     * Uses a non-recursive function to populate the members list with the
     * positions in the nodes[] array of those
     * nodes that belong to the saem component as v.
     * Returns the number of nodes in g component.
     * Oct 18 2010 - was not populating members list, so corrected.
     */
    int GraphUtil::find_component(Graph *g, int v, list<int> *members){
        // g assumes that components vector is initially set to all -1's (-1's)
        if(g->nodes[v].label == -1){
            fatal_error("%s:  Tried to find component for -1 position??\n",
                        __FUNCTION__);
        }

        members->clear();
        int j, cnt = 0;
        vector<bool> visited(g->capacity, false);
        list<int> S;
        list<int>::iterator ii;
        // Put v on the stack
        S.push_front(v);
        while(!(S.empty())){
            j = S.front();
            S.pop_front();
            if(visited[j] == false){
                visited[j] = true;
                cnt++;
                members->push_back(j);
                for(ii = g->nodes[j].nbrs.begin(); ii != g->nodes[j].nbrs.end();
                    ++ii){
                    if(visited[*ii] == false){
                        print_message(100, "Pushing %d\n", *ii);
                        S.push_front(*ii);
                    }
                }
            }
        }

        return cnt;
    } // find_component

    /**
     * Populate the xadj and adjncy vectors with the graph data.
     */
    void GraphUtil::populate_CRS(Graph *g){
        // Only works if capacity=num_nodes
        if(g->capacity != g->num_nodes){
            fatal_error("%s:  Requires capacity=num_nodes!\n", __FUNCTION__);
        }

        g->xadj.resize(g->num_nodes + 1);
        g->adjncy.resize(2 * g->num_edges);
        print_message(10, "%s:  %d nodes, %d edges\n", __FUNCTION__, g->num_nodes,
                      g->num_edges);
        // Nodes are numbered 0,1,...,n-1
        // The neighbors of node i are adjncy[xadj[i]],adjncy[xadj[i]+1],...adjncy[xadj[i+1]-1]
        int i, j = 0;
        for(i = 0; i < g->num_nodes; i++){
            g->xadj[i] = j;
            for(list<int>::iterator ii = g->nodes[i].nbrs.begin();
                ii != g->nodes[i].nbrs.end(); ++ii){
                g->adjncy[j] = *ii;
                j++;
            }
        }

        g->xadj[g->num_nodes] = j;
        print_message(1, "Set xadj[%d]=%d (nedges=%d)\n", g->num_nodes, j,
                      g->num_edges);
        return;
    } // populate_CRS

    /**
     * Free the memory used to create the CRS format.
     */
    void GraphUtil::free_CRS(Graph *g){
        g->xadj.clear();
        g->adjncy.clear();
    }

    #ifdef HAS_BOOST
    /**
     * Populates the boost_graph member of g
     * \param[in] g the input graph
     */
    void GraphUtil::populate_boost(Graph *g){
        const int n = g->get_num_nodes();
        int v;
        boost::graph_traits < BoostUndirected >::vertex_descriptor a, b;
        list<int>::const_iterator it;
        g->boost_graph = new BoostUndirected(n);
        BoostUndirected *bg = g->boost_graph;
        for(v = 0; v < n; v++){
            const list<int> &nbrs = g->get_node(v)->get_nbrs_ref();
            for(it = nbrs.begin(); it != nbrs.end(); ++it){
                if(v < *it){
                    a = boost::vertex(v, *bg);
                    b = boost::vertex(*it, *bg);
                    boost::add_edge(a, b, 1, *bg);
                }
            }
        }
    } // populate_boost

    #endif //HAS_BOOST

    #ifdef HAS_PETSC
    /**
     * Populates the PetscMat member of graph g. This function will also call populate_CRS() if it hasn't been already.
     * Durring cleanup a free_CRS() also needs to be done.
     * \param[in] g the input graph
     */
    void GraphUtil::populate_PetscMat(Graph *g){
        if(g->adjncy.size() == 0){
            populate_CRS(g);
        }
        int start, end;

        //Matrix is g->PetscMat
        PetscErrorCode ierror;
        ierror = MatCreate(PETSC_COMM_WORLD,&(g->PetscMat));
        CHKERRABORT(PETSC_COMM_WORLD,ierror);
        ierror = MatSetSizes(g->PetscMat,PETSC_DECIDE,PETSC_DECIDE,g->num_nodes,g->num_nodes);
        CHKERRABORT(PETSC_COMM_WORLD,ierror);
        ierror = MatSetFromOptions(g->PetscMat);
        CHKERRABORT(PETSC_COMM_WORLD,ierror);
        ierror = MatSetUp(g->PetscMat);
        CHKERRABORT(PETSC_COMM_WORLD,ierror);
        int row, col, size;
        const PetscScalar one_val = 1;
        for(int node = 0; node < g->num_nodes; node++){
            start = g->xadj[node];
            end = g->xadj[node + 1];
            size = 1 + end - start;
            row = node;
            PetscScalar ones[size];
            PetscInt columns[size];
            for(int idx = 0; idx < size; idx++){
                ones[idx] = 1.0;
                columns[idx] = g->adjncy[start + idx];
            }
            ierror = MatSetValues(g->PetscMat, 1, &row, size, columns, ones, INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD,ierror);
        }

        ierror = MatAssemblyBegin(g->PetscMat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_WORLD,ierror);
        ierror = MatAssemblyEnd(g->PetscMat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_WORLD,ierror);
    } // populate_PetscMat

    /**
     * Frees the PetscMat member of g. This function will not cleanup CRS, even if populate_PetscMat() filled it.
     * You will have to also call free_CRS();
     * \param[in] g the input graph
     */
    void GraphUtil::free_PetscMat(Graph *g){
        PetscErrorCode ierror = MatDestroy(&(g->PetscMat));
        CHKERRABORT(PETSC_COMM_WORLD,ierror);
    }

    #endif // ifdef HAS_PETSC

    /**
     * Non-recursive function that fills the members vector with the
     * lists of nodes belonging to the components.  The members vector
     * is resized so that upon completion it has num_connected_components
     * entries and members[i] is a pointer to a list of the nodes in the
     * i-th component.
     */
    int GraphUtil::find_all_components(Graph *g,
                                       vector<list<int> *> *members){
        int i = 0, j = 0;
        // Find the first "active" node
        while(g->nodes[j].label == GD_UNDEFINED){
            j++;
        }

        //list<int> *L;
        list<int> S;
        list<int>::iterator ii;
        // Use S as a stack - use only push_front, pop_front
        bool *visited = new bool[g->capacity];
        int *origin = new int[g->capacity];
        int *component_labels = new int[g->capacity];
        for(i = 0; i < g->capacity; i++){
            visited[i] = false;
            origin[i] = component_labels[i] = -1;
        }
        int numc = 0;
        for(i = 0; i < g->capacity; i++){
            if(!(visited[i]) && (g->nodes[i].label != GD_UNDEFINED) ){
                S.push_front(i);
                while(!(S.empty())){
                    j = S.front();
                    S.pop_front();
                    if(!visited[j]){
                        visited[j] = true;
                        origin[j] = i;
                        if(origin[j] == j){
                            // g indicates a new component
                            // Any node u in g component will have origin[u]=j
                            // Keep track of these "special" values by marking position
                            // j in the component_labels[] array and incrementing numc
                            component_labels[j] = numc;

                            //L = new list<int>;
                            //members->resize(numc + 1);
                            //members->at(numc) = L;
                            members->push_back(new list<int>);
                            numc++;
                        }
                        // Check out j's neighbors - assume these are not GD_UNDEFINED nodes!
                        for(ii = g->nodes[j].nbrs.begin();
                            ii != g->nodes[j].nbrs.end(); ++ii){
                            if(!visited[*ii]){
                                S.push_front(*ii);
                            }
                        }
                    }
                }
            }
            if(g->nodes[i].label != GD_UNDEFINED){
                // i is visited once we get here, so we know that i belongs in list
                // component_labels[origin[i]]
                members->at(component_labels[origin[i]])->push_back(i);                 // Obviously...
            }
        }
        // Clean up
        delete[] visited;
        delete[] origin;
        delete[] component_labels;
        // Sort the lists for convenience
        vector<list<int> *>::iterator jj;
        for(jj = members->begin(); jj != members->end(); ++jj){
            (*jj)->sort();
        }

        // Set the value of and return the # of connected components
        g->num_connected_components = numc;
        return g->num_connected_components;
    } // find_all_components

    int GraphUtil::vertex_separator(Graph *g, list<int> *V,
                                    vector<list<int> *> *members){
        //create a copy of G
        //Changed from = *g to copy constructor July 19 - BDS
        Graph H(*g);
        //remove the vertices in V
        list<int>::iterator i;
        for(i = V->begin(); i != V->end(); ++i){
            H.remove_vertex(*i);
        }
        //find the components of H
        int components = find_all_components(&H, members);
        return components;
    }

    //runs a BFS from start using
    //only vertices with allowed[v] = true. If a vertex u  is reachable,
    //Returns an array of booleans which is true for vertices which are found in the search.
    //num_reached is the number of vertices which are reachable. g does not include
    //disallowed vertices, and does include the start.
    bool *GraphUtil::bfs(Graph *g, int start, bool *allowed,
                         int *num_reached){
        int *dists = bfs_dist(g, start, allowed, num_reached);
        int j;
        bool *visited = new bool[g->capacity];
        for(j = 0; j < g->capacity; j++){
            if(dists[j] == GD_INFINITY
               ){
                visited[j] = false;
            }
            else {
                visited[j] = true;
            }
        }
        delete[] dists;
        return visited;
    }

    /**
     * Finds the degree 0 nodes in the graph and fills the isolated_nodes
     * list with their positions in the nodes[] array.  Returns the
     * number of isolated nodes.
     */
    int GraphUtil::find_isolated_nodes(Graph *g, list<int> *isolated_nodes){
        int i, k = 0;

        for(i = 0; i < g->capacity; i++){
            if((g->nodes[i].label != -1) && (g->nodes[i].nbrs.size() == 0) ){
                isolated_nodes->push_back(i);
                k++;
            }
        }
        return k;
    }

    //runs a BFS from start using only vertices with allowed[v] = true.
    //Returns an array of integers giving distance from the source (0 for source).
    //Unreachable vertices have value GD_INFINITY.
    //num_reached is the number of vertices which are reachable. g does not include
    //disallowed vertices, and does include the start.
    //Returns the maximum distance (vertex eccentricity of start) in ecc.

    int *GraphUtil::bfs_dist(Graph *g, int start, bool *allowed,
                             int *num_reached, int *ecc){
        // We need the graph to be symmetric, as we're going to walk through a neighbor list
        if(!g->canonical){
            fatal_error("%s:  must be in canonical format\n", __FUNCTION__);
        }
        int num_found = 0;

        int j;
        int *dists = new int[g->capacity];
        for(j = 0; j < g->capacity; j++){
            dists[j] = GD_INFINITY;
        }

        // S is a stack that will contain the nodes we have visited -
        // take the first one off, and then add its unvisited neighbors to the end
        list<int> S;
        list<int>::iterator ii;

        int curr_level = 0;
        int left_in_level = 1;
        int next_level = 0;

        // Put start on the stack
        S.clear();
        S.push_front(start);
        while(!(S.empty())){
            if(left_in_level == 0){
                curr_level++;
                left_in_level = next_level;
                next_level = 0;
            }

            // Remove the oldest element from the stack
            j = S.front();
            S.pop_front();
            left_in_level--;
            if(dists[j] == GD_INFINITY){
                // We have now visited j
                dists[j] = curr_level;
                num_found++;

                // Check j's neighbors
                for(ii = g->nodes[j].nbrs.begin(); ii != g->nodes[j].nbrs.end();
                    ++ii){
                    if((dists[*ii] == GD_INFINITY) && (allowed[*ii] == true) ){
                        // Note - g is the only place that we refer to the allowed[] vector
                        // We haven't seen *ii before and it is an "acceptable" vertex,
                        // so it is a candidate to be in the path - add it to the Stack
                        S.push_back(*ii);                         // must be push_back since we pop_front and need FIFO.
                        next_level++;
                    }
                }
            }
        }

        *ecc = curr_level - 1;
        *num_reached = num_found;
        // We emptied the stack.
        return dists;
    } // bfs_dist

    int *GraphUtil::bfs_dist(Graph *g, int start, bool *allowed,
                             int *num_reached){
        int ecc;
        return bfs_dist(g, start, allowed, num_reached, &ecc);
    }

    //Find the eccentricity of each vertex and store it in ecc, which is resized appropriately within this function.
    void GraphUtil::find_ecc(Graph *g, vector<int> *ecc){
        ecc->resize(g->capacity);

        //all nodes are allowed
        // CSG fixing this
        //bool allowed[g->capacity];
        bool *allowed;
        allowed = new bool[g->capacity];
        for(int i = 0; i < g->capacity; i++){
            allowed[i] = true;
        }
        int num_reached;
        int e;
        int *dists;
        for(int i = 0; i < g->capacity; i++){
            //for valid nodes
            if(g->nodes[i].label != -1){
                dists = this->bfs_dist(g, (g->nodes[i]).get_label() - 1, allowed, &num_reached, &e);
                (*ecc)[i] = e;
            }
        }
        delete[] dists;
        delete [] allowed;
        return;
    } // find_ecc

    //Calculate the maximum distance between nodes within a subset of vertices
    //given as a list
    int GraphUtil::subset_max_dist(Graph *g,  vector<int> subset){
        int max = 0;

        //all nodes are allowed
        //bool allowed[g->capacity];

        bool *allowed;
        allowed = new bool[g->capacity];
        for(int i = 0; i < g->capacity; i++){
            allowed[i] = true;
        }
        int num_reached;
        vector<int>::iterator it, it2;
        int *dists;
        for(it = subset.begin(); it != subset.end(); ++it){
            dists = this->bfs_dist(g, *it, allowed, &num_reached);
            for(it2 = subset.begin(); it2 != subset.end(); ++it2){
                if(max < dists[*it2]){
                    max = dists[*it2];
                }
            }
            delete[] dists;
        }
        delete [] allowed;
        return max;
    } // subset_max_dist

    /**
     * Calculates the k-core number of every vertex and stores in kcore vector
     * Return value is degeneracy (maximum non-empty k-core).
     * Uses algorithm of Batagelj and Zaversnik (2003)
     * Implemented by Timothy Goodrich and Matthew Farrell (2013)
     */
    int GraphUtil::find_kcore(Graph *g, vector<int> *kcore){
        int n = g->num_nodes;
        int k = 0;
        int last_deg = 0;
        kcore->resize(n);

        vector<int> deg_lookup(n);
        int max_deg = 0;
        int deg_size = deg_lookup.size();
        for(int i = 0; i < deg_size; i++){
            deg_lookup[i] = g->degree[i];
            if(deg_lookup[i] > max_deg){
                max_deg = deg_lookup[i];
            }
        }
        vector<vector<int> > D;
        D.resize(max_deg + 1);

        int *depth;
        depth = new int[n];
        //int depth[n];
        //can also create an L output list for coloring number optimal ordering

        for(int i = 0; i < n; i++){
            deg_lookup[i] = g->degree[i];
            D[deg_lookup[i]].push_back(i);
            depth[i] = D[deg_lookup[i]].size() - 1;
        }

        for(int i = 0; i < n; i++){
            int v;
            for(int j = last_deg; j <= max_deg; j++){
                if(D[j].size() != 0){
                    v = D[j].back();
                    D[j].pop_back();
                    break;
                }
            }

            (*kcore)[v] = k = max(k,deg_lookup[v]);
            last_deg = max(0, deg_lookup[v] - 1);
            deg_lookup[v] = -1;

            //shift each of v's neighbors down one degree
            list<int> *nbrs = g->nodes[v].get_nbrs_ptr();
            list<int>::const_iterator it;
            for(it = nbrs->begin(); it != nbrs->end(); it++){
                if(deg_lookup[*it] != -1){
                    int it_deg = deg_lookup[*it]--;
                    D[it_deg][depth[*it]] = D[it_deg][D[it_deg].size() - 1];
                    depth[D[it_deg][depth[*it]]] = depth[*it];
                    D[it_deg].pop_back();
                    D[it_deg - 1].push_back(*it);
                    depth[*it] = D[it_deg - 1].size() - 1;
                }
            }
        }
        delete [] depth;
        return k;
    } // find_degen

    /**
     * \input g input graph
     * \output new graph consisting of largest component, or NULL if the input graph is connected
     * NOTE: this function *will* canonicalize your graph!!!
     **/
    Graph *GraphUtil::get_largest_component_graph(Graph *g){
        GraphProperties properties;
        GraphCreatorFile creator;
        Graph *largecomp;
        list<int> *maxlist = NULL;
        size_t maxlen = 0;

        if(!properties.is_connected(g)){
            properties.make_canonical(g);
            vector<list<int> *> members;
            vector<list<int> *>::iterator vlit;
            find_all_components(g, &members);
            // traverse the list and find the one with the most members
            for(vlit = members.begin(); vlit != members.end(); ++vlit){
                //cout << "Got list with size: " << ((*vlit)->size()) << endl;
                if(((*vlit)->size()) >= maxlen){
                    maxlen = (*vlit)->size();
                    //cout << "New list with maxlen > " << maxlen << endl;
                    maxlist = *vlit;
                    //cout << maxlist << endl;
                }
            }

            // found the max component, let's create the graph
            largecomp = creator.create_component((Graph *)g, maxlist, true);
            return largecomp;
        }
        else {
            return NULL;
        }
    } // get_largest_component_graph
}
using namespace std;
/**
 * Populates the provided Graph structure with the necessary information to do TD, visualizations, analysis.
 * Assumes graph_file is in DIMACS format.
 * If needed, stores the largest connected components in DIMACS format at graph_file.giantcomp (and populates Graph from this).
 */
void Graph::create_largestcomponent_graph(char *graph_file, VertexWeightedGraph *&G){
    GraphCreatorFile creator;
    creator.set_file_name(graph_file);
    creator.set_graph_type("DIMACS");
    G = creator.create_weighted_mutable_graph();
    GraphProperties properties;

    if(!G){
        throw(GraphException("Failed to read graph from specified file.\n"));
    }

    if(!properties.is_connected(G)){
        char *bigcompfile = (char *) malloc(100);
        sprintf(bigcompfile, "%s.giantcomp", graph_file);
        G->write_largest_component("DIMACS", "temp_max_comp.dimacs");
        normalize_DIMACS_file("temp_max_comp.dimacs", bigcompfile );
        delete G;
        remove("temp_max_comp.dimacs");
        creator.set_file_name(bigcompfile);
        G = creator.create_weighted_mutable_graph();
        G->set_input_file(bigcompfile);
        free(bigcompfile);
    }
} // create_largestcomponent_graph

