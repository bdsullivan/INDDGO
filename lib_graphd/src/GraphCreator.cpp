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

#include "GraphCreator.h"
#include "RndNumGen.h"
#include "Debug.h"
#include "GraphUtil.h"
#include "VertexWeightedGraph.h"

#include <typeinfo>

namespace Graph {
    GraphCreator::GraphCreator(){
    }

    /**
     * Initialization of an R-MAT graph with 2^l vertices and <= m edges.
     * Four probabilities in *probs used for quadrants, must sum to 1.
     * lc rng initialized with seed. Setting self_loop true throws another
     * edge if a self loop occurs.
     */
    Graph *GraphCreator::initialize_rmat(int l, int m, double *probs, int seed, bool self_loop){
        try {
            //check neighborhood to account for float errors
            if((probs[0] + probs[1] + probs[2] + probs[3] > 1.001) || (probs[0] + probs[1] + probs[2] + probs[3] < 0.999)){
                fatal_error("initialize_rmat: Probabilities did not sum to 1: %f, %f, %f, %f", probs[0], probs[1], probs[2], probs[3]);
            }
        }
        catch(...){
            fatal_error("initialize_rmat: Error occured when checking probabilities");
        }

        //set up quadrants
        const double q1 = probs[0];
        const double q2 = probs[0] + probs[1];
        const double q3 = probs[0] + probs[1] + probs[2];
        set< pair<unsigned long, unsigned long> > edges;

        //set up rng
        init_lcgrand(0,seed);

        //throw m edges
        for(int i = 0; i < m; i++){
            unsigned long x = 0;
            unsigned long y = 0;
            for(int j = l - 1; j >= 0; j--){
                double r = lcgrand(0);
                if(r < q1){
                }
                else if(r < q2){
                    // pow(2,j)?? 1<<j might be slightly faster...
                    x += (int)pow((double)2,j);
                }
                else if(r < q3){
                    y += (int)pow((double)2,j);
                }
                else {
                    x += (int)pow((double)2,j);
                    y += (int)pow((double)2,j);
                }
            }

            //keep edges in order so we don't have a (x,y) and (y,x) duplicate
            if(x == y){
                if(self_loop){
                    i--;
                }
            }
            else if(y > x){
                edges.insert(make_pair(x,y));
            }
            else {
                edges.insert(make_pair(y,x));
            }
        }

        Graph *g = new Graph((int)pow((double)2,l));
        set< pair<unsigned long, unsigned long> >::const_iterator it;
        for(it = edges.begin(); it != edges.end(); it++){
            g->add_edge(it->first,it->second);
        }

        return g;
    } //initialize_rmat

    /**
     * Initialization of a random intersection graph with n vertices, probs->size() attributes,
     * with probabilities on the attributes. lc rng initialized with seed.
     */
    Graph *GraphCreator::initialize_rig(int n, int seed, list<double> *probs){
        set< pair<unsigned long, unsigned long> > edges;
        init_lcgrand(0, seed);

        //for all the attributes, construct a clique in G
        list<double>::const_iterator it;
        for(it = probs->begin(); it != probs->end(); it++){
            //compute which actors have a given attribute
            vector<int> connect;
            for(unsigned long i = 0; i < n; i++){
                if(lcgrand(0) < *it){
                    connect.push_back(i);
                }
            }
            //attach these actors to the graph
            int conn_size = connect.size();
            for(int i = 0; i < conn_size; i++){
                for(int j = i + 1; j < conn_size; j++){
                    edges.insert(make_pair(connect[i],connect[j]));
                }
            }
        }

        Graph *g = new Graph(n);
        set< pair<unsigned long, unsigned long> >::const_iterator edge_it;
        for(edge_it = edges.begin(); edge_it != edges.end(); edge_it++){
            g->add_edge(edge_it->first,edge_it->second);
        }

        return g;
    } //initialize_rig

    /**
     * Initialization of a k-tree with n vertices. RNG seeded with seed.
     * mg starts with a clique on vertices 0 through k and iteratively
     * adds vertices k+1 up to n-1 with each being adjacent to a
     * randomly chosen k-clique in the existing graph at that step.
     */
    VertexWeightedGraph *GraphCreator::initialize_ktree(int n, int k){
        VertexWeightedGraph *mg = new VertexWeightedGraph(n);

        int i, j;
        int c, s;
        list<int>::iterator it;

        //initialize the first k+1 clique in the adjacency list of vertex k
        //also add the appropriate adjacencies for the vertices < k
        for(i = 0; i <= k; i++){
            for(j = 0; j < i; j++){
                mg->add_edge(i, j);
            }
        }

        //adding vertex i for i= k+1 to n, we randomly choose a k+1 clique
        //(given by an index in [k,i-1]), then choose a vertex to omit
        //(by randomly selecting an index in [0,k]), then add the edges from
        //vertex i to the other k vertices.

        for(i = k + 1; i < n; i++){
            c = rand_int(k, i - 1);
            s = rand_int(0, k);
            print_message(1, "Chose color %d and omission index %d\n", c, s);

            it = mg->nodes[c].nbrs.begin();
            for(j = 0; j < k; j++){
                if(j != s){
                    //mg->nodes[i].nbrs.push_back(*it);
                    mg->add_edge(i, *it);
                }
                ++it;
            }
            //if you skipped a list element, we need to add vertex c
            if(s != k){
                //mg->nodes[i].nbrs.push_back(c);
                mg->add_edge(i, c);
            }
        }         //end for i loop

        GraphProperties prop;
        prop.make_simple(mg);

        // this is not necessary, make_simple already recomputes degrees
        //GraphUtil util;
        //util.recompute_degrees(mg);
        //mg->input_file=NULL;
        return mg;
    } // initialize_ktree

    /**
     * Helper function for Graph(Graph, int) constructor. Assumes wmg-> has already been initialized
     * with appropriate number of nodes, and other variables.  Creates a subgraph of H by selecting
     * random edges for removal and then adding them to the current Graph object.
     */
    VertexWeightedGraph *GraphCreator::create_random_edge_subgraph(VertexWeightedGraph *H, int percent_edges){
        //int n = H->Graph::Graph::capacity;
        int n = H->get_capacity();

        VertexWeightedGraph *wmg = new VertexWeightedGraph(n);
        //set up the Graph allocations and default variables.
        wmg->num_edges = 0;
        wmg->capacity = n;
        wmg->num_nodes = n;

        //if(wmg->capacity > H->weight.size()){
        //   FERROR("%s: too much capacity in wmg for H", __FUNCTION__);
        //  throw GraphException("element is out of bounds\n");
        //}

        for(int i = 0; i < wmg->capacity; i++){
            wmg->weight[i] = H->weight[i];
        }

        //wmg will be simple if H was. Otherwise, it will be unknown.
        if(H->simple == true){
            wmg->simple = true;
        }

        else {
            wmg->simple = false;
        }

        int i;
        //find the number of edges to delete
        int m = H->num_edges;
        int num_removed = m - (int) (floor(
                                         ((double) (percent_edges) * (double) (m)) / (double) (100)));
        DEBUG("Removing %d edges out of %d\n", num_removed, m);
        int selected = 0;
        // An array to keep track of which ones we've selected
        vector<int> chosen_edges(m, 0);
        int curr;
        while(selected < num_removed){
            // An occupancy problem
            curr = rand_int(0, m - 1);
            if(chosen_edges[curr] == 0){
                chosen_edges[curr] = 1;
                selected++;
            }
        }
        // Loop over adjacency list, copying the non-removed edges. We only consider the edges of the form (i,j) with i >=j.
        list<int>::iterator it;
        int e = 0;
        wmg->next_label = H->next_label;
        for(i = 0; i < n; i++){
            wmg->nodes[i].label = H->nodes[i].label;
            if(H->nodes[i].label != -1){
                it = H->nodes[i].nbrs.begin();
                while(it != H->nodes[i].nbrs.end() && *it <= i){
                    //if it's not removed, copy it into G
                    if(chosen_edges[e] == 0){
                        wmg->add_edge(i, *it);
                        //wmg->nodes[i].nbrs.push_back(*it);
                        //wmg->num_edges++;
                    }
                    //increment the iterator and edge counter
                    ++it;
                    e++;
                }
            }
        }

        GraphUtil util;
        util.recompute_degrees(wmg);
        wmg->num_connected_components = -1;
        wmg->graph_type = -1;
        // Not canonical since asymmetric - why is wmg the case?
        wmg->canonical = false;
        return wmg;
    } // create_random_edge_subgraph

    /**
     * Creates a graph G with vertex set equal to the set in the members list.
     * The entries in the members list are the positions in the nodes[] array.
     * The edges in G are created by constructing the subgraph induced
     * by the members list.  If make_simple is true, then the resulting graph
     * is guaranteed to be simple.
     */
    VertexWeightedGraph *GraphCreator::create_induced_subgraph(VertexWeightedGraph *g, list<int> *members, bool make_simple){
        VertexWeightedGraph *wmg = new VertexWeightedGraph(members->size());
        // The caller of g function is responsible for allocating a Graph
        // of the appropriate size - make sure!

        if(wmg->capacity != (int) (((members->size())))){
            fatal_error(
                "%s:  Graph of members->size() must be allocated before calling\ncreate_component\n",
                __FUNCTION__);
        }

        // Sort the members - not necessary, but convenient?
        members->sort();
        // Create a subgraph_labels[] array and a graph_indices[]
        // array

        int *subgraph_labels = new int[members->size()];
        memset(subgraph_labels, GD_UNDEFINED, members->size() * sizeof(int));
        int *graph_indices = new int[g->capacity];
        memset(graph_indices, GD_UNDEFINED, g->capacity * sizeof(int));
        list<int>::iterator i, j;
        int k = 0;
        for(i = members->begin(); i != members->end(); ++i){
            subgraph_labels[k] = g->nodes[*i].label;
            graph_indices[*i] = k;
            k++;
        }

        k = 0;
        int new_i, new_j;
        // Use g for make_simple checking
        char *neighbor_vec = new char[g->capacity];
        bool multiple_edges = false;
        wmg->num_edges = 0;
        for(i = members->begin(); i != members->end(); ++i){
            memset(neighbor_vec, 0, g->capacity * sizeof(char));
            new_i = graph_indices[*i];
            if(new_i == GD_UNDEFINED){
                fatal_error(
                    "%s: Encountered GD_UNDEFINED entry for new_i at position %d in graph_indices array\n",
                    __FUNCTION__, *i);
            }

            wmg->nodes[new_i].label = subgraph_labels[k];
            // CSG - 2/24: added g as we are losing weights when
            // creating components!
            wmg->weight[new_i] = g->weight[*i];
            for(j = g->nodes[*i].nbrs.begin(); j != g->nodes[*i].nbrs.end(); ++j){
                // We have an edge between *i and *j in the original graph
                // g should be stored as an edge between
                // graph_indices[*i] and graph_indices[*j] in the subgraph

                new_j = graph_indices[*j];
                print_message(1, "Translating old edge %d-%d to new edge %d-%d\n",
                              *i, *j, new_i, new_j);

                if(new_j == GD_UNDEFINED){
                    print_message(
                        1,
                        "%s: Encountered GD_UNDEFINED entry for new_j at position %d in graph_indices array\n",
                        __FUNCTION__, *j);
                }
                else {
                    if(new_i <= new_j){
                        if(!make_simple || (neighbor_vec[*j] == 0) ){
                            if(neighbor_vec[*j] == 1){
                                // We are adding a duplicate edge to G - G will not be simple
                                multiple_edges = true;
                            }
                            // CSG added make_simple functionality - just check if we have seen g edge before by checking
                            // the neighbor_vec
                            wmg->num_edges++;
                            wmg->nodes[new_i].nbrs.push_back(new_j);
                            if(new_i != new_j){                             // added check Dec 29
                                wmg->nodes[new_j].nbrs.push_back(new_i);
                            }
                        }
                    }
                }
                neighbor_vec[*j] = 1;
            }
            // k is used for the subgraph labels
            k++;
        }
        // If the original was simple, G will be. If make_simple, then it had better be simple!
        // Otherwise, unknown.
        if((g->simple == GD_TRUE) || make_simple){
            wmg->simple = GD_TRUE;
        }
        else {
            wmg->simple = false;
            wmg->canonical = false;
        }
        if(multiple_edges){
            wmg->simple = false;
            wmg->canonical = false;
        }
        // We don't know number of components
        wmg->num_connected_components = GD_UNDEFINED;
        GraphUtil util;
        util.recompute_degrees(wmg);

        delete[] subgraph_labels;
        delete[] graph_indices;
        delete[] neighbor_vec;

        return wmg;
    } // create_induced_subgraph

    /**
     * Creates the subgraph G induced by the members list where it is known
     * that g subgraph is a connected component.  Sets G->num_components=1
     * although g is not verified!  If make_simple is true, then the component
     * is guaranteed to be simple.
     */
    VertexWeightedGraph *GraphCreator::create_component(VertexWeightedGraph *g,
                                                        list<int> *members, bool make_simple){
        VertexWeightedGraph *wmg;
        wmg = create_induced_subgraph(g, members, make_simple);
        // Set num_connected_components to 1 since g is known - not verified - should it be??
        wmg->num_connected_components = 1;
        return wmg;
    }

    /**
     * Creates a graph G with vertex set equal to the set in the members list.
     * The entries in the members list are the positions in the nodes[] array.
     * The edges in G are created by constructing the subgraph induced
     * by the members list.  If make_simple is true, then the resulting graph
     * is guaranteed to be simple.
     */
    Graph *GraphCreator::create_induced_subgraph(Graph *g, list<int> *members, bool make_simple){
        Graph *wmg = new Graph(members->size());
        // The caller of g function is responsible for allocating a Graph
        // of the appropriate size - make sure!

        if(wmg->capacity != (int) (((members->size())))){
            fatal_error(
                "%s:  Graph of members->size() must be allocated before calling\ncreate_component\n",
                __FUNCTION__);
        }

        // Sort the members - not necessary, but convenient?
        members->sort();
        // Create a subgraph_labels[] array and a graph_indices[]
        // array

        int *subgraph_labels = new int[members->size()];
        memset(subgraph_labels, GD_UNDEFINED, members->size() * sizeof(int));
        int *graph_indices = new int[g->capacity];
        memset(graph_indices, GD_UNDEFINED, g->capacity * sizeof(int));
        list<int>::iterator i, j;
        int k = 0;
        for(i = members->begin(); i != members->end(); ++i){
            subgraph_labels[k] = g->nodes[*i].label;
            graph_indices[*i] = k;
            k++;
        }

        k = 0;
        int new_i, new_j;
        // Use g for make_simple checking
        char *neighbor_vec = new char[g->capacity];
        bool multiple_edges = false;
        wmg->num_edges = 0;
        for(i = members->begin(); i != members->end(); ++i){
            memset(neighbor_vec, 0, g->capacity * sizeof(char));
            new_i = graph_indices[*i];
            if(new_i == GD_UNDEFINED){
                fatal_error(
                    "%s: Encountered GD_UNDEFINED entry for new_i at position %d in graph_indices array\n",
                    __FUNCTION__, *i);
            }

            wmg->nodes[new_i].label = subgraph_labels[k];
            for(j = g->nodes[*i].nbrs.begin(); j != g->nodes[*i].nbrs.end(); ++j){
                // We have an edge between *i and *j in the original graph
                // g should be stored as an edge between
                // graph_indices[*i] and graph_indices[*j] in the subgraph

                new_j = graph_indices[*j];
                print_message(1, "Translating old edge %d-%d to new edge %d-%d\n",
                              *i, *j, new_i, new_j);

                if(new_j == GD_UNDEFINED){
                    print_message(
                        1,
                        "%s: Encountered GD_UNDEFINED entry for new_j at position %d in graph_indices array\n",
                        __FUNCTION__, *j);
                }
                else {
                    if(new_i <= new_j){
                        if(!make_simple || (neighbor_vec[*j] == 0) ){
                            if(neighbor_vec[*j] == 1){
                                // We are adding a duplicate edge to G - G will not be simple
                                multiple_edges = true;
                            }
                            // CSG added make_simple functionality - just check if we have seen g edge before by checking
                            // the neighbor_vec
                            wmg->num_edges++;
                            wmg->nodes[new_i].nbrs.push_back(new_j);
                            if(new_i != new_j){                             // added check Dec 29
                                wmg->nodes[new_j].nbrs.push_back(new_i);
                            }
                        }
                    }
                }
                neighbor_vec[*j] = 1;
            }
            // k is used for the subgraph labels
            k++;
        }
        // If the original was simple, G will be. If make_simple, then it had better be simple!
        // Otherwise, unknown.
        if((g->simple == GD_TRUE) || make_simple){
            wmg->simple = GD_TRUE;
        }
        else {
            wmg->simple = false;
            wmg->canonical = false;
        }
        if(multiple_edges){
            wmg->simple = false;
            wmg->canonical = false;
        }
        // We don't know number of components
        wmg->num_connected_components = GD_UNDEFINED;
        GraphUtil util;
        util.recompute_degrees(wmg);

        delete[] subgraph_labels;
        delete[] graph_indices;
        delete[] neighbor_vec;

        return wmg;
    } // create_induced_subgraph

    /**
     * Creates the subgraph G induced by the members list where it is known
     * that g subgraph is a connected component.  Sets G->num_components=1
     * although g is not verified!  If make_simple is true, then the component
     * is guaranteed to be simple.
     */
    Graph *GraphCreator::create_component(Graph *g, list<int> *members, bool make_simple){
        Graph *comp;
        comp = create_induced_subgraph(g, members, make_simple);
        // Set num_connected_components to 1 since g is known - not verified - should it be??
        comp->num_connected_components = 1;
        return comp;
    }

    /**
     * Fills the provided list with graphs that correspond to the connected
     * components of the graph.  Sets num_connected_components to the correct value
     * and returns g value as well.  For large graphs with many components, it takes
     * an extremely long time to populate the list of Graphs, although finding the
     * components is generally very fast (~3 seconds for 2^20 nodes).
     * If make_simple is true, then each component is guaranteed to be simple.
     */
    list<VertexWeightedGraph *> GraphCreator::create_all_components(VertexWeightedGraph *g, bool make_simple){
        int i;
        list<VertexWeightedGraph *> C;
        // Use label_all_components
        // First find/label the components
        // Now create num_components lists - one for each component
        // NOTE!! Currently we add a copy of g to the list even if the graph has only one component
        // Make a single pass through the components vector and populate the lists in the
        // members[] array
        // Added Dec 28
        // i is in component cnum
        // Now allocate room for each of the graphs
        // Allocate memory for a new graph of the appropriate size
        // Create the component for the i-th list in members
        // Add the address of H to the input list C
        // return the # of components
        // Use find_all_components
        // First find the components - members[i] is a list of the node indices in
        // component i
        GraphUtil util;
        vector<list<int> *> members;
        util.find_all_components(g, &members);
        // Now allocate room for each of the graphs

        VertexWeightedGraph *H;
        for(i = 0; i < g->num_connected_components; i++){
            // Create the component for the i-th list in members
            // Respect the make_simple flag here!
            H = create_component(g, members[i], make_simple);

            // Above is done by create_component function already
            // Add the address of H to the input list of graphs, C
            C.push_back(H);
        }

        // Need to delete the lists that comprise members[]
        vector<list<int> *>::iterator jj;
        for(jj = members.begin(); jj != members.end(); ++jj){
            delete (*jj);
        }

        // return the # of components
        return C;
    } // create_all_components

    /**
     * Fills the provided list with graphs that correspond to the connected
     * components of the graph.  Does so by internally calling recursive functions.
     * Sets num_connected_components to the correct value
     * and returns g value as well.  If make_simple is true, then each component
     * is guaranteed to be simple.
     */
    list<VertexWeightedGraph *> GraphCreator::create_rec_all_components(VertexWeightedGraph *g, bool make_simple){
        int i, cnum = -1;
        list<VertexWeightedGraph *> C;
        GraphUtil util;
        VertexWeightedGraph *H;

        // First find the components - populate a components[] vector containing the labels
        vector<int> components(g->capacity, GD_UNDEFINED);
        g->num_connected_components = util.rec_label_all_components(g, &components);

        // Now create num_components lists - one for each component
        // NOTE!! Currently we add a copy of g to the list even if the graph has only one component
        list<int> *members;
        members = new list<int> [g->num_connected_components];
        for(i = 0; i < g->num_connected_components; i++){
            members[i].clear();
        }

        // Make a single pass through the components vector and populate the lists in the
        // members[] array
        for(i = 0; i < g->capacity; i++){
            if(g->nodes[i].label != GD_UNDEFINED){             // Added Dec 28
                // i is in component cnum
                cnum = components.at(i);
                print_message(100, "cnum =%d\n", cnum);
                print_message(100, "members[%d].size()=%d\n", cnum,
                              members[cnum].size());
                members[cnum].push_back(i);
            }
        }

        // Now allocate room for each of the graphs
        for(i = 0; i < g->num_connected_components; i++){
            //H = new VertexWeightedGraph();
            // Create the component for the i-th list in members
            // Respect the make_simple flag here!
            H = create_component(g, members + i, make_simple);

            //H->num_connected_components=1;
            // Add the address of H to the input list C
            C.push_back(H);
        }

        // return the # of components
        delete[] members;
        return C;
    } // create_rec_all_components

    GraphCreator::~GraphCreator(){
    }
}
