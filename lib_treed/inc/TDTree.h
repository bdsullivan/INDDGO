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

#ifndef _TD_TREE_H_
#define _TD_TREE_H_
#ifdef __PARALLEL__
#include <mpi.h>
#endif
#include "WeightedMutableGraph.h"
#include "TDTreeNode.h"

/**
* The TDTree class is used to represent a tree decomposition for a Graph.
*/
class TDTree
{
    friend ostream &operator<<(ostream &, const TDTree &);
	//friend class Graph::WeightedMutableGraph;

public:
    TDTree(Graph::WeightedMutableGraph *);                  // Constructor given an input graph
    TDTree();                          // Default constructor
    ~TDTree();                         // Destructor

    // Copy constructor
    TDTree(const TDTree& rhs);
    
	// Assignment operator
    TDTree& operator=(const TDTree& rhs);

    /**
    * The index of the root node in the tree.
    */
    int root_node;

    /**
    * Boolean flag to indicated whether tree is rooted or not.
    */
    bool rooted;

    /**
    * Roots the tree at node k. Reorders adj. lists so parent is first, followed by children.
    */
    void root(int k);

    /**
     * Finds & sets a root node to minimize the maximum number of tables
     * required to be stored at any given time (assuming no reconstruction is required).
     * Reference: "Memory Requirements for Table Computations" by Aspvall et al.
     * algorithm needs to be one of the following: 
     * TD_ROOT_ASPVALL - original tabreq() function assumes child tables are deleted sequentially due to DFS
     * TD_ROOT_ALLCHILDREN - tries to account for all child tables needing to exist simultaneously to create parent.
     */
    int root_mintables(int algorithm);


    /**
    * Pointer to the graph being decomposed
    */
    Graph::WeightedMutableGraph *G;

    /**
    * Place for the file containing the graph G.
    */
    char *graph_file;                  

    /** 
    * A vector of pointers to TDTreeNodes that make up the graph.
    * We use a vector instead of a list since we need
    * quick access to node j in the tree. These pointers may be NULL if tree has been modified
    * after creation, and tree_nodes.size() should never be used as a replacement for num_tree_nodes.
    */
    vector<TDTreeNode *> tree_nodes;   

    /**
    * Vector of G->n elements where "element i" is a list of the treenode
    * indices whose bags contain node i from G.
    */
    list<int> *node_locations;

    /**
    * The number of nodes in the tree.
    */
    int num_tree_nodes;

	/**
	* Use this to keep track of DP progress.
	*/
	int num_tree_nodes_processed;

	/**
    * The number of leafs in the tree.
    */
	int num_leafs;

    /**
    * True if the tree is nice, false otherwise.
    */
    bool nice;

    // Fills up the node_locations array by searching each tree node's bag 
    // for vertices
    void record_node_locations();

    // Constructs a tree decomposition using Gavril's algorithm
    // with the provided elimination order
    void construct_gavril(vector<int> *elim_order); 

    // Constructs a nice tree decomposition using Kloks' algorithm
    // with the provided elimination order and tree-width
    // Boolean indicates whether or not to push new node creation down the tree when 
    // possible when encountering a one-child node.
    void construct_knice(vector<int> *elim_order, int k, bool descend_one); 
    
    //Converts a given rooted tree decomposition into a nice TD of equal width.
    void make_nice();
    //helper functions for make_nice
    int replace_by_join_tree(int old_node, int max_id);
    int expand_symdiff(int old_node, int max_id);
    //Tree refinement functions
    void remove_duplicate_bags();
    void remove_subset_bags();
    void refine();

    // Constructs a tree decomposition using Bodlaender-Koster's "Algorithm 2"
    // with the provided elimination order
    void construct_BK(vector<int> *elim_order);

    // Finds the path between start and end in the tree
    void find_path(int start, int end, list<int> *path);
    void find_path(int start, int end, int *pred);

    // Verifies that the decomposition is valid
    bool verify();

    // Verifies that the decomposition is nice
    bool is_nice();

    // Writes the decomposition to a file in (an extended) DIMACS format
    void write_DIMACS_file(const char *DIMACS_file);
    void write_DIMACS_file(const char *DIMACS_file, bool preordered);

    // Reads a decomposition from the file - returns the width
    int read_DIMACS_file(const char *DIMACS_file);

    // Writes the decomposition to a file for viewing via graphviz
    void write_graphviz_file(bool spline, const char *GVIZ_file, int style);


    void write_scored_graphviz_file(bool spline, const char *GVIZ_file, const vector<double>& color_vector, const double& max_color, const double& min_color,int style, const bool log_range=false);
    void write_scored_graphviz_file(bool spline, const char *GVIZ_file, const vector<double>& color_vector, int style, const bool log_range=false);

    // same as above, but highlights the subtree containing v.
    void highlight_subtree_gviz(int v, const char *GVIZ_file, int style);
    void highlight_subtree_scored_graphviz_file(int v, const char *GVIZ_file, const vector<double>& color_vector);
    void highlight_subtree_scored_graphviz_file(int v, const char *GVIZ_file, const vector<double>& color_vector, const double max_color, const double min_color);


    int compute_width();

    int width, refined_width;

    void pre_order_walk(vector<int> *walk);
    void post_order_walk(vector<int> *walk);
	void post_order_walk_DFS(vector<int> *walk);

	int compute_table(int(*table_function)(TDTree *T, int k), int k);
    int get_node_type(int k);


    // Fills the TDTree S with the subtree containing the nodes whose indices are contained
    // in the indices list.  Returns false if there is an error (non-connected nodes, etc.), and
    // returns true otherwise
    bool create_subtree(list<int> *indices, TDTree *S);

    // Resets to an empty tree
    void clear();

    
    // The number of 64-bit words in the masks contained in the TDTable
    int num_mask_words;

    /*
     * For rooted tree decompositions only (fatal error if called on unrooted).
     * Populates the vector B (which will be sized to equal the number of bags in the TD)
     * with B[i] = the number of vertices in bag i which are also in the bags of at least 2 children of i. 
     * If i is a leaf, B[i] is the size of the bag (aka the width).
     */
    int find_child_boundaries(vector<int> *B);
    
    /*
     * Deletes all the memory associated with the children of tree node k and sets
     * the relevant pointers in the tree_node vector to NULL.
     */
    void free_children(int k);

	void free_table(int k);

	void fill_bag_vecs();

	/*
	 * Goes through the tree and removes all nodes from the list S from the bags
	 * of all tree nodes */
	void refine_tree(list<int> *S);

	class DP_info *info;
#ifdef __PARALLEL__
    int size;
    int rank;
    int requester_rank;
    int request_size;
    int pool_size;

    int msk_send_count;
    int wgh_send_count;

    int tnum;
    int flag_recv_done;
    int entries;

    vector<int> *tbl_loc;
    vector<int> *tbl_entries;

    vector<int> *request_count;
    vector<int> *rcount;

    vector<int> *request_expect;
    vector<int> *node_sizes;
    vector<int> *iter_count;
    vector<int> *nodes_with_tbl;

    vector<MPI_Request> *requests;

    vector<MPI_Request> *msk_send_requests;
    vector<MPI_Request> *msk_recv_requests;
    vector<MPI_Request> *wgh_send_requests;
    vector<MPI_Request> *wgh_recv_requests;

    vector<MPI_Request> *tbl_loc_requests;
    vector<MPI_Request> *tbl_size_requests;
    vector<MPI_Request> *compute_term_requests;
    vector<MPI_Request> *storage_term_requests;
    vector<MPI_Request> *child_free_requests;

    vector< vector<unsigned long long> *> msk_send_pool;
    vector< vector<unsigned long long> *> msk_recv_pool;
    vector< vector<int> *> wgh_send_pool;
    vector< vector<int> *> wgh_recv_pool;

    vector<pthread_t> *threads;

    vector<bigint_t *> *start_masks;
    vector<bigint_t *> *end_masks;
    vector<bigint_t> *indset_masks;
    vector<int> *ind_wgh;
    vector<int> *tbl_size;
    bigint_t *iteration_counter;
    
    int head_node;
    int tmp;
    int tblloc;
    
    vector<int> compute_nodes;
    vector<int> storage_nodes;
    vector<int> *allnodes;
    vector<int> *all_table_sizes;
    int storage_nodes_size;
    int compute_nodes_size;
    int del_ch;
    
    MPI_Comm *compute_nodes_comm;

#endif
	void compute_parent_child_intersections(vector<int> *intersection_sizes);	
    //Find subtree "below" a vertex v in a rooted TD.
    void compute_subtree(int v, list<int> *subtree);
    void compute_subtree_vertices(int v, list<int> *vertex_list);

protected:

    // Remove and add edges from the tree
    bool remove_tree_edge(int u, int v);
    bool add_tree_edge(int u, int v);



	
};

#endif



