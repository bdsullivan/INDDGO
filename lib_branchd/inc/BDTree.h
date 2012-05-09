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

#ifndef _BD_TREE_H_
#define _BD_TREE_H_

#define NODE_SPLIT  1
#define EDGE_SPLIT  2


class BDTree
{
    friend ostream &operator<<(ostream &, const BDTree &);

public:
    BDTree(Graph *G);       // Constructor where G has n nodes, m edges
    ~BDTree();              // Destructor

    // Flag to state whether we actually have a valid BD (requires ternary tree)
    bool is_valid;
    bool is_rooted;

    int num_nodes;
    int num_edges;
    int num_interior_nodes;
    
    BDTreeNode *nodes; // We know the eventual tree will have 2m-2 nodes and 2m-3 edges
    BDTreeEdge *edges;     

    Graph *G;   // Pointer to the original graph being decomposed

    // Current value of the largest middle set - will represent width when BD is valid
    int max_middle_set_size; 

    void write_graphviz_file(bool spline, const char *gviz_file, bool edge_labels, bool thickness);

    // Create an initial "star" configuration where we have a single central node, 
    //and then m leaf nodes, corresponding to each edge in the graph G
    void create_star();

    // Splits node u into two new nodes and updates the tree accordingly
    // The list A is a subset of indices of neighboring nodes or edges that are
    // to remain "on the same side" of the split.  The value of type should be either
    // NODE_SPLIT or EDGE_SPLIT depending on what A represents.  Returns the index of the
    // new node added to the tree.   Sets the value of new_ms_size to be the size of the 
    // middle set of the newly added edge.
    // Returns the index of the newly created node. Sets new_e_index to index of newly added tree edge.
    int split_node(int u, list<int> *A, int type, int *new_ms_size);
    int split_node(int u, list<int> *A, int type, int *new_ms_size, int *new_e_index);

    // Takes the edges adjacent to node u and creates two new nodes, v and w
    // that are both adjacent to u.  The list A is a list of neighboring LEAF nodes or edges
    // leading to LEAF nodes that will now be adjacent to v.  The remaining LEAF nodes or edges
    // not in A will be adjacent to v.  It is assumed that u is incident to exactly 
    // one edge that does not lead to a leaf.  The value of type should be either
    // NODE_SPLIT or EDGE_SPLIT depending on what A represents.  Returns the index of the
    // new node added to the tree.  Sets the value of new_ms_size1 and new_ms_size2 to be
    // the size of the middle set of the newly added edges u-v and u-w. n1 and n2 are set to
    // be the indices of the two new nodes
    int spawn_nodes(int u, list<int> *A, int type, int *new_ms_size1, int *new_ms_size2, 
        int *n1, int *n2);

    // Returns the size of the middle set if the provided split were to be performed
    int split_size(int u, list<int> *A, int type);

    // Assumes that you want to split u with A and B being edge lists on opposite sides
    // of the split.  Letting t=|nbrs(u)|-|A|-|B|, the function exhausts all 2^t possible
    // ways of adding the remaining edges into the split.  Returns the size of the lowest
    // middle set found, and fills the list C with the edge indices of the corresponding
    // set.
    int best_partition(int u, list<int> *A, list<int> *B, list<int> *C);

    // Considers the edges incident to u and then repeatedly tries to 
    // push new separations until none can be found.
    int push(int u);
   
    //find & run all bfpushes arising if you do a push at u.
    //Modified Nov 2 to be more efficient at checking nodes which may have a push.
    //if newe_only is set to true, we only consider pushes involving the newly created edge at that vertex, and 
    // u_edge should be set to the desired edge for the initial vertex u.
    //Returns the number of splits (pushes) enacted.
    int doall_bfpushes(int u, bool newe_only, int u_edge);
    
    //brute force version of pushing - this just looks for a single separation 
    // at u which separates 2 edges from the rest, satisfies the pushing inequality, and makes that split. 
    // If e1 is specified, only pairs of edges including e1 are considered for the split.
    // If a split is found, returns true and fills in newv, newe with indices of newly created tree vertex/edge, respectively.
    // Otherwise, returns false, sets values to -1.
    // User is responsible for looping over tree repeatedly looking for more.
    bool bf_push(int u, int *newv, int *newe); 
    bool bf_push(int u, int e1, int *newv, int *newe);

    // Returns true if the pushing ineq. is satisfied by edges e1 and e2 incident to u
    bool verify_push(int u, int e1, int e2);
   
    // Useful for small graphs only - does brute force search...
    int find_greedy_split(int u, list<int> *X);

    // Constructs the graph representing the middle sets of all edges incident to v
    // as described in Cook-Seymour TSP paper
    Graph *construct_middle_set_graph(int v, vector<int> *perm, vector<int> *iperm, bool find_CRS);

    // Fills A and B with the indices of the nodes implied by the edge partition created
    // via the Cook-Seymour eigenvector split.  A and B each have \ceil{|D|*p} elements where
    // |D| is the # of edges incident to v, and p is some number in (0,.5].
    void eigenvector_split(int v, list<int> *A, list<int> *B, double p, bool fill_extra);
    void eigenvector_leaf_split(int v, list<int> *A, list<int> *B, double p, bool fill_extra);

    // Just fills the A with a sorting of the BDTreeEdges adjacent to v
    void eigenvector_list(int v, list<int> *A);

    /**
     * Functions to partition edges at a vertex of the branch 
     * decomposition into (X,Y) after using eigenvector (or other) heuristic
     * to classify subsets (A,B) so that A \subset X and B \subset Y. 
     * Specify the node v in the branch decomposition and 
     * lists of the indices of the edges [in BDTreeEdge array] already classified in A,B. 
     * X = A \cup CA and Y = B \cup CB are the desired partition, and 
     * (A \cup B) \cap (CA \cup CB) is empty. 
     */
    void spawn_maxflow_partition(int v, list<int> *A, list<int> *B, list<int>*CA, list<int> *CB, 
        goblinController *CT);

    int split_maxflow_partition(int v, list<int> *A, list<int> *B, list<int>*CA, list<int> *CB, 
        goblinController *CT);

    /**
     * Helper function to check for a three separation in the branch decomposition
     * at a specified vertex v. If one is found, X and Y are populated 
     * with edge indices so that |M(X,Y)| = 3.
     * returns true if one is found, false otherwise.
     */ 
    bool three_separation(int v, list<int> *X, list<int> *Y, goblinController *CT);

    /**
     * Helper function to check for a two separation in the branch decomposition
     * at a specified vertex v. If one is found, X and Y are populated 
     * with edge indices so that |M(X,Y)| = 2.
     * returns true if one is found, false otherwise.
     */ 
    bool two_separation(int v, list<int> *X, list<int> *Y, goblinController *CT);
    
    /**
     * Removes all pairs of vertices u,v from the graph and returns true
     * if the removal of such a pair results in 2 or more components, and false 
     * otherwise. This is finding all 2-separations at the initial star.
     */
    bool brute_force_two_sep();

    /**
     * Removes all triples of vertices u,v,w from the graph and returns true
     * if the removal of such a triple results in 2 or more components, and false 
     * otherwise. This is finding all 3-separations at the initial star.
     */
    bool brute_force_three_sep();
    
    /**
     * specialized middle_set_union allows you to record what new vertices are in the union
     * that you have not seen before (where prior viewings are recorded by W)
     * E is a list of indices into the BDTreeEdge array.
     * V will be filled with a sorted list of the vertices that are
     * in the middle set of at least one edge e\in E.
     * If W != NULL, W[v] will be set to true for all v in V. 
     * newv counts the number of these found vertices for which W[v] 
     * was false at the start of this function.
     */
    int middle_set_union(list<int> *E, list<int> *V, vector<bool>* W, int *newv);

    // E is a list of indices into the BDTreeEdge array
    // V will be filled with a sorted list of the vertices that are
    // in the middle set of at least one edge e\in E.
    int middle_set_union(list<int> *E, list<int> *V);
    
    

    /**
     * Roots the tree by picking some edge a-b and then adding new nodes 
     * to the tree r and s, and new edges a-s, s-b, s-r, and removing the 
     * edge a-b.  Returns the index of the new root node r
     *           r
     *           |
     * a-b ==> a-s-b
     *     
     * Furthermore, we set the edge endpoints throughout the tree so that 
     * the 'end' is nearer the root, effectively forcing all edges to point 
     * towards the root
     */
    int root(int edge_index);
    int root_node;

    // Assuming the tree is rooted, this fills the subtree list with the indices of the
    // nodes beneath u in the tree.
    void find_subtree(int u, list<int> *subtree);

    // Calculate # of leaf neighbors of u
    int calculate_num_leaf_nbrs(int u);

    // Reads a branch decomnposition from a file created by Bill Cook's bwidth binary
    bool read_bwidth_file(const char *infile);
};    

#endif

