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

#ifndef _TD_TREE_NODE_H_
#define _TD_TREE_NODE_H_
#include <string>
#include "bigint.h"
#include "Log.h"
#include <iostream>
#include <list>
#include "TDSolution.h"

using namespace std;

class bigintp_map_compare
{
public:
	bigintp_map_compare()
	{
	}
	;
	bool operator()(const bigint_t *x, const bigint_t *y)
	{
		for (int i = x->S - 1; i >= 0; i--)
			if (x->words[i] < y->words[i])
				return true;
		// x must be >=y if we get here
		return false;
	}
};

class bigint_map_compare
{
public:
	bigint_map_compare()
	{
	}
	;
	bool operator()(const bigint_t x, const bigint_t y)
	{
		for (int i = x.S - 1; i >= 0; i--)
			if (x.words[i] < y.words[i])
				return true;
		// x must be >=y if we get here
		return false;
	}
};

/**
 * The TDTreeNode class is used to represent the individual nodes
 * in the tree decomposition.
 */
class TDTreeNode
{
	friend ostream &operator<<(ostream &, const TDTreeNode &);

private:
	/**
	 * Flag to keep track whether total independent set table is already computed.
	 */
	bool compute_tbl;

	/**
	 * Keep track of completed children.
	 */
	int completed_children;

	/**
	 * keep num_mask_words in TDTreeNode, want to copy data without access TDTree
	 */
	int num_mask_words;

	/**
	 * flag to mark root node
	 */
	bool root;

public:
	TDTreeNode(); // Constructor - no fixed size arrays here so don't need a parameter
	~TDTreeNode(); // Destructor

	// Copy constructor
	TDTreeNode(const TDTreeNode& rhs);
	// Assignment operator
	TDTreeNode& operator=(const TDTreeNode& rhs);

	/**
	 * Unique id for this node in the tree decomposition.
	 */
	int id;

	/**
	 * List of the indices of adjacent TDTreeNodes in this array.
	 * For rooted decompositions, we adopt the convention that the
	 * first node in this list is the parent and the remainder are
	 * the children.  The unique root node of the tree has itself
	 * as a parent.
	 */
	list<int> adj;

	/**
	 *The bag containing the vertices represented by this TreeNode.
	 */
	list<int> bag;
	// Populate the vector with the sorted bag once constructed
	vector<int> bag_vec;

	/**
	 * Returns true if the tree node is a leaf, false otherwise.
	 */
	bool is_leaf();

	/**
	 * Returns the type of node - TD_LEAF_NODE, TD_INTRODUCE_NODE, TD_FORGET_NODE,
	 * TD_JOIN_NODE, or TD_OTHER_NODE
	 */
	int node_type(int k);

	/**
	 * Hash table to store partial solutions in DP
	 */
	TDSolution *hash_table;

	/**
	 * List to store all independent sets
	 */
	list<bigint_t *> all_ind_sets;

	/**
	 * List to keep values for independent sets
	 */
	list<int> all_ind_set_values;

	/**
	 * Place to store best obj_val found in DP at this tree node
	 */
	int obj_val;

	/**
	 * # of edges in the subgraph induced by this tree node's bag.
	 */
	int num_subgraph_edges;

	/**
	 * A place to store the adjacency matrix in a convenient format
	 * for subsequent processing. Set to NULL by the constructor and only
	 * allocated when needed.
	 */
	vector<int_bigint *> *nbr_mask_vec;

	/** 
	 * A place to store the intersection of the children's bags with this treenode's
	 * bag.  In particular, child_intersection[i] has bit j set if the j-th entry in
	 * the bag of the i-th child is contained in this tree node's bag.
	 */
	vector<bigint_t *> *child_intersection;

	/**
	 * A place to store the position in the parent bag of every vertex in this
	 * tree node's bag. For example, if vertex 93 is the 3rd vertex in this tree node's
	 * bag, and vertex 93 is the 5th vertex in the parent's bag, then
	 * parent_position[3]=5;
	 */
	vector<int> *parent_position;

	/**
	 * A place to store weights for the vertices in the bag.
	 */
	vector<int> *vertex_weights;

	/** 
	 * Wrapper to delete the adj mat
	 */
	void free_data();

	/**
	 *  counting completed children for asynchronous table update
	 */
	int get_completed_children() const;

	void set_completed_children(int i);

	void increase_completed_children(int i = 1);

	/**
	 *  set num_mask_words in TDTreeNodes
	 */
	void set_num_mask_words(int i);

	int get_num_mask_words() const;

	/**
	 * set root node
	 */
	void set_root(bool flag);

	/**
	 * check for root node
	 */
	bool is_root();
};

#endif

