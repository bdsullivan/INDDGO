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
#include "TreeDecomposition.h"

/**
* The constructor for the TDTreeNode class.
*/
TDTreeNode::TDTreeNode()
{
	this->id = GD_UNDEFINED;
	// Initialize hash_table to a NULL pointer
	this->hash_table = NULL;
	this->num_subgraph_edges = 0;
	this->completed_children = 1;

	// Set pointers to DP-related data structures to NULL
	this->nbr_mask_vec = NULL;
	this->child_intersection = NULL;
	this->parent_position = NULL;
	this->vertex_weights = NULL;
	this->num_mask_words = 0;
	this->root = false;

	return;
}

/**
* Destructor for TDTreeNode class.
*/
TDTreeNode::~TDTreeNode()
{
	// Delete the hash_table if it exists
#ifndef __MADNESS__
	if (this->hash_table)
	{
		TDSolution *current_set, *tmp;
		HASH_ITER(hh, this->hash_table, current_set, tmp)
		{
			HASH_DEL(this->hash_table,current_set);
			delete current_set;
		}
		this->hash_table = NULL;
	}
#endif

	if (this->nbr_mask_vec)
	{
		delete this->nbr_mask_vec;
		this->nbr_mask_vec = NULL;
	}

	return;
}

/**
* Assignment operator for TDTreeNode.
*/
TDTreeNode& TDTreeNode::operator=(const TDTreeNode& rhs)
{
	//if (!(&rhs) && this != &rhs)
	if (this != &rhs)
	{
		list<int>::const_iterator ii;
		int i = 0;
		// CSG added this after bag vector removal
		this->bag.clear();
		this->bag_vec.clear();
		for (ii = rhs.bag.begin(); ii != rhs.bag.end(); ++ii)
		{
			this->bag.push_back(*ii);
			this->bag_vec.push_back(*ii);
		}

		this->adj.clear();
		for (ii = rhs.adj.begin(); ii != rhs.adj.end(); ++ii)
			this->adj.push_back(*ii);

		this->id = rhs.id;
		this->completed_children = rhs.completed_children;
		this->num_mask_words = rhs.num_mask_words;
		this->root = rhs.root;

		this->hash_table = rhs.hash_table;
		// copy nbr_mask_vec
		this->nbr_mask_vec = NULL;
		if (rhs.nbr_mask_vec)
		{
			int w = rhs.bag.size();
			if (w > 0)
			{
				this->nbr_mask_vec = new vector<int_bigint *> (w);
				for (i = 0; i < w; ++i)
				{
					this->nbr_mask_vec->at(i) = new int_bigint(
						rhs.num_mask_words);
					this->nbr_mask_vec->at(i) = rhs.nbr_mask_vec->at(i);
					this->nbr_mask_vec->at(i)->k = rhs.nbr_mask_vec->at(i)->k;
				}
			}
		}

		//copy child intersection
		this->child_intersection = NULL;
		if (rhs.child_intersection)
		{
			int num_children = rhs.adj.size() - 1;
			if (num_children > 0)
			{
				this->child_intersection
					= new vector<bigint_t *> (num_children);

				for (i = 0; i < num_children; i++)
				{
					this->child_intersection->at(i) = new bigint_t(
						rhs.num_mask_words);
					this->child_intersection->at(i)
						= rhs.child_intersection->at(i);
				}
			}
		}

		//copy parent position
		this->parent_position = NULL;
		if (rhs.parent_position)
		{
			this->parent_position = new vector<int> ();
			this->parent_position = rhs.parent_position;
		}

		//copy vertices weights
		this->vertex_weights = NULL;
		if (rhs.vertex_weights)
		{
			this->vertex_weights = new vector<int> ();
			this->vertex_weights = rhs.vertex_weights;
		}
	}
	return *this;
}

/**
* Copy constructor for TDTreeNode.
*/
TDTreeNode::TDTreeNode(const TDTreeNode& rhs)
{
	list<int>::const_iterator ii;
	int i = 0;

	this->bag.clear();
	this->bag_vec.clear();

	for (ii = rhs.bag.begin(); ii != rhs.bag.end(); ++ii)
	{
		this->bag.push_back(*ii);
		this->bag_vec.push_back(*ii);
	}

	// Copy the actual list
	this->adj.clear();
	for (ii = rhs.adj.begin(); ii != rhs.adj.end(); ++ii)
		this->adj.push_back(*ii);

	// Id's are the same
	this->id = rhs.id;
	this->completed_children = rhs.completed_children;
	this->num_mask_words = rhs.num_mask_words;
	this->root = rhs.root;

	this->hash_table = rhs.hash_table;

	// copy nbr_mask_vec
	this->nbr_mask_vec = NULL;
	if (rhs.nbr_mask_vec)
	{
		int w = rhs.bag.size();
		if (w > 0)
		{
			this->nbr_mask_vec = new vector<int_bigint *> (w);
			for (i = 0; i < w; ++i)
			{
				this->nbr_mask_vec->at(i) = new int_bigint(rhs.num_mask_words);
				this->nbr_mask_vec->at(i) = rhs.nbr_mask_vec->at(i);
				this->nbr_mask_vec->at(i)->k = rhs.nbr_mask_vec->at(i)->k;
			}
		}
	}

	//copy child intersection
	this->child_intersection = NULL;
	if (rhs.child_intersection)
	{
		int num_children = rhs.adj.size() - 1;
		if (num_children > 0)
		{
			this->child_intersection = new vector<bigint_t *> (num_children);

			for (i = 0; i < num_children; i++)
			{
				this->child_intersection->at(i) = new bigint_t(
					rhs.num_mask_words);
				this->child_intersection->at(i) = rhs.child_intersection->at(i);
			}
		}
	}

	//copy parent position
	this->parent_position = NULL;
	if (rhs.parent_position)
	{
		this->parent_position = new vector<int> ();
		this->parent_position = rhs.parent_position;
	}

	//copy vertices weights
	this->vertex_weights = NULL;
	if (rhs.vertex_weights)
	{
		this->vertex_weights = new vector<int> ();
		this->vertex_weights = rhs.vertex_weights;
	}
}

/**
* Overloading the << operator for the TDTreeNode class.
*/
ostream &operator<<(ostream &output, const TDTreeNode &T)
{
	list<int>::const_iterator ii;
	output << "TreeNode " << T.id << endl;
	// Show the bag
	output << "Bag size: " << T.bag.size() << " ( ";
	for (ii = T.bag.begin(); ii != T.bag.end(); ++ii)
		output << *ii << " ";
	output << ")" << endl;

	// Show the edges in the tree adjacent to this TDTreeNode
	output << "Num_neighbors: " << T.adj.size() << " ( ";
	for (ii = T.adj.begin(); ii != T.adj.end(); ++ii)
		output << *ii << " ";
	output << ")" << endl;
	return output;
}

/**
* Simple function to determine if a tree node is a leaf.
*/
bool TDTreeNode::is_leaf()
{
	if (this->adj.size() == 1)
		return true;
	return false;
}

/**
* Returns the number of children nodes for whom DP is complete.
*/
int TDTreeNode::get_completed_children() const
{
	return this->completed_children;
}


/**
* Sets the number of children nodes for whom DP is complete.
*/
void TDTreeNode::set_completed_children(int i)
{
	this->completed_children = i;
}


/**
* Increases the number of children nodes for whom DP is complete by the input i.
*/
void TDTreeNode::increase_completed_children(int i)
{
	this->completed_children += i;
}

/**
* Sets the number of words required in the TDTreeNode's DP masks. This amount is based
* on the width of the tree since all masks must have the same number of words.
*/
void TDTreeNode::set_num_mask_words(int i)
{
	this->num_mask_words = i;
}

/**
* Returns the number of words in the tree node's mask.
*/
int TDTreeNode::get_num_mask_words() const
{
	return num_mask_words;
}

/**
* Sets a node's root status according to the input flag.
*/
void TDTreeNode::set_root(bool flag)
{
	this->root = flag;
}

/**
* Returns true if a node's root flag is true, false otherwise.
*/
bool TDTreeNode::is_root()
{
	return this->root;
}

