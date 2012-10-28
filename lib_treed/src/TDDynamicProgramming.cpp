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
* Constructor for DP_info. Sets all options to false and strings to NULL.
*/
DP_info::DP_info()
{
	this->DIMACS_file = this->model_file = this->ord_file = this->gviz_file
	  = this->scotch_ord_file = this->tree_infile = this->tree_outfile
	  = this->sol_file = this->eorder = NULL;
	this->pbag = false;
	this->parmetis = false;
	this->fix_DIMACS=false;
	this->has_graph = false;
	this->read_ordering = false;
	this->read_scotch_ordering = false;
	this->write_mod = false;
	this->solve_mip = false;
	this->gviz = false;
	this->nice = false;
	this->gavril = false;
	this->superetree = false;
	this->BK = false;
	this->check = false;
	this->nonniceDP = false;
	this->read_tree = false;
	this->write_tree = false;
	this->make_histogram = false;
	this->decompose_only = false;
	this->parent_child = false;
	this->make_nice = false;
	this->verbose = false;
	this->free_children = false;
	this->very_verbose = false;
	this->no_reconstruct = false;
	this->split_bag = false;
	this->asp_root = false;
	this->child_root = false;
	this->DFS = false;
	this->refine_td = false;
	this->write_ordered_tree = false;
	this->async_tbl = false;
	this->mem_est=false;
	this->start_v = GD_UNDEFINED;
	this->elim_order_type = GD_UNDEFINED;
	this->root_node = 0;
	this->total_table_entries = this->total_pc_table_entries = 0;
	this->width = false;
	// set timers to 0
	this->leaf_time = this->introduce_time = this->forget_time
			= this->join_time = this->nonleaf_time = this->reconstruct_time = 0;
	this->mem_est=0;
	this->max_width=999999;
	this->lower_bounds=false;
	// set obj val to 0
	opt_obj = 0;
	aux_info = NULL;
};

/**
* DP_info destructor.
*/
DP_info::~DP_info()
{
	if (this->sol_file)
		delete[] sol_file;
	if (this->aux_info)
		delete[] aux_info;
};


/**
* Function to prepare a node for dynamic programming by constructing the
* adjacency matrix for the induced subgraph and calculating information about
* tree node k's interaction with its children and parent tree nodes.
*/
void DP_prepare_node(TDTree *T, int k)
{
	DP_construct_adj_matrix(T, k);
	DP_create_children_intersections(T, k);
	DP_create_parent_positions(T, k);
	DP_create_vertex_weights(T, k);

	T->tree_nodes[k]->obj_val = 0;
	T->tree_nodes[k]->set_num_mask_words(T->num_mask_words);
}

/**
* Function to create the adjacency matrix for the subgraph induced by
* tree node k's bag.
*/
int DP_construct_adj_matrix(TDTree *T, int k)
{
	int i;
	list<int>::iterator ii, jj;

	// Already constructed - return -1;
	if (T->tree_nodes[k]->nbr_mask_vec)
		return -1;
	// Create a list of neighbor masks (effectively creating a w x w symmetric adj. matrix)
	int w = (int) T->tree_nodes[k]->bag.size();
	// Don't use resize
	T->tree_nodes[k]->nbr_mask_vec = new vector<int_bigint *> (w);
	// Create a tiny hash table for the bag
	key_value *my_table = NULL, *hash_ptr, *temp;
	for (i = 0; i < w; ++i)
	{
		T->tree_nodes[k]->nbr_mask_vec->at(i) = new int_bigint(T->num_mask_words);
		T->tree_nodes[k]->nbr_mask_vec->at(i)->k = i;
		key_value *next_entry = (key_value *) malloc(sizeof(key_value));
		next_entry->key = T->tree_nodes[k]->bag_vec[i];
		next_entry->value = i;
		HASH_ADD_INT(my_table,key,next_entry);
	}

	vector<int> I(w);
	vector<int>::iterator uu, vv;
	Graph::Node *n1;
	list<int> *nbr_list;
	for (i = 0; i < w; i++)
	{
		// Find the set of nbr's of the i-th node in this bag who are also in this bag

		// This is the set_intersection call with the or'ing done on the fly instead of in a
		// separate loop
		uu = T->tree_nodes[k]->bag_vec.begin();
		vv = T->tree_nodes[k]->bag_vec.end();
		int max_in_bag = T->tree_nodes[k]->bag_vec.back();
		n1=T->G->get_node(T->tree_nodes[k]->bag_vec[i]);
		nbr_list=n1->get_nbrs_ptr();
		ii = nbr_list->begin();
		jj = nbr_list->end();
		int max_nbr = nbr_list->back();
		
		while (uu != vv && ii != jj)
		{
			// Break if an additional entry in the intersection is impossible
			if (*ii > max_in_bag || *uu > max_nbr)
				break;
			if (*uu < *ii)
				++uu;
			else
			{
				if (*ii < *uu)
					++ii;
				else
				{
					HASH_FIND_INT(my_table,&(*uu),hash_ptr);
					if (!hash_ptr)
						fatal_error("%s:  WTF\n", __FUNCTION__);
					T->tree_nodes[k]->nbr_mask_vec->at(i)->w->or_bit(
							hash_ptr->value);
					T->tree_nodes[k]->num_subgraph_edges++;
					uu++;
					ii++;
				}
			}
		}
	}

	T->tree_nodes[k]->num_subgraph_edges /= 2;
	// Free the table
	HASH_ITER(hh,my_table,hash_ptr, temp)
	{
		HASH_DEL(my_table,hash_ptr);
		free(hash_ptr);
	}

	return T->tree_nodes[k]->num_subgraph_edges;
}

/**
* Function to free the dynamic programming-specific structures inside
* tree node k.
*/
void DP_free_data(TDTree *T, int k)
{
	if (T->tree_nodes[k]->nbr_mask_vec)
	{
		for (unsigned int i = 0; i < T->tree_nodes[k]->nbr_mask_vec->size(); i++)
			delete T->tree_nodes[k]->nbr_mask_vec->at(i);
		delete T->tree_nodes[k]->nbr_mask_vec;
		T->tree_nodes[k]->nbr_mask_vec = NULL;
	}

	if (T->tree_nodes[k]->child_intersection)
	{
		int i;
		for (i = 0; i < (int) T->tree_nodes[k]->child_intersection->size(); i++)
		{
			delete T->tree_nodes[k]->child_intersection->at(i);
			T->tree_nodes[k]->child_intersection->at(i) = NULL;
		}
		delete T->tree_nodes[k]->child_intersection;
		T->tree_nodes[k]->child_intersection = NULL;
	}

	if (T->tree_nodes[k]->parent_position)
	{
		delete T->tree_nodes[k]->parent_position;
		T->tree_nodes[k]->parent_position = NULL;
	}

	if (T->tree_nodes[k]->vertex_weights)
	{
		delete T->tree_nodes[k]->vertex_weights;
		T->tree_nodes[k]->vertex_weights = NULL;
	}
	return;
}

/**
* Function to create masks that describe the intersection between
* tree node k and his children. Stored in bigint_t's.
*/
int DP_create_children_intersections(TDTree *T, int k)
{
	if (T->tree_nodes[k]->adj.size() == 1)
		// child
		return 0;

	if (T->tree_nodes[k]->child_intersection)
		// Already allocated
		return -1;

	int i = 0, j, w = (int) T->tree_nodes[k]->bag.size();
	vector<int> intersection(w);
	vector<int>::iterator uu, vv;
	list<int>::iterator ii, jj;
	int num_children = T->tree_nodes[k]->adj.size() - 1;
	T->tree_nodes[k]->child_intersection
			= new vector<bigint_t *> (num_children);

	ii = T->tree_nodes[k]->adj.begin();
	++ii;
	for (; ii != T->tree_nodes[k]->adj.end(); ++ii)
	{
		T->tree_nodes[k]->child_intersection->at(i) = new bigint_t(
				T->num_mask_words);
		vv = set_intersection(T->tree_nodes[k]->bag.begin(),
				T->tree_nodes[k]->bag.end(), T->tree_nodes[*ii]->bag.begin(),
				T->tree_nodes[*ii]->bag.end(), intersection.begin());
		// Now go through the intersection and k's bag simultaneously to populate the positions
		jj = T->tree_nodes[k]->bag.begin();
		j = 0;
		for (uu = intersection.begin(); uu != vv; ++uu)
		{
			while (*uu != *jj)
			{
				++jj;
				j++;
			}
			T->tree_nodes[k]->child_intersection->at(i)->set_bit(j);
		}
		i++;
	}
	return num_children;
}

/**
* Function to populate the parent_postion array that gives the position
* in the parent's bag of every vertex j in k's bag.
*/
int DP_create_parent_positions(TDTree *T, int k)
{
	if (T->tree_nodes[k]->parent_position)
	{
		//fprintf(stderr,"Already allocated\n");
		// Already allocated
		return -1;
	}

	int i = 0, j, m, w = (int) T->tree_nodes[k]->bag.size(), parent =
			T->tree_nodes[k]->adj.front();
	vector<int> intersection(w);
	vector<int>::iterator uu, vv;
	list<int>::iterator ii, jj;

	T->tree_nodes[k]->parent_position = new vector<int> (w);
	for (i = 0; i < w; i++)
		T->tree_nodes[k]->parent_position->at(i) = -1;

	vv = set_intersection(T->tree_nodes[k]->bag.begin(),
			T->tree_nodes[k]->bag.end(), T->tree_nodes[parent]->bag.begin(),
			T->tree_nodes[parent]->bag.end(), intersection.begin());
	// Now go through the intersection and k's bag simultaneously to populate the positions
	ii = T->tree_nodes[k]->bag.begin();
	jj = T->tree_nodes[parent]->bag.begin();
	j = m = 0;
	for (uu = intersection.begin(); uu != vv; ++uu)
	{
		// Find *uu in tree node k bag
		while (*uu != *ii)
		{
			++ii;
			j++;
		}
		// Find *uu in parent tree node bag
		while (*uu != *jj)
		{
			++jj;
			m++;
		}
		// Record results of search
		T->tree_nodes[k]->parent_position->at(j) = m;
	}
	return 1;


}

/**
* Function to populate the tree node's internal vertex_weights vector.
* vertex_weights[i] is the weight of the i'th vertex in tree node k's
* sorted bag.
*/
int DP_create_vertex_weights(TDTree *T, int k)
{
	if (T->tree_nodes[k]->vertex_weights)
	{
		// Already allocated
		return -1;
	}

	int i = 0, w = (int) T->tree_nodes[k]->bag.size();
	T->tree_nodes[k]->vertex_weights = new vector<int> (w, 0);

	vector<int> weights = T->G->get_weight();
	for (i = 0; i < w; i++)
	{
		T->tree_nodes[k]->vertex_weights->at(i)
				= weights[T->tree_nodes[k]->bag_vec[i]];
	}

	return 1;
}

