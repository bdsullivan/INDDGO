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

#ifndef __TD_DP_H__
#define __TD_DP_H__

#ifdef __PARALLEL__
#include <mpi.h>
#endif

#include "uthash.h"

struct key_value
{
	int key, value;
	UT_hash_handle hh;
};


// Data structure to keep track of options for TD-based Dynamic Programming
class DP_info
{
public:
	DP_info();
	~DP_info();

	char *DIMACS_file, *model_file, *ord_file, *gviz_file, 
		*scotch_ord_file, *tree_infile, *tree_outfile, *sol_file, *eorder;

	bool lower_bounds, has_graph, read_ordering, read_scotch_ordering,
		write_mod, solve_mip, gviz, nice,
		superetree, gavril, BK, check, nonniceDP, read_tree, write_tree, 
		write_ordered_tree,
		make_histogram, decompose_only, parent_child, make_nice,
		very_verbose, free_children, verbose, no_reconstruct, split_bag,
		asp_root, child_root, refine_td, DFS, async_tbl, mem_est, width, fix_DIMACS, 
		parmetis, pbag;

	int elim_order_type, root_node, start_v, max_width;
	// CSG adding orig_ fields to handle stats for refined tree
	unsigned long long total_table_entries, total_pc_table_entries, 
		orig_total_table_entries, orig_total_pc_table_entries;
	double avg_pc_proportion;
	// A vector to store additional information related to the DP and the tree nodes
	int *aux_info;

	// statistics on the DP computation
	double leaf_time, introduce_time, forget_time, join_time, nonleaf_time, reconstruct_time; 
	clock_t start, stop;
	int opt_obj;
	int process_DP_info(int num_args, char **args);
	double mem_estimate;
};

void DP_prepare_node(TDTree *T, int k);
int DP_construct_adj_matrix(TDTree *T, int k);
void DP_free_data(TDTree *T, int k);
int DP_create_children_intersections(TDTree *T, int k);
int DP_create_parent_positions(TDTree *T, int k);
int DP_create_vertex_weights(TDTree *T, int k);



#endif

