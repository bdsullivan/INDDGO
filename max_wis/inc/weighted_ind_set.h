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

#ifndef __WIS_H__
#define __WIS_H__
#ifdef __PARALLEL__
#include <mpi.h>
#endif
#if HAS_GMP
#include "gmp.h"
#endif

void write_ind_set_model(const char *DIMACS_file, const char *model_file, 
	Graph::WeightedMutableGraph *G);
int create_subgraph_stats(TDTree *T, vector<int_int> *stats);
int compute_weighted_ind_set_table(TDTree *T, int k);
int compute_nonnice_table_parent_child(TDTree *T, int k);
int compute_nonnice_table_standard(TDTree *T, int k);
int compute_nonnice_table_async(TDTree *T, int k);
int compute_nonnice_table_split_bag(TDTree *T, int k);
int compute_nonnice_table(TDTree *T, int k);
int compute_ind_sets_nice(TDTree *T, int k);
int compute_introduce_table(TDTree *T, int k);
int compute_forget_table(TDTree *T, int k);
int compute_join_table(TDTree *T, int k);
int list_ind_sets(TDTree *T, int k, list<int> *bag_subset, list<TDSolution *> *ind_sets);
int hamming_weight(unsigned long long xx);
bool hamming_compare(const int_bigint *x,const int_bigint *y);
void usage(const char *s);
void create_WIS_graph(DP_info *info, Graph::WeightedMutableGraph *&G);
void create_tree_decomposition(DP_info *info, Graph::WeightedMutableGraph *G, TDTree **T, bool suppress_timing);
void create_tree_decomposition(DP_info *info, Graph::WeightedMutableGraph *G, TDTree **T);
int reconstruct_solution(TDTree *T, list<int> *optimal_solution, const char *sol_file);
void print_WIS_results(FILE *stream, TDTree *T, DP_info *info);
int get_optimal_obj(TDTree *T, list<int> *root_intersection, list<int> *root_difference, 
	DP_info *info);
double estimate_memory_usage(TDTree *T, vector<int> *walk, const char *outfile);
double estimate_memory_usage(TDTree *T, vector<int> *walk, bool parent_intersection);
double expected_num_ind_sets(TDTree *T, int k, bool parent_intersection);

#endif  /* __WIS_H__ */
