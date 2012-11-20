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

#ifndef __VIZ_H__
#define __VIZ_H__
#ifdef __PARALLEL__
#include <mpi.h>
#endif
#if HAS_GMP
#include "gmp.h"
#endif


// Data structure to keep track of options for TD-based Dynamic Programming
class TD_info
{
public:
	TD_info();
	~TD_info();


	char *DIMACS_file, *ord_file, *gviz_outfile, 
	  *scotch_ord_file, *tree_infile, *tree_outfile, *score_infile; 

	bool has_graph, read_ordering, read_scotch_ordering,
	  read_tree, write_tree, write_ordered_tree,
	  make_histogram, make_nice, use_scores, subtree, log_range;

	clock_t start, stop;
	int td_alg, elim_order_type, start_v, viz_style, highlight_node;

	void process_TD_info(int num_args, char **args, void (*usage)(const char *str));
};

void usage(const char *s);
void create_tree_decomposition(TD_info *info, Graph::WeightedMutableGraph *G, TDTree **T);

/*These are generic and don't really belong in viz.*/

void form_td(int td_alg, TDTree**T, vector<int>* ordering);

#endif  /* __VIZ_H__ */
