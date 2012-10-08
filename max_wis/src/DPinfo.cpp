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
* Function to process the command line args and populate the
* DP_info structure.
*/
int DP_info::process_DP_info(int num_args, char **args)
{
  int num_td_methods = 0;
	if (num_args == 1 || (num_args == 2 && strcmp(args[1], "-h") == 0)
			|| (num_args == 2 && strcmp(args[1], "--help") == 0) || (num_args
										 == 2 && strcmp(args[1], "--h") == 0) || (num_args == 2 && strcmp(args[1], "-help")==0))
	{
	  return -1;//calling function needs to call usage() function.
	}

	for (int i = 0; i < num_args; i++)
	{
		if (strcmp(args[i], "-s") == 0)
			this->start_v = atoi(args[i + 1]);
		if (strcmp(args[i], "-decompose_only") == 0)
			this->decompose_only = true;
		if (strcmp(args[i], "-fix_DIMACS")==0)
			this->fix_DIMACS=true;
		// ELIM ORDER LOWER BOUNDS
		if (strcmp(args[i], "-lower_bounds")==0)
		        this->lower_bounds=true;
		// ELIM ORDER TYPES
		if (strcmp(args[i], "-lexm") == 0)
			this->elim_order_type = GD_LEXM_BFS;
		if (strcmp(args[i], "-mcsm") == 0)
			this->elim_order_type = GD_MCSM;
		if (strcmp(args[i], "-pktsort") == 0)
			this->elim_order_type = GD_PKT_SORT;
		if (strcmp(args[i], "-mind") == 0)
			this->elim_order_type = GD_MIN_DEGREE;
		if (strcmp(args[i], "-minmaxd") == 0)
			this->elim_order_type = GD_MINMAX_DEGREE;
		if (strcmp(args[i], "-mcs") == 0)
			this->elim_order_type = GD_MCS;
		if (strcmp(args[i], "-lexp") == 0)
			this->elim_order_type = GD_LEXP_BFS;
		if (strcmp(args[i], "-minf") == 0)
			this->elim_order_type = GD_MIN_FILL;
		if (strcmp(args[i], "-bmf") == 0)
			this->elim_order_type = GD_BATCH_MF;
		if (strcmp(args[i], "-mmd") == 0)
			this->elim_order_type = GD_MUL_MIN_DEGREE;
		if (strcmp(args[i], "-beta") == 0)
			this->elim_order_type = GD_BETA;
		if (strcmp(args[i], "-metmmd") == 0)
			this->elim_order_type = GD_METIS_MMD;
		if (strcmp(args[i], "-metnnd") == 0)
			this->elim_order_type = GD_METIS_NODE_ND;
		if (strcmp(args[i], "-amd") == 0)
			this->elim_order_type = GD_AMD;
		if (strcmp(args[i], "-parmetis") == 0)
		  this->parmetis = true;
		if (strcmp(args[i], "-gavril") == 0)
		  {
			this->gavril = true;
			num_td_methods++;
		  }
		if (strcmp(args[i], "-superetree") == 0)
		  {
		    this->superetree = true;
		    num_td_methods++;
		  }
		if (strcmp(args[i], "-pbag") == 0)
			this->pbag = true;
		if (strcmp(args[i], "-bk") == 0)
		  {
		    this->BK = true;
		    num_td_methods++;
		  }
		if (strcmp(args[i], "-nice") == 0)
		  {
			this->nice = true;
			num_td_methods++;
		  }
		if (strcmp(args[i], "-check") == 0)
			this->check = true;
		if (strcmp(args[i], "-nonniceDP") == 0)
			this->nonniceDP = true;
		if (strcmp(args[i], "-asp_root") == 0)
			this->asp_root = true;
		if (strcmp(args[i], "-child_root") == 0)
			this->child_root = true;
		if (strcmp(args[i], "-v") == 0)
			this->verbose = true;
		if (strcmp(args[i], "-vv") == 0)
		{
			this->very_verbose = true;
			// Assume you also want verbose output!
			this->verbose = true;
		}
		if (strcmp(args[i], "-ord") == 0)
		{
			this->read_ordering = true;
			this->ord_file = args[i + 1];
		}
		if (strcmp(args[i], "-sord") == 0)
		{
			this->read_scotch_ordering = true;
			this->scotch_ord_file = args[i + 1];
		}
		if (strcmp(args[i], "-mod") == 0)
		{
			this->write_mod = true;
			this->model_file = args[i + 1];
		}
		if (strcmp(args[i], "-root") == 0)
			this->root_node = atoi(args[i + 1]);
		if (strcmp(args[i], "-gviz") == 0)
		{
			this->gviz = true;
			this->gviz_file = args[i + 1];
		}
		if (strcmp(args[i], "-mip") == 0)
			this->solve_mip = true;
		if (strcmp(args[i], "-t") == 0)
		{
			this->tree_infile = args[i + 1];
			this->read_tree = true;
			num_td_methods++;
		}
		if (strcmp(args[i], "-w") == 0)
		{
			this->tree_outfile = args[i + 1];
			this->write_tree = true;
		}
		if (strcmp(args[i], "-tpreord") == 0)
		{
			this->write_ordered_tree = true;
		}
		if (strcmp(args[i], "-eorder") == 0)
		  this->eorder = args[i + 1];
		if (strcmp(args[i], "-width") == 0)
		  this->width = true;
		if (strcmp(args[i], "-hist") == 0)
			this->make_histogram = true;
		if (strcmp(args[i], "-pc") == 0)
			this->parent_child = true;
		if (strcmp(args[i], "-make_nice") == 0)
			this->make_nice = true;
		if (strcmp(args[i], "-refine_td") == 0)
			this->refine_td = true;
		if (strcmp(args[i], "-del_ch") == 0)
			this->free_children = true;
		if (strcmp(args[i], "-no_reconstruct") == 0)
			this->no_reconstruct = true;
		if (strcmp(args[i], "-sol") == 0)
			this->sol_file = args[i + 1];
		if (strcmp(args[i], "-split_bag") == 0)
			this->split_bag = true;
		if (strcmp(args[i], "-async_tbl") == 0)
			this->async_tbl = true;
		if (strcmp(args[i], "-dfs") == 0)
			this->DFS = true;
		if (strcmp(args[i], "-f") == 0)
		{
			this->DIMACS_file = args[i + 1];
			this->has_graph = true;
		}
		if (strcmp(args[i],"-mem_est")==0)
			this->mem_est=true;
		if (strcmp(args[i],"-max_width")==0)
			this->max_width=atoi(args[i+1]);
	}
	// Check options
	if(this->fix_DIMACS)
		// This is ok since the file is just processed and that's it
		return 0;
	if (this->nice && this->make_nice)
		fatal_error("A nice tree doesn't need to be made nice!\n");

	if (!this->sol_file)
	{
		this->sol_file = new char[100];
		sprintf(sol_file, "%s.WIS.sol", this->DIMACS_file);
	}
	// Make sure we have a graph
	if (!this->has_graph)
	{

		fatal_error("%s:  No graph loaded!\n", __FUNCTION__);
		exit(-1);
	}

	// Make sure gviz file given if gviz option set
	if (this->gviz)
		if (this->gviz_file == NULL)
			fatal_error("%s: graphviz option requires an output file name\n",
					__FUNCTION__);

	// Make sure file given if tree option set
	if (this->write_tree || this->write_ordered_tree)
		if (this->tree_outfile == NULL)
			fatal_error(
					"%s: -w option requires an output file name. -tpreord must be used with -w <filename>.\n",
					__FUNCTION__);

	// Make sure either nice gavril or bk or superetree was selected AND only one tree decomposition option
	if (num_td_methods != 1 && !this->write_mod)
	  //!this->nice && !this->gavril && !this->superetree && !this->BK && !this->write_mod && !this->read_tree && !this->pbag)
		fatal_error(
				"Have to choose exactly one of -superetree, -nice, -gavril or -BK for the TD type or read tree from file\n");


	// Make sure we have at least one method of selecting the elimination ordering
	if (this->elim_order_type == GD_UNDEFINED && !this->read_ordering
			&& !this->read_scotch_ordering && !this->read_tree
			&& !this->write_mod && !this->parmetis)
		fatal_error(
				"Must choose at least one of the elimination ordering routines\n");

	// Only one rooting alg
	if (this->asp_root == true && this->child_root == true)
		fatal_error("Only pick one rooting algorithm\n");

	return 0;

}
