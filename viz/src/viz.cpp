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
#include "viz.h" 

/** 
* Constructor for TD_info. Sets all options to false and strings to NULL.
*/
TD_info::TD_info()
{
  /*IO options*/
  this->DIMACS_file = this->tree_infile = this->tree_outfile = NULL;
  this->gviz_outfile  = (char *)"tree.dot";
  this->has_graph = false;
  this->read_tree = false;
  this->write_tree = false;
  this->make_histogram = false;  
  this->write_ordered_tree = false;

  /*TD construction options*/
  this->td_alg = GD_UNDEFINED;
  this->make_nice = false;
  
  /*EO options*/
  this->read_ordering = false;
  this->read_scotch_ordering = false;
  this->ord_file = this->scotch_ord_file = NULL; 
  this->start_v = GD_UNDEFINED;
  this->elim_order_type = GD_UNDEFINED;
  
  /*visualization options*/
  this->use_scores = false; 
  this->score_infile = NULL; 
  this->viz_style = GV_TREE_ONLY; //default black & white tree.
  this->subtree = false; 
  this->highlight_node = GD_UNDEFINED;
  this->log_range = false; 
};

/**
* TD_info destructor.
*/
TD_info::~TD_info()
{
  //nothing to see here.
};

/**
* Function to process the command line args and populate the
* TD_info structure.
*/
void TD_info::process_TD_info(int num_args, char **args,
		void(*usage)(const char *str))
{
	if (num_args == 1 || (num_args == 2 && strcmp(args[1], "-h") == 0)
			|| (num_args == 2 && strcmp(args[1], "--help") == 0) || (num_args
			== 2 && strcmp(args[1], "--h") == 0))
	{
		usage(args[0]);
		exit(-1);
	}

	for (int i = 0; i < num_args; i++)
	{
		if (strcmp(args[i], "-s") == 0)
			this->start_v = atoi(args[i + 1]);
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
		if(strcmp(args[i], "-superetree") == 0)
		  this->td_alg = TD_SUPERETREE; 
		if (strcmp(args[i], "-gavril") == 0)
		  this->td_alg = TD_GAVRIL;
		if (strcmp(args[i], "-bk") == 0)
			this->td_alg = TD_BK;
		if (strcmp(args[i], "-nice") == 0)
			this->td_alg = TD_NICE;
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
		if (strcmp(args[i], "-t") == 0)
		{
			this->tree_infile = args[i + 1];
			this->read_tree = true;
		}
		if (strcmp(args[i], "-w") == 0)
		{
			this->tree_outfile = args[i + 1];
			this->write_tree = true;
		}
		if (strcmp(args[i], "-scores") == 0)
		{
			this->score_infile = args[i + 1];
			this->use_scores = true;
		}

		if (strcmp(args[i], "-log_range") == 0)
		  this->log_range = true;

		if (strcmp(args[i], "-subtree") == 0)
		{
		  this->highlight_node = atoi(args[i + 1]);
		  this->subtree = true;
		}
		if (strcmp(args[i], "-gviz_file") == 0)
		{
		  this->gviz_outfile = args[i + 1];
		    
		}
		if (strcmp(args[i], "-hist") == 0)
			this->make_histogram = true;
		if (strcmp(args[i], "-make_nice") == 0)
			this->make_nice = true;
		if (strcmp(args[i], "-f") == 0)
		{
			this->DIMACS_file = args[i + 1];
			this->has_graph = true;
		}
		if (strcmp(args[i], "-shaded") == 0)
		  this->viz_style = GV_COLORS;
		if (strcmp(args[i], "-scaled") == 0)
		  {
		    if(this->viz_style != GV_COLORS)//implies scaled, resetting removes scoring colors.
		      this->viz_style = GV_SCALE_BAGS;
		  }
		if (strcmp(args[i], "-labels") == 0)
		  this->viz_style = GV_BAG_LABELS;
		
	}
	// Check options
	if (this->td_alg == TD_NICE && this->make_nice)
		fatal_error("A nice tree doesn't need to be made nice!\n");

	// Make sure we have a graph
	if (!this->has_graph)
	{

		fatal_error("%s:  No graph loaded!\n", __FUNCTION__);
		exit(-1);
	}

	//Make sure scores file given if scores option set
	if (this->use_scores)
	 	if (this->score_infile == NULL)
		  fatal_error("%s: scores option requires an input file name\n",
			      __FUNCTION__);

	//Make sure highlight node given if subtree option set
	if (this->subtree)
	 	if (this->highlight_node < 0 || this->highlight_node == GD_UNDEFINED)
		  fatal_error("%s: subtree option requires an integer for the node to be highlighted\n",
			      __FUNCTION__);

	// Make sure file given if tree option set
	if (this->write_tree)
		if (this->tree_outfile == NULL)
			fatal_error(
					"%s: -w option requires an output file name. -tpreord must be used with -w <filename>.\n",
					__FUNCTION__);

	//check gviz_outfile non null
		if (this->gviz_outfile == NULL)
			fatal_error(
					"%s: -gviz_file option requires an output file name.\n",
					__FUNCTION__);


	// Make sure either nice gavril or bk was selected
	if (this->td_alg == GD_UNDEFINED && !this->read_tree)
		fatal_error(
				"Have to choose either -superetree, -nice, -gavril, -BK for the TD type or read tree from file\n");

	// Make sure we have at least one method of selecting the elimination ordering
	if (this->elim_order_type == GD_UNDEFINED && !this->read_ordering
			&& !this->read_scotch_ordering && !this->read_tree)
		fatal_error(
				"Must choose at least one of the elimination ordering routines\n");

}


void usage(const char *s)
{
	fprintf(stderr, "Usage: %s [options]\n", s);
	fprintf(
			stderr,
			        "\t --- Input/Output Options ---\n"
			        "\t -f DIMACS_file : loads in the graph from file with symmetric adj lists. MANDATORY.\n"
				"\t -t <tree_infile> will read a pre-existing tree from a file (modified DIMACS)\n"
				"\t -w <tree_outfile> will write the tree to a file (modified DIMACS)\n"
				"\t -gviz_file <gviz_file>: writes a graphviz version of the tree decomposition\n"
			        "\t        to the gviz_file. If not specified, the default will be tree.dot\n"
			        "\t -scores <score_file>: reads a score vector from the specified file; scores should be" 
                                "\t separated by spaces, one per graph vertex, in the order associated with the vertex labels.\n"
			        "\t --- Visualization Specific Options ---\n" 
			        "\t -labels : labels each tree node with bag's graph vertices, colored by scores (if specified)\n"
			        "\t -shaded : labels each tree node with the graph vertices its bag contains, colored by average score (if specified), and size of bag otherwise\n"
			        "\t -scaled : scale bag size relative to number of vertices it contains. NOTE: When scores are enabled, this always implies -shaded.\n"
			        "\t -subtree <nodeid> : Highlight the subtree containing the specified vertex (given by node id). Implies -labels style.\n"
			        "\t *NOTE* If none of the above are specified, the tree will be rendered with no bag labels and no score information  - just a black & white tree structure\n"
			        "\t --- Decomposition Construction Options ---\n" 
				"\t -superetree : constructs a non-nice TD using supernodal elimination trees and CHOLMOD\n"
				"\t -gavril : constructs a non-nice TD using Gavril's algorithm\n"
				"\t -bk : constructs a non-nice TD using BK algorithm\n"
				"\t -nice : constructs a nice TD\n"
 			        "\t -make_nice will take a non-nice tree and niceify it by adding new\n"
				"\t        tree nodes\n"
			        "\t --- Elimination Order Options ---\n"
				"\t -s start_v : uses start_v as first vertex in ordering\n"
				"\t -mind : generates an elim. ordering using min degree heuristic\n"
				"\t -mmd: generates an elim. ordering using multiple min degree heuristic\n"
				"\t -minf : generates an elim. ordering using min fill heuristic\n"
				"\t -bmf: generates an elim. ordering using a batched min fill heuristic\n"
				"\t -beta: generates an elim. ordering using the beta heuristic\n"
				"\t -metmmd: generates an elim. ordering using METIS mmd heuristic\n"
				"\t -metnnd: generates an elim. ordering using METS node ND heuristic\n"
				"\t -mcsm : generates an elim. ordering using mcsm euristic\n"
				"\t -mcs  : generates an elim. ordering using mcs\n"
				"\t -lexm : generates an elim. ordering using lex-m bfs heuristic\n"
				"\t -pktsort : generates the elim ordering corresponding to partial\n"
				"\t            ktree generation\n"
				"\t -amd: generates the elim. ordering using UFs approximate min\n"
				"\t       degree heuristic (AMD)\n"
				"\t -minmaxd: generates the elim. ordering using the minimum maximum\n"
				"\t           degree heuristic.\n"
				"\t -ord <ordering_file> reads in an elimination ordering (1-based as in \n"
				"\t      DIMACS file)\n"
				"\t -sord <scotch ordering file> reads in an elimination ordering produced\n"
				"\t                              by Scotch\n"
			        "\t --- Miscellaneous Options ---\n"
			"\t -hist will print out a histogram of bag sizes to std out\n");
}

