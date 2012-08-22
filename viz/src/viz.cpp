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
#include "DIMACSGraphWriter.h"

/** 
* Constructor for TD_info. Sets all options to false and strings to NULL.
*/
TD_info::TD_info()
{
  this->DIMACS_file = this->ord_file =
    this->score_infile 
    = this->scotch_ord_file = this->tree_infile = this->tree_outfile = NULL;
  this->gviz_outfile  = (char *)"tree.dot";
  this->has_graph = false;
  this->read_ordering = false;
  this->read_scotch_ordering = false;
    this->nice = false;
    this->superetree = false;
  this->gavril = false;
  this->BK = false;
  this->read_tree = false;
  this->write_tree = false;
  this->make_histogram = false;
  this->make_nice = false;
  this->write_ordered_tree = false;
  this->start_v = GD_UNDEFINED;
  this->elim_order_type = GD_UNDEFINED;
  this->use_scores = false; //Aaron - by default, this is turned off.
  this->viz_style = GV_TREE_ONLY; //default black & white tree.
  this->subtree = false; 
  this->highlight_node = GD_UNDEFINED;
  this->log_range = false; //Blair, this is turned off by default
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
		  this->superetree = true;
		if (strcmp(args[i], "-gavril") == 0)
		  this->gavril = true;
		if (strcmp(args[i], "-bk") == 0)
			this->BK = true;
		if (strcmp(args[i], "-nice") == 0)
			this->nice = true;
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
		  this->viz_style = GV_SCALE_BAGS;
		if (strcmp(args[i], "-labels") == 0)
		  this->viz_style = GV_BAG_LABELS;
		
	}
	// Check options
	if (this->nice && this->make_nice)
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
	if (!this->nice && !this->gavril && !this->BK && !this->superetree && !this->read_tree)
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


/**
* Populates the provided Graph structure with the necessary information to do TD and visualizations.
*/
void create_TDviz_graph(TD_info *info, Graph::WeightedMutableGraph *&G)
{
	int i;

	Graph::GraphCreatorFile creator;
	creator.set_file_name(info->DIMACS_file);
	creator.set_graph_type("DIMACS");
	G = creator.create_weighted_mutable_graph();

    if (!G)
    {
        return;
    }

    Graph::GraphProperties properties;
    Graph::GraphEOUtil eoutil;
    Graph::GraphWriter *writer;
    Graph::GraphReaderWriterFactory factory;
    
    writer = factory.create_writer("DIMACS");
    // Put the input graph in canonical form for the tests
    properties.make_canonical(G);

// Make sure there is only 1 component in the graph!
	if (!properties.is_connected(G))
	{
		print_message(0,
				"Visualizations run only on connected graphs.  This graph has more than one component!\n");
		list<Graph::WeightedMutableGraph *> C;

		C = creator.create_all_components(G, true);
		print_message(0, "Found %d components\n", C.size());
		list<Graph::WeightedMutableGraph *>::iterator cc;
		int comp_number = 0;
		char temp_file[100], comp_file[100];
		int max_size = -1;
		char *big_file = NULL;
		for (cc = C.begin(); cc != C.end(); ++cc)
		{
			if ((*cc)->get_num_nodes() > 1 && (*cc)->get_num_edges() > 1)
			{
				sprintf(temp_file, "%s_%d_node_component_%d.temp",
						info->DIMACS_file, (*cc)->get_num_nodes(), comp_number);
				sprintf(comp_file, "%s_%d_node_component_%d.comp",
						info->DIMACS_file, (*cc)->get_num_nodes(), comp_number);
				writer->set_out_file_name(temp_file);
				writer->write_graph(*cc);

				normalize_DIMACS_file(temp_file, comp_file);
				remove(temp_file);
				fprintf(stderr, "Wrote component with %d nodes to file %s\n",
						(*cc)->get_num_nodes(), comp_file);
				comp_number++;
				if ((*cc)->get_num_nodes() > max_size)
				{
					if (!big_file)
						big_file = (char *) malloc(100);

					max_size = (*cc)->get_num_nodes();
					sprintf(big_file, "%s", comp_file);
				}
			}
			else
			{
				fprintf(
						stderr,
						"Tiny component with %d nodes & %d edges not written\n",
						(*cc)->get_num_nodes(), (*cc)->get_num_edges());
			}
		}

        
		if (big_file)
		  {
		    info->DIMACS_file = (char *) malloc(100);//Blair - Mar 20, 2012
		    sprintf(info->DIMACS_file, "%s", big_file);
		    free(big_file);
		    print_message(0,
				  "Will run WIS on largest component %s with %d nodes\n",
				  info->DIMACS_file, max_size);
		    delete G;
		    create_TDviz_graph(info, G);
		  }
		else
		  {
		    print_message(0, "No big file found. Exiting\n");
		    exit(-1);
		  }

        int s = C.size();
        for (i = 0; i < s; i++)
        {
            delete C.front();
            C.pop_front();
        }
	}

    delete writer;
}

/**
* Creates a tree decomposition for the provided tree and the DP options
* contained in DP_info.
*/
void create_tree_decomposition(TD_info *info, Graph::WeightedMutableGraph *G,
		TDTree **T)
{
  // Create a copy of G to do the triangulation and elim order
  Graph::GraphEOUtil eoutil;
  
  if (!G)
    {
      return;
    }
  
  Graph::WeightedMutableGraph H = *G;
  (*T) = new TDTree(&H);
  // Set the name of T's graph file
  sprintf((*T)->graph_file, "%s", info->DIMACS_file);
  
  // create a vector for the elimination ordering
  vector<int> ordering(H.get_num_nodes(), GD_UNDEFINED);
  int i;
  
  // Set the info pointer so we have access to the options
  //(*T)->info = info;
  
  //	creating an instance of eoutils class
  
  // Now read in the tree file if we have one
  if (info->read_tree)
    {
      // Read in the decomposition from a file
      (*T)->width = (*T)->read_DIMACS_file(info->tree_infile);
      // Assume not nice - could check this in read_DIMACS_file actually...
      (*T)->nice = false;
    }
  else
    {
      // Figure out how to create the tree
      if (info->read_ordering)
	{
	  read_ordering_file(info->ord_file, &ordering);
	  print_message(0, "Read in ordering\n");
	  print(1, ordering);
	}
      else
	{
	  if (info->read_scotch_ordering)
	    {
	      read_SCOTCH_ordering_file(info->scotch_ord_file, &ordering);
	      print_message(0, "Read in SCOTCH ordering\n");
	      print(1, ordering);
	    }
	  else
	    {
	      // Create the ordering via a heuristic
	      // Create an ordering - if start_v not provided, find a good candidate
	      if (info->start_v == GD_UNDEFINED)
		eoutil.find_elimination_ordering(&H, &ordering,
						 info->elim_order_type, false);
	      //					H.find_elimination_ordering(&ordering,
	      //							info->elim_order_type, false);
	      else
		eoutil.find_elimination_ordering(&H, &ordering,
						 info->elim_order_type, info->start_v, false);
	      //					H.find_elimination_ordering(&ordering,
	      //							info->elim_order_type, info->start_v, false);
	    }
	}
      

			      


      // Triangulate the graph
      clock_t tri_start = clock(), tri_stop;
      if(info->superetree){
	//no need to triangulate
	tri_stop = tri_start;	
      }
      
#if HAS_METIS
      (*T)->width = eoutil.METIS_triangulate(&H, &ordering);
#else
      (*T)->width = eoutil.triangulate(&H, &ordering);
#endif
      tri_stop = clock();

      
      // Now create the tree
      info->start=clock();
      if(info->superetree){
	(*T)->construct_superetree(&ordering);
      }
      if (info->gavril)
	{
	  (*T)->construct_gavril(&ordering);
	}
      
      if (info->BK)
	{
	  (*T)->construct_BK(&ordering);
	}
      
      if (info->nice)
	{
	  print_message(1, "nice\n");
	  (*T)->construct_knice(&ordering, (*T)->width, false);
	}
      info->stop=clock();
      if (0)
	{
	  print_message(0, "Tree construction took %f secs\n",
			(double) (info->stop - info->start) / CLOCKS_PER_SEC);
	}
    }

  
  // Sort the bags
  int num_tree_nodes=(int) (*T)->tree_nodes.size();
  for (i = 0; i < num_tree_nodes; i++)
    {
      if ((*T)->tree_nodes[i])
	(*T)->tree_nodes[i]->bag.sort();
    }
  
  // Reset (*T)'s graph is the original, non-triangulated graph!
  (*T)->G = G;
  
  // Make the tree nice now if requested
  if (info->make_nice)
    (*T)->make_nice();

  // Write out the tree if desired
  if (info->write_tree)
    (*T)->write_DIMACS_file(info->tree_outfile);
  
  // Print out a histogram if desired
  if (info->make_histogram)
    {
      vector<int> counts((*T)->width + 2, 0);
      for (i = 0; i < (int) (*T)->tree_nodes.size(); i++)
	{
	  if ((*T)->tree_nodes[i] != NULL)
	    counts[(*T)->tree_nodes[i]->bag.size()]++;
	}
      printf("Histogram of bag sizes\n");
      for (i = 1; i <= (*T)->width + 1; i++)
	{
	  printf("%02d ", counts[i]);
	}
      printf("\n");
      fflush(stdout);
    }
  
}



