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

/*
 * Helper function for computing a tree decomposition. Valid algorithms include TD_SUPERETREE, TD_GAVRIL, TD_BK, TD_NICE.
 */
void form_td(int td_alg, TDTree**T, vector<int>* ordering)
{
  switch(td_alg)
    {
    case TD_SUPERETREE: 
      (*T)->construct_superetree(ordering);
      break;
    case TD_GAVRIL: 
      (*T)->construct_gavril(ordering);
      break;
    case TD_BK: 
      (*T)->construct_BK(ordering);
      break; 
    case TD_NICE: 
      (*T)->construct_knice(ordering, (*T)->width, false);
      break;
    default: 
      throw(Graph::GraphException("Unrecognized td construction algorithm.\n"));
    }
  /*root this at 0 to make sure writing works correctly. Can always re-root later*/
  (*T)->root(0);
}

void td_size_histogram(TDTree *T, FILE *stream)
{
  int i;
  vector<int> counts(T->width + 2, 0);
  for (i = 0; i < (int) T->tree_nodes.size(); i++)
    {
      if (T->tree_nodes[i] != NULL)
	counts[T->tree_nodes[i]->bag.size()]++;
    }
  fprintf(stream, "Histogram of bag sizes\n");
  for (i = 1; i <= T->width + 1; i++)
    {
      fprintf(stream, "%02d ", counts[i]);
    }
  fprintf(stream, "\n");
  fflush(stream);
}

/**
* Creates a tree decomposition stored in T using graph G.
* To read a td from file, set read_tree to true and specify tree_infile. 
* To use an elimination-ordering based routine for generating the TD: 
*    1) If EO is stored in file, set read_ordering and specify ord_file. Set scotch to true if in SCOTCH format.
*    2) Otherwise, specify elim_order_type as one of GD_XXX from GraphDecomposition.h. Specify start_v or 
*       set to -1/GD_UNDEFINED to use default. 
*    3) Specify a td generation algorithm in td_alg as one of TD_XXX from TreeDecomposition.h
* Set make_nice to true if you want a nice tree decomposition and chose an algorithm such as Gavril initially.
* Set timings to true to print elapsed time (in seconds) for triangulation and tree decomposition construction.
* Note, no timings will be produced if read_tree is used.   
*/
void create_tree_decomposition(Graph::WeightedMutableGraph *G, TDTree **T, bool read_tree, char *tree_infile, bool read_ordering, bool scotch, char *ord_file, int elim_order_type, int start_v, int td_alg, bool make_nice, bool timings)
{
  Graph::GraphEOUtil eoutil;
  if (!G)
    throw(Graph::GraphException("Null graph passed to create_tree_decomposition\n"));
  
  // Create a copy of G to do the triangulation and elim order 
  Graph::WeightedMutableGraph H = *G;
  (*T) = new TDTree(&H);
  // Set the name of T's graph file
  char *cstr = new char[G->get_input_file().length()+1];
  strcpy(cstr, (G->get_input_file()).c_str());
  (*T)->graph_file = cstr;
  int i;
  
  // create a vector for the elimination ordering
  vector<int> ordering(H.get_num_nodes(), GD_UNDEFINED);
  
  if (read_tree)
    {
      // Read in the decomposition from a file
      (*T)->width = (*T)->read_DIMACS_file(tree_infile);
      // Assume not nice - could check this in read_DIMACS_file actually...
      (*T)->nice = false;
      /*root this at 0 to make sure writing works correctly. Can always re-root later*/
      (*T)->root(0);
      print_message(0,"Read tree.\n");
    }
  else
    {
	Graph::form_eo(read_ordering, scotch, ord_file, elim_order_type, start_v, &H, &ordering);
	
	// Triangulate the graph
	clock_t tri_start = clock(), tri_stop;
	if(td_alg == TD_SUPERETREE)
	  {
	    //no need to triangulate
	    tri_stop = tri_start;	
	  }
	else
	  {
#if HAS_METIS
	    (*T)->width = eoutil.METIS_triangulate(&H, &ordering);
#else
	    (*T)->width = eoutil.triangulate(&H, &ordering);
#endif
	    tri_stop = clock();
	  }
	if(timings)
	  {
	    print_message(0, "Triangulation took %f secs\n",
			  (double) (tri_stop - tri_start) / CLOCKS_PER_SEC);
	  }
	
	// Now create the tree
	clock_t td_start = clock(), td_stop;
	form_td(td_alg, T, &ordering);
	td_stop=clock();
	if (timings)
	  {
	    print_message(0, "Tree construction took %f secs\n",
			  (double) (td_stop - td_start) / CLOCKS_PER_SEC);
	  }  
    }//end of creation from ordering

  // Sort the bags
  if(td_alg != TD_SUPERETREE)//bags already sorted
      (*T)->sort_bags();

  // Reset (*T)'s graph is the original, non-triangulated graph!
  (*T)->G = G;
  
  // Make the tree nice now if requested
  if (make_nice)
    (*T)->make_nice();
  
}

/*
 * Calculates a statistic for each bag in the tree decomposition and stores
 * to the stats vector so that stats[i] = statistic(bag assoc w/ node i)
 * stat_flag should be one of GD_XXX defined in Util.h
 * Note this will fail if T's bag vectors have not been filled (e.g. with fill_bag_vecs()).
 */
void bag_statistics(TDTree *T, const vector<double> &scores, vector<double> *stats, int stat_flag)
{
  int size = T->tree_nodes.size();
  TDTreeNode *curr;
  int curr_node = 0;
  stats->resize(T->num_tree_nodes);
  for(int i; i < size; i++)
    {
      curr = T->tree_nodes[i];
      if(curr != NULL)
	{
	  (*stats)[curr_node] = get_statistics(curr->bag_vec, scores, stat_flag);
	  curr_node++;
	}
      //no else clause; we don't record stats for 'missing' nodes.
    }

  
}

/*
 * Calculates the "length" of each bag (maximum shortest path distance in the 
 * original graph between pairs of vertices in that bag).
 * Note this will fail if T's bag vectors have not been filled (e.g. with fill_bag_vecs()).
 */
void bag_lengths(TDTree *T, vector<int> *lengths)
{
  int size = T->tree_nodes.size();
  TDTreeNode *curr;
  Graph::GraphUtil util;
  int curr_node = 0;
  lengths->resize(T->num_tree_nodes);
  for(int i; i < size; i++)
    {
      curr = T->tree_nodes[i];
      if(curr != NULL)
	{
	  (*lengths)[curr_node] = util.subset_max_dist(T->G, curr->bag_vec);
	  curr_node++;
	}
      //no else clause; we don't record stats for 'missing' nodes.
    }
}
