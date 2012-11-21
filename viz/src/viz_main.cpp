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
#include <fstream>
#include <unistd.h>
#include "viz.h"

int main(int argc, char **argv)
{
  int i;
  Graph::WeightedMutableGraph *G;
  TDTree *T = NULL;
  TD_info info; //note that there is no TD_info in T, so you need to use this one!

  /*initialize log files if needed*/
  	int pid;
#if WIN32 || _WIN32
	pid=_getpid();
#else
	pid=(int)getpid();
#endif
	char lfname[100];
	char efname[100];
	sprintf(lfname, "viz-%d.log",pid);
	sprintf(efname, "viz-%d.log",pid);
	//0: debug, 5: critical
	LOG_INIT(lfname, efname, 0);
  
  try
    {
      // Process command line options      
      info.process_TD_info(argc, argv, usage);
      print_message(0, "Options processed\n");

      // Create the graph for WIS
      Graph::create_largestcomponent_graph(info.DIMACS_file, G);
      print_message(0, "Graph structure populated\n");

      // Create the tree decomposition using the options
      if(info.read_scotch_ordering)
      	create_tree_decomposition(G, &T, info.read_tree, info.tree_infile, info.read_ordering, info.read_scotch_ordering, 
      				  info.scotch_ord_file, info.elim_order_type, info.start_v, info.td_alg, info.make_nice, 
      				  true);      
      else
      	create_tree_decomposition(G, &T, info.read_tree, info.tree_infile, info.read_ordering, info.read_scotch_ordering, 
      				  info.ord_file, info.elim_order_type, info.start_v, info.td_alg, info.make_nice, 
      				  true);      

      // Write out the tree if desired
      if (info.write_tree)
	T->write_DIMACS_file(info.tree_outfile);
      
      // Print out a histogram if desired
      if (info.make_histogram)
	td_size_histogram(T, stdout);  

      print_message(0, "Tree Decomposition created\n");

      //fill the bag vectors
      T->fill_bag_vecs();

      if(info.use_scores)
	{
	  vector<double> color_vector;
	  double max_color=0;
	  double min_color=0;
	  cout<<info.score_infile << "\n";
	  bool range = read_color_file(info.score_infile,max_color,min_color,color_vector);

	if(range)
	  {
	    cout<<"Max color: "<<max_color<<"\n";
	    cout<<"Min color: "<<min_color<<"\n";
	    cout<<"Size of color: "<<color_vector.size()<<"\n";
	    if(info.subtree)
	      T->highlight_subtree_scored_graphviz_file(info.highlight_node, info.gviz_outfile, color_vector, max_color, min_color);
	    else
	      {

		if(info.viz_style == GV_BAG_LABELS || info.viz_style==GV_COLORS)
		  {
		    if(info.log_range)
		      T->write_scored_graphviz_file(false,info.gviz_outfile,color_vector,max_color,min_color, info.viz_style, true);
		    else
		      T->write_scored_graphviz_file(false,info.gviz_outfile,color_vector,max_color,min_color, info.viz_style);
		  }
		else
		  {
		    T->write_graphviz_file(false, info.gviz_outfile, GV_TREE_ONLY);
		  }
	      }
	  }
	else
	  {
	    cout<<"No Range provided in file.  Will be computed\n";
	    cout<<"Size of color: "<<color_vector.size()<<"\n";
	    if(info.subtree)
	      T->highlight_subtree_scored_graphviz_file(info.highlight_node, info.gviz_outfile, color_vector);
	    else
	      {
		if(info.viz_style == GV_BAG_LABELS || info.viz_style==GV_COLORS)
		  {
		    if(info.log_range)
		      T->write_scored_graphviz_file(false,info.gviz_outfile,color_vector, info.viz_style, true);
		    else
		      T->write_scored_graphviz_file(false,info.gviz_outfile,color_vector, info.viz_style);
		  }		  
		else
		  {
		    T->write_graphviz_file(false, info.gviz_outfile, GV_TREE_ONLY);
		  }
	      }
	  }	
      }
      else{
	if(info.viz_style == GV_BAG_LABELS)
	  T->write_graphviz_file(false, info.gviz_outfile, GV_BAG_LABELS);
	else
	  T->write_graphviz_file(false, info.gviz_outfile, GV_TREE_ONLY);
	cout<<"No scores were provided\n";
      }


      printf("%s: Treewidth %d\n", info.DIMACS_file, T->width);
      delete G;
      delete T;
      LOG_CLOSE();
      return 1;
	
    }
  catch(Graph::GraphException& e)
    {
      cerr << "exception caught: " << e.what() << endl;
      	LOG_CLOSE();
      return -1;
    }
        

}


