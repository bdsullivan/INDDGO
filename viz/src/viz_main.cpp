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
#include "Log.h"
#include "GraphException.h"
#include <fstream>
#include <string>
#include <unistd.h>

#include "WeightedMutableGraph.h"
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
      //Aaron - add the processing of any options here - should include a -color flag which requires an input file name
      //as well as any options you want in terms of the type of graphviz output you want to be able to create - you'll 
      //probably consider putting things like -avg_scaled or -label_colors or the like to determine the output. 
      info.process_TD_info(argc, argv, usage);

      print_message(0, "Options processed\n");

      // Create the graph for WIS
      create_TDviz_graph(&info, G);

      print_message(0, "Graph structure populated\n");

      // Create the tree decomposition using the options
      create_tree_decomposition(&info, G, &T);

      print_message(0, "Tree Decomposition created\n");

      //fill the bag vectors
      T->fill_bag_vecs();

      //Aaron - here you'll want to make a call to a function that reads and populates your color/score vector.
      //I would place such a function in viz.cpp/viz.h. You'll want to pass in the info struct so it can get out 
      //the color_file name. 
      if(info.use_scores){
	vector<double> color_vector;
	double max_color=0;
	double min_color=0;
	cout<<info.score_infile<<"\n";
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


