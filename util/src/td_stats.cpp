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
#include <sstream> 

/*
 * This file generates one or more tree decompositions and 
 * calculates a variety of statistics (width, bag size distribution, 
 * average bag score), as well as optional treewidth lower bounds. 
 * It generates a single output file containing all requested 
 * information in a csv format. 
 * As usual, this operates on the largest connected component.
 */
int main(int argc, char **argv)
{
  /*current required arguments (in this order, no flags):
   * graph_file (dimacs)
   * kcore_file (scorefile format)
   * output_file_basename
   * tree_outfile_basename
   */
  vector<double> *kcore_score = new vector<double>();
  vector<double> *degree_score = new vector<double>();
  int tdt[] = {TD_SUPERETREE, TD_GAVRIL, TD_BK, TD_NICE}; 
  int et[] = {GD_AMD, GD_METIS_NODE_ND, GD_METIS_MMD};
  const char* tdtype[] = {"Super-E-Tree", "Gavril", "Bodlaender-Koster", "Nice"};
  const char* elimtype[] = {"AMD", "MetisNodeND", "MetisMultMinD"};

  vector<int> td_types(tdt, tdt+sizeof(tdt)/sizeof(tdt[0]));
  vector<int> elim_types(et, et+sizeof(et)/sizeof(et[0]));
  vector<double>* st[] = {kcore_score, degree_score};
  vector<vector<double> *> scores(st, st+sizeof(st)/sizeof(st[0]));

  /*File names; this needs some checking to catch NULLS/miscalls*/
  char *graph_file = argv[1];
  char *kcore_file = argv[2]; 
  /*we'll print to stdout otherwise*/
  char *out_file_base = NULL;
  if(argc > 3)
    out_file_base = argv[3];
  char *tree_file_base = NULL; 
  if(argc > 4)
    tree_file_base = argv[4];
  
  int t, e, i, s;
  Graph::WeightedMutableGraph *G;
  TDTree *T;
  int treewidth, treelength, minecc;
  double kcore_max, kcore_min; 
  int degree_max, degree_min;
  std::ofstream out;
  std::ostream outStream(cout.rdbuf());
  std::stringstream sstm;
  int td_id;

  /*initialize log files if needed*/
  int pid;
#if WIN32 || _WIN32
  pid=_getpid();
#else
  pid=(int)getpid();
#endif
  char lfname[100];
  char efname[100];
  sprintf(lfname, "stats-%d.log",pid);
  sprintf(efname, "stats-%d.log",pid);
  //0: debug, 5: critical
  LOG_INIT(lfname, efname, 0);
  
  try
    {
      if(kcore_file == NULL || graph_file == NULL)
	throw(Graph::GraphException("Call with two arguments: graph_file kcore_file\n"));
     
      /*populate the graph*/
      Graph::create_largestcomponent_graph(graph_file, G);      
     
      /*populate appropriate score vectors*/
      bool range = read_color_file(kcore_file,kcore_max,kcore_min,*kcore_score);
      Graph::GraphUtil gutil; 
      gutil.recompute_degrees(G);
      vector<int> idegree_score = G->get_degree();
      (*degree_score).resize(idegree_score.size());
      for(int i = 0; i < idegree_score.size(); i++)
	(*degree_score)[i] = idegree_score[i];

      
      /*loop over tree decomposition algorithms*/
      for(t = 0; t < td_types.size(); t++)
	{
	  /*loop over elimination order heuristics*/
	  for(e = 0; e < elim_types.size(); e++)
	    {
	      td_id = 10*(t+1) + e; //gives a unique number for the decomposition
	      
	      /*form the tree decomposition*/
	      create_tree_decomposition(G, &T, false, NULL, false, 
					false, NULL, elim_types[e], 
					GD_UNDEFINED, td_types[t], 
					false, true);
	      /*create a place to store all the statistics for this particular decomposition*/
	      vector<vector<double> > stats(T->num_tree_nodes);
	      vector<double> mystats;
	      vector<double>::iterator it;
	      
	      //fill the bag vectors
	      T->fill_bag_vecs();
	      cout << "T has " << T->num_tree_nodes << " tree nodes\n";
	      

	      //Non-score-specific statistics - width, length, eccentricity

	      /* Width - uses degree scores for non-negativity check*/
	      /* We define width = |B| - 1 to coincide with treewidth*/
	      bag_statistics(T, *(scores[0]), &mystats, GD_STAT_COUNT);
	      treewidth = 0;
	      for(int i = 0; i < mystats.size(); i++)
		{
		  if(mystats[i]-1 > treewidth)
		    treewidth = (int)(mystats[i])-1;
		  stats[i].push_back(mystats[i]-1);
		}

	      /*Length*/
	      vector<int> lengths;
	      treelength = 0;
	      bag_lengths(T, &lengths);
	      for(int i = 0; i < lengths.size(); i++)
		{
		  if(lengths[i] > treelength)
		    treelength = lengths[i];
		  stats[i].push_back((double)lengths[i]);
		}

	      /*Eccentricity*/
	      Graph::MutableGraph *GT = T->export_tree();
	      vector<int> ecc;
	      minecc = INT_MAX;
	      gutil.find_ecc(GT, &ecc);
	      for(int i = 0; i < ecc.size(); i++)
		{
		  if(ecc[i] < minecc)
		    minecc = ecc[i];
		  stats[i].push_back(ecc[i]);
		}

	      /*loop over scores and calculate mean med stddev*/
	      for(s = 0; s < scores.size(); s++)
		{
		  /*Mean*/
		  bag_statistics(T, *(scores[s]), &mystats, GD_STAT_MEAN);
		  for(int i = 0; i < mystats.size(); i++)
		    stats[i].push_back(mystats[i]);
		  
		  /*Median*/
		  bag_statistics(T, *(scores[s]), &mystats, GD_STAT_MED);
		  for(int i = 0; i < mystats.size(); i++)
		    stats[i].push_back(mystats[i]);

		  /*Standard Deviation*/
		  bag_statistics(T, *(scores[s]), &mystats, GD_STAT_STD);
		  for(int i = 0; i < mystats.size(); i++)
		    stats[i].push_back(mystats[i]);
		}
	      
	      /*
	       * Open output file for writing results
	       */
	      if(out_file_base == NULL)
		outStream.rdbuf(std::cout.rdbuf());
	      else 
		{
		  sstm << out_file_base << td_id; 
		  out.open((sstm.str()).c_str(), fstream::out);
		  if(out.fail())
		    fatal_error("%s:  Error opening file %s for writing graphviz output\n", __FUNCTION__, (sstm.str()).c_str());
		  outStream.rdbuf(out.rdbuf());
		  sstm.str("");//clears the stream
		}
	      
	      /*
	       * Print the file header
	       */
	      outStream << "# " << graph_file << endl;
	      outStream << "# " << tdtype[t] << endl;
	      outStream << "# " << elimtype[e] << endl;
	      outStream << "# " << "Width " << (treewidth-1) << endl;
	      outStream << "# " << "Length " << treelength << endl;
	      outStream << "# " << "MinEcc " << minecc << endl;
	      outStream << "# Bag\t Cardinality\t Length\t Eccentricity\t MeanDegree\t MedianDegree\t StdDevDegree\t MeanScore\t MedianScore\t StdDevScore" << endl;

	      /*
	       * Print the statistics for each bag
	       */
	      for(int i = 0; i < mystats.size(); i++)
		{
		  outStream << i << "\t";
		  for(it = stats[i].begin(); it != stats[i].end(); ++it)
		    outStream << *it << "\t";
		  outStream << endl;
		}

	      /*
	       * Close the output file.
	       */
	      out.close();

	      /*
	       * Write the tree decomposition to file, if required.
	       */
	      if(tree_file_base != NULL)
		{
		  sstm << tree_file_base << td_id; 
		  T->write_DIMACS_file((sstm.str()).c_str());
		  sstm.str("");//clears the stream
		}

	      /*delete the tree decomposition*/
	      delete T;
	      
	    }
	}

      delete G;
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


