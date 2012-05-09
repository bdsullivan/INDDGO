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

using namespace std;
int main(int argc, char **argv)
{

	// //0: debug, 5: critical
	// LOG_INIT("GraphInfo.log", NULL, 0);
	// DEBUG("######################################\n");
	// DEBUG("Creating DIMACS Graph Reader\n");

	// string graphType("DIMACS");
	// //const char *path = "1dc.128.txt";

	// try
	// {
	// 	int s;
	// 	Graph::GraphCreatorFile creator("1dc.128.txt", graphType);

	// 	Graph::Graph g = creator.create_graph();
	// 	cout << g.get_graph_type() << endl;
	// 	cout << g.get_num_edges() << endl;

	// 	vector<Graph::Node> nodeVec = g.get_nodes();
	// 	s = g.get_num_nodes();
	// 	cout << "Num nodes " << s << endl;
	// 	cout << g.get_graph_type() << endl;

	// 	creator.set_file_name("x.245.txt");
	// 	Graph::WeightedGraph wg = creator.create_weighted_graph();
	// 	cout << g.get_graph_type() << endl;

	// 	s = wg.get_num_nodes();
	// 	vector<int> weights = wg.get_weight();
	// 	vector<Graph::Node> nodes = wg.get_nodes();

	// 	DEBUG("Creating mutable graph\n");
	// 	creator.set_file_name("house.txt");
	// 	Graph::MutableGraph mg = creator.create_mutable_graph();
	// 	cout << mg.get_num_edges() << endl;

	// 	mg.add_edge_advance(1,1);
	// 	cout << mg.get_num_edges() << endl;

	// 	cout << mg.is_simple() << endl;
	// 	Graph::GraphProperties gp;
	// 	cout << mg.is_simple() << endl;
	// 	gp.make_simple(&mg);
	// 	gp.make_canonical(&mg);
	// 	cout << gp.check_simple(&mg) << endl;


	// } catch (Graph::GraphException &ge)
	// {
	// 	cout << ge.what() << endl;
	// }

	// LOG_CLOSE();
	return 0;
}
