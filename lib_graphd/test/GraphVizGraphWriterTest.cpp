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

#include "Log.h"
#include "GraphCreatorFile.h"
#include "Graph.h"
#include "GraphVizGraphWriter.h"
#include <gtest/gtest.h>

using namespace std;
class GraphVizGraphWriterTest: public testing::Test {
public:
	Graph::GraphCreatorFile creator;
    Graph::Graph *g;
    Graph::GraphWriter writer;

    virtual void SetUp(){
 		//SetUp is called before every test
        LOG_INIT("test.log", NULL, 0);
        creator.set_file_name("data/1dc.128.txt");
        creator.set_graph_type("DIMACS");
		g = creator.create_graph();
    }

    virtual void TearDown(){
		//TearDown is called before every test
        LOG_CLOSE();
    }
};



TEST_F(GraphVizGraphWriterTest, testWriteGraph)
{
	string out_file("data/1dc.128.gviz");
	writer.write_graph(g, out_file, string("GRAPHVIZ"), false);
}


// Graph::Node *n;
// n = g->get_node(108);
// list<int> nbrs = n->get_nbrs();
// vector<int> nbrsvec(nbrs.begin(), nbrs.end());
// for (int i = 0; i < nbrsvec.size(); i++)
//     cout << nbrsvec[i] << 
