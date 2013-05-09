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
#include "DIMACSGraphWriter.h"
#include <string>
#include <gtest/gtest.h>

using namespace std;

class DIMACSGraphWriterTest: public testing::Test
{
public:
	Graph::GraphCreatorFile creator;
	Graph::Graph *g;
	Graph::DIMACSGraphWriter *dim_writer;

	virtual void SetUp()
	{
		//SetUp is called before every test
		LOG_INIT("test.log", NULL, 0);
		creator.set_file_name("data/1dc.128.adj");
		creator.set_graph_type("adjmatrix");
		g = creator.create_graph();

		string out_file("data/1dc.128.out");
		dim_writer = new Graph::DIMACSGraphWriter(out_file);
	}

	virtual void TearDown()
	{
		//TearDown is called before every test
		LOG_CLOSE();
	}
};

TEST_F(DIMACSGraphWriterTest, testNumNodes)
{
	dim_writer->write_graph(g);

	creator.set_file_name("data/1dc.128.out");
	creator.set_graph_type("DIMACS");
	g = creator.create_graph();

	EXPECT_EQ(128, g->get_num_nodes());
}

TEST_F(DIMACSGraphWriterTest, testGetDegrees)
{
	dim_writer->write_graph(g);

	creator.set_file_name("data/1dc.128.out");
	creator.set_graph_type("DIMACS");
	g = creator.create_graph();

	vector<int> degree = g->get_degree();
	EXPECT_EQ(31, degree[86]);
}

TEST_F(DIMACSGraphWriterTest, testGetNode)
{
	dim_writer->write_graph(g);

	creator.set_file_name("data/1dc.128.out");
	creator.set_graph_type("DIMACS");
	g = creator.create_graph();

	Graph::Node *n;
	n = g->get_node(108);
	list<int> nbrs = n->get_nbrs();
	EXPECT_EQ(24, nbrs.size());
}

TEST_F(DIMACSGraphWriterTest, testIsEdge)
{
	dim_writer->write_graph(g);

	creator.set_file_name("data/1dc.128.out");
	creator.set_graph_type("DIMACS");
	g = creator.create_graph();

	Graph::Node *n;
	n = g->get_node(108);
	list<int> nbrs = n->get_nbrs();

	EXPECT_EQ(24, nbrs.size());
	EXPECT_TRUE(g->is_edge(108, 100));

}

TEST_F(DIMACSGraphWriterTest, testshuffle)
{
    string out_file("data/1dc.128.out.shuf");
    dim_writer->set_out_file_name(out_file);
    dim_writer->shuffle(g, time(NULL));
	dim_writer->write_graph(g);

	creator.set_file_name("data/1dc.128.out.shuf");
	creator.set_graph_type("DIMACS");
	g = creator.create_graph();

	Graph::Node *n;
	n = g->get_node(108);
	list<int> nbrs = n->get_nbrs();
	//BDS - problem if randomization happens to put a vertex of 
	//degree 24 in this location. Temporary fix.
	//EXPECT_NE(24, nbrs.size());
	EXPECT_NE(0, nbrs.size());
}

// Graph::Node *n;
// n = g->get_node(108);
// list<int> nbrs = n->get_nbrs();
// vector<int> nbrsvec(nbrs.begin(), nbrs.end());
// for (int i = 0; i < nbrsvec.size(); i++)
//     cout << nbrsvec[i] << 
