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
#include "Node.h"
#include "GraphReaderWriterFactory.h"
#include "GraphCreatorFile.h"
#include <gtest/gtest.h>

using namespace std;
// A new one of these is created for each test
class GraphReaderWriterFactoryTest: public testing::Test
{

public:
	Graph::GraphReaderWriterFactory factory;
	Graph::GraphCreatorFile creator;
	Graph::GraphReader *gr;
	Graph::GraphWriter *gw;
	Graph::Graph *g;

	virtual void SetUp()
	{
		LOG_INIT("test.log", NULL, 0);
		gr = factory.create_reader("DIMACS");
		gr->read_graph("../data/1dc.128.txt");

		creator.set_file_name("../data/1dc.128.txt");
		creator.set_graph_type("DIMACS");
		g = creator.create_graph();
	}

	virtual void TearDown()
	{
		LOG_CLOSE();
	}
};

TEST_F(GraphReaderWriterFactoryTest, testCapacity)
{
	EXPECT_EQ(128, gr->get_capacity());
}

TEST_F(GraphReaderWriterFactoryTest, testNumEdges)
{
	EXPECT_EQ(1471, gr->get_num_edges());
}

TEST_F(GraphReaderWriterFactoryTest, testNodeNbrs)
{
	vector<Graph::Node> n;
	int i;
	n = gr->get_nodes();

	EXPECT_EQ(128, n.size());

	list<int> nbrlist = n[35].get_nbrs();
	vector<int> nbrs(nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(75, nbrs[20]);

	nbrlist = n[127].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(111, nbrs[2]);

	nbrlist = n[0].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(1, nbrs[0]);
	EXPECT_EQ(64, nbrs[6l]);
	EXPECT_EQ(8, nbrs[3l]);
}

TEST_F(GraphReaderWriterFactoryTest, testAdjMatrixGraphWriter)
{
	gw = factory.create_writer("adjmatrix");
	gw->set_out_file_name("../data/adjmatrix.out");
	gw->write_graph(g);

	creator.set_file_name("../data/adjmatrix.out");
	creator.set_graph_type("adjmatrix");
	g = creator.create_graph();

	EXPECT_EQ(1471, g->get_num_edges());
}

TEST_F(GraphReaderWriterFactoryTest, testAdjMatrixGraphWriterNodeNbrs)
{
	gw = factory.create_writer("adjmatrix");
	gw->set_out_file_name("../data/adjmatrix.out");
	gw->write_graph(g);

	creator.set_file_name("../data/adjmatrix.out");
	creator.set_graph_type("adjmatrix");
	g = creator.create_graph();

	vector<Graph::Node> n;
	int i;
	n = g->get_nodes();

	EXPECT_EQ(128, n.size());

	list<int> nbrlist = n[35].get_nbrs();
	vector<int> nbrs(nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(75, nbrs[20]);

	nbrlist = n[127].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(111, nbrs[2]);

	nbrlist = n[0].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());

	EXPECT_EQ(1, nbrs[0]);
	EXPECT_EQ(64, nbrs[6l]);
	EXPECT_EQ(8, nbrs[3l]);
}

TEST_F(GraphReaderWriterFactoryTest, testDimacsGraphWriterGetNumEdges)
{
	gw = factory.create_writer("dimacs");
	gw->set_out_file_name("../data/dimacsout.out");
	gw->write_graph(g);

	creator.set_file_name("../data/dimacsout.out");
	creator.set_graph_type("dimacs");
	g = creator.create_graph();

	EXPECT_EQ(1471, g->get_num_edges());
}

TEST_F(GraphReaderWriterFactoryTest, testDimacsGraphWriterGetNode)
{
	gw = factory.create_writer("dimacs", "../data/dimacsout.out");
	gw->write_graph(g);

	creator.set_file_name("../data/dimacsout.out");
	creator.set_graph_type("dimacs");
	g = creator.create_graph();
	Graph::Node *n;
	n = g->get_node(108);
	list<int> nbrs = n->get_nbrs();
	EXPECT_EQ(24, nbrs.size());
}

TEST_F(GraphReaderWriterFactoryTest, testGraphVizWriter)
{
	gw = factory.create_writer("graphviz");
	gw->set_out_file_name("../data/graphviz.out");
	gw->write_graph(g);
	EXPECT_TRUE(true);
}

