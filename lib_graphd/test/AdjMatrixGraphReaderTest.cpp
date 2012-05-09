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
#include "GraphException.h"
#include <gtest/gtest.h>

using namespace std;
class AdjMatrixGraphReaderTest: public testing::Test
{
public:
	Graph::GraphCreatorFile creator;
	Graph::Graph *g;

	virtual void SetUp()
	{
		//SetUp is called before every test
		LOG_INIT("test.log", NULL, 0);
		try
		{
			creator.set_file_name("../data/1dc.128.adj");
			creator.set_graph_type("AdjMatrix");
			g = creator.create_graph();

		} catch (Graph::GraphException &e)
		{
			cout << e.what() << endl;
		}
	}

	virtual void TearDown()
	{
		//TearDown is called before every test
		LOG_CLOSE();
	}
};

TEST_F(AdjMatrixGraphReaderTest, testNumNodes)
{
	EXPECT_EQ(128, g->get_num_nodes());
}

TEST_F(AdjMatrixGraphReaderTest,testNumEdges)
{
	EXPECT_EQ(1471, g->get_num_edges());
}

TEST_F(AdjMatrixGraphReaderTest,testGraphType)
{
	string graph_type = g->get_graph_type();
    EXPECT_STRCASEEQ("adjmatrix", graph_type.c_str());
}

TEST_F(AdjMatrixGraphReaderTest, testGetNode)
{
    Graph::Node *n;
    n = g->get_node(108);
    list<int> nbrs = n->get_nbrs();
    EXPECT_EQ(24, nbrs.size());
}

TEST_F(AdjMatrixGraphReaderTest, testIsEdge)
{
    Graph::Node *n;
    n = g->get_node(108);
    list<int> nbrs = n->get_nbrs();

    EXPECT_EQ(24, nbrs.size());
    EXPECT_TRUE(g->is_edge(108, 100));

}

TEST_F(AdjMatrixGraphReaderTest, testGetDegree)
{
    EXPECT_EQ(19, g->get_degree(28));
}


TEST_F(AdjMatrixGraphReaderTest, testNodeNbrs)
{
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


TEST_F(AdjMatrixGraphReaderTest,testGetDegrees)
{
    vector<int> degree = g->get_degree();
    EXPECT_EQ(31, degree[86]);
}

TEST_F(AdjMatrixGraphReaderTest, testGetNumComponents)
{
    EXPECT_EQ(1, g->get_num_components());
}

TEST_F(AdjMatrixGraphReaderTest, testGetNextLabel)
{
    EXPECT_EQ(129, g->get_next_label());
}

TEST_F(AdjMatrixGraphReaderTest, testSetNextLabel)
{
    EXPECT_EQ(129, g->get_next_label());
    g->set_next_label(140);
    EXPECT_EQ(140, g->get_next_label());
}

