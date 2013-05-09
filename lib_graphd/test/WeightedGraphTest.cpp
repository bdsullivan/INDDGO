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
#include "WeightedGraph.h"
#include <gtest/gtest.h>

using namespace std;
class WeightedGraphTest: public testing::Test
{
public:
	Graph::GraphCreatorFile creator;
	Graph::WeightedGraph *wg;

	virtual void SetUp()
	{
		//SetUp is called before every test
		LOG_INIT("test.log", NULL, 0);
		creator.set_file_name("data/x.245.txt");
		creator.set_graph_type("DIMACS");
		wg = creator.create_weighted_graph();
	}

	virtual void TearDown()
	{
		//TearDown is called before every test
		LOG_CLOSE();
	}
};

TEST_F(WeightedGraphTest, testNumNodes)
{
	EXPECT_EQ(23, wg->get_num_nodes());
}

TEST_F(WeightedGraphTest, testNumEdges)
{
	EXPECT_EQ(71, wg->get_num_edges());
}

TEST_F(WeightedGraphTest, testGraphType)
{
	EXPECT_EQ("DIMACS", wg->get_graph_type());
}

TEST_F(WeightedGraphTest, testGetNode)
{
	Graph::Node *n;
	n = wg->get_node(19);
	list<int> nbrs = n->get_nbrs();
	EXPECT_EQ(4, nbrs.size());
}

TEST_F(WeightedGraphTest, testIsEdge)
{
	Graph::Node *n;
	n = wg->get_node(19);
	list<int> nbrs = n->get_nbrs();

	EXPECT_EQ(4, nbrs.size());
	EXPECT_TRUE(wg->is_edge(19, 10));
	EXPECT_TRUE(wg->is_edge(19, 4));
	EXPECT_FALSE(wg->is_edge(19, 1));

}

TEST_F(WeightedGraphTest, testGetDegree)
{
	EXPECT_EQ(6, wg->get_degree(6));
}

TEST_F(WeightedGraphTest, testNodeNbrs)
{
	vector<Graph::Node> n;
	int i;
	n = wg->get_nodes();

	EXPECT_EQ(23, n.size());

	list<int> nbrlist = n[10].get_nbrs();
	vector<int> nbrs(nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(7, nbrs[2]);

	nbrlist = n[22].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(13, nbrs[2]);

	nbrlist = n[0].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());

	EXPECT_EQ(1, nbrs[0]);
	EXPECT_EQ(17, nbrs[6]);
	EXPECT_EQ(8, nbrs[3]);
}

TEST_F(WeightedGraphTest, testGetDegrees)
{
	vector<int> degree = wg->get_degree();
	EXPECT_EQ(4, degree[19]);
}

TEST_F(WeightedGraphTest, testGetNumComponents)
{
	EXPECT_EQ(1, wg->get_num_components());
}

TEST_F(WeightedGraphTest, testGetNextLabel)
{
	EXPECT_EQ(24, wg->get_next_label());
}

TEST_F(WeightedGraphTest, testSetNextLabel)
{
	EXPECT_EQ(24, wg->get_next_label());
	wg->set_next_label(140);
	EXPECT_EQ(140, wg->get_next_label());
}

TEST_F(WeightedGraphTest, testGetWeight)
{
	vector<int> weights = wg->get_weight();

	EXPECT_EQ(8, weights[0]);
	EXPECT_EQ(8, weights[22]);
	EXPECT_EQ(9, weights[14]);
	EXPECT_EQ(3, weights[5]);

}

TEST_F(WeightedGraphTest, testSetWeight)
{
	vector<int> weights = wg->get_weight();

	EXPECT_EQ(8, weights[0]);
	EXPECT_EQ(8, weights[22]);
	EXPECT_EQ(9, weights[14]);
	EXPECT_EQ(3, weights[5]);

	weights[5] = 100;
	weights[10] = 200;

	wg->set_weight(weights);

	weights = wg->get_weight();

	EXPECT_EQ(8, weights[0]);
	EXPECT_EQ(8, weights[22]);
	EXPECT_EQ(9, weights[14]);
	EXPECT_EQ(100, weights[5]);
	EXPECT_EQ(200, weights[10]);

}
