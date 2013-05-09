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
#include "MutableGraph.h"
#include <gtest/gtest.h>

using namespace std;
class MutableGraphTest: public testing::Test
{
public:
	Graph::GraphCreatorFile creator;
	Graph::MutableGraph *mg;

	virtual void SetUp()
	{
		//SetUp is called before every test
		LOG_INIT("test.log", NULL, 0);
		creator.set_file_name("data/1dc.128.txt");
		creator.set_graph_type("DIMACS");
		mg = creator.create_mutable_graph();
	}

	virtual void TearDown()
	{
		//TearDown is called before every test
		LOG_CLOSE();
	}
};

TEST_F(MutableGraphTest, testNumNodes)
{
	EXPECT_EQ(128, mg->get_num_nodes());
}

TEST_F(MutableGraphTest, testNumEdges)
{
	EXPECT_EQ(1471, mg->get_num_edges());
}

TEST_F(MutableGraphTest, testGraphType)
{
	EXPECT_EQ("DIMACS", mg->get_graph_type());
}

TEST_F(MutableGraphTest, testGetNode)
{
	Graph::Node *n;
	n = mg->get_node(108);
	list<int> nbrs = n->get_nbrs();
	EXPECT_EQ(24, nbrs.size());
}

TEST_F(MutableGraphTest, testIsEdge)
{
	Graph::Node *n;
	n = mg->get_node(108);
	list<int> nbrs = n->get_nbrs();

	EXPECT_EQ(24, nbrs.size());
	EXPECT_TRUE(mg->is_edge(108, 100));

}

TEST_F(MutableGraphTest, testGetDegree)
{
	EXPECT_EQ(19, mg->get_degree(28));
}

TEST_F(MutableGraphTest, testNodeNbrs)
{
	vector<Graph::Node> n;
	int i;
	n = mg->get_nodes();

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

TEST_F(MutableGraphTest, testGetDegrees)
{
	vector<int> degree = mg->get_degree();
	EXPECT_EQ(31, degree[86]);
}

TEST_F(MutableGraphTest, testGetNumComponents)
{
	EXPECT_EQ(1, mg->get_num_components());
}

TEST_F(MutableGraphTest, testGetNextLabel)
{
	EXPECT_EQ(129, mg->get_next_label());
}

TEST_F(MutableGraphTest, testSetNextLabel)
{
	EXPECT_EQ(129, mg->get_next_label());
	mg->set_next_label(140);
	EXPECT_EQ(140, mg->get_next_label());
}

TEST_F(MutableGraphTest, testAddEdge)
{
	EXPECT_FALSE(mg->is_edge(80, 12));
	EXPECT_EQ(1471, mg->get_num_edges());
	EXPECT_EQ(22, mg->get_degree(80));
	EXPECT_EQ(19, mg->get_degree(12));
	mg->add_edge(80, 12);
	EXPECT_TRUE(mg->is_edge(80, 12));
	EXPECT_TRUE(mg->is_edge(12, 80));
	EXPECT_EQ(1472, mg->get_num_edges());
	EXPECT_EQ(23, mg->get_degree(80));
	EXPECT_EQ(20, mg->get_degree(12));
}

TEST_F(MutableGraphTest, testAddEdgeAdvance)
{
	EXPECT_FALSE(mg->is_edge(80, 12));
	EXPECT_EQ(1471, mg->get_num_edges());
	EXPECT_EQ(22, mg->get_degree(80));
	EXPECT_EQ(19, mg->get_degree(12));
	mg->add_edge_advance(80, 12);
	EXPECT_TRUE(mg->is_edge(80, 12));
	EXPECT_TRUE(mg->is_edge(12, 80));
	EXPECT_EQ(1472, mg->get_num_edges());
	EXPECT_EQ(23, mg->get_degree(80));
	EXPECT_EQ(20, mg->get_degree(12));
}

TEST_F(MutableGraphTest, testRemoveEdge)
{
	EXPECT_TRUE(mg->is_edge(108, 100));
	EXPECT_EQ(1471, mg->get_num_edges());
	EXPECT_EQ(24, mg->get_degree(108));
	EXPECT_EQ(24, mg->get_degree(100));

	mg->remove_edge(108, 100);

	EXPECT_FALSE(mg->is_edge(108, 100));
	EXPECT_FALSE(mg->is_edge(100, 108));
	EXPECT_EQ(1470, mg->get_num_edges());
	EXPECT_EQ(23, mg->get_degree(100));
	EXPECT_EQ(23, mg->get_degree(108));
}

TEST_F(MutableGraphTest, testRemoveVertex)
{
	EXPECT_TRUE(mg->is_edge(108, 100));
	EXPECT_EQ(1471, mg->get_num_edges());
	EXPECT_EQ(24, mg->get_degree(108));
	EXPECT_EQ(24, mg->get_degree(100));

	mg->remove_vertex(108);

	EXPECT_FALSE(mg->is_edge(108, 100));
	EXPECT_FALSE(mg->is_edge(100, 108));
	EXPECT_EQ(1447, mg->get_num_edges());
	EXPECT_EQ(23, mg->get_degree(100));
	EXPECT_EQ(0, mg->get_degree(108));
}

TEST_F(MutableGraphTest, testContractEdge)
{
	EXPECT_TRUE(mg->is_edge(108, 100));
	EXPECT_EQ(1471, mg->get_num_edges());
	EXPECT_EQ(24, mg->get_degree(108));
	EXPECT_EQ(24, mg->get_degree(100));

	Graph::Node *na;
	na = mg->get_node(108);
	list<int> nbrs_a = na->get_nbrs();
	vector<int> nbrsvec_a(nbrs_a.begin(), nbrs_a.end());

	Graph::Node *nb;
	nb = mg->get_node(100);
	list<int> nbrs_b = nb->get_nbrs();
	vector<int> nbrsvec_b(nbrs_b.begin(), nbrs_b.end());

	int common_nbrs = 0;
	int final_degree = 0;
	int new_edges = 0;

	for (int ia = 0; ia < nbrsvec_a.size(); ia++)
	{
		for (int ib = 0; ib < nbrsvec_b.size(); ib++)
		{
			if (nbrsvec_a[ia] == nbrsvec_b[ib])
			{
				common_nbrs++;
			}
		}
	}

	final_degree = mg->get_degree(108) + mg->get_degree(100) - common_nbrs - 2;
	new_edges = mg->get_num_edges() - common_nbrs - 1;

	int x = mg->contract_edge(108, 100);

	EXPECT_FALSE(mg->is_edge(108, 100));
	EXPECT_FALSE(mg->is_edge(100, 108));
	EXPECT_EQ(0, mg->get_degree(100));
	EXPECT_EQ(new_edges, mg->get_num_edges());
	EXPECT_EQ(final_degree, mg->get_degree(108));
}

TEST_F(MutableGraphTest, testEliminateVertex)
{
	EXPECT_TRUE(mg->is_edge(108, 100));
	EXPECT_EQ(1471, mg->get_num_edges());
	EXPECT_EQ(24, mg->get_degree(108));
	EXPECT_EQ(24, mg->get_degree(100));

	Graph::Node *na;
	na = mg->get_node(108);
	list<int> nbrs_a = na->get_nbrs();
	vector<int> nbrsvec_a(nbrs_a.begin(), nbrs_a.end());

	Graph::Node *nb;
	nb = mg->get_node(100);
	list<int> nbrs_b = nb->get_nbrs();
	vector<int> nbrsvec_b(nbrs_b.begin(), nbrs_b.end());

	int common_nbrs = 0;
	int final_degree = 0;
	int new_edges = 0;

	for (int ia = 0; ia < nbrsvec_a.size(); ia++)
	{
		for (int ib = 0; ib < nbrsvec_b.size(); ib++)
		{
			if (nbrsvec_a[ia] == nbrsvec_b[ib])
			{
				common_nbrs++;
			}
		}
	}

	new_edges = mg->get_degree(108) - common_nbrs - 1;
	new_edges = new_edges + mg->get_degree(100) - 1;


	mg->eliminate_vertex(108, NULL, true);

	EXPECT_FALSE(mg->is_edge(108, 100));
	EXPECT_FALSE(mg->is_edge(100, 108));
	EXPECT_EQ(1597, mg->get_num_edges());
	EXPECT_EQ(new_edges, mg->get_degree(100));
	EXPECT_EQ(0, mg->get_degree(108));
}

// Graph::Node *n;
// n = mg->get_node(108);
// list<int> nbrs = n->get_nbrs();
// vector<int> nbrsvec(nbrs.begin(), nbrs.end());
// for (int i = 0; i < nbrsvec.size(); i++)
//     cout << nbrsvec[i] << 
