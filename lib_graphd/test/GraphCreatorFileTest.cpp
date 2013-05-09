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

#include "GraphCreatorFile.h"
#include "Log.h"
#include "Node.h"
#include "WeightedMutableGraph.h"
#include "GraphUtil.h"
#include <gtest/gtest.h>

using namespace std;
// A new one of these is created for each test
class GraphCreatorFileTest: public testing::Test
{
public:
	Graph::GraphCreatorFile *creator;
	virtual void SetUp()
        {
            LOG_INIT("test.log", NULL, 0);
            creator = new Graph::GraphCreatorFile("data/1dc.128.txt", "DIMACS");
        }

	virtual void TearDown()
        {
            LOG_CLOSE();
        }
};

TEST_F(GraphCreatorFileTest, testGraph)
{
	Graph::Graph *g = creator->create_graph();
	EXPECT_EQ(128, g->get_num_nodes());
	EXPECT_EQ(1471, g->get_num_edges());

	vector<Graph::Node> n;
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

TEST_F(GraphCreatorFileTest, testMutableGraph)
{
	Graph::MutableGraph *mg = creator->create_mutable_graph();
	EXPECT_EQ(128, mg->get_num_nodes())
		;
	EXPECT_EQ(1471, mg->get_num_edges())
		;

	vector<Graph::Node> n;
	n = mg->get_nodes();

	EXPECT_EQ(128, n.size())
		;

	list<int> nbrlist = n[35].get_nbrs();
	vector<int> nbrs(nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(75, nbrs[20])
		;

	nbrlist = n[127].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(111, nbrs[2])
		;

	nbrlist = n[0].get_nbrs();
	nbrs.clear();
	nbrs.insert(nbrs.end(), nbrlist.begin(), nbrlist.end());
	EXPECT_EQ(1, nbrs[0])
		;
	EXPECT_EQ(64, nbrs[6l])
		;
	EXPECT_EQ(8, nbrs[3l])
		;

	EXPECT_EQ(24, mg->get_degree(35))
		;
	EXPECT_EQ(28, mg->get_degree(75))
		;
	mg->remove_edge(35, 75);
	EXPECT_EQ(23, mg->get_degree(35))
		;
	EXPECT_EQ(27, mg->get_degree(75))
		;

	mg->remove_vertex(35);
	EXPECT_EQ(0, mg->get_degree(35))
		;
	EXPECT_EQ(27, mg->get_degree(75))
		;

}

TEST_F(GraphCreatorFileTest, testGetRandomEdgeSubgraph)
{
	Graph::WeightedMutableGraph *wmg;
	Graph::WeightedMutableGraph *wmg_sub;
	creator->set_file_name("data/1dc.128.txt");
	creator->set_graph_type("DIMACS");
	wmg = creator->create_weighted_mutable_graph();
	wmg_sub = creator->create_random_edge_subgraph(wmg, 70);

	EXPECT_EQ(1471, wmg->get_num_edges())
		;
	EXPECT_EQ(1029, wmg_sub->get_num_edges())
		;
}

TEST_F(GraphCreatorFileTest, testCreateKTree)
{
	Graph::MutableGraph *mg;
	mg = creator->initialize_ktree(100, 50);
	EXPECT_EQ(100, mg->get_num_nodes())
		;
}

TEST_F(GraphCreatorFileTest, testCreateInducedSubGraph)
{
	Graph::WeightedMutableGraph *wmg;
	Graph::WeightedMutableGraph *wmg_sub;
	Graph::GraphUtil util;

	creator->set_file_name("data/1et.64.txt");
	creator->set_graph_type("DIMACS");

	wmg = creator->create_weighted_mutable_graph();

	list<int> f;
	list<int>::iterator it;
	list<int>::iterator jt;
	list<int>::iterator kt;
	util.find_component(wmg, 25, &f);
    
	wmg_sub = creator->create_induced_subgraph(wmg, &f, true);

	EXPECT_EQ(f.size(), wmg_sub->get_num_nodes());

	vector<Graph::Node> wmg_nodes = wmg->get_nodes();
	vector<Graph::Node> wmg_sub_nodes = wmg_sub->get_nodes();
	int i = 0;
	int j = 0;
	int n = wmg_sub->get_num_nodes();

	it = f.begin();
	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(wmg_nodes[*it].get_label(), wmg_sub_nodes[i].get_label())
			;
		EXPECT_EQ(wmg_nodes[*it].get_degree(), wmg_sub_nodes[i].get_degree())
			;
		list<int> wmg_nbrs = wmg_nodes[*it].get_nbrs();
		list<int> wmg_sub_nbrs = wmg_sub_nodes[i].get_nbrs();
		wmg_nbrs.sort();
		wmg_sub_nbrs.sort();

		for (jt = wmg_nbrs.begin(), kt = wmg_sub_nbrs.begin(); jt
                 != wmg_nbrs.end(); ++jt, ++kt)
		{
			ASSERT_EQ(wmg_nodes[*jt].get_label(), wmg_sub_nodes[*kt].get_label())
				;
		}
		++it;
	}

}

TEST_F(GraphCreatorFileTest, testCreateComponent)
{
	Graph::WeightedMutableGraph *wmg;
	Graph::WeightedMutableGraph *wmg_sub;
	Graph::GraphUtil util;

	creator->set_file_name("data/1et.64.txt");
	creator->set_graph_type("DIMACS");

	wmg = creator->create_weighted_mutable_graph();

	list<int> f;
	list<int>::iterator it;
	list<int>::iterator jt;
	list<int>::iterator kt;
	util.find_component(wmg, 25, &f);

	wmg_sub = creator->create_component(wmg, &f, true);

	EXPECT_EQ(f.size(), wmg_sub->get_num_nodes())
		;

	vector<Graph::Node> wmg_nodes = wmg->get_nodes();
	vector<Graph::Node> wmg_sub_nodes = wmg_sub->get_nodes();
	int i = 0;
	int j = 0;
	int n = wmg_sub->get_num_nodes();

	it = f.begin();
	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(wmg_nodes[*it].get_label(), wmg_sub_nodes[i].get_label())
			;
		EXPECT_EQ(wmg_nodes[*it].get_degree(), wmg_sub_nodes[i].get_degree())
			;
		list<int> wmg_nbrs = wmg_nodes[*it].get_nbrs();
		list<int> wmg_sub_nbrs = wmg_sub_nodes[i].get_nbrs();
		wmg_nbrs.sort();
		wmg_sub_nbrs.sort();

		for (jt = wmg_nbrs.begin(), kt = wmg_sub_nbrs.begin(); jt
                 != wmg_nbrs.end(); ++jt, ++kt)
		{
			ASSERT_EQ(wmg_nodes[*jt].get_label(), wmg_sub_nodes[*kt].get_label())
				;
		}
		++it;
	}

}

TEST_F(GraphCreatorFileTest, testCreateAllComponents)
{
	Graph::GraphUtil util;
	Graph::WeightedMutableGraph *wmg;

	creator->set_file_name("data/1et.64.txt");
	creator->set_graph_type("DIMACS");
	wmg = creator->create_weighted_mutable_graph();

	vector<list<int> *> members;
	int x = util.find_all_components(wmg, &members);
	EXPECT_EQ(x, members.size())
		;
	list<Graph::WeightedMutableGraph *> cmembers;
	cmembers = creator->create_all_components(wmg, true);

	EXPECT_EQ(cmembers.size(), members.size())
		;

	int wmg_size = 0;
	int mem_size = members[3]->size();
	list<Graph::WeightedMutableGraph *>::iterator giter;
	Graph::WeightedMutableGraph cmg;
	vector<Graph::Node> nodes;

	for (giter = cmembers.begin(); giter != cmembers.end(); ++giter)
	{
		if (mem_size == (*giter)->get_num_nodes())
		{
			nodes = (*giter)->get_nodes();
			break;
		}
	}

	list<int>::iterator it;
	int i = 0;
	for (it = members[3]->begin(); it != members[3]->end(); ++it)
	{
		EXPECT_EQ((*it), nodes[i].get_label() - 1);
		i++;
	}
}


TEST_F(GraphCreatorFileTest, testCreateRecAllComponents)
{
	Graph::GraphUtil util;
	Graph::WeightedMutableGraph *wmg;

	creator->set_file_name("data/1et.64.txt");
	creator->set_graph_type("DIMACS");
	wmg = creator->create_weighted_mutable_graph();

	vector<list<int> *> members;
	int x = util.find_all_components(wmg, &members);
	EXPECT_EQ(x, members.size())
		;
	list<Graph::WeightedMutableGraph *> cmembers;
	cmembers = creator->create_rec_all_components(wmg, true);

	EXPECT_EQ(cmembers.size(), members.size())
		;

	int wmg_size = 0;
	int mem_size = members[3]->size();
	list<Graph::WeightedMutableGraph *>::iterator giter;
	Graph::WeightedMutableGraph cmg;
	vector<Graph::Node> nodes;

	for (giter = cmembers.begin(); giter != cmembers.end(); ++giter)
	{
		if (mem_size == (*giter)->get_num_nodes())
		{
			nodes = (*giter)->get_nodes();
			break;
		}
	}

	list<int>::iterator it;
	int i = 0;
	for (it = members[3]->begin(); it != members[3]->end(); ++it)
	{
		EXPECT_EQ((*it), nodes[i].get_label() - 1);
		i++;
	}
}
