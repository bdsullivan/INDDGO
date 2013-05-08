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
#include "Graph.h"
#include "GraphCreator.h"
#include "GraphCreatorFile.h"
#include <gtest/gtest.h>

using namespace std;
class GraphTest: public testing::Test {
public:
    Graph::GraphCreatorFile creator;
    Graph::Graph *g;

    virtual void SetUp(){
        LOG_INIT("test.log", NULL, 0);
        //creator.set_file_name("../data/1dc.128.txt");
        creator.set_file_name("data/1dc.128.txt");
        creator.set_graph_type("DIMACS");
		g = creator.create_graph();
    }

    virtual void TearDown(){
        LOG_CLOSE();
    }
};


TEST_F(GraphTest, testNumNodes)
{
    EXPECT_EQ(128, g->get_num_nodes());
}

TEST_F(GraphTest, testNumEdges)
{
    EXPECT_EQ(1471, g->get_num_edges());
}

TEST_F(GraphTest, testGraphType)
{
    EXPECT_EQ("DIMACS", g->get_graph_type());
}

TEST_F(GraphTest, testGetNode)
{
    Graph::Node *n;
    n = g->get_node(108);
    list<int> nbrs = n->get_nbrs();
    EXPECT_EQ(24, nbrs.size());
}

TEST_F(GraphTest, testIsEdge)
{
    Graph::Node *n;
    n = g->get_node(108);
    list<int> nbrs = n->get_nbrs();

    EXPECT_EQ(24, nbrs.size());
    EXPECT_TRUE(g->is_edge(108, 100));

}

TEST_F(GraphTest, testGetDegree)
{
    EXPECT_EQ(19, g->get_degree(28));
}


TEST_F(GraphTest, testNodeNbrs)
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


TEST_F(GraphTest, testGetDegrees)
{
    vector<int> degree = g->get_degree();
    EXPECT_EQ(31, degree[86]);
}

TEST_F(GraphTest, testGetNumComponents)
{
    EXPECT_EQ(1, g->get_num_components());
}

TEST_F(GraphTest, testGetNextLabel)
{
    EXPECT_EQ(129, g->get_next_label());
}

TEST_F(GraphTest, testSetNextLabel)
{
    EXPECT_EQ(129, g->get_next_label());
    g->set_next_label(140);
    EXPECT_EQ(140, g->get_next_label());
}

TEST_F(GraphTest, testComplement)
{
    int d = g->get_degree(86);
    int n = g->get_num_nodes();
    int comp = n*(n-1)/2;
    EXPECT_EQ(31, d);
    EXPECT_TRUE(g->is_edge(108, 100));
    EXPECT_EQ(128, g->get_num_nodes());
    EXPECT_EQ(1471, g->get_num_edges());

    g->complement();
    
    EXPECT_EQ(128, g->get_num_nodes());
    EXPECT_EQ(comp - 1471, g->get_num_edges());

    d = g->get_degree(86);
    EXPECT_EQ(96, d);
    EXPECT_FALSE(g->is_edge(108, 100));
}
