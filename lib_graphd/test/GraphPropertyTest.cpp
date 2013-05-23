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
#include "GraphProperties.h"
#include <gtest/gtest.h>

using namespace std;
class GraphPropertyTest : public testing::Test
{
public:
Graph::GraphCreatorFile creator;
Graph::Graph *mg;
Graph::GraphProperties properties;

virtual void SetUp(){
    //SetUp is called before every test
    LOG_INIT("test.log", NULL, 0);
    creator.set_file_name("data/1dc.128.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();
}

virtual void TearDown(){
    //TearDown is called before every test
    LOG_CLOSE();
}

void remove_all_edges(int i){
    list<int>::iterator it;
    Graph::Node *n = mg->get_node(i);
    list<int> nbrs = n->get_nbrs();
    for(it = nbrs.begin(); it != nbrs.end(); ++it){
        mg->remove_edge(i, *it);
    }
}
};

TEST_F(GraphPropertyTest, testSimple)
{
    Graph::Node *n;
    n = mg->get_node(108);
    list<int> nbrs = n->get_nbrs();

    EXPECT_EQ(24, nbrs.size());
    EXPECT_TRUE(mg->is_edge(108, 100));

    EXPECT_EQ(24, mg->get_degree(108));
    mg->add_edge(108, 100);
    mg->add_edge(108, 108);

    EXPECT_EQ(1473, mg->get_num_edges());
    EXPECT_EQ(26, mg->get_degree(108));

    EXPECT_FALSE(properties.check_simple(mg));
    properties.make_simple(mg);
    EXPECT_TRUE(properties.check_simple(mg));
    EXPECT_EQ(24, mg->get_degree(108));
}

TEST_F(GraphPropertyTest, testClique)
{
    Graph::Node *n;
    n = mg->get_node(108);
    list<int> nbrs = n->get_nbrs();

    EXPECT_EQ(24, nbrs.size());
    EXPECT_FALSE(properties.is_clique(mg, &nbrs));
    nbrs.push_back(108);
    properties.make_clique(mg, &nbrs);
    EXPECT_TRUE(properties.is_clique(mg, &nbrs));
}

TEST_F(GraphPropertyTest, testIsConnected)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    EXPECT_FALSE(properties.is_connected(mg));
}

TEST_F(GraphPropertyTest, testIsIndependentSet)
{
    Graph::Node *n = mg->get_node(20);
    list<int> nbrs = n->get_nbrs();
    EXPECT_FALSE(properties.is_independent_set(mg, &nbrs));

    list<int>::iterator it;
    list<int>::iterator jt;
    for(it = nbrs.begin(); it != nbrs.end(); ++it){
        // Remove all edges between them
        for(jt = nbrs.begin(); jt != nbrs.end(); ++jt){
            if(*it != *jt){
                mg->remove_edge(*it, *jt);
            }
        }
    }

    EXPECT_TRUE(properties.is_independent_set(mg, &nbrs));
    EXPECT_GT(mg->get_num_edges(), 0);
}

TEST_F(GraphPropertyTest, testIsPath)
{
    Graph::Node *n = mg->get_node(20);
    list<int> nbrs = n->get_nbrs();

    EXPECT_GT(nbrs.size(), 0);
    EXPECT_TRUE(properties.is_path(mg, nbrs.front(), nbrs.back()));

    remove_all_edges(nbrs.front());

    EXPECT_FALSE(properties.is_path(mg, nbrs.front(), nbrs.back()));
    EXPECT_GT(mg->get_num_edges(), 0);
}

TEST_F(GraphPropertyTest, testEdgeDensity){
    float ed, val;

    val = (2.0 * 1471) / (128 * (128 - 1.0));

    properties.edgeDensity(mg, ed);
    EXPECT_EQ(val, ed);
}

TEST_F(GraphPropertyTest, testAvgDegree){
    float ad, val;

    val = (2.0 * 1471) / 128.0;

    properties.avgDegree(mg, ad);
    EXPECT_EQ(val, ad);
}

TEST_F(GraphPropertyTest, testDegDist){
    vector<int> dist;
    properties.degDist(mg, dist);
    EXPECT_EQ(4, dist[33]);
}
