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

#define GTEST_HAS_TR1_TUPLE 0

#include <cmath>

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
Graph::Graph *mg_dis;
Graph::GraphProperties properties;

virtual void SetUp(){
    //SetUp is called before every test
    LOG_INIT("test.log", NULL, 0);
    creator.set_file_name("data/1dc.128.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    //disconnected graph for testing
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg_dis = creator.create_mutable_graph();
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
    EXPECT_FALSE(properties.is_connected(mg_dis));
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

TEST_F(GraphPropertyTest, testDiameter){
    int diam, val;

    properties.diameter(mg, diam);
    EXPECT_EQ(7, diam);
}

TEST_F(GraphPropertyTest, testEffectiveDiameter){
    float ediam, val;

    val = 2.0 + (8128.0 * 0.9 - 5173) / (2188.0);

    properties.effective_diameter(mg, ediam);
    EXPECT_EQ(val, ediam);
}

TEST_F(GraphPropertyTest, testEdgeDensity){
    float ed, val;

    val = (2.0 * 1471) / (128 * (128 - 1.0));

    properties.edge_density(mg, ed);
    EXPECT_EQ(val, ed);
}

TEST_F(GraphPropertyTest, testAvgDegree){
    float ad, val;

    val = (2.0 * 1471) / 128.0;

    properties.avg_degree(mg, ad);
    EXPECT_EQ(val, ad);
}

TEST_F(GraphPropertyTest, testDegDist){
    vector<int> dist;

    properties.deg_dist(mg, dist);
    EXPECT_EQ(4, dist[33]);
}

#ifdef HAS_BOOST
TEST_F(GraphPropertyTest, testPowerLaw){
    int xmin;
    double alpha, KS;
    properties.powerlaw(mg, xmin, alpha, KS);
    EXPECT_EQ(18, xmin);
    EXPECT_NEAR(3.5, alpha, 0.1);
    EXPECT_NEAR(0.470626, KS, 0.000001);
}
#endif // ifdef HAS_BOOST

TEST_F(GraphPropertyTest, testDeltaHyperbolicity){
    double max_delta;
    vector<vector<double> > delta;

    properties.delta_hyperbolicity(mg, max_delta, delta);
    EXPECT_NEAR(2.0, max_delta, 0.1);
    EXPECT_NEAR(204.0, delta[3][5], 0.1);

    delta.clear();
    properties.delta_hyperbolicity(mg_dis, max_delta, delta);
    EXPECT_NEAR(1.0, max_delta, 0.1);
    EXPECT_NEAR(126.0, delta[1][3], 0.1);
}

TEST_F(GraphPropertyTest, testClustering){
    double global_cc;
    double avg_cc;
    vector<double> local_ccs;

    properties.clustering_coefficients(mg, global_cc, avg_cc, local_ccs);
    EXPECT_NEAR(0.56209150326797386, local_ccs[32], 0.0000001);
    EXPECT_NEAR(0.491997272348, avg_cc, 0.00000001);
    EXPECT_NEAR(0.460570029725, global_cc, 0.0000001);
}
