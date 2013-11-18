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

#include <fstream>

using namespace std;
class GraphTest : public testing::Test {
public:
Graph::GraphCreatorFile creator;
Graph::Graph *g;
Graph::Graph *g_five;

virtual void SetUp(){
    LOG_INIT("test.log", NULL, 0);
    //creator.set_file_name("../data/1dc.128.txt");
    creator.set_file_name("data/1dc.128.txt");
    creator.set_graph_type("DIMACS");
    g = creator.create_graph();
    g_five = new Graph::Graph(5);
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

TEST_F(GraphTest, testNextLabel){
    int n = g_five->get_next_label();
    EXPECT_EQ(6, n);
}

TEST_F(GraphTest, testIsEdge)
{
    Graph::Node *n;
    n = g->get_node(108);
    list<int> nbrs = n->get_nbrs();

    EXPECT_EQ(24, nbrs.size());
    EXPECT_TRUE(g->is_edge(108, 100));

    n = g->get_node(80);
    list<int> nbrs2 = n->get_nbrs();
    EXPECT_EQ(22, nbrs2.size());
    EXPECT_FALSE(g->is_edge(80, 12));
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

    ofstream outFile;
    outFile.open("degrees.txt");
    for(int i = 0; i < degree.size(); i++){
        outFile << i << " " << degree[i] << endl;
    }

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
    int comp = n * (n - 1) / 2;
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

TEST_F(GraphTest, testAddEdge)
{
    EXPECT_FALSE(g->is_edge(80, 12));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(22, g->get_degree(80));
    EXPECT_EQ(19, g->get_degree(12));
    g->add_edge(80, 12);
    EXPECT_TRUE(g->is_edge(80, 12));
    EXPECT_TRUE(g->is_edge(12, 80));
    EXPECT_EQ(1472, g->get_num_edges());
    EXPECT_EQ(23, g->get_degree(80));
    EXPECT_EQ(20, g->get_degree(12));
}

TEST_F(GraphTest, testAddEdgeAdvance)
{
    EXPECT_FALSE(g->is_edge(80, 12));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(22, g->get_degree(80));
    EXPECT_EQ(19, g->get_degree(12));
    g->add_edge_advance(80, 12);
    EXPECT_TRUE(g->is_edge(80, 12));
    EXPECT_TRUE(g->is_edge(12, 80));
    EXPECT_EQ(1472, g->get_num_edges());
    EXPECT_EQ(23, g->get_degree(80));
    EXPECT_EQ(20, g->get_degree(12));
}

TEST_F(GraphTest, testRemoveEdge)
{
    EXPECT_TRUE(g->is_edge(108, 100));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(24, g->get_degree(108));
    EXPECT_EQ(24, g->get_degree(100));

    g->remove_edge(108, 100);

    EXPECT_FALSE(g->is_edge(108, 100));
    EXPECT_FALSE(g->is_edge(100, 108));
    EXPECT_EQ(1470, g->get_num_edges());
    EXPECT_EQ(23, g->get_degree(100));
    EXPECT_EQ(23, g->get_degree(108));
}

TEST_F(GraphTest, testRemoveVertex)
{
    EXPECT_TRUE(g->is_edge(108, 100));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(24, g->get_degree(108));
    EXPECT_EQ(24, g->get_degree(100));

    g->remove_vertex(108);

    EXPECT_FALSE(g->is_edge(108, 100));
    EXPECT_FALSE(g->is_edge(100, 108));
    EXPECT_EQ(1447, g->get_num_edges());
    EXPECT_EQ(23, g->get_degree(100));
    EXPECT_EQ(0, g->get_degree(108));
}

TEST_F(GraphTest, testContractEdge)
{
    EXPECT_TRUE(g->is_edge(108, 100));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(24, g->get_degree(108));
    EXPECT_EQ(24, g->get_degree(100));

    Graph::Node *na;
    na = g->get_node(108);
    list<int> nbrs_a = na->get_nbrs();
    vector<int> nbrsvec_a(nbrs_a.begin(), nbrs_a.end());

    Graph::Node *nb;
    nb = g->get_node(100);
    list<int> nbrs_b = nb->get_nbrs();
    vector<int> nbrsvec_b(nbrs_b.begin(), nbrs_b.end());

    int common_nbrs = 0;
    int final_degree = 0;
    int new_edges = 0;

    for(int ia = 0; ia < nbrsvec_a.size(); ia++){
        for(int ib = 0; ib < nbrsvec_b.size(); ib++){
            if(nbrsvec_a[ia] == nbrsvec_b[ib]){
                common_nbrs++;
            }
        }
    }

    final_degree = g->get_degree(108) + g->get_degree(100) - common_nbrs - 2;
    new_edges = g->get_num_edges() - common_nbrs - 1;

    int x = g->contract_edge(108, 100);

    EXPECT_FALSE(g->is_edge(108, 100));
    EXPECT_FALSE(g->is_edge(100, 108));
    EXPECT_EQ(0, g->get_degree(100));
    EXPECT_EQ(new_edges, g->get_num_edges());
    EXPECT_EQ(final_degree, g->get_degree(108));
}

TEST_F(GraphTest, testFuseVertices)
{
    EXPECT_FALSE(g->is_edge(80, 12));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(22, g->get_degree(80));
    EXPECT_EQ(19, g->get_degree(12));

    Graph::Node *na;
    na = g->get_node(80);
    list<int> nbrs_a = na->get_nbrs();
    vector<int> nbrsvec_a(nbrs_a.begin(), nbrs_a.end());

    Graph::Node *nb;
    nb = g->get_node(12);
    list<int> nbrs_b = nb->get_nbrs();
    vector<int> nbrsvec_b(nbrs_b.begin(), nbrs_b.end());

    int common_nbrs = 0;
    int final_degree = 0;
    int new_edges = 0;

    for(int ia = 0; ia < nbrsvec_a.size(); ia++){
        for(int ib = 0; ib < nbrsvec_b.size(); ib++){
            if(nbrsvec_a[ia] == nbrsvec_b[ib]){
                common_nbrs++;
            }
        }
    }

    final_degree = g->get_degree(80) + g->get_degree(12) - common_nbrs;
    new_edges = g->get_num_edges() - common_nbrs;

    int x = g->fuse_vertices(80, 12, false);

    EXPECT_FALSE(g->is_edge(80, 12));
    EXPECT_FALSE(g->is_edge(12, 80));
    EXPECT_EQ(0, g->get_degree(12));
    EXPECT_EQ(new_edges, g->get_num_edges());
    EXPECT_EQ(final_degree, g->get_degree(80));
}

TEST_F(GraphTest, testEliminateVertex)
{
    EXPECT_TRUE(g->is_edge(108, 100));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(24, g->get_degree(108));
    EXPECT_EQ(24, g->get_degree(100));

    Graph::Node *na;
    na = g->get_node(108);
    list<int> nbrs_a = na->get_nbrs();
    vector<int> nbrsvec_a(nbrs_a.begin(), nbrs_a.end());

    Graph::Node *nb;
    nb = g->get_node(100);
    list<int> nbrs_b = nb->get_nbrs();
    vector<int> nbrsvec_b(nbrs_b.begin(), nbrs_b.end());

    int common_nbrs = 0;
    int final_degree = 0;
    int new_edges = 0;

    for(int ia = 0; ia < nbrsvec_a.size(); ia++){
        for(int ib = 0; ib < nbrsvec_b.size(); ib++){
            if(nbrsvec_a[ia] == nbrsvec_b[ib]){
                common_nbrs++;
            }
        }
    }

    new_edges = g->get_degree(108) - common_nbrs - 1;
    new_edges = new_edges + g->get_degree(100) - 1;

    g->eliminate_vertex(108, NULL, true);

    EXPECT_FALSE(g->is_edge(108, 100));
    EXPECT_FALSE(g->is_edge(100, 108));
    EXPECT_EQ(1597, g->get_num_edges());
    EXPECT_EQ(new_edges, g->get_degree(100));
    EXPECT_EQ(0, g->get_degree(108));
}

TEST_F(GraphTest, testEdgeSubdivision)
{
    EXPECT_TRUE(g->is_edge(108, 100));
    EXPECT_EQ(1471, g->get_num_edges());
    EXPECT_EQ(24, g->get_degree(108));
    EXPECT_EQ(24, g->get_degree(100));

    //We are creating a new vertex w of degree 2, with neighbors 108 and 100. The neighbors of
    //108 and 100 should change by losing each other and adding w.
    int w = 130; //out of current range of vertice
    int u_old = g->get_degree(108);
    int v_old = g->get_degree(100);

    int x = g->edge_subdivision(108, 100, w);

    //Check existence of correct edges
    EXPECT_FALSE(g->is_edge(108, 100));
    EXPECT_FALSE(g->is_edge(100, 108));
    EXPECT_TRUE(g->is_edge(100, w));
    EXPECT_TRUE(g->is_edge(108, w));
    EXPECT_TRUE(g->is_edge(w, 100));
    EXPECT_TRUE(g->is_edge(w, 108));

    //Check vertex degrees
    EXPECT_EQ(2, g->get_degree(w));
    EXPECT_EQ(u_old, g->get_degree(108));
    EXPECT_EQ(v_old, g->get_degree(100));
}
