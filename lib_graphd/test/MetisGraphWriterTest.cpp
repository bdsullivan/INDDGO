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
#include "MetisGraphWriter.h"
#include <string>
#include <gtest/gtest.h>

using namespace std;

class MetisGraphWriterTest : public testing::Test
{
public:
Graph::GraphCreatorFile creator;
Graph::Graph *g;
Graph::GraphWriter writer;
string out_file;

virtual void SetUp(){
    //SetUp is called before every test
    LOG_INIT("test.log", NULL, 0);
    creator.set_file_name("data/1dc.128.adj");
    creator.set_graph_type("adjmatrix");
    g = creator.create_graph();

    out_file = string("data/1dc.128.met");
}

virtual void TearDown(){
    //TearDown is called before every test
    LOG_CLOSE();
}
};

TEST_F(MetisGraphWriterTest, testNumNodes)
{
    writer.write_graph(g, out_file, string("METIS"), false);

    creator.set_file_name("data/1dc.128.met");
    creator.set_graph_type("Metis");
    g = creator.create_graph();

    EXPECT_EQ(128, g->get_num_nodes());
}

TEST_F(MetisGraphWriterTest, testGetDegrees)
{
    writer.write_graph(g, out_file, string("METIS"), false);

    creator.set_file_name("data/1dc.128.met");
    creator.set_graph_type("Metis");
    g = creator.create_graph();

    vector<int> degree = g->get_degree();
    EXPECT_EQ(31, degree[86]);
}

TEST_F(MetisGraphWriterTest, testGetNode)
{
    writer.write_graph(g, out_file, string("METIS"), false);

    creator.set_file_name("data/1dc.128.met");
    creator.set_graph_type("Metis");
    g = creator.create_graph();

    Graph::Node *n;
    n = g->get_node(108);
    list<int> nbrs = n->get_nbrs();
    EXPECT_EQ(24, nbrs.size());
}

TEST_F(MetisGraphWriterTest, testIsEdge)
{
    writer.write_graph(g, out_file, string("METIS"), false);

    creator.set_file_name("data/1dc.128.met");
    creator.set_graph_type("Metis");
    g = creator.create_graph();

    Graph::Node *n;
    n = g->get_node(108);
    list<int> nbrs = n->get_nbrs();

    EXPECT_EQ(24, nbrs.size());
    EXPECT_TRUE(g->is_edge(108, 100));
}
