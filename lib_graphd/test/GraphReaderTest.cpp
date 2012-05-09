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

#include "DIMACSGraphReader.h"
#include "GraphReader.h"
#include "Log.h"
#include "Node.h"
#include <gtest/gtest.h>

using namespace std;
// A new one of these is created for each test
class GraphReaderTest : public testing::Test {
public:
    Graph::GraphReader *gr;


    virtual void SetUp(){
        LOG_INIT("test.log", NULL, 0);
        gr = new Graph::DIMACSGraphReader ();
        gr->read_graph("../data/1dc.128.txt");
    }

    virtual void TearDown(){
        LOG_CLOSE();
    }
};


TEST_F(GraphReaderTest, testCapacity)
{
    EXPECT_EQ(128, gr->get_capacity());
}

TEST_F(GraphReaderTest, testNumEdges)
{
    EXPECT_EQ(1471, gr->get_num_edges());
}


TEST_F(GraphReaderTest, testNodeNbrs)
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





