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
#include "GraphUtil.h"
#include "Node.h"
#include "VertexWeightedGraph.h"
#include <gtest/gtest.h>

using namespace std;
class GraphUtilTest : public testing::Test
{
public:
Graph::GraphCreatorFile creator;
Graph::Graph *mg;
Graph::GraphUtil util;

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

int get_highest_degree(){
    int n = mg->get_num_nodes();
    int max = 0;
    int i = 0;
    for(i = 0; i < n; i++){
        if(max < mg->get_degree(i)){
            max = mg->get_degree(i);
        }
    }
    return max;
}

int get_lowest_degree(){
    int n = mg->get_num_nodes();
    int min = INT_MAX;
    int i = 0;
    for(i = 0; i < n; i++){
        if(min > mg->get_degree(i)){
            min = mg->get_degree(i);
        }
    }
    return min;
}
};

TEST_F(GraphUtilTest, testFindAllComponent)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    vector<list<int> *> members;
    int x = util.find_all_components(mg, &members);
    EXPECT_EQ(x, members.size());
}

TEST_F(GraphUtilTest, testFindComponent)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    vector<list<int> *>members;
    int x = util.find_all_components(mg, &members);
    EXPECT_EQ(x, members.size());

    list<int> f;
    list<int>::iterator it;
    list<int>::iterator jt;
    util.find_component(mg, 25, &f);
    bool a = true;
    int i = 0;

    (members[3])->sort();
    f.sort();
//	Test whether both components have same elements
    for(it = (members[3])->begin(), jt = f.begin();
        it != (members[3])->end(); ++it, ++jt){
        if(*it != *jt){
            cout << *it << " " << *jt << endl;
            a = false;
        }
        i++;
    }
    EXPECT_TRUE(a);
}

TEST_F(GraphUtilTest, testFindIsolatedNodes)
{
    list<int> isolated_nodes;
    util.find_isolated_nodes(mg, &isolated_nodes);
    EXPECT_EQ(0, isolated_nodes.size());

    remove_all_edges(108);
    remove_all_edges(100);

    util.find_isolated_nodes(mg, &isolated_nodes);
    EXPECT_EQ(2, isolated_nodes.size());
}

TEST_F(GraphUtilTest, testRecomputeDegree)
{
    Graph::Node *n = mg->get_node(36);
    list<int> nbrs = n->get_nbrs();

    EXPECT_EQ(29, nbrs.size());
    EXPECT_FALSE(mg->is_edge(36, 3));

    mg->add_edge(36, 3);
    EXPECT_EQ(30, mg->get_degree(36));
    mg->remove_edge(36, 3);

    nbrs.push_back(3);
    n->set_nbr(nbrs);
    nbrs = n->get_nbrs();
    EXPECT_EQ(30, nbrs.size());
    EXPECT_EQ(29, mg->get_degree(36));

    util.recompute_degrees(mg);
    EXPECT_EQ(30, mg->get_degree(36));
}

TEST_F(GraphUtilTest, testGetRandomHighDegreeVertex)
{
    int i = util.get_random_high_degree_vertex(mg);
    int x = get_highest_degree();
    EXPECT_EQ(x, mg->get_degree(i));
    mg->remove_vertex(i);
    i = util.get_random_high_degree_vertex(mg);
    x = get_highest_degree();
    EXPECT_EQ(x, mg->get_degree(i));
}

TEST_F(GraphUtilTest, testGetRandomLowDegreeVertex)
{
    int i = util.get_random_low_degree_vertex(mg);
    int x = get_lowest_degree();
    EXPECT_EQ(x, mg->get_degree(i));
    mg->remove_vertex(i);
    i = util.get_random_low_degree_vertex(mg);
    x = get_lowest_degree();
    EXPECT_EQ(x, mg->get_degree(i));
}

TEST_F(GraphUtilTest, testLabelComponent)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    vector<int> components(mg->get_capacity(), -1);
    util.label_component(mg, 25, 3, &components);

    list<int> members;
    int x = util.find_component(mg, 25, &members);

    list<int>::iterator it;
    for(it = members.begin(); it != members.end(); ++it){
        EXPECT_EQ(3, components[*it]);
    }
}

TEST_F(GraphUtilTest, testRecLabelComponent)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    vector<int> components(mg->get_capacity(), -1);
    util.rec_label_component(mg, 25, 3, &components);

    list<int> members;
    int x = util.find_component(mg, 25, &members);

    list<int>::iterator it;
    for(it = members.begin(); it != members.end(); ++it){
        EXPECT_EQ(3, components[*it]);
    }
}

TEST_F(GraphUtilTest, testLabelAllComponent)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    vector<int> components(mg->get_capacity(), -1);
    util.label_all_components(mg, &components);

    list<int> members;
    int x = util.find_component(mg, 12, &members);

    list<int>::iterator it;
    for(it = members.begin(); it != members.end(); ++it){
        EXPECT_EQ(components[12], components[*it]);
    }
}

TEST_F(GraphUtilTest, testRecLabelAllComponent)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    vector<int> components(mg->get_capacity(), -1);
    util.rec_label_all_components(mg, &components);

    list<int> members;
    int x = util.find_component(mg, 30, &members);

    list<int>::iterator it;
    for(it = members.begin(); it != members.end(); ++it){
        EXPECT_EQ(components[30], components[*it]);
    }
}

TEST_F(GraphUtilTest, testPopulateCrs)
{
    util.populate_CRS(mg);
    int i = 0;
    vector<int> xvec = mg->get_xadj();
    vector<int> adj = mg->get_adjncy();

    list<int> nbrs = mg->get_node(40)->get_nbrs();
    list<int>::iterator it;
    i = 0;
    for(it = nbrs.begin(); it != nbrs.end(); ++it){
        EXPECT_EQ(*it, adj[xvec[40] + i]);
        i++;
    }

    nbrs = mg->get_node(119)->get_nbrs();
    i = 0;
    for(it = nbrs.begin(); it != nbrs.end(); ++it){
        EXPECT_EQ(*it, adj[xvec[119] + i]);
        i++;
    }
}

TEST_F(GraphUtilTest, testVertexSeparator)
{
    creator.set_file_name("data/1et.64.txt");
    creator.set_graph_type("DIMACS");
    mg = creator.create_mutable_graph();

    list<int> members;
    list<int> cmembers;
    int x = util.find_component(mg, 30, &members);

    cmembers = members;
    vector<list<int> *> newmembers;
    util.vertex_separator(mg, &cmembers, &newmembers);

    int i = 0;
    list<int>::iterator it;
    list<int>::iterator jt;
    for(i = 0; i < newmembers.size(); i++){
        for(it = (newmembers[i])->begin(); it != (newmembers[i])->end();
            ++it){
            for(jt = members.begin(); jt != members.end(); ++jt){
                EXPECT_NE(*jt, *it);
            }
        }
    }
}
