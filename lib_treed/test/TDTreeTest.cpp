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
#include "WeightedMutableGraph.h"
#include "GraphCreator.h"
#include "GraphCreatorFile.h"
#include "GraphEOUtil.h"
#include "TDTree.h"
#include "TDTreeNode.h"
#include "TDDynamicProgramming.h"
#include "weighted_ind_set.h"
#include <gtest/gtest.h>

using namespace std;
class TDTreeTest: public testing::Test
{
public:
	Graph::GraphCreatorFile creator;
	Graph::WeightedMutableGraph *wmg;
	TDTree *T;
	vector<int> *ordering;
	Graph::GraphEOUtil eoutil;
	DP_info *info;
	vector<int> *walk;

	virtual void SetUp()
	{
		LOG_INIT("test.log", NULL, 0);
		creator.set_file_name("../data/1dc.128.txt");
		creator.set_graph_type("DIMACS");
		wmg = creator.create_weighted_mutable_graph();
		T = new TDTree(wmg);

		ordering = new vector<int> (wmg->get_num_nodes(), 0);
		info = new DP_info();

		string filename = "../data/1dc.128.txt";
		info->DIMACS_file = (char *) filename.c_str();
	}

	virtual void TearDown()
	{
		LOG_CLOSE();
	}

	void doDP(TDTree *T)
	{
		int i;
		walk = new vector<int> (T->num_tree_nodes);
		T->post_order_walk(walk);

		// Do the DP using the provided options
		for (i = 0; i < T->num_tree_nodes; i++)
			T->compute_table(compute_weighted_ind_set_table, (*walk)[i]);

		list<int> root_intersection, root_difference;
		get_optimal_obj(T, &root_intersection, &root_difference, T->info);
		//print_WIS_results(stdout, T, T->info);
	}

	void freeTables(TDTree *T)
	{
		int i = 0;
		for (i = 0; i < T->num_tree_nodes; i++)
			T->free_table((*walk)[i]);
	}
};

TEST_F(TDTreeTest, testWidth)
{
	eoutil.find_elimination_ordering(wmg, ordering, GD_MIN_DEGREE, false);
	T->width = eoutil.triangulate(wmg, ordering);
	EXPECT_EQ(67, T->width)
		;
}

#if HAS_METIS
TEST_F(TDTreeTest, testMetisWidth)
{
	eoutil.find_elimination_ordering(wmg, ordering, GD_MIN_DEGREE, false);
	T->width = eoutil.METIS_triangulate(wmg, ordering);
	EXPECT_EQ(67, T->width)
		;
}
#endif

TEST_F(TDTreeTest, testGavril)
{
	eoutil.find_elimination_ordering(wmg, ordering, GD_MIN_DEGREE, false);
	T->width = eoutil.triangulate(wmg, ordering);
	EXPECT_EQ(67, T->width)
		;

	EXPECT_EQ(0, T->num_leafs)
		;
	EXPECT_EQ(0, T->num_tree_nodes)
		;
	T->construct_gavril(ordering);
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;
}

TEST_F(TDTreeTest, testBK)
{
	eoutil.find_elimination_ordering(wmg, ordering, GD_MIN_DEGREE, false);
	T->width = eoutil.triangulate(wmg, ordering);
	EXPECT_EQ(67, T->width)
		;

	EXPECT_EQ(0, T->num_leafs)
		;
	EXPECT_EQ(0, T->num_tree_nodes)
		;
	T->construct_BK(ordering);
	EXPECT_EQ(12, T->num_leafs)
		;
	EXPECT_EQ(128, T->num_tree_nodes)
		;
}

TEST_F(TDTreeTest, testNice)
{
	eoutil.find_elimination_ordering(wmg, ordering, GD_MIN_DEGREE, false);
	T->width = eoutil.triangulate(wmg, ordering);
	EXPECT_EQ(67, T->width)
		;

	EXPECT_EQ(0, T->num_leafs)
		;
	EXPECT_EQ(0, T->num_tree_nodes)
		;

	T->construct_knice(ordering, T->width, false);
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(141, T->num_tree_nodes)
		;
}

TEST_F(TDTreeTest, testGavrilMind)
{
	/**
	 * This test is equivalent to -gavril -mind options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MIN_DEGREE;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	T->free_table(T->root_node);
	T->free_children(T->root_node);
	EXPECT_EQ(39, T->num_tree_nodes)
		;
	EXPECT_EQ(0, T->info->opt_obj)
		;
	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;
	freeTables(T);
}

TEST_F(TDTreeTest, testGavrilMindDfs)
{
	/**
	 * This test is equivalent to -gavril -mind -dfs options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MIN_DEGREE;
	info->DFS = true;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;

	EXPECT_EQ(0, T->info->opt_obj)
		;
	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;

	freeTables(T);
}

TEST_F(TDTreeTest, testGavrilMindDelch)
{
	/**
	 * This test is equivalent to -gavril -mind -dfs options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MIN_DEGREE;
	info->free_children = true;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;
	EXPECT_EQ(0, T->info->opt_obj)
		;

	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;

	T->free_table(T->root_node);
}

TEST_F(TDTreeTest, testGavrilMindDfsDelch)
{
	/**
	 * This test is equivalent to -gavril -mind -dfs -del_ch options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MIN_DEGREE;
	info->DFS = true;
	info->free_children = true;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;

	EXPECT_EQ(0, T->info->opt_obj)
		;

	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;
	T->free_table(T->root_node);
}

/**
 *    -mind : generates an elim. ordering using min degree heuristic
 -mmd: generates an elim. ordering using multiple min degree heuristic
 -minf : generates an elim. ordering using min fill heuristic
 -bmf: generates an elim. ordering using a batched min fill heuristic
 -beta: generates an elim. ordering using the beta heuristic
 -metmmd: generates an elim. ordering using METIS mmd heuristic
 -metnnd: generates an elim. ordering using METS node ND heuristic
 -mcsm : generates an elim. ordering using mcsm euristic
 -mcs  : generates an elim. ordering using mcs
 -lexm : generates an elim. ordering using lex-m bfs heuristic
 -pktsort : generates the elim ordering corresponding to partial
 ktree generation
 -amd: generates the elim. ordering using UFs approximate min
 degree heuristic (AMD)
 -minmaxd: generates the elim. ordering using the minimum maximum
 degree heuristic.
 */

TEST_F(TDTreeTest, testGavrilMmd)
{
	/**
	 * This test is equivalent to -gavril -mmd options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MUL_MIN_DEGREE;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;
	EXPECT_EQ(0, T->info->opt_obj)
		;

	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;

	freeTables(T);
}

TEST_F(TDTreeTest, testGavrilMmdDfs)
{
	/**
	 * This test is equivalent to -gavril -mmd -dfs options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MUL_MIN_DEGREE;
	info->DFS = true;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;
	EXPECT_EQ(0, T->info->opt_obj)
		;

	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;

	freeTables(T);
}

TEST_F(TDTreeTest, testGavrilMmdDelch)
{
	/**
	 * This test is equivalent to -gavril -mmd -del_ch options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MUL_MIN_DEGREE;
	info->free_children = true;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;
	EXPECT_EQ(0, T->info->opt_obj)
		;

	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;

	T->free_table(T->root_node);
}

TEST_F(TDTreeTest, testGavrilMmdDfsDelch)
{
	/**
	 * This test is equivalent to -gavril -mmd -dfs -del_ch options
	 */
	int i = 0;

	create_WIS_graph(info, wmg);

	info->gavril = true;
	info->elim_order_type = GD_MUL_MIN_DEGREE;
	info->free_children = true;

	T->info = info;

	create_tree_decomposition(T->info, wmg, &T);

	EXPECT_EQ(67, T->width)
		;
	EXPECT_EQ(11, T->num_leafs)
		;
	EXPECT_EQ(39, T->num_tree_nodes)
		;

	EXPECT_EQ(0, T->info->opt_obj)
		;

	doDP(T);

	EXPECT_EQ(16, T->info->opt_obj)
		;

	T->free_table(T->root_node);
}
