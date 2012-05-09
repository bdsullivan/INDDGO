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
#include "GraphEOUtil.h"
#include <gtest/gtest.h>

using namespace std;
class GraphEOUtilTest: public testing::Test
{
public:
	Graph::GraphCreatorFile creator;
	Graph::GraphEOUtil eoutil;
	Graph::MutableGraph *mg;

	virtual void SetUp()
	{
		//SetUp is called before every test
		LOG_INIT("test.log", NULL, 0);
		creator.set_file_name("../data/1dc.128.txt");
		creator.set_graph_type("DIMACS");
		mg = creator.create_mutable_graph();
	}

	virtual void TearDown()
	{
		//TearDown is called before every test
		LOG_CLOSE();
	}
};

TEST_F(GraphEOUtilTest, testFindMcsOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MCS, false);
	int mcs_output[] =
	{ 0, 64, 1, 96, 32, 16, 8, 4, 2, 3, 112, 120, 124, 7, 15, 65, 48, 97, 80,
			67, 33, 66, 6, 5, 113, 104, 56, 72, 24, 17, 40, 88, 99, 98, 81, 49,
			68, 34, 35, 71, 11, 9, 69, 12, 14, 10, 70, 20, 36, 100, 116, 28,
			60, 121, 76, 108, 92, 84, 52, 13, 18, 114, 105, 73, 19, 25, 44, 57,
			41, 115, 89, 82, 50, 101, 83, 37, 51, 102, 38, 74, 103, 39, 21, 75,
			22, 78, 26, 77, 23, 79, 30, 29, 27, 42, 106, 85, 43, 86, 53, 90,
			58, 46, 45, 54, 122, 117, 107, 118, 110, 109, 87, 94, 93, 91, 31,
			47, 62, 61, 59, 55, 126, 125, 123, 119, 111, 95, 63, 127 };
	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}

}

TEST_F(GraphEOUtilTest, testFindMinDegreeOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//eoutil.find_min_degree_ordering(mg, &ordering);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MIN_DEGREE, false);
	int mcs_output[] =
	{ 127, 0, 1, 63, 7, 96, 120, 2, 95, 15, 124, 4, 24, 99, 111, 80, 119, 123,
			11, 57, 78, 116, 23, 42, 32, 64, 73, 108, 13, 91, 8, 16, 48, 14,
			39, 60, 101, 20, 40, 34, 59, 61, 62, 118, 121, 122, 125, 126, 17,
			18, 33, 36, 65, 66, 68, 72, 81, 82, 97, 98, 3, 5, 6, 9, 10, 12, 19,
			21, 22, 25, 26, 27, 28, 29, 30, 31, 35, 37, 38, 41, 43, 44, 45, 46,
			47, 49, 50, 51, 52, 53, 54, 55, 56, 58, 67, 69, 70, 71, 74, 75, 76,
			77, 79, 83, 84, 85, 86, 87, 88, 89, 90, 92, 93, 94, 100, 102, 103,
			104, 105, 106, 107, 109, 110, 112, 113, 114, 115, 117 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindMulMinDegreeOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	// eoutil.find_mul_min_degree_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MUL_MIN_DEGREE, 0,
			false);
	int mcs_output[] =
	{ 0, 127, 126, 64, 120, 31, 7, 125, 32, 112, 3, 123, 103, 28, 16, 47, 8, 4,
			116, 70, 49, 11, 104, 85, 95, 63, 54, 19, 114, 36, 119, 111, 79,
			113, 88, 67, 26, 107, 87, 93, 68, 66, 65, 9, 6, 5, 2, 1, 110, 109,
			94, 91, 62, 61, 59, 55, 46, 45, 30, 29, 124, 122, 121, 118, 117,
			115, 108, 106, 105, 102, 101, 100, 99, 98, 97, 96, 92, 90, 89, 86,
			84, 83, 82, 81, 80, 78, 77, 76, 75, 74, 73, 72, 71, 69, 60, 58, 57,
			56, 53, 52, 51, 50, 48, 44, 43, 42, 41, 40, 39, 38, 37, 35, 34, 33,
			27, 25, 24, 23, 22, 21, 20, 18, 17, 15, 14, 13, 12, 10 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindPktSortOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//eoutil.find_pkt_sort_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_PKT_SORT, 0, false);
	int mcs_output[] =
	{ 0, 127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 117, 116, 115, 114,
			113, 112, 111, 110, 109, 108, 107, 106, 105, 104, 103, 102, 101,
			100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88, 87, 86, 85,
			84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68,
			67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51,
			50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34,
			33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,
			16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindLexmBfsOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//	eoutil.find_lexm_bfs_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_LEXM_BFS, 0, false);
	int mcs_output[] =
	{ 127, 63, 95, 111, 119, 123, 125, 126, 31, 47, 55, 59, 61, 62, 79, 87, 91,
			93, 94, 103, 107, 109, 110, 115, 117, 118, 121, 122, 124, 15, 23,
			27, 29, 30, 39, 43, 45, 46, 51, 53, 54, 57, 58, 60, 71, 75, 77, 78,
			83, 85, 86, 89, 90, 92, 99, 101, 102, 105, 106, 108, 113, 114, 116,
			120, 7, 11, 13, 14, 19, 21, 22, 25, 26, 28, 35, 37, 38, 41, 42, 44,
			49, 50, 52, 56, 67, 69, 70, 73, 74, 76, 81, 82, 84, 88, 97, 98,
			100, 104, 112, 3, 5, 6, 9, 10, 12, 17, 18, 20, 24, 33, 34, 36, 40,
			48, 65, 66, 68, 72, 80, 96, 1, 2, 4, 8, 16, 32, 64, 0 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindLexpBfsOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	// eoutil.find_lexm_bfs_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_LEXP_BFS, 0, false);
	int mcs_output[] =
	{ 127, 63, 95, 126, 111, 125, 119, 123, 31, 62, 47, 55, 61, 59, 79, 124,
			94, 87, 93, 91, 103, 110, 122, 118, 121, 115, 107, 109, 117, 15,
			30, 23, 29, 27, 39, 60, 46, 58, 54, 43, 45, 51, 57, 53, 71, 78, 75,
			77, 120, 99, 113, 92, 116, 108, 83, 86, 90, 85, 89, 102, 114, 101,
			106, 105, 7, 14, 11, 13, 28, 19, 22, 26, 21, 25, 35, 56, 49, 38,
			37, 44, 52, 42, 50, 41, 67, 70, 69, 112, 97, 98, 76, 74, 73, 88,
			104, 81, 84, 100, 82, 3, 6, 5, 12, 10, 9, 24, 17, 20, 18, 48, 33,
			34, 40, 36, 96, 65, 66, 68, 80, 72, 1, 2, 4, 8, 16, 32, 64, 0 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindMcsmOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//eoutil.find_mcsm_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MCSM, 0, false);
	int mcs_output[] =
	{ 127, 63, 95, 111, 119, 123, 125, 126, 31, 47, 55, 59, 61, 62, 79, 87, 91,
			93, 94, 103, 107, 109, 110, 115, 117, 118, 121, 122, 124, 15, 23,
			27, 29, 30, 39, 43, 45, 46, 51, 53, 54, 57, 58, 60, 71, 75, 77, 78,
			83, 85, 86, 89, 90, 92, 99, 101, 102, 105, 106, 108, 113, 114, 116,
			120, 7, 11, 13, 14, 19, 21, 22, 25, 26, 28, 35, 37, 38, 41, 42, 44,
			49, 50, 52, 56, 67, 69, 70, 73, 74, 76, 81, 82, 84, 88, 97, 98,
			100, 104, 112, 3, 5, 6, 9, 10, 12, 17, 18, 20, 24, 33, 34, 36, 40,
			48, 65, 66, 68, 72, 80, 96, 1, 2, 4, 8, 16, 32, 64, 0 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindMinFillOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//eoutil.find_min_fill_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MIN_FILL, 0, false);
	int mcs_output[] =
	{ 0, 127, 126, 64, 120, 31, 7, 125, 32, 123, 119, 111, 16, 8, 4, 2, 1, 112,
			95, 63, 62, 61, 59, 55, 47, 65, 66, 68, 72, 11, 33, 17, 34, 18, 36,
			40, 24, 20, 12, 94, 93, 91, 87, 79, 110, 118, 109, 117, 115, 107,
			103, 30, 46, 54, 58, 29, 45, 57, 53, 51, 78, 86, 90, 102, 106, 77,
			89, 85, 105, 101, 99, 83, 124, 122, 121, 116, 114, 113, 108, 92,
			60, 104, 100, 98, 97, 96, 88, 84, 82, 81, 80, 76, 75, 74, 73, 71,
			70, 69, 67, 56, 52, 50, 49, 48, 44, 43, 42, 41, 39, 38, 37, 35, 28,
			27, 26, 25, 23, 22, 21, 19, 15, 14, 13, 10, 9, 6, 5, 3 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindBatchMinFillOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//	eoutil.find_batch_min_fill_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_BATCH_MF, 0, false);
	int mcs_output[] =
	{ 0, 127, 126, 64, 120, 31, 7, 125, 32, 123, 16, 119, 8, 111, 4, 2, 1, 112,
			95, 63, 62, 61, 59, 55, 47, 65, 66, 68, 72, 11, 33, 17, 34, 18, 36,
			40, 24, 20, 12, 94, 93, 91, 87, 79, 110, 118, 109, 117, 115, 107,
			103, 30, 46, 54, 58, 29, 45, 57, 53, 51, 78, 86, 90, 102, 106, 77,
			89, 85, 105, 101, 99, 83, 124, 122, 121, 116, 114, 113, 108, 92,
			60, 104, 100, 98, 97, 96, 88, 84, 82, 81, 80, 76, 75, 74, 73, 71,
			70, 69, 67, 56, 52, 50, 49, 48, 44, 43, 42, 41, 39, 38, 37, 35, 28,
			27, 26, 25, 23, 22, 21, 19, 15, 14, 13, 10, 9, 6, 5, 3 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindBetaOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//	eoutil.find_beta_ordering(mg, &ordering);
	eoutil.find_elimination_ordering(mg, &ordering, GD_BETA, false);
	int mcs_output[] =
	{ 0, 127, 1, 63, 64, 126, 3, 7, 15, 31, 96, 112, 120, 124, 2, 32, 95, 125,
			4, 8, 16, 111, 119, 123, 62, 65, 14, 30, 56, 60, 67, 71, 97, 113,
			6, 28, 48, 79, 99, 121, 12, 24, 103, 115, 5, 47, 80, 122, 11, 23,
			104, 116, 9, 33, 55, 61, 66, 72, 94, 118, 17, 59, 68, 110, 19, 27,
			29, 35, 92, 98, 100, 108, 13, 39, 88, 114, 10, 40, 87, 117, 20,
			107, 49, 57, 70, 78, 25, 51, 76, 102, 21, 43, 84, 106, 42, 85, 46,
			58, 69, 81, 22, 26, 44, 52, 75, 83, 101, 105, 18, 36, 91, 109, 34,
			93, 54, 73, 38, 50, 77, 89, 41, 53, 74, 86, 37, 45, 82, 90 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindMetisMmdOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//eoutil.find_metis_mmd_ordering(mg, &ordering);
	eoutil.find_elimination_ordering(mg, &ordering, GD_METIS_MMD, false);
	int mcs_output[] =
	{ 127, 0, 63, 1, 124, 112, 15, 95, 2, 96, 7, 115, 28, 111, 4, 8, 122, 16,
			70, 49, 104, 11, 55, 85, 100, 26, 19, 121, 118, 80, 32, 64, 72, 44,
			65, 119, 123, 125, 126, 62, 61, 59, 94, 110, 103, 87, 91, 93, 107,
			109, 117, 33, 48, 97, 98, 66, 68, 47, 23, 31, 39, 43, 79, 67, 3, 5,
			6, 9, 17, 34, 35, 37, 51, 69, 71, 73, 75, 81, 82, 83, 99, 101, 102,
			105, 106, 113, 114, 116, 120, 20, 10, 12, 13, 14, 18, 21, 22, 24,
			25, 27, 29, 30, 36, 38, 40, 41, 42, 45, 46, 50, 52, 53, 54, 56, 57,
			58, 60, 74, 76, 77, 78, 84, 86, 88, 89, 90, 92, 108 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindMetisNodeNdOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//	eoutil.find_metis_node_nd_ordering(mg, &ordering);
	eoutil.find_elimination_ordering(mg, &ordering, GD_METIS_NODE_ND, false);
	int mcs_output[] =
	{ 0, 96, 3, 24, 80, 5, 20, 6, 18, 9, 10, 12, 17, 64, 32, 8, 1, 2, 4, 16,
			33, 34, 36, 40, 48, 65, 66, 68, 72, 127, 120, 15, 99, 63, 116, 23,
			95, 78, 57, 83, 124, 111, 27, 122, 119, 101, 108, 85, 121, 47, 31,
			55, 58, 60, 54, 59, 61, 62, 118, 123, 125, 126, 43, 45, 51, 53, 89,
			91, 105, 106, 107, 109, 113, 114, 115, 117, 71, 29, 30, 39, 46, 75,
			77, 79, 86, 87, 90, 92, 93, 94, 102, 103, 110, 112, 28, 100, 104,
			69, 73, 38, 74, 84, 11, 76, 52, 22, 82, 98, 35, 14, 50, 37, 42, 67,
			21, 44, 41, 88, 56, 25, 97, 13, 26, 70, 19, 81, 7, 49 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindMetisEdgeNdOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//	eoutil.find_metis_edge_nd_ordering(mg, &ordering);
	eoutil.find_elimination_ordering(mg, &ordering, GD_METIS_EDGE_ND, false);
	int mcs_output[] =
	{ 0, 112, 7, 28, 1, 96, 11, 2, 70, 49, 26, 84, 4, 19, 8, 16, 104, 100, 21,
			14, 13, 22, 88, 37, 25, 35, 38, 41, 42, 44, 50, 52, 56, 67, 69, 73,
			74, 76, 81, 82, 97, 98, 12, 3, 5, 6, 9, 10, 17, 18, 20, 24, 32, 33,
			34, 36, 40, 48, 64, 65, 66, 68, 72, 80, 127, 124, 31, 115, 122, 47,
			107, 79, 91, 55, 59, 87, 103, 126, 125, 119, 61, 62, 63, 93, 94,
			95, 109, 110, 111, 117, 118, 121, 123, 120, 108, 114, 15, 78, 102,
			101, 106, 99, 105, 90, 89, 86, 71, 83, 85, 46, 75, 77, 92, 58, 57,
			51, 53, 54, 39, 23, 43, 45, 30, 116, 27, 29, 113, 60 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindAmdOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//	eoutil.find_amd_ordering(mg, &ordering);
	eoutil.find_elimination_ordering(mg, &ordering, GD_AMD, false);
	int mcs_output[] =
	{ 0, 64, 32, 16, 8, 4, 36, 70, 7, 3, 11, 19, 67, 1, 2, 5, 6, 9, 65, 66, 68,
			49, 120, 112, 116, 104, 114, 113, 88, 28, 54, 26, 85, 103, 127,
			126, 125, 123, 31, 47, 63, 95, 79, 111, 119, 87, 107, 93, 29, 30,
			45, 46, 55, 59, 61, 62, 91, 94, 110, 109, 10, 12, 13, 14, 15, 17,
			18, 21, 22, 23, 24, 25, 27, 33, 34, 35, 37, 38, 39, 40, 41, 42, 43,
			44, 48, 50, 51, 52, 53, 56, 57, 58, 60, 69, 71, 72, 73, 74, 75, 76,
			77, 78, 80, 81, 82, 83, 84, 86, 89, 90, 92, 96, 97, 98, 99, 100,
			101, 102, 105, 106, 108, 115, 117, 118, 121, 122, 124, 20 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testFindMinMaxDegreeOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	//	eoutil.find_minmaxdegree_ordering(mg, &ordering, 0);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MINMAX_DEGREE, 0, false);
	int mcs_output[] =
	{ 0, 127, 63, 64, 3, 15, 112, 124, 32, 95, 12, 16, 99, 111, 5, 23, 104,
			119, 122, 57, 78, 27, 42, 88, 1, 2, 73, 123, 125, 126, 62, 118, 4,
			6, 8, 109, 56, 39, 9, 65, 17, 33, 66, 68, 107, 10, 18, 20, 34, 36,
			30, 31, 47, 55, 71, 79, 87, 91, 94, 103, 110, 102, 7, 11, 13, 14,
			19, 21, 22, 24, 25, 26, 28, 29, 35, 37, 38, 40, 41, 43, 44, 45, 46,
			48, 49, 50, 51, 52, 53, 54, 58, 59, 60, 61, 67, 69, 70, 72, 74, 75,
			76, 77, 80, 81, 82, 83, 84, 85, 86, 89, 90, 92, 93, 96, 97, 98,
			100, 101, 105, 106, 108, 113, 114, 115, 116, 117, 120, 121 };

	int i = 0;
	int n = mg->get_num_nodes();

	for (i = 0; i < n; i++)
	{
		EXPECT_EQ(mcs_output[i], ordering[i])
			;
	}
}

TEST_F(GraphEOUtilTest, testTriangulate)
{
	int width = 0;
	vector<int> ordering(mg->get_num_nodes(), -1);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MCS, false);
	width = eoutil.triangulate(mg, &ordering);
	EXPECT_EQ(73, width)
		;
}

TEST_F(GraphEOUtilTest, testMETISTriangulate)
{
	int width = 0;
	vector<int> ordering(mg->get_num_nodes(), -1);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MCS, false);
	width = eoutil.METIS_triangulate(mg, &ordering);
	EXPECT_EQ(73, width)
		;
}

TEST_F(GraphEOUtilTest, testGetTwLowerBound)
{
	EXPECT_EQ(16, eoutil.get_tw_lower_bound(mg, GD_MCS_LB, 0))
		;
	EXPECT_EQ(15, eoutil.get_tw_lower_bound(mg, GD_MAX_MIN_DEGREE_LB, 0))
		;
}

TEST_F(GraphEOUtilTest, testIsPerfectOrdering)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	eoutil.find_elimination_ordering(mg, &ordering, GD_MCS, false);
	EXPECT_FALSE(eoutil.is_perfect_ordering(mg, &ordering))
		;
	eoutil.triangulate(mg, &ordering);
	EXPECT_TRUE(eoutil.is_perfect_ordering(mg, &ordering))
		;
}

TEST_F(GraphEOUtilTest, testForwardNeighbors)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	list<int> nbrs;
	int x;
	vector<int>::iterator it;

	eoutil.find_elimination_ordering(mg, &ordering, GD_MCS, false);
	eoutil.find_forward_neighbors(mg, ordering.back(), &ordering,
			mg->get_num_nodes() - 1, &nbrs, &x);

	EXPECT_EQ(0, nbrs.size())
		;

	it = find(ordering.begin(), ordering.end(), 25);
	eoutil.find_forward_neighbors(mg, 25, &ordering,
				*it, &nbrs, &x);

}

TEST_F(GraphEOUtilTest, testBackwardNeighbors)
{
	vector<int> ordering(mg->get_num_nodes(), -1);
	list<int> nbrs;
	int x;

	eoutil.find_elimination_ordering(mg, &ordering, GD_MCS, false);
	eoutil.find_backward_neighbors(mg, 25, &ordering, 0, &nbrs,
			&x);
}
