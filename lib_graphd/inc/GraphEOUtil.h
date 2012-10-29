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

#ifndef GRAPHEOUTIL_H_
#define GRAPHEOUTIL_H_
#include "GraphDecomposition.h"
#ifdef HAS_PARMETIS
#include <mpi.h>
#include <parmetis.h>
#endif
using namespace std;

namespace Graph
{

	class GraphEOUtil
	{
	private:

		//short int *frontvals;
		//int fs;
		// Private elimination ordering routines
		// The versions without a start_v find a random "good" vertex (just low degree for now)
		// GD_MCS
		int find_mcs_ordering(MutableGraph *mg, vector<int> *ordering, int start_v);
		int find_mcs_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_MIN_DEGREE
		int find_min_degree_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);
		int find_min_degree_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_MUL_MIN_DEGREE
		int find_mul_min_degree_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);
		int find_mul_min_degree_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_PKT_SORT
		int find_pkt_sort_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);
		int find_pkt_sort_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_LEXM
		int find_lexm_bfs_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);
		int find_lexm_bfs_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_LEXP
		int find_lexp_bfs_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);
		int find_lexp_bfs_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_MCSM
		int
			find_mcsm_ordering(MutableGraph *mg, vector<int> *ordering, int start_v);
		int find_mcsm_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_MIN_FILL
		int find_min_fill_ordering(MutableGraph *mg, vector<int> *ordering);
		int find_min_fill_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);
		// GD_BATCH_FILL
		int find_batch_min_fill_ordering(MutableGraph *mg, vector<int> *ordering);
		int find_batch_min_fill_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);
		// GD_BETA
		int find_beta_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_METIS_MMD
		// We currently don't control any of the parameters (including start_v)
		int find_metis_mmd_ordering(MutableGraph *mg, vector<int> *ordering);
		/// GD_METIS_NODE_ND
		int find_metis_node_nd_ordering(MutableGraph *mg, vector<int> *ordering);
		// GD_METIS_EDGE_ND
		//int find_metis_edge_nd_ordering(MutableGraph *mg, vector<int> *ordering);

		int find_amd_ordering(MutableGraph *mg, vector<int> *ordering);
		//GD_MAXMIN_DEGREE
		int find_minmaxdegree_ordering(MutableGraph *mg, vector<int> *ordering);
		int find_minmaxdegree_ordering(MutableGraph *mg, vector<int> *ordering,
			int start_v);

		// Private lower bound alg's
		int get_max_min_degree_lower_bound(MutableGraph *mg);
		int get_mcs_lower_bound(MutableGraph *mg, int start_v);

	public:
		GraphEOUtil();
		virtual ~GraphEOUtil();

		/**
		* Tests whether the provided elimination ordering is a perfect elimination
		* ordering for the graph.  Given the ordering v1,v2,..., vn, this
		* means that the neighbors of v_1 form a clique
		* (i.e. it is simplicial) and that v_i is always simplicial after v_1, v_2,ma
		* ..., v_{i-1} are removed from the graph.
		*/
		bool is_perfect_ordering(MutableGraph *mg, vector<int> *ordering);

		// Triangulation
		int triangulate(MutableGraph *mg, vector<int> *ordering);
		int METIS_triangulate(MutableGraph *mg, vector<int> *ordering);

		int find_forward_neighbors(MutableGraph *mg, int v, vector<int> *W,
			int start_pos, list<int> *neighbors, int *min_pos);
		int find_backward_neighbors(MutableGraph *mg, int v, vector<int> *W,
			int end_pos, list<int> *neighbors, int *min_pos);

		int find_forward_neighbors_2(MutableGraph *mg, int v,
			vector<int> *W, int start_pos, list<int> *neighbors, int *min_pos);

		// Manually computes the width corresponding to the ordering by finding the largest
		// number of higher-numbered neighbors
		int get_tree_width(MutableGraph *mg, vector<int> *ordering);

		// Lower bounds on treewidth
		int get_tw_lower_bound(MutableGraph *mg, int algorithm, int start_v);

		// Elimination ordering routines
		void find_elimination_ordering(MutableGraph *mg, vector<int> *ordering,
			int algorithm, int start_v, bool triangulate);
		void find_elimination_ordering(MutableGraph *mg, vector<int> *ordering,
			int algorithm, bool triangulate);
#ifdef HAS_PARMETIS
		void parmetis_elimination_ordering(WeightedMutableGraph *mg, vector<int> &ordering,
			int algorithm, bool triangulate, MPI_Comm comm);

#endif // HAS_PARMETIS

	};

}

#endif /* GRAPHEOUTIL_H_ */
