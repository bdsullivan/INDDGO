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

#include "GraphDecomposition.h"

//BDS April 3, 2012 replaced idxmalloc with imalloc for metis 5.0.2 compatbility

namespace Graph
{
	GraphEOUtil::GraphEOUtil()
	{
	}

	GraphEOUtil::~GraphEOUtil()
	{
	}

	int GraphEOUtil::triangulate(MutableGraph *mg, vector<int> *ordering)
	{
		int i, j, w, x, num_fill_edges = 0;
		list<int> neighbors;
		list<int>::iterator k1, k2;
		time_t start, stop;
		if (mg->num_connected_components != 1)
			fatal_error("%s:  Must have 1 component\nmg graph has %d components\n",
			__FUNCTION__, mg->num_connected_components);

		if (!mg->canonical)
			fatal_error("%s:  Requires canonical form\n", __FUNCTION__);

		if ((int) (((((((((((((((((((((((((((ordering->size())))))))))))))))))))))))))))
			!= mg->num_nodes)
			fatal_error("%s:  Ordering appears to be of incorrect size? (%d!=%d)\n",
			__FUNCTION__, ordering->size(), mg->num_nodes);

		int *f = new int[mg->num_nodes];
		int *index = new int[mg->num_nodes];
		int *fwd_nbrs = new int[mg->num_nodes];
		int tree_w = 0;
		mg->key = 1;
		start = clock();
		for (i = 0; i < mg->num_nodes; i++)
		{
			w = ordering->at(i);
			f[w] = w;
			index[w] = i;
			neighbors.clear();
			mg->key++;
			find_backward_neighbors(mg, w, ordering, i, &neighbors, &j);
			fwd_nbrs[w] = mg->nodes[w].nbrs.size() - neighbors.size();
			if (fwd_nbrs[w] > tree_w)
				tree_w = fwd_nbrs[w];

			for (k1 = neighbors.begin(); k1 != neighbors.end(); ++k1)
			{
				x = *k1;
				while (index[x] < i)
				{
					index[x] = i;
					if (mg->adj_vec[x] != mg->key)
					{
						num_fill_edges++;
						mg->nodes[w].nbrs.push_back(x);
						mg->nodes[x].nbrs.push_back(w);
						mg->degree[w]++;
						mg->degree[x]++;
						mg->num_edges++;
						fwd_nbrs[x]++;
						if (fwd_nbrs[x] > tree_w)
							tree_w = fwd_nbrs[x];

					}
					x = f[x];
				}

				if (f[x] == x)
					f[x] = w;

			}
		}

		stop = clock();
		print_message(1, "Found Y-T arrays - %d fill edges found. Width is %d\n",
			num_fill_edges, tree_w);
		delete[] f;
		delete[] index;
		delete[] fwd_nbrs;
		return tree_w;
	}

	int GraphEOUtil::find_backward_neighbors(MutableGraph *mg, int v,
		vector<int> *W, int end_pos, list<int> *neighbors, int *min_pos)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		if (v < 0 || v >= mg->capacity)
			fatal_error(
			"%s:  find_forward_neighbors() called with vertex %d and there are %d connected nodes\n",
			__FUNCTION__, v, mg->capacity);

		if (mg->nodes[v].nbrs.size() == 0 && mg->num_nodes > 1)
			print_message(
			0,
			"%s:  Calling find_forward_neighbors for node %d which has no neighbors!\n",
			__FUNCTION__, v);

		int j, cnt;
		GraphProperties properties;
		int current_key = properties.fill_adj_vec(mg, v);
		*min_pos = INT_MAX;
		neighbors->clear();
		cnt = 0;
		print_message(
			2,
			"Checking %d neighbors\n",
			(int) ((((((((((((((((((((((((((W->size())))))))))))))))))))))))))));
		for (j = 0; j <= end_pos; j++)
		{
			if (mg->adj_vec[W->at(j)] == current_key)
			{
				if (j < *min_pos)
					*min_pos = j;

				neighbors->push_back(W->at(j));
				cnt++;
			}
		}

		return cnt;
	}

	int GraphEOUtil::get_max_min_degree_lower_bound(MutableGraph *mg)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int i, j, min_degree, n;
		int max_recorded_degree = -GD_INFINITY;

		// Create a copy of G - otherwise we would have to restore it at the end
		//Changed from = *this to copy constructor July 19 - BDS
		MutableGraph G(*mg);

		//Make sure the degrees are up to date before we start
		GraphUtil util;
		util.recompute_degrees(&G);

		// Use n to keep track of how many non-removed nodes remain in the graph
		n = G.num_nodes;
		// We need the # of active nodes here and not capacity -
		// also, isolated (degree 0) nodes are fine here
		// since they will just be at the start of the ordering
		while (n != 0)
		{
			print_message(1, "mmd:  There are %d non-removed nodes in the graph\n",
				n);
			// Find the current vertex that has minimum degree
			min_degree = GD_INFINITY;
			j = -1;
			for (i = 0; i < G.capacity; i++)
			{
				if (G.degree[i] < min_degree && G.nodes[i].label != GD_UNDEFINED)
				{
					min_degree = G.degree[i];
					j = i;
				}
			}
			if (min_degree > max_recorded_degree)
				max_recorded_degree = min_degree;
			if (j == -1)
			{
				fatal_error("%s:  Did not find min_degree < GD_INFINITY??\n",
					__FUNCTION__);
			}

			print_message(1, "Removing vertex %d with degree %d\n", j, G.degree[j]);
			print_message(1, "max_recorded_degree is %d\n", max_recorded_degree);

			//decrement the degrees of j's neighbors (the ones still in G) and inflate j's degree.
			//remember that G is a copy, so this is okay.
			list<int>::iterator ii;
			for (ii = G.nodes[j].nbrs.begin(); ii != G.nodes[j].nbrs.end(); ++ii)
			{
				if (G.degree[*ii] != GD_INFINITY)
					G.degree[*ii]--;
			}
			G.degree[j] = GD_INFINITY;
			n--;
		}

		// G's destructor called here automatically
		// We don't care about the degree array since we modified a copy
		return max_recorded_degree;
	}
	int GraphEOUtil::get_mcs_lower_bound(MutableGraph *mg, int start_v)
	{
		vector<int> ordering(mg->num_nodes, -1);
		int lb = find_mcs_ordering(mg, &ordering, start_v);
		return lb;
	}

	int GraphEOUtil::find_mcs_ordering(MutableGraph *mg, vector<int> *ordering,
		int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int i, j, m, x, global_max_deg = -1;
		list<int>::iterator k;
		list<int> neighbors;
		char *mcs_labels;
		if ((int) ((((((((((((((((((((((ordering->size()))))))))))))))))))))))
			< mg->capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), mg->capacity);

		m = start_v;
		mcs_labels = new char[mg->capacity];
		memset(mcs_labels, 0, mg->capacity * sizeof(char));
		mcs_labels[m] = 1;
		ordering->at(mg->num_nodes - 1) = m;
		i = 1;
		while (i < mg->num_nodes)
		{
			int max_deg = -1, max_pos = -1;
			for (j = 0; j < mg->capacity; j++)
			{
				if (mg->nodes[j].label != -1 && mcs_labels[j] == 0)
				{
					neighbors.clear();
					find_forward_neighbors(mg, j, ordering, mg->num_nodes - i,
						&neighbors, &x);
					if ((int) ((((((((((((((((((((((neighbors.size()))))))))))))))))))))))
			> max_deg)
					{
						max_deg =
							(int) ((((((((((((((((((((((neighbors.size()))))))))))))))))))))));
						max_pos = j;
						if (max_deg > global_max_deg)
						{
							global_max_deg = max_deg;
							print_message(1, "New global_max_deg=%d\n",
								global_max_deg);
						}
					}

				}

			}

			mcs_labels[max_pos] = 1;
			ordering->at(mg->num_nodes - 1 - i) = max_pos;
			i++;
		}

		delete[] mcs_labels;
		return global_max_deg;
	}

	int GraphEOUtil::find_forward_neighbors(MutableGraph *mg, int v, vector<int> *W,
		int start_pos, list<int> *neighbors, int *min_pos)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		if (v < 0 || v >= mg->capacity)
			fatal_error(
			"%s:  find_forward_neighbors() called with vertex %d and there are %d connected nodes\n",
			__FUNCTION__, v, mg->capacity);

		if (mg->nodes[v].nbrs.size() == 0 && mg->num_nodes > 1)
			print_message(
			0,
			"%s:  Calling find_forward_neighbors for node %d which has no neighbors!\n",
			__FUNCTION__, v);

		int j, cnt;
		list<int>::iterator k;
		vector<int>::iterator i;
		GraphProperties properties;
		int current_key = properties.fill_adj_vec(mg, v);
		*min_pos = INT_MAX;
		neighbors->clear();
		cnt = 0;
		print_message(2, "Checking %d neighbors\n",
			(int) (((((((((((((((((((((W->size()))))))))))))))))))))));
		for (j = start_pos;
			j < (int) (((((((((((((((((((((W->size()))))))))))))))))))))); j++)
		{
			if (mg->adj_vec[W->at(j)] == current_key)
			{
				if (j < *min_pos)
					*min_pos = j;

				neighbors->push_back(W->at(j));
				cnt++;
			}
		}

		return cnt;
	}

	int GraphEOUtil::find_forward_neighbors_2(MutableGraph *mg, int v,
		vector<int> *W, int start_pos, list<int> *neighbors, int *min_pos)
	{
		vector<short int> frontvals(mg->get_num_nodes(), -1);
		neighbors->clear();

		for (int j = start_pos; j < (int) W->size(); j++)
		{
			frontvals[W->at(j)] = j;
		}

		list<int>::iterator it;
		Node *n = mg->get_node(v);
		list<int> *nbrs = n->get_nbrs_ptr();
		it = nbrs->begin();
		int count = 0;
		*min_pos = INT_MAX;

		for (; it != nbrs->end(); ++it)
		{
			if (frontvals[*it] > 0)
			{
				neighbors->push_back(*it);
				if (frontvals[*it] < *min_pos)
				{
					*min_pos = frontvals[*it];
				}
				count ++;
			}
		}

		return count;
	}

	int GraphEOUtil::METIS_triangulate(MutableGraph *mg, vector<int> *ordering)
	{
#if !HAS_METIS
		fatal_error("Called METIS triangulation without HAS_METIS.\n");
		// Not reached but MSVC complains
		return 0;
#else
		int i, j, k, nvtxs; 
		idx_t maxlnz, maxsub;
		idx_t *xadj, *adjncy;//need to be idx_t for METIS, but int* for use elsewhere in Graph; this is a problem on 64-bit.
		idx_t *perm, *xlnz, *xnzsub, *nzsub, *iperm;
		nvtxs = mg->num_nodes;
		//print_message(0,"nvtxs=num_nodes=%d\n",nvtxs);
		GraphUtil util;
		util.populate_CRS(mg);

		/**
		 *This block is to allow compatibility with 64-bit metis where idx_t is not an int, and replaces the original code below it.
		 */
		int sizexa = mg->get_num_nodes() + 1;
		int sizead = mg->xadj[nvtxs]; //2*(mg->get_num_edges());

		maxsub = 4*sizead;
				
		xadj = new idx_t[sizexa];
		adjncy = new idx_t[sizead];
		for(i = 0; i < sizexa;i++)
		  {
		    xadj[i] = ((mg->xadj)[i])+1;
		  }
		for(i = 0; i < sizead; i++)
		  {
		    adjncy[i] = (mg->adjncy)[i]+1;
		  } 
		// xadj = &(mg->xadj[0]);
		// adjncy = &(mg->adjncy[0]);
		// maxsub = 4 * xadj[nvtxs];
		// k = xadj[nvtxs];
		// for (i = 0; i < k; i++)
		// 	adjncy[i]++;

		// for (i = 0; i < nvtxs + 1; i++)
		// 	xadj[i]++;

		perm = imalloc(nvtxs + 1, (char *) "ComputeFillIn: perm");
		xlnz = imalloc(nvtxs + 1, (char *) "ComputeFillIn: xlnz");
		xnzsub = imalloc(nvtxs + 1, (char *) "ComputeFillIn: xnzsub");
		nzsub = imalloc(maxsub, (char *) "ComputeFillIn: nzsub");
		iperm = imalloc(nvtxs + 1, (char *) "ComputeFillIn: perm");
		memcpy(perm, &(ordering->at(0)), nvtxs * sizeof(int));
		for (i = 0; i < nvtxs; i++)
			iperm[perm[i]] = i;

		for (i = 0; i < nvtxs; i++)
		{
			iperm[i]++;
			perm[i]++;
		}
		if (smbfct(nvtxs, xadj, adjncy, perm, iperm, xlnz, &maxlnz, xnzsub, nzsub,
			&maxsub))
		{
			free(nzsub);

			maxsub = 4 * maxsub;
			nzsub = imalloc(maxsub, (char *) "ComputeFillIn: nzsub");
			if (smbfct(nvtxs, xadj, adjncy, perm, iperm, xlnz, &maxlnz, xnzsub,
				nzsub, &maxsub))
				errexit((char *) "MAXSUB is too small!");
		}
		//print_message(0, "smbfct found %d nonzeros\n", maxlnz);
		//print_message(0, "Clearing graph\n");

		for (i = 0; i < nvtxs; i++)
		{
			xlnz[i]--;
			mg->nodes[i].nbrs.clear();
		}
		mg->num_edges = 0;
		int width = 0;
		for (i = 0; i < nvtxs - 1; i++)
		{
			k = 0;
			for (j = xlnz[i]; j < xlnz[i + 1]; j++)
			{
				//print_message(0, "[%d,%d]: (%d-%d; %d)\n", i,j, perm[i],
				//	perm[nzsub[xnzsub[i] - 1 + j - xlnz[i]] - 1],
				//	perm[nzsub[xnzsub[i] - 1 + j - xlnz[i]] - 1] - 1);
				mg->nodes[perm[i] - 1].nbrs.push_back(
					perm[nzsub[xnzsub[i] - 1 + j - xlnz[i]] - 1] - 1);
				mg->nodes[perm[nzsub[xnzsub[i] - 1 + j - xlnz[i]] - 1] - 1].nbrs.push_back(
					perm[i] - 1);
				mg->num_edges++;
				mg->degree[perm[nzsub[xnzsub[i] - 1 + j - xlnz[i]] - 1] - 1]++;
				mg->degree[perm[i] - 1]++;
				if (iperm[perm[i] - 1]
				< iperm[perm[nzsub[xnzsub[i] - 1 + j - xlnz[i]] - 1] - 1])
					k++;

			}
			if (k > width)
				width = k;

		}

		util.free_CRS(mg);

		delete[] xadj; 
		delete[] adjncy;

		free(perm);
		free(xlnz);
		free(xnzsub);
		free(nzsub);
		free(iperm);
		return width;
#endif
	}

	bool GraphEOUtil::is_perfect_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		if (!mg->canonical)
			fatal_error("%s:  Requires canonical form\n", __FUNCTION__);

		if ((int) (((((((((((((((((((ordering->size())))))))))))))))))))
			!= mg->num_nodes)
			fatal_error("%s:  Ordering appears to be of incorrect size? (%d!=%d)\n",
			__FUNCTION__, ordering->size(), mg->num_nodes);

		int i, v;
		list<int>::iterator jj;
		GraphProperties properties;
		MutableGraph H(*mg);
		for (i = 0; i < H.num_nodes; i++)
		{
			v = ordering->at(i);
			if (properties.is_clique(&H, &(H.nodes[v].nbrs)) == false)
			{
				DEBUG("Neighbors of %d do not form a clique!\n", v);
				for (jj = H.nodes[v].nbrs.begin(); jj != H.nodes[v].nbrs.end();
					++jj)
					DEBUG("%d, ", *jj);

				DEBUG("\n");
				return false;
			}
			H.remove_vertex(v);
		}

		return true;
	}

	int GraphEOUtil::find_mcs_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_low_degree_vertex(mg);
		return find_mcs_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_min_degree_ordering(MutableGraph *mg,
		vector<int> *ordering, int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int i, m, min_deg, min_pos = -1;
		MutableGraph G(*mg);
		if ((int) (((((((((((((((((ordering->size()))))))))))))))))) < G.capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), G.capacity);

		GraphUtil util;
		util.recompute_degrees(mg);
		i = 0;
		ordering->at(i) = start_v;
		G.eliminate_vertex(start_v, NULL, true);
		G.degree[start_v] = INT_MAX;
		int width = 0;
		for (i = 1; i < G.num_nodes; i++)
		{
			// CSG - mg should be ok if num_nodes!=capacity since in the for(m=..) loop
			// below we check for -1 labels and will encounter exactly num_nodes
			// active nodes
			// Find the minimum degree node in the current graph
			// We obviously don't need to search all the nodes here since
			// i have already been eliminated.
			min_deg = INT_MAX;
			for (m = 0; m < G.capacity; m++)
			{
				// Make sure location m is valid
				if (G.nodes[m].label != -1)
				{
					if (G.degree[m] < min_deg)
					{
						min_deg = G.degree[m];
						min_pos = m;
					}
				}

			}

			ordering->at(i) = min_pos;
			// Now "eliminate" node min_pos from the graph
			// Note that mg modifies G.degree[] array!
			if ((int) (G.nodes[min_pos].nbrs.size()) > width)
				width = G.nodes[min_pos].nbrs.size();

			G.eliminate_vertex(min_pos, NULL, true);
			// Artificially inflate degree[min_pos] to infinity - G is a copy so mg abuse is ok
			G.degree[min_pos] = INT_MAX;
		}
		print_message(1, "Width in %s: %d\n", __FUNCTION__, width);
		return 1;
	}
	int GraphEOUtil::find_min_degree_ordering(MutableGraph *mg,
		vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_low_degree_vertex(mg);
		return find_min_degree_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_mul_min_degree_ordering(MutableGraph *mg,
		vector<int> *ordering, int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int i, m, min_deg;

		// Make a copy of this graph
		//Changed from =*mg to copy constructor - BDS July 19
		MutableGraph G(*mg);

		// Make sure ordering[] is big enough for the graph
		if ((int) ordering->size() < G.capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), G.capacity);

		// Make sure the degrees are accurate
		GraphUtil util;

		util.recompute_degrees(&G);

		// Begin with the provided start_v ("best" ordering probably comes
		// by selecting start_v to be a vertex of minimum degree)
		i = 0;
		ordering->at(i) = start_v;
		G.eliminate_vertex(start_v, NULL, true);
		G.degree[start_v] = GD_INFINITY;
		i++;

		int width = 0;
		list<int> min_nodes;
		while (i < G.num_nodes)
		{
			// Find the minimum degree nodes in the current graph

			min_deg = GD_INFINITY;

			min_nodes.clear();
			for (m = 0; m < G.capacity; m++)
			{
				// Make sure location m is valid
				if (G.nodes[m].label != GD_UNDEFINED)
				{
					if (G.degree[m] <= min_deg)
					{
						min_deg = G.degree[m];
						min_nodes.push_front(m);
					}
				}
			}
			list<int>::iterator ii;
			for (ii = min_nodes.begin(); ii != min_nodes.end(); ++ii)
			{
				if (G.degree[*ii] == min_deg)
				{
					ordering->at(i) = *ii;
					i++;
					// CSG - this width may not be accurate?
					if ((int) G.nodes[*ii].nbrs.size() > width)
						width = G.nodes[*ii].nbrs.size();
					G.eliminate_vertex(*ii, NULL, true);
					// Artificially inflate degree[min_pos] to infinity - G is a copy so this abuse is ok
					G.degree[*ii] = GD_INFINITY;
				}
				else
					break;
			}
		}

		print_message(1, "Width in %s: %d\n", __FUNCTION__, width);
		return 1;
	}
	int GraphEOUtil::find_mul_min_degree_ordering(MutableGraph *mg,
		vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_low_degree_vertex(mg);
		return find_mul_min_degree_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_pkt_sort_ordering(MutableGraph *mg, vector<int> *ordering,
		int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		if ((int) (((((((((((((((ordering->size()))))))))))))))) < mg->capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), mg->capacity);

		GraphUtil util;
		util.recompute_degrees(mg);
		int i, j;
		i = 0;
		ordering->at(i) = start_v;
		j = mg->num_nodes - 1;
		for (i = 1; i < mg->num_nodes; i++)
		{
			if (j != start_v)
				ordering->at(i) = j;

			else
			{
				ordering->at(i) = j - 1;
				j--;
			}
			j--;
		}

		return 1;
	}

	int GraphEOUtil::find_pkt_sort_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		int start_v = mg->num_nodes - 1;
		return find_pkt_sort_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_lexm_bfs_ordering(MutableGraph *mg, vector<int> *ordering,
		int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		VertexLabel *labels;
		labels = new VertexLabel[mg->capacity];
		VertexLabel largest_label, current_label;
		bool *S, *t;
		S = new bool[mg->capacity];
		t = new bool[mg->capacity];
		list<int>::iterator ii;
		int i, j, u, v;
		list<int> nodes_to_relabel;
		list<int> path;
		GraphProperties properties;
		for (i = 0; i < mg->capacity; i++)
		{
			if (mg->nodes[i].label != -1)
			{
				S[i] = true;
				labels[i].id = i;
				labels[i].entries.clear();
			}
			else
			{
				S[i] = false;
				labels[i].id = -1;
				labels[i].entries.clear();
			}
		}

		bool started = false;
		for (i = mg->capacity - 1; i >= 0; i--)
		{
			if (mg->nodes[i].label != -1)
			{
				u = -1;
				if (!started)
				{
					started = true;
					v = start_v;
				}
				else
				{
					v = -1;
					largest_label.id = -1;
					largest_label.entries.clear();
					for (j = mg->capacity - 1; j >= 0; j--)
					{
						if (S[j] && largest_label < labels[j])
						{
							largest_label.id = labels[j].id;
							largest_label.entries.clear();
							for (ii = labels[j].entries.begin();
								ii != labels[j].entries.end(); ++ii)
								largest_label.entries.push_back(*ii);

							v = j;
						}
					}

				}

				if (v == -1)
					fatal_error("%s:  Didn't find valid v in lex search??\n",
					__FUNCTION__);

				ordering->at(i) = v;
				S[v] = false;
				int current_key = properties.fill_adj_vec(mg, v);
				nodes_to_relabel.clear();
				for (u = 0; u < mg->capacity; u++)
				{
					if (S[u] == true && mg->nodes[u].label != -1)
					{
						print_message(1, "Finding path from %d to %d\n", u, v);
						for (int k = 0; k < mg->capacity; k++)
						{
							if (mg->nodes[k].label != -1 && S[k] == true
								&& labels[k] < labels[u])
								t[k] = true;

							else
								t[k] = false;

						}
						if (mg->adj_vec[u] == current_key
							|| properties.is_path(mg, u, v, t))
						{
							if (path.size() > 0)
							{
								print_message(
									0,
									"Found path (length %d) from %d to %d through un-numbered, lower-labeled vertices\n",
									path.size(), u, v);
							}
							nodes_to_relabel.push_back(u);
						}
						else
							print_message(1, "Found no path from %d to %d\n", u, v);

					}
				}

				for (ii = nodes_to_relabel.begin(); ii != nodes_to_relabel.end();
					++ii)
					labels[*ii].entries.push_back(i);

			}
		}

		delete[] labels;
		delete[] S;
		delete[] t;
		return 1;
	}

	int GraphEOUtil::find_lexm_bfs_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_low_degree_vertex(mg);
		return find_lexm_bfs_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_lexp_bfs_ordering(MutableGraph *mg, vector<int> *ordering,
		int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		VertexLabel *labels;
		labels = new VertexLabel[mg->capacity];
		VertexLabel largest_label, current_label;
		bool *S;
		S = new bool[mg->capacity];
		list<int>::iterator ii;
		int i, j, u, v;
		list<int> nodes_to_relabel;
		GraphProperties properties;
		for (i = 0; i < mg->capacity; i++)
		{
			// The set S consists of all nodes to be placed into the ordering
			// Once i is in the ordering we set S[i]=false;
			if (mg->nodes[i].label != GD_UNDEFINED)
			{
				S[i] = true;
				labels[i].id = i;
			}
			else
			{
				S[i] = false;
				labels[i].id = GD_UNDEFINED;
			}
		}
		bool started = false;
		for (i = mg->capacity - 1; i >= 0; i--)
		{
			if (mg->nodes[i].label != GD_UNDEFINED)
			{
				u = GD_UNDEFINED;
				if (!started)
				{
					started = true;
					v = start_v;
				}
				else
				{
					// Let v be one of the un-numbered vertices w/ largest label
					v = GD_UNDEFINED;
					largest_label.id = GD_UNDEFINED;
					largest_label.entries.clear();
					for (j = mg->capacity - 1; j >= 0; j--)
					{
						if (S[j] && largest_label < labels[j])
						{
							largest_label.id = labels[j].id;
							largest_label.entries.clear();
							for (ii = labels[j].entries.begin();
								ii != labels[j].entries.end(); ++ii)
								largest_label.entries.push_back(*ii);

							v = j;
						}
					}

				}

				// Node v is selected
				if (v == GD_UNDEFINED)
					fatal_error("%s:  Didn't find valid v in lex search??\n",
					__FUNCTION__);

				ordering->at(i) = v;
				S[v] = false;
				// Now look for all unlabeled vertices u such that there is a path
				// u=x0-x1-...-xk=v where L(x_j)<L(u)
				int current_key = properties.fill_adj_vec(mg, v);
				// Clear the list of nodes whose weight we will increment
				nodes_to_relabel.clear();
				for (u = 0; u < mg->capacity; u++)
				{
					if (S[u] == true && mg->nodes[u].label != GD_UNDEFINED)
					{
						if (mg->adj_vec[u] == current_key)
						{
							// We will augment the label of u
							nodes_to_relabel.push_back(u);
						}
					}

				}

				// Now augment the labels
				for (ii = nodes_to_relabel.begin(); ii != nodes_to_relabel.end();
					++ii)
					// Note that since i goes from n->0, we keep the labels in
					// sorted decreasing order by pushing back
					labels[*ii].entries.push_back(i);

			}
		}
		delete[] labels;
		delete[] S;
		return 1;
	}
	int GraphEOUtil::find_lexp_bfs_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_low_degree_vertex(mg);
		return find_lexp_bfs_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_mcsm_ordering(MutableGraph *mg, vector<int> *ordering,
		int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int *mcs_weight;
		bool *S, *t;
		list<int>::iterator ii;
		list<int> path;
		list<int> nodes_to_increment;
		GraphProperties properties;
		S = new bool[mg->capacity];
		t = new bool[mg->capacity];
		mcs_weight = new int[mg->capacity];
		int i, j, u, v;
		for (i = 0; i < mg->capacity; i++)
		{
			S[i] = true;
			mcs_weight[i] = 0;
			if (mg->nodes[i].label == GD_UNDEFINED)
			{
				mcs_weight[i] = -GD_INFINITY;
				S[i] = false;
			}
		}
		bool started = false;
		int *perm = new int[mg->capacity];
		for (i = 0; i < mg->capacity; i++)
			perm[i] = i;

		random_permutation(perm, mg->capacity);
		for (i = 0; i <= mg->capacity - 1; i++)
		{
			if (mg->nodes[i].label == GD_UNDEFINED)
				// Treat mg as fatal for now
				fatal_error("Undefined label in MCSM !\n");

			v = GD_UNDEFINED;
			// Find a node in S with largest mcs_weight
			if (!started)
			{
				v = start_v;
				started = true;
			}
			else
			{
				int max_val = -1;
				for (int l = 0; l < mg->capacity; l++)
				{
					j = l;
					print_message(10, "MCS_M: j=%d\n", j);
					// There could be numerous tiebreaking strategies here
					if (mg->nodes[j].label != GD_UNDEFINED
						&& mcs_weight[j] >= max_val && S[j] == true)
					{
						max_val = mcs_weight[j];
						v = j;
					}
				}

			}

			if (v == GD_UNDEFINED)
				fatal_error("%s:  Didn't find valid v in mcs??\n", __FUNCTION__);

			// Add v to the ordering
			ordering->at(mg->num_nodes - 1 - i) = v;
			// Remove v from S
			S[v] = false;
			print_message(
				1,
				"At step %d, the vertex with max mcs_weight is %d (mcs_weight=%d)\n",
				i, v, mcs_weight[v]);
			// Find all unordered nodes u such that there is a path from u to v:
			// u=x0-x1-x2-...-xk=v where mcs_weight[x_j] < mcs_weight[u] for
			// j=1,2,...k-1
			// Create the adj_vec for v
			int current_key = properties.fill_adj_vec(mg, v);
			// Clear the list of nodes whose mcs_weight we will increment
			nodes_to_increment.clear();
			for (u = 0; u < mg->capacity; u++)
			{
				if (S[u] == true && mg->nodes[u].label != GD_UNDEFINED)
				{
					// u is a valid node index and not yet in the ordering
					// Try and find a path from u to v that traverses
					// only nodes with lower mcs_weight values
					print_message(1,
						"Searching for low weight path from %d to %d\n", u, v);
					// Construct t[] to be the vector where t[m]=true for all vertices m with lower weight than u
					// We could be much more efficient here - maintain a sorted list of the weights
					// and then just know where in the list we need to start for a particular u, etc.
					for (int m = 0; m < mg->capacity; m++)
					{
						if (mg->nodes[m].label != GD_UNDEFINED && S[m] == true
							&& mcs_weight[m] < mcs_weight[u])
							t[m] = true;

						else
							t[m] = false;

					}
					if (mg->adj_vec[u] == current_key
						|| properties.is_path(mg, u, v, t))
					{
						print_message(
							1,
							"Found path from %d[%d] to %d[%d] through un-numbered, lower-weighted vertices\n",
							u, mcs_weight[u], v, mcs_weight[v]);
						// We will increment the mcs_weight of u
						nodes_to_increment.push_back(u);
					}
					else
						print_message(1, "Found no path from %d to %d\n", u, v);

				}
			}

			print_message(1,
				"Incrementing the weights of %d vertices (%d possible)\n",
				nodes_to_increment.size(), mg->capacity - i);
			// Now increment the weights
			for (ii = nodes_to_increment.begin(); ii != nodes_to_increment.end();
				++ii)
				mcs_weight[*ii]++;

			for (u = 0; u < mg->capacity; u++)
				print_message(1, "Weight[%d]=%d\n", u, mcs_weight[u]);

		}
		delete[] S;
		delete[] t;
		delete[] mcs_weight;
		delete[] perm;
		return 1;
	}
	int GraphEOUtil::find_mcsm_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_high_degree_vertex(mg);
		return find_mcsm_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_min_fill_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_low_degree_vertex(mg);
		return find_min_fill_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_min_fill_ordering(MutableGraph *mg, vector<int> *ordering,
		int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int i, m, min_fill, min_pos = -1, my_num_edges;
		int u;
		list<int>::iterator nbr_it, twohop_it;
		GraphUtil util;
		int *fill = new int[mg->num_nodes];
		for (i = 0; i < mg->num_nodes; i++)
			fill[i] = GD_UNDEFINED;
		MutableGraph G(*mg);
		if ((int) ((((((((((ordering->size())))))))))) < G.capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), G.capacity);

		util.recompute_degrees(&G);
		i = 0;
		ordering->at(i) = start_v;
		G.eliminate_vertex(start_v, NULL, true);
		fill[start_v] = GD_INFINITY;
		int h, j;
		int *vn = new int[G.num_nodes];
		list<int>::iterator it;
		for (i = 1; i < G.num_nodes; i++)
		{
			// Find the minimum fill-in node in the current graph, computing fill where needed.
			min_fill = GD_INFINITY;
			//loop over the fill array
			for (m = 0; m < G.num_nodes; m++)
			{
				//skip nodes that have already been ordered
				if (fill[m] != GD_INFINITY)
				{
					if (fill[m] == GD_UNDEFINED)
					{
						// Count the number of edges among the neighbors
						for (j = 0; j < G.num_nodes; j++)
							vn[j] = 0;

						it = G.nodes[m].nbrs.begin();
						while (it != G.nodes[m].nbrs.end())
						{
							vn[*it] = 1;
							++it;
						}
						my_num_edges = 0;
						for (j = 0; j < G.num_nodes; j++)
						{
							if (vn[j] == 1)
							{
								for (nbr_it = G.nodes[j].nbrs.begin();
									nbr_it != G.nodes[j].nbrs.end(); ++nbr_it)
								{
									if (vn[*nbr_it] == 1 && j <= *nbr_it)
										// mg is a common neighbor - if we have self-loops mg could cause a problem.
										my_num_edges++;

								}
							}

						}

						h = G.degree[m];
						fill[m] = (h * (h - 1)) / 2 - my_num_edges;
					}

					//check for a new minimum - change <= to < to get first one encountered.
					if (fill[m] <= min_fill)
					{
						min_pos = m;
						min_fill = fill[m];
					}
				} //end of if statement checking if node has been eliminated

			} //end of loop over fill array

			ordering->at(i) = min_pos;
			// set fill to undefined for all vertices within 2 hops of v
			for (nbr_it = G.nodes[min_pos].nbrs.begin();
				nbr_it != G.nodes[min_pos].nbrs.end(); ++nbr_it)
			{
				u = *nbr_it;
				for (twohop_it = G.nodes[u].nbrs.begin();
					twohop_it != G.nodes[u].nbrs.end(); ++twohop_it)
				{
					fill[*twohop_it] = GD_UNDEFINED;
				}
			}

			// Now "eliminate" node min_pos from the graph
			G.eliminate_vertex(min_pos, NULL, true);
			//set the fill array of the eliminated vertex to GD_INFINITY
			fill[min_pos] = GD_INFINITY;
			//continue to next step in finding ordering
		}
		delete[] vn;
		delete[] fill;
		return 1;
	}
	int GraphEOUtil::find_batch_min_fill_ordering(MutableGraph *mg,
		vector<int> *ordering, int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int i, m, min_fill, min_pos = -1, my_num_edges;
		int u;
		list<int>::iterator nbr_it, twohop_it;
		int *fill = new int[mg->num_nodes];
		for (i = 0; i < mg->num_nodes; i++)
			fill[i] = GD_UNDEFINED;
		GraphUtil util;
		list<int> min_list;
		MutableGraph G(*mg);
		if ((int) (((((((((ordering->size()))))))))) < G.capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), G.capacity);

		util.recompute_degrees(&G);
		i = 0;
		ordering->at(i) = start_v;
		G.eliminate_vertex(start_v, NULL, true);
		fill[start_v] = GD_INFINITY;
		int h, j;
		int *vn = new int[G.num_nodes];
		list<int>::iterator it;
		int eliminated;
		for (i = 1; i < G.num_nodes; i++)
		{
			// Find the minimum fill-in node in the current graph, computing fill where needed.
			min_fill = GD_INFINITY;
			min_list.clear(); //clear all the old vertices from the list
			//loop over the fill array
			for (m = 0; m < G.num_nodes; m++)
			{
				//skip nodes that have already been ordered
				if (fill[m] != GD_INFINITY)
				{
					if (fill[m] == GD_UNDEFINED)
					{
						// Count the number of edges among the neighbors
						for (j = 0; j < G.num_nodes; j++)
							vn[j] = 0;

						it = G.nodes[m].nbrs.begin();
						while (it != G.nodes[m].nbrs.end())
						{
							vn[*it] = 1;
							++it;
						}
						my_num_edges = 0;
						for (j = 0; j < G.num_nodes; j++)
						{
							if (vn[j] == 1)
							{
								for (nbr_it = G.nodes[j].nbrs.begin();
									nbr_it != G.nodes[j].nbrs.end(); ++nbr_it)
								{
									if (vn[*nbr_it] == 1 && j <= *nbr_it)
										// mg is a common neighbor - if we have self-loops mg could cause a problem.
										my_num_edges++;

								}
							}

						}

						h = G.degree[m];
						fill[m] = (h * (h - 1)) / 2 - my_num_edges;
					}

					//check for a new minimum
					if (fill[m] <= min_fill)
					{
						//edited to clear the list so we can push from the front or back to reverse
						//order vertices are being considered in - Mar 29
						if (fill[m] < min_fill)
							min_list.clear();

						min_list.push_front(m);
						min_fill = fill[m];
					}
				} //end of if statement checking if node has been eliminated

			} //end of loop over fill array

			print(10, min_list);
			eliminated = 0;
			while (min_list.size() > 0 && fill[min_list.front()] == min_fill)
			{
				eliminated++;
				min_pos = min_list.front();
				min_list.pop_front();
				ordering->at(i) = min_pos;
				// set fill to undefined for all vertices within 2 hops of v
				for (nbr_it = G.nodes[min_pos].nbrs.begin();
					nbr_it != G.nodes[min_pos].nbrs.end(); ++nbr_it)
				{
					u = *nbr_it;
					for (twohop_it = G.nodes[u].nbrs.begin();
						twohop_it != G.nodes[u].nbrs.end(); ++twohop_it)
					{
						fill[*twohop_it] = GD_UNDEFINED;
					}
				}

				i++; // the next run of the loop will fill in the next part of the order
				// Now "eliminate" node min_pos from the graph
				G.eliminate_vertex(min_pos, NULL, true);
				//set the fill array of the eliminated vertex to GD_INFINITY
				fill[min_pos] = GD_INFINITY;
			}

			print_message(1, "Eliminated %d vertices in mg batch\n", eliminated);
			i--; // correct the index term
			//continue to next step in finding ordering
		}
		delete[] vn;
		delete[] fill;
		return 1;
	}
	int GraphEOUtil::find_batch_min_fill_ordering(MutableGraph *mg,
		vector<int> *ordering)
	{
		GraphUtil util;
		int start_v = util.get_random_low_degree_vertex(mg);
		return find_batch_min_fill_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::find_beta_ordering(MutableGraph *mg, vector<int> *ordering)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int m, my_num_edges;
		list<int>::iterator nbr_it, twohop_it;
		list<int_int> fill;
		list<int_int> ties;
		GraphUtil util;
		if ((int) ((((((((ordering->size())))))))) < mg->capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), mg->capacity);

		util.recompute_degrees(mg);
		int *vn = new int[mg->num_nodes];
		int h, j;
		list<int>::iterator it;
		for (m = 0; m < mg->num_nodes; m++)
		{
			int_int curr;
			for (j = 0; j < mg->num_nodes; j++)
				vn[j] = 0;

			it = mg->nodes[m].nbrs.begin();
			while (it != mg->nodes[m].nbrs.end())
			{
				vn[*it] = 1;
				++it;
			}
			my_num_edges = 0;
			for (j = 0; j < mg->num_nodes; j++)
			{
				if (vn[j] == 1)
				{
					for (nbr_it = mg->nodes[j].nbrs.begin();
						nbr_it != mg->nodes[j].nbrs.end(); ++nbr_it)
					{
						if (vn[*nbr_it] == 1 && j <= *nbr_it)
							my_num_edges++;

					}
				}

			}

			h = mg->degree[m];
			curr.p1 = (h * (h - 1)) / 2 - my_num_edges;
			curr.p2 = m;
			fill.push_front(curr);
		}

		fill.sort();
		int i = 0;
		int min;
		list<int_int>::iterator tit;
		while (fill.size() > 0)
		{
			ties.clear();
			min = fill.front().p1;
			while (fill.size() > 0 && fill.front().p1 == min)
			{
				ties.push_front(fill.front());
				ties.front().p1 = mg->degree[ties.front().p2];
				fill.pop_front();
			}
			ties.sort();
			for (tit = ties.begin(); tit != ties.end(); ++tit)
			{
				ordering->at(i) = (*tit).p2;
				i++;
			}
		}

		delete[] vn;
		return 1;
	}

	int GraphEOUtil::find_metis_mmd_ordering(MutableGraph *mg,
		vector<int> *ordering)
	{
#if !HAS_METIS
		fatal_error("Called METIS MMD without HAS_METIS.\n");
		// Not reached but MSVC complains
		return 0;
#else
		int i;
		GraphUtil util;
		util.populate_CRS(mg);
		//BDS April 3, 2012 switched from GraphType, ControlType with metis4->metis5 
		graph_t graph;
		ctrl_t *ctrl;

		graph.nvtxs = mg->get_num_nodes();
		graph.nedges = mg->get_num_edges();
		graph.label = imalloc(graph.nvtxs, (char *) "TD Graph: graph label");
		for (i = 0; i < graph.nvtxs; i++)
			graph.label[i] = i;

		/*
		 * for 64-bit compatibility with METIS
		 * Concern: in triangulation, we increment all vertex numbers by 1 in 
		 * xadj and adjncy; this wasn't being done in this function. Are we getting correct results in both?
		 */
		int sizexa = mg->get_num_nodes() + 1;
		int sizead = 2*(mg->get_num_edges());
		idx_t *xadj = new idx_t[sizexa];
		idx_t *adjncy = new idx_t[sizead];
		for(i = 0; i < sizexa;i++)
		  {
		    xadj[i] = (mg->xadj)[i];//++;
		  }
		for(i = 0; i < sizead; i++)
		  {
		    adjncy[i] = (mg->adjncy)[i];//++;
		  } 

		graph.adjncy = adjncy;//&(mg->adjncy[0]);
		graph.xadj = xadj;//&(mg->xadj[0]);
		//workspace mem alloc uses graph->ncon, which we were not setting before
		graph.ncon = 1;


		//cannot leave ctrl completely uninitialized
		ctrl = SetupCtrl(METIS_OP_OMETIS, NULL, 1, 3, NULL, NULL);
		/* allocate workspace memory */
		AllocateWorkSpace(ctrl, &graph);

		vector<idx_t> metis_order(mg->get_num_nodes(), GD_UNDEFINED);
		MMDOrder(ctrl, &graph, &(metis_order[0]), graph.nvtxs);
		util.free_CRS(mg);
		for (i = 0; i < mg->get_num_nodes(); i++)
			ordering->at(metis_order[i]) = i;

		delete[] xadj; 
		delete[] adjncy; 

		return 1;
#endif
	}



	int GraphEOUtil::find_metis_node_nd_ordering(MutableGraph *mg,
		vector<int> *ordering)
	{
#if !HAS_METIS
		fatal_error("Called METIS node_nd without HAS_METIS.\n");
		// Not reached but MSVC complains
		return 0;
#else
		GraphUtil util;
		util.populate_CRS(mg);
		vector<idx_t> metis_order(mg->get_num_nodes(), GD_UNDEFINED);
		idx_t nvtxs = (idx_t)mg->num_nodes;
		int numflag = 0;

		//int options[10];
		//options[0] = 0;
		//METIS_NodeND(&nvtxs, &(mg->xadj[0]), &(mg->adjncy[0]), &numflag, options,
		//		&(ordering->at(0)), &(metis_order[0]));

		//BDS 04/03/11 - testing with Metis5.0.2 numflag is now set within ctrl structure
		//via options array, and default is 0, so no changes. New NULL parameter 
		// is for vertex weights.

		//WARNING -- 
		//Tried new default options (NULL options array) - does not reproduce old behaviour. 
		//04/04 These options were chosen to (1) try and match what was being set 
		//as default for ONMETIS before and (2) satisfy the checkParams function
		//in 5.0.2. They do NOT currently exactly replicate old behaviour.
		idx_t options[METIS_NOPTIONS];
		METIS_SetDefaultOptions(options);
		options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_NODE; 
		options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; 
		options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE; 
		options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP1SIDED; 
		options[METIS_OPTION_DBGLVL] = 0; 
		options[METIS_OPTION_COMPRESS] = 1;
		options[METIS_OPTION_PFACTOR] = -1; 
		options[METIS_OPTION_NSEPS] = 1; 
		options[METIS_OPTION_SEED] = -1; 
		options[METIS_OPTION_CCORDER] = 0; 
		options[METIS_OPTION_NITER] = 10;

		/*
		 * for 64-bit compatibility with METIS
		 * Concern: in triangulation, we increment all vertex numbers by 1 in 
		 * xadj and adjncy; this wasn't being done in this function. Are we getting correct results in both?
		 */
		int i;
		int sizexa = mg->get_num_nodes() + 1;
		int sizead = 2*(mg->get_num_edges());
		idx_t *xadj = new idx_t[sizexa];
		idx_t *adjncy = new idx_t[sizead];
		for(i = 0; i < sizexa;i++)
		  {
		    xadj[i] = (mg->xadj)[i];//++;
		  }
		for(i = 0; i < sizead; i++)
		  {
		    adjncy[i] = (mg->adjncy)[i];//++;
		  } 
		vector<idx_t> my_order(mg->get_num_nodes(), GD_UNDEFINED);
		METIS_NodeND(&nvtxs, xadj, adjncy, NULL, options,
			     &(my_order[0]), &(metis_order[0]));
		util.free_CRS(mg);

		for (i = 0; i < mg->get_num_nodes(); i++)
		  ordering->at(i) = (int) my_order[i];

		delete[] xadj; 
		delete[] adjncy; 
		
		return 1;

#endif
	}


	//DEPRECATED in Metis 5.0.2 as per discussion with G. Karypis. 
	//Removed from INDDGO Aug 17, 2012.
	//     int GraphEOUtil::find_metis_edge_nd_ordering(MutableGraph *mg,
	//                                                  vector<int> *ordering)
	//     {
	// #if !HAS_METIS
	//       fatal_error("Called METIS edge_nd without HAS_METIS.\n");
	// #else
	//         GraphUtil util;
	//         util.populate_CRS(mg);
	//         vector<int> metis_order(mg->get_num_nodes(), GD_UNDEFINED);
	//         int nvtxs = mg->num_nodes;
	//         int numflag = 0;
	//         //	int options[10];
	//         //options[0] = 0;
	//         //METIS_EdgeND(&nvtxs, &(mg->xadj[0]), &(mg->adjncy[0]), &numflag, options,
	//         //		&(ordering->at(0)), &(metis_order[0]));

	//         //BDS 04/03/11 - testing with Metis5.0.2. numflag is now set within ctrl structure
	//         //via options array, and default is 0, so no changes. New NULL parameter 
	//         // is for vertex weights.

	//         //WARNING -- 
	//         //4/04 These options were chosen to (1) try and match what was being set 
	//         //as default for OEMETIS before and (2) satisfy the checkParams function
	//         //in 5.0.2. They do NOT currently exactly replicate old behaviour.
	//         int options[METIS_NOPTIONS];
	//         options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_NODE; 
	//         options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; 
	//         options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_EDGE; 
	//         options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP1SIDED; 
	//         options[METIS_OPTION_DBGLVL] = 0; 
	//         options[METIS_OPTION_NSEPS] = 1;
	//         options[METIS_OPTION_SEED] = -1; 
	//         options[METIS_OPTION_CCORDER] = 0; 
	//         options[METIS_OPTION_NITER] = 10; 
	//         options[METIS_OPTION_PFACTOR] = 0; 
	//         options[METIS_OPTION_COMPRESS] = 1;
	//         METIS_NodeND(&nvtxs, &(mg->xadj[0]), &(mg->adjncy[0]), NULL, options,
	//                      &(ordering->at(0)), &(metis_order[0]));
	//         util.free_CRS(mg);
	//         return 1;
	// #endif
	//     }

	int GraphEOUtil::find_amd_ordering(MutableGraph *mg, vector<int> *ordering)
	{
#if !HAS_SUITESPARSE
		return 0;
#else
		GraphUtil util;
		util.populate_CRS(mg);
		amd_order(mg->num_nodes, &(mg->xadj[0]), &(mg->adjncy[0]), &ordering->at(0),
			(double *) NULL, (double *) NULL);
		util.free_CRS(mg);
		return 1;
#endif
	}

	int GraphEOUtil::find_minmaxdegree_ordering(MutableGraph *mg,
		vector<int> *ordering, int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int i, m, min_maxd, min_pos = -1;
		int u;
		list<int>::iterator nbr_it;
		GraphUtil util;
		GraphProperties properties;
		int *maxd = new int[mg->num_nodes];
		for (i = 0; i < mg->num_nodes; i++)
			maxd[i] = GD_UNDEFINED;
		MutableGraph G(*mg);
		if ((int) ((((ordering->size())))) < G.capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), G.capacity);

		util.recompute_degrees(&G);
		i = 0;
		if (start_v != -1)
		{
			ordering->at(i) = start_v;
			G.eliminate_vertex(start_v, NULL, true);
			// set maxd to infinity for v (it's been eliminated)
			maxd[start_v] = GD_INFINITY;
			i++;
		}
		int delta;
		list<int>::iterator it, jt;
		for (; i < G.num_nodes; i++)
		{
			// Find the minimum maxd node in the current graph, computing fill where needed.
			min_maxd = GD_INFINITY;
			//loop over the maxd array
			for (m = 0; m < G.num_nodes; m++)
			{
				//skip nodes that have already been ordered
				if (maxd[m] != GD_INFINITY)
				{
					if (maxd[m] == GD_UNDEFINED)
					{
						maxd[m] = 0;
						//compute the actual max degree when m is eliminated
						//create a vector of which vertices are neighbors
						properties.fill_adj_vec(&G, m); // fills G.adj_vec
						//loop over the neighbors of m, counting their number of neighbors in N(m)
						//their new degree is their current degree + N(m) - currMnbrs
						for (it = G.nodes[m].nbrs.begin();
							it != G.nodes[m].nbrs.end(); ++it)
						{
							delta = G.get_degree(m);
							for (jt = G.nodes[*it].nbrs.begin();
								jt != G.nodes[*it].nbrs.end(); ++jt)
							{
								if (G.adj_vec[*jt])
									//m is also adjacent to mg neighbor of *it
									delta--; //we gain one less edge from elimination

							}
							if (maxd[m] < G.get_degree(*it) + delta)
								//we found a new worst case vertex in G eliminate m
								maxd[m] = G.get_degree(*it) + delta;

						}
					}

					//check for a new minimum - if there's  a tie, let min degree decide.
					if (maxd[m] < min_maxd
						|| (maxd[m] == min_maxd
						&& G.get_degree(m) < G.get_degree(min_pos)))
					{
						min_pos = m;
						min_maxd = maxd[m];
					}
				} //end of if statement checking if node has been eliminated

			} //end of loop over fill array

			ordering->at(i) = min_pos;
			// set maxd to undefined for all vertices within 2 hops of v
			for (nbr_it = G.nodes[min_pos].nbrs.begin();
				nbr_it != G.nodes[min_pos].nbrs.end(); ++nbr_it)
			{
				u = *nbr_it;
				for (jt = G.nodes[u].nbrs.begin(); jt != G.nodes[u].nbrs.end();
					++jt)
				{
					maxd[*jt] = GD_UNDEFINED;
				}
			}

			// Now "eliminate" node min_pos from the graph
			G.eliminate_vertex(min_pos, NULL, true);
			//set the fill array of the eliminated vertex to GD_INFINITY
			maxd[min_pos] = GD_INFINITY;
			//continue to next step in finding ordering
		}
		delete[] maxd;
		return 1;
	}
	int GraphEOUtil::find_minmaxdegree_ordering(MutableGraph *mg,
		vector<int> *ordering)
	{
		int start_v = -1;
		return find_minmaxdegree_ordering(mg, ordering, start_v);
	}

	int GraphEOUtil::get_tree_width(MutableGraph *mg, vector<int> *ordering)
	{
		int width = 0;
		int i, k;
		vector<int> inv_order(ordering->size());
		for (i = 0; i < (int) (((ordering->size()))); i++)
			inv_order[ordering->at(i)] = i;

		for (i = 0; i < mg->capacity; i++)
		{
			// Note that degree[i] must be greater than the current width in order for mg
			// node to possibly have more than width higher #'ed nbrs
			if (mg->nodes[i].label != GD_UNDEFINED && mg->degree[i] > width)
			{
				k = 0;
				for (list<int>::iterator ii = mg->nodes[i].nbrs.begin();
					ii != mg->nodes[i].nbrs.end(); ++ii)
				{
					if (inv_order[*ii] > inv_order[i])
						k++;
				}
				if (k > width)
					width = k;
			}
		}
		return width;
	}
	int GraphEOUtil::get_tw_lower_bound(MutableGraph *mg, int algorithm,
		int start_v)
	{
		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		int lb = GD_UNDEFINED;
		switch (algorithm)
		{
		case GD_MAX_MIN_DEGREE_LB:
			// Ignoring start_v here - is there a better way to handle mg??
			lb = get_max_min_degree_lower_bound(mg);
			break;

		case GD_MCS_LB:
			lb = get_mcs_lower_bound(mg, start_v);
			break;

		default:
			fatal_error("%s:  Didn't understand provided algorithm of %d\n",
				__FUNCTION__, algorithm);
			break;
		}
		if (lb == GD_UNDEFINED)
			fatal_error("%s:  Didn't find valid lower bound?\n");
		return lb;
	}

	/**
	* Fills the ordering[] array with a permutation of (0,1,2,...,num_nodes-1).  The
	* permutation is computed by using the specified algorithm (such as GD_MIN_DEGREE)
	* and then running the appropriate algorithm.  Assumes sufficient space has been
	* allocated to the ordering[] array.  The elimination ordering algorithm begins with
	* the provided start_v vertex.
	*/
	void GraphEOUtil::find_elimination_ordering(MutableGraph *mg,
		vector<int> *ordering, int algorithm, int start_v, bool triangulate)
	{
		// Valid algorithms: GD_MIN_DEGREE, GD_MCS, GD_MCSM, GD_LEXM_BFS, GD_LEXP_BFS, GD_MIN_FILL, PKT_SORT, GD_BATCH_MF, GD_MINMAX_DEGREE

		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		// Make sure ordering[] is big enough for the graph
		if ((int) ordering->size() < mg->capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), mg->capacity);

		switch (algorithm)
		{
		case GD_MIN_DEGREE:
			this->find_min_degree_ordering(mg, ordering, start_v);
			break;

		case GD_MINMAX_DEGREE:
			this->find_minmaxdegree_ordering(mg, ordering, start_v);
			break;

		case GD_MUL_MIN_DEGREE:
			this->find_mul_min_degree_ordering(mg, ordering, start_v);
			break;

		case GD_MCS:
			this->find_mcs_ordering(mg, ordering, start_v);
			break;

		case GD_MCSM:
			this->find_mcsm_ordering(mg, ordering, start_v);
			break;

		case GD_LEXM_BFS:
			this->find_lexm_bfs_ordering(mg, ordering, start_v);
			break;

		case GD_LEXP_BFS:
			this->find_lexp_bfs_ordering(mg, ordering, start_v);
			break;

		case GD_MIN_FILL:
			this->find_min_fill_ordering(mg, ordering, start_v);
			break;

		case GD_BATCH_MF:
			this->find_batch_min_fill_ordering(mg, ordering, start_v);
			break;

		case GD_PKT_SORT:
			this->find_pkt_sort_ordering(mg, ordering, start_v);
			break;
		default:
			fatal_error("%s:  Didn't understand provided algorithm of %d\n",
				__FUNCTION__, algorithm);
			break;
		}

		// Triangulate the graph if desired
		if (triangulate)
			this->triangulate(mg, ordering);

		return;

	}

	/**
	* Fills the ordering[] array with a permutation of (0,1,2,...,num_nodes-1).  The
	* permutation is computed by using the specified algorithm (such as GD_MIN_DEGREE)
	* and then running the appropriate algorithm.  Assumes sufficient space has been
	* allocated to the ordering[] array.  Routines begin with a randomly chosen low
	* degree vertex.
	*/
	void GraphEOUtil::find_elimination_ordering(MutableGraph *mg,
		vector<int> *ordering, int algorithm, bool triangulate)
	{
		// Valid algorithms: GD_MIN_DEGREE, GD_MCS, GD_MCSM, GD_LEXM_BFS, GD_LEXP_BFS, GD_MIN_FILL, PKT_SORT,
		// GD_BATCH_MF, GD_METIS_MMD, GD_METIS_NODE_ND, GD_BETA, GD_MINMAX_DEGREE
		// Removed GD_METIS_EDGE_ND Aug 17

		if (!mg->canonical)
			print_message(0, "%s:  Graph must be in canonical form\n",
			__FUNCTION__);

		// Make sure ordering[] is big enough for the graph
		if ((int) ordering->size() < mg->capacity)
			fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
			__FUNCTION__, ordering->size(), mg->capacity);

		switch (algorithm)
		{
		case GD_MIN_DEGREE:
			this->find_min_degree_ordering(mg, ordering);
			break;

		case GD_MINMAX_DEGREE:
			this->find_minmaxdegree_ordering(mg, ordering);
			break;

		case GD_MUL_MIN_DEGREE:
			this->find_mul_min_degree_ordering(mg, ordering);
			break;

		case GD_MCS:
			this->find_mcs_ordering(mg, ordering);
			break;

		case GD_MCSM:
			this->find_mcsm_ordering(mg, ordering);
			break;

		case GD_LEXM_BFS:
			this->find_lexm_bfs_ordering(mg, ordering);
			break;

		case GD_LEXP_BFS:
			this->find_lexp_bfs_ordering(mg, ordering);
			break;

		case GD_MIN_FILL:
			this->find_min_fill_ordering(mg, ordering);
			break;

		case GD_BATCH_MF:
			this->find_batch_min_fill_ordering(mg, ordering);
			break;

		case GD_BETA:
			this->find_beta_ordering(mg, ordering);
			break;

		case GD_PKT_SORT:
			this->find_pkt_sort_ordering(mg, ordering);
			break;

		case GD_METIS_MMD:
			this->find_metis_mmd_ordering(mg, ordering);
			break;

		case GD_METIS_NODE_ND:
			this->find_metis_node_nd_ordering(mg, ordering);
			break;

			//case GD_METIS_EDGE_ND:
			//this->find_metis_edge_nd_ordering(mg, ordering);
			//break;

		case GD_AMD:
			this->find_amd_ordering(mg, ordering);
			break;

		default:
			fatal_error("%s:  Didn't understand provided algorithm of %d\n",
				__FUNCTION__, algorithm);
			break;
		}

		// Triangulate the graph if desired
		if (triangulate)
			this->triangulate(mg, ordering);

		return;
	}

#ifdef HAS_PARMETIS
	void GraphEOUtil::parmetis_elimination_ordering(WeightedMutableGraph *mg, vector<int> &orderingout,
		int algorithm, bool triangulate, MPI_Comm comm)
	{
		int ws;
		int wr;

		char eoname[512];
		char eoname_other[512];

		// Get size and rank from the communicator
		MPI_Comm_size(comm, &ws);
		MPI_Comm_rank(comm, &wr);

		WeightedMutableGraph G(*mg);
		double xtime = MPI_Wtime();
		sprintf(eoname, "%s.order.%d", mg->get_input_file().c_str(), ws);
		sprintf(eoname_other, "%s.order_other.%d", mg->get_input_file().c_str(), ws);

		DEBUG("size: %d, rank %d \n", ws, wr);
		int n = G.get_num_nodes();
		int x = n/ws;
		int xm = n%ws;
		int i = 0;
		DEBUG("n: %d x: %d xm: %d \n", n, x, xm);

		vector<int> xadj;
		vector<int> adjncy;

		vector<int> vtxdist(ws + 1, 0);
		vector<int> sizes(2*ws,0);
		vector<int> ordering(x+1, 0);
		vector<int> recvcnt(ws, 0);
		vector<int> displ(ws, 0);

		int numflag = 0;
		int options[10];

		options[0] = 0;
		vtxdist[0] = 0;
		for (i = 1; i <= ws; i++)
		{
			vtxdist[i] = vtxdist[i - 1] + x;
			if (i <= xm)
				vtxdist[i]++;
		}

		// prepareing displacement and receive counts to use with MPI_Gatherv
		for (i = 0; i < ws; i++)
		{
			recvcnt[i] = x;
			if (i < xm)
				recvcnt[i] ++;

			if (i > 0)
				displ[i] += displ[i - 1] + recvcnt[i - 1];
		}

		DEBUG("range: %d, %d\n", vtxdist[wr], vtxdist[wr + 1]);
		int j = 0;
		xadj.push_back(0);
		for (i = vtxdist[wr]; i < vtxdist[wr + 1]; i++)
		{
			Node *no = G.get_node(i);
			list<int> *l = no->get_nbrs_ptr();
			list<int>::iterator it = l->begin();

			for (; it != l->end(); ++it)
			{
				adjncy.push_back(*it);
				j++;
			}
			xadj.push_back(j);
		}

		if (METIS_OK != ParMETIS_SerialNodeND(&vtxdist.front(), &xadj.front(), &adjncy.front(), &numflag, options, &ordering.front(), &sizes.front(), &comm))
		{
			FERROR("error occured while processing parmetis, aborting\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}


		double parmet_time = MPI_Wtime() - xtime;
		n = G.get_num_nodes();
		vector<int> recvbuf(n, 0);
		vector<int> eo(n, 0);


		if (algorithm > 0)
		{
			if (MPI_SUCCESS !=
				MPI_Gatherv((void *)&ordering.front(), recvcnt[wr], MPI_INT,
				(void *)&recvbuf.front(), &recvcnt.front(), &displ.front(), MPI_INT,
				0, comm))
			{
				FERROR("MPI error occured at Gatherv, Abort!\n");
				MPI_Abort(comm, -1);
			}

			if (wr == 0)
			{
				for (int i = 0; i < n; i++)
					eo[recvbuf[i]] = i;
			}
		}
		else
		{

			if (MPI_SUCCESS !=
				MPI_Allgatherv((void *)&ordering.front(), recvcnt[wr], MPI_INT,
				(void *)&recvbuf.front(), &recvcnt.front(), &displ.front(), MPI_INT, comm))
			{
				FERROR("MPI error occured at Gatherv, Abort!\n");
				MPI_Abort(comm, -1);
			}

			for (int i = 0; i < n; i++)
				orderingout[recvbuf[i]] = i;

			return;
		}

		ordering.clear();
		ordering.resize(recvcnt[wr], 0);

		if (MPI_SUCCESS !=
			MPI_Scatterv ((void *)&eo.front(), &recvcnt.front(), &displ.front(), MPI_INT,
			(void *)&ordering.front(), recvcnt[wr], MPI_INT,
			0, comm))
		{
			FERROR("MPI error occured at Scatterv, Abort! \n");
			MPI_Abort(comm, -1);
		}

		GraphCreatorFile gf;
		WeightedMutableGraph *wg;
		GraphEOUtil eoutil;
		GraphProperties prop;
		list<int>members(ordering.begin(), ordering.end());

		wg = gf.create_component(&G, &members, false);
		prop.make_canonical(wg);

		vector<int> ord(recvcnt[wr], 0);
		vector<int> ordsend(recvcnt[wr, 0]);
		double xxtime = MPI_Wtime();
		eoutil.find_elimination_ordering(wg, &ord, algorithm, false);
		DEBUG("eo time : %f\n", MPI_Wtime() - xxtime);

		int sz = recvcnt[wr];

		for (int i = 0; i < sz; i++)
			ordsend[i] = wg->get_node(ord[i])->get_label() - 1;

		recvbuf.assign(n, -1);
		if (MPI_SUCCESS !=
			MPI_Allgatherv((void *)&ordsend.front(), recvcnt[wr], MPI_INT, (void *)&recvbuf.front(), &recvcnt.front(), &displ.front(), MPI_INT, comm))
		{
			FERROR("MPI error occured at Gatherv, Abort!\n");
			MPI_Abort(comm, -1);
		}

		orderingout = recvbuf;

		double p_amd_time = MPI_Wtime() - xtime;

		DEBUG("ordering is written into %s\n", eoname);
		DEBUG("%f,%f\n", parmet_time, p_amd_time);
	}
#endif // HAS_PARMETIS

}

/**
 * Helper functions for creating elimination orderings. Allow reading from file (normal or SCOTCH), or generating internally (with or without start vertex)
 */
void Graph::form_eo(bool read_order, bool scotch, char* ord_file, int elim_order_type, int start_v, MutableGraph *G, vector<int> *ordering)
{
    GraphEOUtil eoutil;
    if(read_order)
      {
  	if(scotch)
  	  read_SCOTCH_ordering_file(ord_file, ordering);      
  	else
  	  read_ordering_file(ord_file, ordering);      
      }
    else
      {
  	if (start_v == GD_UNDEFINED || start_v < 0)
	  eoutil.find_elimination_ordering(G, ordering, elim_order_type, false);
  	else
	  eoutil.find_elimination_ordering(G, ordering, elim_order_type, start_v, false);
      }
  }
  
void Graph::form_eo(bool read_order, bool scotch, char* ord_file, MutableGraph *G,   vector<int> *ordering)
  {
    form_eo(read_order, scotch, ord_file, 0, 0, G, ordering);
  }
  
void Graph::form_eo(int elim_order_type, int start_v, MutableGraph *G, vector<int> *ordering)
  {
    form_eo(false, false, NULL, elim_order_type, start_v, G, ordering);
  }
  
void Graph::form_eo(int elim_order_type, MutableGraph *G, vector<int> *ordering)
  {
    form_eo(false, false, NULL, elim_order_type, GD_UNDEFINED, G, ordering);
  }




