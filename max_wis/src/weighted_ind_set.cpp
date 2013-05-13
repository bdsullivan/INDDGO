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
#include "TreeDecomposition.h"
#include "weighted_ind_set.h"
#include "Debug.h"
#include "DIMACSGraphWriter.h"

void usage(const char *s)
{
	fprintf(stderr, "Usage: %s [options]\n", s);
	fprintf(
		stderr,
		"\t --- Input/Output Options ---\n"
		"\t -f <DIMACS_file> : loads in the file with symmetric adj lists\n"
		"\t -t <tree_infile> will read in the tree from a file\n"
		"\t -w <tree_outfile> will write the tree to a file\n"
		"\t -gviz <gviz_file>: writes a graphviz version of the tree decomposition\n"
		"\t        to gviz_file, and the DIMACS form to gviz_file.dimacs\n"
		"\t -tpreord may be used in conjunction with -w to write the tree so bags \n"
		"\t        are labeled with a preordering (1 is root, i > j means i is \n"
		"\t        farther from root than j\n"
		"\t -eorder <eo_file> writes elimination ordering into specified file.\n"
		"\t -sol <sol_file> : will write the weighted independent set solution to\n"
		"\t                   the given file.  If not given, then solution is\n"
		"\t                   written to <DIMACS_file>.WIS.sol\n"
		"\t -fix_DIMACS <DIMACS_file> will attempt to repair an existing file\n"
		"\t         and write DIMACS_file.fixed\n"	
		"\t --- Decomposition Construction Options ---\n" 
		"\t -superetree : constructs a non-nice TD using CHOLMOD and supernodal\n"
		"\t               etrees (fast; no triangulation)\n"
		"\t -gavril : constructs a non-nice TD using Gavril's algorithm\n"
		"\t -bk : constructs a non-nice TD using BK algorithm\n"
		"\t -pbag : parallel bag generation, non-nice TD using Gavril's algorithm\n"
		"\t -nice : constructs a nice TD\n"
		"\t -make_nice : will take a non-nice tree and niceify it by adding new\n"
		"\t              tree nodes\n"
		"\t -refine_td : will try to reduce the sizes of the bag intersections by \n"
		"\t              post-processing the tree decomposition\n"
		"\t -check : will manually verify the TD (can be slow!)\n"
		"\t --- Elimination Order Options ---\n"
		"\t -lower_bounds : Runs two lower bound heuristics and prints results.\n"
		"\t                 Can be used in combination with other EO options.\n"
		"\t -s start_v : uses start_v as first vertex in ordering\n"
		"\t -mind : generates an elim. ordering using min degree heuristic\n"
		"\t -mmd : generates an elim. ordering using multiple min degree heuristic\n"
		"\t -minf : generates an elim. ordering using min fill heuristic\n"
		"\t -bmf : generates an elim. ordering using a batched min fill heuristic\n"
		"\t -beta : generates an elim. ordering using the beta heuristic\n"
		"\t -metmmd : generates an elim. ordering using METIS mmd heuristic\n"
		"\t -metnnd : generates an elim. ordering using METS node ND heuristic\n"
		"\t -mcsm : generates an elim. ordering using mcsm euristic\n"
		"\t -mcs  : generates an elim. ordering using mcs\n"
		"\t -lexm : generates an elim. ordering using lex-m bfs heuristic\n"
		"\t -pktsort : generates the elim ordering corresponding to partial\n"
		"\t            ktree generation\n"
		"\t -amd : generates the elim. ordering using UFs approximate min\n"
		"\t       degree heuristic (AMD)\n"
		"\t -minmaxd : generates the elim. ordering using the minimum maximum\n"
		"\t           degree heuristic.\n"
		"\t -parmetis : generates an elim. ordering using ParMETIS\n"
		"\t -ord <ordering_file> : reads in an elimination ordering (1-based as\n"
		"\t                        in DIMACS file)\n"
		"\t -sord <scotch ordering file> : reads in an elimination ordering\n"
		"\t                                produced by Scotch\n"
		"\t --- Dynamic Programming Options ---\n"
		"\t -nonniceDP : will use the non-nice Dynamic Programming routines\n"
		"\t -root <root_node> : will set the TD's root node (default is 0)\n"
		"\t                     (Only works for Gavril and BK trees)\n"
		"\t -asp_root : will use the Aspvall algorithm to find a root node\n"
		"\t -child_root : will use All My Children algorithm to find root node\n"
		"\t -pc : will use the child sets to help construct the parent\n\n"			
		"\t -del_ch : will delete a tree node's table once its parent is processed\n"
		"\t           NOTE: you cannot reconstruct the solution if you do this!\n"
		"\t -no_reconstruct : will not reconstruct the solution\n"
		"\t -split_bag : will search for solutions by splitting the bag in half\n"
		"\t              and then merging the two halves\n"
		"\t -async_tbl : will do the DP one child at a time rather than looping\n"
		"\t              over all children at once\n"
		"\t -dfs : will generate a post-order walk using a DFS (BFS is default)\n"
		"\t        DIMACS_file.mem_est\n"
		"\t -max_width <w> : will not attempt to run the DP if the decomposition\n"
		"\t                  has width larger than w.\n"
		"\t --- Miscellaneous Options ---\n"
		"\t -decompose_only : will just generate the tree decomposition and will\n"
		"\t                  *not* run MWIS\n"
		"\t -width : will compute and print width of the specified decomposition\n"
		"\t          and will *not* proceed\n"
		"\t -hist : will print out a histogram of bag sizes prior to running DP\n"
		"\t -mem_est : will estimate the memory usage at each treenode\n"
		"\t -mod <model_file> : writes a GMPL MIP formulation of the problem\n"
		"\t                     to model_file.\n"
		"\t -mip : solves the problem by runing GLPK's glpsol solver\n"
		"\t        via a system() call to glpsol -m model_file\n"
		"\t -v : runs in verbose mode\n"
		"\t -vv : runs in very verbose mode - lots of output about the ind. sets\n"
		"\t -noheader : will not print header for output data to stdout\n"
		);
		

}
;

/**
* Function to quickly get the number of bits set
* in an unsigned long long.
*/
int hamming_weight(unsigned long long xx)
{
	unsigned long long x = xx;
	int ans = 0;
	while (x)
	{
		x &= (x - 1);
		ans++;
	}
	return ans;
}

/**
* Function to compare two int_bigint's based on their hamming weight.
* The idea is that when we check for independent sets, we should try to and with
* the densest masks first.
*/
bool hamming_compare(const int_bigint *x, const int_bigint *y)
{
	int x_ans = 0, y_ans = 0;
	for (int i = 0; i < x->w->S; i++)
	{
		x_ans += hamming_weight(x->w->words[i]);
		y_ans += hamming_weight(y->w->words[i]);
	}

	if (x_ans > y_ans)
		return true;
	return false;
}

/**
* Takes in a sorted list of vertices taken from k's bag.  Fills the ind_sets
* list with all independent sets from bag_subset and translates these into masks that
* pick off entries in k's actual bag (so the masks are stretched in other words).
*/
int list_ind_sets(TDTree *T, int k, list<int> *bag_subset,
	list<TDSolution *> *ind_sets)
{
	// check for idiotic argument
	if (k < 0 || !(T->tree_nodes[k]))
		fatal_error("%s: tried to compute table for impossible tree node\n",
		__FUNCTION__);

	int i, j;
	TDSolution *new_set;
	list<int>::iterator ii, jj, kk;
	bigint_t current_mask(T->num_mask_words);
	list<bigint_t*>::iterator LLit;

	list<TDSolution*>::iterator ss, tt;
	vector<int> pos_vec(bag_subset->size(), -1);

	// For each position i in the list bag_subset, pos_vec[i] is the position of the ith entry
	// in bag_subset in tree node k's bag
	jj = T->tree_nodes[k]->bag.begin();
	i = j = 0;
	for (ii = bag_subset->begin(); ii != bag_subset->end(); ++ii)
	{
		while (*ii != *jj)
		{
			j++;
			++jj;
		}
		pos_vec[i] = j;
		i++;
	}

	// Clear the ind_sets list
	ind_sets->clear();
	int w = (int) bag_subset->size();
	// Copy the bag (stored as a list) to a w-long bag_subset_vec
	vector<int> bag_subset_vec(w, 0);
	i = 0;
	for (ii = bag_subset->begin(); ii != bag_subset->end(); ++ii)
	{
		bag_subset_vec[i] = *ii;
		i++;
	}

	// Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
	vector<int_bigint *> nbr_mask_vec(w);
	Graph::Node *n1;
	list<int> *nbrs;
	for (i = 0; i < w; i++)
	{
		nbr_mask_vec[i] = new int_bigint(T->num_mask_words);
		nbr_mask_vec[i]->k = i;
		for (j = 0; j < w; j++)
		{
			// Consider all poss. edges bag_subset_vec[i]-bag_subset_vec[j]
			n1=T->G->get_node(bag_subset_vec[i]);
			nbrs=n1->get_nbrs_ptr();
			//for (ii = T->G->nodes[bag_subset_vec[i]].nbrs.begin(); ii
			//		!= T->G->nodes[bag_subset_vec[i]].nbrs.end(); ++ii)
			for(ii=nbrs->begin();ii!=nbrs->end();++ii)
			{
				if (i != j && bag_subset_vec[j] == *ii)
				{
					nbr_mask_vec[i]->w->or_bit(j);
					T->tree_nodes[k]->num_subgraph_edges++;
					// No point looking further, so break out
					break;
				}
			}
		}
	}
	T->tree_nodes[k]->num_subgraph_edges = T->tree_nodes[k]->num_subgraph_edges
		/ 2;
	// Sort the masks and put the densest first 
	sort(nbr_mask_vec.begin(), nbr_mask_vec.end(), hamming_compare);

	bool is_independent;
	current_mask = 0;
	bigint_t two_pow_w(T->num_mask_words);
	bigint_t one(T->num_mask_words);
	one = BIGINT_ONE;
	two_pow_w.two_x(w);
	int two_w_word = 1, num_sets = 0;
	while (w >= BIGINT_WORD_SIZE * two_w_word)
		two_w_word++;
	two_w_word--;

	bool has_bit;

	while (current_mask.words[two_w_word] < two_pow_w.words[two_w_word])
	{
		is_independent = true;
		for (i = 0; i < w; i++)
		{
			if (current_mask.test_bit(nbr_mask_vec[i]->k))
			{
				// current_mask has the bit corresponding to nbr_mask_vec[i] set - so and with the appropriate vector for a hit
				has_bit = false;
				// Do the loop below manually rather than calling the overloaded AND
				for (j = 0; j < current_mask.S; j++)
				{
					if (current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
					{
						// Some nbr is found in word j
						has_bit = true;
						break;
					}
				}

				if (has_bit)
				{
					is_independent = false;
#if 0
					// Find the rightmost bit in the current_mask's word where we had the intersection
					// with the edge mask and then add this to the current mask - a nice one liner!
					current_mask.words[j] += (current_mask.words[j] & ~(current_mask.words[j]-1));		
#else
					int rpos = 0;
					// Find the rightmost bit in current_mask - this could be sped up
					while (!(current_mask.test_bit(rpos)))
						rpos++;

					// Do the add w/ only bitwise op's
					// Would be faster here to figure out which word we are in and 
					// do this w/o test_bit
					int b = 0;
					while (current_mask.test_bit(rpos + b))
						b++;

					for (int a = 0; a < b; a++)
						current_mask.xor_bit(rpos + a);

					current_mask.or_bit(rpos + b);
#endif
					// Break out of the for loop since current_mask cannot represent
					// an independent set
					break;
				}
			}
		}

		if (is_independent)
		{
			num_sets++;
			new_set = new TDSolution(T->num_mask_words);
			// Now translate the mask by translating the bits from current_mask into bits in the parent's full bag
			*(new_set->mask) = 0;
			new_set->value = 0;
			for (j = 0; j < w; j++)
			{
				if (current_mask.test_bit(j))
				{
					// bit j is set in current_mask - that corresponds to bit pos_vec[j] in a full mask
					new_set->mask->set_bit(pos_vec[j]);
					new_set->value += T->G->get_vertex_weight(bag_subset_vec[j]);
					//T->G->weight[bag_subset_vec[j]];
				}
			}
			// Add the new set to the list here -
			ind_sets->push_back(new_set);

			// Advance to next candidate ind. set
			++current_mask;
		}
	}

	for (i = 0; i < w; i++)
		delete nbr_mask_vec[i];
	// Return the size of the list
	return num_sets;
}

/**
* Computes the independent sets at a non-leaf node using information from the children.
* The idea is that any independent set in the parent is the union of a set that lives
* only in the parent and a set that lives in parent \cap child.
* Then looks "up" to the parent to find the
* ind. set in this bag that has the largest weight with a given intersection
* with the parent.  Letting M denote the mask of this intersection in the
* parent's bag, weight_table[M]=w where w is the weight of this largest set.
*/
int compute_nonnice_table_parent_child(TDTree *T, int k)
{
	// if k is a leaf, then just compute it the normal way
	if (T->tree_nodes[k]->adj.size() == 1)
	{
		// k is a leaf!
		int retval = compute_nonnice_table_standard(T, k);
		return retval;
	}

	int i, j;
	list<int>::iterator ii, jj, kk;
	bigint_t current_mask(T->num_mask_words);
	list<bigint_t*>::iterator current_mask_it;
	int w = (int) T->tree_nodes[k]->bag.size();
	vector<int> I(w);
	vector<int> D(w);
	vector<int>::iterator uu, vv;
	list<TDSolution*>::iterator ss, tt;

	vector<int> pos_vec(T->G->get_num_nodes(), -1);
	i = 0;
	for (ii = T->tree_nodes[k]->bag.begin(); ii != T->tree_nodes[k]->bag.end(); ++ii)
	{
		pos_vec[*ii] = i;
		i++;
	}

	// Some vectors to compute intersection with parent - for hash table
	int parent = T->tree_nodes[k]->adj.front();
	bool is_independent;
	T->tree_nodes[k]->obj_val = 0;
	// Use hash_ptr to look up things in the table (no allocation)
	TDSolution *hash_ptr = NULL;
	int max_intersection = -1, best_child = -1, diff_size;
	T->tree_nodes[k]->bag.sort();
	// k is not a leaf - find the child whose bag has the largest intersection
	ii = T->tree_nodes[k]->adj.begin();
	++ii;
	for (; ii != T->tree_nodes[k]->adj.end(); ++ii)
	{
		T->tree_nodes[*ii]->bag.sort();
		vv = set_intersection(T->tree_nodes[k]->bag.begin(),
			T->tree_nodes[k]->bag.end(), T->tree_nodes[*ii]->bag.begin(),
			T->tree_nodes[*ii]->bag.end(), I.begin());
		if (vv - I.begin() > max_intersection)
		{
			max_intersection = vv - I.begin();
			best_child = *ii;
		}
	}

	// Compute parent - best_child
	vv = set_difference(T->tree_nodes[k]->bag.begin(),
		T->tree_nodes[k]->bag.end(),
		T->tree_nodes[best_child]->bag.begin(),
		T->tree_nodes[best_child]->bag.end(), D.begin());
	diff_size = vv - D.begin();
	if (diff_size + max_intersection != w)
		fatal_error("%s:  Incorrect bag size computation %d!=%d??\n",
		__FUNCTION__, diff_size + max_intersection, w);
	list<int> bag_subset;
	list<TDSolution *> ind_sets;
	for (uu = D.begin(); uu != vv; ++uu)
		bag_subset.push_back(*uu);

	vector<bool> subset_pos(T->G->get_num_nodes(), false);
	for (ii = bag_subset.begin(); ii != bag_subset.end(); ++ii)
		subset_pos[*ii] = true;

	// The cvec[]'s allow us to quickly do the DP with the children by getting their
	// contribution to each ind. set in k's bag
	vector<int> children(T->tree_nodes[k]->adj.size() - 1);
	ii = T->tree_nodes[k]->adj.begin();
	++ii;
	i = 0;
	for (; ii != T->tree_nodes[k]->adj.end(); ++ii)
	{
		children[i] = *ii;
		i++;
	}
	vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
	for (i = 0; i < (int) T->tree_nodes[k]->adj.size() - 1; i++)
	{
		cvec[i] = new bool[T->G->get_capacity()];
		for (j = 0; j < T->G->get_capacity(); j++)
			cvec[i][j] = false;
		for (ii = T->tree_nodes[children[i]]->bag.begin(); ii
			!= T->tree_nodes[children[i]]->bag.end(); ++ii)
			cvec[i][*ii] = true;
	}

	// Create all independent sets that live only in the difference between the parent and best_child
	list_ind_sets(T, k, &bag_subset, &ind_sets);

	TDSolution *s, *tmp, temp_sol(T->num_mask_words);
	bigint_t orig_mask(T->num_mask_words);
	bool is_best;
	int num_in_hash = 0, num_hash_ind_found = 0, num_ind_sets = 0;
	int new_value = 0;
	// Use HASH_ITER to loop through best_child's hash table    
	HASH_ITER(hh,T->tree_nodes[best_child]->hash_table,s,tmp)
	{
		// s is the mask of the intersection in the best_child's hash_table
		for (ss = ind_sets.begin(); ss != ind_sets.end(); ++ss)
		{
			// ss is the ind. set that lives only in the parent
			// Need to see if s | ss represents an independent set - the only possible edges
			// are between s and ss 
			is_independent = true;
			for (i = 0; i < w; i++)
			{
				if (subset_pos[T->tree_nodes[k]->bag_vec[T->tree_nodes[k]->nbr_mask_vec->at(i)->k]]
				&& (*ss)->mask->test_bit(T->tree_nodes[k]->nbr_mask_vec->at(i)->k))
				{
					// Do the loop below manually rather than calling the overloaded AND
					for (j = 0; j < s->mask->S; j++)
					{
						if (s->mask->words[j] & T->tree_nodes[k]->nbr_mask_vec->at(i)->w->words[j])
							//if (s->mask->words[j] & nbr_mask_vec[i]->w->words[j])
						{
							// Some nbr is found in word j
							is_independent = false;
							break;
						}
					}
				}
			}
			if (is_independent)
			{
				num_ind_sets++;
				// increment total_table_entries
				T->info->total_table_entries++;
				// The real ind. set is s | ss - place it in orig_mask
				orig_mask.zeroize();
				orig_mask = *(s->mask);
				orig_mask |= *((*ss)->mask);
				is_best = false;
				int parent_weight = 0;
				new_value = 0;
				for (j = 0; j < w; j++)
				{
					if (s->mask->test_bit(j) || (*ss)->mask->test_bit(j))
						new_value += T->G->get_vertex_weight(T->tree_nodes[k]->bag_vec[j]);
					//T->G->weight[T->tree_nodes[k]->bag_vec[j]];
				}
				// Check for new max_obj
				if (new_value > T->tree_nodes[k]->obj_val)
				{
					T->tree_nodes[k]->obj_val = new_value;
					is_best = true;
				}
				for (i = 0; i < (int) T->tree_nodes[k]->adj.size() - 1; i++)
				{
					temp_sol.mask->zeroize();
					for (int z = 0; z < w; z++)
					{
						// This part could be sped up by considering the words and then doing some xoring, I think
						if (cvec[i][T->tree_nodes[k]->bag_vec[z]])
						{
							if (s->mask->test_bit(z)
								|| (*ss)->mask->test_bit(z))
								// Bit z is set, and the corresponding vertex is also in the child
								temp_sol.mask->set_bit(z);
						}
					}
					hash_ptr = NULL;
					HASH_FIND(hh, T->tree_nodes[children[i]]->hash_table, temp_sol.mask->words,
						T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
					if (!hash_ptr)
						fatal_error(
						"%s: Did not find set in child!!! Cannot continue\n",
						__FUNCTION__);
					new_value += hash_ptr->value;
					// Check for new max obj
					if (new_value > T->tree_nodes[k]->obj_val)
					{
						T->tree_nodes[k]->obj_val = new_value;
						is_best = true;
					}
				}
				// Now look up to the parent
				temp_sol.mask->zeroize();
				parent_weight = 0;
				for (j = 0; j < w; j++)
				{
					if (s->mask->test_bit(j) || (*ss)->mask->test_bit(j))
					{
						//if (parent_location_vec[T->tree_nodes[k]->bag_vec[j]]!= -1)
						if(T->tree_nodes[k]->parent_position->at(j)!=-1)
						{
							temp_sol.mask->set_bit(T->tree_nodes[k]->parent_position->at(j));
							//parent_location_vec[T->tree_nodes[k]->bag_vec[j]]);
							parent_weight
								+=  T->G->get_vertex_weight(T->tree_nodes[k]->bag_vec[j]);
						}
					}
				}
				// temp_mask is a mask of this set \cap parent_bag
				hash_ptr = NULL;
				temp_sol.value = new_value;
				// Only add the best set to the root!
				if (k == T->root_node)
				{
					if (is_best)
					{
						TDSolution *added_set;
						added_set = new TDSolution(T->num_mask_words);
						*(added_set->mask) = *(temp_sol.mask);
						// Set the orig_mask
						*(added_set->orig_mask) = orig_mask;
						added_set->value = temp_sol.value;
						HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
							T->num_mask_words*sizeof(BIGINT_WORD), added_set);
					}
				}
				else
				{
					// k is not the root - so do the usual
					HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_sol.mask->words,
						T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
					if (hash_ptr)
					{
						if (temp_sol.value - parent_weight > hash_ptr->value)
						{
							// We found an entry for this mask and this mask has a better value 
							// so update the entry in the table
							hash_ptr->value = temp_sol.value - parent_weight;
							// Set the orig_mask
							*(hash_ptr->orig_mask) = orig_mask;
						}
					}
					else

					{
						// There is no entry for this mask - add it to the table
						TDSolution *added_set;
						added_set = new TDSolution(T->num_mask_words);
						*(added_set->mask) = *(temp_sol.mask);
						// Set the orig_mask
						*(added_set->orig_mask) = orig_mask;
						added_set->value = temp_sol.value - parent_weight;
						HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
							T->num_mask_words*sizeof(BIGINT_WORD), added_set);
					}
				}
				num_hash_ind_found++;
			}
			// end if(is_independent)
		}
		// end for(ss...)
		num_in_hash++;
	}

	// Done with ind_sets now - so clear it
	for (ss = ind_sets.begin(); ss != ind_sets.end(); ++ss)
		delete *ss;
	ind_sets.clear();
	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
	T->info->total_pc_table_entries += (unsigned long long) table_size;
	T->info->aux_info[k] = num_ind_sets;

	return table_size;
}

/**
* Computes the table for a leaf or non-leaf node in a non-nice TD. Uses only
* the hash_table and only stores an entry for the unique intersections of node
* k independent sets with the parent's.
*/
int compute_nonnice_table_standard(TDTree *T, int k)
{
	int i, j, w = (int) T->tree_nodes[k]->bag.size();
	list<int>::iterator ii, jj;

	bigint_t current_mask(T->num_mask_words);
	list<bigint_t*>::iterator current_mask_it;
	list<TDSolution*>::iterator ss, tt;

	// We already have all the necessary information inside the TDTreeNode
	current_mask = 0;
	// Find out when to stop the loop
	bigint_t two_pow_w(T->num_mask_words), one(T->num_mask_words), temp(T->num_mask_words);
	one = BIGINT_ONE;
	two_pow_w.two_x(w);
	int two_w_word = 1;
	while (w >= BIGINT_WORD_SIZE * two_w_word)
		two_w_word++;
	two_w_word--;
	bool is_independent, has_bit;
	int num_ind_sets = 0;
	TDSolution *hash_ptr, temp_set(T->num_mask_words);

	clock_t loop_start = clock();
	while (current_mask.words[two_w_word] < two_pow_w.words[two_w_word])
	{
		is_independent = true;
		for (i = 0; i < w; i++)
		{
			if (current_mask.test_bit(T->tree_nodes[k]->nbr_mask_vec->at(i)->k))
			{
				// current_mask has the bit corresponding to nbr_mask_vec[i] set - so and with the appropriate vector for a hit
				has_bit = false;
				// Do the loop below manually rather than calling the overloaded AND
				for (j = 0; j < current_mask.S; j++)
				{
					if ( (current_mask.words[j]
					& T->tree_nodes[k]->nbr_mask_vec->at(i)->w->words[j]) )
					{
						// Some nbr is found in word j so current_mask cannot represent an ind. set
						has_bit = true;
						break;
					}
				}
				// breaks to here
				if (has_bit)
				{
					// Set represented by current_mask is not independent
					is_independent = false;
					int rpos = 0;

#if 0
					// Unfortunately this won't work since bigint doesn't know the + operator
					// I should have used GMP from day one, esp. now that I'm using it for the memory
					// estimation
					// Find the rightmost bit in the current_mask's word where we had the intersection
					// with the edge mask and then add this to the current mask - a nice one liner!
					current_mask.words[j] += (current_mask.words[j] & ~(current_mask.words[j]-1));		
#else
					// Find the rightmost bit in current_mask - this could be sped up by looking at the words
					// SPEEDUP
					while (!(current_mask.test_bit(rpos)))
						rpos++;

					// Now advance current_mask by a potentially massive amount
					// Would be faster here to figure out which word we are in and 
					// do this w/o test_bit
					int b = 0;
					while (current_mask.test_bit(rpos + b))
						b++;
					for (int a = 0; a < b; a++)
					{
						current_mask.xor_bit(rpos + a);
					}
					current_mask.or_bit(rpos + b);
#endif
					break;
				}
			}
		}
		if (is_independent)
		{
			// increment total_table_entries
			T->info->total_table_entries++;
			num_ind_sets++;
			bool is_best = false;
			int new_value = 0, parent_weight = 0;
			for (j = 0; j < w; j++)
			{
				if (current_mask.test_bit(j))
				{
					// Compute the actual weight in G of this ind. set
					new_value += T->G->get_vertex_weight(T->tree_nodes[k]->bag_vec[j]);
					//T->G->weight[T->tree_nodes[k]->bag_vec[j]];
				}
			}

			// Check for new max obj
			if (new_value > T->tree_nodes[k]->obj_val)
			{
				T->tree_nodes[k]->obj_val = new_value;
				is_best = true;
			}

			// Now do the DP over the children for this independent set - if this is a leaf, then
			// the body of the loop is never executed
			ii = T->tree_nodes[k]->adj.begin();
			++ii;
			i = 0;
			for (; ii != T->tree_nodes[k]->adj.end(); ++ii)
			{
				temp_set.mask->zeroize();
				*(temp_set.mask) = current_mask;
				// Get the intersection of this set with child i
				*(temp_set.mask) &= *(T->tree_nodes[k]->child_intersection->at(i));
				HASH_FIND(hh, T->tree_nodes[*ii]->hash_table, temp_set.mask->words,
					T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
				if (!hash_ptr)
				{
					fatal_error(
						"%s: Did not find set in child!!! Cannot continue\n",
						__FUNCTION__);

				}
				new_value += hash_ptr->value;
				// Check for new max obj
				if (new_value > T->tree_nodes[k]->obj_val)
				{
					T->tree_nodes[k]->obj_val = new_value;
					is_best = true;
				}
				i++;
			}
			// Done with child DP

			// Now look up to the parent to translate the masks
			temp_set.mask->zeroize();
			parent_weight = 0;
			for (j = 0; j < w; j++)
			{
				if (current_mask.test_bit(j))
				{
					if (T->tree_nodes[k]->parent_position->at(j) != -1)
					{
						temp_set.mask->set_bit(
							T->tree_nodes[k]->parent_position->at(j));
						parent_weight
							+= T->G->get_vertex_weight(T->tree_nodes[k]->bag_vec[j]);
					}
				}
			}
			// temp_mask is a mask of this set \cap parent_bag
			hash_ptr = NULL;
			temp_set.value = new_value;
			// Only add the best set to the root!
			if (k == T->root_node)
			{
				if (is_best)
				{
					TDSolution *added_set;
					added_set = new TDSolution(T->num_mask_words);
					*(added_set->mask) = *(temp_set.mask);
					*(added_set->orig_mask) = current_mask;
					added_set->value = temp_set.value;
					HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
						T->num_mask_words*sizeof(BIGINT_WORD), added_set);
				}
			}
			else
			{
				// k is not the root - so do the usual
				HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
					T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
				if (hash_ptr)
				{
					if (temp_set.value - parent_weight > hash_ptr->value)
					{
						// We found an entry for this mask and this mask has a better value 
						// so update the entry in the table
						hash_ptr->value = temp_set.value - parent_weight;
						*(hash_ptr->orig_mask) = current_mask;
					}
				}
				else
				{
					// There is no entry for this mask - add it to the table
					TDSolution *added_set;
					added_set = new TDSolution(T->num_mask_words);
					*(added_set->mask) = *(temp_set.mask);
					*(added_set->orig_mask) = current_mask;
					added_set->value = temp_set.value - parent_weight;
					HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
						T->num_mask_words*sizeof(BIGINT_WORD), added_set);
				}
			}
			// Advance to next candidate ind. set       
			++current_mask;
		}
	}
	//loop_end = clock();

	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
	T->info->total_pc_table_entries += (unsigned long long) table_size;
	T->info->aux_info[k] = num_ind_sets;
	return table_size;
}

/**
* Computes the table for tree node k by splitting the bag into
* two parts and finding the independent sets in each half. Then
* considers all pairs on independent sets from each half to construct
* the table for tree node k.
*/
int compute_nonnice_table_split_bag(TDTree *T, int k)
{
	int i, j, w = (int) T->tree_nodes[k]->bag.size(), parent =
		T->tree_nodes[k]->adj.front(), num_nbrs =
		(int) T->tree_nodes[k]->adj.size();
	int half_w = w / 2;
	list<int>::iterator ii, jj, kk;
	bigint_t current_mask(T->num_mask_words);
	list<bigint_t*>::iterator current_mask_it;
	vector<int> I(w), D(w),	children(num_nbrs - 1);

	list<TDSolution*>::iterator ss, tt;

	// Create a vector of the children
	ii = T->tree_nodes[k]->adj.begin();
	++ii;
	i = 0;
	for (; ii != T->tree_nodes[k]->adj.end(); ++ii)
	{
		children[i] = *ii;
		i++;
	}


	T->tree_nodes[k]->obj_val = 0;
	TDSolution *hash_ptr = NULL;

	// Here is where the splitting begins - first list all the ind. sets in the two halves
	list<int> first_half, second_half;
	ii = T->tree_nodes[k]->bag.begin();

	// Create the two halves 
	for (i = 0; i < half_w; i++)
	{
		first_half.push_back(*ii);
		++ii;
	}
	for (; i < w; i++)
	{
		second_half.push_back(*ii);
		++ii;
	}
	clock_t half_stop, half_start = clock();
	// Find the ind. sets in each half
	list<TDSolution *> first_sets, second_sets;
	list_ind_sets(T, k, &first_half, &first_sets);
	list_ind_sets(T, k, &second_half, &second_sets);
	int ss_size = (int) second_sets.size();

	// Need a vector of second sets
	vector<TDSolution *> second_vec(ss_size);
	i = 0;
	for (ss = second_sets.begin(); ss != second_sets.end(); ++ss)
	{
		second_vec[i] = *ss;
		i++;
	}

	bool **V;
	V = new bool*[half_w];
	V[0] = new bool[half_w * ss_size];
	for (i = 0; i < half_w * ss_size; i++)
		V[0][i] = false;
	for (i = 1; i < half_w; i++)
		V[i] = V[i - 1] + ss_size;

	// Now set V[i][j] if (1<<i)|S[j] is an ind. set
	for (i = 0; i < half_w; i++)
	{
		for (j = 0; j < ss_size; j++)
		{
			bool good = true;
			for (int m = 0; m < second_vec[j]->mask->S; m++)
			{
				if (second_vec[j]->mask->words[m] & T->tree_nodes[k]->nbr_mask_vec->at(i)->w->words[m])
				{
					good = false;
					break;
				}
			}
			// breaks to here
			if (good)
			{
				V[i][j] = true;
			}
		}
	}
	half_stop = clock();

	// Now go through ff and use the V[i]'s to determine independence
	list<TDSolution *>::iterator ff;
	list<int> L;
	int num_ind_sets = 0;
	TDSolution temp_set(T->num_mask_words);
	clock_t DP_time = 0, ind_time = 0;
	clock_t ind_start, ind_stop;
	bool *Lvec = new bool[ss_size];
	int ctr = 0;

	for (ff = first_sets.begin(); ff != first_sets.end(); ++ff)
	{
		ind_start = clock();
		// Set Lvec to true and then go through V to determine which second_sets
		// actually work
		for (i = 0; i < ss_size; i++)
			Lvec[i] = true;
		if (ctr != 0)
		{
			for (i = 0; i < half_w; i++)
			{
				if ((*ff)->mask->test_bit(i))
				{
					// And V[i] to Lvec
					for (j = 0; j < ss_size; j++)
					{
						if (!V[i][j])
							// Equivalent to Lvec[j] &= V[i][j] but this operator 
							// doesn't exist
							Lvec[j] = false;
					}
				}
			}
		}
		ind_stop = clock();
		ind_time += (ind_stop - ind_start);

		clock_t DP_start = clock();
		clock_t DP_end;
		// Now we know that ff|second_vec[j] is an independent set for every j with Lvec[j]=true
		// Do the DP for all such j now
		for (int x = 0; x < ss_size; x++)
		{
			while (!Lvec[x] && x < ss_size)
				x++;
			if (x == ss_size)
				break;

			// This is a new ind. set, so increment total_table_entries
			T->info->total_table_entries++;
			num_ind_sets++;

			int new_value, parent_weight = 0;
			current_mask = *((*ff)->mask);
			current_mask |= (*(second_vec[x]->mask));
			new_value = (*ff)->value + second_vec[x]->value;
			bool is_best = false;
			// Check for new max obj
			if (new_value > T->tree_nodes[k]->obj_val)
			{
				T->tree_nodes[k]->obj_val = new_value;
				is_best = true;
			}

			// Now do the DP over the children for this independent set - if this is a leaf, then
			// the body of the loop is never executed
			for (i = 0; i < num_nbrs - 1; i++)
			{
				temp_set.mask->zeroize();
				for (int z = 0; z < w; z++)
				{
					// This part could be sped up by considering the words and then doing some xoring, I think
					if(T->tree_nodes[k]->child_intersection->at(i)->test_bit(z))
					{
						if (current_mask.test_bit(z))
							// Bit z is set, and the corresponding vertex is also in the child
							temp_set.mask->set_bit(z);
					}
				}

				HASH_FIND(hh, T->tree_nodes[children[i]]->hash_table, temp_set.mask->words,
					T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
				if (!hash_ptr)
					fatal_error(
					"%s: Did not find set in child!!! Cannot continue\n",
					__FUNCTION__);
				new_value += hash_ptr->value;
				// Check for new max obj
				if (new_value > T->tree_nodes[k]->obj_val)
				{
					T->tree_nodes[k]->obj_val = new_value;
					is_best = true;
				}
			}
			// Done with child DP

			// Now look up to the parent
			temp_set.mask->zeroize();
			parent_weight = 0;
			for (j = 0; j < w; j++)
			{
				if (current_mask.test_bit(j))
				{
					//if (parent_location_vec[T->tree_nodes[k]->bag_vec[j]] != -1)
					if(T->tree_nodes[k]->parent_position->at(j)!=-1)
					{
						temp_set.mask->set_bit(T->tree_nodes[k]->parent_position->at(j));
						//		parent_location_vec[T->tree_nodes[k]->bag_vec[j]]);
						parent_weight+=T->G->get_vertex_weight(T->tree_nodes[k]->bag_vec[j]);
						//		+= T->G->weight[T->tree_nodes[k]->bag_vec[j]];
					}
				}
			}
			// temp_mask is a mask of this set \cap parent_bag
			hash_ptr = NULL;
			temp_set.value = new_value;
			// Only add the best set to the root!
			if (k == T->root_node)
			{
				if (is_best)
				{
					TDSolution *added_set;
					added_set = new TDSolution(T->num_mask_words);
					*(added_set->mask) = *(temp_set.mask);
					*(added_set->orig_mask) = current_mask;
					added_set->value = temp_set.value;
					HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
						T->num_mask_words*sizeof(BIGINT_WORD), added_set);
				}
			}
			else
			{
				// k is not the root - so do the usual
				HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
					T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
				if (hash_ptr)
				{
					if (temp_set.value - parent_weight > hash_ptr->value)
					{
						// We found an entry for this mask and this mask has a better value 
						// so update the entry in the table
						hash_ptr->value = temp_set.value - parent_weight;
						*(hash_ptr->orig_mask) = current_mask;
					}
				}
				else
				{
					// There is no entry for this mask - add it to the table
					TDSolution *added_set;
					added_set = new TDSolution(T->num_mask_words);
					*(added_set->mask) = *(temp_set.mask);
					*(added_set->orig_mask) = current_mask;
					added_set->value = temp_set.value - parent_weight;
					HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
						T->num_mask_words*sizeof(BIGINT_WORD), added_set);
				}
			}
		}
		// end loop over entries in L
		DP_end = clock();
		DP_time += (DP_end - DP_start);
		ctr++;
	}

	if (T->info->very_verbose)
		// Print out some timing info for the DP
		print_message(0,
		"\n\tHALF: %f secs\n\tIND : %f secs\n\tDP  : %f secs\n",
		(double) (half_stop - half_start) / CLOCKS_PER_SEC,
		(double) ind_time / CLOCKS_PER_SEC,
		(double) DP_time / CLOCKS_PER_SEC);

	// Need to delete first and second lists - iterator has to be const_ to do this!!
	for (list<TDSolution *>::const_iterator dd = first_sets.begin(); dd
		!= first_sets.end(); ++dd)
		delete *dd;
	first_sets.clear();
	for (list<TDSolution *>::const_iterator dd = second_sets.begin(); dd
		!= second_sets.end(); ++dd)
		delete *dd;
	second_sets.clear();


	delete[] V[0];
	delete[] V;
	delete[] Lvec;

	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
	T->info->total_pc_table_entries += (unsigned long long) table_size;
	T->info->aux_info[k] = num_ind_sets;

	return table_size;

}

/**
* Function to add given mask into parent's hash table.
*/
void add_to_parent_intersection_table(TDTree *T, int k, int w, bigint_t *mask,
	int value)
{
	TDSolution temp_set(T->num_mask_words);
	bigint_t current_mask(T->num_mask_words);
	TDSolution *hash_ptr;
	int j = 0;

	temp_set.mask->zeroize();
	current_mask = *mask;
	for (j = 0; j < w; j++)
	{
		if (current_mask.test_bit(j))
		{
			if (T->tree_nodes[k]->parent_position->at(j) != -1)
			{
				temp_set.mask->set_bit(T->tree_nodes[k]->parent_position->at(j));
			}
		}
	}

	temp_set.value = value;
	// k is not the root - so do the usual
	HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
		T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

	if (hash_ptr)
	{
		if (temp_set.value > hash_ptr->value)
		{
			// We found an entry for this mask and this mask has a better value
			// so update the entry in the table
			hash_ptr->value = temp_set.value;
			*(hash_ptr->orig_mask) = current_mask;
		}
	}
	else
	{
		// There is no entry for this mask - add it to the table
		TDSolution *added_set;
		added_set = new TDSolution(T->num_mask_words);
		*(added_set->mask) = *(temp_set.mask);
		*(added_set->orig_mask) = current_mask;
		added_set->value = temp_set.value;
		HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
			T->num_mask_words*sizeof(BIGINT_WORD), added_set);
	}
}

/**
* Computes the table for a leaf or non-leaf node in a non-nice TD. 
* Parent does not have to wait for all children to finish.
* As soon as one child's table is computed, the parent can compute 
* its table and free that child table.
*/
int compute_nonnice_table_async(TDTree *T, int k)
{
	/**
	* TODO: In parallel implementation to dynamic programming intersection with the children should be found.
	* Therefore mask which represents child intersection should also be sent.
	*/

	bool is_leaf_node = false;
	int table_size = 0;

	if (T->tree_nodes[k]->adj.size() <= 1)
	{
		//compute parent child intersection for leaf nodes.
		table_size = compute_nonnice_table_standard(T, k);
		is_leaf_node = true;
	}

	if (T->root_node == k)
	{
		table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
		return table_size;
	}

	bigint_t tmp_set(T->num_mask_words);
	//compute complete table for parent.
	int parent = T->tree_nodes[k]->adj.front();
	int i, j;
	int w = (int) T->tree_nodes[parent]->bag.size();
	TDSolution *tmp_sol;
	bigint_t current_mask(T->num_mask_words);
	TDSolution *hash_ptr = NULL;

	DP_prepare_node(T, parent);

	//DEBUG("k: %d parent: %d\n", k, parent);
	//DEBUG("hash count %d:%d\n", HASH_COUNT(T->tree_nodes[k]->hash_table), HASH_COUNT(T->tree_nodes[k]->complete_table));

	vector<int> intersection(w);
	bigint_t child_intersection(T->num_mask_words);

	list<int>::iterator ii, jj;
	vector<int>::iterator uu, vv;
	list<bigint_t*>::iterator current_mask_it;
	list<TDSolution*>::iterator ss, tt;
	int num_ind_sets = 0;

	int inter_weight = 0;

	vv = set_intersection(T->tree_nodes[parent]->bag.begin(),
		T->tree_nodes[parent]->bag.end(), T->tree_nodes[k]->bag.begin(),
		T->tree_nodes[k]->bag.end(), intersection.begin());

	jj = T->tree_nodes[parent]->bag.begin();
	j = 0;
	for (uu = intersection.begin(); uu != vv; ++uu)
	{
		while (*uu != *jj)
		{
			++jj;
			j++;
		}
		child_intersection.set_bit(j);
	}

	/**
	* check whether parent has a table to store all independent sets.
	*/
	if (T->tree_nodes[parent]->all_ind_sets.empty())
	{
		// We already have all the necessary information inside the TDTreeNode
		current_mask = 0;
		// Find out when to stop the loop
		bigint_t two_pow_w(T->num_mask_words), one(T->num_mask_words);
		one = BIGINT_ONE;
		two_pow_w.two_x(w);
		int two_w_word = 1;
		while (w >= BIGINT_WORD_SIZE * two_w_word)
			two_w_word++;
		two_w_word--;
		bool is_independent, has_bit;

		TDSolution *hash_ptr, temp_set(T->num_mask_words);
		clock_t loop_start = clock();

		while (current_mask.words[two_w_word] < two_pow_w.words[two_w_word])
		{
			is_independent = true;
			for (i = 0; i < w; i++)
			{
				if (current_mask.test_bit(
					T->tree_nodes[parent]->nbr_mask_vec->at(i)->k))
				{
					// current_mask has the bit corresponding to nbr_mask_vec[i] set - so and with the appropriate vector for a hit
					has_bit = false;
					// Do the loop below manually rather than calling the overloaded AND
					for (j = 0; j < current_mask.S; j++)
					{
						if (current_mask.words[j]
						& T->tree_nodes[parent]->nbr_mask_vec->at(i)->w->words[j])
						{
							// Some nbr is found in word j so current_mask cannot represent an ind. set
							has_bit = true;
							break;
						}
					}
					// breaks to here
					if (has_bit)
					{
						// Set represented by current_mask is not independent
						is_independent = false;
#if 0
						// Find the rightmost bit in the current_mask's word where we had the intersection
						// with the edge mask and then add this to the current mask - a nice one liner!
						current_mask.words[j] += (current_mask.words[j] & ~(current_mask.words[j]-1));		
#else
						int rpos = 0;
						// Find the rightmost bit in current_mask - this could be sped up by looking at the words
						// SPEEDUP
						while (!(current_mask.test_bit(rpos)))
							rpos++;

						// Now advance current_mask by a potentially massive amount
						// Would be faster here to figure out which word we are in and
						// do this w/o test_bit
						int b = 0;
						while (current_mask.test_bit(rpos + b))
							b++;
						for (int a = 0; a < b; a++)
							current_mask.xor_bit(rpos + a);
						current_mask.or_bit(rpos + b);
#endif
						break;
					}
				}
			}

			if (is_independent)
			{
				// increment total_table_entries
				T->info->total_table_entries++;
				num_ind_sets++;
				int new_value = 0;
				for (j = 0; j < w; j++)
				{
					if (current_mask.test_bit(j))
					{
						// Compute the actual weight in G of this ind. set
						new_value+=T->G->get_vertex_weight(T->tree_nodes[parent]->bag_vec[j]);
						//		+= T->G->weight[T->tree_nodes[parent]->bag_vec[j]];
					}
				}

				// Now do the DP for the intersection of current_mask using tree node k's table.
				// Compute mask wrt k's bag.
				tmp_set.zeroize();
				tmp_set = current_mask;
				tmp_set &= child_intersection;
				hash_ptr = NULL;

				inter_weight = 0;
				if (!is_leaf_node)
				{
					int i = 0;
					for (i = 0; i < w; i++)
					{
						if (tmp_set.test_bit(i))
						{
							inter_weight+=T->G->get_vertex_weight(T->tree_nodes[parent]->bag_vec[i]);
							//		+= T->G->weight[T->tree_nodes[parent]->bag_vec[i]];
						}
					}
				}
				new_value -= inter_weight;

				HASH_FIND(hh,T->tree_nodes[k]->hash_table, tmp_set.words,
					T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
				if (hash_ptr)
					new_value += hash_ptr->value;
				else
				{
					fatal_error(
						"can not find intersection for current mask, abort !!\n");
					current_mask.print(__error_file__);
				}

				if (T->tree_nodes[parent]->adj.size() > 2)
				{
					bigint_t *cmask = new bigint_t(T->num_mask_words);
					*cmask = current_mask;
					T->tree_nodes[parent]->all_ind_sets.push_back(cmask);
					T->tree_nodes[parent]->all_ind_set_values.push_back(
						new_value);
				}
				else
				{
					add_to_parent_intersection_table(T, parent, w,
						&current_mask, new_value);
				}

				// Done with child DP - advance to look for the next independent set in the parent's bag
				++current_mask;
			}
		}
		//loop_end = clock();

		T->info->aux_info[parent] = num_ind_sets;
		// finish complete table creation and computation for k's parent.
	}
	else
	{
		// The independent sets for k's parent have already been computed, so just loop over them
		// and do the DP treating tree node k as the child
		// Update k's parent complete table
		list<bigint_t *>::iterator bi_iter;
		list<int>::iterator i_iter;

		for (bi_iter = T->tree_nodes[parent]->all_ind_sets.begin(), i_iter
			= T->tree_nodes[parent]->all_ind_set_values.begin(); bi_iter
			!= T->tree_nodes[parent]->all_ind_sets.end(); ++bi_iter, ++i_iter)
		{
			//HASH_ITER(hh, T->tree_nodes[parent]->complete_table, current_set, tmp_sol)
			//{
			tmp_set.zeroize();
			tmp_set = *(*bi_iter);
			tmp_set &= child_intersection;
			hash_ptr = NULL;

			//calculating weights for the nodes in the intersection
			inter_weight = 0;
			if (!is_leaf_node)
			{
				int i = 0;
				for (i = 0; i < w; i++)
				{
					if (tmp_set.test_bit(i))
					{
						inter_weight
							+=T->G->get_vertex_weight(T->tree_nodes[parent]->bag_vec[i]);
						//+= T->G->weight[T->tree_nodes[parent]->bag_vec[i]];
					}
				}
			}

			HASH_FIND(hh,T->tree_nodes[k]->hash_table, tmp_set.words, T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
			if (hash_ptr)
			{
				*i_iter += hash_ptr->value;
				*i_iter -= inter_weight;
			}
			else
			{
				fatal_error(
					"can not find intersection for current mask, abort !!\n");
				current_mask.print(__error_file__);
			}
		}
		//finish update table for k's parent
	}

	// Free k's hash table since we updated his parent's table using this information
	table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
	T->info->total_pc_table_entries += (unsigned long long) table_size;
	if (T->info->free_children)
	{
		TDSolution *added_set;
		HASH_ITER(hh, T->tree_nodes[k]->hash_table, added_set, tmp_sol)
		{
			HASH_DEL(T->tree_nodes[k]->hash_table, added_set);
			delete added_set;
		}
		T->tree_nodes[k]->hash_table = NULL;
	}

	// increase completed children count.
	T->tree_nodes[parent]->increase_completed_children();
	if (!(T->tree_nodes[parent]->get_completed_children()
		< (int) T->tree_nodes[parent]->adj.size()))
	{
		// All the children have been processed so now go through the complete table
		// and create the hash table for the parent node. Since the masks in the hash table
		// are stored in terms of the parent's bag, this requires you to now look up to the
		// parent's parent using the parent_position information!
		//Save parent child intersection in parent hash and clean completed_table
		int max = 0;
		list<bigint_t *>::iterator bi_iter;
		list<int>::iterator i_iter;

		for (bi_iter = T->tree_nodes[parent]->all_ind_sets.begin(), i_iter
			= T->tree_nodes[parent]->all_ind_set_values.begin(); bi_iter
			!= T->tree_nodes[parent]->all_ind_sets.end(); ++bi_iter, ++i_iter)
		{
			if (parent == T->root_node)
			{
				TDSolution *added_set;
				if (*i_iter > max)
				{
					max = *i_iter;
					added_set = new TDSolution(T->num_mask_words);
					*(added_set->mask) = *(*bi_iter);
					*(added_set->orig_mask) = *(*bi_iter);
					added_set->value = *i_iter;
					HASH_ADD_KEYPTR(hh, T->tree_nodes[parent]->hash_table, added_set->mask->words,
						T->num_mask_words*sizeof(BIGINT_WORD), added_set);
				}
				delete *bi_iter;
				continue;
			}

			add_to_parent_intersection_table(T, parent, w, *bi_iter, *i_iter);
			delete *bi_iter;
		}
		T->tree_nodes[parent]->all_ind_set_values.clear();
		T->tree_nodes[parent]->all_ind_sets.clear();
	}
	return table_size;
}

/** 
* Wrapper for the various nonnice table functions. The computation varies depending
* on the options to the dynamic programming.
*/
int compute_nonnice_table(TDTree *T, int k)
{
	int retval;
	if (T->info->parent_child)
	{
		retval = compute_nonnice_table_parent_child(T, k);
		return retval;
	}

	if (T->info->split_bag)
	{
		retval = compute_nonnice_table_split_bag(T, k);
		return retval;
	}

	if (T->info->async_tbl)
	{
		retval = compute_nonnice_table_async(T, k);
		return retval;
	}

	// Default
	retval = compute_nonnice_table_standard(T, k);
	return retval;
}

// Functions specific to a NICE Tree Decomposition

/**
* Computes the independent sets at a leaf node in a nice TD.  Does not look up to the
* parent and worry about the intersection.
*/
int compute_ind_sets_nice(TDTree *T, int k)
{
	int i, j;
	int w = (int) T->tree_nodes[k]->bag.size();
	TDSolution *new_set;
	list<int>::iterator ii, jj, kk;
	bigint_t current_mask(T->num_mask_words);
	list<bigint_t*>::iterator current_mask_it;
	vector<int> I(w);
	vector<int> D(w);

	list<TDSolution*>::iterator ss, tt;

	int parent = T->tree_nodes[k]->adj.front();
	bool is_independent;
	T->tree_nodes[k]->obj_val = 0;
	int num_added = 0;

	// Add the empty set to the table
	new_set = new TDSolution(T->num_mask_words);
	*(new_set->mask) = (BIGINT_WORD) 0;
	*(new_set->orig_mask) = (BIGINT_WORD) 0;
	new_set->value = 0;
	HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, new_set->mask->words,
		T->num_mask_words*sizeof(BIGINT_WORD), new_set);
	T->info->total_pc_table_entries++;
	T->info->total_table_entries++;
	num_added++;
	// k is a leaf - have to do it from scratch
	current_mask = 1;
	bigint_t two_pow_w(T->num_mask_words);
	bigint_t one(T->num_mask_words);
	one = BIGINT_ONE;
	two_pow_w.two_x(w);
	int two_w_word = 1;
	while (w >= BIGINT_WORD_SIZE * two_w_word)
		two_w_word++;
	two_w_word--;
	bool has_bit;
	while (current_mask.words[two_w_word] < two_pow_w.words[two_w_word])
	{
		is_independent = true;
		for (i = 0; i < w; i++)
		{
			// This needs to be changed if we sort the nbr_mask_vec
			if (current_mask.test_bit(T->tree_nodes[k]->nbr_mask_vec->at(i)->k))
			{
				// current_mask has the bit corresponding to nbr_mask_vec[i] set - so and with the appropriate vector for a hit
				has_bit = false;
				// Do the loop below manually rather than calling the overloaded AND
				for (j = 0; j < current_mask.S; j++)
				{
					if (current_mask.words[j] & T->tree_nodes[k]->nbr_mask_vec->at(i)->w->words[j])
					{
						// Some nbr is found in word j
						has_bit = true;
						break;
					}
				}

				if (has_bit)
				{
					is_independent = false;
#if 0
					// Find the rightmost bit in the current_mask's word where we had the intersection
					// with the edge mask and then add this to the current mask - a nice one liner!
					current_mask.words[j] += (current_mask.words[j] & ~(current_mask.words[j]-1));		
#else
					int rpos = 0;
					// Find the rightmost bit in current_mask - this could be sped up
					while (!(current_mask.test_bit(rpos)))
						rpos++;

					// Do the add w/ only bitwise op's
					// Would be faster here to figure out which word we are in and
					// do this w/o test_bit
					int b = 0;
					while (current_mask.test_bit(rpos + b))
						b++;

					for (int a = 0; a < b; a++)
						current_mask.xor_bit(rpos + a);

					current_mask.or_bit(rpos + b);
#endif
					// Break out of the for loop since current_mask cannot represent
					// an independent set
					break;
				}
			}
		}

		if (is_independent)
		{
			// new_set is the actual subset in this bag that is independent
			new_set = new TDSolution(T->num_mask_words);
			*(new_set->mask) = current_mask;
			// Set the orig_mask
			*(new_set->orig_mask) = current_mask;
			new_set->value = 0;
			for (j = 0; j < w; j++)
			{
				if (current_mask.test_bit(j))
					new_set->value
					+= T->G->get_vertex_weight(T->tree_nodes[k]->bag_vec[j]);
				//weight[T->tree_nodes[k]->bag_vec[j]];
			}
			// There is no point searching the hash table for nice leaves since every set
			// is added to the hash table
			HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, new_set->mask->words,
				T->num_mask_words*sizeof(BIGINT_WORD), new_set);
			T->info->total_pc_table_entries++;
			T->info->total_table_entries++;
			num_added++;
			if (new_set->value > T->tree_nodes[k]->obj_val)
				T->tree_nodes[k]->obj_val = new_set->value;
			// Advance to next candidate ind. set
			++current_mask;
		}
	}

	// Return the size of the table
	return num_added;
}

/**
* Computes the table for an introduce node (node has 1 child whose bag is smaller by 1)
* in a nice tree decomposition.
* Fatal error if the tree is not nice.
*/
int compute_introduce_table(TDTree *T, int k)
{
	if (!T->nice)
		fatal_error(
		"%s: only know how to compute table for introduce node in nice tree!\n",
		__FUNCTION__);
	if (T->tree_nodes[k]->adj.size() != 2)
		fatal_error("Tree node %d does not appear to be an introduce node!\n",
		k);

	int i, j, v = GD_UNDEFINED, v_pos;
	int child = GD_UNDEFINED;
	TDSolution *new_set = NULL, *new_v_set = NULL;
	list<int>::iterator ii, jj, kk;

	list<bigint_t>::iterator current_mask_it;
	vector<int> I(T->G->get_capacity());
	vector<int>::iterator vv;
	list<TDSolution *>::iterator ss;

	// First find the parents and children and the difference in the bags
	//parent = T->tree_nodes[k]->adj.front();
	child = T->tree_nodes[k]->adj.back();

	if (T->tree_nodes[k]->bag.size() <= T->tree_nodes[child]->bag.size())
		fatal_error(
		"%s: Node %d has a smaller or same size bag than its child --> not an introduce node\n",
		__FUNCTION__, k);
	// Find v
	vv = set_difference(T->tree_nodes[k]->bag.begin(),
		T->tree_nodes[k]->bag.end(), T->tree_nodes[child]->bag.begin(),
		T->tree_nodes[child]->bag.end(), I.begin());
	if (int(vv - I.begin()) != 1)
		fatal_error("%s:  Child/parent differ by something other than 1!\n",
		__FUNCTION__);
	v = *(I.begin());
	// Find v_pos, the position of the node that differentiates parent and child
	// and create a mask that has a 1 in position k when bag[i]-v is an edge in the graph
	Graph::GraphProperties properties;
	int current_key = properties.fill_adj_vec(T->G, v);
	//int current_key=T->G->fill_adj_vec(v);
	bigint_t v_nbr_mask(T->num_mask_words);
	i = 0;
	v_pos = GD_UNDEFINED;
	vector<int> *adj_vec = T->G->get_adj_vec_ptr();
	for (ii = T->tree_nodes[k]->bag.begin(); ii != T->tree_nodes[k]->bag.end(); ++ii)
	{
		if ((*ii) == v)
			v_pos = i;
		if (adj_vec->at(*ii) == current_key)
			v_nbr_mask.set_bit(i);
		i++;
	}
	if (v_pos == GD_UNDEFINED)
		fatal_error("%s:  Didn't find v=%d!\n", __FUNCTION__, v);

	T->tree_nodes[k]->obj_val = -1;
	bool is_independent;
	TDSolution *tt, *tmp;
	int num_added = 0;

	HASH_ITER(hh,T->tree_nodes[child]->hash_table, tt, tmp)
	{
		// CSG, 8/16/2011 - this is really wasteful. Maybe this could be fixed by
		// using HASH_ADD_PTR and then just modifying the data/mask that is pointed to??
		// Not sure but this is less than ideal to just repeatedly copy the data

		// Loop over every set tt in the child's hash table
		// Copy tt into new_set
		new_set = new TDSolution(T->num_mask_words);
		*(new_set->mask) = *(tt->mask);
		new_set->mask->insert_zero_bit(v_pos);
		*(new_set->orig_mask) = *(new_set->mask);
		new_set->value = tt->value;
		HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, new_set->mask->words,
			T->num_mask_words*sizeof(BIGINT_WORD), new_set);
		if (new_set->value > T->tree_nodes[k]->obj_val)
			T->tree_nodes[k]->obj_val = new_set->value;
		num_added++;
		T->info->total_table_entries++;
		// Now consider adding vertex v to this set - need to construct the set from
		// the mask and the bag
		is_independent = true;
		for (j = 0; j < new_set->mask->S; j++)
		{
			if (new_set->mask->words[j] & v_nbr_mask.words[j])
			{
				// Some nbr is found in word j
				is_independent = false;
				break;
			}
		}
		if (is_independent)
		{
			T->info->total_table_entries++;
			// The set tt \cup {v} is independent
			new_v_set = new TDSolution(T->num_mask_words);
			*(new_v_set->mask) = *(new_set->mask);
			// set the bit for v!
			new_v_set->mask->set_bit(v_pos);
			new_v_set->value = new_set->value + T->G->get_vertex_weight(v);
			//T->G->weight[v];
			*(new_v_set->orig_mask) = *(new_v_set->mask);
			// Add to the hash table
			HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, new_v_set->mask->words,
				T->num_mask_words*sizeof(BIGINT_WORD), new_v_set);
			num_added++;
			if (new_v_set->value > T->tree_nodes[k]->obj_val)
				T->tree_nodes[k]->obj_val = new_v_set->value;
		}
	}

	return num_added;
}

/**
* Computes the table for a forget node (node has 1 child whose bag is larger by 1) in a 
* nice tree decomposition.  Fatal error if the tree is not nice.
*/
int compute_forget_table(TDTree *T, int k)
{
	if (!T->nice)
		fatal_error(
		"%s: only know how to compute table for forget node in nice tree!\n",
		__FUNCTION__);

	if (T->tree_nodes[k]->adj.size() != 2)
	{
		cerr << *(T->tree_nodes[k]);
		fatal_error("Tree node %d does not appear to be a forget node!\n", k);
	}
	int i, v = GD_UNDEFINED;
	int child = GD_UNDEFINED;
	list<int>::iterator ii, jj, kk;
	vector<int> I(T->G->get_capacity());
	vector<int>::iterator vv;
	list<TDSolution *>::iterator ss;
	// First find the parents and children and the difference in the bags
	//parent = T->tree_nodes[k]->adj.front();
	child = T->tree_nodes[k]->adj.back();

	// Compute the child difference
	if (T->tree_nodes[k]->bag.size() >= T->tree_nodes[child]->bag.size())
	{
		print_message(0, "Current tree node %d\n", k);
		cerr << *(T->tree_nodes[k]);
		print_message(0, "Child tree node %d\n", child);
		cerr << *(T->tree_nodes[child]);
		fatal_error(
			"Node %d has a larger or same size bag than its child --> not a forget node\n",
			k);
	}

	// Find v
	vv = set_difference(T->tree_nodes[child]->bag.begin(),
		T->tree_nodes[child]->bag.end(), T->tree_nodes[k]->bag.begin(),
		T->tree_nodes[k]->bag.end(), I.begin());
	if (int(vv - I.begin()) != 1)
		fatal_error("%s:  Child/parent differ by something other than 1!\n",
		__FUNCTION__);
	v = *(I.begin());
	// v is the differentiating node
	// Find v_pos
	int v_pos = GD_UNDEFINED;
	i = 0;
	for (ii = T->tree_nodes[child]->bag.begin(); ii
		!= T->tree_nodes[child]->bag.end(); ++ii)
	{
		if ((*ii) == v)
		{
			v_pos = i;
			break;
		}
		i++;
	}
	// Sanity check
	if (v_pos == GD_UNDEFINED)
		fatal_error("%s:  Didn't find v=%d!\n", __FUNCTION__, v);

	TDSolution *hash_ptr = NULL;
	bigint_t Scupv(T->num_mask_words);
	TDSolution *new_sol = NULL;
	T->tree_nodes[k]->obj_val = -1;

	TDSolution *tt, *tmp;
	int num_added = 0;

	HASH_ITER(hh,T->tree_nodes[child]->hash_table, tt, tmp)
	{
		// Loop over every set tt in the child's hash table
		// Check if the bit for v is set
		if (!(tt->mask->test_bit(v_pos)))
		{
			// The set tt (aka S) does not contain v, so this is an ind. set in the parent
			// Look for S \cup v in the child
			// Then parent set's value will be the max table value of S or S\cup v
			Scupv = *(tt->mask);
			Scupv.set_bit(v_pos);
			// Now look for Scupv in the table of the child
			hash_ptr = NULL;
			HASH_FIND(hh,T->tree_nodes[child]->hash_table,Scupv.words,
				T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
			int new_val = tt->value;
			if (hash_ptr)
			{
				// both S and S cup v are in the child - pick the max to stick in the parent's table
				if (hash_ptr->value > new_val)
					new_val = hash_ptr->value;
				// else we don't touch the value
			}
			new_sol = new TDSolution(T->num_mask_words);
			new_sol->value = new_val;
			*(new_sol->mask) = *(tt->mask);
			new_sol->mask->remove_bit(v_pos);
			*(new_sol->orig_mask) = *(new_sol->mask);
			T->info->total_table_entries++;
			HASH_ADD_KEYPTR(hh,T->tree_nodes[k]->hash_table,new_sol->mask->words,
				T->num_mask_words*sizeof(BIGINT_WORD), new_sol);
			num_added++;
			if (new_sol->value > T->tree_nodes[k]->obj_val)
				T->tree_nodes[k]->obj_val = new_sol->value;
		}
	}
	print_message(1, "FORGET done\n");
	return num_added;
}

/**
* Computes the table for a join node (node has 2 children with identical bags) in a nice tree decomposition.
* Fatal error if the tree is not nice.
*/
int compute_join_table(TDTree *T, int k)
{
	if (!T->nice)
		fatal_error(
		"%s: only know how to compute table for join node in nice tree!\n",
		__FUNCTION__);

	if (T->tree_nodes[k]->adj.size() != 3)
		fatal_error("Tree node %d does not appear to be a join node!\n", k);

	int j, child1 = GD_UNDEFINED, child2 = GD_UNDEFINED;
	list<int>::iterator ii;

	list<TDSolution *>::iterator ss;

	// First get the 2 children
	ii = T->tree_nodes[k]->adj.begin();
	++ii;
	child1 = *ii;
	++ii;
	child2 = *ii;
	// Make sure bag sizes are all the same as a sanity check
	int w = T->tree_nodes[k]->bag.size();
	if ((int) T->tree_nodes[child1]->bag.size() != w
		|| (int) T->tree_nodes[child2]->bag.size() != w)
		fatal_error("%s:  child bags not the same size as parents!\n");

	// Go through child1 and look up in child2's hash table
	// We only need to worry about the intersection since a subset is assumed
	// to have cost = -\infty if it is not in the list!
	T->tree_nodes[k]->obj_val = -1;
	TDSolution *hash_ptr, *tt, *tmp, *new_sol = NULL;
	int num_added = 0;

	HASH_ITER(hh,T->tree_nodes[child1]->hash_table, tt, tmp)
	{
		// If free_children is set, then this table should be "stolen" and not
		// just duplicated
		// Loop over every set tt in child1 - note that child2 will have the exact same
		// sets in his hash table, just with potentially different values (I think)
		// Find tt in child2 - it has to live there
		hash_ptr = NULL;
		HASH_FIND(hh,T->tree_nodes[child2]->hash_table,tt->mask->words,
			T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
		if (!hash_ptr)
			// This is a fatal error - should not happen
			fatal_error("%s:  Didn't find set in child2=%d's table!!!!!\n",
			__FUNCTION__, child2);
		new_sol = new TDSolution(T->num_mask_words);
		*(new_sol->mask) = *(tt->mask);
		*(new_sol->orig_mask) = *(new_sol->mask);
		new_sol->value = tt->value + hash_ptr->value;

		for (j = 0; j < w; j++)
		{
			if (new_sol->mask->test_bit(j))
			{
				// Bit j is set in i - so bag_vec[j] is in this set --> subtract the weight
				new_sol->value -= T->G->get_vertex_weight(T->tree_nodes[k]->bag_vec[j]);
				//weight[T->tree_nodes[k]->bag_vec[j]];
			}
		}
		// Update k's hash table
		HASH_ADD_KEYPTR(hh,T->tree_nodes[k]->hash_table,new_sol->mask->words,
			T->num_mask_words*sizeof(BIGINT_WORD), new_sol);
		T->info->total_table_entries++;
		num_added++;
		if (new_sol->value > T->tree_nodes[k]->obj_val)
			T->tree_nodes[k]->obj_val = new_sol->value;
	}

	return num_added;
}

/**
* Computes the number of nodes and edges in the graph induced by the vertices in the
* intersection of each tree node's bag with its parent's bag. Stores results in
* stats vector that contains n and m for every node in the tree.
*/
int create_subgraph_stats(TDTree *T, vector<int_int> *stats)
{
	if((int)stats->size() < T->num_tree_nodes)
		fatal_error("%s:  stats array is too small\n",__FUNCTION__);

	vector<int> parent_intersection_vec(T->G->get_num_nodes());
	vector<int>::iterator vv,ww;
	list<int> parent_intersection;
	int parent;
	for(int i=0;i<=T->num_tree_nodes;i++)
	{
		if(T->tree_nodes[i])
		{
			// Create the intersection of i with his parent since this is what we actually
			// store in the DP tables
			parent=T->tree_nodes[i]->adj.front();
			parent_intersection.clear();
			print_message(10,"Parent of %d is %d\n",i,parent);
			if(parent!=i)
			{
				ww=set_intersection(T->tree_nodes[i]->bag.begin(), T->tree_nodes[i]->bag.end(),
					T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end(),
					parent_intersection_vec.begin());
				print_message(10,"intersection has %d entries\n",ww-parent_intersection_vec.begin());

				for(vv=parent_intersection_vec.begin(); vv!=ww;++vv)
					parent_intersection.push_back(*vv);
			}
			else
				parent_intersection=T->tree_nodes[i]->bag;

#if 1
			// Compute the stats on the subgraph that is the intersection of this bag with his
			// parent since that is what we need to store in the tables!
			stats->at(i).p1 = (int)parent_intersection.size();
			stats->at(i).p2 = T->G->get_num_edges_in_subgraph(&parent_intersection);
#else
			stats->at(i).p1 = (int)T->tree_nodes[i]->bag.size();
			// DP_construct_adj_matrix returns the # of edges in the subgraph
			stats->at(i).p2 = DP_construct_adj_matrix(T,i);
#endif
			print_message(10,"%d\n",stats->at(i).p2);
		}	
	}
	return 1;
}


/**
* A user-defined table_function that computes a table for
* the weighted independent set dynamic programming.
*/
int compute_weighted_ind_set_table(TDTree *T, int k)
{
	if (T->info->verbose)
	{
		if (T->num_tree_nodes_processed == 0)
			// print a header line for verbose output
			print_message(0,
			"id #left #bag #child DP_prep_time DP_run_time #table #sets max_obj #sg_edges memHW\n");
		print_message(0, "%04d %04d %03d %d ", k,
			T->num_tree_nodes - T->num_tree_nodes_processed,
			T->tree_nodes[k]->bag.size(), T->tree_nodes[k]->adj.size() - 1);
	}
	// start the timer
	T->info->start = clock();

	// Prepare the node for the DP
	DP_prepare_node(T, k);
	T->info->stop = clock();
	// print out the "prep time"
	if (T->info->verbose)
		fprintf(stderr, "%1.5f ",
		(double) (T->info->stop - T->info->start) / CLOCKS_PER_SEC);

	int table_size = GD_UNDEFINED;
	int node_type = T->get_node_type(k);
	switch (node_type)
	{
	case TD_NICE_LEAF_NODE:
		table_size = compute_ind_sets_nice(T, k);
		T->info->stop = clock();
		T->info->leaf_time += (double) (T->info->stop - T->info->start)
			/ CLOCKS_PER_SEC;
		break;
	case TD_NONNICE_LEAF_NODE:
		table_size = compute_nonnice_table(T, k);
		T->info->stop = clock();
		T->info->leaf_time += (double) (T->info->stop - T->info->start)
			/ CLOCKS_PER_SEC;
		break;
	case TD_NONLEAF_NODE:
		// The choice of algorithm to use is contained inside T's DP_info struct
		table_size = compute_nonnice_table(T, k);
		T->info->stop = clock();
		T->info->nonleaf_time += (double) (T->info->stop - T->info->start)
			/ CLOCKS_PER_SEC;
		break;
	case TD_INTRODUCE_NODE:
		table_size = compute_introduce_table(T, k);
		T->info->stop = clock();
		T->info->introduce_time += (double) (T->info->stop - T->info->start)
			/ CLOCKS_PER_SEC;
		break;
	case TD_FORGET_NODE:
		table_size = compute_forget_table(T, k);
		T->info->stop = clock();
		T->info->forget_time += (double) (T->info->stop - T->info->start)
			/ CLOCKS_PER_SEC;
		break;
	case TD_JOIN_NODE:
		table_size = compute_join_table(T, k);
		T->info->stop = clock();
		T->info->join_time += (double) (T->info->stop - T->info->start)
			/ CLOCKS_PER_SEC;
		break;
	default:
		fatal_error(
			"Unknown tree node type %d encountered when computing table of tree node %d??\n",
			node_type, k);
		break;
	}
	T->num_tree_nodes_processed++;
	DP_free_data(T, k);

	// Free the children memory if requested
	if (T->info->free_children)
		T->free_children(k);

	// Print results on the node if verbose
	if (T->info->verbose)
	{
		int x, y;
		if (T->nice)
			x = -1;
		else
			x = T->info->aux_info[k];
#if !(__CYGWIN__ || WIN32 || _WIN32)
		y = getHWmem();
#else
		y=-1;
#endif

		// PRINT: table_size num_sets(-1 if nice) max_weight num_secs num_subgraph_edges mem_usage(-1 if Windows)
		print_message(0, "%1.5f %06d %06d %05d %05d %d\n", 
			(double) (T->info->stop - T->info->start) / CLOCKS_PER_SEC,table_size, x,
			T->tree_nodes[k]->obj_val,

			T->tree_nodes[k]->num_subgraph_edges, y);

	}

	if (T->info->very_verbose)
	{
		GEN("%d:%d:%d:", k, table_size, T->info->aux_info[k]);
		list<int>::iterator ii;
		ii = T->tree_nodes[k]->adj.begin();
		GEN("%d:", *ii);
		++ii;
		for (; ii != T->tree_nodes[k]->adj.end(); ++ii)
			GEN("%d,", *ii);
		GEN("\n");
	}
	return table_size;
}

/**
* Writes a GMPL model file to the file using information from the graph G
* If solved, the solution is written to <DIMACS_file>.IND_SET.sol.
*/
void write_ind_set_model(const char *DIMACS_file, const char *model_file,
	Graph::VertexWeightedGraph *G)
{
	int i, j;
	char sol_file[200];

	FILE *out = fopen(model_file, "w");
	if (!out)
		fatal_error("Can't open %s for writing\n", model_file);

	sprintf(sol_file, "%s.MIP.WIS.sol", DIMACS_file);

	// Write out the model
	fprintf(out,
		"# Automatically generated model file for max independent set\n"
		"# Assumes 1-based node names with no holes!!\n"
		"# Original graph file is %s\n"
		"# Model section\n\n"
		"param num_nodes ;\n"
		"set Nodes := 1..num_nodes ;\n"
		"set E within Nodes cross Nodes ;\n"
		"param w{Nodes};\n"
		"var x{i in Nodes} binary;\n\n"
		// Change to minimization for archaic MPS reasons
		"minimize Cost: sum{i in Nodes} -x[i]*w[i];\n\n"
		"subject to Independent{i in Nodes, j in Nodes: (i,j) in E}:\n"
		"x[i]+x[j]<=1;\n\n"
		"solve;\n\n", DIMACS_file);
	fprintf(
		out,
		"printf \" \\nMaximum Independent Set: %%d\\n\", sum{i in Nodes} x[i]*w[i] > \"%s\";\n",
		sol_file);
	fprintf(out, "for { i in Nodes }\n"
		"{\n"
		"   printf{0..0: x[i]!=0} \"%%d (%%3.1f)\\n\", i, w[i] >> \"%s\";\n",
		sol_file);
	fprintf(out, "}\n\n");

	fprintf(out, "# Data Section\n\ndata;\n\n");
	// Write out the problem data
	fprintf(out, "param num_nodes:=%d;\n", G->get_num_nodes());
	fprintf(out, "param w:=\n");
	vector<Graph::Node> nodes = G->get_nodes();
	vector<int> weight = G->get_weight();

	for (i = 0; i < G->get_num_nodes(); i++)
	{
		fprintf(out, "%d %d\n", nodes[i].get_label(), weight[i]);
	}
	fprintf(out, ";\n");
	fprintf(out, "set E:=\n");
	Graph::GraphProperties properties;

	vector<int> adj_vec;
	for (i = 0; i < G->get_num_nodes(); i++)
	{
		// This assumes 1-based labeling!
		int current_key = properties.fill_adj_vec(G, i);
		//int current_key=G->fill_adj_vec(i);
		adj_vec = G->get_adj_vec();
		for (j = i + 1; j < G->get_num_nodes(); j++)
		{
			if (adj_vec[j] == current_key)
				fprintf(out, "(%d,%d)\n", i + 1, j + 1);
		}
		fprintf(out, "\n");
	}
	fprintf(out, ";\nend;\n");
	fclose(out);
}

/**
* Used to reconstruct the solution after dynamic programming.  This takes the provided 
* mask (assumed to live in tree node k's hash table),
* and computes the contribution of the corresponding set to the provided optimal_vec.
*/
void compute_contribution(TDTree *T, int k, bigint_t *mask,
	vector<bool> *optimal_vec)
{
	TDSolution *hash_ptr;
	list<int>::iterator ii;
	int i;

	if (T->nice)
	{
		// If the decomposition is nice, then the mask just picks off nodes in k's bag
		for (i = 0; i < T->num_mask_words * BIGINT_WORD_SIZE; i++)
		{
			if (mask->test_bit(i))
				optimal_vec->at(T->tree_nodes[k]->bag_vec[i]) = true;
		}
		return;
	}
	else
	{
		// Non-nice
		// Find the mask in k's table - the mask is in the language of k's parent if non-nice
		HASH_FIND(hh,T->tree_nodes[k]->hash_table, mask->words, T->num_mask_words*sizeof(BIGINT_WORD),
			hash_ptr);
		if (!hash_ptr)
		{
			// This is a fatal errorint
			print_message(
				0,
				"%s:  Didn't find mask in table of tree node %d\n\tMissing mask: ",
				__FUNCTION__, k);
			mask->print(stderr);
			list<int> missing_set;
			int parent = T->tree_nodes[k]->adj.front();
			vector<int> parent_vec(T->tree_nodes[parent]->bag.size(), -1);
			i = 0;
			for (ii = T->tree_nodes[parent]->bag.begin(); ii
				!= T->tree_nodes[parent]->bag.end(); ++ii)
			{
				parent_vec[i] = *ii;
				i++;
			}
			for (i = 0; i < T->num_mask_words * BIGINT_WORD_SIZE; i++)
				if (mask->test_bit(i))
					missing_set.push_back(parent_vec[i]);
			print_message(0, "Missing set is:\n");
			print(0, missing_set);
			int val;
			Graph::GraphProperties properties;

			if (properties.is_independent_set(
				(Graph::VertexWeightedGraph *) (T->G), &missing_set, &val))
				print_message(0, "Missing set is independent with weight %f\n",
				val);
			else
				print_message(0, "Missing set is NOT independent!!!\n");
			fatal_error("Cannot continue.\n");
		}

		// Same computation if tree is/isn't nice
		// Now what we care about is the orig_mask of hash_ptr since that contains
		// the other nodes from k's bag in addition to those found in the intersection
		// with the parent (that set is represented by mask in the language of the parent)
		for (i = 0; i < T->num_mask_words * BIGINT_WORD_SIZE; i++)
		{
			if (hash_ptr->orig_mask->test_bit(i))
				// bit i is set in the orig_mask, so bag_vec[i] is in this set
				// Set position optimal_vec[bag_vec[i]]=true;
				optimal_vec->at(T->tree_nodes[k]->bag_vec[i]) = true;
		}
	}
	return;
}

/**
* Function used to reconstruct solutions. The mask points to a bigint_t that should live 
* in the table of tree node k. If the TD is non-nice, then this function looks at the children
* of tree node k and computes mask_child=mask\cap bag(child) and sets the appropriate pointer in the
* masks_to_process vector for each child.  Then we know to look for the
* mask_child in the child's bag at the
* next step in the reconstruction.
* If the TD is nice, then the mask represents a solution in terms of k's bag, and
* I set the child masks by doing a case-by-case consideration of the TDNode type (forget, introduce, join)
*/
void record_children_solutions(TDTree *T, int k, bigint_t *mask,
	vector<bigint_t *> *masks_to_process)
{
	list<int>::iterator ii, jj;

	// kset will be the actual subset of nodes from k's bag
	list<int> kset;
	vector<int> intersection(T->tree_nodes[k]->bag.size(), -1);
	vector<int>::iterator vv;
	vector<int> pos_vec(T->G->get_capacity() + 1, -1);
	// pos_vec[j]=i means that j is the i'th entry in k's bag

	TDSolution *hash_ptr;
	bigint_t *new_mask, *new_mask2;
	int i = 0;

	if (T->nice)
	{
		// T is nice - mask is in the language of k
		// Need to figure out what mask to look for in the children based on the
		// relationship between k's bag and his children
		ii = T->tree_nodes[k]->adj.begin();
		++ii; // advance to first child

		// Case 1 - join node where we have 2 children whose bags are identical to k's
		if (T->tree_nodes[*ii]->bag.size() == T->tree_nodes[k]->bag.size())
		{
			new_mask = new bigint_t(T->num_mask_words);
			*new_mask = *mask;
			masks_to_process->at(*ii) = new_mask;
			// Same thing for child 2
			++ii;
			new_mask2 = new bigint_t(T->num_mask_words);
			*new_mask2 = *mask;
			masks_to_process->at(*ii) = new_mask2;
			return;
		}

		// Case 2 - child has one less node than k
		if (T->tree_nodes[*ii]->bag.size() < T->tree_nodes[k]->bag.size())
		{
			// find the difference - k's bag contains one node that ii's doesn't
			i = 0;
			list<int>::iterator uu, vv;
			uu = T->tree_nodes[k]->bag.begin();
			vv = T->tree_nodes[*ii]->bag.begin();
			while (*uu == *vv)
			{
				++uu;
				++vv;
				i++;
			}
			int v_pos = i;
			// It doesn't matter if bit v_pos is set or not - just remove it to create the
			// mask for the child!
			new_mask = new bigint_t(T->num_mask_words);
			*new_mask = *mask;
			new_mask->remove_bit(v_pos);
			masks_to_process->at(*ii) = new_mask;
			return;
		}

		// Case 3 - child has one more node than k
		if (T->tree_nodes[*ii]->bag.size() > T->tree_nodes[k]->bag.size())
		{
			// find the difference - child bag contains one node that k doesn't
			i = 0;
			list<int>::iterator uu, vv;
			uu = T->tree_nodes[k]->bag.begin();
			vv = T->tree_nodes[*ii]->bag.begin();
			while (*uu == *vv)
			{
				++uu;
				++vv;
				i++;
			}
			int v_pos = i;
			new_mask = new bigint_t(T->num_mask_words);
			*new_mask = *mask;
			new_mask->insert_zero_bit(v_pos);
			// The mask in k was derived either from the same set in the child or by setting bit v_pos
			new_mask->set_bit(v_pos);
			HASH_FIND(hh,T->tree_nodes[*ii]->hash_table, new_mask->words,T->num_mask_words*sizeof(BIGINT_WORD),
				hash_ptr);
			int val = -1;
			if (hash_ptr)
				val = hash_ptr->value;
			// Now clear v_pos and look
			new_mask->xor_bit(v_pos);
			HASH_FIND(hh,T->tree_nodes[*ii]->hash_table, new_mask->words,T->num_mask_words*sizeof(BIGINT_WORD),
				hash_ptr);
			if (!hash_ptr)
			{
				// This is an error for WIS
				// May not be an error for other DP's, not sure...
				mask->print(stderr);
				new_mask->print(stderr);
				fatal_error("%s:  Didn't find mask with zeroized bit??\n",
					__FUNCTION__);
			}
			if (hash_ptr->value > val)
				masks_to_process->at(*ii) = new_mask;
			else
			{
				// Better off with the v_pos bit set
				new_mask->set_bit(v_pos);
				masks_to_process->at(*ii) = new_mask;
			}
			return;
		}
	}
	else
	{
		// T is not nice
		// mask is in the language of k's parent - need to find the corresponding orig_mask in k!!
		HASH_FIND(hh,T->tree_nodes[k]->hash_table, mask->words,T->num_mask_words*sizeof(BIGINT_WORD),
			hash_ptr);
		if (!hash_ptr)
			fatal_error("%s:  couldn't find mask???\n", __FUNCTION__);

		i = 0;
		for (ii = T->tree_nodes[k]->bag.begin(); ii
			!= T->tree_nodes[k]->bag.end(); ++ii)
		{
			// Use the orig_mask as that references nodes in k's bag
			if (hash_ptr->orig_mask->test_bit(i))
				kset.push_back(*ii);
			pos_vec[*ii] = i;
			i++;
		}
		// kset is now a list of the vertices in k's parent's bag represented by the mask

		ii = T->tree_nodes[k]->adj.begin();
		++ii; // advance to first child
		for (; ii != T->tree_nodes[k]->adj.end(); ++ii)
		{
			// Intersect the child bag with kset
			vv = set_intersection(T->tree_nodes[*ii]->bag.begin(),
				T->tree_nodes[*ii]->bag.end(), kset.begin(), kset.end(),
				intersection.begin());

			// Create the child mask for the intersection - since we want to look this mask up,
			// it is expressed in the language of the parent bag (k in this case)
			new_mask = new bigint_t(T->num_mask_words);
			for (i = 0; i < vv - intersection.begin(); i++)
				new_mask->set_bit(pos_vec[intersection[i]]);
			masks_to_process->at(*ii) = new_mask;
		}
	}

	return;
}

/**
* Function to reconstruct the optimal solution after reaching the root node in a DP by descending
* back down the tree.  If sol_file is non-NULL, solution is written to sol_file
*/
int reconstruct_solution(TDTree *T, list<int> *optimal_solution,
	const char *sol_file)
{
	// Get the best solution from the root
	TDSolution *ss, *temp, *best_sol = NULL;

	int i;
	int max_val = -1;
	// Record the optimal solution here
	vector<bool> optimal_vec(T->G->get_capacity(), false);

	// Compute a pre-order walk of the tree
	vector<int> walk(T->tree_nodes.size());
	T->pre_order_walk(&walk);
	if (walk[0] != T->root_node)
		fatal_error("walk[0]!=root??? %d!=%d???\n", walk[0], T->root_node);

	// Find the best solution in the root_node table
	HASH_ITER(hh,T->tree_nodes[T->root_node]->hash_table,ss,temp)
	{
		if (ss->value > max_val)
		{
			max_val = ss->value;
			best_sol = ss;
		}
	}
	// Note that mask==orig_mask in the root since he is his own parent

	// The masks_to_process vector contains an entry for each node in the Tree.
	// Each entry is a pointer to a mask that represents the contribution from
	// this tree node's bag to the optimal set found in the root node.
	// If the TD is non-nice, then the masks are written in terms of the tree
	// node's parent's bag since this
	// is the key that we use to look up things in the hash table
	// If the TD is nice, then masks_to_process[k] references nodes in k's bag
	vector<bigint_t *> masks_to_process(T->tree_nodes.size(), NULL);

	// Set the mask for the root node
	masks_to_process[T->root_node] = best_sol->mask;

	// Print out the # of vertices in best_sol versus the # in that bag
	int num_bits = 0;
	for (i = 0; i < T->width + 1; i++)
		if (best_sol->mask->test_bit(i))
			num_bits++;
	if(T->info->verbose)
	{
		print_message(0, "Best solution %d has %d/%d nodes from root bag\n",
			max_val, num_bits, (int) T->tree_nodes[T->root_node]->bag.size());
	}
	// Compute the contribution of the best sol in the root to the optimal solution
	compute_contribution(T, walk[0], masks_to_process[walk[0]], &optimal_vec);
	int max_tree_nodes=(int)T->tree_nodes.size();
	for (i = 0; i < T->num_tree_nodes; i++)
	{
		// Look at the solution in masks_to_process[walk[i]] and record the relevant masks
		// for walk[i]'s children in the masks_to_process array
		// Only do this if node walk[i] is a non-leaf
		if (T->tree_nodes[walk[i]]->adj.size() > 1)
			record_children_solutions(T, walk[i], masks_to_process[walk[i]],
			&masks_to_process);

		// Now compute the contribution of the children whose masks we just added
		list<int>::iterator ii = T->tree_nodes[walk[i]]->adj.begin();
		++ii;
		// Advance to the first child
		for (; ii != T->tree_nodes[walk[i]]->adj.end(); ++ii)
			// Compute contribution of each child
			compute_contribution(T, *ii, masks_to_process[*ii], &optimal_vec);
	}

	// Now record the solution in the list
	int optimal_obj_val = 0;
	for (i = 0; i < T->G->get_capacity(); i++)
	{
		if (optimal_vec[i])
		{
			optimal_solution->push_back(i);
			optimal_obj_val += T->G->get_vertex_weight(i);
			//weight[i];
		}
	}

	// Delete the masks_to_process entries
	for (i = 0; i < (int) T->tree_nodes.size(); i++)
	{
		if (masks_to_process[i] && i != T->root_node)
		{
			delete masks_to_process[i];
			masks_to_process[i] = NULL;
		}
	}
	// Check the optimal solution
	Graph::GraphProperties properties;
	if (!properties.is_independent_set(T->G, optimal_solution))
		fprintf(stderr, "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"Claimed optimal solution is not an ind. set!\n"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");

	Graph::Node *n1;
	if (sol_file)
	{
		FILE *out = fopen(sol_file, "w");
		if (!out)
		{
			fprintf(stderr, "Can't write to sol_file %s\n", sol_file);
			exit(-1);
		}
		fprintf(out, "# Optimal solution to %s has obj. function value %d\n",
			T->graph_file, optimal_obj_val);
		fprintf(out,
			"# Solution written in terms of original labels from input file\n");
		for (list<int>::iterator ii = optimal_solution->begin(); ii
			!= optimal_solution->end(); ++ii)
		{
			n1=T->G->get_node(*ii);
			fprintf(out, "%d %d\n",n1->get_label(),T->G->get_vertex_weight(*ii));
			//	T->G->nodes[*ii].label, T->G->weight[*ii]);
		}
		fclose(out);

	}
	if(T->info->verbose)
	{
		fprintf(stderr, "Optimal solution (%d) written to %s\n", optimal_obj_val,
			sol_file);
	}
	return optimal_obj_val;
}

/**
* Populates the provided Graph structure with the necessary information to 
* do the WIS calculation.
*/
void create_WIS_graph(DP_info *info, Graph::VertexWeightedGraph *&G)
{
	int i;

	Graph::GraphCreatorFile creator;
	creator.set_file_name(info->DIMACS_file);
	creator.set_graph_type("DIMACS");

	G = creator.create_weighted_mutable_graph();
	if (!G)
	{
		return;
	}

	// Add weights of 1 to all vertices if we didn't load any from the file
	bool has_weights = false;
	vector<int> weight = (G)->get_weight();
	for (i = 0; i < (G)->get_capacity(); i++)
	{
		if (weight[i] != 0)
		{
			has_weights = true;
			break;
		}
	}

	if (!has_weights)
		for (i = 0; i < (G)->get_capacity(); i++)
			weight[i] = 1;

	Graph::GraphProperties properties;
	Graph::GraphEOUtil eoutil;
	Graph::GraphWriter *writer;
	Graph::GraphReaderWriterFactory factory;

	writer = factory.create_writer("DIMACS");
	// Put the input graph in canonical form for the tests
	properties.make_canonical(G);
	// Make sure there is only 1 component in the graph!
	if (!properties.is_connected(G))
	{
		print_message(0,
			"WIS only runs on connected graphs.  This graph has more than one component!\n");
		list<Graph::VertexWeightedGraph *> C;

		C = creator.create_all_components(G, true);
		print_message(0, "Found %d components\n", C.size());
		list<Graph::VertexWeightedGraph *>::iterator cc;
		int comp_number = 0;
		char temp_file[100], comp_file[100];
		int max_size = -1;
		char *big_file = NULL;
		for (cc = C.begin(); cc != C.end(); ++cc)
		{
			if ((*cc)->get_num_nodes() > 1 && (*cc)->get_num_edges() > 1)
			{
				sprintf(temp_file, "%s_%d_node_component_%d.temp",
					info->DIMACS_file, (*cc)->get_num_nodes(), comp_number);
				sprintf(comp_file, "%s_%d_node_component_%d.comp",
					info->DIMACS_file, (*cc)->get_num_nodes(), comp_number);
				writer->set_out_file_name(temp_file);
				writer->write_graph(*cc);

				normalize_DIMACS_file(temp_file, comp_file);
				remove(temp_file);
				fprintf(stderr, "Wrote component with %d nodes to file %s\n",
					(*cc)->get_num_nodes(), comp_file);
				comp_number++;
				if ((*cc)->get_num_nodes() > max_size)
				{
					if (!big_file)
						big_file = (char *) malloc(100);

					max_size = (*cc)->get_num_nodes();
					sprintf(big_file, "%s", comp_file);
				}
			}
			else
			{
				fprintf(
					stderr,
					"Tiny component with %d nodes & %d edges not written\n",
					(*cc)->get_num_nodes(), (*cc)->get_num_edges());
			}
		}


		if (big_file)
		{
			info->DIMACS_file = (char *) malloc(100);//3-21-2012 - before, we were overflowing the char*.
			sprintf(info->DIMACS_file, "%s", big_file);
			free(big_file);
			print_message(0,
				"Will run WIS on largest component %s with %d nodes\n",
				info->DIMACS_file, max_size);
			delete G;
			create_WIS_graph(info, G);
		}
		else
		{
			print_message(0, "No big file found. Exiting\n");
			exit(-1);
		}

		int s = C.size();
		for (i = 0; i < s; i++)
		{
			delete C.front();
			C.pop_front();
		}
	}

	delete writer;

}


/**
* Creates a tree decomposition for the provided tree and the DP options
* contained in DP_info.
*/
void create_tree_decomposition(DP_info *info, Graph::VertexWeightedGraph *G,
	TDTree **T)
{
	create_tree_decomposition(info, G, T, true);
}

/**
* Creates a tree decomposition for the provided tree and the DP options
* contained in DP_info.
*/
void create_tree_decomposition(DP_info *info, Graph::VertexWeightedGraph *G,
	TDTree **T, bool suppress_timing)
{
	// Create a copy of G to do the triangulation and elim order
	Graph::GraphEOUtil eoutil;

	if (!G)
	{
		return;
	}

	Graph::VertexWeightedGraph H = *G;

	//Run (optional) lower bound heuristics and print information
	if(info->lower_bounds)
	{
		int lb; 

		lb = eoutil.get_tw_lower_bound(G, GD_MAX_MIN_DEGREE_LB, 0); 
		printf("%s: MMD Lower Bound %d\n", info->DIMACS_file, lb);
		lb = eoutil.get_tw_lower_bound(G, GD_MCS_LB, 0); 
		printf("%s: MCS Lower Bound %d\n", info->DIMACS_file, lb);

	}

	(*T) = new TDTree(&H);
	// Set the name of T's graph file
	sprintf((*T)->graph_file, "%s", info->DIMACS_file);

	// create a vector for the elimination ordering
	vector<int> ordering(H.get_num_nodes(), GD_UNDEFINED);
	int i;

	// Set the info pointer so we have access to the options
	(*T)->info = info;

	//	creating an instance of eoutils class

	// Now read in the tree file if we have one
	if (info->read_tree)
	{
		// Read in the decomposition from a file
		(*T)->width = (*T)->read_DIMACS_file(info->tree_infile);
		// Assume not nice - could check this in read_DIMACS_file actually...
		(*T)->nice = false;
	}
	else
	{
		// Figure out how to create the tree
		if (info->read_ordering)
		{
			read_ordering_file(info->ord_file, &ordering);
			print_message(0, "Read in ordering\n");
			print(1, ordering);
		}
		else
		{
			if (info->read_scotch_ordering)
			{
				read_SCOTCH_ordering_file(info->scotch_ord_file, &ordering);
				print_message(0, "Read in SCOTCH ordering\n");
				print(1, ordering);
			}
			else
			{
				// Create the ordering via a heuristic
				// Create an ordering - if start_v not provided, find
				// a good candidate
				double start = clock();
				if (!info->parmetis)
				{
					if (info->start_v == GD_UNDEFINED)
						eoutil.find_elimination_ordering(&H, &ordering,
						info->elim_order_type, false);
					else
						eoutil.find_elimination_ordering(&H, &ordering,
						info->elim_order_type, info->start_v, false);
				}
				else
				{
#ifdef HAS_PARMETIS
					MPI_Comm comm;
					MPI_Comm_dup(MPI_COMM_WORLD, &comm);
					eoutil.parmetis_elimination_ordering(&H, ordering, info->elim_order_type, false, comm);
#endif
				}

				if(!suppress_timing)
					print_message(0, "%.2f:", (double) (clock() - start) / CLOCKS_PER_SEC);
			}

			// Write elimination ordering into a file
			if (info->eorder)
			{
				int os = ordering.size();
				FILE *ef = fopen(info->eorder, "w");
				if (!ef)
				{
					FERROR("can not open %s for writing\n", info->eorder);
					return;
				}

				for (int i = 0; i < os; i++)
				{
					fprintf (ef, "%d\n", ordering[i] + 1);
				}
				fclose (ef);
			}
		}

		if (info->very_verbose)
		{
			print_message(0, "\nORDERING:\n");
			print(0, ordering);
			print_message(0, "\n\n");
		}

		// Triangulate the graph for methods requiring it.
		clock_t tri_start = clock(), tri_stop;
		if(info->superetree){
			//no need to triangulate; width set in construction.
			tri_stop = tri_start;	
		}
		else
		{
#if HAS_METIS
			(*T)->width = eoutil.METIS_triangulate(&H, &ordering);
#else
			(*T)->width = eoutil.triangulate(&H, &ordering);
#endif
			tri_stop = clock();

		}
		if(!suppress_timing)
			print_message(0, "%.2f:", (double) (tri_stop - tri_start) / CLOCKS_PER_SEC);


		// Now create the tree
		info->start=clock();
		if(info->superetree)
		{
			(*T)->construct_superetree(&ordering);
		}
		if (info->gavril)
		{
			(*T)->construct_gavril(&ordering);
		}

		if (info->BK)
		{
			(*T)->construct_BK(&ordering);
		}

		if (info->nice)
		{
			print_message(1, "nice\n");
			(*T)->construct_knice(&ordering, (*T)->width, false);
		}
		info->stop=clock();

		//moved information about triangulation and width below construction since width in superetree is unknown until here.
		if (info->verbose)
		{
			print_message(0, "Triangulation took %f secs\n",
				(double) (tri_stop - tri_start) / CLOCKS_PER_SEC);
			print_message(0, "Width=%d\n", (*T)->width);
		}

		if (info->width)
		{
			DEBUG("%s: Width=%d\n", info->DIMACS_file, (*T)->width);
			return;
		}


		if(!suppress_timing)
			print_message(0, "%.2f:", (double) (info->stop - info->start) / CLOCKS_PER_SEC);

		if (info->verbose)
		{
			print_message(0, "Tree construction took %f secs\n",
				(double) (info->stop - info->start) / CLOCKS_PER_SEC);
		}
	}



	if (!info->nice)
	{
		// Root the tree
		if (!info->asp_root && !info->child_root)
			(*T)->root(info->root_node);
		else
		{
			// A rooting alg was chosen
			if (info->asp_root)
				(*T)->root_mintables(TD_ROOT_ASPVALL);
			if (info->child_root)
				(*T)->root_mintables(TD_ROOT_ALLCHILDREN);
		}
	}

	// Sort the bags
	int num_tree_nodes=(int) (*T)->tree_nodes.size();
	if(!info->superetree)//bags already sorted
	{
		for (i = 0; i < num_tree_nodes; i++)
		{
			if ((*T)->tree_nodes[i])
				(*T)->tree_nodes[i]->bag.sort();
		}
	}
	// Reset (*T)'s graph is the original, non-triangulated graph!
	(*T)->G = G;

	// Make the tree nice now if requested
	if (info->make_nice)
		(*T)->make_nice();

	// Refine the tree decomposition if requested
	if (info->refine_td)
	{
		if (1) // for testing
		{
			(*T)->write_DIMACS_file("pre_refinement.td");
			if (info->make_histogram)
			{
				vector<int> counts((*T)->width + 2, 0);
				for (i = 0; i < (int) (*T)->tree_nodes.size(); i++)
				{
					if ((*T)->tree_nodes[i] != NULL)
						counts[(*T)->tree_nodes[i]->bag.size()]++;
				}
				printf("Histogram of bag sizes\n");
				for (i = 1; i <= (*T)->width + 1; i++)
				{
					printf("%02d ", counts[i]);
				}
				printf("\n");
				fflush(stdout);
			}
		}
		(*T)->refine();
	}

	// Write out the tree if desired
	// Write out the tree if desired with root labelled 1, and rest
	//labelled according to pre-order walk
	if (info->write_ordered_tree)
		(*T)->write_DIMACS_file(info->tree_outfile, true);
	else if (info->write_tree)
		(*T)->write_DIMACS_file(info->tree_outfile);

	// Print out a histogram if desired
	if (info->make_histogram)
	{
		vector<int> counts((*T)->width + 2, 0);
		for (i = 0; i < (int) (*T)->tree_nodes.size(); i++)
		{
			if ((*T)->tree_nodes[i] != NULL)
				counts[(*T)->tree_nodes[i]->bag.size()]++;
		}
		printf("Histogram of bag sizes\n");
		for (i = 1; i <= (*T)->width + 1; i++)
		{
			printf("%02d ", counts[i]);
		}
		printf("\n");
		fflush(stdout);
	}

	if (info->gviz)
	{
		char dimacs_tree_file[100];
		// Create a non-labeled Graphviz output
		sprintf(dimacs_tree_file, "%s.dimacs", info->gviz_file);
		(*T)->write_DIMACS_file(dimacs_tree_file);
		//(*T)->write_graphviz_file(true, info->gviz_file, GV_BAG_LABELS);
		(*T)->write_graphviz_file(true, info->gviz_file, GV_COLORS);
	}

	// Verify the tree if this option was passed in
	if (info->check)
	{
		if ((*T)->nice)
		{
			print_message(0, "Verifying tree\n");
			if ((*T)->is_nice())
			{
				// this calls verify()
				printf("Nice TD is verified\n");
			}
			else
				fatal_error("Nice tree invalid!\n");
		}
		else
		{
			// non-nice
			if ((*T)->verify())
			{
				printf("Non-Nice TD is verified\n");
			}
			else
				fatal_error("Non-Nice tree invalid!\n");
		}
	}

	if (info->nonniceDP)
	{
		print_message(0, "Forcing non-nice DP!\n");
		(*T)->nice = false;
	}

	// Determine how big the masks need to be based on the tree width
	// CSG changing July 6!
	// Width is 1 less than the largest bag! Then we need an extra word
	// to handle a shift of BIGINT_WORD_SIZE!!
	(*T)->num_mask_words = (int) ceil(
		((double) (*T)->width + 2) / BIGINT_WORD_SIZE);
	print_message(1, "width=%d; num_words=%d\n", (*T)->width,
		(*T)->num_mask_words);

	// set T->num_mask_word in nodes, num_mask_words need to copy data such as nbr_mask_vec,
	// child_intersection and parent position.
	for (i = 0; i < (*T)->num_tree_nodes; i++)
	{
		if( (*T)->tree_nodes[i])
			(*T)->tree_nodes[i]->set_num_mask_words((*T)->num_mask_words);
	}

	// Fill the bag vecs
	(*T)->fill_bag_vecs();

	// Create the aux_info vector - for WIS, use this to contain the # of ind. sets
	/// found when processing each tree node (invalid for nice nonleaf nodes since all ind.
	// sets are not actually found in this case)
	info->aux_info = new int[(*T)->num_tree_nodes];
	for (i = 0; i < (*T)->num_tree_nodes; i++)
		info->aux_info[i] = 0;

}

/**
* Function to output the results of the WIS calculation to the provided stream.
*/
void print_WIS_results(FILE *stream, TDTree *T, DP_info *info)
{
	if(info->verbose || info->very_verbose)
	{
		// Print out:
		// filename G_n G_m td_type PC/NPC MN/NMN width refined_width # treenodes # leafs leaf_time
		// nonleaf_time introtime forgettime jointime total_pc_entries total_entries MEM reconsruct_time avg_pc _proportion obj

		// Print header line
		if(!info->noheader)
			fprintf(stream,"file n m td_type PC/NPC MN/NMN FC/NFC width ref_width num_tnodes num_leafs leaf_t "
			"nonleaf_t intro_t forget_t join_t tot_pc_entries tot_entries tot_ref_pc_entries tot_ref_entries mem_est "
			"recon_t avg_pc_perc obj\n");
		char td_type[3];
		if (info->gavril)
			sprintf(td_type, "GV");
		if (info->BK)
			sprintf(td_type, "BK");
		if (info->nice)
			sprintf(td_type, "NC");
		if (info->read_tree)
			sprintf(td_type, "FL");

		fprintf(stream, "%s %d %d %s ", info->DIMACS_file, T->G->get_num_nodes(),
			T->G->get_num_edges(), td_type);
		if (info->parent_child)
			fprintf(stream, "PC ");
		else
			fprintf(stream, "NPC ");

		if (info->make_nice)
			fprintf(stream, "MN ");
		else
			fprintf(stream, "NMN ");

		if (info->free_children)
			fprintf(stream, "FC ");
		else
			fprintf(stream, "NFC ");

		fprintf(stream, "%d %d %d %d ", T->width, T->refined_width,
			T->num_tree_nodes, T->num_leafs);
		fprintf(
			stream,
			"%3.2f %3.2f %3.2f %3.2f %3.2f",
			info->leaf_time,
			info->introduce_time + info->forget_time + info->join_time
			+ info->nonleaf_time, info->introduce_time,
			info->forget_time, info->join_time);
		fprintf(stream, " %lld %lld %lld %lld", info->orig_total_pc_table_entries,
			info->orig_total_table_entries, info->total_pc_table_entries,
			info->total_table_entries);
#if !(__CYGWIN__ || WIN32 || _WIN32)
		fprintf(stream, " %d ", getHWmem());
#else
		fprintf(stream," -1 ");
#endif
		fprintf(stream,"%3.1f ", info->mem_estimate);
		fprintf(stream, "%3.3f ", info->reconstruct_time);
		fprintf(stream, "%0.3f ", info->avg_pc_proportion);
		// obj is last entry
		fprintf(stream, "%d\n", info->opt_obj);
	}
	else
	{

		if(!info->noheader)
		{
			fprintf(stream,"file n m w obj");
			if(info->mem_estimate>0)
				fprintf(stream," num_pc_entries num_entries est_entries\n");
			else
				fprintf(stream,"\n");
		}
		fprintf(stream,"%s %d %d %d %d",info->DIMACS_file, T->G->get_num_nodes(),
			T->G->get_num_edges(), T->width, info->opt_obj);
		if(info->mem_estimate>0)
			fprintf(stream," %lld %lld %lld\n",info->orig_total_pc_table_entries, 
			info->orig_total_table_entries, (unsigned long long)info->mem_estimate);
		else
			fprintf(stream,"\n");
	}

}


/**
* Finds the optimal solution in the given tree by looking in the root node's table.
* Also populates two provided lists with information about how the optimal solution
* intersects with the parent's bag.
*/
int get_optimal_obj(TDTree *T, list<int> *root_intersection,
	list<int> *root_difference, DP_info *info)
{
	TDSolution *hash_ptr, *temp_sol;
	info->opt_obj = 0;

	// Fill the bag vecs
	T->fill_bag_vecs();
	int root_size = (int) T->tree_nodes[T->root_node]->bag.size();
	HASH_ITER(hh,T->tree_nodes[T->root_node]->hash_table,hash_ptr,temp_sol)
	{
		if (hash_ptr->value > info->opt_obj)
		{
			// New best solution
			root_intersection->clear();
			root_difference->clear();
			info->opt_obj = hash_ptr->value;
			for (int i = 0; i < root_size; i++)
				if (hash_ptr->mask->test_bit(i))
					root_intersection->push_back(
					T->tree_nodes[T->root_node]->bag_vec[i]);
				else
					root_difference->push_back(
					T->tree_nodes[T->root_node]->bag_vec[i]);
			root_intersection->sort();
			root_difference->sort();
		}
	}

	return info->opt_obj;

}


/**
* Estimates the memory usage in a given tree decomposition by computing the actual
* sum of the binomial coefficients using exact arithmetic and GMP. Returns the
* estimated number of independent sets to be found in all nodes in the tree.
* Writes details to provided outfile
*/
double estimate_memory_usage(TDTree *T, vector<int> *walk, const char *outfile)
{
#if !HAS_GMP
	fprintf(stderr,"HAS_GMP is set to 0\nCannot run %s!",__FUNCTION__);
	return 0;
#else

	FILE *out=fopen(outfile,"w");
	if(!out)
	{
		fatal_error("%s:  can't open %s for writing\n",__FUNCTION__,outfile);
		exit(-1);
	}

	vector<int_int> stats(T->num_tree_nodes);
	create_subgraph_stats(T,&stats);
	vector<int> inv_walk(T->num_tree_nodes);

	mpz_t sum, w_choose_k, rho_num, rho_denom, sum_num, sum_denom;
	mpq_t q_sum, q_num, q_denom, q_term, grand_total;

	mpz_init(sum); mpz_init(w_choose_k); mpz_init(rho_num); mpz_init(rho_denom);
	mpz_init(sum_num); mpz_init(sum_denom);
	mpq_init(q_sum); mpq_init(q_num); mpq_init(q_denom); mpq_init(q_term);
	mpq_init(grand_total); mpq_set_d(grand_total,0);

	int i,j;
	for(i=0;i<T->num_tree_nodes;i++)
		inv_walk[walk->at(i)]=i;
	for(j=0;j<T->num_tree_nodes;j++)
	{
		i=walk->at(j);
		fprintf(out,"%d %d %d ",i,stats[i].p1,stats[i].p2);
		// calculate the estimate of # ind. sets
		unsigned long ww=(unsigned long)stats[i].p1;
		unsigned long ss=(unsigned long)stats[i].p2;

		// rho=1-s/(w choose 2)=((w choose 2)-s)/(w choose 2)
		// rho_denom=w choose 2
		mpz_set_ui(rho_denom, ww*(ww-1)/2);
		// rho_num= rho_denom-s;
		mpz_sub_ui(rho_num,rho_denom,ss);
		mpq_set_ui(q_sum,1,1);
		unsigned long exponent;
		for(int k=1;k<=stats[i].p1;k++)
		{
			//printf("k=%d\n",k);
			//mpq_out_str(stderr,10,q_sum);
			exponent=(unsigned long)k*((unsigned long)k-1)/2;
			mpz_bin_uiui(w_choose_k,ww,(unsigned long)k);
			mpz_pow_ui(sum_num,rho_num,exponent);
			mpz_pow_ui(sum_denom,rho_denom,exponent);
			mpz_mul(sum_num,sum_num,w_choose_k);
			mpq_set_z(q_num,sum_num);
			mpq_set_z(q_denom,sum_denom);
			mpq_div(q_term,q_num,q_denom);
			mpq_canonicalize(q_term);
			mpq_add(q_sum,q_sum,q_term);
		}
		mpq_add(grand_total, grand_total, q_sum);
		double ans=mpq_get_d(q_sum);
		fprintf(out," %f\n",ans);
	}
	double final_ans=mpq_get_d(grand_total);
	double log2_final=log(final_ans)/log(2.0);
	fprintf(out,"Total %f\nLog_2(Total) %f\n",final_ans,log2_final);

	// Clean up GMP stuff
	mpz_clear(sum); mpz_clear(w_choose_k);mpz_clear(rho_num);
	mpz_clear(rho_denom);mpz_clear(sum_num);mpz_clear(sum_denom);
	mpq_clear(q_sum);mpq_clear(q_num); mpq_clear(q_denom);mpq_clear(q_term);
	mpq_clear(grand_total);

	fclose(out);

	return final_ans;
#endif
}

/**
* Returns a (potentially very large) double that is the estimate of the number of independent
* sets that must be stored in the table for tree node k.  If parent_intersection flag is true,
* then the estimate is for the # of ind. sets in the intersection of k's bag with his parent.
* Otherwise only k's bag is considered.
*/
double expected_num_ind_sets(TDTree *T, int k, bool parent_intersection)
{
#if !HAS_GMP
	fprintf(stderr,"HAS_GMP=0\n");
	return 0.0;
#else
	if(!(T->tree_nodes[k]))
		fatal_error("%s:  no memory for tree node k=%d\n",__FUNCTION__,k);

	list<int> subgraph_v;
	vector<int> subgraph_v_vec(T->tree_nodes[k]->bag.size());
	vector<int>::iterator vv,uu;

	if(parent_intersection)
	{
		// Create the intersection of i with his parent since this is what we actually
		// store in the DP tables
		int parent=T->tree_nodes[k]->adj.front();
		subgraph_v.clear();
		if(parent!=k)
		{
			uu=set_intersection(T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end(),
				T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end(),
				subgraph_v_vec.begin());
			// Copy to list
			for(vv=subgraph_v_vec.begin(); vv!=uu;++vv)
				subgraph_v.push_back(*vv);
		}
		else
			// must be root node
			subgraph_v=T->tree_nodes[k]->bag;
	}
	else
		// Don't worry about the parent - just use k's bag
		subgraph_v=T->tree_nodes[k]->bag;

	unsigned long ww=(unsigned long)subgraph_v.size();
	unsigned long ss=(unsigned long)T->G->get_num_edges_in_subgraph(&subgraph_v);

	mpz_t sum, w_choose_k, rho_num, rho_denom, sum_num, sum_denom;
	mpq_t q_sum, q_num, q_denom, q_term;

	mpz_init(sum); mpz_init(w_choose_k); mpz_init(rho_num); mpz_init(rho_denom);
	mpz_init(sum_num); mpz_init(sum_denom);
	mpq_init(q_sum); mpq_init(q_num); mpq_init(q_denom); mpq_init(q_term);

	// rho=1-s/(w choose 2)=((w choose 2)-s)/(w choose 2)
	// rho_denom=w choose 2
	mpz_set_ui(rho_denom, ww*(ww-1)/2);
	// rho_num= rho_denom-s;
	mpz_sub_ui(rho_num,rho_denom,ss);
	mpq_set_ui(q_sum,1,1);
	unsigned long exponent;
	for(unsigned long j=1;j<=ww;j++)
	{
		//printf("k=%d\n",k);
		//mpq_out_str(stderr,10,q_sum);
		exponent=(unsigned long)j*((unsigned long)j-1)/2;
		mpz_bin_uiui(w_choose_k,ww,(unsigned long)j);
		mpz_pow_ui(sum_num,rho_num,exponent);
		mpz_pow_ui(sum_denom,rho_denom,exponent);
		mpz_mul(sum_num,sum_num,w_choose_k);
		mpq_set_z(q_num,sum_num);
		mpq_set_z(q_denom,sum_denom);
		mpq_div(q_term,q_num,q_denom);
		mpq_canonicalize(q_term);
		mpq_add(q_sum,q_sum,q_term);
	}

	double ans=mpq_get_d(q_sum);

	mpz_clear(sum); mpz_clear(w_choose_k);mpz_clear(rho_num);
	mpz_clear(rho_denom);mpz_clear(sum_num);mpz_clear(sum_denom);
	mpq_clear(q_sum);mpq_clear(q_num); mpq_clear(q_denom);mpq_clear(q_term);
	return ans;
#endif
}


