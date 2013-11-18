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

#ifndef _TREE_DEC_H
#define _TREE_DEC_H

// Node types in a nice or not-nice decomposition
#define TD_UNKNOWN_NODE         0
#define TD_NICE_LEAF_NODE       1
#define TD_NONNICE_LEAF_NODE    2
#define TD_INTRODUCE_NODE       3
#define TD_FORGET_NODE          4
#define TD_JOIN_NODE            5
#define TD_NONLEAF_NODE         6

// Algorithms for estimating cost.
#define TD_WIDTH_ONLY       1
#define TD_BAG_SUM_ALG      2
#define TD_WIS_SUM_ALG      3

// Algorithms for choosing root. 
#define TD_ROOT_ASPVALL     1
#define TD_ROOT_ALLCHILDREN 2

// Algorithms for generation
#define TD_SUPERETREE 1
#define TD_GAVRIL     2
#define TD_BK         3
#define TD_NICE       4
#define TD_BFS        5

// Types of graphviz output that are available for tree decompositions
#define GV_BAG_LABELS  0
#define GV_SCALE_BAGS  1
#define GV_TREE_ONLY   2
#define GV_COLORS      3

// Visualization related constants
#define TD_BAG_SCALE   1.5  //Multiplier used for bag radius in visualization. 

// Include the header files specific to tree decompositions
// Include this first - custom implementation of some big integer math for 
// bit masking, etc.

#include "uthash.h"
#include "bigint.h"
#include <time.h>
#include "GraphDecomposition.h"
#include "TDSolution.h"
#include "TDTreeNode.h"
#include "TDTree.h"
#include "TDDynamicProgramming.h"
#include "TDMadTreeNode.h"

#if HAS_SUITESPARSE
extern "C"{
#include "SuperWrap.h"
}
#endif

void form_td(int td_alg, TDTree**T, vector<int>* ordering);
void td_size_histogram(TDTree *T, FILE *stream);
void create_tree_decomposition(Graph::VertexWeightedGraph *G, TDTree **T, bool read_tree, char *tree_infile, bool read_ordering, bool scotch, char *ord_file, int elim_order_type, int start_v, int td_alg, bool make_nice, bool timings);
void bag_statistics(TDTree *T, const vector<double> &scores, vector<double> *stats, int stat_flag);
void find_eccentricities(TDTree &T, vector<int>* ec);
void bag_lengths(TDTree *T, vector<int> *lengths);
#endif

