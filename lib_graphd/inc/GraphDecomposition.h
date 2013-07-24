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

#ifndef GRAPHDECOMPOSITION_H_
#define GRAPHDECOMPOSITION_H_

#if WIN32 || _WIN32
// Minor hack to deal with apparent inconsistency in VS versions regarding
// WIN32 and _WIN32?? WIN32 is used in metis.h...
  #define WIN32 1
  #define _WIN32 1
// Define these here if using Windows
// Otherwise, they are set in the makefile
  #define HAS_METIS   0
  #define HAS_SUITESPARSE     0
//#define HAS_ARPACK 1
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdarg.h>
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>

#ifdef HAS_BOOST
  #ifdef HAS_GTEST
    #define GTEST_HAS_TR1_TUPLE 0
  #endif
  #include <boost/math/special_functions/zeta.hpp>
#endif

#include <set>
#include <map>
#include <limits.h>

using namespace std;

#if HAS_METIS
extern "C" {
#include <metis.h>
#include <metislib.h>
#include <metisbin.h> //appropriate header for smbfactor in metis programs
}

#endif

#if HAS_SUITESPARSE
extern "C" {
#include <amd.h>
}
#endif

// Some user-defined constants, always prefaced with GD_
// to avoid conflicts w/ other code, etc.
#define    GD_SUCCESS   0
#define    GD_FAILURE   -1
#define    GD_INFINITY  (1 << 30)
#define    GD_MAX(a,b)    (((a) > (b)) ? (a) : (b))
#define    GD_MIN(a,b)    (((a) < (b)) ? (a) : (b))
#define    GD_UNDEFINED    -1
#define    GD_TRUE          1
#define    GD_FALSE         0

// Set this to 0 if you want reproducability and a
// repeatable random stream
#define GD_ADD_ENTROPY      0

// Strategies to compute elimination ordering
#define     GD_MIN_DEGREE           0
#define     GD_MCS                  1
#define     GD_MCSM                 2
#define     GD_LEXM_BFS             3
#define     GD_LEXP_BFS             4
#define     GD_MIN_FILL             5
#define     GD_PKT_SORT             6
#define     GD_BATCH_MF             7
#define     GD_METIS_MMD            8
#define     GD_METIS_NODE_ND        9
//#define     GD_METIS_EDGE_ND        10
#define     GD_BETA                 11
#define     GD_AMD                  12
#define     GD_MUL_MIN_DEGREE       13
#define     GD_MINMAX_DEGREE        14

//array of names for the elimination ordering routines defined above
static const char EO_NAMES[][30] = {
    "MinDegree", "MCS", "MCS-M", "LEX-M", "LEX-P", "MinFill", "PKT_sort", "Batch-MinFill", "MetisMMD", "MetisNodeND", "DEPRECATED" /*"MetisEdgeND"*/, "BetaOrdering", "AMDOrdering", "MultipleMinDegree", "MinMaxDegree", "Unknown"
};

// Update the # of heuristics as new ones are added
#define     GD_NUM_ELIM_HEURISTICS  15

// Lower bound strategies
#define     GD_MAX_MIN_DEGREE_LB    1
#define     GD_MCS_LB               2

// Specialized graph types
#define     GD_KTREE           1

// Include header files
#include "GraphDecomposition.h"
#include "Graph.h"
//#include "Graph.h"
#include "DIMACSGraphReader.h"
#include "DIMACSGraphWriter.h"
#include "GraphDisplay.h"
#include "GraphInterface.h"
#include "GraphUtil.h"
#include "Node.h"
#include "VertexWeightedGraph.h"
#include "GraphCreatorFile.h"
#include "GraphEOUtil.h"
#include "GraphProperties.h"
#include "GraphWriter.h"
#include "RndNumGen.h"
#include "GraphCreator.h"
#include "GraphException.h"
#include "GraphReader.h"

#include "Log.h"
#include "Util.h"
#include "RndNumGen.h"
#include "Debug.h"

#endif /* GRAPHDECOMPOSITION_H_ */
