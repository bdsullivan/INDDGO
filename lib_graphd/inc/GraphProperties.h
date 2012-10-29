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

#ifndef GRAPHPROPERTIES_H_
#define GRAPHPROPERTIES_H_
#include "GraphDecomposition.h"

using namespace std;

namespace Graph
{

	class GraphProperties
	{

	public:
		GraphProperties();
		virtual ~GraphProperties();
		void make_simple(MutableGraph *mg); //Removes all loops and duplicate edges. Sets the simple flag to true.
		void make_canonical(MutableGraph *mg);
		int make_clique(MutableGraph *mg, list<int> *vertices);
		int fill_adj_vec(MutableGraph *mg, int v);
		bool check_simple(MutableGraph *mg);
		// Connectivity functions

		bool is_connected(MutableGraph *mg);

		bool is_clique(MutableGraph *mg, list<int> *vertices);
		bool is_independent_set(MutableGraph *mg, list<int> *vertices);
		bool is_independent_set(WeightedMutableGraph *mg, list<int> *vertices,
			int *val);
		bool is_path(MutableGraph *mg, int start, int end, bool *t); // Determines existence of a path between start
		// and end involving only vertices v with t[v]=true
		bool is_path(MutableGraph *mg, int start, int end);

	};

}

#endif /* GRAPHPROPERTIES_H_ */
