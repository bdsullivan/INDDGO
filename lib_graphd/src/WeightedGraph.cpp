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

#include "WeightedGraph.h"

namespace Graph
{

	WeightedGraph::WeightedGraph()
	{

	}

	vector<int> WeightedGraph::get_weight() const
	{
		return weight;
	}

	int WeightedGraph::get_vertex_weight(int k)
	{
		return this->weight[k];
	}
	WeightedGraph::WeightedGraph(int n) :
	Graph(n)
	{
		weight.resize(n, 0);
	}

	void WeightedGraph::set_weight(vector<int> weight)
	{
		this->weight = weight;
	}

	WeightedGraph::~WeightedGraph()
	{

	}

}
