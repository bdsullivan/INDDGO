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

#include "WeightedMutableGraph.h"

namespace Graph
{

	WeightedMutableGraph::WeightedMutableGraph()
	{


	}

	WeightedMutableGraph::WeightedMutableGraph(int n) :
	Graph(n), MutableGraph(n), WeightedGraph(n)
	{
	}

	WeightedMutableGraph::~WeightedMutableGraph()
	{

	}


	WeightedMutableGraph& WeightedMutableGraph::operator=(
		const WeightedMutableGraph& rhs)
	{
		weight = rhs.weight;
		simple = rhs.simple;
		canonical = rhs.canonical;
		key = rhs.key;

		xadj = rhs.xadj;
		adjncy = rhs.adjncy;

		adj_vec = rhs.adj_vec;

		num_nodes = rhs.num_nodes;
		num_edges = rhs.num_edges;
		capacity = rhs.capacity;
		next_label = rhs.next_label;
		num_connected_components = rhs.num_connected_components;
		graph_type = rhs.graph_type;
		degree = rhs.degree;
		nodes = rhs.nodes;
		input_file = rhs.input_file;

		return *this;
	}


}
