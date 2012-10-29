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

#ifndef GRAPH_H_
#define GRAPH_H_

#include "GraphInterface.h"
#include "Node.h"
#include <string>

using namespace std;

namespace Graph
{

	class Graph: public GraphInterface
	{

	protected:
		int num_nodes;
		int num_edges;
		int capacity;
		int next_label;
		int num_connected_components;
		string graph_type;
		vector<int> degree;
		vector<Node> nodes;
		string input_file;

	public:
		Graph();
		Graph(int n);
		virtual ~Graph();

		void set_degree(vector<int> degree);
		void set_nodes(vector<Node> nodes);
		void set_graph_type(string graphType);
		void set_num_edges(int numEdges);
		void set_num_nodes(int numNodes);
		void set_input_file(string input_file);

		vector<int> get_degree() const;
		vector<Node> get_nodes() const;
		Node *get_node(int i);
		string get_graph_type() const;
		int get_num_edges() const;
		int get_num_nodes() const;
		virtual int get_degree(int v) const;
		virtual bool is_edge(int i, int j) const;
		virtual int *serialize();
		virtual void deserialize(int *buffer);
		virtual int get_capacity() const;
		void set_capacity(int capacity);
		int get_next_label() const;
		void set_next_label(int nextLabel);
		int get_num_components() const;
		void set_num_components(int num_components);
		string get_input_file();
		int get_num_edges_in_subgraph(list<int> *vertices);
		void complement();

		friend class GraphUtil;
		friend class GraphProperties;
		friend class GraphDisplay;
		friend class GraphEOUtil;
		friend class GraphCreator;

	};

}

#endif /* GRAPH_H_ */



