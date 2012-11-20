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

namespace Graph
{

	GraphUtil::GraphUtil()
	{

	}

	GraphUtil::~GraphUtil()
	{

	}
  
 
  
  /**
   * Recomputes the entries in the degree[] array using the size of
   * the adjacency lists. Updates the number of edges in the graph.
   */
  void GraphUtil::recompute_degrees(MutableGraph *wmg)
  {
    wmg->num_edges = 0;
		int i;
		//set the degrees to 0 - should wmg be -1 if the nodes label is -1?
		// CSG - all these loops used to start at i=1?
		fill(wmg->degree.begin(), wmg->degree.end(), 0);

		//not strictly necessary to symmetrize, but much easier - lists are sorted now.
		if (wmg->simple == true)
		{
			for (i = 0; i < wmg->capacity; i++)
			{
				wmg->degree[i] = wmg->nodes[i].nbrs.size();
				wmg->num_edges += wmg->degree[i];
			}
			// We double counted all the edges;
			wmg->num_edges = wmg->num_edges / 2;
		}
		else
		{
			print_message(0,"Not simple??\n");
			// we could have loops or multiedges - the former are the only thing that screw up our accounting.
			size_t origsize, newsize;
			int loops = 0;
			int j;
			for (i = 0; i < wmg->capacity; i++)
			{
				origsize = wmg->nodes[i].nbrs.size();
				wmg->nodes[i].nbrs.remove(i);
				newsize = wmg->nodes[i].nbrs.size();
				if (wmg->nodes[i].nbrs.size() < origsize)
				{
					loops += (int) ((((((((((((((((((origsize
						- newsize))))))))))))))))));
					wmg->num_edges +=
						(int) (((((((((((((((((newsize))))))))))))))))); //all but the loop are normal edges
					for (j = 0;
						j
						< (int) ((((((((((((((((((origsize - newsize))))))))))))))))));
					j++)
					{
						wmg->nodes[i].nbrs.push_back(i); //add it back
					}
					wmg->degree[i] = wmg->nodes[i].nbrs.size();
				}
				else
				{
					wmg->degree[i] = origsize;
					wmg->num_edges += wmg->degree[i]; //all are normal edges
				}
			}

			// We double counted all the normal edges & need to add the loops
			wmg->num_edges = wmg->num_edges / 2 + loops;
		}

		return;
	}

	/**
	* Returns the index of a vertex with the lowest degree. If there are multiple minimum degree vertices,
	* then a random representative is selected.
	*/
	int GraphUtil::get_random_high_degree_vertex(MutableGraph *mg) const
	{
		vector<int> high_degree_vs(mg->capacity, -1);
		int max_deg = -INT_MAX;
		int i, j;
		j = 0;
		for (i = 0; i < mg->capacity; i++)
		{
			if (mg->nodes[i].label != -1
				&& (int) (((((((((((((((mg->nodes[i].nbrs.size())))))))))))))))
				>= max_deg)
			{
				max_deg = mg->nodes[i].nbrs.size();
				high_degree_vs[j] = i;
				//print_message(0, "High[%d]=%d(deg=%d)\n", j, i, max_deg);
				j++;
			}
		}

		// We now have j high degree vertices in positions 0,1,...,j-1
		i = j - 1;
		while (high_degree_vs[i] == high_degree_vs[i - 1])
			i--;

		// Pick a random integer k in {i,i+1,...,j-1} and return
		// high_degree_vs[k]
		int k = rand_int(i, j - 1);
		k = high_degree_vs[k];
		// sanity check
		if (k == -1)
			fatal_error("%s: Didn't find valid low degree vertex!\n", __FUNCTION__);

		return k;
	}
	/**
	* Returns the index of a vertex with the lowest degree. If there are multiple minimum degree vertices,
	* then a random representative is selected.
	*/
	int GraphUtil::get_random_low_degree_vertex(MutableGraph *mg) const
	{
		vector<int> low_degree_vs(mg->capacity, -1);
		int min_deg = INT_MAX;
		int i, j;
		j = 0;
		for (i = 0; i < mg->capacity; i++)
		{
			if (mg->nodes[i].label != -1
				&& (int) ((((((((((((((mg->nodes[i].nbrs.size()))))))))))))))
				<= min_deg)
			{
				min_deg = mg->nodes[i].nbrs.size();
				low_degree_vs[j] = i;
				print_message(1, "Low[%d]=%d(deg=%d)\n", j, i, min_deg);
				j++;
			}
		}

		//BDS - fixed Aug 5.
		if (j == 1)
		{
			return low_degree_vs[0];
		}
		// We now have j low degree vertices in positions 0,1,...,j-1
		i = j - 1;
		while (low_degree_vs[i] == low_degree_vs[i - 1])
			i--;

		// Pick a random integer k in {i,i+1,...,j-1} and return
		// low_degree_vs[k]
		int k = rand_int(i, j - 1);
		k = low_degree_vs[k];
		// sanity check
		if (k == -1)
			fatal_error("%s: Didn't find valid low degree vertex!\n", __FUNCTION__);

		return k;
	}
	/**
	* Given the components vector of size capacity, mg uses a non-recursive
	* procedure to find the positions in the nodes[] array of those nodes that
	* belong to the same component C as v and sets component[i]=label for
	* all i in C.  The value of label should be non-negative and the components
	* vector should be initially set with all entries -1.  Returns the
	* number of nodes in mg component.
	*/
	int GraphUtil::label_component(MutableGraph *mg, int v, int label,
		vector<int> *components)
	{
		// mg assumes that components vector is initially set to all -1's (-1's)
		if (mg->nodes[v].label == -1)
			fatal_error("%s:  Tried to find component for -1 position??\n",
			__FUNCTION__);

		if (label < 0)
			fatal_error("%s:  The value of label should be non-negative\n",
			__FUNCTION__);

		print_message(1, "Labeling component of vertex %d with %d\n", v, label);
		if ((int) (((((((((((((components->size()))))))))))))) < mg->capacity)
			fatal_error("%s:  Component vector does not have enough space.\n",
			__FUNCTION__);

		//components->resize(mg->capacity,-1);
		// Just return 0 if we have already labeled v in the components array
		if (components->at(v) != -1)
			return 0;

		int j, cnt = 0;
		list<int> S;
		list<int>::iterator ii;
		// Put v on the stack
		S.push_front(v);
		while (!(S.empty()))
		{
			j = S.front();
			S.pop_front();
			if (components->at(j) != label)
			{
				components->at(j) = label;
				cnt++;
				for (ii = mg->nodes[j].nbrs.begin(); ii != mg->nodes[j].nbrs.end();
					++ii)
				{
					if (components->at(*ii) != label)
					{
						print_message(100, "Pushing %d\n", *ii);
						S.push_front(*ii);
					}
				}

			}

		}

		return cnt;
	}

	/**
	* Given the components vector of size capacity, mg uses a recursive
	* procedure to find the positions in the nodes[] array of those nodes that
	* belong to the same component C as v and sets component[i]=label for
	* all i in C.  The value of label should be non-negative and the components
	* vector should be initially set with all entries -1.  Returns the
	* number of nodes in mg component.
	*/
	int GraphUtil::rec_label_component(MutableGraph *mg, int v, int label,
		vector<int> *components)
	{
		// mg assumes that components vector is initially set to all -1's (-1's)
		if (mg->nodes[v].label == -1)
			fatal_error("%s:  Tried to find component for -1 position??\n",
			__FUNCTION__);

		if (label < 0)
			fatal_error("%s:  The value of label should be non-negative\n",
			__FUNCTION__);

		print_message(1, "Labeling component of vertex %d with %d\n", v, label);
		if ((int) ((((((((((((components->size())))))))))))) < mg->capacity)
			fatal_error("%s:  Component vector does not have enough space.\n",
			__FUNCTION__);

		//components->resize(mg->capacity,-1);
		// Just return 0 if we have already labeled v in the components array
		if (components->at(v) != -1)
			return 0;

		// Mark v
		components->at(v) = label;
		list<int>::iterator ii;
		// Recursive DFS
		for (ii = mg->nodes[v].nbrs.begin(); ii != mg->nodes[v].nbrs.end(); ++ii)
		{
			print_message(1, "Neighbor %d of %d\n", *ii, v);
			// Put if check here to eliminate unnecessary function calls?
			if (components->at(*ii) == -1)
				rec_label_component(mg, *ii, label, components);

		}
		// Count the # of nodes in mg component
		int cnt = 0;
		for (int i = 0; i < mg->capacity; i++)
			if (components->at(i) == label)
				cnt++;

		return cnt;
	}
	/**
	* Non-recursive function that labels all connected components in the graph
	* by filling the components vector with the component label for each node.
	* Every entry in the components vector is initially set to
	* -1 and is eventually populated with integers in the range
	* 0,1,...,num_components-1.  If components[i]=k, then mg means
	* that the node in position i in the nodes[] array is in component k.
	* Returns the number of components found.
	*/
	int GraphUtil::label_all_components(MutableGraph *mg,
		vector<int> *components)
	{
		// Make sure the components vector is the right size
		if ((int) (((((((((((components->size()))))))))))) < mg->capacity)
			components->resize(mg->capacity, -1);

		else
			// Initialize components vector to all -1's
			fill(components->begin(), components->end(), -1);

		int i = 0, j = 0;
		// Find the first "active" node
		while (mg->nodes[j].label == -1)
			j++;

		list<int> S;
		list<int>::iterator ii;
		// Use S as a stack - use only push_front, pop_front
		bool *visited = new bool[mg->capacity];
		int *origin = new int[mg->capacity];
		int *component_labels = new int[mg->capacity];
		for (i = 0; i < mg->capacity; i++)
		{
			visited[i] = false;
			origin[i] = component_labels[i] = -1;
		}
		int numc = 0;
		for (i = 0; i < mg->capacity; i++)
		{
			if (!(visited[i]) && mg->nodes[i].label != -1)
			{
				S.push_front(i);
				while (!(S.empty()))
				{
					j = S.front();
					S.pop_front();
					if (!visited[j])
					{
						visited[j] = true;
						origin[j] = i;
						if (origin[j] == j)
						{
							// mg indicates a new component
							// Any node u in mg component will have origin[u]=j
							// Keep track of these "special" values by marking position
							// j in the component_labels[] array and incrementing numc
							component_labels[j] = numc;
							numc++;
						}
						// Check out j's neighbors - assume these are not -1 nodes!
						for (ii = mg->nodes[j].nbrs.begin();
							ii != mg->nodes[j].nbrs.end(); ++ii)
						{
							if (!visited[*ii])
								S.push_front(*ii);

						}
					}

				}

			}

			if (mg->nodes[i].label != -1)
				// i is visited once we get here, so we know the value of components[i]
				components->at(i) = component_labels[origin[i]];

		}
		// Clean up
		delete[] visited;
		delete[] origin;
		delete[] component_labels;
		// Set the value of and return the # of connected components
		mg->num_connected_components = numc;
		return mg->num_connected_components;
	}

	/**
	* Recursive function that labels all connected components in the graph
	* by filling the components vector with the component label for each node.
	* Every entry in the components vector is initially set to
	* -1 and is eventually populated with integers in the range
	* 0,1,...,num_components-1.  If components[i]=k, then mg means
	* that the node in position i in the nodes[] array is in component k.
	* Returns the number of components found.
	*/
	int GraphUtil::rec_label_all_components(MutableGraph *mg,
		vector<int> *components)
	{
		// Make sure the components vector is the right size
		if ((int) ((((((((((components->size())))))))))) < mg->capacity)
			components->resize(mg->capacity, -1);

		else
			// Initialize components vector to all -1's
			fill(components->begin(), components->end(), -1);

		int i = 0, label = 0, j = 0;
		// Find the first "active" node
		while (mg->nodes[j].label == -1)
			j++;

		// Label the component containing node i
		rec_label_component(mg, j, label, components);
		for (i = j + 1; i < mg->capacity; i++)
		{
			if (components->at(i) == -1 && mg->nodes[i].label != -1)
			{
				// We found an active node that we haven't seen before
				label++;
				rec_label_component(mg, i, label, components);
			}
		}

		// We now know number of components
		mg->num_connected_components = label + 1;
		// components[] will now contain integers in 0,1,...,label-1
		// where components[v]=c means that node v is in component c
		// If components[v]=-1, then mg node must have an -1 label
		// and is therefore not in a component of the graph
		// Return the # of components we discovered
		return mg->num_connected_components;
	}

	/**
	* Uses a recursive function to populate the members list with the
	* positions in the nodes[] array of those
	* nodes that belong to the saem component as v.
	* Returns the number of nodes in mg component.
	*/
	int GraphUtil::rec_find_component(MutableGraph *mg, int v,
		list<int> *members)
	{
		// Create a vector of -1's and then label v's component w/ 0
		vector<int> components(mg->capacity, -1);
		rec_label_component(mg, v, 0, &components);
		members->clear();
		for (int i = 0; i < mg->capacity; i++)
		{
			if (components.at(i) == 0)
				members->push_back(i);

		}
		// Sort and return the component size
		members->sort();
		return members->size();
	}

	/**
	* Uses a non-recursive function to populate the members list with the
	* positions in the nodes[] array of those
	* nodes that belong to the saem component as v.
	* Returns the number of nodes in mg component.
	* Oct 18 2010 - was not populating members list, so corrected.
	*/
	int GraphUtil::find_component(MutableGraph *mg, int v, list<int> *members)
	{
		// mg assumes that components vector is initially set to all -1's (-1's)
		if (mg->nodes[v].label == -1)
			fatal_error("%s:  Tried to find component for -1 position??\n",
			__FUNCTION__);

		members->clear();
		int j, cnt = 0;
		vector<bool> visited(mg->capacity, false);
		list<int> S;
		list<int>::iterator ii;
		// Put v on the stack
		S.push_front(v);
		while (!(S.empty()))
		{
			j = S.front();
			S.pop_front();
			if (visited[j] == false)
			{
				visited[j] = true;
				cnt++;
				members->push_back(j);
				for (ii = mg->nodes[j].nbrs.begin(); ii != mg->nodes[j].nbrs.end();
					++ii)
				{
					if (visited[*ii] == false)
					{
						print_message(100, "Pushing %d\n", *ii);
						S.push_front(*ii);
					}
				}

			}

		}

		return cnt;
	}

	/**
	* Populate the xadj and adjncy vectors with the graph data.
	*/
	void GraphUtil::populate_CRS(MutableGraph *mg)
	{
		// Only works if capacity=num_nodes
		if (mg->capacity != mg->num_nodes)
			fatal_error("%s:  Requires capacity=num_nodes!\n", __FUNCTION__);

		mg->xadj.resize(mg->num_nodes + 1);
		mg->adjncy.resize(2 * mg->num_edges);
		print_message(10, "%s:  %d nodes, %d edges\n", __FUNCTION__, mg->num_nodes,
			mg->num_edges);
		// Nodes are numbered 0,1,...,n-1
		// The neighbors of node i are adjncy[xadj[i]],adjncy[xadj[i]+1],...adjncy[xadj[i+1]-1]
		int i, j = 0;
		for (i = 0; i < mg->num_nodes; i++)
		{
			mg->xadj[i] = j;
			for (list<int>::iterator ii = mg->nodes[i].nbrs.begin();
				ii != mg->nodes[i].nbrs.end(); ++ii)
			{
				mg->adjncy[j] = *ii;
				j++;
			}
		}

		mg->xadj[mg->num_nodes] = j;
		print_message(1, "Set xadj[%d]=%d (nedges=%d)\n", mg->num_nodes, j,
			mg->num_edges);
		return;
	}

	/**
	* Free the memory used to create the CRS format.
	*/
	void GraphUtil::free_CRS(MutableGraph *mg)
	{
		mg->xadj.clear();
		mg->adjncy.clear();
	}

	/**
	* Non-recursive function that fills the members vector with the
	* lists of nodes belonging to the components.  The members vector
	* is resized so that upon completion it has num_connected_components
	* entries and members[i] is a pointer to a list of the nodes in the
	* i-th component.
	*/
	int GraphUtil::find_all_components(MutableGraph *mg,
		vector<list<int> *> *members)
	{
		int i = 0, j = 0;
		// Find the first "active" node
		while (mg->nodes[j].label == GD_UNDEFINED)
			j++;

		//list<int> *L;
		list<int> S;
		list<int>::iterator ii;
		// Use S as a stack - use only push_front, pop_front
		bool *visited = new bool[mg->capacity];
		int *origin = new int[mg->capacity];
		int *component_labels = new int[mg->capacity];
		for (i = 0; i < mg->capacity; i++)
		{
			visited[i] = false;
			origin[i] = component_labels[i] = -1;
		}
		int numc = 0;
		for (i = 0; i < mg->capacity; i++)
		{
			if (!(visited[i]) && mg->nodes[i].label != GD_UNDEFINED)
			{
				S.push_front(i);
				while (!(S.empty()))
				{
					j = S.front();
					S.pop_front();
					if (!visited[j])
					{
						visited[j] = true;
						origin[j] = i;
						if (origin[j] == j)
						{
							// mg indicates a new component
							// Any node u in mg component will have origin[u]=j
							// Keep track of these "special" values by marking position
							// j in the component_labels[] array and incrementing numc
							component_labels[j] = numc;

							//L = new list<int>;
							//members->resize(numc + 1);
							//members->at(numc) = L;
							members->push_back(new list<int>);
							numc++;
						}
						// Check out j's neighbors - assume these are not GD_UNDEFINED nodes!
						for (ii = mg->nodes[j].nbrs.begin();
							ii != mg->nodes[j].nbrs.end(); ++ii)
						{
							if (!visited[*ii])
								S.push_front(*ii);
						}
					}
				}
			}
			if (mg->nodes[i].label != GD_UNDEFINED)
				// i is visited once we get here, so we know that i belongs in list
				// component_labels[origin[i]]
				members->at(component_labels[origin[i]])->push_back(i); // Obviously...
		}
		// Clean up
		delete[] visited;
		delete[] origin;
		delete[] component_labels;
		// Sort the lists for convenience
		vector<list<int> *>::iterator jj;
		for (jj = members->begin(); jj != members->end(); ++jj)
			(*jj)->sort();

		// Set the value of and return the # of connected components
		mg->num_connected_components = numc;
		return mg->num_connected_components;
	}
	int GraphUtil::vertex_separator(MutableGraph *mg, list<int> *V,
		vector<list<int> *> *members)
	{
		//create a copy of G
		//Changed from = *mg to copy constructor July 19 - BDS
		MutableGraph H(*mg);
		//remove the vertices in V
		list<int>::iterator i;
		for (i = V->begin(); i != V->end(); ++i)
		{
			H.remove_vertex(*i);
		}
		//find the components of H
		int components = find_all_components(&H, members);
		return components;
	}

	//runs a BFS from start using
	//only vertices with allowed[v] = true. If a vertex u  is reachable,
	//Returns an array of booleans which is true for vertices which are found in the search.
	//num_reached is the number of vertices which are reachable. mg does not include
	//disallowed vertices, and does include the start.
	bool *GraphUtil::bfs(MutableGraph *mg, int start, bool *allowed,
		int *num_reached)
	{
		int *dists = bfs_dist(mg, start, allowed, num_reached);
		int j;
		bool *visited = new bool[mg->capacity];
		for (j = 0; j < mg->capacity; j++)
		{
			if (dists[j] == GD_INFINITY
				)
				visited[j] = false;
			else
				visited[j] = true;
		}
		delete[] dists;
		return visited;
	}

	/**
	* Finds the degree 0 nodes in the graph and fills the isolated_nodes
	* list with their positions in the nodes[] array.  Returns the
	* number of isolated nodes.
	*/
	int GraphUtil::find_isolated_nodes(MutableGraph *mg, list<int> *isolated_nodes)
	{
		int i, k = 0;

		for (i = 0; i < mg->capacity; i++)
		{

			if (mg->nodes[i].label != -1 && mg->nodes[i].nbrs.size() == 0)
			{
				isolated_nodes->push_back(i);
				k++;
			}
		}
		return k;
	}

	//runs a BFS from start using only vertices with allowed[v] = true.
	//Returns an array of integers giving distance from the source (0 for source).
	//Unreachable vertices have value GD_INFINITY.
	//num_reached is the number of vertices which are reachable. mg does not include
	//disallowed vertices, and does include the start.

	int *GraphUtil::bfs_dist(MutableGraph *mg, int start, bool *allowed,
		int *num_reached)
	{

		// We need the graph to be symmetric, as we're going to walk through a neighbor list
		if (!mg->canonical)
			fatal_error("%s:  must be in canonical format\n", __FUNCTION__);
		int num_found = 0;

		int j;
		int *dists = new int[mg->capacity];
		for (j = 0; j < mg->capacity; j++)
			dists[j] = GD_INFINITY;

		// S is a stack that will contain the nodes we have visited -
		// take the first one off, and then add its unvisited neighbors to the end
		list<int> S;
		list<int>::iterator ii;

		int curr_level = 0;
		int left_in_level = 1;
		int next_level = 0;

		// Put start on the stack
		S.clear();
		S.push_front(start);
		while (!(S.empty()))
		{
			if (left_in_level == 0)
			{
				curr_level++;
				left_in_level = next_level;
				next_level = 0;
			}

			// Remove the oldest element from the stack
			j = S.front();
			S.pop_front();
			left_in_level--;
			if (dists[j] == GD_INFINITY)
			{
				// We have now visited j
				dists[j] = curr_level;
				num_found++;

				// Check j's neighbors
				for (ii = mg->nodes[j].nbrs.begin(); ii != mg->nodes[j].nbrs.end();
					++ii)
				{
					if (dists[*ii] == GD_INFINITY && allowed[*ii] == true)
					{
						// Note - mg is the only place that we refer to the allowed[] vector
						// We haven't seen *ii before and it is an "acceptable" vertex,
						// so it is a candidate to be in the path - add it to the Stack
						S.push_back(*ii); // must be push_back since we pop_front and need FIFO.
						next_level++;
					}
				}
			}
		}

		*num_reached = num_found;
		// We emptied the stack.
		return dists;
	}
}

using namespace std;
 /**
   * Populates the provided Graph structure with the necessary information to do TD, visualizations, analysis.
   * Assumes graph_file is in DIMACS format.
   * If needed, stores the largest connected components in DIMACS format at graph_file.giantcomp (and populates Graph from this).
   */
void Graph::create_largestcomponent_graph(char* graph_file, WeightedMutableGraph *&G)
  {
    GraphCreatorFile creator;
    creator.set_file_name(graph_file);
    creator.set_graph_type("DIMACS");
    G = creator.create_weighted_mutable_graph();
    GraphProperties properties;
    
    if (!G)
      throw(GraphException("Failed to read graph from specified file.\n"));
    
    if (!properties.is_connected(G))
      {
  	char* bigcompfile = (char*) malloc(100);
  	sprintf(bigcompfile, "%s.giantcomp", graph_file);	
  	G->write_largest_component("DIMACS", "temp_max_comp.dimacs");
  	normalize_DIMACS_file("temp_max_comp.dimacs", bigcompfile );
  	delete G;
  	remove("temp_max_comp.dimacs");
  	creator.set_file_name(bigcompfile);
  	G = creator.create_weighted_mutable_graph();
  	G->set_input_file(bigcompfile);
  	free(bigcompfile);
      }
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




