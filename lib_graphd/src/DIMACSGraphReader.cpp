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

#include "DIMACSGraphReader.h"
#include "Log.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "Debug.h"

namespace Graph
{

	DIMACSGraphReader::DIMACSGraphReader()
	{
		this->num_edges = 0;

	}

	void DIMACSGraphReader::read_graph(const char *filename)
	{
		char line[100], format[100], x;
		int i, j, m, n, retval, id, count;
		int val;
		FILE *in;

		m = n = 0;
		count = 0;
		// Use the above to count # of connected nodes, etc. and make sure that it
		// matches what is in the existing Graph struct.

		if ((in = fopen(filename, "r")) == NULL)
		{
			FERROR("%s:  Error opening %s for reading\n", __FUNCTION__, filename);
			fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
				filename);
		}

		while (!feof(in))
		{

			retval = fscanf(in, "%2c", line);

			if (feof(in))
				break;

			if (retval != 1)
				FERROR("%s:  fscanf read %d char's (expected to read 1!)\n",
				__FUNCTION__, retval);

			switch (line[0])
			{
			case 'p':
				// This is the "problem line" - get n and m here
				// Make sure we don't already know n and m!
				if (n != 0 || m != 0)
				{
					FERROR(
						"%s:  DIMACS read error - are there more than one problem lines?\n",
						__FUNCTION__);
				}
				retval = fscanf(in, "%s %d %d", format, &n, &m);

				DEBUG("DIMACS: read n = %d, m = %d\n", n, m);

				// Simple error checking
				if (n <= 0 || m <= 0 || retval != 3)
				{
					FERROR(
						"%s:  DIMACS read error - didn't understand problem line! p %s %d %d\n",
						__FUNCTION__, format, n, m);
					exit(-1);
				}

				if (strncmp(format, "e", 1) != 0 && strncmp(format, "E", 1) != 0
					&& strncmp(format, "edge", 4) != 0
					&& strncmp(format, "EDGE", 4) != 0)
				{
					FERROR(
						"%s:  DIMACS read error - problem line - FORMAT must be one of\n"
						" e, E, edge, or EDGE\n", __FUNCTION__);
				}

				this->weights.resize(n, 1);
				this->degree.resize(n, 0);
				this->nodes.resize(n);
				this->capacity = n;
				// Make sure that we match the number of total nodes
				//			if (n != this->capacity)
				//				FERROR("%s:  Graph has %d total nodes, file has %d\n",
				//						__FUNCTION__, this->capacity, n);
				break;

			case 'c':
				// Comment line - skip and move on
				break;

			case 'n':
				// Node line - of the form n id VALUE
				retval = fscanf(in, "%d %d", &id, &val);

				// Simple error checking - make sure we know n and m already!
				if (n == 0 || m == 0)
					FERROR("%s:  DIMACS read error - node line found before problem line!\n",
					__FUNCTION__);
				if (retval != 2)
					FERROR( "%s:  DIMACS read error - didn't understand node line\n",
					__FUNCTION__);

				// Check id
				if (id < 1 || id > n)
					FERROR("%s:  DIMACS read error - node id (%d) must be between 1 and n= %d\n",
					__FUNCTION__, id, n);

				// Store value - first allocate the value array if NULL
				//if((int)this->weight.size()!=n)
				//    this->weight.resize(n);

				// Store value
				this->weights[id - 1] = val; // 1->0
				DEBUG("Storing weight[%d]= %d\n", id - 1, this->weights[id - 1]);

				break;

			case 'e':
				// Edge line - of the form e start end
				retval = fscanf(in, "%d %d", &i, &j);

				//DEBUG("Read edge %d-%d\n", i, j);

				// Simple error checking - make sure we know n and m already!
				if (n == 0 || m == 0)
					FERROR(
					"%s:  DIMACS read error - edge line found before problem line!\n",
					__FUNCTION__);

				if (retval != 2)
					FERROR( "%s:  DIMACS read error - didn't understand edge line\n",
					__FUNCTION__);
				// Check start and end values
				if (i < 1 || i > n || j < 1 || j > n)
					FERROR( "%s:  DIMACS read error - edge (u,v) (%d-%d) must have"
					" both u and v between 1 and n=%d, inclusive\n",
					__FUNCTION__, i, j, n);

				// So when we read in an edge (a,b) in the DIMACS file
				// we need to add G.label_index[b] to nodes[G.label_index[a]]
				// When going back at the end, we need to use the label field
				// The edge appears valid - store in adj lists

				this->degree[i - 1]++; // 1->0
				this->degree[j - 1]++; // 1->0

				this->nodes[i - 1].set_label(i);
				this->nodes[j - 1].set_label(j);
				// Add to i's adjlist

				this->nodes[i - 1].add_nbr(j - 1);
				this->nodes[j - 1].add_nbr(i - 1);

				// Increment edge counter
				this->num_edges++;
				//DEBUG("%d:%d\t%d\n", i, j, this->numEdges);

				break;
			case 'x':
				// Shows up in TSP files - we don't need it
				break;

			default:
				// We encountered some character that we don't understand
				FERROR("%s:  DIMACS read error - didn't understand the character %c to start a line\n",
					__FUNCTION__, line[0]);
				break;
			}; // end switch

			// Advance to next line if there is format at the end of it
			while (!feof(in) && (x = getc(in)) != '\n')
				;
		}
		fclose(in);
		DEBUG("Finish reading DIMACS file\n");

		for (i = 0; i < this->capacity; i++)
		{
			if (nodes[i].get_label() == -1)
			{
				nodes[i].set_label(i + 1);
			}
		}
		//#if 1
		//	// Sort the neighbor lists
		//	for (i = 0; i < this->capacity; i++)
		//		if (this->nodes[i].label != GD_UNDEFINED)
		//			this->nodes[i].nbrs.sort();
		//#endif
		return;
	}

	std::vector<int> DIMACSGraphReader::get_degree()
	{
		return degree;
	}

	std::vector<Node> DIMACSGraphReader::get_nodes()
	{
		return nodes;
	}

	int DIMACSGraphReader::get_num_edges() const
	{
		return num_edges;
	}

	int DIMACSGraphReader::get_capacity() const
	{
		return capacity;
	}

	void DIMACSGraphReader::set_capacity(int capacity)
	{
		this->capacity = capacity;
	}

	std::vector<int> DIMACSGraphReader::get_weights()
	{
		return weights;
	}

	DIMACSGraphReader::~DIMACSGraphReader()
	{

	}

}
