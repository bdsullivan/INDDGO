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

#include "EdgelistGraphReader.h"
#include "Log.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "Debug.h"

namespace Graph
{

	EdgelistGraphReader::EdgelistGraphReader()
	{
		this->num_edges = 0;

	}

	void EdgelistGraphReader::read_graph(const char *filename)
	{
		char line[100], line2[100], format[100], x;
		int i, j, m, n, retval, id, count;
        char *retp;
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

        // get the first non-comment lines
        retp = fgets(line, 100, in);
        fprintf(stderr, "retval=%d\n", retval);
        while(line[0] == '#'){
            retp = fgets(line, 100, in);
            fprintf(stderr, "retval=%d\n", retval);
        }
        // it has out number of nodes and number of edges
        retval = sscanf(line, "%d %d", &n, &m);
        fprintf(stderr, "n=%d m=%d\n", n, m);
        this->weights.resize(n, 1);
        this->degree.resize(n, 0);
        this->nodes.resize(n);
        this->capacity = n;


		while(!feof(in)) {

			retp = fgets(line, 100, in);
            //fprintf(stderr, "Got line: %s\n",line);

            if (feof(in))
                break;

            if(line[0] == '#'){
                continue;
            }

            // Edge line - start end
            retval = sscanf(line, "%d %d", &i, &j);

				//DEBUG("Read edge %d-%d\n", i, j);

            if (retval != 2)
            {
                FERROR( "%s:  Edgelist read error - didn't understand edge line\n",
                __FUNCTION__);
            }

				// Check start and end values
            if (i < 0 || i > n-1 || j < 0 || j > n-1)
            {
                FERROR( "%s:  Edgelist read error - edge (u,v) (%d-%d) must have"
                " both u and v between 1 and n=%d, inclusive\n",
                __FUNCTION__, i, j, n);
            }
            //fprintf(stderr, "i=%d j=%d\n", i, j);

            // So when we read in an edge (a,b) in the Edgelist file
            // we need to add G.label_index[b] to nodes[G.label_index[a]]
            // When going back at the end, we need to use the label field
            // The edge appears valid - store in adj lists

            this->nodes[i].set_label(i+1);
            this->nodes[j].set_label(j+1);
            // Add to i's adjlist

            if(i != j){
                this->nodes[j].add_nbr(i);
                this->degree[j]++; // 1->0
            }
            this->nodes[i].add_nbr(j);
            this->degree[i]++; // 1->0


            // Increment edge counter
            fprintf(stderr, "Read line: %s",line);
            this->num_edges++;
				//DEBUG("%d:%d\t%d\n", i, j, this->numEdges);
        };

		DEBUG("Finish reading Edgelist file\n");
        fprintf(stderr, "EdgelistGraphReader: read in %d edges\n", this->num_edges);

		for (i = 0; i < this->capacity; i++)
		{
			if (nodes[i].get_label() == -1)
			{
				nodes[i].set_label(i+1);
			}
		}


		return;
	}

	std::vector<int> EdgelistGraphReader::get_degree()
	{
		return degree;
	}

	std::vector<Node> EdgelistGraphReader::get_nodes()
	{
		return nodes;
	}

	int EdgelistGraphReader::get_num_edges() const
	{
		return num_edges;
	}

	int EdgelistGraphReader::get_capacity() const
	{
		return capacity;
	}

	void EdgelistGraphReader::set_capacity(int capacity)
	{
		this->capacity = capacity;
	}

	std::vector<int> EdgelistGraphReader::get_weights()
	{
		return weights;
	}

	EdgelistGraphReader::~EdgelistGraphReader()
	{

	}

}
