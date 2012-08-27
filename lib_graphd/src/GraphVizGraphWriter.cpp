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

#include "GraphVizGraphWriter.h"
#include "Debug.h"
#include <stdio.h>


namespace Graph
{

	GraphVizGraphWriter::GraphVizGraphWriter()
	{
	}

	GraphVizGraphWriter::GraphVizGraphWriter(string & outfile) :
	GraphWriter(outfile)
	{
	}

	GraphVizGraphWriter::~GraphVizGraphWriter()
	{
	}

	/**
	* Writes the graph to the provide graph_viz file.
	*/
	void GraphVizGraphWriter::write_graph(Graph *g)
	{
		int i, j;
		FILE *out;
		string graphviz_file = out_file_name;
		if ((out = fopen(graphviz_file.c_str(), "w")) == NULL
			)
			fatal_error(
			"%s:  Error opening file %s for writing graphviz output\n",
			__FUNCTION__, graphviz_file.c_str());

		fprintf(out, "Graph G{\n");
		fprintf(out, "overlap=false;\n");
		list<int>::iterator it;
		vector<Node> nodes = g->get_nodes();
		list<int> nbrs;

		for (i = 0; i < nodes.size(); i++)
		{
			nbrs = nodes[i].get_nbrs();
			it = nbrs.begin();
			while (it != nbrs.end())
			{
				j = *it;
				if (i < j)
					// Make this 1-based in true DIMACS spirit
					fprintf(out, "%d -- %d;\n", i + 1, j + 1);
				++it;
			}
		}
		fprintf(out, "}\n");
		fflush(out);
		fclose(out);

		return;
	}

}
