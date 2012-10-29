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

#include "GraphReaderWriterFactory.h"
#include "DIMACSGraphReader.h"
#include "AdjMatrixGraphReader.h"
#include "AdjMatrixGraphWriter.h"
#include "DIMACSGraphWriter.h"
#include "MetisGraphWriter.h"
#include "MetisGraphReader.h"
#include "GraphVizGraphWriter.h"
#include "GraphException.h"

namespace Graph
{

	GraphReaderWriterFactory::GraphReaderWriterFactory()
	{

	}

	GraphReader *GraphReaderWriterFactory::create_reader(std::string name)
	{
		GraphReader *r;
		if (name == "DIMACS" || name == "dimacs")
		{
			r = new DIMACSGraphReader();
		}
		else if (name == "AdjMatrix" || name == "ADJMATRIX" || name == "adjmatrix")
		{
			r = new AdjMatrixGraphReader();
		}
		else if (name == "METIS" || name == "Metis" || name == "metis")
		{
			r = new MetisGraphReader();
		}
		else
		{
			std::string desc = "Reader for file type " + name
				+ " has not been implemented";
			throw GraphException(desc);
		}

		return r;
	}

	GraphWriter *GraphReaderWriterFactory::create_writer(string name,
		string outfile)
	{
		GraphWriter *w;
		if (name == "DIMACS" || name == "dimacs")
		{
			w = new DIMACSGraphWriter(outfile);
		}
		else if (name == "METIS" || name == "Metis" || name == "metis")
		{
			w = new MetisGraphWriter();
		}
		else if (name == "AdjMatrix" || name == "ADJMATRIX" || name == "adjmatrix")
		{
			w = new AdjMatrixGraphWriter(outfile);
		}
		else if (name == "GraphViz" || name == "GRAPHVIZ" || name == "graphviz")
		{
			w = new GraphVizGraphWriter(outfile);
		}
		else
		{
			std::string desc = "Reader for file type " + name
				+ " has not been implemented";
			throw GraphException(desc);
		}

		return w;
	}

	GraphReaderWriterFactory::~GraphReaderWriterFactory()
	{

	}

}
