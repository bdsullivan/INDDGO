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

#include "GraphWriter.h"

namespace Graph
{

	GraphWriter::GraphWriter()
	{

	}

	GraphWriter::GraphWriter(string &filename)
	{
		out_file_name = filename;
	}

    void GraphWriter::shuffle(Graph *g, int seed)
    {
		int capacity = g->get_capacity();
        perm.clear();
        for (int i = 0; i < capacity; i++)
            perm.push_back(i);

        int s = seed % 97;
        for (int i = 0; i < s; i++)
            random_shuffle(perm.begin(), perm.end());
            
    }

	string GraphWriter::get_out_file_name() const
	{
		return out_file_name;
	}

	void GraphWriter::set_out_file_name(string out_file_name)
	{
		this->out_file_name = out_file_name;
	}

	GraphWriter::~GraphWriter()
	{

	}

}
