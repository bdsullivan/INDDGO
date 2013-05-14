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

#ifndef GRAPHREADER_H_
#define GRAPHREADER_H_

#include "Node.h"
#include <vector>

namespace Graph {
    class GraphReader
    {
public:

    GraphReader();
    virtual ~GraphReader();

    virtual void read_graph(const char *filename) = 0;
    virtual std::vector<int> get_degree() = 0;
    virtual std::vector<Node> get_nodes() = 0;
    virtual int get_num_edges() const = 0;
    virtual std::vector<int> get_weights() = 0;
    virtual int get_capacity() const = 0;
    virtual void set_capacity(int capacity) = 0;
    };
}

#endif /* GRAPHREADER_H_ */
