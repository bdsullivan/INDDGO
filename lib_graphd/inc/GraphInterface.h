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

#ifndef GRAPHINTERFACE_H_
#define GRAPHINTERFACE_H_

#include <string>
#include <vector>

namespace Graph {
    class GraphInterface
    {
    virtual bool is_edge(int i, int j) const = 0;
    virtual int get_num_nodes() const = 0;
    virtual int get_num_edges() const = 0;
    virtual std::string get_graph_type() const = 0;
    virtual int get_capacity() const = 0;
    virtual int get_degree(int v) const = 0;
    virtual std::vector <int> get_degree() const = 0;
    virtual int *serialize() = 0;
    virtual void deserialize(int *buffer) = 0;
    };
}

#endif /* GRAPHINTERFACE_H_ */
