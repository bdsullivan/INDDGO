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

#ifndef GRAPHWRITER_H_
#define GRAPHWRITER_H_
#include "Graph.h"
#include "WeightedGraph.h"
#include "MutableGraph.h"
#include "WeightedMutableGraph.h"
#include <string>
#include <algorithm>

using namespace std;

namespace Graph {
    class GraphWriter {
protected:
    string out_file_name;
    vector<int> perm;

public:
    GraphWriter();
    GraphWriter(string &filename);
    virtual ~GraphWriter();
    virtual void write_graph(Graph *g) = 0;
    virtual void shuffle(Graph *g, int seed = 1);
    string get_out_file_name() const;
    void set_out_file_name(string out_file_name);
    };
}

#endif /* GRAPHWRITER_H_ */
