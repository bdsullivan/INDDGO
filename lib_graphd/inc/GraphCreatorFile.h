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

#ifndef GRAPHCREATORFILE_H_
#define GRAPHCREATORFILE_H_

#ifdef HAS_BOOST
  #ifdef HAS_GTEST
    #define GTEST_HAS_TR1_TUPLE 0
  #endif
#endif

#include "GraphCreator.h"
#include "GraphReader.h"
#include <string>
#include <vector>

namespace Graph {
    class GraphCreatorFile : public GraphCreator
    {
private:
    std::string file_name;
    std::string graph_type;

public:
    GraphCreatorFile();
    GraphCreatorFile(std::string fileName, std::string graphType = "DIMACS");

    virtual ~GraphCreatorFile();
    virtual Graph *create_graph();
    virtual VertexWeightedGraph *create_vertex_weighted_graph();
    virtual Graph *create_mutable_graph();
    /**
     * Compatibility function, simply calls create_vertex_weighted_graph under the hood.  Deprecated.
     * \return G new graph
     */
    virtual VertexWeightedGraph *create_weighted_mutable_graph();

    std::string get_file_name() const;
    std::string get_graph_type() const;
    void set_file_name(std::string fname);
    void set_graph_type(std::string gtyep);
    };
}

#endif /* GRAPHCREATORFILE_H_ */
