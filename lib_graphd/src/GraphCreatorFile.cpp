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

#include "GraphCreatorFile.h"
#include "Graph.h"
#include "VertexWeightedGraph.h"
#include "GraphException.h"
#include <iostream>

namespace Graph {
    GraphCreatorFile::GraphCreatorFile(){
        this->factory_rw = new GraphReaderWriterFactory();
    }

    GraphCreatorFile::GraphCreatorFile(string file, string graphType){
        this->factory_rw = new GraphReaderWriterFactory();
        this->file_name = file;
        this->graph_type = graphType;
    }

    GraphCreatorFile::~GraphCreatorFile(){
        delete this->factory_rw;
    }

    string GraphCreatorFile::get_file_name() const {
        return file_name;
    }

    string GraphCreatorFile::get_graph_type() const {
        return graph_type;
    }

    void GraphCreatorFile::set_file_name(string fileName){
        this->file_name = fileName;
    }

    VertexWeightedGraph *GraphCreatorFile::create_vertex_weighted_graph(){
        VertexWeightedGraph *vwg = new VertexWeightedGraph();
        NewGraphReader gr;

        try
        {
            gr.read_graph(vwg, this->file_name, this->graph_type, true);
            vwg->set_input_file(file_name);
        }
        catch(GraphException& e)
        {
            delete vwg;
            cerr << "exception caught: " << e.what() << endl;
            const string desc("Can not create a weighted graph\n");
            throw GraphException(desc);
        }

        return vwg;
    } // create_vertex_weighted_graph

    Graph *GraphCreatorFile::create_graph(){
        Graph *g = new Graph();
        NewGraphReader gr;

        try
        {
            gr.read_graph(g, this->file_name, this->graph_type, false);
            g->set_input_file(file_name);
        }

        catch(GraphException& e)
        {
            delete g;
            cerr << "exception caught: " << e.what() << endl;
            const string desc("Can not create a mutable graph\n");
            throw GraphException(desc);
        }

        return g;
    } // create_graph

    Graph *GraphCreatorFile::create_mutable_graph(){
        Graph *g;

        try
        {
            g = this->create_graph();
        }
        catch(GraphException& e)
        {
            cerr << "exception caught: " << e.what() << endl;
            const string desc("Can not create a weighted mutable graph\n");
            throw GraphException(desc);
        }

        return g;
    } // create_mutable_graph

    VertexWeightedGraph *GraphCreatorFile::create_weighted_mutable_graph(){
        VertexWeightedGraph *vwg;

        try
        {
            vwg = this->create_vertex_weighted_graph();
        }
        catch(GraphException& e)
        {
            cerr << "exception caught: " << e.what() << endl;
            const string desc("Can not create a weighted mutable graph\n");
            throw GraphException(desc);
        }

        return vwg;
    } // create_weighted_mutable_graph

    void GraphCreatorFile::set_graph_type(string graphType){
        this->graph_type = graphType;
    }
}
