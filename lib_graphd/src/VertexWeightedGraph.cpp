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

#include "VertexWeightedGraph.h"
#include "GraphProperties.h"
#include "GraphReaderWriterFactory.h"
#include "GraphCreatorFile.h"
#include <string.h>

namespace Graph {
    VertexWeightedGraph::VertexWeightedGraph(){
    }

    VertexWeightedGraph::VertexWeightedGraph(int n) :
        Graph(n), WeightedGraph(n){
    }

    VertexWeightedGraph::~VertexWeightedGraph(){
    }

    VertexWeightedGraph& VertexWeightedGraph::operator=(const VertexWeightedGraph& rhs){
        weight = rhs.weight;
        simple = rhs.simple;
        canonical = rhs.canonical;
        key = rhs.key;

        xadj = rhs.xadj;
        adjncy = rhs.adjncy;

        adj_vec = rhs.adj_vec;

        num_nodes = rhs.num_nodes;
        num_edges = rhs.num_edges;
        capacity = rhs.capacity;
        next_label = rhs.next_label;
        num_connected_components = rhs.num_connected_components;
        graph_type = rhs.graph_type;
        degree = rhs.degree;
        nodes = rhs.nodes;
        input_file = rhs.input_file;

        return *this;
    } // =

    /* This should really live in Graph, but all the GraphCreator code is currently for WMGs only.
     * Finds all connected components and saves one with largest number of nodes to the specified filename.
     * If you desire normalized DIMACS files, you must run that after the component file is created.
     */
    void VertexWeightedGraph::write_largest_component(std::string filetype, std::string filename){
        GraphProperties properties;
        GraphWriter *writer;
        GraphReaderWriterFactory factory;
        GraphCreatorFile creator;

        writer = factory.create_writer(filetype);
        writer->set_out_file_name(filename);

        // Put the input graph in canonical form for the tests
        properties.make_canonical(this);

        list<VertexWeightedGraph *> C;
        C = creator.create_all_components(this, true);
        print_message(10, "Found %d components\n", C.size());
        list<VertexWeightedGraph *>::iterator cc;
        list<VertexWeightedGraph *>::iterator maxcomp;
        int max_size = -1;

        for(cc = C.begin(); cc != C.end(); ++cc){
            if(( (*cc)->get_num_nodes() > 1) && ( (*cc)->get_num_edges() > 1) ){
                if((*cc)->get_num_nodes() > max_size){
                    /*mark this as potential largest component*/
                    maxcomp = cc;
                    max_size = (*cc)->get_num_nodes();
                }
            }
        }

        writer->write_graph(*maxcomp);

        /*cleanup*/
        int s = C.size();
        for(int i = 0; i < s; i++){
            delete C.front();
            C.pop_front();
        }
        delete writer;
    } // write_largest_component
}

