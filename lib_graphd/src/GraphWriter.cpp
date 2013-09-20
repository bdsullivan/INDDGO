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
#include "Util.h"

using namespace std;

namespace Graph {
    GraphWriter::GraphWriter(){
        this->shuffle = false;
        this->shuffle_seed = 42;
    }

    GraphWriter::~GraphWriter(){
    }

    /**
     * \param[in] shuf bool, shuffle graphs out from this GraphWriter
     */
    void GraphWriter::set_shuffle(bool shuf){
        this->shuffle = shuf;
    }

    /**
     * \param[in] shuf bool, shuffle graphs out from this GraphWriter
     */
    void GraphWriter::set_shuffle_seed(int seed){
        this->shuffle_seed = seed;
    }

    /**
     * \return 0 on success, nonzero on failure
     * \param[in] g Pointer to the graph to be written
     * \param[in] filename full path of the file to be written to
     * \param[in] type format of the file to write
     * \param[in] write_vertex_weights whether to write vertex weights to the output file
     */
    int GraphWriter::write_graph(Graph *g, const string filename, const string type, bool write_vertex_weights, bool shuffle){
        string t = str_to_up(type);
        if("ADJMATRIX" == t){
            write_adjmatrix(g, filename);
        }
        else if("DIMACS" == t){
            write_dimacs(g, filename, write_vertex_weights);
        }
        else if("METIS" == t){
            write_metis(g, filename);
        }
        else if("GRAPHVIZ" == t){
            write_graphviz(g, filename);
        }
        else if("ADJLIST" == t){
            write_adjlist(g, filename);
        }
        else {
            cerr << "unknown type: " << type << endl;
        }
        return 0;
    } // write_graph

    /**
     * \return 0 on success, nonzero on failure
     * \param[in] g Pointer to the graph to be written
     * \param[in] filename full path of the file to be written to
     */
    int GraphWriter::write_adjmatrix(Graph *g, const string filename){
        FILE *fout = fopen(filename.c_str(), "w");
        if(!fout){
            fatal_error("can not open %s for write", filename.c_str());
        }

        int n = g->get_num_nodes();
        int i = 0;
        int j = 0;
        int val = 0;

        fprintf(fout, "%d \n", n);

        for(i = 0; i < n; i++){
            for(j = 0; j <  n; j++){
                val = 0;
                if(g->is_edge(i, j)){
                    val = 1;
                }
                fprintf(fout, "%d ", val);
            }
            fprintf(fout, "\n");
        }

        fclose(fout);

        return 0;
    } // write_adjmatrix

    /**
     * \return 0 on success, nonzero on failure
     * \param[in] g Pointer to the graph to be written
     * \param[in] filename full path of the file to be written to
     */
    int GraphWriter::write_adjlist(Graph *g, const string filename){
        FILE *fout = fopen(filename.c_str(), "w");
        if(!fout){
            fatal_error("can not open %s for write", filename.c_str());
        }

        int n = g->get_num_nodes();
        int i = 0;
        int j = 0;
        int val = 0;
        string line;
        string is;
        ostringstream os;

        for(i = 0; i < n; i++){
            const list<int> &nbrs = g->get_node(i)->get_nbrs_ref();
            os.str("");
            os << i;
            line = is;
            for(list<int>::const_iterator it = nbrs.begin(); it != nbrs.end(); ++it){
                os << " " << *it;
            }
            os << "\n";
            fprintf(fout, "%s", os.str().c_str());
        }

        fclose(fout);

        return 0;
    } // write_adjmatrix

    /**
     * Expects to output to a 1-based METIS file
     * \param[in] g Pointer to the graph to be written
     * \param[in] filename full path of the file to be written to
     * \param[in] write_vertex_weights whether to write vertex weights to the output file
     * \return 0 on success, nonzero on failure
     */
    int GraphWriter::write_metis(Graph *g, const string filename){
        int capacity = g->get_capacity();
        int num_nodes = g->get_num_nodes();
        int num_edges = g->get_num_edges();
        // FIXME: there should be a better way than copying our entire node vector here
        vector<Node> nodes = g->get_nodes();

        FILE *out;
        if((out = fopen(filename.c_str(), "w")) == NULL){
            fatal_error("%s:  Couldn't open %s for writing\n", __FUNCTION__,
                        filename.c_str());
        }

        if(capacity != num_nodes){
            fprintf(stderr, "%s:  Warning! Metis file may have holes.\n",
                    __FUNCTION__);
            // Print a comment warning in the file
            fprintf(
                out,
                "c Warning: capacity=%d, num_nodes=%d. Resulting graph may have holes\n",
                capacity, num_nodes);
        }

        fprintf(out, "%d %d\n", num_nodes, num_edges);
        list<int>::iterator ii;
        list<int> *nbrs;
        for(int i = 0; i < num_nodes; i++){
            nbrs = nodes[i].get_nbrs_ptr();
            for(ii = nbrs->begin(); ii != nbrs->end(); ++ii){
                fprintf(out, "%d ", (*ii) + 1);
            }
            fprintf(out, "\n");
        }

        fflush(out);
        fclose(out);
        return 0;
    } // write_metis

    /**
     * \return 0 on success, nonzero on failure
     * \param[in] g Pointer to the graph to be written
     * \param[in] filename full path of the file to be written to
     */
    int GraphWriter::write_dimacs(Graph *g, const string filename, bool write_vertex_weights){
        //FIXME: figure out if this is actually necessary
        vector<int> perm;
        int capacity = g->get_capacity();
        int num_nodes = g->get_num_nodes();
        int num_edges = g->get_num_edges();
        //vector<Node> nodes = g->get_nodes();

        for(int i = 0; i < capacity; i++){
            perm.push_back(i);
        }

        if(this->shuffle){
            int s = this->shuffle_seed % 97;
            for(int i = 0; i < s; i++){
                random_shuffle(perm.begin(), perm.end());
            }
        }

        FILE *out;
        if((out = fopen(filename.c_str(), "w")) == NULL){
            fatal_error("%s:  Couldn't open %s for writing\n", __FUNCTION__,
                        filename.c_str());
        }

        // Print out a comment line
        //fprintf(out,"c Graph has %d connected components\n",g->num_connected_components);

        // Will write out non-canonical - note that we may have "holes" in the numbering
        // Print a warning to stderr and the file in this case
        if(capacity != num_nodes){
            fprintf(stderr, "%s:  Warning! DIMACS file may have holes.\n",
                    __FUNCTION__);
            // Print a comment warning in the file
            fprintf(
                out,
                "c Warning: capacity=%d, num_nodes=%d. Resulting graph may have holes\n",
                capacity, num_nodes);
        }

        //EDITED TO PRINT p edge instead of p e for compatibility with Koster.
        fprintf(out, "p edge %d %d\n", num_nodes, num_edges);
        // Write edges out using the label field from the nodes
        // Only write out edges u-v with u<=v - so file is at least simple
        list<int>::iterator ii;
        int edges_written = 0;
        list<int> nbrs;
        for(int i = 0; i < num_nodes; i++){
            // CSG- could add sort here although that won't work for the randomized labels
            // g->nodes[perm[i]].nbrs.sort();
            nbrs = g->get_node(i)->get_nbrs();
            for(ii = nbrs.begin(); ii != nbrs.end(); ++ii){
                //if(i<=*ii)
                // CSG - we want the edges in the file to be u-v with u<=v
                if(g->get_node(perm[i])->get_label() <= g->get_node(perm[*ii])->get_label()){
                    // Changing june 21 for testing - not using labels
                    //if(i<=*ii)
                    edges_written++;
                    fprintf(out, "e %d %d\n", g->get_node(perm[i])->get_label(),
                            g->get_node(perm[*ii])->get_label());
                    //,i+1,*ii+1);//
                }
            }
        }
        if(edges_written != num_edges){
            fatal_error("%s:  Didn't write the correct number of edges! (%d!=%d)\n",
                        __FUNCTION__, edges_written, num_edges);
        }

        fflush(out);
        fclose(out);

        return 0;
    } // write_dimacs

    int GraphWriter::write_graphviz(Graph *g, const string filename){
        int i, j;
        FILE *out;
        if((out = fopen(filename.c_str(), "w")) == NULL){
            fatal_error( "%s:  Error opening file %s for writing graphviz output\n", __FUNCTION__, filename.c_str());
        }

        fprintf(out, "Graph G{\n");
        fprintf(out, "overlap=false;\n");
        list<int>::iterator it;
        //FIXME: really shouldn't be copying our vector around
        vector<Node> nodes = g->get_nodes();
        list<int> nbrs;

        for(i = 0; i < g->get_num_nodes(); i++){
            nbrs = nodes[i].get_nbrs();
            it = nbrs.begin();
            while(it != nbrs.end()){
                j = *it;
                if(i < j){
                    // Make this 1-based in true DIMACS spirit
                    fprintf(out, "%d -- %d;\n", i + 1, j + 1);
                }
                ++it;
            }
        }
        fprintf(out, "}\n");
        fflush(out);
        fclose(out);

        return(0);
    } // write_graphviz
}
