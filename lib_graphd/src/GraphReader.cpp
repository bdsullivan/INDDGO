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

#include "GraphReader.h"
#include "Log.h"
#include "Debug.h"
#include "Util.h"
#if !WIN32
  #include <strings.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

namespace Graph {
    GraphReader::GraphReader(){
    }

    GraphReader::~GraphReader(){
    }

/**
 * Public function to read in a graph
 * \param[in,out] g a graph object
 * \param[in] filename full path to the file to read in (relative or absolute)
 * \param[in] type type of file to be read in
 * \param[in] read_vertex_weights read vertex weights from input file
 * \return code 0 on success, nonzero on failure
 */
    int GraphReader::read_graph(Graph *g, const string filename, string type, bool read_vertex_weights){
        string t = str_to_up(type);
        int retcode = -1;
        if("EDGE" == t){
            retcode = GraphReader::read_edgelist(g, filename);
            if(read_vertex_weights){
                VertexWeightedGraph *vg;
                vg = (VertexWeightedGraph *) g;
                vg->weight.resize(vg->num_nodes, 1);
            }
        }
        else if("ADJMATRIX" == t){
            retcode = GraphReader::read_adjmatrix(g, filename);
            if(read_vertex_weights){
                VertexWeightedGraph *vg;
                vg = (VertexWeightedGraph *) g;
                vg->weight.resize(vg->num_nodes, 1);
            }
        }
        else if("ADJLIST" == t){
            retcode = GraphReader::read_adjlist(g, filename);
            if(read_vertex_weights){
                VertexWeightedGraph *vg;
                vg = (VertexWeightedGraph *) g;
                vg->weight.resize(vg->num_nodes, 1);
            }
        }
        else if("METIS" == t){
            retcode = GraphReader::read_metis(g, filename);
            if(read_vertex_weights){
                VertexWeightedGraph *vg;
                vg = (VertexWeightedGraph *) g;
                vg->weight.resize(vg->num_nodes, 1);
            }
        }
        else if("DIMACS" == t){
            retcode = GraphReader::read_dimacs(g, filename, read_vertex_weights);
        }
        else {
            cerr << "Error: type " << type << " not currently supported\n";
            return -1;
        }

        g->set_input_file(filename);
        g->set_graph_type(type);
        g->resize_adj_vec(g->get_capacity());

        return retcode;
    } // read_graph

/**
 * Private function to read in a graph in adjacency list format
 * \param[in,out] g a graph object
 * \param[in] filename full path to the file to read in (relative or absolute)
 * \return code 0 on success, nonzero on failure
 */
    int GraphReader::read_adjlist(Graph *g, const string filename){
        int new_vid;
        int new_nbr;

        string line;
        ifstream input(filename.c_str());
        if(!input.is_open()){
            FERROR("%s:  Error opening %s for reading\n", __FUNCTION__, filename.c_str());
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        filename.c_str());
        }

        while(std::getline(input, line)){
            std::stringstream ss(line);
            string token;
            ss >> token; // the first element
            if('#' != token[0]){
                istringstream(token) >> new_vid;
                if(new_vid > g->get_num_nodes() - 1){
                    g->add_vertices((g->get_num_nodes() - new_vid) + 1);
                }
                while( ss >> token){
                    istringstream(token) >> new_nbr;
                    if(new_nbr > g->get_num_nodes() - 1){
                        g->add_vertices((new_nbr - g->get_num_nodes()) + 1);
                    }
                    g->add_edge(new_vid, new_nbr);
                }
            }
        }
        return 0;
    } // read_adjlist

/**
 * Private function to read in a graph in edgelist format
 * \param[in,out] g a graph object
 * \param[in] filename full path to the file to read in (relative or absolute)
 * \return code 0 on success, nonzero on failure
 */
    int GraphReader::read_edgelist(Graph *g, const string filename){
        char line[100];
        int i, j, m, n, retval;
        char *retp;
        FILE *in;

        m = n = 0;

        if((in = fopen(filename.c_str(), "r")) == NULL){
            FERROR("%s:  Error opening %s for reading\n", __FUNCTION__, filename.c_str());
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        filename.c_str());
        }

        // get the first non-comment lines
        retp = fgets(line, 100, in);
        while(line[0] == '#'){
            retp = fgets(line, 100, in);
            if(NULL == retp){
                FERROR( "%s:  Edgelist read error - got EOF before any useful data\n",
                        __FUNCTION__);
            }
        }
        // it has out number of nodes and number of edges
        retval = sscanf(line, "%d %d", &n, &m);

        if(2 != retval){
            FERROR( "%s:  Edgelist read error - unable to read number of nodes and/or edges\n",
                    __FUNCTION__);
        }

        g->add_vertices(n);

        while(!feof(in)){
            retp = fgets(line, 100, in);

            // not sure why this doesn't trigger correctly in the while condition....
            if(feof(in)){
                break;
            }

            if((line[0] == '#') || (line[0] == '\n') ){
                continue;
            }

            // Edge line - start end
            retval = sscanf(line, "%d %d", &i, &j);

            // TODO: add smome detection of comments, blank lines, etc
            if(retval != 2){
                FERROR( "%s:  Edgelist read error - didn't understand edge line\n",
                        __FUNCTION__);
            }

            // Check start and end values
            if((i < 0) || (i > n - 1) || (j < 0) || (j > n - 1)){
                FERROR( "%s:  Edgelist read error - edge (u,v) (%d-%d) must have"
                        " both u and v between 0 and n=%d, inclusive\n",
                        __FUNCTION__, i, j, n - 1);
            }
            //fprintf(stderr, "i=%d j=%d\n", i, j);

            // So when we read in an edge (a,b) in the Edgelist file
            // we need to add G.label_index[b] to nodes[G.label_index[a]]
            // When going back at the end, we need to use the label field
            // The edge appears valid - store in adj lists

            g->add_edge(i, j);
        }
        ;

        DEBUG("Finish reading Edgelist file\n");
        //fprintf(stderr, "EdgelistGraphReader: read in %d edges\n", this->num_edges);

        // FIXME : how do we label properly???
        /* for(i = 0; i < this->capacity; i++){
            if(nodes[i].get_label() == -1){
                nodes[i].set_label(i + 1);
            }
           }

         */
        fclose(in);
        return 0;
    } // read_edgelist

/**
 * Private function to read in a graph in DIMACS format
 * \param[in,out] g a graph object
 * \param[in] filename full path to the file to read in (relative or absolute)
 * \param[in] read_vertex_weights read vertex weights from input file
 * \return code 0 on success, nonzero on failure
 */
    int GraphReader::read_dimacs(Graph *g, const string filename, bool read_vertex_weights){
        char line[100], format[100];
        int i, j, m, n, retval, id, x;
        int val;
        FILE *in;

        m = n = 0;

        if((in = fopen(filename.c_str(), "r")) == NULL){
            FERROR("%s:  Error opening %s for reading\n", __FUNCTION__, filename.c_str());
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        filename.c_str());
        }

        while(!feof(in)){
            retval = fscanf(in, "%2c", line);

            if(feof(in)){
                break;
            }

            if(retval != 1){
                FERROR("%s:  fscanf read %d char's (expected to read 1!)\n",
                       __FUNCTION__, retval);
            }
            switch(line[0])
            {
            case 'p':
                // This is the "problem line" - get n and m here
                // Make sure we don't already know n and m!
                if((n != 0) || (m != 0) ){
                    FERROR(
                        "%s:  DIMACS read error - are there more than one problem lines?\n",
                        __FUNCTION__);
                }
                retval = fscanf(in, "%s %d %d", format, &n, &m);

                DEBUG("DIMACS: read n = %d, m = %d\n", n, m);

                // Simple error checking
                if((n <= 0) || (m <= 0) || (retval != 3) ){
                    FERROR(
                        "%s:  DIMACS read error - didn't understand problem line! p %s %d %d\n",
                        __FUNCTION__, format, n, m);
                    exit(-1);
                }

                if((strncasecmp(format, "e", 1) != 0) && (strncasecmp(format, "edge", 4) != 0)){
                    FERROR(
                        "%s:  DIMACS read error - problem line - FORMAT must be one of\n"
                        " e, E, edge, or EDGE\n", __FUNCTION__);
                }

                g->add_vertices(n);

                /* here we assume that we initialize all weights to 1 */
                if(read_vertex_weights){
                    VertexWeightedGraph *vg = (VertexWeightedGraph *)g;
                    if((int)vg->weight.size() != n){
                        vg->weight.resize(n, 1);
                    }
                }

                break;

            case 'c':
                // Comment line - skip and move on
                break;

            case 'n':
                // Node line - of the form n id VALUE
                retval = fscanf(in, "%d %d", &id, &val);

                // Simple error checking - make sure we know n and m already!
                if((n == 0) || (m == 0) ){
                    FERROR("%s:  DIMACS read error - node line found before problem line!\n",
                           __FUNCTION__);
                }
                if(retval != 2){
                    FERROR( "%s:  DIMACS read error - didn't understand node line\n",
                            __FUNCTION__);
                }
                // Check id
                if((id < 1) || (id > n) ){
                    FERROR("%s:  DIMACS read error - node id (%d) must be between 1 and n= %d\n",
                           __FUNCTION__, id, n);
                }
                // Store value - first allocate the value array if NULL
                if(read_vertex_weights){
                    VertexWeightedGraph *vg = (VertexWeightedGraph *)g;
                    if((int)vg->weight.size() != n){
                        vg->weight.resize(n);
                    }

                    // Store value
                    // FIXME:: figure out how weighted graphs work
                    vg->weight[id - 1] = val;                 // 1->0
                    // DEBUG("Storing weight[%d]= %d\n", id - 1, this->weights[id - 1]);
                }

                break;

            case 'e':
                // Edge line - of the form e start end
                retval = fscanf(in, "%d %d", &i, &j);

                //DEBUG("Read edge %d-%d\n", i, j);

                // Simple error checking - make sure we know n and m already!
                if((n == 0) || (m == 0) ){
                    FERROR(
                        "%s:  DIMACS read error - edge line found before problem line!\n",
                        __FUNCTION__);
                }
                if(retval != 2){
                    FERROR( "%s:  DIMACS read error - didn't understand edge line\n",
                            __FUNCTION__);
                }
                // Check start and end values
                if((i < 1) || (i > n) || (j < 1) || (j > n) ){
                    FERROR( "%s:  DIMACS read error - edge (u,v) (%d-%d) must have"
                            " both u and v between 1 and n=%d, inclusive\n",
                            __FUNCTION__, i, j, n);
                }

                // So when we read in an edge (a,b) in the DIMACS file
                // we need to add G.label_index[b] to nodes[G.label_index[a]]
                // When going back at the end, we need to use the label field
                // The edge appears valid - store in adj lists

                g->add_edge(j - 1, i - 1);
                break;
            case 'x':
                // Shows up in TSP files - we don't need it
                break;

            default:
                // We encountered some character that we don't understand
                FERROR("%s:  DIMACS read error - didn't understand the character %c to start a line\n",
                       __FUNCTION__, line[0]);
                break;
            };                 // end switch

            // Advance to next line if there is format at the end of it
            while(!feof(in) && (x = getc(in)) != '\n'){
                ;
            }
        }
        fclose(in);
        DEBUG("Finish reading DIMACS file\n");

        return 0;
    }     // read_dimacs

/**
 * Private function to read in a graph in adjacency matrix format.  This is an
 * INDDGO-specific format - no claims are made to compatibility with other files.
 *
 * \param[in,out] g a graph object
 * \param[in] filename full path to the file to read in (relative or absolute)
 * \return code 0 on success, nonzero on failure
 */
    int GraphReader::read_adjmatrix(Graph *g, const string filename){
        int i, j, k, m, n, retval;
        FILE *in;

        if((in = fopen(filename.c_str(), "r")) == NULL){
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        filename.c_str());
        }

        fscanf(in, "%d\n", &n);
        //printf("%s:  Graph has %d nodes\n", __FUNCTION__, n);
        // Each row will have n binary entries, and there will be n rows.

        g->add_vertices(n);

        i = 0;             // Use i to keep track of row
        m = 0;             // To count edges
        while(!feof(in)){
            if(i >= n){
                break;
            }
            for(k = 0; k < n; k++){
                retval = fscanf(in, "%d ", &j);
                //printf("%s:  Read in %d\n", __FUNCTION__, j);
                if(retval != 1){
                    fatal_error("Didn't read digit properly from matrix\n");
                }
                if((j != 0) && (j != 1) ){
                    fatal_error("%s:  Encountered entry %d!=0 or 1!\n",
                                __FUNCTION__, j);
                }
                if(j == 1){
                    //ignore self loops!
                    if(i != k){
                        if(k > i){  // assume we have already seen this edge
                            g->add_edge(i, k);
                        }
                        m++;
                    }
                }
            }
            i++;
        }
        fclose(in);

        if(2 * (m / 2) != m){
            fatal_error(
                "%s: Found an odd number of ones (%d). Non-symmetric matrix. Please correct.\n",
                __FUNCTION__, m);
        }

        return 0;
    }     // read_adjmatrix

/**
 * Private function to read in a graph in METIS format
 * \param[in,out] g a graph object
 * \param[in] filename full path to the file to read in (relative or absolute)
 * \return code 0 on success, nonzero on failure
 */
    int GraphReader::read_metis(Graph *g, const string filename){
        // here we assume metis graphs are based on 1.
        string s;
        ifstream in(filename.c_str());
        vector<int> elem;
        int n, e = 0;
        int count = 0;

        if(in.is_open()){
            elem.clear();
            do {
                getline(in, s);
            } while(s[0] == '%');       // skip comments
            split(s, ' ', elem);
            n = elem[0];
            e = elem[1];

            fprintf(stderr, "n is %d\n", n);

            g->add_vertices(n);

            if((n <= 0) || (e <= 0) ){
                fatal_error("Metis read error - didn't understand problem line!");
            }

            while(in.good()){
                elem.clear();
                getline(in, s);
                if(s[0] != '%'){
                    split(s, ' ', elem);
                    int is = (int) elem.size();
                    for(int i = 0; i < is; i++){
                        if(!g->is_edge(count,elem[i] - 1)){
                            g->add_edge(count,elem[i] - 1);
                        }
                    }
                    count++;
                }
            }

            in.close();
        }
        return 0;
    }     // read_metis
}
