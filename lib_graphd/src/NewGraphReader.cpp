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

#include "NewGraphReader.h"
#include "Log.h"
#include "Debug.h"
#include <strings.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

namespace Graph {
    NewGraphReader::NewGraphReader(){
    }

    NewGraphReader::~NewGraphReader(){
    }

    MutableGraph *NewGraphReader::read_graph(const char *filename, char *type) {
        MutableGraph *g;
        // FIXME: use regexes here
        fprintf(stderr,"- in read_graph type is %s\n", type);
        if(strcasecmp("edge", type) == 0){
            fprintf(stderr,"- calling read_edgelist\n");
            g = NewGraphReader::read_edgelist(filename);
        }
        else if(strncasecmp("adjmatrix",type , 9) == 0){
            fprintf(stderr,"- calling read_adjmatrix\n");
            g = NewGraphReader::read_adjmatrix(filename);
        }
        else if(strncasecmp("metis", type, 5) == 0){
            fprintf(stderr,"- calling read_metis\n");
            g = NewGraphReader::read_metis(filename);
        } else {
            fprintf(stderr, "WTF\n");
        }
        return g;
    }

    MutableGraph *NewGraphReader::read_edgelist(const char *filename){
        char line[100];
        int i, j, m, n, retval, count;
        char *retp;
        FILE *in;
        MutableGraph *g;

        m = n = 0;
        count = 0;
        // Use the above to count # of connected nodes, etc. and make sure that it
        // matches what is in the existing Graph struct.

        if((in = fopen(filename, "r")) == NULL){
            FERROR("%s:  Error opening %s for reading\n", __FUNCTION__, filename);
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        filename);
        }

        // get the first non-comment lines
        retp = fgets(line, 100, in);
        while(line[0] == '#'){
            retp = fgets(line, 100, in);
        }
        // it has out number of nodes and number of edges
        retval = sscanf(line, "%d %d", &n, &m);
        g = new MutableGraph(n);

        while(!feof(in)){
            retp = fgets(line, 100, in);

            // not sure why this doesn't trigger correctly in the while condition....
            if(feof(in)){
                break;
            }

            if(line[0] == '#'){
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
        };

        DEBUG("Finish reading Edgelist file\n");
        //fprintf(stderr, "EdgelistGraphReader: read in %d edges\n", this->num_edges);

        // FIXME : how do we label properly???
        /* for(i = 0; i < this->capacity; i++){
            if(nodes[i].get_label() == -1){
                nodes[i].set_label(i + 1);
            }
        }

        */
        return g;

    }

    MutableGraph *NewGraphReader::read_dimacs(const char *filename){
        char line[100], format[100], x;
        int i, j, m, n, retval, id, count;
        int val;
        FILE *in;
        MutableGraph *g;

        m = n = 0;
        count = 0;
        // Use the above to count # of connected nodes, etc. and make sure that it
        // matches what is in the existing Graph struct.

        if((in = fopen(filename, "r")) == NULL){
            FERROR("%s:  Error opening %s for reading\n", __FUNCTION__, filename);
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        filename);
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


                g = new MutableGraph(n);

                if((strncasecmp(format, "e", 1) != 0) && (strncasecmp(format, "edge", 4) != 0)){
                    FERROR(
                        "%s:  DIMACS read error - problem line - FORMAT must be one of\n"
                        " e, E, edge, or EDGE\n", __FUNCTION__);
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
                //if((int)this->weight.size()!=n)
                //    this->weight.resize(n);

                // Store value
                // FIXME:: figure out how weighted graphs work
                // this->weights[id - 1] = val;                 // 1->0
                // DEBUG("Storing weight[%d]= %d\n", id - 1, this->weights[id - 1]);

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

                g->add_edge(j-1, i-1);
                break;
            case 'x':
                // Shows up in TSP files - we don't need it
                break;

            default:
                // We encountered some character that we don't understand
                FERROR("%s:  DIMACS read error - didn't understand the character %c to start a line\n",
                       __FUNCTION__, line[0]);
                break;
            };             // end switch

            // Advance to next line if there is format at the end of it
            while(!feof(in) && (x = getc(in)) != '\n'){
                ;
            }
        }
        fclose(in);
        DEBUG("Finish reading DIMACS file\n");

        return g;
    } // read_dimacs

    MutableGraph *NewGraphReader::read_adjmatrix(const char *filename){
        int i, j, k, m, n, retval;
        FILE *in;
        MutableGraph *g;

        if((in = fopen(filename, "r")) == NULL){
            fatal_error("%s:  Error opening %s for reading\n", __FUNCTION__,
                        filename);
        }

        fscanf(in, "%d\n", &n);
        //printf("%s:  Graph has %d nodes\n", __FUNCTION__, n);
        // Each row will have n binary entries, and there will be n rows.

        g = new MutableGraph(n);

        i = 0;         // Use i to keep track of row
        m = 0;         // To count edges
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
                        g->add_edge(i, k);
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

        return g;
    } // read_adjmatrix

    void NewGraphReader::split(const string& s, char sep, vector<int>& v){
        string::size_type i = 0;
        stringstream ss(s);
        int k=-1;
        string temp;
        while(std::getline(ss, temp, sep)){
            stringstream convert(temp);
            convert >> k;
            if(k != -1){
                v.push_back(k);
            }
            k=-1;
        }
    } // split

    MutableGraph *NewGraphReader::read_metis(const char *filename){
        // here we assume metis graphs are based on 1.
        string s;
        ifstream in(filename);
        char c = ' ';
        vector<int> elem;
        MutableGraph *g;
        int n, e = 0;
        int count = 0;

        if(in.is_open()){
            elem.clear();
            do {
                getline(in, s);
            } while (s[0] == '%');  // skip comments
            split(s, ' ', elem);
            n = elem[0];
            e = elem[1];

            fprintf(stderr, "n is %d\n", n);
            if((n <= 0) || (e <= 0) ){
                fatal_error("Metis read error - didn't understand problem line!");
            }


            g = new MutableGraph(n);
            while(in.good()){
                elem.clear();
                getline(in, s);
                if(s[0] != '%') {
                    split(s, ' ', elem);
                    int is = (int) elem.size();
                    for(int i = 0; i < is; i++){
                        if(! g->is_edge(count,elem[i] - 1)){
                            g->add_edge(count,elem[i] - 1);
                        }
                    }
                    count++;
                }
            }

            in.close();
        }
        return g;
    } // read_metis
}


