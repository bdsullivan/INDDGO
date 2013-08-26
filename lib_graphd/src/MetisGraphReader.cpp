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

#include "MetisGraphReader.h"
#include "Debug.h"
#include "Log.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

namespace Graph {
    MetisGraphReader::MetisGraphReader(){
    }

    MetisGraphReader::MetisGraphReader(string f) :
        filename(f){
    }

    MetisGraphReader::~MetisGraphReader(){
    }

    void MetisGraphReader::split(const string& s, char c, vector<int>& v){
        string::size_type i = 0;
        string::size_type j = s.find(c);
        stringstream ss;
        int k;
        while(j != string::npos){
            //v.push_back(s.substr(i, j - i));
            k = 0;
            ss.clear();
            ss << s.substr(i, j - i);
            ss >> k;
            v.push_back(k);
            i = ++j;
            j = s.find(c, j);
            if(j == string::npos){
                k = 0;
                //v.push_back(s.substr(i, s.length( )));
                ss.clear();
                ss << s.substr(i, s.length());
                ss >> k;
                if(k > 0){
                    v.push_back(k);
                }
            }
        }
    } // split

    void MetisGraphReader::read_graph(const char *filename){
        // here we assume metis graphs are based on 1.
        string s;
        ifstream in(filename);
        vector<int> elem;
        int n, e = 0;
        int count = 0;

        if(in.is_open()){
            elem.clear();
            getline(in, s);
            split(s, ' ', elem);
            n = elem[0];
            e = elem[1];

            if((n <= 0) || (e <= 0) ){
                fatal_error("Metis read error - didn't understand problem line!");
            }

            this->weights.resize(n, 1);
            this->degree.resize(n, 0);
            this->nodes.resize(n);
            this->capacity = n;

            while(in.good()){
                getline(in, s);
                elem.clear();
                split(s, ' ', elem);
                int is = (int) elem.size();
                for(int i = 0; i < is; i++){
                    this->nodes[count].add_nbr(elem[i] - 1);
                    this->degree[count]++;
                }
                count++;
            }

            in.close();
        }
    } // read_graph

    void MetisGraphReader::read_graph(){
        this->read_graph(filename.c_str());
    }

    vector<int> MetisGraphReader::get_degree(){
        return degree;
    }

    vector<Node> MetisGraphReader::get_nodes(){
        return nodes;
    }

    int MetisGraphReader::get_num_edges() const {
        return num_edges;
    }

    vector<int> MetisGraphReader::get_weights(){
        return weights;
    }

    int MetisGraphReader::get_capacity() const {
        return capacity;
    }

    void MetisGraphReader::set_capacity(int capacity){
        this->capacity = capacity;
    }
}
