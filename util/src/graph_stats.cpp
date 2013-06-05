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

#include <stdlib.h>
#include <getopt.h>
#include <sstream>
#include <map>

#include "Graph.h"
#include "GraphReader.h"
#include "GraphProperties.h"
#include "Debug.h"
#include "Log.h"
#include "Util.h"
#include "GraphException.h"

#include "orbconfig.h"
#include "orbtimer.h"

using namespace std;

void print_time(string prefix, ORB_t start, ORB_t end){
    cout << prefix + ": " << ORB_seconds(end, start) << "\n";
}

const string allowed_methods ("edge_density,avg_degree,degree_dist,global_cc,avg_cc,local_ccs,shortest_paths,assortativity");

/**
 * Creates a map from a comma-separated string
 * \param[in] list comma-separated list of methods
 * \param[out] outmap a std::map ref with our newly separated methods
 */
void create_map(string list, map<string, bool> &outmap){
    string token;
    stringstream convert(list);
    while(getline(convert, token, ',')){
        outmap[token] = true;
    }
}

/**
 * Parse all of our options
 * \param[in] argc argument count
 * \param[in] argv arguments
 * \param[out] infile file to read as input
 * \param[out] intype type of intput file
 * \param[out] outfilename file to write output to
 * \param[out] methods list of methods we want to run.  Valid values currently: edge_density,avg_degree,degree_dist,global_cc, avg_cc, local_ccs
 */
int parse_options(int argc, char **argv, string& infile, string& intype, string& outfilename, string &outprefix, std::map<string, bool>& methods){
    int flags, opt;
    while((opt = getopt(argc, argv, "i:t:o:m:p:")) != -1){
        switch(opt){
        case 'i':
            infile = optarg;
            break;
        case 't':
            intype = optarg;
            break;
        case 'o':
            outfilename = optarg;
            break;
        case 'p':
            outprefix = optarg;
            break;
        case 'm':
            create_map(optarg, methods);
            break;
        }
    }

    return 0;
} // parse_options

int main(int argc, char **argv){
    string infile;
    string outfilename ("graph-stats.txt");
    string outprefix;
    ofstream outfile;
    string intype ("edge");
    std::map<string, bool> req_methods;
    std::map<string, bool> val_methods;
    ORB_t t1, t2;

    create_map(allowed_methods, val_methods);
    parse_options(argc, argv, infile, intype, outfilename, outprefix, req_methods);
    if(outprefix.length() == 0){
        outprefix = infile;
    }

    // we'd like higher precision when printing values
    std::cout.precision(10);

    cout << "done parsing options\n";
    cout << "Input  file: " << infile << "\n";
    cout << "Input  type: " << intype << "\n";
    cout << "Output file: " << outfilename << "\n";
    cout << "Methods    :";
    for(map<string, bool>::iterator it = req_methods.begin(); it != req_methods.end(); ++it){
        cout << " " << it->first;
        if(val_methods[it->first] != true){
            cerr << "Error: " << it->first << " is not a valid method! " << endl;
        }
    }
    cout << "\n";
    cout << "Calibrating timers\n";
    ORB_calibrate();

    // let's do some calculations

    Graph::Graph g;
    Graph::GraphReader gr;
    Graph::GraphProperties gp;

    cout << "Reading graph\n";
    ORB_read(t1);
    if(gr.read_graph(&g, infile, intype, false) == -1){
        exit(1);
    }
    ORB_read(t2);
    print_time("Time(read_graph)", t1, t2);

    double global_cc, avg_cc, assortativity;
    vector<double> local_cc;
    float edge_density, avg_degree;
    vector<int> deg_dist;
    vector< vector<int> > shortest_path_distances;

    outfile.open(outfilename.c_str());
    if(!outfile.is_open()){
        cerr << "Error opening " << outfilename << " for writing, exiting\n";
        exit(1);
    }

    outfile.precision(16);

    cout << "Simplifying graph\n";
    ORB_read(t1);
    gp.make_simple(&g);
    ORB_read(t2);
    print_time("Time(make_simple)", t1, t2);

    outfile << "num_nodes " << g.get_num_nodes() << "\n";
    outfile << "num_edges " << g.get_num_edges() << "\n";

    if(req_methods["edge_density"] == true){
        cout << "Calculating edge density\n";
        ORB_read(t1);
        gp.edge_density(&g, edge_density);
        ORB_read(t2);
        print_time("Time(edge_density)", t1, t2);
        outfile << "edge_density " << edge_density << "\n";
    }
    if(req_methods["avg_degree"] == true){
        cout << "Calculating average degree\n";
        ORB_read(t1);
        gp.avg_degree(&g, avg_degree);
        ORB_read(t2);
        print_time("Time(average_degree)", t1, t2);
        outfile << "avg_degree " << avg_degree << "\n";
    }
    if(req_methods["degree_dist"] == true){
        cout << "Calculating degree distribution\n";
        ORB_read(t1);
        gp.deg_dist(&g, deg_dist);
        ORB_read(t2);
        print_time("Time(degree_distribution)", t1, t2);
        string of = outprefix + ".deg_dist";
        write_degree_distribution(of, deg_dist);
        outfile << "degree_distribution " <<  of << "\n";
    }
    if(req_methods["assortativity"] == true){
        cout << "Calculating degree assortativity\n";
        ORB_read(t1);
        gp.deg_assortativity(&g, assortativity);
        ORB_read(t2);
        print_time("Time(assortativity)", t1, t2);
        outfile << "assortativity " <<  assortativity << "\n";
    }
    if((req_methods["global_cc"] == true) || (req_methods["local_ccs"] == true) || (req_methods["avg_cc"] == true)){
        cout << "Calculating clustering coefficients\n";
        ORB_read(t1);
        gp.clustering_coefficients(&g, global_cc, avg_cc, local_cc);
        ORB_read(t2);
        print_time("Time(clustering_coeffecients)", t1, t2);
        if(req_methods["global_cc"] == true){
            outfile << "global_clustering_coefficient " << global_cc << "\n";
        }
        if(req_methods["avg_cc"] == true){
            outfile << "average_clustering_coefficient " << avg_cc << "\n";
        }
        if(req_methods["local_ccs"] == true){
        }
    }

    if(req_methods["shortest_paths"] == true){
        cout << "Calculating shortest paths\n";
        ORB_read(t1);
        gp.paths_dijkstra_all(&g, shortest_path_distances);
        ORB_read(t2);
        print_time("Time(shortest_paths_dijkstra)", t1, t2);
    }
    outfile.close();
    exit(0);
} // main

