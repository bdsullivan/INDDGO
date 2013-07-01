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
<<<<<<< HEAD
#include <fstream>
=======
#include <map>
>>>>>>> upstream/master

#include "Graph.h"
#include "GraphReader.h"
#include "GraphProperties.h"
#include "Debug.h"
#include "Log.h"
#include "Util.h"
#include "GraphException.h"
#include "GraphUtil.h"
#include "GraphDecomposition.h"
#include "DIMACSGraphReader.h"

#include "orbconfig.h"
#include "orbtimer.h"

using namespace std;

void print_time(string prefix, ORB_t start, ORB_t end){
    cout << prefix + ": " << ORB_seconds(end, start) << "\n";
}

const string allowed_methods ("edge_density,avg_degree,degree_dist,global_cc,avg_cc,local_ccs,shortest_paths,assortativity,eccentricity,eccentricity_dist,expansion");

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

void print_usage(char **argv){
    cerr << "Usage: " << argv[0] << " [-h] -i infile [-t input-type] [-o outfile] [-p output-prefix] [-m methods]\n";
    cerr << "Allowed methods: " << allowed_methods << "\n";
    cerr << "Input type should be one of: edge, adjlist, adjmatrix, dimacs\n";
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
<<<<<<< HEAD
    while((opt = getopt(argc, argv, "i:t:o:m:")) != -1){
	switch(opt){
	case 'i':
	    infile = optarg;
	    break;
	case 't':
	    intype = optarg;
	    break;
	case 'o':
	    outfile = optarg;
	    break;
	case 'm':
	    string tempstr;
	    string method;
	    string token;
	    tempstr = optarg;
	    stringstream convert(tempstr);
	    /* split the string on commas */
	    while(getline(convert, token, ',')){
		methods.push_back(token);
	    }
	    break;
	}
=======
    while((opt = getopt(argc, argv, "hi:t:o:m:p:")) != -1){
        switch(opt){
        case 'h':
            print_usage(argv);
            exit(0);
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
>>>>>>> upstream/master
    }
    
    return 0;
} // parse_options

void write_results(vector<int> *kshell, string infile,  string outfile) {
    // TODO: make more robust
    ofstream myfile;
    myfile.open(outfile.c_str());
    myfile << "#This file contains the color generated from the file " << infile;
    myfile << " the colors are provided in ascending nodeID order\n";
    myfile << "#min " << *min_element(kshell->begin(), kshell->end()) << "\n";
    myfile << "#max " << *max_element(kshell->begin(), kshell->end()) << "\n";
    for(int i = 0; i < kshell->size(); i++)
	myfile << (*kshell)[i] << " "; 
    myfile << "\n";
    myfile.close();
}

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
<<<<<<< HEAD
    for(vector<string>::iterator it = methods.begin(); it != methods.end(); ++it){
	cout << " " << *it;
    }
    cout << "\n";
    
    Graph::GraphReader read;
    Graph::GraphUtil util;  
    Graph::GraphProperties prop;
    
    time_t start, stop;
    
    Graph::Graph temp;
    read.read_graph(&temp, infile, "DIMACS", false);
    cout << "n: " << temp.get_capacity() << endl;
    int max_deg = 0;
    int isolated = 0;
    
    for(int i = 0; i < methods.size(); i++) {
	Graph::Graph g;
	vector<int> kshell(10);  
	int k_degen = -1;
	read.read_graph(&g, infile, "DIMACS", false);
	prop.make_simple(&g);
	max_deg = 0;
	
	for(int j= 0; j < g.get_capacity(); j++) {
	    if(g.get_degree()[j] > max_deg)
		max_deg = g.get_degree()[j];
	    if(g.get_degree()[j] == 0)
		isolated++;
	}
	cout << "max_deg: " << max_deg << endl;
	cout << "isolated: " << isolated << endl;
	//Graph::create_largestcomponent_graph(infile.c_str(), &g);
	
	
	
	if(methods[i].compare("1") == 0) {
	    start = clock();
	    //      k_degen = util.find_kcore(&g, &kshell);  
	    stop = clock();
	}
	else if(methods[i].compare("2") == 0) {
	    start = clock();
	    k_degen = util.find_kcore2(&g, &kshell);  
	    stop = clock();
	}
	else if (methods[i].compare("3") == 0) {
	    start = clock();
	    k_degen = util.find_degen(&g, &kshell);
	    stop = clock();
	}
	ofstream output;
	output.open("results");
	output << "kcore: " << k_degen << endl;
	output << "max_deg: " << max_deg << endl;
	output << "isolated: " << isolated << endl;
	output.close();
	
	cout << "Algorithm " << methods[i] << " runtime: " << (((double)(stop - start )) / CLOCKS_PER_SEC) << endl;
	write_results(&kshell, infile, outfile+methods[i]);
	cout << "k_degen" << methods[i] << ":\t" << k_degen << endl;
    }
}
=======
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
    vector<double> local_cc, freq_ecc, norm_hops;
    float edge_density, avg_degree;
    vector<int> deg_dist, ecc;
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
    if(req_methods["eccentricity"] == true){
        cout << "Calculating eccentricities\n";
        ORB_read(t1);
        gp.eccentricity(&g, ecc);
        ORB_read(t2);
        print_time("Time(eccentricity)",t1,t2);
    }
    if(req_methods["eccentricity_dist"] == true){
        cout << "Calculating distribution of eccentricities\n";
        ORB_read(t1);
        gp.eccentricity_dist(&g, ecc, freq_ecc);
        ORB_read(t2);
        print_time("Time(eccentricity distribution)",t1,t2);
    }
    if(req_methods["expansion"] == true){
        cout << "Calculating normalized expansion (distance distribution) - no self loops allowed\n";
        ORB_read(t1);
        gp.expansion(&g, norm_hops);
        ORB_read(t2);
        print_time("Time(expansion)",t1,t2);
    }

    outfile.close();
    exit(0);
} // main

>>>>>>> upstream/master
