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
#include <fstream>

#include "Graph.h"
#include "Debug.h"
#include "Log.h"
#include "GraphException.h"
#include "GraphUtil.h"
#include "GraphDecomposition.h"
#include "DIMACSGraphReader.h"

using namespace std;

/**
 * Parse all of our options
 * \param[in] argc argument count
 * \param[in] argv arguments
 * \param[out] infile file to read as input
 * \param[out] intype type of intput file
 * \param[out] outfile file to write output to
 * \param[out] methods list of methods we want to run
 */
int parse_options(int argc, char **argv, string& infile, string& intype, string& outfile, vector<string>& methods){
    int flags, opt;
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
    string outfile;
    string intype;
    vector<string> methods;
    parse_options(argc, argv, infile, intype, outfile, methods);
    cout << "done parsing options\n";
    cout << "Input  file: " << infile << "\n";
    cout << "Input  type: " << intype << "\n";
    cout << "Output file: " << outfile << "\n";
    cout << "Methods    :";
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
