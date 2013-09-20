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

#ifdef HAS_SLEPC
  #include <slepceps.h>
#endif

using namespace std;

void print_time(ofstream &of, string prefix, ORB_t start, ORB_t end){
    cout << prefix + ": " << ORB_seconds(end, start) << endl;
    if(of.is_open()){
        of << prefix + ": " << ORB_seconds(end, start) << endl;
    }
}

const string allowed_methods ("edge_density,avg_degree,degree_dist,global_cc,avg_cc,local_ccs,shortest_paths,assortativity,eccentricity,eccentricity_dist,expansion,apsp_output,avg_shortest_path,shortest_paths_boost,eigen_spectrum,k_cores,degeneracy,betweenness,powerlaw,delta_hyperbolicity");

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
    cerr << "Usage: " << argv[0] << " [-h] -i infile [-t input-type] [-o outfile] [-p output-prefix] [-m methods] [-s eigen spectrum size] [-r]" << endl;
    cerr << "Allowed methods: " << allowed_methods << endl;
    cerr << "Input type should be one of: edge, adjlist, adjmatrix, dimacs" << endl;
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

int parse_options(int argc, char **argv, string& infile, string& intype, string& outfilename, string &outprefix, std::map<string, bool>& methods, bool& record_timings, bool &file_append, int *spectrum_spread){
    int flags, opt;
    while((opt = getopt(argc, argv, "hi:t:o:m:p:s:ra")) != -1){
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
        case 'r':
            record_timings = true;
            break;
        case 's':
            *spectrum_spread = atoi(optarg);
            break;
        case 'a':
            file_append = true;
            break;
        }
    }

    return 0;
} // parse_options

void run_all_methods(Graph::Graph *g, ofstream &outfile, ofstream &timing_file, string outprefix, std::map<string, bool> req_methods,  bool &file_append, int spectrum_spread){
    Graph::GraphReader gr;
    Graph::GraphProperties gp;
    Graph::GraphUtil gu;
    ORB_t t1, t2;

    double global_cc, avg_cc, assortativity;

    vector<double> local_cc, freq_ecc, norm_hops, eigen_spectrum;
    float edge_density, avg_degree;
    vector<int> deg_dist, ecc;
    int degeneracy;
    vector<int> k_cores;
    double avg_path_length;
    int xmin;
    double alpha, KS, max_delta;
    vector<vector<double> > delta;
    vector<double> betweenness;

    Graph::Graph *largest_component;

    vector< vector<int> > shortest_path_distances;

    cout << "Simplifying graph" << endl;
    ORB_read(t1);
    gp.make_simple(g);
    ORB_read(t2);
    print_time(timing_file, "Time(make_simple)", t1, t2);
    int num_components = g->get_num_connected_components();
    outfile << "connected_components " << num_components << endl;

    if(req_methods["edge_density"] == true){
        cout << "Calculating edge density" << endl;
        ORB_read(t1);
        gp.edge_density(g, edge_density);
        ORB_read(t2);
        print_time(timing_file, "Time(edge_density)", t1, t2);
        outfile << "edge_density " << edge_density << endl;
    }
    if(req_methods["avg_degree"] == true){
        cout << "Calculating average degree" << endl;
        ORB_read(t1);
        gp.avg_degree(g, avg_degree);
        ORB_read(t2);
        print_time(timing_file, "Time(average_degree)", t1, t2);
        outfile << "avg_degree " << avg_degree << endl;
    }
    if(req_methods["degree_dist"] == true){
        cout << "Calculating degree distribution" << endl;
        ORB_read(t1);
        gp.deg_dist(g, deg_dist);
        ORB_read(t2);
        print_time(timing_file, "Time(degree_distribution)", t1, t2);
        string of = outprefix + ".deg_dist";
        write_degree_distribution(of, deg_dist);
        outfile << "degree_distribution " <<  of << endl;
    }
    if(req_methods["assortativity"] == true){
        cout << "Calculating degree assortativity" << endl;
        ORB_read(t1);
        gp.deg_assortativity(g, assortativity);
        ORB_read(t2);
        print_time(timing_file, "Time(assortativity)", t1, t2);
        outfile << "assortativity " <<  assortativity << endl;
    }
    if((req_methods["degeneracy"] == true) || (req_methods["k_cores"] == true)){
        cout << "Calculating k_cores and degeneracy" << endl;
        ORB_read(t1);
        degeneracy = gu.find_kcore(g, &k_cores);
        ORB_read(t2);
        print_time(timing_file, "Time(find_kcore)", t1, t2);
        outfile << "degeneracy " << degeneracy << endl;
        if(req_methods["k_cores"] == true){
            string of = outprefix + ".kcores";
            outfile << "kcore file " << of << endl;
            write_kcores(of, k_cores);
        }
    }

    if((req_methods["global_cc"] == true) || (req_methods["local_ccs"] == true) || (req_methods["avg_cc"] == true)){
        cout << "Calculating clustering coefficients" << endl;
        ORB_read(t1);
        gp.clustering_coefficients(g, global_cc, avg_cc, local_cc);
        ORB_read(t2);
        print_time(timing_file, "Time(clustering_coeffecients)", t1, t2);
        if(req_methods["global_cc"] == true){
            outfile << "global_clustering_coefficient " << global_cc << endl;
        }
        if(req_methods["avg_cc"] == true){
            outfile << "average_clustering_coefficient " << avg_cc << endl;
        }
        if(req_methods["local_ccs"] == true){
        }
    }
    if(req_methods["shortest_paths"] == true){
        cout << "Calculating shortest paths" << endl;
        ORB_read(t1);
        gp.paths_dijkstra_all(g, shortest_path_distances);
        ORB_read(t2);
        print_time(timing_file, "Time(shortest_paths_dijkstra)", t1, t2);
    }

    #ifdef HAS_BOOST
    if((req_methods["shortest_paths_boost"] == true)){
        cout << "Creating BOOST representation of g" << endl;
        ORB_read(t1);
        gu.populate_boost(g);
        ORB_read(t2);
        print_time(timing_file, "Time(populate_boost)", t1, t2);
        cout << "Calculating shortest paths (boost)" << endl;
        ORB_read(t1);
        gp.paths_dijkstra_boost_all(g, shortest_path_distances);
        ORB_read(t2);
        print_time(timing_file, "Time(shortest_paths_dijkstra_boost)", t1, t2);
        if((req_methods["apsp_output"] == true)){
            string of = outprefix + ".apsp";
            write_apsp_matrix(of, shortest_path_distances);
            ORB_read(t2);
            print_time(timing_file, "Time(write_apsp_matrix)", t1, t2);
        }
    }
    if(req_methods["betweenness"]){
        /* cout << "Creating BOOST representation of g" << endl;
           ORB_read(t1);
           gu.populate_boost(g);
           ORB_read(t2);
           print_time(timing_file, "Time(populate_boost)", t1, t2);
         */cout << "Calculating betweeneess centrality" << endl;
        ORB_read(t1);
        gp.betweenness_centrality(g, betweenness);
        ORB_read(t2);
        print_time(timing_file, "Time(betweenness_centrality",t1,t2);
        string of = outprefix + ".betweenness";
        outfile << "betweenness_file " << of << endl;
        write_betweenness(of, g->get_betweenness_ref());
    }
    #else // ifdef HAS_BOOST
    cerr << "Error: BOOST support was not compiled, cannot run shortest_paths_boost or betweenness" << endl;
    #endif // ifdef HAS_BOOST

    if(num_components == 1){
        if(req_methods["eccentricity"] == true){
            cout << "Calculating eccentricities" << endl;
            ORB_read(t1);
            gp.eccentricity(g, ecc);
            ORB_read(t2);
            print_time(timing_file, "Time(eccentricity)",t1,t2);
        }
        if(req_methods["eccentricity_dist"] == true){
            cout << "Calculating distribution of eccentricities" << endl;
            ORB_read(t1);
            gp.eccentricity_dist(g, ecc, freq_ecc);
            ORB_read(t2);
            print_time(timing_file, "Time(eccentricity distribution)",t1,t2);
        }
    }
    else {
        cout << "Graph is disconnected - not calculating eccentricities" << endl;
    }

    if(req_methods["expansion"] == true){
        cout << "Calculating normalized expansion (distance distribution) - no self loops allowed" << endl;
        ORB_read(t1);
        gp.expansion(g, norm_hops);
        ORB_read(t2);
        print_time(timing_file, "Time(expansion)",t1,t2);
    }
    if(req_methods["avg_shortest_path"] == true){
        cout << "Calculating average shortest path length" << endl;
        ORB_read(t1);
        gp.avg_path_length(g, avg_path_length);
        ORB_read(t2);
        print_time(timing_file, "Time(avg_path_length)", t1, t2);
        outfile << "avg_path_length " << avg_path_length << endl;
    }

    #ifdef HAS_PETSC
    if(req_methods["eigen_spectrum"] == true){
        //If petsc/slepc are present, initalize those.
        //If MPI support is added in the future, init MPI before Petsc. Petsc will do it's own MPI
        //init if MPI isn't already inited.
        #ifdef HAS_SLEPC
        SlepcInitializeNoArguments();
        #elif HAVE_PETSC
        PetscInitializeNoArguments();
        #endif
        if(spectrum_spread == 0){
            spectrum_spread = 3;
        }

        cout << "Calculating adjacency matrix eigen spectrum\n";
        ORB_read(t1);
        gp.eigen_spectrum(g, eigen_spectrum, spectrum_spread);
        ORB_read(t2);
        print_time(timing_file, "Time(eigen spectrum)",t1,t2);
        outfile << "eigen_spectrum ";
        for(int idx = 0; idx < eigen_spectrum.size(); idx++){
            outfile << eigen_spectrum[idx];
        }
        outfile << "\n";
    }
    #endif // ifdef HAS_PETSC

    #ifdef HAS_BOOST
    if(req_methods["powerlaw"] == true){
        cout << "Calculating power law parameters" << endl;
        ORB_read(t1);
        gp.powerlaw(g, xmin, alpha, KS);
        ORB_read(t2);

        print_time(timing_file, "Time(powerlaw)", t1, t2);
        outfile << "powerlaw " << xmin << " " << alpha << " " << KS << endl;
    }
    #else
    cerr << "Error: BOOST support was not compiled, cannot run shortest_paths_boost or betweenness" << endl;
    #endif //HAS_BOOST
    if(num_components == 1){
        if(req_methods["delta_hyperbolicity"] == true){
            cout << "Calculating delta hyperbolicity" << endl;
            ORB_read(t1);
            gp.delta_hyperbolicity(g, max_delta, delta);
            ORB_read(t2);

            print_time(timing_file, "Time(delta_hyperbolicity)", t1, t2);
            //outfile << "delta_hyperbolicity " << max_delta << endl;
            //for(int idx = 0; idx < delta.size(); idx++){
            //    for(int jdx = 0; jdx < delta[idx].size(); jdx++){
            //        outfile << delta[idx][jdx] << " ";
            //    }
            //    outfile << endl;
            //}
            string of = outprefix + ".delta_hyp";
            write_delta_hyperbolicity(of, delta);
            outfile << "max_delta_hyperbolicity " << max_delta;
        }
    }
    else {
        cout << "Graph is disconnected - not calculating delta hyperbolicity" << endl;
    }

    outfile.close();

    #ifdef HAS_SLEPC
    SlepcFinalize();
    #elif HAVE_PETSC
    PetscFinalize();
    #endif
} // run_all_methods

int main(int argc, char **argv){
    string infile;
    string outfilename ("graph-stats.txt");
    string outprefix;
    ofstream outfile;
    ofstream timing_file;
    bool record_timings = false;
    bool file_append = false;
    string intype ("edge");
    std::map<string, bool> req_methods;
    std::map<string, bool> val_methods;
    ORB_t t1, t2;
    int spectrum_spread = 0;
    create_map(allowed_methods, val_methods);
    parse_options(argc, argv, infile, intype, outfilename, outprefix, req_methods, record_timings, file_append, &spectrum_spread);
    if(outprefix.length() == 0){
        outprefix = infile;
    }

    // we'd like higher precision when printing values
    std::cout.precision(10);

    cout << "done parsing options" << endl;
    cout << "Input  file: " << infile << endl;
    cout << "Input  type: " << intype << endl;
    cout << "Output file: " << outfilename << endl;
    cout << "Appending  : ";
    cout << std::boolalpha << file_append << endl;
    cout << "Methods    :";
    for(map<string, bool>::iterator it = req_methods.begin(); it != req_methods.end(); ++it){
        cout << " " << it->first;
        if(val_methods[it->first] != true){
            cerr << "Error: " << it->first << " is not a valid method! " << endl;
        }
    }
    cout << endl;
    cout << "Calibrating timers" << endl;
    ORB_calibrate();

    // let's do some calculations

    Graph::Graph *g = new(Graph::Graph);
    Graph::GraphReader gr;
    Graph::GraphProperties gp;
    Graph::GraphUtil gu;

    // Set up output streams
    if(file_append == true){
        outfile.open(outfilename.c_str());
    }
    else {
        outfile.open(outfilename.c_str(), ios_base::out | ios_base::app);
    }
    if(!outfile.is_open()){
        cerr << "Error opening " << outfilename << " for writing, exiting" << endl;
        exit(1);
    }

    if(outfile.tellp() == 0){
        outfile << "filename" << infile << endl;
        outfile << "num_nodes " << g->get_num_nodes() << endl;
        outfile << "num_edges " << g->get_num_edges() << endl;
    }

    if(record_timings){
        string of = outfilename + ".timings";
        timing_file.open(of.c_str());
        if(!timing_file.is_open()){
            cerr << "Error opening " << timing_file << " for writing, exiting" << endl;
            exit(1);
        }
        outfile << "Recording timings to " << of << endl;
    }

    // Read in the graph and start recording things to output streams
    cout << "Reading graph" << endl;
    ORB_read(t1);
    if(gr.read_graph(g, infile, intype, false) == -1){
        exit(1);
    }
    ORB_read(t2);
    print_time(timing_file, "Time(read_graph)", t1, t2);

    outfile.precision(16);
    outfile << "filename " << infile << endl;
    vector<int> components;
    cout << "Graph is connected?: " << std::boolalpha << gp.is_connected(g) << endl;
    ORB_read(t1);
    cout << "GU.label_all_components says: " << gu.label_all_components(g, &components) << endl;
    ORB_read(t2);
    print_time(timing_file, "Time(label_all_components)", t1, t2);
    bool is_connected = gp.is_connected(g);
    cout << "Connected components: " << g->get_num_connected_components() << endl;
    cout << "Graph is connected?: " << std::boolalpha << is_connected << endl;

    cout << "Creating largest connected component graph" << endl;
    run_all_methods(g, outfile, timing_file, outprefix, req_methods, file_append, spectrum_spread);
    outfile.close();
    timing_file.close();

    // some algorithms only make sense to run on a connected graph/component
    if(not is_connected){  // run everything against the other algorithms
        cout << "Graph is not connected, re-running stats on largest connected component" << endl;
        outfilename += ".largest_component";
        if(file_append == true){
            outfile.open(outfilename.c_str());
        }
        else {
            outfile.open(outfilename.c_str(), ios_base::out | ios_base::app);
        }
        if(!outfile.is_open()){
            cerr << "Error opening " << outfilename << " for writing, exiting" << endl;
            exit(1);
        }
        if(outfile.tellp() == 0){
            outfile << "filename " << infile << endl;
            outfile << "num_nodes " << g->get_num_nodes() << endl;
            outfile << "num_edges " << g->get_num_edges() << endl;
        }
        if(record_timings){
            string of = outfilename + ".timings";
            timing_file.open(of.c_str());
            if(!timing_file.is_open()){
                cerr << "Error opening " << timing_file << " for writing, exiting" << endl;
                exit(1);
            }
            outfile << "Recording timings to " << of << endl;
        }

        outprefix = outprefix + ".largest_component";

        outfile.precision(16);
        outfile << "filename" << infile << endl;
        Graph::Graph *largest_component = gu.get_largest_component_graph(g);
        delete(g);  // delete g here to save on memory

        run_all_methods(largest_component, outfile, timing_file, outprefix, req_methods, file_append, spectrum_spread);
        outfile.close();
        timing_file.close();
    }

    exit(0);
} // main

