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

#include "GraphDecomposition.h"
#include "TreeDecomposition.h"
#include <fstream>
#include <unistd.h>
#include <sstream>

/*
 * Usage message for td_stats function.
 */
void usage(){
    cerr << "/* required arguments:\n* -g graph_file (dimacs)\n* -s  kcore_file (scorefile format)\n* optional arguments:\n* -o output_file_basename (otherwise it writes to stdout)\n* -t tree_outfile_basename (otherwise no tree is written)\n* -l : This computes lower bounds on the treewidth \n* -L : Lower bounds ONLY. Turns off the computation of all upper bounds, bag distributions, etc\n*/\n";
}

/*
 * This file generates one or more tree decompositions and
 * calculates a variety of statistics (width, bag size distribution,
 * average bag score), as well as optional treewidth lower bounds.
 * It generates a single output file containing all requested
 * information in a csv format.
 * As usual, this operates on the largest connected component.
 */
int main(int argc, char **argv){
    /* required arguments:
     * -g graph_file (dimacs)
     * -s  kcore_file (scorefile format)
     * -o output_file_basename
     * -t tree_outfile_basename
     * optional arguments:
     * -l : This computes lower bounds on the treewidth
     * -L : Lower bounds ONLY. Turns off the computation of all upper bounds, bag distributions, etc
     */

    vector<double> *kcore_score = new vector<double>();
    vector<double> *degree_score = new vector<double>();
    int tdt[] = {
        TD_SUPERETREE, TD_GAVRIL
    };                                    //, TD_BK, TD_NICE};
    int et[] = {
        GD_AMD, GD_METIS_NODE_ND, GD_METIS_MMD
    };
    int lbt[] = {
        GD_MAX_MIN_DEGREE_LB, GD_MCS_LB
    };
    const char *tdtype[] = {
        "Super-E-Tree", "Gavril","Bodlaender-Koster", "Nice"
    };
    const char *elimtype[] = {
        "AMD", "MetisNodeND", "MetisMultMinD", "MinFill", "MCS"
    };
    bool lower_bounds = false;
    bool upper_bounds = true;
    bool all_methods = false;

    vector<int> td_types(tdt, tdt + sizeof(tdt) / sizeof(tdt[0]));
    vector<int> elim_types(et, et + sizeof(et) / sizeof(et[0]));
    vector<double> *st[] = {
        degree_score, kcore_score
    };
    vector<vector<double> *> scores(st, st + sizeof(st) / sizeof(st[0]));

    char *graph_file = NULL;
    char *kcore_file = NULL;
    char *out_file_base = NULL;
    char *tree_file_base = NULL;

    //process arguments
    int c;
    while((c = getopt (argc, argv, "hlLag:s:o:t:")) != -1){
        switch(c)
        {
        case 'h':
            usage();
            return 0;
        case 'L':
            if(all_methods){
                cerr << "Error: cannot specify both -L and -a options" << endl;
                exit(1);
            }
            upper_bounds = false;
            lower_bounds = true;
            break;
        case 'a':
            if(!upper_bounds){
                cerr << "Error: cannot specify both -L and -a options" << endl;
                exit(1);
            }
            all_methods = true;
            lower_bounds = true;
            break;
        case 'l':
            lower_bounds = true;
            break;
        case 'g':
            graph_file = optarg;
            break;
        case 's':
            kcore_file = optarg;
            break;
        case 'o':
            out_file_base = optarg;
            break;
        case 't':
            tree_file_base = optarg;
            break;
        case '?':
            if((optopt == 'g') || (optopt == 's') || (optopt == 'o') || (optopt == 't') ){
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            }
            else if(isprint (optopt)){
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            }
            else {
                fprintf (stderr,
                         "Unknown option character `\\x%x'.\n",
                         optopt);
            }
            return 1;
        default:
            abort ();
        }
    }
    for(int index = optind; index < argc; index++){
        printf ("Non-option argument %s\n", argv[index]);
    }

    int t, e, s;
    Graph::VertexWeightedGraph *G;
    TDTree *T;
    int treewidth, treelength, minecc;
    double kcore_max, kcore_min;
    std::ofstream out;
    std::ostream outStream(cout.rdbuf());
    std::stringstream sstm;
    int td_id;

    /* Add methods to the stack if -a was passed */
    if(all_methods){
        td_types.push_back(TD_BK);
        td_types.push_back(TD_NICE);
        elim_types.push_back(GD_MIN_FILL);
        elim_types.push_back(GD_MCS);
    }

    /*initialize log files if needed*/
    int pid;
    #if WIN32 || _WIN32
    pid = _getpid();
    #else
    pid = (int)getpid();
    #endif
    char lfname[100];
    char efname[100];
    sprintf(lfname, "stats-%d.log",pid);
    sprintf(efname, "stats-%d.log",pid);
    //0: debug, 5: critical
    LOG_INIT(lfname, efname, 0);

    try
    {
        if(graph_file == NULL){
            throw(Graph::GraphException("Failed to specify input graph file (-g). Aborting.\n"));
        }

        /*populate the graph*/
        Graph::create_largestcomponent_graph(graph_file, G);
        Graph::GraphUtil gutil;
        gutil.recompute_degrees(G);

        if(lower_bounds){
            /*
             * Open output file for writing results
             */
            if(out_file_base == NULL){
                outStream.rdbuf(std::cout.rdbuf());
            }
            else {
                sstm << out_file_base << "LB";
                out.open((sstm.str()).c_str(), fstream::out);
                if(out.fail()){
                    fatal_error("%s:  Error opening file %s for writing graphviz output\n", __FUNCTION__, (sstm.str()).c_str());
                }
                outStream.rdbuf(out.rdbuf());
                sstm.str(""); //clears the stream
            }

            outStream << "# " << graph_file << endl;
            outStream << "# Lower bounds: min_max_degree and mcs\n";

            Graph::GraphEOUtil eoutil;
            /*
             * Print the lower bounds
             */
            for(int i = 0; i < sizeof(lbt) / sizeof(lbt[0]); i++){
                //currently defaulting start_v to 0. Should probably improve this.
                outStream << eoutil.get_tw_lower_bound(G, lbt[i], 0 ) << "\t";
            }
            outStream << endl;

            /*
             * Close the output file.
             */
            out.close();
        }

        if(!upper_bounds){
            delete G;
            LOG_CLOSE();
            return 1;
        }

        if((kcore_file == NULL) || (tree_file_base == NULL) ){
            throw(Graph::GraphException("TD stats requested, but no score file (-s) specified. Aborting.\n"));
        }

        /*populate appropriate score vectors*/
        bool range = read_color_file(kcore_file,kcore_max,kcore_min,*kcore_score);
        vector<int> idegree_score = G->get_degree();
        (*degree_score).resize(idegree_score.size());
        for(int i = 0; i < idegree_score.size(); i++){
            (*degree_score)[i] = idegree_score[i];
        }

        /*loop over tree decomposition algorithms*/
        for(t = 0; t < td_types.size(); t++){
            /*loop over elimination order heuristics*/
            for(e = 0; e < elim_types.size(); e++){
                td_id = 10 * (t + 1) + e; //gives a unique number for the decomposition

                /*form the tree decomposition*/
                create_tree_decomposition(G, &T, false, NULL, false,
                                          false, NULL, elim_types[e],
                                          GD_UNDEFINED, td_types[t],
                                          false, true);
                /*create a place to store all the statistics for this particular decomposition*/
                vector<vector<double> > stats(T->num_tree_nodes);
                vector<double> mystats;
                vector<double>::iterator it;

                //fill the bag vectors
                T->fill_bag_vecs();
                cout << "T has " << T->num_tree_nodes << " tree nodes\n";

                //Non-score-specific statistics - width, length, eccentricity

                /* Width - uses degree scores for non-negativity check*/
                /* We define width = |B| - 1 to coincide with treewidth*/
                bag_statistics(T, *(scores[0]), &mystats, GD_STAT_COUNT);
                treewidth = 0;
                for(int i = 0; i < mystats.size(); i++){
                    if(mystats[i] - 1 > treewidth){
                        treewidth = (int)(mystats[i]) - 1;
                    }
                    stats[i].push_back(mystats[i] - 1);
                }

                /*Length*/
                vector<int> lengths;
                treelength = 0;
                bag_lengths(T, &lengths);
                for(int i = 0; i < lengths.size(); i++){
                    if(lengths[i] > treelength){
                        treelength = lengths[i];
                    }
                    stats[i].push_back((double)lengths[i]);
                }

                /*Eccentricity*/
                Graph::Graph *GT = T->export_tree();
                vector<int> ecc;
                minecc = INT_MAX;
                gutil.find_ecc(GT, &ecc);
                for(int i = 0; i < ecc.size(); i++){
                    if(ecc[i] < minecc){
                        minecc = ecc[i];
                    }
                    stats[i].push_back(ecc[i]);
                }

                /*loop over scores and calculate mean med stddev*/
                for(s = 0; s < scores.size(); s++){
                    /*Mean*/
                    bag_statistics(T, *(scores[s]), &mystats, GD_STAT_MEAN);
                    for(int i = 0; i < mystats.size(); i++){
                        stats[i].push_back(mystats[i]);
                    }

                    /*Median*/
                    bag_statistics(T, *(scores[s]), &mystats, GD_STAT_MED);
                    for(int i = 0; i < mystats.size(); i++){
                        stats[i].push_back(mystats[i]);
                    }

                    /*Standard Deviation*/
                    bag_statistics(T, *(scores[s]), &mystats, GD_STAT_STD);
                    for(int i = 0; i < mystats.size(); i++){
                        stats[i].push_back(mystats[i]);
                    }
                }

                /*
                 * Open output file for writing results
                 */
                if(out_file_base == NULL){
                    outStream.rdbuf(std::cout.rdbuf());
                }
                else {
                    sstm << out_file_base << td_id;
                    out.open((sstm.str()).c_str(), fstream::out);
                    if(out.fail()){
                        fatal_error("%s:  Error opening file %s for writing graphviz output\n", __FUNCTION__, (sstm.str()).c_str());
                    }
                    outStream.rdbuf(out.rdbuf());
                    sstm.str(""); //clears the stream
                }

                /*
                 * Print the file header
                 */
                outStream << "# " << graph_file << endl;
                outStream << "# " << tdtype[t] << endl;
                outStream << "# " << elimtype[e] << endl;
                outStream << "# " << "Width " << (treewidth - 1) << endl;
                outStream << "# " << "Length " << treelength << endl;
                outStream << "# " << "MinEcc " << minecc << endl;
                outStream << "# Bag\t Cardinality\t Length\t Eccentricity\t MeanDegree\t MedianDegree\t StdDevDegree\t MeanScore\t MedianScore\t StdDevScore" << endl;

                /*
                 * Print the statistics for each bag
                 */
                for(int i = 0; i < mystats.size(); i++){
                    outStream << i << "\t";
                    for(it = stats[i].begin(); it != stats[i].end(); ++it){
                        outStream << *it << "\t";
                    }
                    outStream << endl;
                }

                /*
                 * Close the output file.
                 */
                out.close();

                /*
                 * Write the tree decomposition to file, if required.
                 */
                if(tree_file_base != NULL){
                    sstm << tree_file_base << td_id;
                    T->write_DIMACS_file((sstm.str()).c_str());
                    sstm.str(""); //clears the stream
                }

                /*delete the tree decomposition*/
                delete T;
            }
        }

        delete G;
        LOG_CLOSE();
        return 1;
    }
    catch(Graph::GraphException& e)
    {
        cerr << "exception caught: " << e.what() << endl;
        LOG_CLOSE();
        return -1;
    }
} // main

