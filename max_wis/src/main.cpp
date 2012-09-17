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
#include "Log.h"
#include "GraphException.h"
#include <fstream>
#include <string>

// Include header for getpid and define GLPSOL_EXE if we are trying to solve
// the problem with GLPK via integer programming
#if WIN32 || _WIN32
#define GLPSOL_EXE "glpsol.exe"
#include <process.h>
#else
#define GLPSOL_EXE "glpsol"
#include <unistd.h>
#endif
#define WRAP_UTHASH 0

#include "weighted_ind_set.h"
#include "WeightedMutableGraph.h"

int main(int argc, char **argv)
{
	int i;
	Graph::WeightedMutableGraph *G;
	TDTree *T = NULL;
	DP_info info;
	list<int> optimal_solution;
	int pid;
#if WIN32 || _WIN32
	pid=_getpid();
#else
	pid=(int)getpid();
#endif
	char lfname[100];
	char efname[100];
	sprintf(lfname, "weighted_ind_set-%d.log",pid);
	sprintf(efname, "weighted_ind_set-%d.log",pid);
	//0: debug, 5: critical
	LOG_INIT(lfname, efname, 0);

    try
    {
        // Process command line options
      int help = info.process_DP_info(argc, argv);
      if(help == -1){
	usage(argv[0]); 
	exit(1);
      }

		if(info.fix_DIMACS)
		{
			char fixed_DIMACS_file[100];
			sprintf(fixed_DIMACS_file,"%s.fixed",info.DIMACS_file);
			fprintf(stderr,"Writing fixed DIMACS file to %s\n",fixed_DIMACS_file);
			normalize_DIMACS_file(info.DIMACS_file,fixed_DIMACS_file);
			fprintf(stderr,"%s contains normalized DIMACS file\n",fixed_DIMACS_file);
			exit(0);
		}
        // Create the graph for WIS
        create_WIS_graph(&info, G);

        // See if the MIP is to be written/solved
        if (info.write_mod)
        {
            write_ind_set_model(info.DIMACS_file, info.model_file, G);
            if (info.solve_mip)
            {
                // Use GLPSOL to run the MIP solver
                char command_string[100];
                sprintf(command_string, "%s -m %s", GLPSOL_EXE, info.model_file);
                fprintf(
                    stderr,
                    "Will use %s to solve MIP. Solution will be written to %s.MIP.WIS.sol\n\n",
                    GLPSOL_EXE, info.DIMACS_file);
                system(command_string);
            }
            return 1;
        }

        // Create the tree decomposition using the options
        create_tree_decomposition(&info, G, &T);

        if (info.width)
        {
            printf("%s: Treewidth %d\n", info.DIMACS_file, T->width);
            delete G;
            delete T;
            return 1;
        }

        int max_tree_nodes=T->tree_nodes.size();


	// Check to see if decompose_only and no memory estimation
        if (info.decompose_only && !info.mem_est)
        {
            printf("%s: Treewidth %d\n", info.DIMACS_file, T->width);
	    delete G;
            delete T;
            return 1;
        }


        // Create the post_order walk - gives user the option of doing memory estimate
        // and then not doing the DP - the mem estimate requires the walk
        vector<int> walk(T->num_tree_nodes, GD_UNDEFINED);
        if (T->info->DFS)
            T->post_order_walk_DFS(&walk);
        else
            T->post_order_walk(&walk);

        if (T->info->verbose)
        {
            print_message(0, "WALK:\n");
            int walk_size=(int) walk.size();
            for (i = 0; i < walk_size; i++)
                print_message(0, "%d ", walk[i]);
            print_message(0, "\n");
        }
        // Check to see if we are generating the mem_est
        if(info.mem_est)
        {
            char mem_output[100];
            sprintf(mem_output,"%s.mem_est",info.DIMACS_file);
            info.mem_estimate=estimate_memory_usage(T, &walk, (const char *)mem_output);
        }

        // Check to see if decompose_only (
        if (info.decompose_only)
        {
            printf("%s: Treewidth %d\n", info.DIMACS_file, T->width);
	    delete G;
            delete T;
            return 1;
        }

        // Do the DP using the provided options - this loop is where all the time is!!!
        // Note that T->num_tree_nodes is safe here since the walk
        // handles NULLs
        double dstart = clock();
        for (i = 0; i < T->num_tree_nodes; i++)
            T->compute_table(compute_weighted_ind_set_table, walk[i]);
        //print_message(0, "%.2f\n", (clock() - dstart) / CLOCKS_PER_SEC);
        //printf("%s: Treewidth %d\n", info.DIMACS_file, T->width);

        // Store and reset info's table stats 
        info.orig_total_pc_table_entries = info.total_pc_table_entries;
        info.orig_total_table_entries = info.total_table_entries;
        info.total_pc_table_entries = info.total_table_entries = 0;

        // Compute some stats here in parent-child intersection to analyze 
        // why memory savings varies across different graphs
        vector<int> intersection_sizes(max_tree_nodes, -1);
        T->compute_parent_child_intersections(&intersection_sizes);
        info.avg_pc_proportion = 0;
        int num_in_avg = 0;
        for (i = 0; i < max_tree_nodes; i++)
        {
            if (T->tree_nodes[i])
            {
                num_in_avg++;
                info.avg_pc_proportion += ((double) intersection_sizes[i]
                                           / (double) T->tree_nodes[i]->bag.size());
            }
        }
        info.avg_pc_proportion = info.avg_pc_proportion / (double) num_in_avg;

        // Reconstruct the solution if desired.
        list<int> root_intersection, root_difference;
        T->refined_width = 0;
        if (info.no_reconstruct)
        {
            fprintf(stderr,"Getting optimal_obj\n");
            get_optimal_obj(T, &root_intersection, &root_difference, &info);
            info.reconstruct_time = 0;
            fprintf(stderr,"Not reconstructing\n");
            print_WIS_results(stdout, T, &info);
        }
        else
        {
            clock_t recon_stop, recon_start = clock();
            // Need to reconstruct the optimal solution - if we don't have the 
            // tables, then refine the tree and re-run the DP
            if (!info.free_children)
            {
                // We have the children's tables so just reconstruct the optimal solution
                info.opt_obj = reconstruct_solution(T, &optimal_solution,
                                                    info.sol_file);
                recon_stop = clock();
                info.reconstruct_time = (double) (recon_stop - recon_start)
                    / CLOCKS_PER_SEC;

                // Print out results to stdout
                print_WIS_results(stdout, T, &info);
            }
            else
            {
                // Find the best solution and the root intersection
                get_optimal_obj(T, &root_intersection, &root_difference, &info);
                // Free the root's children and the root's table
                T->free_children(T->root_node);
                T->free_table(T->root_node);
                // Refine the tree by removing the difference between the root bag and the
                // optimal solution with root bag from all bags in the tree 
                T->refine_tree(&root_difference);
                if(T->info->verbose)
                    print_message(0, "Root difference has %d entries\n",
                                  (int) root_difference.size());
                list<int> opt_sol_nbrs;
                list<int>::iterator ii, jj;
                Graph::Node *n1;
                list<int> *nbr_list;
                for (ii = root_intersection.begin(); ii != root_intersection.end(); ++ii)
                {
                    n1=T->G->get_node(*ii);
                    nbr_list=n1->get_nbrs_ptr();
                    for (jj = nbr_list->begin(); jj!= nbr_list->end(); ++jj)
                    {
                        opt_sol_nbrs.push_back(*jj);
                    }
                    opt_sol_nbrs.sort();
                    opt_sol_nbrs.unique();
                }
                if(T->info->verbose)
                    print_message(0, "Removing %d nbrs from tree\n",
                                  (int) opt_sol_nbrs.size());
                T->refine_tree(&opt_sol_nbrs);
                // Run the DP again on the smaller tree
                T->num_tree_nodes_processed = 0;
                T->info->opt_obj = 0;
                T->info->free_children = false;
                T->info->parent_child = false;

                // async table update disable for reconstruction
                T->info->async_tbl = false;
                T->info->verbose = false;

                for (i = 0; i < T->num_tree_nodes; i++)
                {
                    T->compute_table(compute_weighted_ind_set_table, walk[i]);
                    if ((int) T->tree_nodes[walk[i]]->bag.size() > T->refined_width)
                        T->refined_width = (int) T->tree_nodes[walk[i]]->bag.size();
                }
                T->refined_width--;
                // Now reconstruct the optimal solution
                T->info->opt_obj = reconstruct_solution(T, &optimal_solution,
                                                        info.sol_file);

                recon_stop = clock();
                T->info->reconstruct_time = (double) (recon_stop - recon_start)
                    / CLOCKS_PER_SEC;
                // Print out results to stdout
                // Reset info.free_children!
                T->info->free_children = true;
                print_WIS_results(stdout, T, T->info);
            }
        }
    }
    catch(Graph::GraphException& e)
    {
        cerr << "exception caught: " << e.what() << endl;
        LOG_CLOSE();
        return -1;
    }

	LOG_CLOSE();


	//BDS 04/05 - always returning 0 to allow successful calls in make.
	// return the opt obj val if we have it
	int retval;
	//if(T->info->opt_obj>0)
	//		retval=T->info->opt_obj;
	//else
	retval= 0;

	delete T;
	delete G;

	return retval;
}


