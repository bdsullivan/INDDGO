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
#include "TDDynamicProgramming.h"

#include <weighted_ind_set.h>
#include <parallel_wis.h>
#include "Log.h"

#define WRAP_UTHASH 0
#define GLPSOL "glpsol"

#ifdef __PARALLEL__

int main(int argc, char **argv)
{
	int i, j;
	Graph::WeightedMutableGraph *G;
	TDTree *T=NULL;
	DP_info info;
	list<int> optimal_solution;
    int rank, size;
    double stime, etime;

    int storage_nodes = 1;
    
    MPI_Group orig_group, new_group;

    double leaf_time;
    double nonleaf_time;
    long smem_hwm;
    long rmem_hwm;

    

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    // There is a seperate log file for each processor, strictly for
    // debugging and development purposes only.

    char lfname[100];
    char efname[100];
    sprintf (lfname, "parallel_wis-%d.log", rank);
    sprintf (efname, "parallel_wis-error.log", rank);
    //0: debug, 5: critical
    LOG_INIT(lfname, efname, 0);


     if (!(size > storage_nodes + 2))
     {
         FERROR ("No enough processors, please allocate more\n");
         FERROR ("Processor distribution: head_node: %d storage_nodes: %d compute_nodes: %d\n", 1, storage_nodes, size - (storage_nodes + 1));
         MPI_Finalize ();
         exit (0);
     }


	// Process command line options
	int help = info.process_DP_info(argc, argv);
	if(help == -1){
	  if(rank == 0)
	    {
	      usage(argv[0]);
	    } 
	  MPI_Finalize ();
	  exit (1);
	}
	
	// Create the graph for WIS
	create_WIS_graph(&info, G);
	// See if the MIP is to be written/solved
	if(info.write_mod)
	{
		write_ind_set_model(info.DIMACS_file, info.model_file, G);
		if(info.solve_mip)
		{
			// Use GLPSOL to run the MIP solver
			char command_string[100];
			sprintf(command_string,"%s -m %s",GLPSOL,info.model_file);
			fprintf(stderr,"Will use %s to solve MIP. Solution will be written to %s.MIP.WIS.sol",
                    GLPSOL,info.DIMACS_file);
			system(command_string);
		}
		return 1;		
	}
		
	
	// Create the tree decomposition using the options
	create_tree_decomposition(&info, G, &T);
	// Check to see if the DP is to be done
	if(info.decompose_only)
	{
		printf("%s: Treewidth %d\n", info.DIMACS_file, T->width);
		//delete G;
		//delete T;
		return 1;
	}


    T->head_node = 0;
    T->allnodes = new vector<int>(size, -1);
    //T->storage_nodes = new vector<int>(size, -1);
    //T->compute_nodes = new vector<int>(size, -1);

    (*T->allnodes)[0] = HEAD_NODE;
    for (i = 1; i < size; i++)
    {
        if ((storage_nodes > 0) && (i % 4 == 1))
        {
            T->storage_nodes.push_back(i);
            (*T->allnodes)[i] = STORAGE_NODE;
            storage_nodes --;
        }
        else
        {
            T->compute_nodes.push_back(i);
            (*T->allnodes)[i] = COMPUTE_NODE;
        }
    }

    if (MPI_SUCCESS != 
        MPI_Comm_group(MPI_COMM_WORLD, &orig_group))
    {
        FERROR ("MPI Test any failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    if (MPI_SUCCESS != 
        MPI_Group_incl (orig_group, T->compute_nodes.size(), &T->compute_nodes.front(), &new_group))
    {
        FERROR ("MPI Test any failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    MPI_Comm compute_nodes_comm;
    if (MPI_SUCCESS != 
        MPI_Comm_create (MPI_COMM_WORLD, new_group, &compute_nodes_comm))
    {
        FERROR ("MPI Test any failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    T->compute_nodes_comm = &compute_nodes_comm;
    T->storage_nodes_size = T->storage_nodes.size();
    T->compute_nodes_size = T->compute_nodes.size();

    if ((*T->allnodes)[rank] == HEAD_NODE)
    {
        parallel_wis_head_init (T, size, rank);
        parallel_wis_head_start (T);
    }
    else if ((*T->allnodes)[rank] == COMPUTE_NODE)
    {
        parallel_wis_compute_init (T, size, rank);
    }
    else if ((*T->allnodes)[rank] == STORAGE_NODE)
    {
        parallel_wis_storage_init (T, size, rank);
        parallel_wis_storage_start (T);
    }


    if ((*T->allnodes)[T->rank] == COMPUTE_NODE)
    {
        DEBUG("post order walk\n");
        // Create the post_order walk
        vector<int> walk(T->num_tree_nodes,GD_UNDEFINED);
        T->post_order_walk(&walk);

        stime = MPI_Wtime();
        int j;
        DEBUG("starting loop\n");
        for(i = 0; i < T->num_tree_nodes; i++)
        {
            CRIT("processing node : ############# %d\n", walk[i]);
            DEBUG(" %d  more nodes left\n", (T->num_tree_nodes - i));
            T->compute_table(compute_weighted_ind_sets_parallel, walk[i]);
            MPI_Barrier (compute_nodes_comm);

            if(T->info->free_children)
            {
                for(j = 0; (j < T->storage_nodes_size) && (T->rank == T->compute_nodes[0]); j++)
                {
                    if (MPI_SUCCESS !=
                        MPI_Send ((void *)&walk[i], 1, MPI_INT, T->storage_nodes[j], MPI_CHILD_FREE_TAG, MPI_COMM_WORLD))
                    {
                        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }
                }
            }
        }

        // Send termination signal to head node to nofity that
        // computation is over for the given tree
        if (MPI_SUCCESS != 
            MPI_Send (&T->rank, 1, MPI_INT, 0, MPI_COMP_TERM_TAG, MPI_COMM_WORLD))
        {
            FERROR ("MPI Send failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
        
        // Set flag_recv_done when all tree nodes have been finished then
        // receiving threads can exit.
        T->flag_recv_done = 1;
        
        etime = MPI_Wtime ();
        

        // Store and reset info's table stats 
        info.orig_total_pc_table_entries=info.total_pc_table_entries;
        info.orig_total_table_entries=info.total_table_entries;
        
        info.total_pc_table_entries=info.total_table_entries=0;
        
        
        // Compute some stats here in parent-child intersection to analyze 
        // why memory savings varies across different graphs
        // This is less than ideal, but you have to do it this way for something like make_nice
        // since there are "holes" in the tree_nodes[] array.  The above compute_table() loop
        // gets around this by using the post order walk that is filled with only valid tree node
        // indices!
        
        vector<int> intersection_sizes(T->tree_nodes.size(),-1);
        T->compute_parent_child_intersections(&intersection_sizes);
        info.avg_pc_proportion=0;
        int num_in_avg=0;
        for(i=0;i<(int)T->tree_nodes.size();i++)
        {
            if(T->tree_nodes[i])
            {
                num_in_avg++;
                info.avg_pc_proportion += ((double)intersection_sizes[i]/(double)T->tree_nodes[i]->bag.size());
            }
        }
        info.avg_pc_proportion = info.avg_pc_proportion/(double)num_in_avg;
    }

    smem_hwm = getHWmem();
    if ((*T->allnodes)[rank] == COMPUTE_NODE)
    {
    
        if (MPI_SUCCESS != 
            MPI_Reduce ((void *)&T->info->leaf_time, (void *)&leaf_time, 1, MPI_DOUBLE, MPI_MAX, 0, *T->compute_nodes_comm))
        {
            FERROR ("MPI Reduce failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
        T->info->leaf_time = leaf_time;
    
        if (MPI_SUCCESS != 
            MPI_Reduce ((void *)&T->info->nonleaf_time, (void *)&nonleaf_time, 1, MPI_DOUBLE, MPI_MAX, 0, *T->compute_nodes_comm))
        {
            FERROR ("MPI Reduce failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
        T->info->nonleaf_time = nonleaf_time;
    
        if (MPI_SUCCESS != 
            MPI_Reduce ((void *)&smem_hwm, (void *)&rmem_hwm, 1, MPI_LONG, MPI_MAX, 0, *T->compute_nodes_comm))
        {
            FERROR ("MPI Reduce failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
    
        if (T->rank == T->compute_nodes[0])
            fprintf (stderr, "Memory HWM for all compute nodes: %d\n", rmem_hwm);
    }
    
    if (MPI_SUCCESS != 
        MPI_Reduce ((void *)&smem_hwm, (void *)&rmem_hwm, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD))
    {
        FERROR ("MPI Reduce failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    if (T->rank == 0)
        fprintf (stderr, "Memory HWM for all  nodes: %d\n", rmem_hwm);


    
    if ((*T->allnodes)[rank] == HEAD_NODE)
    {
        parallel_wis_head_cleanup (T);
    }
    else if ((*T->allnodes)[rank] == COMPUTE_NODE)
    {
        parallel_wis_compute_cleanup (T);
    }
    else if ((*T->allnodes)[rank] == STORAGE_NODE)
    {
        parallel_wis_storage_cleanup (T);
    }

    LOG_CLOSE();
    MPI_Finalize();

    if (T->rank == T->compute_nodes[0] )
        print_WIS_results(stdout, T, &info);

  	//delete T;
    //delete &G;
 	return 1;
}
#else
int main(int argc, char **argv)
{
    fprintf (stderr, "please compile with -D__PARALLEL__ flag to use parallel program\n");
    exit (0);
}
#endif  // __PARALLEL__





