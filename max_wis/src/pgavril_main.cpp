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

#include <iostream>
#include <string.h>
#include <assert.h>
#include "GraphDecomposition.h"
#include "TreeDecomposition.h"
#include "Log.h"

#include "PTD.h"
#include "weighted_ind_set.h"
#include "WeightedMutableGraph.h"
#define GLPSOL "glpsol"

using namespace std;

int main(int argc, char** argv)
{
	initialize(argc, argv);
	World world(MPI::COMM_WORLD);
	redirectio(world);

	Graph::WeightedMutableGraph *G;
    PTD  *p;

	DP_info info;

    double xtime = MPI_Wtime();
    double otime = MPI_Wtime();

	char lfname[100];
	char efname[100];

    int rank = world.rank();

	sprintf(lfname, "log.%05d", rank);
	//0: debug, 5: critical
	LOG_INIT(lfname, NULL, 0);

	// Process command line options
	int help = info.process_DP_info(argc, argv);
	if(help == -1){
	  if(rank == 0){
	    usage(argv[0]);
	  }
	  finalize();
	  exit(1);
	}

	// Create the graph for WIS
	create_WIS_graph(&info, G);

	// Create the tree decomposition using the options
    xtime = MPI_Wtime();


	p = new PTD(world, *G);
    p->set_file_name(p->get_random_name(info.DIMACS_file).c_str());
    p->compute_elim_order(info.elim_order_type, info.read_ordering, info.ord_file);


    DEBUG("eot: %f\n", MPI_Wtime() - xtime);
    DEBUG("trangulate\n");
    xtime = MPI_Wtime();

    xtime = MPI_Wtime();
    
    // processing bags using MPI and pthreads.
    if (!info.width)
    {
        // trangulation
        int width = p->triangulate();
        DEBUG("tt: %d\n", MPI_Wtime() - xtime);
        DEBUG("graph: %s ord:%s set width: %d \n", info.DIMACS_file, info.ord_file, width);
        
        DEBUG("starting data distribution\n");

        p->process_bag();

        DEBUG("finish data distribution\n");
        
        DEBUG("ptime: %f\n", MPI_Wtime() - xtime);
        
        xtime = MPI_Wtime();
        //barrier
        world.gop.fence();
        
        DEBUG("ftime: %f\n", MPI_Wtime() - xtime);
        DEBUG("total time: %f\n", MPI_Wtime() - otime);
        
        if (p->rank() == 0)
            p->merge_files(info.DIMACS_file);
    }
    else
    {
        int width = p->triangulate();        
        DEBUG("%s: %d\n", info.DIMACS_file, width);
        if (p->rank() == 0)
            fprintf(stdout, "%s: %d\n", info.DIMACS_file, width);
    }

	delete G;
    delete p;
    
	LOG_CLOSE();
    finalize();
	return 0;
}
