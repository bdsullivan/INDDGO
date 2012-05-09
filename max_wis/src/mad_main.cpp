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

#ifdef __MADNESS__
#include <iostream>
#include <string.h>
#include <assert.h>

#include <world/world.h>
#include <world/worlddc.h>

#include "GraphDecomposition.h"
#include "TreeDecomposition.h"
#include "Log.h"
#if WIN32 || _WIN32
#define GLPSOL_EXE "glpsol.exe"
#else
#define GLPSOL_EXE "glpsol"
#endif
#define WRAP_UTHASH 0

#include "weighted_ind_set.h"
#include "WeightedMutableGraph.h"
#define GLPSOL "glpsol"

#include "MadnessTableProcessor.h"
#include <unordered_map>
#include "PTD.h"
#include "GraphException.h" 

using namespace madness;
using namespace std;

static void prepare(TDTree *T, int k, bool x)
{
    if (x)
        DP_prepare_node(T, k);
}

int main(int argc, char** argv)
{
  initialize(argc, argv);
  World world(MPI::COMM_WORLD);
  // For debugging uncomment following line.
  //redirectio(world);
  char *env;
  
  Graph::WeightedMutableGraph *G;
  TDTree *T = NULL;
  DP_info info;
  list<int> optimal_solution;
  
  char lfname[100];
  char efname[100];
  
  int ptd_threads;
  int mad_threads;
  int omp_threads;
  
  int rank = world.rank();
  ostringstream os;

  try{
    // Process command line options
    int help = info.process_DP_info(argc, argv);
    if(help == -1){
      if(rank == 0){
	usage(argv[0]);
      }
      finalize();
      exit(1);
    }
  }
  catch(Graph::GraphException& e)
    {
      cerr << "exception caught: " << e.what() << endl;
      finalize(); 
      exit(-1);
    }


  env = getenv("MAD_NUM_THREADS");
  if(env != NULL){
    mad_threads = atoi(env);
  } else {
    fprintf(stderr,"ERROR: MAD_NUM_THREADS must be set in your environment\n");
    finalize();
    exit(1);
  }
  
  env = getenv("PTD_NUM_THREADS");
  if(env != NULL){
    ptd_threads = atoi(env);
  } else {
    fprintf(stderr,"ERROR: PTD_NUM_THREADS must be set in your environment\n");
    finalize();
    exit(1);
  }
  env = getenv("OMP_NUM_THREADS");
  if(env != NULL){
    omp_threads = atoi(env);
  } else {
    fprintf(stderr,"ERROR: OMP_NUM_THREADS must be set in your environment\n");
    finalize();
    exit(1);
  }
  
  
  //DEBUG("%d:%d:x%d:%f:%f:%d\n", max, world.size(), nthreads,
  //etime, prepare_node_time, y);
  os << world.size() << ":";
  os << "x" << ptd_threads << ":";
  
  sprintf(lfname, "madlog-%04d", rank);
  //0: debug, 5: critical
  LOG_INIT(lfname, NULL, 0);
  
  try{
  DEBUG("create graph\n");
  // Create the graph for WIS
  create_WIS_graph(&info, G);
  
  //world.gop.fence();
  // See if the MIP is to be written/solved
  if (info.write_mod)
    {
      write_ind_set_model(info.DIMACS_file, info.model_file, G);
      if (info.solve_mip)
	{
	  // Use GLPSOL to run the MIP solver
	  char command_string[100];
	  sprintf(command_string, "%s -m %s", GLPSOL, info.model_file);
	  fprintf(
		  stderr,
		  "Will use %s to solve MIP. Solution will be written to %s.MIP.WIS.sol\n\n",
		  GLPSOL, info.DIMACS_file);
	  system(command_string);
	}
      return 1;
    }

    //DEBUG("create td\n");

	// Create the tree decomposition using the options
    double tdtime = cpu_time();
    if (!info.pbag)
      create_tree_decomposition(&info, G, &T, true);
    else
    {
        world.gop.fence();
        PTD *p = new PTD(world, *G);
        p->set_file_name(p->get_random_name(info.DIMACS_file).c_str());
        if (!info.parmetis)
            p->compute_elim_order(info.elim_order_type, info.read_ordering, info.ord_file);
        else
	  p->parmetis_elim_order(info.elim_order_type);
	
	// Write elimination ordering into a file
	if (info.eorder && rank == 0)
	{
            p->write_elim_order(info.eorder);
	}	

        p->triangulate();
        p->create_tree_decomposition(info.DIMACS_file, T);
	// Write out the tree if desired with root labelled 1, and rest
	//labelled according to pre-order walk
        if(rank == 0)
        {
                if (info.write_ordered_tree)
                        T->write_DIMACS_file(info.tree_outfile, true);
                else if (info.write_tree)
                        T->write_DIMACS_file(info.tree_outfile);
        }
    }
	world.gop.fence();

    if (info.decompose_only)
    {
        os << cpu_time() - tdtime << endl;
        if (rank == 0)
            fprintf(stdout, "%s\n", os.str().c_str());

        delete T;
        delete G;

        LOG_CLOSE();
        finalize();
        return 0;
    }

    os << cpu_time() - tdtime << ":";
    os << "x" << omp_threads << ":";

    srand(time(NULL));
    long r = rand();
    world.gop.max(r);

	MadnessTableProcessor *mtp = new MadnessTableProcessor(world, G->get_num_nodes(), T->num_mask_words, r);

    //DEBUG("populate dc\n");

    double prepare_node_time =  0.0;
    double pstart = 0.0;

    pstart = cpu_time();
#pragma omp parallel for
	for (int i = 0; i < T->num_tree_nodes; i++)
	{
        prepare(T, i, (mtp->owner(i) == world.rank()));
    }
	world.gop.fence();
    prepare_node_time += cpu_time() - pstart;
    os << prepare_node_time << ":";

	for (int i = 0; i < T->num_tree_nodes; i++)
	{
		if (mtp->owner(i) == world.rank())
		{
            if (i == T->root_node)
			{
				T->tree_nodes[i]->set_root(true);
			}
			mtp->add_node(i, *(static_cast<TDMadTreeNode*>(T->tree_nodes[i])));
		}
	}

    int root_node = T->root_node;
	delete T;
	delete G;

	world.gop.fence();    

    DEBUG("start computation\n");
    int max = 0;
    double etime = 0.0;
    double ttime = 0.0;
    double tstart = cpu_time();
	if (world.rank() == mtp->owner(root_node))
	{
        Future<TDMadTreeNode> fsol =  mtp->task(mtp->owner(root_node),
                  &MadnessTableProcessor::compute_table, root_node);
        max = mtp->task(mtp->owner(root_node),
                          &MadnessTableProcessor::find_max, fsol);
	}

    world.gop.fence();
    etime = cpu_time() - pstart;

    int y = getHWmem();
    world.gop.max(max);
    world.gop.max(y);

    /*Output remaining solution information*/

    os << "x" << mad_threads << ":";
    os << etime << ":";
    os << y << ":";
    os << max << endl;
    
    GEN("%s:", info.DIMACS_file);
    GEN("%s\n", os.str().c_str());
    if (rank == 0)
    {
        fprintf(stdout,  "%s:", info.DIMACS_file);
        fprintf(stdout, "%s\n", os.str().c_str());
    }
    delete mtp;
    LOG_CLOSE();
    finalize();
    return 0;
  }
  catch(Graph::GraphException& e)
    {
      cerr << "exception caught: " << e.what() << endl;
      LOG_CLOSE();
      finalize(); 
      exit(-1);
    }

}
#else
int main()
{
    printf("To use madness compile with -D__MADNESS__ flag\n");
    return 0;
}
#endif  /* __MADNESS__ */
