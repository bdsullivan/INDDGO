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

#include <mpi.h>
#include "MutableGraph.h"
#include "GraphEOUtil.h"
#include "Log.h"
#include "GraphCreatorFile.h"
#include "assert.h"

#ifdef HAS_PARMETIS
int testcount;
void parmetis_with_amd ()
{
    DEBUG("########## Testing Function : %s #################\n", __FUNCTION__);
    int size, rank;
    MPI_Comm comm;
	Graph::GraphEOUtil eoutil;
    Graph::VertexWeightedGraph *wmg;
	Graph::GraphCreatorFile creator;

	int output[] =
        {58,128,64,96,112,120,124,126,127,63,121,117,102,86,103,87,59,61,62,90,91,93,107,108,109,110,111,115,119,123,125,95,55,79,78,122,114,106,118,116,100,84,32,16,24,28,30,31,40,44,46,47,48,52,54,56,60,72,76,80,92,94,104,88,25,42,101,97,81,1,8,4,2,11,3,5,6,7,9,10,13,17,18,19,21,34,35,37,41,49,65,66,67,69,73,33,15,14,68,36,20,12,43,77,113,105,57,99,22,23,26,27,29,38,39,45,50,51,53,71,74,75,82,83,85,89,98,70};  

    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    creator.set_file_name("../data/1dc.128.txt");
    creator.set_graph_type("DIMACS");
    wmg = creator.create_weighted_mutable_graph();
    
	vector<int> ordering(wmg->get_num_nodes(), -1);
    eoutil.parmetis_elimination_ordering(wmg, ordering, GD_AMD, false, comm);
    if (rank == 0)
    {
        int n = wmg->get_num_nodes();
         for (int i = 0; i < n; i++)
             assert (output[i] == ordering[i]);
        cout << "Parmetis Test "<< ++testcount << " : Passed" << endl;
    }
}

void parmetis_only ()
{
    DEBUG("########## Testing Function : %s #################\n", __FUNCTION__);
    int size, rank;
    MPI_Comm comm;
	Graph::GraphEOUtil eoutil;
    Graph::VertexWeightedGraph *wmg;
	Graph::GraphCreatorFile creator;

	int output[] =
        {128,64,121,103,86,96,112,58,120,117,124,102,126,127,63,95,111,87,108,110,109,91,107,119,123,125,115,59,61,62,90,93,94,106,114,116,118,122,28,72,40,24,16,44,55,52,92,46,54,47,88,31,78,84,56,104,32,79,80,100,30,60,48,76,1,25,34,7,2,21,18,9,10,17,6,3,4,5,11,13,19,97,81,67,37,73,66,65,49,35,41,69,33,8,42,101,27,50,12,113,20,26,83,75,53,105,23,74,71,89,77,15,57,98,70,14,82,22,29,45,99,43,85,68,38,51,39,36};

    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    creator.set_file_name("../data/1dc.128.txt");
    creator.set_graph_type("DIMACS");
    wmg = creator.create_weighted_mutable_graph();
    
	vector<int> ordering(wmg->get_num_nodes(), -1);
    eoutil.parmetis_elimination_ordering(wmg, ordering, -1, false, comm);
    if (rank == 0)
    {
        int n = wmg->get_num_nodes();
         for (int i = 0; i < n; i++)
             assert (output[i] == ordering[i]);
        cout << "Parmetis Test "<< ++testcount << " : Passed" << endl;
    }
}

void parmetis_with_metmmd ()
{
    DEBUG("########## Testing Function : %s #################\n", __FUNCTION__);
    int size, rank;
    MPI_Comm comm;
	Graph::GraphEOUtil eoutil;
    Graph::VertexWeightedGraph *wmg;
	Graph::GraphCreatorFile creator;

	int output[] =
        {102,58,128,121,86,64,96,112,103,117,120,87,124,126,127,63,95,108,110,111,91,107,109,119,123,125,115,59,61,62,90,93,122,114,106,118,32,116,55,100,16,79,84,78,24,48,52,56,44,54,60,80,28,30,31,40,46,47,72,76,88,92,94,104,8,101,1,42,4,25,97,81,6,7,3,2,73,49,13,5,9,11,17,33,65,66,67,69,10,18,19,21,34,35,37,41,113,15,68,105,77,50,14,43,20,57,12,70,53,51,23,22,36,38,82,83,98,99,74,26,27,29,39,45,71,75,85,89};

    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    creator.set_file_name("../data/1dc.128.txt");
    creator.set_graph_type("DIMACS");
    wmg = creator.create_weighted_mutable_graph();
    
	vector<int> ordering(wmg->get_num_nodes(), -1);
    eoutil.parmetis_elimination_ordering(wmg, ordering, GD_METIS_MMD, false, comm);
    if (rank == 0)
    {
        int n = wmg->get_num_nodes();
        for (int i = 0; i < n; i++)
             assert (output[i] == ordering[i]);
        cout << "Parmetis Test "<< ++testcount << " : Passed" << endl;
    }
}
#endif
//for (int i = 0; i < n; i++)
//    GEN("%d,", ordering[i]);
//GEN("\n");


int main(int argc, char** argv)
{
#ifdef HAS_PARMETIS
    int size, rank, rc;
	char lfname[100];
	char efname[100];
    testcount = 0;
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        FERROR("MPI Initialization error\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rank = 0;    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	sprintf(lfname, "parmetislog.%04d", rank);
	//0: debug, 5: critical
	LOG_INIT(lfname, NULL, 0);

    parmetis_with_amd();
    parmetis_only();
    parmetis_with_metmmd();

    MPI_Finalize();
    LOG_CLOSE();
#else
    cout << "Enable HAS_PARMETIS flag and recompile to test parmetis functions" << endl;
#endif // HAS_PARMETIS
	return 0;
}









