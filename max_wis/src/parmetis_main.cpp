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
#include "ORD.h"

static void *write_file(void *a)
{
    Graph::GraphWriter *writer;
    Graph::GraphReaderWriterFactory factory;
    Graph::WeightedMutableGraph *g;

    ORD *ord = reinterpret_cast<ORD*>(a);
    if (!ord->is_original())
    {
        writer = factory.create_writer("DIMACS");
        const char *compname = ord->get_filename();
        g = ord->get_graph();
        writer->set_out_file_name(compname);
        writer->write_graph(g);
        DEBUG("wrote into filesystem: %s\n", compname);
    }
}

int main(int argc, char** argv)
{
    int size, rank, rc;
	char lfname[100];
	char efname[100];
    MPI_Comm comm;

    rank = 0;    
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
        FERROR("MPI Initialization error\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	sprintf(lfname, "madlog.%04d", rank);
	//0: debug, 5: critical
	LOG_INIT(lfname, NULL, 0);

    if (argc < 3)
    {
        fprintf (stderr, "usage: %s -f file\n", argv[0]);
        exit (-1);
    }

    if (strcmp(argv[1], "-f"))
    {
        fprintf (stderr, "usage: %s -f file\n", argv[0]);
        exit (-1);
    }

    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    // creating ORD class generate elimination ordering using PARMetis
    ORD *p = new ORD(comm, argv[2]);
    
    DEBUG("building graph\n");
    p->buildGraph();

    pthread_t pid;
    if (rank == 0)
    {
        int rc = pthread_create(&pid, 0, write_file, p);
        if (rc)
        {
            FERROR("error occured while creating pthread: %d\n", rc);
        }
    }

    rc = MPI_Barrier(MPI_COMM_WORLD);
    if (rc != MPI_SUCCESS)
    {
        FERROR("error with MPI Barrier\n");
    }

    DEBUG("find elimination ordering\n");
    p->find_elim_ordering();

    if (rank == 0)
    {
        rc = pthread_join (pid, NULL);
        if (rc)
        {
            FERROR("error occured while joining pthreads\n");
        }
    }

    rc = MPI_Comm_free(&comm);
    if (rc != MPI_SUCCESS)
    {
        FERROR("error with MPI_Comm_free\n");
    }

    rc = MPI_Finalize();
    if (rc != MPI_SUCCESS)
    {
        FERROR("error with MPI_Finalize\n");
    }

	LOG_CLOSE();
	return 0;
}






