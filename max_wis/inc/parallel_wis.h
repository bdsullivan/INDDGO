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

#ifndef __PARALLEL_WIS__
#define __PARALLEL_WIS__


enum MPI_TAGS {
    MPI_SEND_RECV_TAG = 100000,
    MPI_REQ_WRK_TAG,
    MPI_WRK_TAG,
    MPI_REQ_WGH_TAG,
    MPI_WGH_TAG,
    MPI_REQ_COUNT_TAG,
    MPI_TBL_LOC_TAG,
    MPI_TBL_SIZE_TAG,
    MPI_COMP_TERM_TAG,
    MPI_STOR_TERM_TAG,
    MPI_TBL_SAVE_TAG,
    MPI_CHILD_FREE_TAG,
};

enum NODE_TYPES {
    HEAD_NODE = 100,
    COMPUTE_NODE,
    STORAGE_NODE,
};

#define REQUEST_SIZE 4096
#define THRESHOLD 60
#define CHUNKSIZE 16
#define MINENTRIES 512
#define LOG2SIZE(x) (int)(floor(log((double)x)/log((double)2)))

int compute_weighted_ind_sets_parallel (TDTree *T, int k);
int compute_weighted_ind_sets_test (TDTree *T, int k);
int compute_nonnice_table_test (TDTree *T, int k);
void parallel_wis_init (TDTree *T, int size, int rank);
void parallel_wis_cleanup (TDTree *T);
void *thread_recv_msk_send_wgh(void *v);
void parallel_wis_init_non_thread (TDTree *T, int size, int rank);
void parallel_wis_cleanup_non_thread (TDTree *T);
void parallel_wis_init_9 (TDTree *T, int size, int rank);
void parallel_wis_cleanup_9 (TDTree *T);

void parallel_wis_head_init (TDTree *T, int size, int rank);
void parallel_wis_compute_init (TDTree *T, int size, int rank);
void parallel_wis_storage_init (TDTree *T, int size, int rank);

void parallel_wis_head_cleanup (TDTree *T);
void parallel_wis_compute_cleanup (TDTree *T);
void parallel_wis_storage_cleanup (TDTree *T);

void parallel_wis_head_start (TDTree *T);
void parallel_wis_storage_start (TDTree *T);


#endif  /* __PARALLEL_WIS__ */
