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

#ifndef PTD_H_
#define ORD_H_

#include "GraphDecomposition.h"
#include "Log.h"
#include "mpi.h"

class ORD
{

private:
	vector<int> ordering;
	Graph::WeightedMutableGraph *G;
	int width;
    MPI_Comm comm;
    string filename;
    bool cflag;

public:
	ORD(MPI_Comm comm, char *filename);
	ORD(char *filename);
	virtual ~ORD();
    void buildGraph();
    void set_filename(char *filename);
    const char *get_filename() const;
    int find_elim_ordering();
    void outputMetisFormat();
    Graph::WeightedMutableGraph *get_graph();
    bool is_original();
    void set_original(bool x);

};
#endif /* ORD_H_ */
