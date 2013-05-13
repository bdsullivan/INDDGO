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
#define PTD_H_
#ifdef __MADNESS__
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>
#include <world/worldobj.h>
#include <world/worlddc.h>
#include <iomanip>

#include "GraphDecomposition.h"
#include "TreeDecomposition.h"
#include "Log.h"
#include  "mpi.h"
#include "pthread.h"
using namespace std;
using namespace madness;


using namespace std;

class PTD: public WorldObject<PTD>
{
  private:

	vector<int> ordering;
    vector<int> position;
    Graph::VertexWeightedGraph G;
    Graph::VertexWeightedGraph H;
	int width;
    double xtime;
    Graph::GraphEOUtil eoutil;
    string tmp_filename;
    string tree_name;
    int s;
    int r;

  public:
    MutexFair tlock;
    ofstream fw;

    PTD(World& world, Graph::VertexWeightedGraph g);
	virtual ~PTD();

    int owner(int i);
	int get_width() const;
	void set_width(int w);
	int triangulate();
    int size();
    int rank();

    int compute_elim_order(int order_type, bool t, const char *s);
    int parmetis_elim_order(int order_type);
	vector<int> *get_elim_order();
    void set_elim_order(vector<int> &order);
    int write_elim_order(char *filename);

    int create_bag(int v, int i, list<string> &ls);
    void process_bag();
    string get_random_name(char *x);
    void set_file_name(string s);
    string get_file_name() const;
    void set_tree_name (string t);
    string get_tree_name() const;

    void merge_files(char *x);
    void refine_td(string x);
    
    int get_int(char *x);

    void create_tree_decomposition(char *filename, TDTree *&T);
};

// TPTD class is used to pass data into threaded environment.
class TPTD {
  private:
    PTD *p;
    int id;
    vector<int> range_vector;

  public:
    
    TPTD(PTD *p, vector<int> range, int id): 
        p(p), range_vector(range), id(id)
    {};
        
        virtual ~TPTD() 
        {};
        
        int get_id() 
        {
            return id;
        };

        vector<int> get_range_vector() 
        {
            return range_vector;
        };

        PTD* get_ptd() 
        {
            return p;
        };
};
#endif
#endif /* PTD_H_ */

















