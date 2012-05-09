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

#ifndef MADNESSTABLEPROCESSOR_H_
#define MADNESSTABLEPROCESSOR_H_
#ifdef __MADNESS__
#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>
#include <world/worldobj.h>
#include <world/worlddc.h>
#include <iomanip>

#include "GraphDecomposition.h"
#include "TreeDecomposition.h"
#include "Log.h"
using namespace std;
using namespace madness;

struct TDKey
{
	int k;

	TDKey() :
		k(-1)
        {
        }

	TDKey(int k) :
		k(k)
        {
        }

	hashT hash() const
        {
            return k;
        }

	int get_val() const
        {
            return k;
        }

	template<typename Archive>
	void serialize(const Archive& ar)
        {
            ar & k;
        }

	bool operator==(const TDKey& b) const
        {
            return k == b.k;
        }
};

class TDPmap : public WorldDCPmapInterface<TDKey> {
  private:
    const int nproc;
    const int seed;
    vector<int> node_map;
  public:
    TDPmap(World& world, int seed)
        : nproc(world.mpi.nproc()), seed(seed)
    { 
        int i = 0;
        int j = 0;
        srand(seed);
        for (i = 0; i < nproc; i++)
            node_map.push_back(i);
        random_shuffle(node_map.begin(), node_map.end());
        random_shuffle(node_map.begin(), node_map.end());
    };

        ProcessID owner(const TDKey& key) const {
            int k = key.get_val();
            if (nproc == 1) return 0;
            return node_map[(k + k/nproc)%nproc];
        };
};


class MadnessTableProcessor: public WorldObject<MadnessTableProcessor>
{
  private:
	WorldContainer<TDKey, TDMadTreeNode> *node_container;

    // d8w:
    // replaced tr1::unordered_map with a map to fix an error with 2
    // million node graph. 
    //
    // Note: I think, with low width graph DP finished so fast that
    // unordered map can not keep up with the speed and return no
    // table error. Frequent addition and deletion of entries may
    // cause the slowness in unordered_map. On the other hand map can 
    // keep up with the speed since each key has only one entry.

    // d8w: 04/20/2012
    // Had occasional segfaults with larger graphs since map
    // operations are not thread safe therefore replaced maps with arrays.
	//map<int, list<pair<bigint_t, int>>*> ind_sets;
    //map<int, MutexFair *> nodelocks;
    //map<int, TDSolution *> hashtbl_ptr;

    MutexFair **vec_lock;
    TDSolution **vec_hashtbl;
    list<pair<bigint_t,int>> **vec_ind_sets;

    MutexFair lock;
	int num_mask_words;
    int count;
    double t;

	TDMadTreeNode compute_table_leaf(int inode, int p);
	int update_table(int inode, TDMadTreeNode cnode, int cindex, int child);
	int create_new_table(int inode);
	Future<TDMadTreeNode> check_updates(vector<Future<int>> status, int inode);
	TDMadTreeNode create_intersection_table(int inode, int p);
	vector<bigint_t *> generate_masks(int w, int nchunks);
	vector<pair<bigint_t, int>> search_region(int inode, bigint_t& start, bigint_t &end);
	vector<pair<bigint_t, int>> get_final_list(vector<Future<vector<pair<bigint_t, int>>> > futlist);

  public:
    MadnessTableProcessor(World& world, int num_nodes, int n, int seed);
    virtual ~MadnessTableProcessor();
    void add_node(int i, TDMadTreeNode& node) const;
    TDMadTreeNode& get_node(int i);
    Void delete_node(int i);
    ProcessID owner(int i) const;
    Future<TDMadTreeNode> compute_table(int i);
    int find_max(TDMadTreeNode node);
    Void delete_table(int i);
};

#endif /* MADNESSTABLEPROCESSOR_H_ */
#endif /* __MADNESS__ */
