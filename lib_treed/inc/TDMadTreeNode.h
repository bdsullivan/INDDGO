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

#ifndef TDMADTREENODE_H_
#define TDMADTREENODE_H_
#ifdef __MADNESS__
#include "TDTreeNode.h"
#include "TDSolution.h"
#include <world/archive.h>
#include <world/worldmutex.h>

class TDMadTreeNode: public TDTreeNode
{

  public:
	TDMadTreeNode();
	TDMadTreeNode(const TDMadTreeNode& rhs);
	TDMadTreeNode& operator=(const TDMadTreeNode& rhs);
	virtual ~TDMadTreeNode();
    
};

/*
 * Asymetric store and load implementation using madness
 */

namespace madness
{
namespace archive
{

template<class Archive>
struct ArchiveStoreImpl<Archive, TDMadTreeNode>
{
	static void store(const Archive& ar, const TDMadTreeNode& node)
	{
		int i = 0;
		int adjsize = 0;
		int w = node.bag_vec.size();
		int num_ch = 0;
		int num_mask_words;
		bool sflag = false;
		int completed_children = 0;
		int count = 0;
		int masklen = 0;

		std::list<int>::const_iterator it;

		ar & node.id;

		num_mask_words = node.get_num_mask_words();
		ar & num_mask_words;

		completed_children = node.get_completed_children();
		ar & completed_children;

		// store bag vec
		ar & node.bag_vec;

		// store adj list
		i = 0;
		adjsize = node.adj.size();
		ar & adjsize;
		it = node.adj.begin();
		while (i < adjsize)
		{
			ar & (*it);
			++it;
			i++;
		}

		//store parent position vector
		if (node.parent_position)
		{
			sflag = true;
			ar & sflag;
			ar & *node.parent_position;
		}
		else
		{
			sflag = false;
			ar & sflag;
		}

		if (node.nbr_mask_vec)
		{
			sflag = true;
			ar & sflag;
			//store nbr_mask_vec a.k.a adj matrix
			ar & node.nbr_mask_vec->at(0)->w->S;
			for (i = 0; i < w; i++)
				ar & *(node.nbr_mask_vec->at(i));
		}
		else
		{
			sflag = false;
			ar & sflag;
		}

		if (node.child_intersection)
		{
			sflag = true;
			ar & sflag;
			//store children intersection
			num_ch = node.adj.size() - 1;
			if (num_ch > 0)
			{
				ar & node.child_intersection->at(0)->S;
				for (i = 0; i < num_ch; i++)
					ar & *(node.child_intersection->at(i));
			}
		}
		else
		{
			sflag = false;
			ar & sflag;
		}

		//store vertices weights vector
		if (node.vertex_weights)
		{
			sflag = true;
			ar & sflag;
			ar & *node.vertex_weights;
		}
		else
		{
			sflag = false;
			ar & sflag;
		}

		//serialize solution hashtable
		if (node.hash_table)
		{
			count = HASH_COUNT(node.hash_table);
			ar & count;
			if (count > 0)
			{
				masklen = node.hash_table->get_mask_length();
				ar & masklen;
				TDSolution *added_set;
				TDSolution *tmp_sol;
				HASH_ITER(hh, node.hash_table, added_set, tmp_sol)
				{
					ar & *(added_set->mask);
					ar & added_set->value;
				}
			}
		}
		else
		{
			ar & count;
		}
	}
	;
};

template<class Archive>
struct ArchiveLoadImpl<Archive, TDMadTreeNode>
{
	static void load(const Archive& ar, TDMadTreeNode& node)
	{
		int value;
		int adjsize;
		int i = 0;
		int w = 0;
		int ws = 0;
		int num_ch = 0;
		int num_mask_words = 8;
		bool flag = false;
		int completed_children = 0;
		int count = 0;
		int masklen = 0;

		ar & node.id;

		ar & num_mask_words;
		node.set_num_mask_words(num_mask_words);

		ar & completed_children;
		node.set_completed_children(completed_children);

		// load bag vec
		ar & node.bag_vec;

		// restoring bag list
		w = node.bag_vec.size();
		for (i = 0; i < w; i++)
			node.bag.push_back(node.bag_vec[i]);

		// load adj list
		i = 0;
		ar & adjsize;
		while (i < adjsize)
		{
			ar & value;
			node.adj.push_back(value);
			i++;
		}

		// load parent position vector
		flag = false;
		ar & flag;
		if (flag)
		{
			node.parent_position = new vector<int> ();
			ar & *node.parent_position;
		}

		flag = false;
		ar & flag;

		if (flag)
		{
			//load nbr_mask_vec
			ar & ws;
			w = node.bag_vec.size();
			node.nbr_mask_vec = new vector<int_bigint *> (w);
			for (i = 0; i < w; i++)
			{
				node.nbr_mask_vec->at(i) = new int_bigint(ws);
				ar & *(node.nbr_mask_vec->at(i));
			}
		}

		flag = false;
		ar & flag;

		if (flag)
		{
			//load children intersection
			num_ch = node.adj.size() - 1;
			if (num_ch > 0)
			{
				ar & ws;

				node.child_intersection = new vector<bigint_t *> (num_ch);
				for (i = 0; i < num_ch; i++)
				{
					node.child_intersection->at(i) = new bigint_t(ws);
					ar & *(node.child_intersection->at(i));
				}
			}
		}

		// load vertices weights vector
		flag = false;
		ar & flag;
		if (flag)
		{
			node.vertex_weights = new vector<int> ();
			ar & *node.vertex_weights;
		}

		ar & count;
		if (count > 0)
		{
			TDSolution *added_set;
			ar & masklen;
			i = 0;
			while (i < count)
			{
				added_set = new TDSolution(masklen);
				ar & *(added_set->mask);
				ar & added_set->value;
				HASH_ADD_KEYPTR(hh, node.hash_table, added_set->mask->words,
						masklen*sizeof(BIGINT_WORD), added_set);
				i++;
			}
		}

	}
	;
};
}
}
#endif /* __MADNESS__ */

#endif /* TDMADTREENODE_H_ */
