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
#include "MadnessTableProcessor.h"
#include "UTHashWrapper.h"
#include "GraphException.h"
#ifdef VTRACE
#include "vt_user.h"
#endif


MadnessTableProcessor::MadnessTableProcessor(World& world, int num_nodes, int n, int seed) :
	count(0), t(cpu_time()), num_mask_words(n), WorldObject<MadnessTableProcessor> (world) 
{
    // For random assignment:
    std::shared_ptr< WorldDCPmapInterface<TDKey> > pmap(new TDPmap(world, seed));
	node_container = new WorldContainer<TDKey, TDMadTreeNode> (world, pmap);
    // For static assignment:
    //node_container = new WorldContainer<TDKey, TDMadTreeNode> (world, pmap);
    vec_lock = (MutexFair **) malloc(sizeof(MutexFair *) * num_nodes);
    vec_hashtbl = (TDSolution **)malloc(sizeof(TDSolution *) * num_nodes);
    vec_ind_sets = (list<pair<bigint_t, int>> **) malloc (sizeof(list<pair<bigint_t,int>> *)*num_nodes);

	process_pending();
}

MadnessTableProcessor::~MadnessTableProcessor()
{
	delete node_container;
    free (vec_lock);
    free (vec_hashtbl);
    free (vec_ind_sets);
}

void MadnessTableProcessor::add_node(int i, TDMadTreeNode& node) const
{
	TDKey key(i);
	node_container->replace(key, node);
}

TDMadTreeNode& MadnessTableProcessor::get_node(int i)
{
	TDKey key(i);
	return node_container->find(key).get()->second;

}

Void MadnessTableProcessor::delete_node(int i)
{

    ScopedMutex<MutexFair>m(lock);
	TDKey key(i);
	node_container->erase(key);
    return None;
}


ProcessID MadnessTableProcessor::owner(int i) const
{
	TDKey key(i);
	return node_container->owner(key);
}

int MadnessTableProcessor::find_max(TDMadTreeNode node)
{
	TDSolution *added_set;
	TDSolution *tmp_sol;
	int max = -1;
	HASH_ITER(hh, node.hash_table, added_set, tmp_sol)
	{
		if (max < added_set->value)
			max = added_set->value;
	}

	return max;
}

Void MadnessTableProcessor::delete_table(int i)
{
    //if (hashtbl_ptr.count(i) > 0)
    if (vec_hashtbl[i])
    {
        //UTHashWrapper uth(hashtbl_ptr[i], num_mask_words);
        UTHashWrapper uth(vec_hashtbl[i], num_mask_words);
        uth.hash_delete();
        //free(vec_hashtbl[i]);
    }
    return None;
}


TDMadTreeNode MadnessTableProcessor::compute_table_leaf(int inode, int p)
{
#ifdef VTRACE
    VT_TRACER("comp_leaf");
#endif 

	TDMadTreeNode node = get_node(inode);

    int i, j, w = (int) node.bag.size();
    list<int>::iterator ii, jj;
    vector<int>::iterator uu, vv;
    bigint_t current_mask(num_mask_words);
    list<bigint_t*>::iterator current_mask_it;
    list<TDSolution*>::iterator ss, tt;

    // We already have all the necessary information inside the TDTreeNode
    current_mask = 0;
    // Find out when to stop the loop
    bigint_t two_pow_w(num_mask_words), one(num_mask_words);
    one = BIGINT_ONE;
    two_pow_w.two_x(w);
    int two_w_word = 1;
    while (w >= BIGINT_WORD_SIZE * two_w_word)
        two_w_word++;
    two_w_word--;
    bool is_independent, has_bit;
    int num_ind_sets = 0;
    TDSolution temp_set(num_mask_words);

    UTHashWrapper uth(node.hash_table, num_mask_words);
    while (current_mask.words[two_w_word] < two_pow_w.words[two_w_word])
    {
        is_independent = true;
        for (i = 0; i < w; i++)
        {
            if (current_mask.test_bit(node.nbr_mask_vec->at(i)->k))
            {
                // current_mask has the bit corresponding to nbr_mask_vec[i] set - so and with the appropriate vector for a hit
                has_bit = false;
                // Do the loop below manually rather than calling the overloaded AND
                for (j = 0; j < current_mask.S; j++)
                {
                    if (current_mask.words[j]
                        & node.nbr_mask_vec->at(i)->w->words[j])
                    {
                        // Some nbr is found in word j so current_mask cannot represent an ind. set
                        has_bit = true;
                        break;
                    }
                }
                // breaks to here
                if (has_bit)
                {
                    // Set represented by current_mask is not independent
                    is_independent = false;
                    int rpos = 0;
                    // Find the rightmost bit in current_mask - this could be sped up by looking at the words
                    // SPEEDUP
                    while (!(current_mask.test_bit(rpos)))
                        rpos++;

                    // Now advance current_mask by a potentially massive amount
                    // Would be faster here to figure out which word we are in and
                    // do this w/o test_bit
                    int b = 0;
                    while (current_mask.test_bit(rpos + b))
                        b++;
                    for (int a = 0; a < b; a++)
                        current_mask.xor_bit(rpos + a);
                    current_mask.or_bit(rpos + b);
                    break;
                }
            }
        }

        if (is_independent)
        {
            // increment total_table_entries
            num_ind_sets++;
            bool is_best = false;
            int new_value = 0, parent_weight = 0;
            for (j = 0; j < w; j++)
            {
                if (current_mask.test_bit(j))
                {
                    // Compute the actual weight in G of this ind. set
                    new_value += node.vertex_weights->at(j);
                }
            }

            // Now look up to the parent to translate the masks
            temp_set.mask->zeroize();
            parent_weight = 0;
            for (j = 0; j < w; j++)
            {
                if (current_mask.test_bit(j))
                {
                    if (node.parent_position->at(j) != -1)
                    {
                        temp_set.mask->set_bit(node.parent_position->at(j));
                        parent_weight += node.vertex_weights->at(j);

                    }
                }
            }
            // temp_mask is a mask of this set \cap parent_bag
            temp_set.value = new_value;

            TDSolution *hash_ptr;

            if (node.is_root())
            {
                if (is_best)
                    uth.hash_add(*(temp_set.mask), temp_set.value);
            }
            else
            {
                hash_ptr = uth.hash_find(*temp_set.mask);
                if (hash_ptr)
                {
                    if (temp_set.value - parent_weight > hash_ptr->value)
                    {
                        // We found an entry for this mask and this mask has a better value
                        // so update the entry in the table
                        hash_ptr->value = temp_set.value - parent_weight;
                        *(hash_ptr->orig_mask) = current_mask;
                    }
                }
                else
                {
                    // There is no entry for this mask - add it to the
                    // table
                    uth.hash_add(*(temp_set.mask), (temp_set.value - parent_weight));
                }
            }

            // Advance to next candidate ind. set
            ++current_mask;
        }

    }

    vec_hashtbl[node.id] = node.hash_table;
    return node;
}

int MadnessTableProcessor::create_new_table(int inode)
{
#ifdef VTRACE
    VT_TRACER("create_new_tbl");
#endif 

    // Following scoped mutex will make sure only one thread will work
    // on critical section .
    ScopedMutex<MutexFair> mu(vec_lock[inode]);

    //if (ind_sets.count(inode) > 0)
    if (vec_ind_sets[inode] != NULL)
    {
        return 1;
    }

    TDMadTreeNode node = get_node(inode);

    bigint_t tmp_set(num_mask_words);
    int i, j;
    int w = (int) node.bag.size();

    bigint_t current_mask(num_mask_words);
    long num_ind_sets;

    list<pair<bigint_t, int>> *mvs = new list<pair<bigint_t, int>>();
    // We already have all the necessary information inside the TDTreeNode
    current_mask = 0;
    // Find out when to stop the loop
    bigint_t two_pow_w(num_mask_words), one(num_mask_words);
    one = BIGINT_ONE;
    two_pow_w.two_x(w);
    int two_w_word = 1;
    while (w >= BIGINT_WORD_SIZE * two_w_word)
        two_w_word++;
    two_w_word--;
    bool is_independent, has_bit;

    TDSolution temp_set(num_mask_words);

    while (current_mask.words[two_w_word] < two_pow_w.words[two_w_word])
    {
        is_independent = true;
        for (i = 0; i < w; i++)
        {
            if (current_mask.test_bit(node.nbr_mask_vec->at(i)->k))
            {
                // current_mask has the bit corresponding to nbr_mask_vec[i] set - so and with the appropriate vector for a hit
                has_bit = false;
                // Do the loop below manually rather than calling the overloaded AND
                for (j = 0; j < current_mask.S; j++)
                {
                    if (current_mask.words[j]
                        & node.nbr_mask_vec->at(i)->w->words[j])
                    {
                        // Some nbr is found in word j so current_mask cannot represent an ind. set
                        has_bit = true;
                        break;
                    }
                }
                // breaks to here
                if (has_bit)
                {
                    // Set represented by current_mask is not independent
                    is_independent = false;
                    int rpos = 0;
                    // Find the rightmost bit in current_mask - this could be sped up by looking at the words
                    // SPEEDUP
                    while (!(current_mask.test_bit(rpos)))
                        rpos++;

                    // Now advance current_mask by a potentially massive amount
                    // Would be faster here to figure out which word we are in and
                    // do this w/o test_bit
                    int b = 0;
                    while (current_mask.test_bit(rpos + b))
                        b++;
                    for (int a = 0; a < b; a++)
                        current_mask.xor_bit(rpos + a);
                    current_mask.or_bit(rpos + b);
                    break;
                }
            }
        }

        if (is_independent)
        {
            // increment total_table_entries
            num_ind_sets++;
            int new_value = 0;
            for (j = 0; j < w; j++)
            {
                if (current_mask.test_bit(j))
                {
                    // Compute the actual weight in G of this ind. set
                    new_value += node.vertex_weights->at(j);
                }
            }

            mvs->push_back(make_pair(current_mask, new_value));
            // Done with child DP - advance to look for the next independent set in the k's bag
            ++current_mask;
        }
    }

    vec_ind_sets[inode] = mvs;
    return 1;
}

/*
 * If all independent set table is exists already, just update values.
 */
int MadnessTableProcessor::update_table(int inode, TDMadTreeNode cnode,
                                        int cindex, int child)
{
#ifdef VTRACE
    VT_TRACER("update_tbl");
#endif 

    // For debugging, uncomment following block to see the progress
    //if (world.rank() == 0)
    //{
    //    lock.lock();
    //    count ++;
    //    if ((count % 1000) == 0)
    //        print(cpu_time() - t, ":", count, ":", inode, ":", child, ":", cindex);
    //    lock.unlock();
    //}

    // Following scoped mutex will make sure only one thread will work
    // on critical section .
    ScopedMutex<MutexFair> mu(vec_lock[inode]);

	TDMadTreeNode node = get_node(inode);

    if (!node.is_root())
        task(owner(node.adj.front()), &MadnessTableProcessor::create_new_table, node.adj.front());    

    bigint_t tmp_set(num_mask_words);
    int w = (int) node.bag.size();
    bigint_t current_mask(num_mask_words);
    int inter_weight = 0;

    list<bigint_t *>::iterator bi_iter;
    list<int>::iterator i_iter;
    bigint_t *child_intersection = node.child_intersection->at(cindex);


    //while(!(ind_sets.count(inode) > 0)) usleep(100);

    list<pair<bigint_t, int> >::iterator itl;
    list<pair<bigint_t, int>> *mvs = vec_ind_sets[inode];

    if (!mvs)
    {
        print("[error] no independent set table : ", inode, child, cindex, world.rank());
        string dsc("[error] no independent set table : " +  inode);
        throw Graph::GraphException(dsc);
    }

    itl = mvs->begin();

    
    // class wrapper for uthash
    UTHashWrapper uth(cnode.hash_table, num_mask_words);

    for (; itl != mvs->end(); ++itl)
    {
        tmp_set.zeroize();
        tmp_set = itl->first;
        tmp_set &= *child_intersection;

        //calculating weights for the nodes in the intersection
        inter_weight = 0;

        int i = 0;
        for (i = 0; i < w; i++)
        {
            if (tmp_set.test_bit(i))
            {
                inter_weight += node.vertex_weights->at(i);
            }
        }


        TDSolution *hash_ptr;
        hash_ptr = uth.hash_find(tmp_set);
            
        if (hash_ptr)
        {
            itl->second += hash_ptr->value;
        }
        else
        {
            fatal_error(
                "can not find intersection for current mask, abort !!\n");
            current_mask.print(__error_file__);
        }
    }

    // deletes local hash table
    uth.hash_delete();
    if (owner(child) != world.rank())
    {
        // send a message to delete hashtable in remote machine
        task(owner(child), &MadnessTableProcessor::delete_table, child);
    }

    return 1;
}

/**
 * After processing all children, all independent set table will be reduced to intersection (with parent) table.
 */
TDMadTreeNode MadnessTableProcessor::create_intersection_table(int inode, int p)
{
#ifdef VTRACE
    VT_TRACER("pc_inter_tbl");
#endif 

    // ScopedMutex<MutexFair> mu(nodelocks[inode]);

	TDMadTreeNode node = get_node(inode);

    // All the children have been processed so now go through the complete table
    // and create the hash table for the parent node. Since the masks in the hash table
    // are stored in terms of the parent's bag, this requires you to now look up to the
    // parent's parent using the parent_position information!
    //Save parent child intersection in parent hash and clean completed_table
    int max = 0;
    TDSolution temp_set(num_mask_words);
    bigint_t current_mask(num_mask_words);
    int w = (int) node.bag.size();

    //while(!(ind_sets.count(inode) > 0)) usleep(100);

    list<pair<bigint_t, int> >::iterator itl;
    list<pair<bigint_t, int>> *mvs = vec_ind_sets[inode];
    UTHashWrapper uth(node.hash_table, num_mask_words);

    if (!mvs)
    {
        print("[error] intersection no independent set table : ", inode, world.rank());
        string dsc("[error] no independent set table : " +  inode);
        throw Graph::GraphException(dsc);
    }
    
    itl = mvs->begin();
    while (itl != mvs->end())
    {
        if (node.is_root())
        {
            if (itl->second > max)
            {
                max = itl->second;
                uth.hash_add(itl->first, itl->second);
            }
            
            //itl = mvs->erase(itl);
            ++itl;
            continue;
        }

        int j = 0;
        int parent_weight = 0;

        temp_set.mask->zeroize();
        current_mask = itl->first;

        for (j = 0; j < w; j++)
        {
            if (current_mask.test_bit(j))
            {
                if (node.parent_position->at(j) != -1)
                {
                    temp_set.mask->set_bit(node.parent_position->at(j));
                    parent_weight += node.vertex_weights->at(j);
                }
            }
        }

        temp_set.value = itl->second;

        TDSolution *hash_ptr;
        // k is not the root - so do the usual
        hash_ptr = uth.hash_find(*(temp_set.mask));

        if (hash_ptr)
        {
            if (temp_set.value - parent_weight > hash_ptr->value)
            {
                // We found an entry for this mask and this mask has a better value
                // so update the entry in the table
                hash_ptr->value = temp_set.value - parent_weight;
                *(hash_ptr->orig_mask) = current_mask;
            }
        }
        else
            uth.hash_add(*(temp_set.mask), (temp_set.value - parent_weight));

        //itl = mvs->erase(itl);
        ++itl;
    }

    mvs->clear();
    // Deletes mutex for current node.
    //nodelocks.erase(inode);
    //delete vec_lock[inode];
    //vec_lock[inode] = NULL;
    //hashtbl_ptr[inode] = node.hash_table;
    vec_hashtbl[inode] = node.hash_table;
    task(owner(inode), &MadnessTableProcessor::delete_node, inode);
    return node;
} 

Future<TDMadTreeNode> MadnessTableProcessor::check_updates(
    vector<Future<int>> status, int inode)
{
	int ssize = status.size();
	int i = 0;
    int a = 0;
    Future<int> p;

	TDMadTreeNode node = get_node(inode);

	for (i = 0; i < ssize; i++)
		a = status[i].get();

    p = task(owner(node.adj.front()), &MadnessTableProcessor::create_new_table, node.adj.front());    
	return task(owner(inode),
                &MadnessTableProcessor::create_intersection_table, inode, p);
}

Future<TDMadTreeNode> MadnessTableProcessor::compute_table(int inode)
{
#ifdef VTRACE
    VT_TRACER("comp_tbl");
#endif 

	// process local
	TDMadTreeNode node = get_node(inode);

	if (node.is_leaf())
    {
        Future<int> p;
        p = task(owner(node.adj.front()), &MadnessTableProcessor::create_new_table, node.adj.front());    
		return task(owner(inode), &MadnessTableProcessor::compute_table_leaf, inode, p);
    }
	int i = 0;
    Future<int> a;
	vector<Future<TDMadTreeNode> > ch_sol;
	vector<Future<int>> status_vec;
	
    list<int>::iterator it = node.adj.begin();
	++it; // first entry is the parent

    // create a mutex for current node, this mutex will make sure
    // threads will update all independent set table without any
    // conflicts. 
    vec_lock[inode] = new MutexFair();
    vec_hashtbl[inode] = NULL;
    vec_ind_sets[inode] = NULL;
    
    // tasks for compute_table on children.
	for (i = 0; it != node.adj.end(); ++it, i++)
		ch_sol.push_back(
            task(owner(*it), &MadnessTableProcessor::compute_table, *it));
    
    // launch tasks which updates all independent sets table as
    // children finish.
    it = node.adj.begin();
    ++it;
    for (i = 0; it != node.adj.end(); ++it, i++)
		status_vec.push_back(
            task(owner(inode), &MadnessTableProcessor::update_table, inode,
                 ch_sol[i], i, *it));

    // After all updates are finished, create intersection table and
    // return TDMadTreeNode to parent.
	return task(owner(inode), &MadnessTableProcessor::check_updates,
                status_vec, inode);
}

#endif /* __MADNESS__ */
