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

#include "GraphDecomposition.h"
#include "TreeDecomposition.h"
#include <iostream>
#include <cmath>



#include "Log.h"
#include <algorithm>
#include <numeric>

#include <weighted_ind_set.h>
#include <cmath>

#define WRAP_UTHASH 0

#ifdef __PARALLEL__
#include <omp.h>
#include <mpi.h>
#include <parallel_wis.h>

static unsigned long long *encode_weights(vector<double> *, int *, int *);
static vector<int> decode_weights(unsigned long long *, double *, int *);
static int parallel_mask_gen (TDTree *T, int k, int bagsize, vector<bigint_t *> *start_mask_vec, vector<bigint_t *> *end_mask_vec, int nchunks);
static int parallel_mask_gen_1 (TDTree *T, int k, int bagsize, vector<bigint_t *> *start_mask_vec, vector<bigint_t *> *end_mask_vec, int nchunks);
static void mask_split (TDTree *T, int k,  bigint_t *start_mask, bigint_t *end_mask, bigint_t *resultsvec, int nchunks);
static void populate_count_disp (TDTree *T, int total_size, vector<int> *count, vector<int> *disp);


static int compute_nonnice_table_static_lb (TDTree *T, int k);
static int compute_nonnice_table_parallel_0 (TDTree *T, int k);
static int compute_nonnice_table_parallel_1 (TDTree *T, int k);
static int compute_nonnice_table_parallel_2 (TDTree *T, int k);
static int compute_nonnice_table_parallel_3 (TDTree *T, int k);
static int compute_nonnice_table_parallel_4 (TDTree *T, int k);
static int compute_nonnice_table_parallel_5 (TDTree *T, int k);
static int compute_nonnice_table_parallel_6 (TDTree *T, int k);
static int compute_nonnice_table_parallel_7 (TDTree *T, int k);
static int compute_nonnice_table_parallel_8 (TDTree *T, int k);
static int compute_nonnice_table_parallel_9 (TDTree *T, int k);
static int compute_nonnice_table_parallel_10 (TDTree *T, int k);

static int list_ind_sets_parallel (TDTree *T, int k, list<int> *bag_subset, list<TDSolution *> *ind_sets, bigint_t *start_mask, bigint_t *end_mask);
static void generate_ind_sets (TDTree *T, int k, list<int> *bag_subset, vector<TDSolution *> *ind_sets);
static int get_msk_send_buffer_index(TDTree *T);
static int get_wgh_send_buffer_index(TDTree *T);
static int send_masks (TDTree *T, int k, int startpos, int iteration, 
                       vector<int> *k_bag_vec, vector<bool *> *cvec, vector<int> *children, vector<bigint_t> *indset_masks);

static int send_masks_9 (TDTree *T, int k, int startpos, int iteration, 
                       vector<int> *k_bag_vec, vector<bool *> *cvec, vector<int> *children);

static void recv_msk_send_wgh(TDTree *T);
static bool recv_msk_send_wgh_9(TDTree *T);

static int recv_wgh(TDTree *T, vector<int> *children, vector<int> *ind_wgh);
static bool recv_wgh_9(TDTree *T, int k, vector<int> *children);


static void mask_split (TDTree *T, int k,  bigint_t *start_mask, bigint_t *end_mask,  bigint_t *resultsvec, int nchunks)
{
    int w;
    int i;
    int msize;
    int log2_block_size;

    resultsvec->zeroize();

    for (i = 0; i < (int)T->tree_nodes[k]->bag.size(); i++)
    {
        if (end_mask->test_bit(i))
            w = i;
    }
    
    msize = LOG2SIZE(nchunks);
    log2_block_size = w - msize;
    msize --;
    
    for (i = 0; i < msize; i++)
    {
        if (msize & (1 << i))
            resultsvec->set_bit(log2_block_size + i);
    }

    //resultsvec->print(__log_file__);
    //end_mask->print(__log_file__);
}

static inline void generate_ind_sets (TDTree *T, int k, list<int> *bag_subset, vector<TDSolution *> *ind_sets)
{
    DEBUG("In function : %s \n", __FUNCTION__);
    int i, n, w;
    vector<bigint_t *> start_masks (CHUNKSIZE);
    vector<bigint_t *> end_masks (CHUNKSIZE);
    vector< list<TDSolution  *> > generated_sets(CHUNKSIZE);
    vector<TDSolution *> ind_sets_vec;

    w = (int) bag_subset->size();
    parallel_mask_gen (T, k, w, &start_masks, &end_masks, CHUNKSIZE);
/*
  for (i = 0; i < CHUNKSIZE; i++)
  {
  DEBUG("start mask 1 : ");
  start_masks_1[i]->print(__log_file__);
  DEBUG("end mask 1 : ");
  end_masks_1[i]->print (__log_file__);
  DEBUG("############\n");
  }
*/


#pragma omp parallel for
    for (i = 0; i < CHUNKSIZE; i++)
    {
        if (i < 1)
        {
            DEBUG("number of threads %d \n", omp_get_num_threads());
        }

        list_ind_sets_parallel (T, k, bag_subset, &(generated_sets[i]), start_masks[i], end_masks[i]);
        //DEBUG("i : %d ind size %d\n", i, generated_sets[i].size());
    }

    for (i = 0; i < CHUNKSIZE; i++)
    {
        ind_sets_vec.insert (ind_sets_vec.end(), generated_sets[i].begin(), generated_sets[i].end());
    }

    *ind_sets = ind_sets_vec;
    DEBUG("return from %s \n", __FUNCTION__);
}

static void populate_count_disp (TDTree *T, int total_size, vector<int> *count, vector<int> *disp)
{
    int modval;
    int blocksize;
    int i;
    vector<int> icount(T->size);
    vector<int> idisp(T->size, 0);

    modval = total_size % T->size;
    blocksize = (total_size - modval)/T->size;
    idisp[0] = 0;
    for (i = 0; i < T->size; i++)
    {
        icount[i] = blocksize;
        if (i < modval)
            icount[i]++;

        if (i > 0)
        {
            idisp[i] = idisp[i - 1] + icount[i - 1];
        }
    }

    *count = icount;
    *disp = idisp;
}

inline static bool mask_compare (bigint_t *a, bigint_t *b)
{
    int i = 0;
    for(i = 0; i < a->S; i++)
    {
        if (a->words[i] & b->words[i])
            return true;
    }

    return false;
}

static void print_bag (TDTree *T, int k)
{
    GEN("bag : ");
    list<int>::iterator itv;
    for (itv = T->tree_nodes[k]->bag.begin(); itv != T->tree_nodes[k]->bag.end(); ++itv)
    {
        GEN("%d ", *itv);
    }
    GEN("\n");
}



static int parallel_mask_gen (TDTree *T, int k, int bagsize, 
                              vector<bigint_t *> *start_mask_vec, vector<bigint_t *> *end_mask_vec, int nchunks)
{
    int w;
    int i;
    int j;
    int msize;
    int log2_block_size;
    int minsize;
    int ln_nchunks;
    bigint_t mask(T->num_mask_words);

    w = bagsize;
    for (i = 0; i < nchunks; i++)
    {
        (*start_mask_vec)[i]->zeroize();
        (*end_mask_vec)[i]->zeroize();
    }

    minsize = w - LOG2SIZE(MINENTRIES);
    ln_nchunks = LOG2SIZE(nchunks);

    if (minsize < 0)
    {
        (*end_mask_vec)[0]->two_x(w);
        return nchunks;
    }

    if (ln_nchunks > minsize)
    {
        nchunks = 0;
        nchunks = 1 << minsize;
        msize = minsize;
        log2_block_size = w - msize;
    }
    else
    {
        msize = LOG2SIZE(nchunks);
        log2_block_size = w - msize;
    }

    // DEBUG("nchunks %d msize %d log2_block_size %d w %d \n", nchunks, msize, log2_block_size, w);
    mask.zeroize();
    for (i = 0; i < nchunks; i++)
    {
        for (j = 0; j < msize; j++)
        {
            if (mask.test_bit(j))
            {
                (*start_mask_vec)[i]->set_bit(log2_block_size + j);
            }
        }
        ++mask;
        if (i > 0)
        {
            (*end_mask_vec)[i - 1] = (*start_mask_vec)[i];
        }
    }

    (*end_mask_vec)[i - 1]->two_x(w);
    return nchunks;
}


static int parallel_mask_gen_1 (TDTree *T, int k, int bagsize, 
                                vector<bigint_t *> *start_mask_vec, vector<bigint_t *> *end_mask_vec, int nchunks)
{
    int w;
    int i;
    int j;
    int msize;
    int log2_block_size;
    int minsize;
    int kk = 0;

    bigint_t mask(T->num_mask_words);

    w = bagsize;
    for (i = 0; i < nchunks; i++)
    {
        (*start_mask_vec)[i] = new bigint_t (T->num_mask_words);
        (*end_mask_vec)[i] = new bigint_t (T->num_mask_words);

        (*start_mask_vec)[i]->zeroize();
        (*end_mask_vec)[i]->zeroize();
    }


    minsize = w - LOG2SIZE(MINENTRIES);
    if (minsize < 0)
    {
        (*end_mask_vec)[0]->two_x(w);
        return nchunks;
    }

    if (nchunks > minsize)
    {
        nchunks = minsize;
        msize = minsize;
        log2_block_size = w - msize;
    }
    else
    {

        kk = w - CHUNKSIZE;
        log2_block_size = kk;
        msize = CHUNKSIZE;
    }


    mask.zeroize();
    for (i = 0; i < nchunks; i++)
    {
        mask.zeroize();
        for (j = 0; 2*j < i*(i+1); j++)
        {
            ++mask;
        }

        for (j = 0; j < msize; j++)
        {
            if (mask.test_bit(j))
            {
                (*start_mask_vec)[i]->set_bit(log2_block_size + j);
            }
        }

        if (i > 0)
        {
            (*end_mask_vec)[i - 1] = (*start_mask_vec)[i];
        }
    }

    (*end_mask_vec)[nchunks - 1]->two_x(w);

    return nchunks;
}


static unsigned long long *encode_weights(vector<double> *weights_vec, int *trows, int *ewords)
{
    int encode_words;
    int outgoing_size;
    int entries = (int) weights_vec->size();
    double max_elem = *(max_element(weights_vec->begin(), weights_vec->end()));
    int table_rows = 0;
    unsigned long long *table_words;
    long long value = 0;
    int i, j = 0;

    if ((int)max_elem > 0)
        table_rows = (int) floor(log(max_elem)/log(2.0)) + 1;
    encode_words = 0;
    if (entries >= 0)
        encode_words = (int)ceil(((double)entries)/BIGINT_WORD_SIZE);

    *trows = table_rows;
    *ewords = encode_words;

    vector<bigint_t *> encode_table_value(table_rows);
    outgoing_size = sizeof(unsigned long long)*encode_words*table_rows;
    table_words = (unsigned long long *) malloc (outgoing_size);
    //bzero(table_words, outgoing_size);
	memset(table_words,0, outgoing_size*sizeof(unsigned long long));
    // DEBUG("table_rows %d encode_words %d\n", table_rows, encode_words);

    for(i = 0; i < table_rows; i++)
    {
        value = (long long)1 << i;
        encode_table_value[i]=new bigint_t (encode_words);
        encode_table_value[i]->zeroize();
        for (j = 0; j < entries; j++)
        {
            if ((int)(*weights_vec)[j] & (value))
            {
                encode_table_value[i]->set_bit(j);
            }
        }
    }

    for (i = 0; i < table_rows; i++)
    {
        for (j = 0; j < encode_words; j++)
        {
            table_words[i*encode_words + j] = encode_table_value[i]->words[j];
        }
    }

    return table_words;
}


static vector<int> decode_weights(unsigned long long *table_words, double *max_elem, int *entries)
{
    int encode_words;
    int table_rows = 0;
    long long value = 0;
    int i, j = 0;

    table_rows = 0;
    if ((int)*max_elem > 0)
        table_rows = (int) floor(log(*max_elem)/log(2.0)) + 1;

    encode_words = 0;
    if (*entries > 0)
        encode_words = (int)ceil(((double)*entries)/BIGINT_WORD_SIZE);

    vector<bigint_t *> decode_table_value(table_rows);
    vector<int> dweights(*entries, 0);
  
    for(i = 0; i < table_rows; i++)
    {
        decode_table_value[i] = new bigint_t (encode_words);
        for (j = 0; j < encode_words; j++)
        {
            decode_table_value[i]->words[j] = table_words[i*encode_words + j];
        }
    }
    
    for(i = 0; i < table_rows; i++)
    {
        value = 0;
        value = (long long)1 << i;
        for (j = 0; j < *entries; j++)
        {
            if (decode_table_value[i]->test_bit(j))
            {
                dweights[j] |= value;
            }
        }
    }

    return dweights;
}


static int list_ind_sets_parallel (TDTree *T, int k, list<int> *bag_subset, list<TDSolution *> *ind_sets, bigint_t *start_mask, bigint_t *end_mask)
{

    //DEBUG("in list ind sets parallel method\n");
    //start_mask->print(__log_file__);
    //end_mask->print(__log_file__);

    if( k<0 || !(T->tree_nodes[k]))
    {
        FERROR("%s: tried to compute table for impossible tree node\n",__FUNCTION__);
        MPI_Finalize();
        exit (-1);
    }

    int i,j;
    TDSolution *new_set;
    list<int>::iterator ii,jj,kk;
    bigint_t current_mask(T->num_mask_words);
    vector<int> pos_vec(bag_subset->size(),-1), k_vec(T->tree_nodes[k]->bag.size(),-1);

    i=0;
    for(ii=T->tree_nodes[k]->bag.begin();ii!=T->tree_nodes[k]->bag.end();++ii)
    {
        k_vec[i]=*ii;
        i++;
    }

    // For each position i in the list bag_subset, pos_vec[i] is the position of the ith entry
    // in bag_subset in tree node k's bag
    jj=T->tree_nodes[k]->bag.begin();
    i=j=0;
    for(ii=bag_subset->begin();ii!=bag_subset->end();++ii)
    {
        while(*ii!=*jj)
        {
            j++;
            ++jj;
        }
        pos_vec[i]=j;
        i++;
    }

    // Clear the ind_sets list
    ind_sets->clear();
    int w= (int)bag_subset->size();
    // Copy the bag (stored as a list) to a w-long bag_vec
    vector<int> bag_vec(w,0);
    i=0;
    for(ii=bag_subset->begin();ii!=bag_subset->end();++ii)
    {
        bag_vec[i]=*ii;
        i++;
    }

    // Add the empty set to the list
    /* new_set=new TDSolution(T->num_mask_words);
     *(new_set->mask)=(BIGINT_WORD)0;
     new_set->value=0;
     ind_sets->push_back(new_set);*/

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i=0;i<w;i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j=0;j<w;j++)
        {
            // Consider all poss. edges bag_vec[i]-bag_vec[j]
            nbrs = nodes[bag_vec[i]].get_nbrs();
            for(ii=nbrs.begin(); ii!=nbrs.end();++ii)
            {
                if(i!=j && bag_vec[j]==*ii)
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    // No point looking further, so break out
                    break;
                }
            }
        }
    }

    // Sort the masks and put the densest first 
    //sort(nbr_mask_vec.begin(),nbr_mask_vec.end(),hamming_compare);
    long long value = 0;
    bool is_independent;
    current_mask=0;

    bigint_t two_pow_w(T->num_mask_words);
    two_pow_w = *end_mask;
    current_mask = *start_mask;    

    bool has_bit;
    while(current_mask < two_pow_w)
    {
        is_independent=true;
        for(i=0; i<w ; i++) 
        {            
            if(current_mask.test_bit(nbr_mask_vec[i]->k))
            {
                // current_mask has the bit corresponding to nbr_mask_vec[i] set - so and with the appropriate vector for a hit
                has_bit=false;
                // Do the loop below manually rather than calling the overloaded AND
                for(j=0;j<current_mask.S;j++)
                {
                    if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                    {
                        // Some nbr is found in word j
                        has_bit=true;
                        break;
                    }
                }

                if( has_bit )
                {                    
                    is_independent=false;
                    int rpos=0;
                    // Find the rightmost bit in current_mask - this could be sped up
                    while(!(current_mask.test_bit(rpos)))
                        rpos++;

                    // Do the add w/ only bitwise op's
                    // Would be faster here to figure out which word we are in and 
                    // do this w/o test_bit
                    int b=0;
                    while(current_mask.test_bit(rpos+b))
                        b++;

                    for(int a=0;a<b;a++)
                        current_mask.xor_bit(rpos+a);

                    current_mask.or_bit(rpos+b); 

                    // Break out of the for loop since current_mask cannot represent
                    // an independent set
                    break;                   
                }
            }
        }


        if(is_independent)
        {            
            new_set=new TDSolution(T->num_mask_words);
            // Now translate the mask by translating the bits from current_mask into bits in the parent's full bag
            *(new_set->mask)=0;
            new_set->value=0;
            for(j=0;j<(int)T->tree_nodes[k]->bag.size();j++)
            {
                if(current_mask.test_bit(j))
                {
                    // bit j is set in current_mask - that corresponds to bit pos_vec[j] in a full mask
                    new_set->mask->set_bit(pos_vec[j]);
                }
            }
            // Add the new set to the list here -
            ind_sets->push_back(new_set);

            // Advance to next candidate ind. set
            ++current_mask;
        }            
    }

	for(i=0;i<w;i++)
		delete nbr_mask_vec[i];

	// Return the size of the list
	return ind_sets->size();

}

static int compute_nonnice_table_static_lb (TDTree *T, int k)
{

    vector<int>::iterator pci;
    vector<int>::iterator itv;
    vector<int> I (T->tree_nodes[k]->bag.size());

    list<int> pc_isec;
    list<int> pc_diff;
    list<TDSolution *> ind_sets;
    list<TDSolution *> ind_sets_parallel;
    list<TDSolution *> co_ind_sets;
    list<TDSolution *>::iterator itl;
    list<int>::iterator it_li;

    int i, j, kk;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();
    int indw = 0;
    int chunksize;

    int new_value = 0;
    int iweight = 0;


    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes());

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();


    pci = set_intersection(T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end(), 
                           T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end(), I.begin());

    // save parent child intersection into a list so that we can use
    // it to compute all IS.
    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_isec.push_back(*itv);
    }

    // computing all IS in parent child intersection.
    pc_isec.sort ();

    I.clear();
    pci = set_difference (T->tree_nodes[k]->bag.begin (), T->tree_nodes[k]->bag.end (),
                          pc_isec.begin (), pc_isec.end(), I.begin());

    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_diff.push_back(*itv);
    }


    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < (int)T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }



    // All worker nodes require to have all IS only in child, then
    // they can extend IS in the intersection to findout maximum
    // weight. 
    DEBUG("computing IS only in the bag\n");
    DEBUG("%d elements in the set\n", (int)pc_diff.size());
    if (!pc_diff.empty())
    {
        list_ind_sets (T, k, &pc_diff, &co_ind_sets);
    }
    else
    {
        TDSolution *new_set=new TDSolution(T->num_mask_words);
        *(new_set->mask)=0;
        new_set->value=0;
        co_ind_sets.push_back(new_set);
    }
    DEBUG("finished IS computation\n");
    // Design Idea:
    // ============ 
    // We can  let all process to compute independent sets in the
    // intersection and later determine subset which belong to them OR
    // we can compute IS in intersection at head node and distribute
    // other values using MPI_Scatter to workers and later collect
    // results using MPI_Gather.
    // At the moment we do not know which method is going to perform
    // better we have to implement both methods and see how it
    // performs. 

    // Important to note that these ind_sets mask is computed wrt
    // child's bag. 


    // DEBUG("computing IS in the  intersection\n");
    //list_ind_sets (T, k, &pc_isec, &ind_sets);
    //int ind_size = (int) ind_sets.size();


    vector<bigint_t *> start_masks;
    vector<bigint_t *> end_masks;
    if (LOG2SIZE(CHUNKSIZE) > (int)pc_isec.size())
        chunksize = T->size;
    else
        chunksize = CHUNKSIZE;


    parallel_mask_gen (T, k, (int)pc_isec.size(), &start_masks, &end_masks, chunksize);

    int imask = 0;
    vector<int> weights_collection_vec;
    int wcv_size;
    DEBUG("start IS DP loop\n");
    for (imask = 0; imask < chunksize; imask++)
    {
        if (imask % T->size != T->rank)
            continue;

        ind_sets.clear();
        list_ind_sets_parallel (T, k, &pc_isec, &ind_sets, start_masks[imask], end_masks[imask]);
        int ind_size = (int) ind_sets.size();

        //if (ind_size > 0)
        //DEBUG("%d : size %d\n", imask, ind_size);

        // populating ind_sets_vec using elements in the ind_sets list
        // because vector provide better direct access to elements
        // compared to list.
        //vector<TDSolution *> ind_sets_par_vec (ind_sets_parallel.begin(), ind_sets_parallel.end());
        vector<TDSolution *> ind_sets_vec (ind_sets.begin(), ind_sets.end());
        vector<TDSolution *> co_ind_sets_vec (co_ind_sets.begin(), co_ind_sets.end());
        vector<TDSolution *> parent_masks_vec(ind_sets.begin(), ind_sets.end());

        coi_size = (int) co_ind_sets_vec.size ();
        vector<int> ind_weights (ind_size, 0);
        //DEBUG("processor: %d start: %d end : %d block size: %d\n", T->rank, rstart, rend, block_size);

        // Sub set of IS's in intersection of parent and child which
        // processor (worker) i going to work on.
        indw = 0;
        bigint_t tmp_set(T->num_mask_words);
        bigint_t c_set(T->num_mask_words);
        //print_bag (T, k);
        vector<int> weight = T->G->get_weight();
        for (i = 0; i < ind_size; i++)
        {
            c_set.zeroize ();
            new_value = 0;
            ind_weights[indw] = -1;

            if (k == T->root_node)
            {
                for (j = 0; j < w; j++)
                {
                    if (ind_sets_vec[i]->mask->test_bit(j))
                    {
                        new_value += weight[k_bag_vec[j]];
                    }
                
                }
            }
        
            // Iterate through the ISs' only in child bag to extend IS.
            iweight = new_value;
            //ind_sets_vec[i]->mask->print(__log_file__);
            //DEBUG("Extensions: iweight %f\n", iweight);

            for (coi = 0; coi < coi_size; coi++)
            {
                new_value = iweight;
                tmp_set.zeroize ();
                tmp_set = *(ind_sets_vec[i]->mask);
                for (j = 0; j < w; j++)
                {
                    // Fist we check whether jth node is in the child only
                    // IS. If it is in the child only IS, we get that node
                    // and check whether is adjacent to any node in the
                    // current IS. if not add it to the current IS.

                    if (co_ind_sets_vec[coi]->mask->test_bit(j))
                    {
                        if (!(mask_compare(ind_sets_vec[i]->mask, nbr_mask_vec[j]->w)))
                        {
                            //GEN("%d ", k_bag_vec[j]);
                            new_value += weight[k_bag_vec[j]];
                            tmp_set.set_bit(j);
                        }
                    }
                }
                //GEN("\n");

                // Now tmp_set has the mask for complete IS, iterate
                // through children and find the maximum weight and then
                // substract weight from the parent before updating
                // maximum weight.

                // Iterate through children and get the maximum value.
                // if this is leaf node, will not execute the body.
                TDSolution *hash_ptr;
                //DEBUG("k: %d # of children : %d \n", k,
                //T->tree_nodes[k]->adj.size() - 1);

                for (kk = 0; kk < (int)T->tree_nodes[k]->adj.size() - 1; kk++)
                {
                    c_set.zeroize();
                    for (j = 0; j < w; j++)
                    {
                        if (cvec[kk][k_bag_vec[j]] && tmp_set.test_bit(j))
                        {
                            c_set.set_bit(j);
                        }
                    }
                
                    // Now we have IS which is in child hash table.
                    //DEBUG("Looking for  mask ");
                    //c_set.print(__log_file__);
                    HASH_FIND (hh, T->tree_nodes[children[kk]]->hash_table, c_set.words, 
                               T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                    if (!hash_ptr)
                        FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                    new_value += hash_ptr->value;
                }

                if (new_value > ind_weights[indw])
                {
                    ind_weights[indw] = new_value;
                }
            }
            //GEN("\n");

            indw ++;
        }
        weights_collection_vec.insert(weights_collection_vec.end(), ind_weights.begin(), ind_weights.end());
    }

    wcv_size = weights_collection_vec.size();
    // Finish loop ind set weight calculation
    DEBUG("finish IS DP loop\n");
    T->free_children(k);


    //print_bag(T, k);
    // encoding table into smaller space.
    //vector<double> max_elem_size_vec(T->size*2);
    //vector<double> max_elem_size(2, 0.0);
    //max_elem_size[0] = *(max_element(ind_weights.begin(), ind_weights.end()));
    //max_elem_size[1] = (double)ind_size;
    vector<int> n_elem_vec(T->size);
        
    //unsigned long long *table_words = encode_weights (&ind_weights, &table_rows, &encode_words);
    //DEBUG("table rows %d encode words %d\n", table_rows, encode_words);
    //int outgoing_size = sizeof(unsigned long long)*table_rows*encode_words;
    //DEBUG("entries %d outgoing size %d bytes\n", ind_size,
    //outgoing_size);

    DEBUG("processor %d waiting on gather\n", T->rank);
    MPI_Allgather((void *)&wcv_size, 1, MPI_INT, (void *)&n_elem_vec.front(), 1, MPI_INT, 
                  MPI_COMM_WORLD);
    
    //DEBUG("send weights i %d\n", T->rank);
    //DEBUG("mpi send bytes : %d instead of :%d bytes\n", table_rows*encode_words*sizeof(unsigned long long), 
    //     sizeof(int)*ind_weights.size());

    // Now preparing for receiving data. Once we recieve data we can
    // update the hash table straigh away.

    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elem_vec[i];
        w_recv_counts[i] = n_elem_vec[i];
        //DEBUG("i : %d count %d\n", i, n_elem_vec[i]);
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }


    //total_data_size = accumulate (w_recv_counts.begin(), w_recv_counts.end(), 0);
    //DEBUG("total data size : %d\n", total_data_size);
    //DEBUG("total entries : %d\n", total_entries);
    // vector<unsigned long long> complete_table(total_data_size,0);
    /* unsigned long long *complete_table = (unsigned long long *) malloc (sizeof(unsigned long long) * total_data_size);
       if (!complete_table)
       {
       FERROR ("no enough memory \n");
       MPI_Finalize ();
       exit (-1);
       }*/

    vector<int> ind_weights_recv(total_entries, -1);

    // Note
    // =====
    // Given that weights table become huge we should pursue other
    // methods (which is a big question right now) to distribute
    // weights to processors.  

    // Now parent collects all weights into ind_weights_recv vector,
    // since it is in order we can update before inserting into hash
    // table. 

    // ind_weights.front () is very important as it gives a pointer to
    // the first element. Spent lot of time to figure this out.

    //MPI_Allgatherv((void *)&ind_weights.front(), recv_counts[T->rank], MPI_DOUBLE,
    //              (void *)&ind_weights_recv.front(), &recv_counts.front(), &recv_disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);
    //MPI_Barrier (MPI_COMM_WORLD);


    DEBUG("waiting on all gatherV\n");
    MPI_Allgatherv((void *)&weights_collection_vec.front(), w_recv_counts[T->rank], MPI_DOUBLE,
                   (void *)&ind_weights_recv.front(), &w_recv_counts.front(), &w_recv_disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);
    DEBUG("all gather finished\n");

    list<TDSolution *>::iterator it_td;
    vector<TDSolution *> all_ind_sets_vec (total_entries);
    list<TDSolution *> tmp_ind_sets;
    vector <int> offset_counter(T->size,  0);
    int mod_val;
    
    if (k == T->root_node)
    {
        GEN("maximum independent set value : %f\n", *(max_element(ind_weights_recv.begin(), ind_weights_recv.end())));
        CRIT("maximum independent set value : %f\n", *(max_element(ind_weights_recv.begin(), ind_weights_recv.end())));
        fprintf(stderr, "maximum independent set value : %f\n", *(max_element(ind_weights_recv.begin(), ind_weights_recv.end())));
    }
    else
    {
        
        for (imask = 0; imask < chunksize; imask++)
        {
            tmp_ind_sets.clear();
            list_ind_sets_parallel (T, k, &pc_isec, &tmp_ind_sets, start_masks[imask], end_masks[imask]);
            mod_val = imask % T->size;
            
            //start_masks[imask]->print(__log_file__);
            //end_masks[imask]->print(__log_file__);
            //DEBUG("processor %d : \t size %d \n", mod_val, tmp_ind_sets.size());
            //DEBUG("start offset counter %d\n", offset_counter[mod_val]);
            //DEBUG("location in all weights start from %d\n", w_recv_disp[mod_val] + offset_counter[mod_val]);
            for (it_td = tmp_ind_sets.begin(); it_td != tmp_ind_sets.end(); it_td++)
            {
                all_ind_sets_vec[w_recv_disp[mod_val] + offset_counter[mod_val]] = *it_td;
                offset_counter[mod_val] ++;
            }
            //DEBUG("end offset counter %d\n", offset_counter[mod_val]);
        }

        DEBUG("Finished:mask generation for other processors\n");
    }

    //list<TDSolution *> test_ind;
    //list_ind_sets(T, k, &pc_isec, &test_ind);
    //DEBUG("size of all masks %d total_entries %d\n", all_ind_sets_vec.size(), total_entries);
    //DEBUG("recevied weights %d \n", ind_weights_recv.size());


    // At this point ind_weights_recv vector contains best values for
    // each independent set in the intersection, therefore we can
    // directly update the hash table for the node (in the decomposed
    // tree).

    if (k != T->root_node)
    {
        TDSolution *add_set;
        bigint_t pmask(T->num_mask_words);
            
        DEBUG("updating hash table\n");
        for (i = 0; i < (int)all_ind_sets_vec.size(); i++)
        {
            pmask.zeroize ();
            for (j = 0; j < w; j++)
            {
                if (all_ind_sets_vec[i]->mask->test_bit(j))
                {
                    pmask.set_bit(parent_location_vec[k_bag_vec[j]]);
                }
            }
            add_set = new TDSolution (T->num_mask_words);
            *(add_set->mask) = pmask;
            add_set->value = ind_weights_recv[i];
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, add_set->mask->words, 
                            T->num_mask_words*sizeof(BIGINT_WORD), add_set);
        }
    }
    
    return 0;
}
            
int compute_nonnice_table_parallel_0 (TDTree *T, int k)
{
    double stat_usec = MPI_Wtime ();
    vector<int>::iterator pci;
    vector<int>::iterator itv;
    vector<int> I (T->tree_nodes[k]->bag.size());

    list<int> pc_isec;
    list<int> pc_diff;
    list<TDSolution *> ind_sets;
    list<TDSolution *> ind_sets_parallel;
    list<TDSolution *> co_ind_sets;
    list<TDSolution *>::iterator itl;
    list<int>::iterator it_li;

    int i, j, kk;
    int w, pw;
    int coi, coi_size;
    int is_size;
    int cblock_size, fblock_size, block_size;
    int rstart, rend;
    int parent = T->tree_nodes[k]->adj.front();
    int indw = 0;
    int position = 0;

    int max_weight = -1, new_value = 0;
    int parent_weight;
    int iweight = 0;


    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes());

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;
    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }


    pci = set_intersection(T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end(), 
                           T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end(), I.begin());

    // save parent child intersection into a list so that we can use
    // it to compute all IS.
    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_isec.push_back(*itv);
    }

    // computing all IS in parent child intersection.
    pc_isec.sort ();

    I.clear();
    pci = set_difference (T->tree_nodes[k]->bag.begin (), T->tree_nodes[k]->bag.end (),
                          pc_isec.begin (), pc_isec.end(), I.begin());

    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_diff.push_back(*itv);
    }

    // All worker nodes require to have all IS only in child, then
    // they can extend IS in the intersection to findout maximum
    // weight. 
    DEBUG("computing IS only in the bag\n");
    DEBUG("%d elements in the set\n", pc_diff.size());
    if (!pc_diff.empty())
    {
        list_ind_sets (T, k, &pc_diff, &co_ind_sets);
    }
    else
    {
        TDSolution *new_set=new TDSolution(T->num_mask_words);
        *(new_set->mask)=0;
        new_set->value=0;
        co_ind_sets.push_back(new_set);
    }
    DEBUG("finished IS computation\n");
    // Design Idea:
    // ============ 
    // We can  let all process to compute independent sets in the
    // intersection and later determine subset which belong to them OR
    // we can compute IS in intersection at head node and distribute
    // other values using MPI_Scatter to workers and later collect
    // results using MPI_Gather.
    // At the moment we do not know which method is going to perform
    // better we have to implement both methods and see how it
    // performs. 

    // Important to note that these ind_sets mask is computed wrt
    // child's bag. 


    // DEBUG("computing IS in the  intersection\n");
    //list_ind_sets (T, k, &pc_isec, &ind_sets);
    //int ind_size = (int) ind_sets.size();


    vector<bigint_t *> start_masks (T->size);
    vector<bigint_t *> end_masks (T->size);
    int ims;

    ims = parallel_mask_gen (T, k, (int)pc_isec.size(), &start_masks, &end_masks, T->size);
    DEBUG("finshed parallel mask generation\n");

    /*DEBUG("pc_isec.size %d \n", pc_isec.size());
      for (i = 0 ; i < T->size; i++)
      {
      start_masks[i]->print(__log_file__);
      end_masks[i]->print(__log_file__);
      GEN("xxxxxxxxxxxxxxxxx\n");
      }*/
    //DEBUG("processor %d\n", T->rank);
    //start_masks[T->rank]->print(__log_file__);
    //end_masks[T->rank]->print(__log_file__);

    list_ind_sets_parallel (T, k, &pc_isec, &ind_sets, start_masks[T->rank], end_masks[T->rank]);
    int ind_size = (int) ind_sets.size();
    DEBUG("total ind sets in the region %d \n", ind_size);

    // populating ind_sets_vec using elements in the ind_sets list
    // because vector provide better direct access to elements
    // compared to list.
    //vector<TDSolution *> ind_sets_par_vec (ind_sets_parallel.begin(), ind_sets_parallel.end());
    vector<TDSolution *> ind_sets_vec (ind_sets.begin(), ind_sets.end());
    vector<TDSolution *> co_ind_sets_vec (co_ind_sets.begin(), co_ind_sets.end());
    vector<TDSolution *> parent_masks_vec(ind_sets.begin(), ind_sets.end());

    coi_size = (int) co_ind_sets_vec.size ();
    vector<int> ind_weights (ind_size, 0);
    //DEBUG("processor: %d start: %d end : %d block size: %d\n", T->rank, rstart, rend, block_size);

    // Sub set of IS's in intersection of parent and child which
    // processor (worker) i going to work on.
    indw = 0;
    bigint_t tmp_set(T->num_mask_words);
    bigint_t c_set(T->num_mask_words);
    //print_bag (T, k);
    bool is_best = false;
    int ib = 0;
    DEBUG("start IS DP loop\n");
    vector<int> weight;
    weight = T->G->get_weight();

    for (i = 0; i < ind_size; i++)
    {
        c_set.zeroize ();
        new_value = 0;
        ind_weights[indw] = -1;

        if (k == T->root_node)
        {
            for (j = 0; j < w; j++)
            {
                if (ind_sets_vec[i]->mask->test_bit(j))
                {
                    //GEN("%d ", k_bag_vec[j]);
                    new_value += weight[k_bag_vec[j]];
                }
                
            }
        }
        
        // Iterate through the ISs' only in child bag to extend IS.
        iweight = new_value;

        for (coi = 0; coi < coi_size; coi++)
        {
            new_value = iweight;
            tmp_set.zeroize ();
            tmp_set = *(ind_sets_vec[i]->mask);
            for (j = 0; j < w; j++)
            {
                // Fist we check whether jth node is in the child only
                // IS. If it is in the child only IS, we get that node
                // and check whether is adjacent to any node in the
                // current IS. if not add it to the current IS.

                if (co_ind_sets_vec[coi]->mask->test_bit(j))
                {
                    if (!(mask_compare(ind_sets_vec[i]->mask, nbr_mask_vec[j]->w)))
                    {
                        //GEN("%d ", k_bag_vec[j]);
                        new_value += weight[k_bag_vec[j]];
                        tmp_set.set_bit(j);
                    }
                }
            }
            
            //GEN("\n");

            // Now tmp_set has the mask for complete IS, iterate
            // through children and find the maximum weight and then
            // substract weight from the parent before updating
            // maximum weight.

            // Iterate through children and get the maximum value.
            // if this is leaf node, will not execute the body.
            TDSolution *hash_ptr;
            //DEBUG("k: %d # of children : %d \n", k,
            //T->tree_nodes[k]->adj.size() - 1);

            
            for (kk = 0; kk < T->tree_nodes[k]->adj.size() - 1; kk++)
            {
                c_set.zeroize();
                for (j = 0; j < w; j++)
                {
                    if (cvec[kk][k_bag_vec[j]] && tmp_set.test_bit(j))
                    {
                        c_set.set_bit(j);
                    }
                }
                
                // Now we have IS which is in child hash table.
                //DEBUG("Looking for  mask ");
                //c_set.print(__log_file__);

                HASH_FIND (hh, T->tree_nodes[children[kk]]->hash_table, c_set.words, 
                           T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                if (!hash_ptr)
                {
                    FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                    FERROR("error occured at child %d  kk value %d \n", children[kk], kk);
                    FERROR("mask looked for : ");
                    c_set.print(__log_file__);
                    MPI_Finalize ();
                    exit (-1);
                }
                new_value += hash_ptr->value;
            }

            if (new_value > ind_weights[indw])
            {
                ind_weights[indw] = new_value;
            }
        }
        //GEN("\n");

        indw ++;
    }

    // Finish loop ind set weight calculation
    DEBUG("finish IS DP loop\n");
    T->free_children(k);


    vector<int> n_elem_vec(T->size);
     DEBUG("processor %d waiting on gather\n", T->rank);
    //DEBUG("ind size %d\n", ind_size);
    //DEBUG("T->size %d\n", T->size);

    MPI_Allgather((void *)&ind_size, 1, MPI_INT, (void *)&n_elem_vec.front(), 1, MPI_INT, MPI_COMM_WORLD);

     // Now preparing for receiving data. Once we recieve data we can
    // update the hash table straigh away.

    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elem_vec[i];
        w_recv_counts[i] = n_elem_vec[i];
        //DEBUG("i : %d count %d\n", i, n_elem_vec[i]);
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

 
    DEBUG("total entries %d\n", total_entries);
    vector<int> ind_weights_recv(total_entries, -1);

    // Note
    // =====
    // Given that weights table become huge we should pursue other
    // methods (which is a big question right now) to distribute
    // weights to processors.  

    // Now parent collects all weights into ind_weights_recv vector,
    // since it is in order we can update before inserting into hash
    // table. 

    // ind_weights.front () is very important as it gives a pointer to
    // the first element. Spent lot of time to figure this out.

    //MPI_Allgatherv((void *)&ind_weights.front(), recv_counts[T->rank], MPI_DOUBLE,
    //              (void *)&ind_weights_recv.front(), &recv_counts.front(), &recv_disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);
    //MPI_Barrier (MPI_COMM_WORLD);

    list<TDSolution *>::iterator it_td;
    vector<TDSolution *> all_ind_sets_vec (total_entries);
    //vector<TDSolution *> all_ind_sets_vec;
    list<TDSolution *> tmp_ind_sets;

#pragma omp parallel for private(i, j, it_td, tmp_ind_sets)
    for (i = -1; i < T->size; i++)
    {
        if (i < 0)
        {
            DEBUG("waiting on all gatherV\n");
            MPI_Allgatherv((void *)&ind_weights.front(), w_recv_counts[T->rank], MPI_DOUBLE,
                           (void *)&ind_weights_recv.front(), &w_recv_counts.front(), &w_recv_disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);
            DEBUG("all gather finished\n");
            if (k == T->root_node)
            {
                GEN("maximum independent set value : %f\n", *(max_element(ind_weights_recv.begin(), ind_weights_recv.end())));
                CRIT("maximum independent set value : %f\n", *(max_element(ind_weights_recv.begin(), ind_weights_recv.end())));
                fprintf(stderr, "maximum independent set value : %f\n", *(max_element(ind_weights_recv.begin(), ind_weights_recv.end())));
            }
        }
        else
        {

            j = 0;
            tmp_ind_sets.clear();
            list_ind_sets_parallel (T, k, &pc_isec, &tmp_ind_sets, start_masks[i], end_masks[i]);
            for (it_td = tmp_ind_sets.begin(); it_td != tmp_ind_sets.end(); it_td++)
            {
                all_ind_sets_vec[w_recv_disp[i] + j] = *it_td;
                j++;
            }
        }
    }

    DEBUG("Finished:mask generation for other processors\n");


    //list<TDSolution *> test_ind;
    //list_ind_sets(T, k, &pc_isec, &test_ind);
    //DEBUG("size of all masks %d total_entries %d\n", all_ind_sets_vec.size(), total_entries);
    //DEBUG("recevied weights %d \n", ind_weights_recv.size());


    // At this point ind_weights_recv vector contains best values for
    // each independent set in the intersection, therefore we can
    // directly update the hash table for the node (in the decomposed
    // tree).

    if (k != T->root_node)
    {
        TDSolution *add_set;
        bigint_t pmask(T->num_mask_words);

        DEBUG("updating hash table\n");
        for (i = 0; i < (int)all_ind_sets_vec.size(); i++)
        {
            pmask.zeroize ();

            for (j = 0; j < w; j++)
            {
                if (all_ind_sets_vec[i]->mask->test_bit(j))
                {
                    pmask.set_bit(parent_location_vec[k_bag_vec[j]]);
                }
            }
            add_set = new TDSolution (T->num_mask_words);
            *(add_set->mask) = pmask;
            add_set->value = ind_weights_recv[i];
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, add_set->mask->words, 
                            T->num_mask_words*sizeof(BIGINT_WORD), add_set);

        }
    }

    DEBUG("function : %s execution time %lf\n", __FUNCTION__, MPI_Wtime() - stat_usec);   
    return 0;
}

int compute_weighted_ind_sets_parallel (TDTree *T, int k)
{
    int table_size=GD_UNDEFINED;
    int node_type=T->get_node_type(k);
    double stime = 0.0;
    switch(node_type)
    {
        case TD_NICE_LEAF_NODE:
            break;
        case TD_NONNICE_LEAF_NODE:
            stime = MPI_Wtime();
            table_size = compute_nonnice_table_parallel_10 (T,k);
            T->info->leaf_time += MPI_Wtime() - stime;
            break;
        case TD_NONLEAF_NODE:
            // Use parent_child if option given
            stime = MPI_Wtime();
            table_size = compute_nonnice_table_parallel_10 (T,k);
            T->info->nonleaf_time += MPI_Wtime() - stime;
            break;
        case TD_INTRODUCE_NODE:
            break;
        case TD_FORGET_NODE:
            break;
        case TD_JOIN_NODE:
            break;
        default:
            FERROR("Unknown tree node type %d encountered when computing table of tree node %d??\n",node_type,k);
    }
    return table_size;
}

static int compute_nonnice_table_parallel_1 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
    double stat_usec = MPI_Wtime();
    vector<int>::iterator pci;
    vector<int>::iterator itv;
    vector<int> I (T->tree_nodes[k]->bag.size());

    list<int> pc_isec;
    list<int> pc_diff;
    list<TDSolution *> ind_sets;
    list<TDSolution *> ind_sets_parallel;
    list<TDSolution *> co_ind_sets;
    list<TDSolution *>::iterator itl;
    list<int>::iterator it_li;

    int i, j, kk;
    int w, pw;
    int coi, coi_size;
    int is_size;
    int cblock_size, fblock_size, block_size;
    int rstart, rend;
    int parent = T->tree_nodes[k]->adj.front();
    int indw = 0;
    int position = 0;

    int max_weight = -1, new_value = 0;
    int parent_weight;
    int iweight = 0;


    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes());

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }


    pci = set_intersection(T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end(), 
                           T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end(), I.begin());

    // save parent child intersection into a list so that we can use
    // it to compute all IS.
    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_isec.push_back(*itv);
    }

    // computing all IS in parent child intersection.
    pc_isec.sort ();

    //I.clear();
	for(i=0;i<T->tree_nodes[k]->bag.size();i++)
		I[i]=-1;

    pci = set_difference (T->tree_nodes[k]->bag.begin (), T->tree_nodes[k]->bag.end (),
                          pc_isec.begin (), pc_isec.end(), I.begin());

    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_diff.push_back(*itv);
    }

    // All worker nodes require to have all IS only in child, then
    // they can extend IS in the intersection to findout maximum
    // weight. 
    DEBUG("computing ind sets only in the bag\n");
    DEBUG("%d elements in the set\n", pc_diff.size());
    //vector<TDSolution *> co_ind_sets_vec;
    if (!pc_diff.empty())
    {
        list_ind_sets (T, k, &pc_diff, &co_ind_sets);
    }
    else
    {
        TDSolution *new_set=new TDSolution(T->num_mask_words);
        *(new_set->mask)=0;
        new_set->value=0;
        co_ind_sets.push_back(new_set);
        //co_ind_sets_vec.push_back(new_set);
    }
    DEBUG("finished IS computation\n");

    // Design Idea:
    // ============ 
    // We can  let all process to compute independent sets in the
    // intersection and later determine subset which belong to them OR
    // we can compute IS in intersection at head node and distribute
    // other values using MPI_Scatter to workers and later collect
    // results using MPI_Gather.
    // At the moment we do not know which method is going to perform
    // better we have to implement both methods and see how it
    // performs. 

    // Important to note that these ind_sets mask is computed wrt
    // child's bag. 

    vector<int> counts(T->size);
    vector<int> disp(T->size);

    // populating ind_sets_vec using elements in the ind_sets list
    // because vector provide better direct access to elements
    // compared to list.
    vector<TDSolution *> co_ind_sets_vec (co_ind_sets.begin(), co_ind_sets.end());
    //vector<TDSolution *> parent_masks_vec(ind_sets.begin(), ind_sets.end());
    //vector<TDSolution *> ind_sets_vec;

    DEBUG("computing ind sets in the intesection using list_ind_sets\n");
    list_ind_sets(T, k, &pc_isec,  &ind_sets);
    vector<TDSolution *> ind_sets_vec(ind_sets.begin(), ind_sets.end());
    int ind_size = (int)ind_sets_vec.size();
    DEBUG("total sets %d\n", ind_size);

    // populate_cond_disp function computes size of the workload each
    // processor get, and disp array contains starting position
    // (displacement from the beginning) for each processor.  
    populate_count_disp (T, ind_size, &counts, &disp);

    coi_size = (int) co_ind_sets_vec.size ();
    vector<int> ind_weights (ind_size, 0);

    // Sub set of IS's in intersection of parent and child which
    // processor (worker) i going to work on.
    indw = 0;
    bigint_t tmp_set(T->num_mask_words);
    bigint_t c_set(T->num_mask_words);
    int ib = 0;
    int localmax = -1;

    DEBUG("start IS DP loop\n");
    vector<int> weight;
    weight = T->G->get_weight();

    for (i = disp[T->rank]; i < counts[T->rank] + disp[T->rank]; i++)
    {
        c_set.zeroize ();
        new_value = 0;
        ind_weights[indw] = -1;

        if (k == T->root_node)
        {
            for (j = 0; j < w; j++)
            {
                if (ind_sets_vec[i]->mask->test_bit(j))
                {
                    new_value += weight[k_bag_vec[j]];
                }
                
            }
        }
        
        // Iterate through the ISs' only in child bag to extend IS.
        iweight = new_value;

        for (coi = 0; coi < coi_size; coi++)
        {
            new_value = iweight;
            tmp_set.zeroize ();
            tmp_set = *(ind_sets_vec[i]->mask);
            for (j = 0; j < w; j++)
            {
                // Fist we check whether jth node is in the child only
                // IS. If it is in the child only IS, we get that node
                // and check whether is adjacent to any node in the
                // current IS. if not add it to the current IS.

                if (co_ind_sets_vec[coi]->mask->test_bit(j))
                {
                    if (!(mask_compare(ind_sets_vec[i]->mask, nbr_mask_vec[j]->w)))
                    {
                        new_value += weight[k_bag_vec[j]];
                        tmp_set.set_bit(j);
                    }
                }
            }

            // Now tmp_set has the mask for complete IS, iterate
            // through children and find the maximum weight and then
            // substract weight from the parent before updating
            // maximum weight.

            // Iterate through children and get the maximum value.
            // if this is leaf node, will not execute the body.
            TDSolution *hash_ptr;
            
            for (kk = 0; kk < T->tree_nodes[k]->adj.size() - 1; kk++)
            {
                c_set.zeroize();
                for (j = 0; j < w; j++)
                {
                    if (cvec[kk][k_bag_vec[j]] && tmp_set.test_bit(j))
                    {
                        c_set.set_bit(j);
                    }
                }

                // Now we have IS which is in child hash table.
                HASH_FIND (hh, T->tree_nodes[children[kk]]->hash_table, c_set.words, 
                           T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                if (!hash_ptr)
                {
                    FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                    FERROR("error occured at child %d  kk value %d \n", children[kk], kk);
                    FERROR("mask looked for : ");
                    c_set.print(__log_file__);
                    MPI_Finalize ();
                    exit (-1);
                }
                new_value += hash_ptr->value;
            }

            
            if (k == T->root_node)
            {
                // At root node we just need maximum value. Each
                // processor send their local  maximum value to root
                // and root picks the global maximum.
                if (new_value > localmax)
                {
                    localmax = new_value;
                }
            }
            else
            {
                if (new_value > ind_weights[indw])
                {
                    ind_weights[indw] = new_value;
                }
            }
        }
        //GEN("\n");
        indw ++;
    }

    // Finish loop ind set weight calculation
    DEBUG("finish IS DP loop\n");
    T->free_children(k);


    vector<int> ind_weights_recv(ind_size, -1);

    // Note
    // =====
    // Given that weights table become huge we should pursue other
    // methods (which is a big question right now) to distribute
    // weights to processors.  

    // Now parent collects all weights into ind_weights_recv vector,
    // since it is in order we can update before inserting into hash
    // table. 

    // ind_weights.front () is very important as it gives a pointer to
    // the first element. Spent lot of time to figure this out.

    if (k == T->root_node)
    {
        DEBUG("preparing for reduce operation\n");
        int max;
        MPI_Reduce ((void *)&localmax, (void *) &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (T->rank == 0)
        {
            CRIT("Maximum independent set value : %f \n", max);
            fprintf(stderr, "maximum independent set value : %f\n", max);
        }
    }
    else
    {
        DEBUG("waiting on all gatherV\n");
        MPI_Allgatherv((void *)&ind_weights.front(), counts[T->rank], MPI_DOUBLE,
                       (void *)&ind_weights_recv.front(), &counts.front(), &disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);
        DEBUG("all gather finished\n");
    }
    


    // each independent set in the intersection, therefore we can
    // directly update the hash table for the node (in the decomposed
    // tree).

    if (k != T->root_node)
    {
        TDSolution *add_set;
        bigint_t pmask(T->num_mask_words);

        DEBUG("updating hash table\n");
        for (i = 0; i < (int)ind_sets_vec.size(); i++)
        {
            pmask.zeroize ();

            for (j = 0; j < w; j++)
            {
                if (ind_sets_vec[i]->mask->test_bit(j))
                {
                    pmask.set_bit(parent_location_vec[k_bag_vec[j]]);
                }
            }
            add_set = new TDSolution (T->num_mask_words);
            *(add_set->mask) = pmask;
            add_set->value = ind_weights_recv[i];
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, add_set->mask->words, 
                            T->num_mask_words*sizeof(BIGINT_WORD), add_set);
        }
    }

    DEBUG("function : %s execution time %lf\n", __FUNCTION__, MPI_Wtime() - stat_usec);
    return 0;
}

static int compute_nonnice_table_parallel_2 (TDTree *T, int k)
{

    DEBUG("In function : %s \n", __FUNCTION__);
    double stat_usec = MPI_Wtime();
    vector<int>::iterator pci;
    vector<int>::iterator itv;
    vector<int> I (T->tree_nodes[k]->bag.size());

    list<int> pc_isec;
    list<int> pc_diff;
    list<TDSolution *> ind_sets;
    list<TDSolution *> ind_sets_parallel;
    list<TDSolution *> co_ind_sets;
    list<TDSolution *>::iterator itl;
    list<int>::iterator it_li;

    int i, j, kk;
    int w, pw;
    int coi, coi_size;
    int is_size;
    int cblock_size, fblock_size, block_size;
    int rstart, rend;
    int parent = T->tree_nodes[k]->adj.front();
    int indw = 0;
    int position = 0;

    int max_weight = -1, new_value = 0;
    int parent_weight;
    int iweight = 0;


    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes());

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }


    pci = set_intersection(T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end(), 
                           T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end(), I.begin());

    // save parent child intersection into a list so that we can use
    // it to compute all IS.
    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_isec.push_back(*itv);
    }

    // computing all IS in parent child intersection.
    pc_isec.sort ();

    I.clear();
    pci = set_difference (T->tree_nodes[k]->bag.begin (), T->tree_nodes[k]->bag.end (),
                          pc_isec.begin (), pc_isec.end(), I.begin());

    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_diff.push_back(*itv);
    }

    // All worker nodes require to have all IS only in child, then
    // they can extend IS in the intersection to findout maximum
    // weight. 
    DEBUG("computing ind sets only in the bag\n");
    DEBUG("%d elements in the set\n", pc_diff.size());

    if (!pc_diff.empty())
    {
        list_ind_sets (T, k, &pc_diff, &co_ind_sets);
    }
    else
    {
        TDSolution *new_set=new TDSolution(T->num_mask_words);
        *(new_set->mask)=0;
        new_set->value=0;
        co_ind_sets.push_back(new_set);
    }
    DEBUG("finished IS computation\n");


    // Design Idea:
    // ============ 
    // We can  let all process to compute independent sets in the
    // intersection and later determine subset which belong to them OR
    // we can compute IS in intersection at head node and distribute
    // other values using MPI_Scatter to workers and later collect
    // results using MPI_Gather.
    // At the moment we do not know which method is going to perform
    // better we have to implement both methods and see how it
    // performs. 

    // Important to note that these ind_sets mask is computed wrt
    // child's bag. 

    vector<int> counts(T->size);
    vector<int> disp(T->size);

    // populating ind_sets_vec using elements in the ind_sets list
    // because vector provide better direct access to elements
    // compared to list.
    vector<TDSolution *> co_ind_sets_vec (co_ind_sets.begin(), co_ind_sets.end());

    DEBUG("computing ind sets in the intesection using generate_ind_sets\n");

    vector<TDSolution *> ind_sets_vec;
    generate_ind_sets (T, k, &pc_isec, &ind_sets_vec);
    int ind_size = (int)ind_sets_vec.size();
    DEBUG("total sets %d\n", ind_size);

    // populate_cond_disp function computes size of the workload each
    // processor get, and disp array contains starting position
    // (displacement from the beginning) for each processor.  
    populate_count_disp (T, ind_size, &counts, &disp);

    coi_size = (int) co_ind_sets_vec.size ();
    vector<int> ind_weights (ind_size, 0);

    // Sub set of IS's in intersection of parent and child which
    // processor (worker) i going to work on.
    indw = 0;
    bigint_t tmp_set(T->num_mask_words);
    bigint_t c_set(T->num_mask_words);
    int ib = 0;
    int localmax = -1;

    DEBUG("start IS DP loop\n");
    vector<int>weight;
    weight = T->G->get_weight();
    for (i = disp[T->rank]; i < counts[T->rank] + disp[T->rank]; i++)
    {
        c_set.zeroize ();
        new_value = 0;
        ind_weights[indw] = -1;

        if (k == T->root_node)
        {
            for (j = 0; j < w; j++)
            {
                if (ind_sets_vec[i]->mask->test_bit(j))
                {
                    new_value += weight[k_bag_vec[j]];
                }
                
            }
        }
        
        // Iterate through the ISs' only in child bag to extend IS.
        iweight = new_value;

        for (coi = 0; coi < coi_size; coi++)
        {
            new_value = iweight;
            tmp_set.zeroize ();
            tmp_set = *(ind_sets_vec[i]->mask);
            for (j = 0; j < w; j++)
            {
                // Fist we check whether jth node is in the child only
                // IS. If it is in the child only IS, we get that node
                // and check whether is adjacent to any node in the
                // current IS. if not add it to the current IS.

                if (co_ind_sets_vec[coi]->mask->test_bit(j))
                {
                    if (!(mask_compare(ind_sets_vec[i]->mask, nbr_mask_vec[j]->w)))
                    {
                        new_value += weight[k_bag_vec[j]];
                        tmp_set.set_bit(j);
                    }
                }
            }

            // Now tmp_set has the mask for complete IS, iterate
            // through children and find the maximum weight and then
            // substract weight from the parent before updating
            // maximum weight.

            // Iterate through children and get the maximum value.
            // if this is leaf node, will not execute the body.
            TDSolution *hash_ptr;
            
            for (kk = 0; kk < T->tree_nodes[k]->adj.size() - 1; kk++)
            {
                c_set.zeroize();
                for (j = 0; j < w; j++)
                {
                    if (cvec[kk][k_bag_vec[j]] && tmp_set.test_bit(j))
                    {
                        c_set.set_bit(j);
                    }
                }
                
                // Now we have IS which is in child hash table.

                HASH_FIND (hh, T->tree_nodes[children[kk]]->hash_table, c_set.words, 
                           T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                if (!hash_ptr)
                {
                    FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                    FERROR("error occured at child %d  kk value %d \n", children[kk], kk);
                    FERROR("mask looked for : ");
                    c_set.print(__log_file__);
                    MPI_Finalize ();
                    exit (-1);
                }
                new_value += hash_ptr->value;
            }

            if (k == T->root_node)
            {
                if (new_value > localmax)
                {
                    localmax = new_value;
                }
            }
            else
            {
                if (new_value > ind_weights[indw])
                {
                    ind_weights[indw] = new_value;
                }
            }
        }
        indw ++;
    }

    // Finish loop ind set weight calculation
    DEBUG("finish IS DP loop\n");
    T->free_children(k);

    vector<int> ind_weights_recv(ind_size, -1);

    // Note
    // =====
    // Given that weights table become huge we should pursue other
    // methods (which is a big question right now) to distribute
    // weights to processors.  

    // Now parent collects all weights into ind_weights_recv vector,
    // since it is in order we can update before inserting into hash
    // table. 

    // ind_weights.front () is very important as it gives a pointer to
    // the first element. Spent lot of time to figure this out.

    if (k == T->root_node)
    {
        DEBUG("preparing for reduce operation\n");
        int max;
        MPI_Reduce ((void *)&localmax, (void *) &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (T->rank == 0)
        {
            CRIT("Maximum independent set value : %f \n", max);
            fprintf(stderr, "maximum independent set value : %f\n", max);
        }
    }
    else
    {
        DEBUG("waiting on all gatherV\n");
        MPI_Allgatherv((void *)&ind_weights.front(), counts[T->rank], MPI_DOUBLE,
                       (void *)&ind_weights_recv.front(), &counts.front(), &disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);
        DEBUG("all gather finished\n");
    }
    


    // each independent set in the intersection, therefore we can
    // directly update the hash table for the node (in the decomposed
    // tree).

    if (k != T->root_node)
    {
        TDSolution *add_set;
        bigint_t pmask(T->num_mask_words);

        DEBUG("updating hash table\n");
        for (i = 0; i < (int)ind_sets_vec.size(); i++)
        {
            pmask.zeroize ();

            for (j = 0; j < w; j++)
            {
                if (ind_sets_vec[i]->mask->test_bit(j))
                {
                    pmask.set_bit(parent_location_vec[k_bag_vec[j]]);
                }
            }

            add_set = new TDSolution (T->num_mask_words);
            *(add_set->mask) = pmask;
            add_set->value = ind_weights_recv[i];
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, add_set->mask->words, 
                            T->num_mask_words*sizeof(BIGINT_WORD), add_set);
        }
    }
    
    DEBUG("function : %s execution time %lf\n", __FUNCTION__, MPI_Wtime() - stat_usec);
    return 0;
}

int compute_nonnice_table_test (TDTree *T, int k)
{
    double stat_usec;
    vector<int>::iterator pci;
    vector<int>::iterator itv;
    vector<int> I (T->tree_nodes[k]->bag.size());

    list<int> pc_isec;
    list<int> pc_diff;
    list<TDSolution *> ind_sets;
    list<TDSolution *> ind_sets_parallel;
    list<TDSolution *> co_ind_sets;
    list<TDSolution *>::iterator itl;
    list<int>::iterator it_li;

    int i, j, kk;
    int w, pw;
    int coi, coi_size;
    int is_size;
    int cblock_size, fblock_size, block_size;
    int rstart, rend;
    int parent = T->tree_nodes[k]->adj.front();
    int indw = 0;
    int position = 0;

    int max_weight = -1, new_value = 0;
    int parent_weight;
    int iweight = 0;


    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes());

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }


    pci = set_intersection(T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end(), 
                           T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end(), I.begin());

    // save parent child intersection into a list so that we can use
    // it to compute all IS.
    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_isec.push_back(*itv);
    }

    // computing all IS in parent child intersection.
    pc_isec.sort ();

    I.clear();
    pci = set_difference (T->tree_nodes[k]->bag.begin (), T->tree_nodes[k]->bag.end (),
                          pc_isec.begin (), pc_isec.end(), I.begin());

    for (itv = I.begin(); itv != pci; ++itv)
    {
        pc_diff.push_back(*itv);
    }

    // All worker nodes require to have all IS only in child, then
    // they can extend IS in the intersection to findout maximum
    // weight. 
    //DEBUG("computing ind sets only in the bag\n");
    //DEBUG("%d elements in the set\n", pc_diff.size());
    //vector<TDSolution *> co_ind_sets_vec;
    if (!pc_diff.empty())
    {
        list_ind_sets (T, k, &pc_diff, &co_ind_sets);
        //generate_ind_sets (T, k, &pc_diff, &co_ind_sets_vec);
    }
    else
    {
        TDSolution *new_set=new TDSolution(T->num_mask_words);
        *(new_set->mask)=0;
        new_set->value=0;
        co_ind_sets.push_back(new_set);
        //co_ind_sets_vec.push_back(new_set);
    }
    //DEBUG("finished IS computation\n");


    // Design Idea:
    // ============ 
    // We can  let all process to compute independent sets in the
    // intersection and later determine subset which belong to them OR
    // we can compute IS in intersection at head node and distribute
    // other values using MPIScatter to workers and later collect
    // results using MPI_Gather.
    // At the moment we do not know which method is going to perform
    // better we have to implement both methods and see how it
    // performs. 

    // Important to note that these ind_sets mask is computed wrt
    // child's bag. 

    vector<int> counts(T->size);
    vector<int> disp(T->size);

    // populating ind_sets_vec using elements in the ind_sets list
    // because vector provide better direct access to elements
    // compared to list.
    vector<TDSolution *> co_ind_sets_vec (co_ind_sets.begin(), co_ind_sets.end());
    //vector<TDSolution *> parent_masks_vec(ind_sets.begin(), ind_sets.end());
    //vector<TDSolution *> ind_sets_vec;

    //DEBUG("computing ind sets in the intesection with %d chunks\n", CHUNKSIZE);
    //list<TDSolution *> ind_sets;
    vector<TDSolution *> ind_sets_vec;
    stat_usec = MPI_Wtime();
    generate_ind_sets (T, k, &pc_isec, &ind_sets_vec);
    //DEBUG("Finished ind sets generation\n");
    DEBUG("execution time : %lf\n", MPI_Wtime() - stat_usec);

    DEBUG("start sequential computation\n");
    stat_usec = MPI_Wtime();
    list_ind_sets(T, k, &pc_isec,  &ind_sets);
    DEBUG("execution time : %lf \n", MPI_Wtime() - stat_usec);
    DEBUG("finish sequential computation\n");
    //vector<TDSolution *> ind_sets_vec(ind_sets.begin(), ind_sets.end());
    int ind_size = (int)ind_sets_vec.size();
    DEBUG("total sets %d\n", ind_size);
    if (ind_size != (int)ind_sets.size())
    {
        FERROR("sequential total sets %d\n", ind_sets.size());
        MPI_Finalize();
        exit(-1);
    }
    //list<TDSolution *> seqlist;
    //list_ind_sets(T, k, &pc_isec,  &seqlist);
    //DEBUG("total sets sequential %d\n", seqlist.size());

    populate_count_disp (T, ind_size, &counts, &disp);

    coi_size = (int) co_ind_sets_vec.size ();
    vector<int> ind_weights (ind_size, 0);

    // Sub set of IS's in intersection of parent and child which
    // processor (worker) i going to work on.
    indw = 0;
    bigint_t tmp_set(T->num_mask_words);
    bigint_t c_set(T->num_mask_words);
    int ib = 0;
    int localmax = -1;

    // DEBUG("start IS DP loop\n");
    vector<int> weight = T->G->get_weight();

    for (i = disp[T->rank]; i < counts[T->rank] + disp[T->rank]; i++)
    {
        c_set.zeroize ();
        new_value = 0;
        ind_weights[indw] = -1;

        if (k == T->root_node)
        {
            for (j = 0; j < w; j++)
            {
                if (ind_sets_vec[i]->mask->test_bit(j))
                {
                    //GEN("%d ", k_bag_vec[j]);
                    new_value += weight[k_bag_vec[j]];
                }
                
            }
        }
        
        // Iterate through the ISs' only in child bag to extend IS.
        iweight = new_value;

        for (coi = 0; coi < coi_size; coi++)
        {
            new_value = iweight;
            tmp_set.zeroize ();
            tmp_set = *(ind_sets_vec[i]->mask);
            for (j = 0; j < w; j++)
            {
                // Fist we check whether jth node is in the child only
                // IS. If it is in the child only IS, we get that node
                // and check whether is adjacent to any node in the
                // current IS. if not add it to the current IS.

                if (co_ind_sets_vec[coi]->mask->test_bit(j))
                {
                    if (!(mask_compare(ind_sets_vec[i]->mask, nbr_mask_vec[j]->w)))
                    {
                        new_value += weight[k_bag_vec[j]];
                        tmp_set.set_bit(j);
                    }
                }
            }
            
            // Now tmp_set has the mask for complete IS, iterate
            // through children and find the maximum weight and then
            // substract weight from the parent before updating
            // maximum weight.

            // Iterate through children and get the maximum value.
            // if this is leaf node, will not execute the body.
            TDSolution *hash_ptr;
            
            for (kk = 0; kk < T->tree_nodes[k]->adj.size() - 1; kk++)
            {
                c_set.zeroize();
                for (j = 0; j < w; j++)
                {
                    if (cvec[kk][k_bag_vec[j]] && tmp_set.test_bit(j))
                    {
                        c_set.set_bit(j);
                    }
                }
                
                // Now we have IS which is in child hash table.

                HASH_FIND (hh, T->tree_nodes[children[kk]]->hash_table, c_set.words, 
                           T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                if (!hash_ptr)
                {
                    FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                    FERROR("error occured at child %d  kk value %d \n", children[kk], kk);
                    FERROR("mask looked for : ");
                    c_set.print(__log_file__);
                    MPI_Finalize ();
                    exit (-1);
                }
                new_value += hash_ptr->value;
            }

            if (k == T->root_node)
            {
                if (new_value > localmax)
                {
                    localmax = new_value;
                }
            }
            else
            {
                if (new_value > ind_weights[indw])
                {
                    ind_weights[indw] = new_value;
                }
            }
        }
        indw ++;
    }

    // Finish loop ind set weight calculation
    T->free_children(k);


    vector<int> ind_weights_recv(ind_size, -1);

    // Note
    // =====
    // Given that weights table become huge we should pursue other
    // methods (which is a big question right now) to distribute
    // weights to processors.  

    // Now parent collects all weights into ind_weights_recv vector,
    // since it is in order we can update before inserting into hash
    // table. 

    // ind_weights.front () is very important as it gives a pointer to
    // the first element. Spent lot of time to figure this out.

    if (k == T->root_node)
    {
        int max;
        MPI_Reduce ((void *)&localmax, (void *) &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (T->rank == 0)
        {
            CRIT("Maximum independent set value : %f \n", max);
            fprintf(stderr, "maximum independent set value : %f\n", max);
        }
    }
    else
    {
        // DEBUG("waiting on all gatherV\n");
        MPI_Allgatherv((void *)&ind_weights.front(), counts[T->rank], MPI_DOUBLE,
                       (void *)&ind_weights_recv.front(), &counts.front(), &disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);
        //DEBUG("all gather finished\n");
    }
    


    // each independent set in the intersection, therefore we can
    // directly update the hash table for the node (in the decomposed
    // tree).

    if (k != T->root_node)
    {
        TDSolution *add_set;
        bigint_t pmask(T->num_mask_words);

        //DEBUG("updating hash table\n");
        for (i = 0; i < (int)ind_sets_vec.size(); i++)
        {
            pmask.zeroize ();

            for (j = 0; j < w; j++)
            {
                if (ind_sets_vec[i]->mask->test_bit(j))
                {
                    pmask.set_bit(parent_location_vec[k_bag_vec[j]]);
                }
            }
            add_set = new TDSolution (T->num_mask_words);
            *(add_set->mask) = pmask;
            add_set->value = ind_weights_recv[i];
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, add_set->mask->words, 
                            T->num_mask_words*sizeof(BIGINT_WORD), add_set);
        }
    }
    

    return 0;
}

static int compute_nonnice_table_parallel_3 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
	list<int>::iterator ii,jj,kk;
	bigint_t current_mask(T->num_mask_words);
	vector<int> I(T->tree_nodes[k]->bag.size());
    list<int>::iterator it_li;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    vector<bigint_t *> start_masks (T->size);
    vector<bigint_t *> end_masks (T->size);

    DEBUG("mask generation\n");
    parallel_mask_gen (T, k, w, &start_masks, &end_masks, T->size);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;
	bigint_t end_mask(T->num_mask_words);

	current_mask = *(start_masks[T->rank]);
	end_mask = *(end_masks[T->rank]);

	TDSolution temp_set(T->num_mask_words);
    DEBUG("DP loop\n");

    double stime = MPI_Wtime();
	while(current_mask < end_mask)
	{
		is_independent=true;
		for(i = 0; i < w ; i++) 
		{            
			if(current_mask.test_bit(nbr_mask_vec[i]->k))
			{
				// current_mask has the bit corresponding to
				// nbr_mask_vec[i] set - so and with the appropriate
				// vector for a hit 
				has_bit=false;
				// Do the loop below manually rather than calling the overloaded AND
				for(j=0; j<current_mask.S; j++)
				{
					if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
					{
						// Some nbr is found in word j
						has_bit=true;
						break;
					}
				}
				// breaks to here
				if( has_bit )
				{                
					// Set represented by current_mask is not independent
					is_independent=false;
					int rpos=0;
					// Find the rightmost bit in current_mask - this
					// could be sped up by looking at the words 
					// SPEEDUP
					while(!(current_mask.test_bit(rpos)))
						rpos++;

					// Now advance current_mask by a potentially massive amount
					// Would be faster here to figure out which word we are in and 
					// do this w/o test_bit
					int b=0;
					while(current_mask.test_bit(rpos+b))
						b++;
					for(int a=0;a<b;a++)
						current_mask.xor_bit(rpos+a);
					current_mask.or_bit(rpos+b);
					break;                   
				}
			}
		}


		if(is_independent)
		{   
			bool is_best=false;
			int new_value=0, parent_weight=0;
            vector<int> weight = T->G->get_weight();
			for(j=0;j<w;j++)
			{                
				if(current_mask.test_bit(j))
				{
					// Compute the actual weight in G of this ind. set
					new_value += weight[k_bag_vec[j]];
				}
			}
			// Check for new max obj
			if(new_value>T->tree_nodes[k]->obj_val)
			{
				T->tree_nodes[k]->obj_val=new_value;
				is_best=true;
			}

			// Now do the DP over the children for this independent set - if this is a leaf, then
			// the body of the loop is never executed
			for(i=0;i<(int)T->tree_nodes[k]->adj.size()-1;i++)
			{
				temp_set.mask->zeroize();
				for(int z=0;z<w;z++)
				{
					// This part could be sped up by considering the words and then doing some xoring, I think
					if( cvec[i][k_bag_vec[z]] )
					{
						if(current_mask.test_bit(z))
							// Bit z is set, and the corresponding vertex is also in the child
							temp_set.mask->set_bit(z);                        
					}
				}

				HASH_FIND(hh, T->tree_nodes[children[i]]->hash_table, temp_set.mask->words,
                          T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
                
                if (!hash_ptr)
                {
                    FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                    FERROR("error occured at child %d  i value %d \n", children[i], i);
                    FERROR("mask looked for : ");
                    temp_set.mask->print(__log_file__);
                    MPI_Finalize ();
                    exit (-1);
                }
				new_value += hash_ptr->value;
				// Check for new max obj
				if(new_value>T->tree_nodes[k]->obj_val)
				{
					T->tree_nodes[k]->obj_val=new_value ;
					is_best=true;
				}
			}
			// Done with child DP

			// Now look up to the parent
			temp_set.mask->zeroize();
			parent_weight=0;
            
			for(j=0;j<w;j++)
			{                
				if(current_mask.test_bit(j))
				{       
					if(parent_location_vec[k_bag_vec[j]]!=-1)
					{
						temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
						parent_weight += weight[k_bag_vec[j]];
					}
				}
			}

			// temp_mask is a mask of this set \cap parent_bag
			hash_ptr=NULL;
			temp_set.value=new_value;
			if(k==T->root_node)
			{
				if(is_best)
				{
					TDSolution *added_set;
					added_set=new TDSolution(T->num_mask_words);
					*(added_set->mask)=*(temp_set.mask);
					*(added_set->orig_mask)=current_mask;
					added_set->value=temp_set.value;
					HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                    T->num_mask_words*sizeof(BIGINT_WORD), added_set);
				}
			}
			else
			{
				// k is not the root - so do the usual
				HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                          T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
				if(hash_ptr)
				{
					if(temp_set.value - parent_weight > hash_ptr->value)
					{
						// We found an entry for this mask and this mask has a better value 
						// so update the entry in the table
						hash_ptr->value = temp_set.value - parent_weight;
						*(hash_ptr->orig_mask)=current_mask;
                    }
				}
				else
				{
					// There is no entry for this mask - add it to the table
					TDSolution *added_set;
					added_set=new TDSolution(T->num_mask_words);
					*(added_set->mask)=*(temp_set.mask);
					*(added_set->orig_mask)=current_mask;
					added_set->value=temp_set.value - parent_weight;
					HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                    T->num_mask_words*sizeof(BIGINT_WORD), added_set);
				}
			}
			// Advance to next candidate ind. set       
			++current_mask;
		}            
	}

    T->free_children(k);
	int table_size=HASH_COUNT(T->tree_nodes[k]->hash_table);
    DEBUG("size of the table : %d\n", table_size);
    DEBUG("Time for DP loop : %lf\n", MPI_Wtime() - stime);
    hash_ptr = NULL;
    i = 0;
    TDSolution *temp_sol;
    int masksize = T->num_mask_words + 1;
    vector<unsigned long long> maskweight_array(table_size*(masksize));

    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        for (j = 0; j < T->num_mask_words; j++)
        {
            maskweight_array[i*(masksize) + j] = hash_ptr->mask->words[j];
        }
        maskweight_array[i*(masksize) + j] = (unsigned long long)ceil(hash_ptr->value);
        i++;
    }

    vector<int> n_elements_vec(T->size, 0);
    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    MPI_Allgather(&(table_size), 1, MPI_INT, &n_elements_vec.front(), 1, MPI_INT, MPI_COMM_WORLD);

    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elements_vec[i];
        w_recv_counts[i] = n_elements_vec[i];
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

    // Deriving datatype
    MPI_Datatype mask_type;
    MPI_Type_contiguous (T->num_mask_words + 1, MPI_UNSIGNED_LONG_LONG, &mask_type);
    MPI_Type_commit(&mask_type);
    vector<unsigned long long> maskweights_vec(total_entries*masksize);    

    DEBUG("starting allgather v\n");
    MPI_Allgatherv((void *)&maskweight_array.front(), w_recv_counts[T->rank], mask_type,
                   (void *)&maskweights_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), mask_type, MPI_COMM_WORLD);

    maskweight_array.clear();
    DEBUG("complete allgather v\n");

    unsigned long long max = 0;
    DEBUG("total entries : %d\n", total_entries);
    DEBUG("entering loop to update hash table \n");
    stime = MPI_Wtime ();
    // Loop to  update hash table after gathering values at other
    // processors 
    for (i = 0; i < total_entries; i++)
    {
        if (T->root_node == k)
        {
            unsigned long long cval = maskweights_vec[i*masksize + T->num_mask_words];
            if ( cval > max)
            {
                max = cval;
            }
            continue;
        }

        TDSolution *added_set;
        added_set=new TDSolution(T->num_mask_words);
        if (!added_set)
        {
            FERROR("no enough memory \n");
            MPI_Finalize();
            exit(-1);
        }

        for (j = 0; j < T->num_mask_words; j++)
        {
            added_set->mask->words[j] = maskweights_vec[i*masksize + j];
        }
        added_set->value= maskweights_vec[i*masksize + j];

        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
        }
        else
        {
            if (added_set->value > hash_ptr->value)
            {
                hash_ptr->value = added_set->value;
            }
        }
    }

    if (T->root_node == k)
    {
        DEBUG("maximum value %d\n", max);
		// CSG adding just to check answer
		fprintf(stderr,"Max value is %d\n",max);
    }
    DEBUG("Finish hash table update\n");
    DEBUG("Time for hash table update : %lf\n", MPI_Wtime() - stime);


	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	// Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}


static int compute_nonnice_table_parallel_4 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
	list<int>::iterator ii,jj,kk;
	vector<int> I(T->tree_nodes[k]->bag.size());
    list<int>::iterator it_li;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    //vector<bigint_t *> start_masks (entries);
    //vector<bigint_t *> end_masks (entries);
 
    DEBUG("mask generation\n");
    parallel_mask_gen (T, k, w, T->start_masks, T->end_masks, T->entries);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;

	bigint_t current_mask(T->num_mask_words);
	bigint_t end_mask(T->num_mask_words);

	TDSolution temp_set(T->num_mask_words);

    double stime = MPI_Wtime();
    int request_flag = 0;
    int request_index = -1;

    unsigned long long start_pos, end_pos, wchunk, counter, stwork;
    wchunk = T->entries/T->size;
    start_pos = (wchunk * T->rank);
    end_pos = (wchunk * (T->rank + 1));

    vector<int> requests(T->size, 0);
    DEBUG("DP loop\n");
    counter = start_pos;

    while (1)
    {
        
        // DEBUG("counter %d end_pos %d\n", counter, end_pos);
        if (counter < end_pos)
        {
            current_mask = *(*T->start_masks)[counter];
            end_mask = *(*T->end_masks)[counter];
            counter ++;
        }

        while(current_mask < end_mask)
        {

            if ((T->rank == 0) && (counter < end_pos))
            {

                MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                while (request_flag)
                {
                    MPI_Irecv (&(T->requester_rank), 1, MPI_INT, request_index, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, &((*T->requests)[request_index]));
                    requests[request_index] = 1;
                    
                    request_flag = 0;
                    request_index = -1;
                    T->requester_rank = -1;
                    MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                }
                
                    // Send work to the requester.
                for (i = 1; i < T->size && (end_pos > counter); i++)
                {
                    if (requests[i])
                    {
                        end_pos --;
                        MPI_Send (&end_pos, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                        requests[i] = 0;
                    }
                }
            }
            
            is_independent=true;

            for(i = 0; i < w ; i++) 
            {            
                if(current_mask.test_bit(nbr_mask_vec[i]->k))
                {
                    // current_mask has the bit corresponding to
                    // nbr_mask_vec[i] set - so and with the appropriate
                    // vector for a hit 
                    has_bit=false;
                    // Do the loop below manually rather than calling the overloaded AND
                    for(j=0; j<current_mask.S; j++)
                    {
                        if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                        {
                            // Some nbr is found in word j
                            has_bit=true;
                            break;
                        }
                    }
                    // breaks to here
                    if( has_bit )
                    {                
                        // Set represented by current_mask is not independent
                        is_independent=false;
                        int rpos=0;
                        // Find the rightmost bit in current_mask - this
                        // could be sped up by looking at the words 
                        // SPEEDUP
                        while(!(current_mask.test_bit(rpos)))
                            rpos++;

                        // Now advance current_mask by a potentially massive amount
                        // Would be faster here to figure out which word we are in and 
                        // do this w/o test_bit
                        int b=0;
                        while(current_mask.test_bit(rpos+b))
                            b++;
                        for(int a=0;a<b;a++)
                            current_mask.xor_bit(rpos+a);
                        current_mask.or_bit(rpos+b);
                        break;                   
                    }
                }
            }


            if(is_independent)
            {   
                T->info->total_table_entries++;
                bool is_best=false;
                int new_value=0, parent_weight=0;
                vector<int> weight = T->G->get_weight();

                for(j=0;j<w;j++)
                {                
                    if(current_mask.test_bit(j))
                    {
                        // Compute the actual weight in G of this ind. set
                        new_value += weight[k_bag_vec[j]];
                    }
                }
                // Check for new max obj
                if(new_value>T->tree_nodes[k]->obj_val)
                {
                    T->tree_nodes[k]->obj_val=new_value;
                    is_best=true;
                }

                // Now do the DP over the children for this independent set - if this is a leaf, then
                // the body of the loop is never executed
                for(i=0;i<(int)T->tree_nodes[k]->adj.size()-1;i++)
                {
                    temp_set.mask->zeroize();
                    for(int z=0;z<w;z++)
                    {
                        // This part could be sped up by considering the words and then doing some xoring, I think
                        if( cvec[i][k_bag_vec[z]] )
                        {
                            if(current_mask.test_bit(z))
                                // Bit z is set, and the corresponding vertex is also in the child
                                temp_set.mask->set_bit(z);                        
                        }
                    }

                    HASH_FIND(hh, T->tree_nodes[children[i]]->hash_table, temp_set.mask->words,
                              T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
                
                    if (!hash_ptr)
                    {
                        FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                        FERROR("error occured at child %d  i value %d \n", children[i], i);
                        FERROR("mask looked for : ");
                        temp_set.mask->print(__log_file__);
                        MPI_Finalize ();
                        exit (-1);
                    }
                    new_value += hash_ptr->value;
                    // Check for new max obj
                    if(new_value>T->tree_nodes[k]->obj_val)
                    {
                        T->tree_nodes[k]->obj_val=new_value ;
                        is_best=true;
                    }
                }
                // Done with child DP

                // Now look up to the parent
                temp_set.mask->zeroize();
                parent_weight=0;
            
                for(j=0;j<w;j++)
                {                
                    if(current_mask.test_bit(j))
                    {       
                        if(parent_location_vec[k_bag_vec[j]]!=-1)
                        {
                            temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
                            parent_weight += weight[k_bag_vec[j]];
                        }
                    }
                }

                // temp_mask is a mask of this set \cap parent_bag
                hash_ptr=NULL;
                temp_set.value=new_value;
                if(k==T->root_node)
                {
                    if(is_best)
                    {
                        TDSolution *added_set;
                        added_set=new TDSolution(T->num_mask_words);
                        *(added_set->mask)=*(temp_set.mask);
                        *(added_set->orig_mask)=current_mask;
                        added_set->value=temp_set.value;
                        HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                        T->num_mask_words*sizeof(BIGINT_WORD), added_set);
                    }
                }
                else
                {
                    // k is not the root - so do the usual
                    HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                              T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
                    if(hash_ptr)
                    {
                        if(temp_set.value - parent_weight > hash_ptr->value)
                        {
                            // We found an entry for this mask and this mask has a better value 
                            // so update the entry in the table
                            hash_ptr->value = temp_set.value - parent_weight;
                            *(hash_ptr->orig_mask)=current_mask;
                        }
                    }
                    else
                    {
                        // There is no entry for this mask - add it to the table
                        TDSolution *added_set;
                        added_set=new TDSolution(T->num_mask_words);
                        *(added_set->mask)=*(temp_set.mask);
                        *(added_set->orig_mask)=current_mask;
                        added_set->value=temp_set.value - parent_weight;
                        HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                        T->num_mask_words*sizeof(BIGINT_WORD), added_set);
                    }
                }
                // Advance to next candidate ind. set       
                ++current_mask;
            }
        }

        if (!(counter < end_pos))
        {
            if (T->rank > 0)
            {
                MPI_Request request;
                MPI_Status status;
                int wflag = 0;
                int count = 0;
            
                //DEBUG("process : %d ask for work from : %d\n", T->rank, 0);
                MPI_Send(&(T->rank), 1, MPI_INT, 0, MPI_REQ_WRK_TAG, MPI_COMM_WORLD);
                MPI_Recv (&stwork, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_WRK_TAG, MPI_COMM_WORLD, &status);

                count = 0;
                MPI_Get_count (&status, MPI_UNSIGNED_LONG_LONG, &count);

                if (!(count > 0))
                    break;

                current_mask.zeroize ();
                end_mask.zeroize ();
                current_mask = *(*T->start_masks)[stwork];
                end_mask = *(*T->end_masks)[stwork];
            }
            else
            {
                //Sending terminate signal
                for (i = 0; i < T->size; i++)
                    MPI_Send (&end_pos, 0, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                break;
            }
        }
    }
    
    T->free_children(k);
	int table_size=HASH_COUNT(T->tree_nodes[k]->hash_table);
    DEBUG("size of the table : %d\n", table_size);
    DEBUG("Time for DP loop : %lf\n", MPI_Wtime() - stime);

    hash_ptr = NULL;
    i = 0;
    TDSolution *temp_sol;
    int masksize = T->num_mask_words + 1;
    vector<unsigned long long> maskweight_array(table_size*(masksize));

    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        for (j = 0; j < T->num_mask_words; j++)
        {
            maskweight_array[i*(masksize) + j] = hash_ptr->mask->words[j];
        }
        maskweight_array[i*(masksize) + j] = hash_ptr->value;
        i++;
    }

    vector<int> n_elements_vec(T->size, 0);
    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    MPI_Allgather(&(table_size), 1, MPI_INT, &n_elements_vec.front(), 1, MPI_INT, MPI_COMM_WORLD);

    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elements_vec[i];
        w_recv_counts[i] = n_elements_vec[i];
        // DEBUG("node %d: work %d\n", i, n_elements_vec[i]);
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

    // Deriving datatype
    MPI_Datatype mask_type;
    MPI_Type_contiguous (T->num_mask_words + 1, MPI_UNSIGNED_LONG_LONG, &mask_type);
    MPI_Type_commit(&mask_type);
    vector<unsigned long long> maskweights_vec(total_entries*masksize);    


    DEBUG("starting allgather v\n");
    MPI_Allgatherv((void *)&maskweight_array.front(), w_recv_counts[T->rank], mask_type,
                   (void *)&maskweights_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), mask_type, MPI_COMM_WORLD);
    DEBUG("complete allgather v\n");

    unsigned long long max = 0;
    DEBUG("total entries : %d\n", total_entries);
    DEBUG("entering loop to update hash table \n");

    stime = MPI_Wtime ();
    // Loop to  update hash table after gathering values at other
    // processors 
    for (i = 0; i < total_entries; i++)
    {
        if (T->root_node == k)
        {
            unsigned long long cval = maskweights_vec[i*masksize + T->num_mask_words];
            if ( cval > max)
            {
                max = cval;
                T->info->opt_obj = max;
            }
            continue;
        }

        TDSolution *added_set;
        added_set=new TDSolution(T->num_mask_words);
        if (!added_set)
        {
            FERROR("no enough memory \n");
            MPI_Finalize();
            exit(-1);
        }

        for (j = 0; j < T->num_mask_words; j++)
        {
            added_set->mask->words[j] = maskweights_vec[i*masksize + j];
        }
        added_set->value= maskweights_vec[i*masksize + j];

        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
        }
        else
        {
            if (added_set->value > hash_ptr->value)
            {
                hash_ptr->value = added_set->value;
            }
        }
    }

	T->info->total_pc_table_entries += (unsigned long long)HASH_COUNT(T->tree_nodes[k]->hash_table);
    DEBUG("completed table size : %d\n", (int)HASH_COUNT(T->tree_nodes[k]->hash_table));
    // if (T->root_node == k)
    // {
    //     DEBUG("maximum value %d\n", max);
	// 	// CSG adding just to check answer
    //     if (T->rank == 0)
    //         fprintf(stderr,"Max value is %d\n",max);
    // }
    DEBUG("Finish hash table update\n");
    DEBUG("Time for hash table update : %lf\n", MPI_Wtime() - stime);

    // clear some vectors to gain more memory.
    maskweights_vec.clear();
    maskweight_array.clear();
    w_recv_disp.clear();
    w_recv_counts.clear();

	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	// Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}

static int compute_nonnice_table_parallel_5 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
	list<int>::iterator ii,jj,kk;
	vector<int> I(T->tree_nodes[k]->bag.size());
    list<int>::iterator it_li;

    int entries = 32768;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    vector<bigint_t *> start_masks (entries);
    vector<bigint_t *> end_masks (entries);

    DEBUG("mask generation\n");
    parallel_mask_gen (T, k, w, &start_masks, &end_masks, entries);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;

	bigint_t current_mask(T->num_mask_words);
	bigint_t end_mask(T->num_mask_words);

	TDSolution temp_set(T->num_mask_words);

    double stime = MPI_Wtime();
    int request_flag = 0;
    int request_index = -1;

    unsigned long long start_pos, end_pos, wchunk, counter, stwork;
    wchunk = entries/T->size;
    start_pos = (wchunk * T->rank);
    end_pos = (wchunk * (T->rank + 1));

    vector<int> requests(T->size, 0);
    DEBUG("DP loop\n");
    counter = start_pos;

    while (1)
    {
        
        // DEBUG("counter %d end_pos %d\n", counter, end_pos);
        if (counter < end_pos)
        {
            current_mask = *(start_masks[counter]);
            end_mask = *(end_masks[counter]);
            counter ++;
        }

        while(current_mask < end_mask)
        {

            if ((T->rank == 0) && (counter < end_pos))
            {

                MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                while (request_flag)
                {
                    MPI_Irecv (&(T->requester_rank), 1, MPI_INT, request_index, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, &((*T->requests)[request_index]));
                    requests[request_index] = 1;
                    
                    request_flag = 0;
                    request_index = -1;
                    T->requester_rank = -1;
                    MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                }
                
                    // Send work to the requester.
                for (i = 1; i < T->size && (end_pos > counter); i++)
                {
                    if (requests[i])
                    {
                        end_pos --;
                        MPI_Send (&end_pos, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                        requests[i] = 0;
                    }
                }
            }
            
            is_independent=true;

            for(i = 0; i < w ; i++) 
            {            
                if(current_mask.test_bit(nbr_mask_vec[i]->k))
                {
                    // current_mask has the bit corresponding to
                    // nbr_mask_vec[i] set - so and with the appropriate
                    // vector for a hit 
                    has_bit=false;
                    // Do the loop below manually rather than calling the overloaded AND
                    for(j=0; j<current_mask.S; j++)
                    {
                        if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                        {
                            // Some nbr is found in word j
                            has_bit=true;
                            break;
                        }
                    }
                    // breaks to here
                    if( has_bit )
                    {                
                        // Set represented by current_mask is not independent
                        is_independent=false;
                        int rpos=0;
                        // Find the rightmost bit in current_mask - this
                        // could be sped up by looking at the words 
                        // SPEEDUP
                        while(!(current_mask.test_bit(rpos)))
                            rpos++;

                        // Now advance current_mask by a potentially massive amount
                        // Would be faster here to figure out which word we are in and 
                        // do this w/o test_bit
                        int b=0;
                        while(current_mask.test_bit(rpos+b))
                            b++;
                        for(int a=0;a<b;a++)
                            current_mask.xor_bit(rpos+a);
                        current_mask.or_bit(rpos+b);
                        break;                   
                    }
                }
            }


            if(is_independent)
            {   
                bool is_best=false;
                int new_value=0, parent_weight=0;
                vector<int> weight = T->G->get_weight();
                for(j=0;j<w;j++)
                {                
                    if(current_mask.test_bit(j))
                    {
                        // Compute the actual weight in G of this ind. set
                        new_value += weight[k_bag_vec[j]];
                    }
                }
                // Check for new max obj
                if(new_value>T->tree_nodes[k]->obj_val)
                {
                    T->tree_nodes[k]->obj_val=new_value;
                    is_best=true;
                }

                // Now do the DP over the children for this independent set - if this is a leaf, then
                // the body of the loop is never executed
                for(i=0;i<(int)T->tree_nodes[k]->adj.size()-1;i++)
                {
                    temp_set.mask->zeroize();
                    for(int z=0;z<w;z++)
                    {
                        // This part could be sped up by considering the words and then doing some xoring, I think
                        if( cvec[i][k_bag_vec[z]] )
                        {
                            if(current_mask.test_bit(z))
                                // Bit z is set, and the corresponding vertex is also in the child
                                temp_set.mask->set_bit(z);                        
                        }
                    }

                    HASH_FIND(hh, T->tree_nodes[children[i]]->hash_table, temp_set.mask->words,
                              T->num_mask_words*sizeof(BIGINT_WORD),hash_ptr);
                
                    if (!hash_ptr)
                    {
                        FERROR("%s : Set not available in the child, can not proceed\n", __FUNCTION__);
                        FERROR("error occured at child %d  i value %d \n", children[i], i);
                        FERROR("mask looked for : ");
                        temp_set.mask->print(__log_file__);
                        MPI_Finalize ();
                        exit (-1);
                    }
                    new_value += hash_ptr->value;
                    // Check for new max obj
                    if(new_value>T->tree_nodes[k]->obj_val)
                    {
                        T->tree_nodes[k]->obj_val=new_value ;
                        is_best=true;
                    }
                }
                // Done with child DP

                // Now look up to the parent
                temp_set.mask->zeroize();
                parent_weight=0;
            
                for(j=0;j<w;j++)
                {                
                    if(current_mask.test_bit(j))
                    {       
                        if(parent_location_vec[k_bag_vec[j]]!=-1)
                        {
                            temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
                            parent_weight += weight[k_bag_vec[j]];
                        }
                    }
                }

                // temp_mask is a mask of this set \cap parent_bag
                hash_ptr=NULL;
                temp_set.value=new_value;
                if(k==T->root_node)
                {
                    if(is_best)
                    {
                        TDSolution *added_set;
                        added_set=new TDSolution(T->num_mask_words);
                        *(added_set->mask)=*(temp_set.mask);
                        *(added_set->orig_mask)=current_mask;
                        added_set->value=temp_set.value;
                        HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                        T->num_mask_words*sizeof(BIGINT_WORD), added_set);
                    }
                }
                else
                {
                    // k is not the root - so do the usual
                    HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                              T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);
                    if(hash_ptr)
                    {
                        if(temp_set.value - parent_weight > hash_ptr->value)
                        {
                            // We found an entry for this mask and this mask has a better value 
                            // so update the entry in the table
                            hash_ptr->value = temp_set.value - parent_weight;
                            *(hash_ptr->orig_mask)=current_mask;
                        }
                    }
                    else
                    {
                        // There is no entry for this mask - add it to the table
                        TDSolution *added_set;
                        added_set=new TDSolution(T->num_mask_words);
                        *(added_set->mask)=*(temp_set.mask);
                        *(added_set->orig_mask)=current_mask;
                        added_set->value=temp_set.value - parent_weight;
                        HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                        T->num_mask_words*sizeof(BIGINT_WORD), added_set);
                    }
                }
                // Advance to next candidate ind. set       
                ++current_mask;
            }
        }

        if (!(counter < end_pos))
        {
            if (T->rank > 0)
            {
                MPI_Request request;
                MPI_Status status;
                int wflag = 0;
                int count = 0;
            
                //DEBUG("process : %d ask for work from : %d\n", T->rank, 0);
                MPI_Send(&(T->rank), 1, MPI_INT, 0, MPI_REQ_WRK_TAG, MPI_COMM_WORLD);
                MPI_Recv (&stwork, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_WRK_TAG, MPI_COMM_WORLD, &status);

                count = 0;
                MPI_Get_count (&status, MPI_UNSIGNED_LONG_LONG, &count);

                if (!(count > 0))
                    break;

                current_mask.zeroize ();
                end_mask.zeroize ();
                current_mask = *(start_masks[stwork]);
                end_mask = *(end_masks[stwork]);
            }
            else
            {
                //Sending terminate signal
                for (i = 0; i < T->size; i++)
                    MPI_Send (&end_pos, 0, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                break;
            }
        }
    }
    
    T->free_children(k);
	int table_size=HASH_COUNT(T->tree_nodes[k]->hash_table);
    DEBUG("size of the table : %d\n", table_size);
    DEBUG("Time for DP loop : %lf\n", MPI_Wtime() - stime);

    hash_ptr = NULL;
    i = 0;
    TDSolution *temp_sol;
    vector<unsigned long long> masks_array (table_size*T->num_mask_words);
    vector<int> weights_array(table_size);
  
    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        for (j = 0; j < T->num_mask_words; j++)
            masks_array[i*(T->num_mask_words) + j] = hash_ptr->mask->words[j];
        weights_array[i] = hash_ptr->value;
        i++;
    }


    vector<int> n_elements_vec(T->size, 0);
    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    MPI_Allgather(&(table_size), 1, MPI_INT, &n_elements_vec.front(), 1, MPI_INT, MPI_COMM_WORLD);

    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elements_vec[i];
        w_recv_counts[i] = n_elements_vec[i];
        // DEBUG("node %d: work %d\n", i, n_elements_vec[i]);
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

    // Deriving datatype
    MPI_Datatype mask_type;
    MPI_Type_contiguous (T->num_mask_words, MPI_UNSIGNED_LONG_LONG, &mask_type);
    MPI_Type_commit(&mask_type);

    vector<unsigned long long> masks_vec(total_entries*T->num_mask_words);    
    vector<int> weights_vec(total_entries);

    DEBUG("starting allgather v\n");

    MPI_Allgatherv((void *)&masks_array.front(), w_recv_counts[T->rank], mask_type,
                   (void *)&masks_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), mask_type, MPI_COMM_WORLD);

    MPI_Allgatherv((void *)&weights_array.front(), w_recv_counts[T->rank], MPI_DOUBLE,
                   (void *)&weights_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), MPI_DOUBLE, MPI_COMM_WORLD);

    int max = 0;
    DEBUG("total entries : %d\n", total_entries);
    DEBUG("entering loop to update hash table \n");

    stime = MPI_Wtime ();
    // Loop to  update hash table after gathering values at other
    // processors 
    for (i = 0; i < total_entries; i++)
    {
        if (T->root_node == k)
        {
            int cval = weights_vec[i];
            if ( cval > max)
            {
                max = cval;
            }
            continue;
        }

        TDSolution *added_set;
        added_set=new TDSolution(T->num_mask_words);
        if (!added_set)
        {
            FERROR("no enough memory \n");
            MPI_Finalize();
            exit(-1);
        }

        for (j = 0; j < T->num_mask_words; j++)
            added_set->mask->words[j] = masks_vec[i*T->num_mask_words + j];
        added_set->value = weights_vec[i];

        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
        }
        else
        {
            if (added_set->value > hash_ptr->value)
            {
                hash_ptr->value = added_set->value;
            }
        }
    }

    DEBUG("completed table size : %d\n", (int)HASH_COUNT(T->tree_nodes[k]->hash_table));
    if (T->root_node == k)
    {
        DEBUG("maximum value %f\n", max);
		// CSG adding just to check answer
        if (T->rank == 0)
            fprintf(stderr,"Max value is %f\n",max);
    }
    DEBUG("Finish hash table update\n");
    DEBUG("Time for hash table update : %lf\n", MPI_Wtime() - stime);

    // clear some vectors to gain more memory.
    weights_vec.clear();
    masks_vec.clear();
    w_recv_disp.clear();
    w_recv_counts.clear();
    weights_array.clear();
    masks_array.clear();

	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	// Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}

void parallel_wis_init (TDTree *T, int size, int rank)
{
    int i = 0;

    T->size = size;
    T->rank = rank;
    T->tnum = 1;
    T->request_size = 2048;
    T->pool_size = 2048;
    T->flag_recv_done = 0;
    T->msk_send_count = 0;
    T->wgh_send_count = 0;

    if (!(T->size > 0))
    {
        FERROR("Environment size can not be zero\n");
        MPI_Finalize();
        exit (-1);
    }

    if (T->rank < 0)
    {
        FERROR("Ranks can not be negative\n");
        MPI_Finalize();
        exit (-1);
    }

    try 
    {
        T->tbl_entries = new vector<int> (T->size,0);
        T->tbl_loc = new vector<int> (T->num_tree_nodes, -1);
        T->node_sizes = new vector<int> (T->num_tree_nodes, 0);
        T->iter_count = new vector<int> (T->size, 0);
        T->threads = new vector<pthread_t>(T->tnum);

        T->requests = new vector<MPI_Request>(T->size);
        if (T->rank == 0)
        {
            for (i = 0; i < T->size; i++)
            {
                MPI_Irecv (&(T->requester_rank), 1, MPI_INT, i, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, &((*T->requests)[i]));
            }
        }

        // request_count vector keep track of # of requests each node has
        // to receive.
        //T->request_count = new vector<int>(T->size, -1);

        // Allocating buffers to receive weight requests. Masks are sent
        // over the network along with the child number and iteration #.
        T->msk_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->msk_recv_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            MPI_Irecv((void *)&(*T->msk_recv_pool[i]).front(), T->request_size*T->num_mask_words + 2, 
                      MPI_UNSIGNED_LONG_LONG, i, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_recv_requests)[i]);
        }


        // Response with weights for masks comes with one additional
        // element, iteration #
        // Iteration # is used to find the correct position in original
        // vector for the weights received. 
        T->wgh_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->wgh_recv_pool.push_back(new vector<int>(T->request_size + 1, 0));
            MPI_Irecv ((void *)&T->wgh_recv_pool[i]->front(), T->request_size + 1, 
                       MPI_INT, i, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[i]);
        }

        T->msk_send_requests = new vector<MPI_Request>(T->pool_size);
        T->wgh_send_requests = new vector<MPI_Request>(T->pool_size);
        for (i = 0; i < T->pool_size; i++)
        {
            T->msk_send_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            T->wgh_send_pool.push_back(new vector<int>(T->request_size + 1, 0));
        }

        DEBUG("Creating threads : %d\n", T->tnum);
        for (i = 0; i < T->tnum; i++)
        {
            if(pthread_create(&(*T->threads)[i], NULL, thread_recv_msk_send_wgh, (void *)T))
            {
                FERROR("Error while creating threads, processor: %d Abort\n", rank);
                MPI_Finalize();
                exit (-1);
            }
        }
    }
    catch (std::bad_alloc e)
    {
        cerr << e.what() << endl;
        FERROR("Bad allocation: %s\n", e.what());
        MPI_Finalize();
        exit (-1);
    }
}

void parallel_wis_init_non_thread (TDTree *T, int size, int rank)
{
    int i = 0;

    T->size = size;
    T->rank = rank;
    T->request_size = 2048;
    T->pool_size = 2048;
    T->msk_send_count = 0;
    T->wgh_send_count = 0;
    T->entries = 65536;

    if (!(T->size > 0))
    {
        FERROR("Environment size can not be zero\n");
        MPI_Finalize();
        exit (-1);
    }

    if (T->rank < 0)
    {
        FERROR("Ranks can not be negative\n");
        MPI_Finalize();
        exit (-1);
    }

    try 
    {

        T->start_masks = new vector<bigint_t *> (T->entries);
        T->end_masks = new vector<bigint_t *> (T->entries);
        for (i = 0; i < T->entries; i++)
        {
            (*T->start_masks)[i] = new bigint_t (T->num_mask_words);
            (*T->end_masks)[i] = new bigint_t (T->num_mask_words);
        }


        T->tbl_entries = new vector<int> (T->size,0);
        T->tbl_loc = new vector<int> (T->num_tree_nodes, -1);
        T->node_sizes = new vector<int> (T->num_tree_nodes, 0);
        T->iter_count = new vector<int> (T->size, 0);
        T->rcount = new vector<int> (T->size, 0);
        // request_count vector keep track of # of requests each node has
        // to receive.
        T->request_count = new vector<int>(T->size, -1);

    
        T->requests = new vector<MPI_Request>(T->size);
        if (T->rank == 0)
        {
            for (i = 0; i < T->size; i++)
            {
                MPI_Irecv (&(T->requester_rank), 1, MPI_INT, i, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, &((*T->requests)[i]));
            }
        }

        // Allocating buffers to receive weight requests. Masks are sent
        // over the network along with the child number and iteration #.
        T->msk_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->msk_recv_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            MPI_Irecv((void *)&(*T->msk_recv_pool[i]).front(), T->request_size*T->num_mask_words + 2, 
                      MPI_UNSIGNED_LONG_LONG, i, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_recv_requests)[i]);
        }


        // Response with weights for masks comes with one additional
        // element, iteration #
        // Iteration # is used to find the correct position in original
        // vector for the weights received. 
        T->wgh_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->wgh_recv_pool.push_back(new vector<int>(T->request_size + 1, 0));
            MPI_Irecv ((void *)&T->wgh_recv_pool[i]->front(), T->request_size + 1, 
                       MPI_INT, i, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[i]);
        }

        T->msk_send_requests = new vector<MPI_Request>(T->pool_size);
        T->wgh_send_requests = new vector<MPI_Request>(T->pool_size);
        for (i = 0; i < T->pool_size; i++)
        {
            T->msk_send_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            T->wgh_send_pool.push_back(new vector<int>(T->request_size + 1, 0));
        }
    }
    catch (std::bad_alloc e)
    {
        cerr << e.what() << endl;
        FERROR("Bad allocation: %s\n", e.what());
        MPI_Finalize();
        exit (-1);
    }
}


void parallel_wis_cleanup_non_thread (TDTree *T)
{
    int i;

     // joining threads
    DEBUG("join threads\n");

    // Cancel waiting requests.
    for (i = 0; i < T->size; i++)
    {
        if (T->rank == 0)
        {
            MPI_Cancel (&(*T->requests)[i]);
        }
        MPI_Cancel(&(*T->wgh_recv_requests)[i]);
        MPI_Cancel(&(*T->msk_recv_requests)[i]);
    }

    delete T->tbl_entries;
    delete T->tbl_loc;
    delete T->node_sizes;
    delete T->iter_count;
    delete T->requests;
    delete T->start_masks;
    delete T->end_masks;
//    delete T->request_count;
    //  delete T->rcount;

    delete T->msk_recv_requests;
    for (i = 0; i < T->size; i++)
    {
        delete T->msk_recv_pool[i];
        delete T->wgh_recv_pool[i];
        
    }
    delete T->wgh_recv_requests;

    delete T->msk_send_requests;
    delete T->wgh_send_requests;

    for (i = 0; i < T->pool_size; i++)
    {
        delete  T->msk_send_pool[i];
        delete  T->wgh_send_pool[i];
    }
}

void parallel_wis_cleanup (TDTree *T)
{
    int i;


     // joining threads
    DEBUG("join threads\n");

    // Cancel waiting requests.
    for (i = 0; i < T->size; i++)
    {
        if (T->rank == 0)
        {
            MPI_Cancel (&(*T->requests)[i]);
        }
        MPI_Cancel(&(*T->wgh_recv_requests)[i]);
        MPI_Cancel(&(*T->msk_recv_requests)[i]);
    }

    for (i = 0; i < T->tnum; i++)
        pthread_join((*T->threads)[i], NULL);

    delete T->tbl_entries;
    delete T->tbl_loc;
    delete T->node_sizes;
    delete T->iter_count;
    delete T->requests;
    delete T->request_count;

    delete T->msk_recv_requests;
    for (i = 0; i < T->size; i++)
    {
        delete T->msk_recv_pool[i];
        delete T->wgh_recv_pool[i];
        
    }
    delete T->wgh_recv_requests;

    delete T->msk_send_requests;
    delete T->wgh_send_requests;

    for (i = 0; i < T->pool_size; i++)
    {
        delete  T->msk_send_pool[i];
        delete  T->wgh_send_pool[i];
    }
}

static int get_msk_send_buffer_index(TDTree *T)
{
    int msk_send_index = 0;
    MPI_Status msk_send_status;
    if (!(T->msk_send_count < T->pool_size))
    {
        MPI_Waitany(T->pool_size, &T->msk_send_requests->front(),
        &msk_send_index, &msk_send_status); 
        return msk_send_index;
    }
    return T->msk_send_count++;
}

static int get_wgh_send_buffer_index(TDTree *T)
{
//     int wgh_send_index = 0;
//     MPI_Status wgh_send_status;
//     if (!(T->wgh_send_count < T->pool_size))
//     {
//         MPI_Waitany(T->pool_size, &T->wgh_send_requests->front(), &wgh_send_index, &wgh_send_status);
//         return wgh_send_index;
//     }
//     return T->wgh_send_count++;

    int wgh_send_index = 0;
    int wgh_send_flag = 0;
    MPI_Status wgh_send_status;
    if (!(T->wgh_send_count < T->pool_size))
    {
        MPI_Testany(T->pool_size, &T->wgh_send_requests->front(), &wgh_send_index, &wgh_send_flag, &wgh_send_status);
        if (wgh_send_flag)
            return wgh_send_index;
        else
            return -1;
    }
    return T->wgh_send_count++;

}


static int compute_nonnice_table_parallel_6 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
	vector<int> I(T->tree_nodes[k]->bag.size());
    list<int>::iterator it_li;
    int kk;
    int entries = 32768;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());
    int children_size = (int) children.size();

    int child_tbl_flag = 0;
    // If a child table is stored in current node, child_tbl_flag will
    // be set to 1.  Which means current node should be ready to serve
    // weight requests. 
    for (i = 0; i < children_size; i++)
    {
        if (T->rank == (*T->tbl_loc)[children[i]])
        {
            child_tbl_flag = 1;
            break;
        }
    }


    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    vector<bigint_t *> start_masks (entries);
    vector<bigint_t *> end_masks (entries);

    DEBUG("mask generation\n");
    parallel_mask_gen (T, k, w, &start_masks, &end_masks, entries);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;

    vector<bigint_t> indset_masks;
    vector<int> ind_wgh;
    unsigned long long ind_counter = 0;
    
 
	bigint_t current_mask(T->num_mask_words);
	bigint_t end_mask(T->num_mask_words);

	TDSolution temp_set(T->num_mask_words);

    double stime = MPI_Wtime();
    int request_flag = 0;
    int request_index = -1;

    unsigned long long start_pos, end_pos, wchunk, counter, stwork;
    wchunk = entries/T->size;
    start_pos = (wchunk * T->rank);
    end_pos = (wchunk * (T->rank + 1));

    vector<int> requests(T->size, 0);
    DEBUG("DP loop\n");
    counter = start_pos;
    while (1)
    {
        DEBUG("counter %d end_pos %d\n", counter, end_pos);
        if (counter < end_pos)
        {
            current_mask = *(start_masks[counter]);
            end_mask = *(end_masks[counter]);
            counter ++;
        }

        while(current_mask < end_mask)
        {

            if ((T->rank == 0) && (counter < end_pos))
            {

                MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                while (request_flag)
                {
                    MPI_Irecv (&(T->requester_rank), 1, MPI_INT, request_index, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, 
                               &((*T->requests)[request_index]));
                    requests[request_index] = 1;
                    
                    request_flag = 0;
                    request_index = -1;
                    T->requester_rank = -1;
                    MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                }
                
                // Send work to the requester.
                for (i = 1; i < T->size && (end_pos > counter); i++)
                {
                    if (requests[i])
                    {
                        end_pos --;
                        MPI_Send (&end_pos, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                        requests[i] = 0;
                    }
                }
            }
            
            is_independent=true;

            for(i = 0; i < w ; i++) 
            {            
                if(current_mask.test_bit(nbr_mask_vec[i]->k))
                {
                    // current_mask has the bit corresponding to
                    // nbr_mask_vec[i] set - so and with the appropriate
                    // vector for a hit 
                    has_bit=false;
                    // Do the loop below manually rather than calling the overloaded AND
                    for(j=0; j<current_mask.S; j++)
                    {
                        if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                        {
                            // Some nbr is found in word j
                            has_bit=true;
                            break;
                        }
                    }
                    // breaks to here
                    if( has_bit )
                    {                
                        // Set represented by current_mask is not independent
                        is_independent=false;
                        int rpos=0;
                        // Find the rightmost bit in current_mask - this
                        // could be sped up by looking at the words 
                        // SPEEDUP
                        while(!(current_mask.test_bit(rpos)))
                            rpos++;

                        // Now advance current_mask by a potentially massive amount
                        // Would be faster here to figure out which word we are in and 
                        // do this w/o test_bit
                        int b=0;
                        while(current_mask.test_bit(rpos+b))
                            b++;
                        for(int a=0;a<b;a++)
                            current_mask.xor_bit(rpos+a);
                        current_mask.or_bit(rpos+b);
                        break;                   
                    }
                }
            }


            if(is_independent)
            {   
                // ind_counter keep track of how many independent sets
                // are in the indset_mask vector
                ind_counter ++;
                indset_masks.push_back(current_mask);
                ind_wgh.push_back(0);
                ++current_mask;
            }
        }


        // Once assigned work is finished, ask for more work
        if (!(counter < end_pos))
        {
            if (T->rank > 0)
            {
                MPI_Request request;
                MPI_Status status;
                int wflag = 0;
                int count = 0;
            
                //DEBUG("process : %d ask for work from : %d\n", T->rank, 0);
                MPI_Send(&(T->rank), 1, MPI_INT, 0, MPI_REQ_WRK_TAG, MPI_COMM_WORLD);
                MPI_Recv (&stwork, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_WRK_TAG, MPI_COMM_WORLD, &status);

                count = 0;
                MPI_Get_count (&status, MPI_UNSIGNED_LONG_LONG, &count);

                if (!(count > 0))
                    break;

                current_mask.zeroize ();
                end_mask.zeroize ();
                current_mask = *(start_masks[stwork]);
                end_mask = *(end_masks[stwork]);
            }
            else
            {
                //Sending terminate signal
                for (i = 0; i < T->size; i++)
                    MPI_Send (&end_pos, 0, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                break;
            }
        }
    }

    // Now do the DP over the children for this independent set - if this is a leaf, then
    // the body of the loop is never executed
    int ind_set_tbl_size = (int)indset_masks.size();
    MPI_Request wgh_s_request;
    MPI_Request wgh_r_request;
    int imaskcount = 0;
    int request_size = T->request_size;
    MPI_Status wstatus[T->size];
    request_index = 0;
    request_flag = 0;
    //vector<vector<unsigned long long> *>chd_masks_recv;
    //vector<vector<unsigned long long> *> chd_masks_array(children_size);
    //vector<int> request_count(T->size, 0);
    vector<int> *wgh_response;

    

    DEBUG("computing mask weights\n");
    int iteration = 0;
    DEBUG("mask count %d tble size %d\n", imaskcount, ind_set_tbl_size);
    if (!(ind_set_tbl_size > 0))
    {
        for(i = 0; i < (int)children_size; i++)
            MPI_Isend ((void *)NULL, 0, MPI_UNSIGNED_LONG_LONG, 
                       (*T->tbl_loc)[children[i]], MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &wgh_s_request);
    }

    //DEBUG("sending out requests\n");
    DEBUG("# of nbrs %d\n", children_size);
    DEBUG("child:location: ");
    for ( i = 0; i < children_size; i++)
    {
        GEN("%d:%d ", children[i], (*T->tbl_loc)[children[i]]);
    }
    GEN("\n");

    // Sending request of size request_size.
    int count_send = 0;
    int destination = 0;
    MPI_Request send_request[children_size];
    int msk_send_index = 0;
    //MPI_Request msk_send_request;
    vector<unsigned long long> *msk_send_buffer;

    while (imaskcount < ind_set_tbl_size)
    {
        for(i = 0; i < (int)children_size; i++)
        {
            msk_send_index = get_msk_send_buffer_index (T);
            //DEBUG("send index: %d\n", msk_send_index);
            msk_send_buffer =  T->msk_send_pool[msk_send_index];
            //msk_send_request = (*T->msk_send_requests)[msk_send_index];
            for (j = 0; (j < request_size) && (imaskcount + j < ind_set_tbl_size); j++)
            {
                temp_set.mask->zeroize();
                for(int z = 0; z < w; z++)
                {
                    // This part could be sped up by considering the words and then doing some xoring, I think
                    if( cvec[i][k_bag_vec[z]] )
                    {
                        if(indset_masks[imaskcount + j].test_bit(z))
                        {
                            // Bit z is set, and the corresponding vertex is also in the child
                            temp_set.mask->set_bit(z);                        
                        }
                            
                    }
                }

                for (kk = 0; kk < T->num_mask_words; kk++)
                    (*msk_send_buffer)[j*T->num_mask_words + kk] = temp_set.mask->words[kk];
            }
            // Find the location of the child (k_adj_vec[i + 1]) using
         // T->tbl_loc vector.Then send the masks to get the
         // weights, we do this asynchonously therefore we can move
         // forward with other requests too.  
         // Note: k_adj_vec[0] is the parent.

         // Tags:
         // Since there will be bunch of requests, different
         // requests are identified using the tag. 
         // Tag is used in following way,
         // i*MPI_REQ_WGH_TAG + iteration
         // where i is the position of the child in the bag and
         // iteration is the request number. This is important when
         // indset_masks table is larger than request_size. Because
         // there has to be multiple requests. 
         //DEBUG("request_size: %d j : %d iteration : %d\n",
         // T->request_size, j, iteration);
            (*msk_send_buffer)[j*T->num_mask_words] = children[i];
            (*msk_send_buffer)[j*T->num_mask_words + 1] = iteration;
            destination = (*T->tbl_loc)[children[i]];
            //DEBUG("child: %d iteration : %d destination : %d\n",  children[i], iteration, destination);
            MPI_Isend ((void *)&(*msk_send_buffer).front(), j*T->num_mask_words + 2, MPI_UNSIGNED_LONG_LONG, 
                       destination, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_send_requests)[msk_send_index]);
        }
        iteration ++;
        imaskcount += request_size;
    }

    DEBUG("processor : %d iterations %d \n", T->rank, iteration);
    int local_children = 0;
    // This extra messages will help to make sure that any request for
    // weights get lost during the communication. Whoever has table
    // should receive this many requests.
    for(i = 0; i < (int)children_size; i++)
    {
        // DEBUG("send iteration: %d %d\n", iteration, (*T->tbl_loc)[children[i]]);
        MPI_Send ((void *)&iteration, 1, MPI_INT, 
                   (*T->tbl_loc)[children[i]], MPI_REQ_COUNT_TAG, MPI_COMM_WORLD);

        if ((*T->tbl_loc)[children[i]] == T->rank) local_children ++;
    }

    DEBUG("processor : %d has local children : %d\n", T->rank, local_children);
    // Current compute node ready to serve weight requests, 
    // Only if,
    // 1. If current tree node has children
    // 2. If current compute node has at least one child table 
    //MPI_Status status;
    request_flag = 0;
    vector<int> rcount(T->size, 0);
    int countflag = 0;
    int countindex = 0;
    int rcomplete_flag = 0;

    int msk_recv_count_elements = 0;
    int child_node = 0;
    int wgh_elem = ((msk_recv_count_elements - 1)/T->num_mask_words);
    MPI_Status recvstat;
    //MPI_Status wreq_status;
    MPI_Status msk_recv_status;
    int msk_recv_index = -1;
    int msk_recv_flag = 0;
    MPI_Status wgh_send_status;
    vector<unsigned long long> msk_recv_buffer;
    int wgh_send_index = -1;

    if (child_tbl_flag)
    {
        while (*(min_element(T->request_count->begin(), T->request_count->end())) < 0)
        {
            int tmp = -1;
            MPI_Recv((void *)&tmp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_REQ_COUNT_TAG, MPI_COMM_WORLD, &recvstat);
            (*T->request_count)[recvstat.MPI_SOURCE] = tmp;
        }
        
        // for (i = 0; i < T->size; i++)
        //     DEBUG("i : %d received : %d \n", i, (*T->request_count)[i]);

        while (1)
        {
            MPI_Waitany(T->size, &T->msk_recv_requests->front(), &msk_recv_index, &msk_recv_status);
            msk_recv_flag = 1;

            while (msk_recv_flag)
            {

                wgh_send_index = get_wgh_send_buffer_index(T);
                wgh_response = T->wgh_send_pool[wgh_send_index];

                wgh_response->clear();
                MPI_Get_count (&msk_recv_status, MPI_UNSIGNED_LONG_LONG, &msk_recv_count_elements);
                if (msk_recv_count_elements < 0)
                {
                    FERROR("Received elements can not be negative, Abort !\n");
                    MPI_Finalize();
                    exit (-1);
                }

                wgh_elem = 0;
                if (msk_recv_count_elements > 2)
                {
                    // Computing neighbor and iteration using last element
                    // of the request.
                    // Get masks and populate wgh_response vector. 
                    // 1st element is the iteration number.
                    
                    //(*wgh_response)[0] =
                    //(*T->chd_masks_recv[request_index])[msk_recv_count_elements
                    //- 1];
                    msk_recv_buffer = (*T->msk_recv_pool[msk_recv_index]);
                    wgh_response->push_back(msk_recv_buffer[msk_recv_count_elements - 1]);
                    child_node = msk_recv_buffer[msk_recv_count_elements - 2];
                    wgh_elem = ((msk_recv_count_elements - 2)/T->num_mask_words);

                    if (T->rank != (*T->tbl_loc)[child_node])
                    {
                        if ((*T->tbl_loc)[child_node] > -1)
                        {
                            FERROR("received wrong request for child : %d at node : %d, but should be at node : %d\n",
                                  child_node, T->rank, (*T->tbl_loc)[child_node]); 
                            FERROR("source is : %d iteration : %d\n", msk_recv_status.MPI_SOURCE, (*wgh_response)[0]);
                            MPI_Finalize();
                            exit (-1);
                        }
                    }
                }

                //DEBUG("wgh elem: %d\n", wgh_elem);
                for (i = 0; i < wgh_elem; i++)
                {
                    hash_ptr = NULL;
                    temp_set.mask->zeroize ();

                    for (kk = 0; kk < T->num_mask_words; kk++)
                        temp_set.mask->words[kk] = msk_recv_buffer[i*T->num_mask_words + kk];

                    HASH_FIND(hh, T->tree_nodes[child_node]->hash_table, temp_set.mask->words,
                              T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                    if (!hash_ptr)
                    {
                        temp_set.mask->print(__log_file__);
                        FERROR("i : %d can not find the mask in %d table, Abort !!", i, child_node);
                        MPI_Finalize();
                        exit (-1);
                    }

                    wgh_response->push_back(hash_ptr->value);
                }

                // Send weights back to the requester
                // if (!(T->wgh_send_count < pool_size))
                //     DEBUG("weights : %d:%d iteration: %d ready to send to : %d reuse: %d\n", wgh_response->size(), wgh_elem, 
                //           (*wgh_response)[0], msk_recv_status.MPI_SOURCE, wgh_send_index);
                
                if (wgh_elem > 0)
                {
                    MPI_Isend ((void *)&(wgh_response->front()), (wgh_elem + 1), MPI_INT,
                               msk_recv_status.MPI_SOURCE, MPI_WGH_TAG, MPI_COMM_WORLD, 
                               &(*T->wgh_send_requests)[wgh_send_index]);
                }
                else
                {
                    MPI_Isend ((void *)&(wgh_response->front()), 0, MPI_INT, 
                               msk_recv_status.MPI_SOURCE, MPI_WGH_TAG, MPI_COMM_WORLD, 
                               &(*T->wgh_send_requests)[wgh_send_index]);
                }

                T->msk_recv_pool[msk_recv_index]->assign(request_size*T->num_mask_words + 2, 0);
                rcount[msk_recv_index] ++;
                
                MPI_Irecv((void *)&T->msk_recv_pool[msk_recv_index]->front(), 
                          request_size*T->num_mask_words + 2, 
                          MPI_UNSIGNED_LONG_LONG, msk_recv_index, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, 
                          &(*T->msk_recv_requests)[msk_recv_index]);

                msk_recv_index = -1;
                msk_recv_flag = 0;
                MPI_Testany (T->size, &T->msk_recv_requests->front(), &msk_recv_index, &msk_recv_flag, &msk_recv_status);
            }
            
            for (i = 0; i < T->size; i++)
            {
                if ((rcount[i] < (*T->request_count)[i]*local_children))
                {
                    rcomplete_flag = 0;
                    break;
                }
                rcomplete_flag = 1;
            }

            if (rcomplete_flag) break;

        }
        T->request_count->assign(T->size, -1);
    }
    DEBUG("finished serving requests\n");


    // child_iterations vector make sure that, we have recieved
    // intended number of requests from compute nodes where children
    // tables are located.
    vector<int> child_iterations(T->size, 0);
    for (i = 0; i < children_size; i++)
    {
        child_iterations[(*T->tbl_loc)[children[i]]] += iteration;
    }

    // tbl_loc vector make sures that nodes recieves requests from
    // intended sources.
    vector<int> tbl_loc(T->size, 0);
    for (i = 0; i < children_size; i++)
    {
        //GEN("%d:%d  ", children[i], (*T->tbl_loc)[children[i]]);
        tbl_loc[(*T->tbl_loc)[children[i]]] = 1;
    }
    //GEN("\n");
    vector<int> wgh_recv_buffer(request_size + 1, 0);
    T->iter_count->assign(T->size, 0);
    
    // Receiving weights for the masks sent over
    if (children_size*iteration > 0)
    {
        int wgh_recv_index;
        int wgh_recv_flag = 0;
        int start_pos = 0;
        int recv_count = 0;
        int iter = 0;
        j = 0;

        DEBUG("waiting\n");
        while (1)
        {
            MPI_Status wgh_recv_status;
            wgh_recv_index = 0;
            MPI_Waitany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_status);
            //DEBUG("wgh_index: %d\n", wgh_recv_index);
            while(!tbl_loc[wgh_recv_index])
            {
              T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
              MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
                         wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);
              MPI_Waitany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_status);
            }
            
            recv_count = 0;
            MPI_Get_count (&wgh_recv_status, MPI_INT, &recv_count);
            //DEBUG("source : %d iteration : %d recv_count: %d\n",
            //status.MPI_SOURCE, (*T->wgh_recv[wgh_recv_index])[0],
            //recv_count);
            wgh_recv_buffer = (*T->wgh_recv_pool[wgh_recv_index]);
            T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
            MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
                       wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);


            iter = wgh_recv_buffer[0];
            start_pos = request_size*iter;
            (*T->iter_count)[wgh_recv_index]++;
            //DEBUG("recv count: %d\n", recv_count);
            for (i = 0; i < recv_count - 1; i++)
            {
                if (start_pos + i > ind_wgh.size())
                {
                    FERROR("location : %d is out of boundaries of the array: %d mask array size: %d \n", 
                          start_pos + i, ind_wgh.size(), indset_masks.size());
                    FERROR("source : %d iteration : %d recv_count: %d\n", wgh_recv_status.MPI_SOURCE, 
                          wgh_recv_buffer[0], recv_count);
                    MPI_Finalize();
                    exit (-1);
                }
                ind_wgh[start_pos + i] += wgh_recv_buffer[i + 1];
            }


            j = 1;
            for (i = 0; i < T->size; i++)
            {
                if ((*T->iter_count)[i] < child_iterations[i])
                {
                    //GEN("%d:%d:%d ", i, (*T->iter_count)[i], child_iterations[i]); 
                    j = 0;
                    break;
                }
                //GEN("\n");
            }
            if (j) break;
        }
        
        DEBUG("iteration counts \n");
         for (i = 0; i < T->size; i++)
             GEN("%d:%d ", i, (*T->iter_count)[i]); 
         GEN("\n");
    }

    DEBUG("finished recv weights\n");
    // Done with child DP

    // Free children table and  update tbl_entries values.
    T->free_children(k);
    for (i = 0; i < children_size; i++)
    {
        (*T->tbl_entries)[(*T->tbl_loc)[children[i]]] -= (*T->node_sizes)[children[i]];
        (*T->node_sizes)[children[i]] = 0;
        (*T->tbl_loc)[children[i]] = -1;
    }


    TDSolution *temp_sol;
    int parent_weight = 0;

    // Updating weight for intersection with the parent.
    int maxvalue = 0;
    vector<int> weight = T->G->get_weight();
    for (i = 0; i < ind_set_tbl_size; i++)
    {
        temp_set.mask->zeroize();
        parent_weight=0;
        current_mask = indset_masks[i];
    
        for(j = 0; j < w; j++)
        {                
            if(current_mask.test_bit(j))
            { 
                ind_wgh[i] += weight[k_bag_vec[j]];
                if(parent_location_vec[k_bag_vec[j]]!=-1)
                {
                    temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
                    parent_weight += weight[k_bag_vec[j]];
                }
            }
        }
        
        // temp_mask is a mask of this set \cap parent_bag
        hash_ptr=NULL;
        
        // k is not the root - so do the usual
        HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (k == T->root_node)
        {
            if (ind_wgh[i] > maxvalue)
            {
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
                *(added_set->mask)=*(temp_set.mask);
                *(added_set->orig_mask)=current_mask;
                
                added_set->value = ind_wgh[i];
                maxvalue = ind_wgh[i];
                
                HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
        else
        {
            
            if(hash_ptr)
            {
                if (ind_wgh[i] - parent_weight > hash_ptr->value)
                {
                    hash_ptr->value = ind_wgh[i] - parent_weight;
                }
            }
            else
            {
                // There is no entry for this mask - add it to the table
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
            *(added_set->mask)=*(temp_set.mask);
            *(added_set->orig_mask)=current_mask;
            
            added_set->value = ind_wgh[i] - parent_weight;
            
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
    }
    
    DEBUG("Hash tables ready to save\n");

    // for (i = 0; i < T->size; i++)
    //   GEN("%d:%d ", i, (*T->tbl_entries)[i]);
    //GEN("\n");
    
    // Find a location to save the table.
    int minloc = -1;
    int minent = INT_MAX;
    for (i = 0; i < T->size; i++)
    {
        if ((*T->tbl_entries)[i] < minent)
        {
            minent = (*T->tbl_entries)[i];
            minloc = i;
        }
    }
    (*T->tbl_loc)[k] = minloc;
    DEBUG("table for node : %d will be saved in : %d\n", k, minloc);

    // Before reaching this point we have calculate weights for the
    // masks in indset_masks vector.
	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
    DEBUG("size of the table : %d\n", table_size);
    DEBUG("Time for DP loop : %lf\n", MPI_Wtime() - stime);

    hash_ptr = NULL;
    i = 0;
    int masksize = T->num_mask_words + 1;
    vector<unsigned long long> maskweight_array(table_size*(masksize));

    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        for (j = 0; j < T->num_mask_words; j++)
        {
            maskweight_array[i*(masksize) + j] = hash_ptr->mask->words[j];
        }
        maskweight_array[i*(masksize) + j] = hash_ptr->value;
        i++;
    }

    vector<int> n_elements_vec(T->size, 0);
    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    DEBUG("In gather phase\n");
    
    MPI_Allgather(&table_size, 1, MPI_INT, &n_elements_vec.front(), 1, MPI_INT, MPI_COMM_WORLD);
    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elements_vec[i];
        w_recv_counts[i] = n_elements_vec[i];
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

    (*T->tbl_entries)[minloc] += total_entries;
    (*T->node_sizes)[k] = total_entries;

    // Deriving datatype
    MPI_Datatype mask_type;
    MPI_Type_contiguous (T->num_mask_words + 1, MPI_UNSIGNED_LONG_LONG, &mask_type);
    MPI_Type_commit(&mask_type);
    vector<unsigned long long> maskweights_vec(total_entries*masksize, 0);    

    DEBUG("starting Gather v\n");
    MPI_Gatherv((void *)&maskweight_array.front(), w_recv_counts[T->rank], mask_type,
                (void *)&maskweights_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), mask_type, minloc, MPI_COMM_WORLD);
    DEBUG("complete allgather v\n");

    unsigned long long max = 0;
    DEBUG("total entries : %d\n", total_entries);
    DEBUG("entering loop to update hash table \n");

    stime = MPI_Wtime ();
    // Loop to  update hash table after gathering values at other
    // processors 
    DEBUG("gather all results to %d\n", minloc);
    for (i = 0; i < total_entries && (T->rank == minloc); i++)
    {

       if (T->root_node == k)
        {
            unsigned long long cval = maskweights_vec[i*masksize + T->num_mask_words];
            if ( cval > max)
            {
                max = cval;
            }
            continue;
        }

        TDSolution *added_set;
        added_set=new TDSolution(T->num_mask_words);
        if (!added_set)
        {
            FERROR("no enough memory \n");
            MPI_Finalize();
            exit(-1);
        }

        for (j = 0; j < T->num_mask_words; j++)
        {
            added_set->mask->words[j] = maskweights_vec[i*masksize + j];
        }
        added_set->value= maskweights_vec[i*masksize + j];

        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
        }
        else
        {
            if (added_set->value > hash_ptr->value)
            {
                hash_ptr->value = added_set->value;
            }
        }
    }

    DEBUG("completed table size : %d\n", (int)HASH_COUNT(T->tree_nodes[k]->hash_table));
    if (T->root_node == k)
    {
        DEBUG("maximum value %d\n", max);
		// CSG adding just to check answer
        if (T->rank == minloc)
            fprintf(stderr,"Max value is %d\n",max);
    }
    DEBUG("Finish hash table update\n");
    DEBUG("Time for hash table update : %lf\n", MPI_Wtime() - stime);

    // clear some vectors to gain more memory.
     maskweights_vec.clear();
     maskweight_array.clear();
     w_recv_disp.clear();
     w_recv_counts.clear();
     indset_masks.clear();
     //  wgh_request_pool.clear();

    
    // for (i = 0; i < pool_size; i++)
    //     if (children_size > 0)
    //         delete wgh_vec_pool[i];

    for (i = 0; i < entries; i++)
        delete start_masks[i];

	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	
    // Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}

void recv_msk_send_wgh(TDTree *T)
{
    int i;
    int kk;

	TDSolution *hash_ptr=NULL;
	TDSolution temp_set(T->num_mask_words);

    int msk_recv_count_elements = 0;
    int child_node = 0;
    int wgh_elem = 0;
    MPI_Status recvstat;

    MPI_Status msk_recv_status;
    int msk_recv_index = -1;
    int msk_recv_flag = 0;

    MPI_Status wgh_send_status;
    vector<unsigned long long> msk_recv_buffer;
    int wgh_send_index = -1;
    vector<int> *wgh_response;
    int static threadcount;

    MPI_Testany(T->size, &T->msk_recv_requests->front(), &msk_recv_index, &msk_recv_flag, &msk_recv_status);

    while (msk_recv_flag)
    {
        wgh_send_index = get_wgh_send_buffer_index(T);
        if (wgh_send_index < 0)
            break;

        wgh_response = T->wgh_send_pool[wgh_send_index];

        wgh_response->clear();
        MPI_Get_count (&msk_recv_status, MPI_UNSIGNED_LONG_LONG, &msk_recv_count_elements);
        if (msk_recv_count_elements < 0)
        {

            FERROR("Received elements can not be negative, Abort !\n");
            MPI_Finalize();
            exit (-1);
        }

        wgh_elem = 0;
        if (msk_recv_count_elements > 2)
        {
            // Computing neighbor and iteration using last element
            // of the request.
            // Get masks and populate wgh_response vector. 
            // 1st element is the iteration number.

            msk_recv_buffer = (*T->msk_recv_pool[msk_recv_index]);
            wgh_response->push_back(msk_recv_buffer[msk_recv_count_elements - 1]);
            child_node = msk_recv_buffer[msk_recv_count_elements - 2];
            wgh_elem = ((msk_recv_count_elements - 2)/T->num_mask_words);

            if (T->rank != (*T->tbl_loc)[child_node])
            {
                if ((*T->tbl_loc)[child_node] > -1)
                {
                    FERROR("received wrong request for child : %d at node : %d, but should be at node : %d\n",
                          child_node, T->rank, (*T->tbl_loc)[child_node]); 
                    FERROR("source is : %d iteration : %d\n", msk_recv_status.MPI_SOURCE, (*wgh_response)[0]);
                    MPI_Finalize();
                    exit (-1);
                }
            }
        }

        for (i = 0; i < wgh_elem; i++)
        {
            hash_ptr = NULL;
            temp_set.mask->zeroize ();

            for (kk = 0; kk < T->num_mask_words; kk++)
                temp_set.mask->words[kk] = msk_recv_buffer[i*T->num_mask_words + kk];

            HASH_FIND(hh, T->tree_nodes[child_node]->hash_table, temp_set.mask->words,
                      T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

            if (!hash_ptr)
            {
                temp_set.mask->print(__log_file__);
                FERROR("i : %d can not find the mask in %d table, Abort !!", i, child_node);
                MPI_Finalize();
                exit (-1);
            }

            wgh_response->push_back(hash_ptr->value);
        }


        // counting received and responded requests.
        (*T->rcount)[msk_recv_index]++;
        if (wgh_elem > 0)
        {
            MPI_Isend ((void *)&(wgh_response->front()), (wgh_elem + 1), MPI_INT,
                       msk_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, 
                       &(*T->wgh_send_requests)[wgh_send_index]);
        }
        else
        {
            MPI_Isend ((void *)&(wgh_response->front()), 0, MPI_INT, 
                       msk_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, 
                       &(*T->wgh_send_requests)[wgh_send_index]);
        }
            
        // we can post/repost this recv as soon as we get the content
        // into a local buffer.
        T->msk_recv_pool[msk_recv_index]->assign(T->request_size*T->num_mask_words + 2, 0);
        MPI_Irecv((void *)&T->msk_recv_pool[msk_recv_index]->front(), 
                  T->request_size*T->num_mask_words + 2, 
                  MPI_UNSIGNED_LONG_LONG, msk_recv_index, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, 
                  &(*T->msk_recv_requests)[msk_recv_index]);

        msk_recv_index = -1;
        msk_recv_flag = 0;

        MPI_Testany (T->size, &T->msk_recv_requests->front(), &msk_recv_index, &msk_recv_flag, &msk_recv_status);
    }
}

bool recv_msk_send_wgh_9 (TDTree *T)
{
    int i;
    int kk;

	TDSolution *hash_ptr=NULL;
	TDSolution temp_set(T->num_mask_words);

    int msk_recv_count_elements = 0;
    int child_node = 0;
    int wgh_elem = 0;
    MPI_Status recvstat;

    MPI_Status msk_recv_status;
    int msk_recv_index = -1;
    int msk_recv_flag = 0;

    MPI_Status wgh_send_status;
    vector<unsigned long long> msk_recv_buffer;
    //int wgh_send_index = -1;
    vector<int> *wgh_response;
    int static threadcount;

    if (MPI_SUCCESS != 
        MPI_Testany(T->size, &T->msk_recv_requests->front(), &msk_recv_index, &msk_recv_flag, &msk_recv_status))
    {
        FERROR ("MPI Testany failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
        return false;
    }

    if (!msk_recv_flag)
        return false;

    wgh_response = T->wgh_send_pool[msk_recv_index];
    wgh_response->clear();

    if (MPI_SUCCESS != 
        MPI_Get_count (&msk_recv_status, MPI_UNSIGNED_LONG_LONG, &msk_recv_count_elements))
    {
        FERROR ("MPI Get count failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
        return false;
    }

    if (msk_recv_count_elements < 0)
    {

        FERROR("Received elements can not be negative, Abort !\n");
        MPI_Finalize();
        exit (0);
    }

    wgh_elem = 0;
    if (msk_recv_count_elements >= 2)
    {
        // Computing neighbor and iteration using last element
        // of the request.
        // Get masks and populate wgh_response vector. 
        // 1st element is the iteration number.
        //DEBUG("msk recv count : %d index : %d\n", msk_recv_count_elements, msk_recv_index);
        msk_recv_buffer = (*T->msk_recv_pool[msk_recv_index]);
        wgh_response->push_back(msk_recv_buffer[msk_recv_count_elements - 1]);
        child_node = msk_recv_buffer[msk_recv_count_elements - 2];

        wgh_response->push_back(child_node);
        wgh_elem = ((msk_recv_count_elements - 2)/T->num_mask_words);
    }
    else
    {
        FERROR("Wrong element count in the request\n");
        MPI_Finalize ();
        return false;
    }

    for (i = 0; i < wgh_elem; i++)
    {
        hash_ptr = NULL;
        temp_set.mask->zeroize ();

        for (kk = 0; kk < T->num_mask_words; kk++)
            temp_set.mask->words[kk] = msk_recv_buffer[i*T->num_mask_words + kk];

        HASH_FIND(hh, T->tree_nodes[child_node]->hash_table, temp_set.mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            temp_set.mask->print(__log_file__);
            FERROR("i : %d can not find the mask in %d table at processor : %d, Abort !!\n", i, child_node, T->rank);
            MPI_Finalize();
            return false;
        }

        wgh_response->push_back(hash_ptr->value);
    }

    // counting received and responded requests.
    //DEBUG("Send %d elements of child : %d to : %d\n", wgh_elem + 2, child_node, msk_recv_index);
    
    if (MPI_SUCCESS != 
        MPI_Send ((void *)&(wgh_response->front()), (wgh_elem + 2), MPI_INT,
                   msk_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD))
    {
        FERROR ("MPI Isend failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
        return false;
    }

    // we can post/repost this recv as soon as we get the content
    // into a local buffer.
    T->msk_recv_pool[msk_recv_index]->assign(T->request_size*T->num_mask_words + 2, 0);
    if (MPI_SUCCESS != MPI_Irecv((void *)&T->msk_recv_pool[msk_recv_index]->front(), 
                                 T->request_size*T->num_mask_words + 2, 
                                 MPI_UNSIGNED_LONG_LONG, msk_recv_index, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, 
                                 &(*T->msk_recv_requests)[msk_recv_index]))
    {
        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
        return false;
    }

    msk_recv_index = -1;
    msk_recv_flag = 0;
    return true;
}

void *thread_recv_msk_send_wgh(void *v)
{
    int i;
    int kk;

    TDTree *T = (TDTree *)v;

	TDSolution *hash_ptr=NULL;
	TDSolution temp_set(T->num_mask_words);

    int msk_recv_count_elements = 0;
    int child_node = 0;
    int wgh_elem = 0;
    MPI_Status recvstat;

    MPI_Status msk_recv_status;
    int msk_recv_index = -1;
    int msk_recv_flag = 0;

    MPI_Status wgh_send_status;
    vector<unsigned long long> msk_recv_buffer;
    int wgh_send_index = -1;
    vector<int> *wgh_response;
    int static threadcount;
    while (1)
    {

        while (!msk_recv_flag)
        {
            if (T->flag_recv_done)
                break;
            MPI_Testany(T->size, &T->msk_recv_requests->front(), &msk_recv_index, &msk_recv_flag, &msk_recv_status);
        }
        
        // When flag_recv_done is set, thread is preaparing for join.
        // T->flag_recv_done set in the main thread after all nodes finished.
        if (T->flag_recv_done)
            break;

        while (msk_recv_flag)
        {

            wgh_send_index = get_wgh_send_buffer_index(T);
            wgh_response = T->wgh_send_pool[wgh_send_index];

            wgh_response->clear();
            MPI_Get_count (&msk_recv_status, MPI_UNSIGNED_LONG_LONG, &msk_recv_count_elements);
            if (msk_recv_count_elements < 0)
            {

                FERROR("Received elements can not be negative, Abort !\n");
                MPI_Finalize();
                exit (-1);
            }

            wgh_elem = 0;
            if (msk_recv_count_elements > 2)
            {
                // Computing neighbor and iteration using last element
                // of the request.
                // Get masks and populate wgh_response vector. 
                // 1st element is the iteration number.

                msk_recv_buffer = (*T->msk_recv_pool[msk_recv_index]);
                wgh_response->push_back(msk_recv_buffer[msk_recv_count_elements - 1]);
                child_node = msk_recv_buffer[msk_recv_count_elements - 2];
                wgh_elem = ((msk_recv_count_elements - 2)/T->num_mask_words);

                if (T->rank != (*T->tbl_loc)[child_node])
                {
                    if ((*T->tbl_loc)[child_node] > -1)
                    {
                        FERROR("received wrong request for child : %d at node : %d, but should be at node : %d\n",
                              child_node, T->rank, (*T->tbl_loc)[child_node]); 
                        FERROR("source is : %d iteration : %d\n", msk_recv_status.MPI_SOURCE, (*wgh_response)[0]);
                        MPI_Finalize();
                        exit (-1);
                    }
                }
            }

            for (i = 0; i < wgh_elem; i++)
            {
                hash_ptr = NULL;
                temp_set.mask->zeroize ();

                for (kk = 0; kk < T->num_mask_words; kk++)
                    temp_set.mask->words[kk] = msk_recv_buffer[i*T->num_mask_words + kk];

                HASH_FIND(hh, T->tree_nodes[child_node]->hash_table, temp_set.mask->words,
                          T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                if (!hash_ptr)
                {
                    temp_set.mask->print(__log_file__);
                    FERROR("i : %d can not find the mask in %d table, Abort !!", i, child_node);
                    MPI_Finalize();
                    exit (-1);
                }

                wgh_response->push_back(hash_ptr->value);
            }

            if (wgh_elem > 0)
            {
                MPI_Isend ((void *)&(wgh_response->front()), (wgh_elem + 1), MPI_INT,
                           msk_recv_status.MPI_SOURCE, MPI_WGH_TAG, MPI_COMM_WORLD, 
                           &(*T->wgh_send_requests)[wgh_send_index]);
            }
            else
            {
                MPI_Isend ((void *)&(wgh_response->front()), 0, MPI_INT, 
                           msk_recv_status.MPI_SOURCE, MPI_WGH_TAG, MPI_COMM_WORLD, 
                           &(*T->wgh_send_requests)[wgh_send_index]);
            }
            
            // we can post/repost this recv as soon as we get the content
            // into a local buffer.
            T->msk_recv_pool[msk_recv_index]->assign(T->request_size*T->num_mask_words + 2, 0);
            MPI_Irecv((void *)&T->msk_recv_pool[msk_recv_index]->front(), 
                      T->request_size*T->num_mask_words + 2, 
                      MPI_UNSIGNED_LONG_LONG, msk_recv_index, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, 
                      &(*T->msk_recv_requests)[msk_recv_index]);

            msk_recv_index = -1;
            msk_recv_flag = 0;

            MPI_Testany (T->size, &T->msk_recv_requests->front(), &msk_recv_index, &msk_recv_flag, &msk_recv_status);
        }
    }
}



static int compute_nonnice_table_parallel_7 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
	vector<int> I(T->tree_nodes[k]->bag.size());
    list<int>::iterator it_li;
    int kk;
    int entries = 32768;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());
    int children_size = (int) children.size();

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    vector<bigint_t *> start_masks (entries);
    vector<bigint_t *> end_masks (entries);

    DEBUG("mask generation\n");
    parallel_mask_gen (T, k, w, &start_masks, &end_masks, entries);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;

    vector<bigint_t> indset_masks;
    vector<int> ind_wgh;
    unsigned long long ind_counter = 0;
    
 
	bigint_t current_mask(T->num_mask_words);
	bigint_t end_mask(T->num_mask_words);

	TDSolution temp_set(T->num_mask_words);

    double stime = MPI_Wtime();
    int request_flag = 0;
    int request_index = -1;

    unsigned long long start_pos, end_pos, wchunk, counter, stwork;
    wchunk = entries/T->size;
    start_pos = (wchunk * T->rank);
    end_pos = (wchunk * (T->rank + 1));

    vector<int> requests(T->size, 0);
    vector<int> terminate(T->size, -1);

    DEBUG("DP loop\n");
    counter = start_pos;
    while (1)
    {
        
        // DEBUG("counter %d end_pos %d\n", counter, end_pos);
        if (counter < end_pos)
        {
            current_mask = *(start_masks[counter]);
            end_mask = *(end_masks[counter]);
            counter ++;
        }

        while(current_mask < end_mask)
        {

            if ((T->rank == 0) && (counter < end_pos))
            {

                MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                while (request_flag)
                {
                    MPI_Irecv (&(T->requester_rank), 1, MPI_INT, request_index, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, 
                               &((*T->requests)[request_index]));
                    requests[request_index] = 1;
                    
                    request_flag = 0;
                    request_index = -1;
                    T->requester_rank = -1;
                    MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                }
                
                // Send work to the requester.
                for (i = 1; i < T->size && (end_pos > counter); i++)
                {
                    if (requests[i])
                    {
                        end_pos --;
                        MPI_Send (&end_pos, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                        requests[i] = 0;
                    }
                }
            }
            
            is_independent=true;

            for(i = 0; i < w ; i++) 
            {            
                if(current_mask.test_bit(nbr_mask_vec[i]->k))
                {
                    // current_mask has the bit corresponding to
                    // nbr_mask_vec[i] set - so and with the appropriate
                    // vector for a hit 
                    has_bit=false;
                    // Do the loop below manually rather than calling the overloaded AND
                    for(j=0; j<current_mask.S; j++)
                    {
                        if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                        {
                            // Some nbr is found in word j
                            has_bit=true;
                            break;
                        }
                    }
                    // breaks to here
                    if( has_bit )
                    {                
                        // Set represented by current_mask is not independent
                        is_independent=false;
                        int rpos=0;
                        // Find the rightmost bit in current_mask - this
                        // could be sped up by looking at the words 
                        // SPEEDUP
                        while(!(current_mask.test_bit(rpos)))
                            rpos++;

                        // Now advance current_mask by a potentially massive amount
                        // Would be faster here to figure out which word we are in and 
                        // do this w/o test_bit
                        int b=0;
                        while(current_mask.test_bit(rpos+b))
                            b++;
                        for(int a=0;a<b;a++)
                            current_mask.xor_bit(rpos+a);
                        current_mask.or_bit(rpos+b);
                        break;                   
                    }
                }
            }


            if(is_independent)
            {   
                // ind_counter keep track of how many independent sets
                // are in the indset_mask vector
                ind_counter ++;
                indset_masks.push_back(current_mask);
                ind_wgh.push_back(0);
                ++current_mask;
            }
        }


        // Once assigned work is finished, ask for more work
        if (!(counter < end_pos))
        {
            if (T->rank > 0)
            {
                MPI_Request maskrequest;
                MPI_Status maskstatus;
                int wflag = 0;
                int count = 0;
            
                //DEBUG("process : %d ask for work from : %d\n", T->rank, 0);
                MPI_Send(&(T->rank), 1, MPI_INT, 0, MPI_REQ_WRK_TAG, MPI_COMM_WORLD);
                MPI_Recv (&stwork, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_WRK_TAG, MPI_COMM_WORLD, &maskstatus);
                count = 0;
                MPI_Get_count (&maskstatus, MPI_UNSIGNED_LONG_LONG, &count);

                if (!(count > 0))
                    break;

                current_mask.zeroize ();
                end_mask.zeroize ();
                current_mask = *(start_masks[stwork]);
                end_mask = *(end_masks[stwork]);
            }
            else
            {
                //Sending terminate signal
                for (i = 1; i < T->size; i++)
                    MPI_Send (&end_pos, 0, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                break;
            }
        }
    }

    // Now do the DP over the children for this independent set - if this is a leaf, then
    // the body of the loop is never executed
    int ind_set_tbl_size = (int)indset_masks.size();
    MPI_Request wgh_s_request;
    MPI_Request wgh_r_request;
    int imaskcount = 0;
    int request_size = T->request_size;
    MPI_Status wstatus[T->size];
    request_index = 0;
    request_flag = 0;
    int child_tbl_flag = 0;
    vector<int> *wgh_response;

    // If a child table is stored in current node, child_tbl_flag will
    // be set to 1.  Which means current node should be ready to serve
    // weight requests. 
    for (i = 0; i < children_size; i++)
    {
        if (T->rank == (*T->tbl_loc)[children[i]])
        {
            child_tbl_flag = 1;
            break;
        }
    }
    

    DEBUG("computing mask weights\n");
    int iteration = 0;
    DEBUG("mask count %d tble size %d\n", imaskcount, ind_set_tbl_size);
    if (!(ind_set_tbl_size > 0))
    {
        for(i = 0; i < (int)children_size; i++)
            MPI_Isend ((void *)NULL, 0, MPI_UNSIGNED_LONG_LONG, 
                       (*T->tbl_loc)[children[i]], MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &wgh_s_request);
    }

    //DEBUG("sending out requests\n");
    DEBUG("# of nbrs %d\n", children_size);
    DEBUG("child:location: ");
    for ( i = 0; i < children_size; i++)
    {
        GEN("%d:%d ", children[i], (*T->tbl_loc)[children[i]]);
    }
    GEN("\n");

    // Sending request of size request_size.
    int count_send = 0;
    int destination = 0;
    int msk_send_index = 0;
    vector<unsigned long long> *msk_send_buffer;

    while (imaskcount < ind_set_tbl_size)
    {
        for(i = 0; i < (int)children_size; i++)
        {
            msk_send_index = get_msk_send_buffer_index (T);
            msk_send_buffer =  T->msk_send_pool[msk_send_index];

            for (j = 0; (j < request_size) && (imaskcount + j < ind_set_tbl_size); j++)
            {
                temp_set.mask->zeroize();
                for(int z = 0; z < w; z++)
                {
                    // This part could be sped up by considering the words and then doing some xoring, I think
                    if( cvec[i][k_bag_vec[z]] )
                    {
                        if(indset_masks[imaskcount + j].test_bit(z))
                        {
                            // Bit z is set, and the corresponding vertex is also in the child
                            temp_set.mask->set_bit(z);                        
                        }
                            
                    }
                }

                for (kk = 0; kk < T->num_mask_words; kk++)
                    (*msk_send_buffer)[j*T->num_mask_words + kk] = temp_set.mask->words[kk];
            }
            // Find the location of the child (k_adj_vec[i + 1]) using
         // T->tbl_loc vector.Then send the masks to get the
         // weights, we do this asynchonously therefore we can move
         // forward with other requests too.  
         // Note: k_adj_vec[0] is the parent.

         // Tags:
         // Since there will be bunch of requests, different
         // requests are identified using the tag. 
         // Tag is used in following way,
         // i*MPI_REQ_WGH_TAG + iteration
         // where i is the position of the child in the bag and
         // iteration is the request number. This is important when
         // indset_masks table is larger than request_size. Because
         // there has to be multiple requests. 
         //DEBUG("request_size: %d j : %d iteration : %d\n",
         // T->request_size, j, iteration);
            (*msk_send_buffer)[j*T->num_mask_words] = children[i];
            (*msk_send_buffer)[j*T->num_mask_words + 1] = iteration;
            destination = (*T->tbl_loc)[children[i]];
            MPI_Isend ((void *)&(*msk_send_buffer).front(), j*T->num_mask_words + 2, MPI_UNSIGNED_LONG_LONG, 
                       destination, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_send_requests)[msk_send_index]);
        }
        iteration ++;
        imaskcount += request_size;
    }

    DEBUG("processor : %d iterations %d \n", T->rank, iteration);
    int local_children = 0;
    // This extra messages will help to make sure that any request for
    // weights get lost during the communication. Whoever has table
    // should receive this many requests.
    for(i = 0; i < (int)children_size; i++)
    {
        MPI_Send ((void *)&iteration, 1, MPI_INT, 
                   (*T->tbl_loc)[children[i]], MPI_REQ_COUNT_TAG, MPI_COMM_WORLD);

        if ((*T->tbl_loc)[children[i]] == T->rank) local_children ++;
    }

    DEBUG("processor : %d has local children : %d\n", T->rank, local_children);
    // Current compute node ready to serve weight requests, 
    // Only if,
    // 1. If current tree node has children
    // 2. If current compute node has at least one child table 
    //MPI_Status status;
    request_flag = 0;

    // child_iterations vector make sure that, we have recieved
    // intended number of requests from compute nodes where children
    // tables are located.
    vector<int> child_iterations(T->size, 0);
    for (i = 0; i < children_size; i++)
    {
        child_iterations[(*T->tbl_loc)[children[i]]] += iteration;
    }

    // tbl_loc vector make sures that nodes recieves requests from
    // intended sources.
    vector<int> tbl_loc(T->size, 0);
    for (i = 0; i < children_size; i++)
        tbl_loc[(*T->tbl_loc)[children[i]]] = 1;

    vector<int> wgh_recv_buffer(request_size + 1, 0);
    
    // Receiving weights for the masks sent over
    if (children_size*iteration > 0)
    {
        int wgh_recv_index;
        int wgh_recv_flag = 0;
        int start_pos = 0;
        int recv_count = 0;
        int iter = 0;
        j = 0;

        DEBUG("waiting\n");
        while (1)
        {
            MPI_Status wgh_recv_status;
            wgh_recv_index = 0;
            MPI_Waitany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_status);

            while(!tbl_loc[wgh_recv_index])
            {
              T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
              MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
                         wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);
              MPI_Waitany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_status);
            }
            
            recv_count = 0;
            MPI_Get_count (&wgh_recv_status, MPI_INT, &recv_count);

            wgh_recv_buffer = (*T->wgh_recv_pool[wgh_recv_index]);
            T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
            MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
                       wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);


            iter = wgh_recv_buffer[0];
            start_pos = request_size*iter;
            (*T->iter_count)[wgh_recv_index]++;

            for (i = 0; i < recv_count - 1; i++)
            {
                if (start_pos + i > ind_wgh.size())
                {

                    FERROR("location : %d is out of boundaries of the array: %d mask array size: %d \n", 
                          start_pos + i, ind_wgh.size(), indset_masks.size());
                    FERROR("source : %d iteration : %d recv_count: %d\n", wgh_recv_status.MPI_SOURCE, 
                          wgh_recv_buffer[0], recv_count);
                    MPI_Finalize();
                    exit (-1);
                }
                ind_wgh[start_pos + i] += wgh_recv_buffer[i + 1];
            }


            j = 1;
            for (i = 0; i < T->size; i++)
            {
                if ((*T->iter_count)[i] < child_iterations[i])
                {
                    j = 0;
                    break;
                }
            }
            if (j) break;
        }
        
        DEBUG("iteration counts \n");
         for (i = 0; i < T->size; i++)
             GEN("%d:%d ", i, (*T->iter_count)[i]); 
         GEN("\n");
    }

    DEBUG("finished recv weights\n");
    // Done with child DP


    TDSolution *temp_sol;
    int parent_weight = 0;

    // Updating weight for intersection with the parent.
    int maxvalue = 0;
    vector<int> weight = T->G->get_weight();
    for (i = 0; i < ind_set_tbl_size; i++)
    {
        temp_set.mask->zeroize();
        parent_weight=0;
        current_mask = indset_masks[i];
    
        for(j = 0; j < w; j++)
        {                
            if(current_mask.test_bit(j))
            { 
                ind_wgh[i] += weight[k_bag_vec[j]];
                if(parent_location_vec[k_bag_vec[j]]!=-1)
                {
                    temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
                    parent_weight += weight[k_bag_vec[j]];
                }
            }
        }
        
        // temp_mask is a mask of this set \cap parent_bag
        hash_ptr=NULL;
        
        // k is not the root - so do the usual
        HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (k == T->root_node)
        {
            if (ind_wgh[i] > maxvalue)
            {
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
                *(added_set->mask)=*(temp_set.mask);
                *(added_set->orig_mask)=current_mask;
                
                added_set->value = ind_wgh[i];
                maxvalue = ind_wgh[i];
                
                HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
        else
        {
            
            if(hash_ptr)
            {
                if (ind_wgh[i] - parent_weight > hash_ptr->value)
                {
                    hash_ptr->value = ind_wgh[i] - parent_weight;
                }
            }
            else
            {
                // There is no entry for this mask - add it to the table
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
            *(added_set->mask)=*(temp_set.mask);
            *(added_set->orig_mask)=current_mask;
            
            added_set->value = ind_wgh[i] - parent_weight;
            
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
    }
    
    DEBUG("Hash tables ready to save\n");

    
    // Find a location to save the table.
    int minloc = -1;
    int minent = INT_MAX;
    for (i = 0; i < T->size; i++)
    {
        if ((*T->tbl_entries)[i] < minent)
        {
            minent = (*T->tbl_entries)[i];
            minloc = i;
        }
    }
    (*T->tbl_loc)[k] = minloc;
    DEBUG("table for node : %d will be saved in : %d\n", k, minloc);

    // Before reaching this point we have calculate weights for the
    // masks in indset_masks vector.
	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
    DEBUG("size of the table : %d\n", table_size);
    DEBUG("Time for DP loop : %lf\n", MPI_Wtime() - stime);

    hash_ptr = NULL;
    i = 0;
    int masksize = T->num_mask_words + 1;
    vector<unsigned long long> maskweight_array(table_size*(masksize));

    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        for (j = 0; j < T->num_mask_words; j++)
        {
            maskweight_array[i*(masksize) + j] = hash_ptr->mask->words[j];
        }
        maskweight_array[i*(masksize) + j] = hash_ptr->value;
        i++;
    }

    vector<int> n_elements_vec(T->size, 0);
    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    DEBUG("In gather phase\n");
    
    MPI_Allgather(&table_size, 1, MPI_INT, &n_elements_vec.front(), 1, MPI_INT, MPI_COMM_WORLD);
    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elements_vec[i];
        w_recv_counts[i] = n_elements_vec[i];
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

    (*T->tbl_entries)[minloc] += total_entries;
    (*T->node_sizes)[k] = total_entries;

    // Deriving datatype
    MPI_Datatype mask_type;
    MPI_Type_contiguous (T->num_mask_words + 1, MPI_UNSIGNED_LONG_LONG, &mask_type);
    MPI_Type_commit(&mask_type);
    vector<unsigned long long> maskweights_vec(total_entries*masksize, 0);    

    DEBUG("starting Gather v\n");
    MPI_Gatherv((void *)&maskweight_array.front(), w_recv_counts[T->rank], mask_type,
                (void *)&maskweights_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), mask_type, minloc, MPI_COMM_WORLD);

    unsigned long long max = 0;
    DEBUG("total entries : %d\n", total_entries);
    DEBUG("entering loop to update hash table \n");

    stime = MPI_Wtime ();
    // Loop to  update hash table after gathering values at other
    // processors 
    DEBUG("gather all results to %d\n", minloc);
    for (i = 0; i < total_entries && (T->rank == minloc); i++)
    {

       if (T->root_node == k)
        {
            unsigned long long cval = maskweights_vec[i*masksize + T->num_mask_words];
            if ( cval > max)
            {
                max = cval;
            }
            continue;
        }

        TDSolution *added_set;
        added_set=new TDSolution(T->num_mask_words);
        if (!added_set)
        {
            FERROR("no enough memory \n");
            MPI_Finalize();
            exit(-1);
        }

        for (j = 0; j < T->num_mask_words; j++)
        {
            added_set->mask->words[j] = maskweights_vec[i*masksize + j];
        }
        added_set->value= maskweights_vec[i*masksize + j];

        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
        }
        else
        {
            if (added_set->value > hash_ptr->value)
            {
                hash_ptr->value = added_set->value;
            }
        }
    }

    DEBUG("completed table size : %d\n", (int)HASH_COUNT(T->tree_nodes[k]->hash_table));
    if (T->root_node == k)
    {
        DEBUG("maximum value %d\n", max);
		// CSG adding just to check answer
        if (T->rank == minloc)
            fprintf(stderr,"Max value is %d\n",max);
    }
    DEBUG("Finish hash table update\n");
    DEBUG("Time for hash table update : %lf\n", MPI_Wtime() - stime);

    // clear some vectors to gain more memory.
     maskweights_vec.clear();
     maskweight_array.clear();
     w_recv_disp.clear();
     w_recv_counts.clear();
     indset_masks.clear();

    for (i = 0; i < children_size; i++)
    {
        (*T->tbl_entries)[(*T->tbl_loc)[children[i]]] -= (*T->node_sizes)[children[i]];
        (*T->node_sizes)[children[i]] = 0;
        (*T->tbl_loc)[children[i]] = -1;
    }

    
    for (i = 0; i < entries; i++)
        delete start_masks[i];

	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	
    // Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}

static int compute_nonnice_table_parallel_8 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
	vector<int> I(T->tree_nodes[k]->bag.size());
    list<int>::iterator it_li;
    int kk;
    int entries = 32768;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());
    int children_size = (int) children.size();

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    vector<bigint_t *> start_masks (entries);
    vector<bigint_t *> end_masks (entries);

    DEBUG("mask generation\n");
    parallel_mask_gen (T, k, w, &start_masks, &end_masks, entries);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;

    vector<bigint_t> indset_masks;
    vector<int> ind_wgh;
    unsigned long long ind_counter = 0;
    
 
	bigint_t current_mask(T->num_mask_words);
	bigint_t end_mask(T->num_mask_words);

	TDSolution temp_set(T->num_mask_words);

    double stime = MPI_Wtime();
    int request_flag = 0;
    int request_index = -1;

    unsigned long long start_pos, end_pos, wchunk, counter, stwork;
    wchunk = entries/T->size;
    start_pos = (wchunk * T->rank);
    end_pos = (wchunk * (T->rank + 1));

    vector<int> requests(T->size, 0);
    vector<int> terminate(T->size, -1);

    DEBUG("DP loop\n");
    counter = start_pos;
    while (1)
    {
        
        // DEBUG("counter %d end_pos %d\n", counter, end_pos);
        if (counter < end_pos)
        {
            current_mask = *(start_masks[counter]);
            end_mask = *(end_masks[counter]);
            counter ++;
        }

        while(current_mask < end_mask)
        {

            if ((T->rank == 0) && (counter < end_pos))
            {

                MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                while (request_flag)
                {
                    MPI_Irecv (&(T->requester_rank), 1, MPI_INT, request_index, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, 
                               &((*T->requests)[request_index]));
                    requests[request_index] = 1;
                    
                    request_flag = 0;
                    request_index = -1;
                    T->requester_rank = -1;
                    MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE);
                }
                
                // Send work to the requester.
                for (i = 1; i < T->size && (end_pos > counter); i++)
                {
                    if (requests[i])
                    {
                        end_pos --;
                        MPI_Send (&end_pos, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                        requests[i] = 0;
                    }
                }
            }
            
            is_independent=true;

            for(i = 0; i < w ; i++) 
            {            
                if(current_mask.test_bit(nbr_mask_vec[i]->k))
                {
                    // current_mask has the bit corresponding to
                    // nbr_mask_vec[i] set - so and with the appropriate
                    // vector for a hit 
                    has_bit=false;
                    // Do the loop below manually rather than calling the overloaded AND
                    for(j=0; j<current_mask.S; j++)
                    {
                        if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                        {
                            // Some nbr is found in word j
                            has_bit=true;
                            break;
                        }
                    }
                    // breaks to here
                    if( has_bit )
                    {                
                        // Set represented by current_mask is not independent
                        is_independent=false;
                        int rpos=0;
                        // Find the rightmost bit in current_mask - this
                        // could be sped up by looking at the words 
                        // SPEEDUP
                        while(!(current_mask.test_bit(rpos)))
                            rpos++;

                        // Now advance current_mask by a potentially massive amount
                        // Would be faster here to figure out which word we are in and 
                        // do this w/o test_bit
                        int b=0;
                        while(current_mask.test_bit(rpos+b))
                            b++;
                        for(int a=0;a<b;a++)
                            current_mask.xor_bit(rpos+a);
                        current_mask.or_bit(rpos+b);
                        break;                   
                    }
                }
            }


            if(is_independent)
            {   
                // ind_counter keep track of how many independent sets
                // are in the indset_mask vector
                ind_counter ++;
                indset_masks.push_back(current_mask);
                ind_wgh.push_back(0);
                ++current_mask;
            }
        }


        // Once assigned work is finished, ask for more work
        if (!(counter < end_pos))
        {
            if (T->rank > 0)
            {
                MPI_Request maskrequest;
                MPI_Status maskstatus;
                int wflag = 0;
                int count = 0;
            
                //DEBUG("process : %d ask for work from : %d\n", T->rank, 0);
                MPI_Send(&(T->rank), 1, MPI_INT, 0, MPI_REQ_WRK_TAG, MPI_COMM_WORLD);
                MPI_Recv (&stwork, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_WRK_TAG, MPI_COMM_WORLD, &maskstatus);
                count = 0;
                MPI_Get_count (&maskstatus, MPI_UNSIGNED_LONG_LONG, &count);

                if (!(count > 0))
                    break;

                current_mask.zeroize ();
                end_mask.zeroize ();
                current_mask = *(start_masks[stwork]);
                end_mask = *(end_masks[stwork]);
            }
            else
            {
                //Sending terminate signal
                for (i = 1; i < T->size; i++)
                    MPI_Send (&end_pos, 0, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD);
                break;
            }
        }
    }

    // Now do the DP over the children for this independent set - if this is a leaf, then
    // the body of the loop is never executed
    int ind_set_tbl_size = (int)indset_masks.size();
    MPI_Request wgh_s_request;
    MPI_Request wgh_r_request;
    int imaskcount = 0;
    int request_size = T->request_size;
    MPI_Status wstatus[T->size];
    request_index = 0;
    request_flag = 0;
    int child_tbl_flag = 0;
    vector<int> *wgh_response;

    T->rcount->assign (T->size, 0);
    T->request_count->assign (T->size, -1);
    // If a child table is stored in current node, child_tbl_flag will
    // be set to 1.  Which means current node should be ready to serve
    // weight requests. 
    for (i = 0; i < children_size; i++)
    {
        if (T->rank == (*T->tbl_loc)[children[i]])
        {
            child_tbl_flag = 1;
            break;
        }
    }
    

    DEBUG("computing mask weights\n");
    int iteration = 0;
    DEBUG("mask count %d tble size %d\n", imaskcount, ind_set_tbl_size);
    if (!(ind_set_tbl_size > 0))
    {
        for(i = 0; i < (int)children_size; i++)
            MPI_Isend ((void *)NULL, 0, MPI_UNSIGNED_LONG_LONG, 
                       (*T->tbl_loc)[children[i]], MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &wgh_s_request);
    }

    //DEBUG("sending out requests\n");
    DEBUG("# of nbrs %d\n", children_size);
    DEBUG("child:location: ");
    for ( i = 0; i < children_size; i++)
    {
        GEN("%d:%d ", children[i], (*T->tbl_loc)[children[i]]);
    }
    GEN("\n");

    // Sending request of size request_size.
    int count_send = 0;
    int destination = 0;
    int msk_send_index = 0;
    vector<unsigned long long> *msk_send_buffer;

    while (imaskcount < ind_set_tbl_size)
    {
        for(i = 0; i < (int)children_size; i++)
        {
            msk_send_index = get_msk_send_buffer_index (T);
            msk_send_buffer =  T->msk_send_pool[msk_send_index];

            for (j = 0; (j < request_size) && (imaskcount + j < ind_set_tbl_size); j++)
            {
                temp_set.mask->zeroize();
                for(int z = 0; z < w; z++)
                {
                    // This part could be sped up by considering the words and then doing some xoring, I think
                    if( cvec[i][k_bag_vec[z]] )
                    {
                        if(indset_masks[imaskcount + j].test_bit(z))
                        {
                            // Bit z is set, and the corresponding vertex is also in the child
                            temp_set.mask->set_bit(z);                        
                        }
                            
                    }
                }

                for (kk = 0; kk < T->num_mask_words; kk++)
                    (*msk_send_buffer)[j*T->num_mask_words + kk] = temp_set.mask->words[kk];
            }
            // Find the location of the child (k_adj_vec[i + 1]) using
         // T->tbl_loc vector.Then send the masks to get the
         // weights, we do this asynchonously therefore we can move
         // forward with other requests too.  
         // Note: k_adj_vec[0] is the parent.

         // Tags:
         // Since there will be bunch of requests, different
         // requests are identified using the tag. 
         // Tag is used in following way,
         // i*MPI_REQ_WGH_TAG + iteration
         // where i is the position of the child in the bag and
         // iteration is the request number. This is important when
         // indset_masks table is larger than request_size. Because
         // there has to be multiple requests. 
         //DEBUG("request_size: %d j : %d iteration : %d\n",
         // T->request_size, j, iteration);
            (*msk_send_buffer)[j*T->num_mask_words] = children[i];
            (*msk_send_buffer)[j*T->num_mask_words + 1] = iteration;
            destination = (*T->tbl_loc)[children[i]];
            MPI_Isend ((void *)&(*msk_send_buffer).front(), j*T->num_mask_words + 2, MPI_UNSIGNED_LONG_LONG, 
                       destination, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_send_requests)[msk_send_index]);
        }
        iteration ++;
        imaskcount += request_size;

        if (child_tbl_flag)
            recv_msk_send_wgh(T);
    }

    DEBUG("processor : %d iterations %d \n", T->rank, iteration);
    int local_children = 0;

    // This extra messages will help to make sure that any request for
    // weights get lost during the communication. Whoever has table
    // should receive this many requests.
    for(i = 0; i < (int)children_size; i++)
    {
        MPI_Request sendmskcountreq;
        MPI_Isend ((void *)&iteration, 1, MPI_INT, 
                   (*T->tbl_loc)[children[i]], MPI_REQ_COUNT_TAG, MPI_COMM_WORLD, &sendmskcountreq);

        if ((*T->tbl_loc)[children[i]] == T->rank) local_children ++;
    }

    if (child_tbl_flag)
    {
        while (*(min_element(T->request_count->begin(), T->request_count->end())) < 0)
        {
            MPI_Status recvstat;
            int tmp = -1;
            MPI_Recv((void *)&tmp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_REQ_COUNT_TAG, MPI_COMM_WORLD, &recvstat);
            (*T->request_count)[recvstat.MPI_SOURCE] = tmp;
        }
        recv_msk_send_wgh(T);

        int rcomplete_flag = 1;
        while(1)
        {
            recv_msk_send_wgh(T);
            for (i = 0; i < T->size; i++)
            {
                DEBUG("i : %d rcount : %d request count : %d \n", i, (*T->rcount)[i], (*T->request_count)[i]);
                if (((*T->rcount)[i] < (*T->request_count)[i]*local_children))
                {
                    rcomplete_flag = 0;
                    break;
                }
                rcomplete_flag = 1;
            }
            
            if (rcomplete_flag) break;
        }
    }



    DEBUG("processor : %d has local children : %d\n", T->rank, local_children);
    // Current compute node ready to serve weight requests, 
    // Only if,
    // 1. If current tree node has children
    // 2. If current compute node has at least one child table 
    //MPI_Status status;
    request_flag = 0;

    // child_iterations vector make sure that, we have recieved
    // intended number of requests from compute nodes where children
    // tables are located.
    vector<int> child_iterations(T->size, 0);
    for (i = 0; i < children_size; i++)
    {
        child_iterations[(*T->tbl_loc)[children[i]]] += iteration;
    }

    // tbl_loc vector make sures that nodes recieves requests from
    // intended sources.
    vector<int> tbl_loc(T->size, 0);
    for (i = 0; i < children_size; i++)
        tbl_loc[(*T->tbl_loc)[children[i]]] = 1;

    vector<int> wgh_recv_buffer(request_size + 1, 0);
    T->iter_count->assign(T->size, 0);
    
    // Receiving weights for the masks sent over
    if (children_size*iteration > 0)
    {
        int wgh_recv_index;
        int wgh_recv_flag = 0;
        int start_pos = 0;
        int recv_count = 0;
        int iter = 0;
        j = 0;

        DEBUG("waiting\n");
        while (1)
        {

            MPI_Status wgh_recv_status;
            wgh_recv_index = 0;
            MPI_Waitany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_status);

            while(!tbl_loc[wgh_recv_index])
            {
              T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
              MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
                         wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);
              MPI_Waitany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_status);
            }
            
            recv_count = 0;
            MPI_Get_count (&wgh_recv_status, MPI_INT, &recv_count);

            wgh_recv_buffer = (*T->wgh_recv_pool[wgh_recv_index]);
            T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
            MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
                       wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);


            iter = wgh_recv_buffer[0];
            start_pos = request_size*iter;
            (*T->iter_count)[wgh_recv_index]++;

            for (i = 0; i < recv_count - 1; i++)
            {
                if (start_pos + i > ind_wgh.size())
                {

                    FERROR("location : %d is out of boundaries of the array: %d mask array size: %d \n", 
                          start_pos + i, ind_wgh.size(), indset_masks.size());
                    FERROR("source : %d iteration : %d recv_count: %d\n", wgh_recv_status.MPI_SOURCE, 
                          wgh_recv_buffer[0], recv_count);
                    MPI_Finalize();
                    exit (-1);
                }
                ind_wgh[start_pos + i] += wgh_recv_buffer[i + 1];
            }


            j = 1;
            for (i = 0; i < T->size; i++)
            {
                if ((*T->iter_count)[i] < child_iterations[i])
                {
                    j = 0;
                    break;
                }
            }
            if (j) break;

            if (child_tbl_flag)
                recv_msk_send_wgh(T);
        }
    }
    DEBUG("finished recv weights\n");
    // Done with child DP
    TDSolution *temp_sol;
    int parent_weight = 0;

    // Updating weight for intersection with the parent.
    int maxvalue = 0;
    vector<int> weight = T->G->get_weight();
    for (i = 0; i < ind_set_tbl_size; i++)
    {
        temp_set.mask->zeroize();
        parent_weight=0;
        current_mask = indset_masks[i];
    
        for(j = 0; j < w; j++)
        {                
            if(current_mask.test_bit(j))
            { 
                ind_wgh[i] += weight[k_bag_vec[j]];
                if(parent_location_vec[k_bag_vec[j]]!=-1)
                {
                    temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
                    parent_weight += weight[k_bag_vec[j]];
                }
            }
        }
        
        // temp_mask is a mask of this set \cap parent_bag
        hash_ptr=NULL;
        
        // k is not the root - so do the usual
        HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (k == T->root_node)
        {
            if (ind_wgh[i] > maxvalue)
            {
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
                *(added_set->mask)=*(temp_set.mask);
                *(added_set->orig_mask)=current_mask;
                
                added_set->value = ind_wgh[i];
                maxvalue = ind_wgh[i];
                
                HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
        else
        {
            
            if(hash_ptr)
            {
                if (ind_wgh[i] - parent_weight > hash_ptr->value)
                {
                    hash_ptr->value = ind_wgh[i] - parent_weight;
                }
            }
            else
            {
                // There is no entry for this mask - add it to the table
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
            *(added_set->mask)=*(temp_set.mask);
            *(added_set->orig_mask)=current_mask;
            
            added_set->value = ind_wgh[i] - parent_weight;
            
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
    }
    
    DEBUG("Hash tables ready to save\n");

    
    // Find a location to save the table.
    int minloc = -1;
    int minent = INT_MAX;
    for (i = 0; i < T->size; i++)
    {
        if ((*T->tbl_entries)[i] < minent)
        {
            minent = (*T->tbl_entries)[i];
            minloc = i;
        }
    }
    (*T->tbl_loc)[k] = minloc;
    DEBUG("table for node : %d will be saved in : %d\n", k, minloc);

    // Before reaching this point we have calculate weights for the
    // masks in indset_masks vector.
	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
    DEBUG("size of the table : %d\n", table_size);
    DEBUG("Time for DP loop : %lf\n", MPI_Wtime() - stime);

    hash_ptr = NULL;
    i = 0;
    int masksize = T->num_mask_words + 1;
    vector<unsigned long long> maskweight_array(table_size*(masksize));

    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        for (j = 0; j < T->num_mask_words; j++)
        {
            maskweight_array[i*(masksize) + j] = hash_ptr->mask->words[j];
        }
        maskweight_array[i*(masksize) + j] = hash_ptr->value;
        i++;
    }

    vector<int> n_elements_vec(T->size, 0);
    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    DEBUG("In gather phase\n");
    
    MPI_Allgather(&table_size, 1, MPI_INT, &n_elements_vec.front(), 1, MPI_INT, MPI_COMM_WORLD);
    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elements_vec[i];
        w_recv_counts[i] = n_elements_vec[i];
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

    (*T->tbl_entries)[minloc] += total_entries;
    (*T->node_sizes)[k] = total_entries;

    // Deriving datatype
    MPI_Datatype mask_type;
    MPI_Type_contiguous (T->num_mask_words + 1, MPI_UNSIGNED_LONG_LONG, &mask_type);
    MPI_Type_commit(&mask_type);
    vector<unsigned long long> maskweights_vec(total_entries*masksize, 0);    

    DEBUG("starting Gather v\n");
    MPI_Gatherv((void *)&maskweight_array.front(), w_recv_counts[T->rank], mask_type,
                (void *)&maskweights_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), mask_type, minloc, MPI_COMM_WORLD);

    unsigned long long max = 0;
    DEBUG("total entries : %d\n", total_entries);
    DEBUG("entering loop to update hash table \n");

    stime = MPI_Wtime ();
    // Loop to  update hash table after gathering values at other
    // processors 
    DEBUG("gather all results to %d\n", minloc);
    for (i = 0; i < total_entries && (T->rank == minloc); i++)
    {

       if (T->root_node == k)
        {
            unsigned long long cval = maskweights_vec[i*masksize + T->num_mask_words];
            if ( cval > max)
            {
                max = cval;
            }
            continue;
        }

        TDSolution *added_set;
        added_set=new TDSolution(T->num_mask_words);
        if (!added_set)
        {
            FERROR("no enough memory \n");
            MPI_Finalize();
            exit(-1);
        }

        for (j = 0; j < T->num_mask_words; j++)
        {
            added_set->mask->words[j] = maskweights_vec[i*masksize + j];
        }
        added_set->value= maskweights_vec[i*masksize + j];

        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
        }
        else
        {
            if (added_set->value > hash_ptr->value)
            {
                hash_ptr->value = added_set->value;
            }
        }
    }

    DEBUG("completed table size : %d\n", (int)HASH_COUNT(T->tree_nodes[k]->hash_table));
    if (T->root_node == k)
    {
        DEBUG("maximum value %d\n", max);
		// CSG adding just to check answer
        if (T->rank == minloc)
            fprintf(stderr,"Max value is %d\n",max);
    }
    DEBUG("Finish hash table update\n");
    DEBUG("Time for hash table update : %lf\n", MPI_Wtime() - stime);

    // clear some vectors to gain more memory.
     maskweights_vec.clear();
     maskweight_array.clear();
     w_recv_disp.clear();
     w_recv_counts.clear();
     indset_masks.clear();

    for (i = 0; i < children_size; i++)
    {
        (*T->tbl_entries)[(*T->tbl_loc)[children[i]]] -= (*T->node_sizes)[children[i]];
        (*T->node_sizes)[children[i]] = 0;
        (*T->tbl_loc)[children[i]] = -1;
    }

    
    for (i = 0; i < entries; i++)
        delete start_masks[i];

	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	
    // Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}


static int send_masks (TDTree *T, int k, int startpos, int iteration, 
                       vector<int> *k_bag_vec, vector<bool *> *cvec, vector<int> *children, vector<bigint_t> *indset_masks)
{
    int i;
    int j;
    int z;
    int kk;
    int ind_set_tbl_size = indset_masks->size();
    int destination = 0;
    int msk_send_index = 0;
    vector<unsigned long long> *msk_send_buffer;
    int children_size = children->size();
	TDSolution temp_set(T->num_mask_words);
    int w = T->tree_nodes[k]->bag.size();

    for(i = 0; i < (int)children_size; i++)
    {
        msk_send_index = get_msk_send_buffer_index (T);
        if (msk_send_index < 0)
            return msk_send_index;

        msk_send_buffer =  T->msk_send_pool[msk_send_index];
        
        for (j = startpos; (j < T->request_size) && (startpos + j < ind_set_tbl_size); j++)
        {
            temp_set.mask->zeroize();
            for(int z = 0; z < w; z++)
            {
                // This part could be sped up by considering the words and then doing some xoring, I think
                if( (*cvec)[i][(*k_bag_vec)[z]] )
                {
                    if((*indset_masks)[j].test_bit(z))
                    {
                        // Bit z is set, and the corresponding vertex is also in the child
                        temp_set.mask->set_bit(z);                        
                    }
                    
                }
            }
            
            for (kk = 0; kk < T->num_mask_words; kk++)
                (*msk_send_buffer)[j*T->num_mask_words + kk] = temp_set.mask->words[kk];
        }
            // Find the location of the child (k_adj_vec[i + 1]) using
        // T->tbl_loc vector.Then send the masks to get the
        // weights, we do this asynchonously therefore we can move
        // forward with other requests too.  
        // Note: k_adj_vec[0] is the parent.
        
        // Tags:
        // Since there will be bunch of requests, different
        // requests are identified using the tag. 
        // Tag is used in following way,
        // i*MPI_REQ_WGH_TAG + iteration
        // where i is the position of the child in the bag and
        // iteration is the request number. This is important when
        // indset_masks table is larger than request_size. Because
        // there has to be multiple requests. 
        //DEBUG("request_size: %d j : %d iteration : %d\n",
        // T->request_size, j, iteration);
        (*msk_send_buffer)[j*T->num_mask_words] = (*children)[i];
        (*msk_send_buffer)[j*T->num_mask_words + 1] = iteration;
        destination = (*T->tbl_loc)[(*children)[i]];
        
        MPI_Isend ((void *)&(*msk_send_buffer).front(), j*T->num_mask_words + 2, MPI_UNSIGNED_LONG_LONG, 
                   destination, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_send_requests)[msk_send_index]);
    }
    return 1;
}

static int send_masks_9 (TDTree *T, int k, int startpos, int iteration, 
                       vector<int> *k_bag_vec, vector<bool *> *cvec, vector<int> *children)
{
    int i;
    int j;
    int z;
    int kk;
    int ind_set_tbl_size = T->indset_masks->size();
    int destination = 0;
    int msk_send_index = 0;
    vector<unsigned long long> *msk_send_buffer;
    int children_size = children->size();
	TDSolution temp_set(T->num_mask_words);
    int w = T->tree_nodes[k]->bag.size();
	TDSolution *hash_ptr=NULL;

    //DEBUG("startpos : %d iteration : %d\n", startpos, iteration);
    for(i = 0; i < (int)children_size; i++)
    {
        destination = (*T->tbl_loc)[(*children)[i]];
        msk_send_buffer =  T->msk_send_pool[destination];
        for (j = 0; (j < T->request_size) && (startpos + j < ind_set_tbl_size); j++)
        {
            temp_set.mask->zeroize();
            for(int z = 0; z < w; z++)
            {
                // This part could be sped up by considering the words and then doing some xoring, I think
                if( (*cvec)[i][(*k_bag_vec)[z]] )
                {
                    if((*T->indset_masks)[startpos + j].test_bit(z))
                    {
                        // Bit z is set, and the corresponding vertex is also in the child
                        temp_set.mask->set_bit(z);                        
                    }
                       
                }
            }
               
            for (kk = 0; kk < T->num_mask_words; kk++)
                (*msk_send_buffer)[j*T->num_mask_words + kk] = temp_set.mask->words[kk];
        }
        //DEBUG("j : %d startpos + j : %d ind_set_tbl_size : %d\n", j, startpos + j, ind_set_tbl_size);
            
        // Find the location of the child (k_adj_vec[i + 1]) using
        // T->tbl_loc vector.Then send the masks to get the
        // weights, we do this asynchonously therefore we can move
        // forward with other requests too.  
        // Note: k_adj_vec[0] is the parent.
        
        // Tags:
        // Since there will be bunch of requests, different
        // requests are identified using the tag. 
        // Tag is used in following way,
        // i*MPI_REQ_WGH_TAG + iteration
        // where i is the position of the child in the bag and
        // iteration is the request number. This is important when
        // indset_masks table is larger than request_size. Because
        // there has to be multiple requests. 
        //DEBUG("request_size: %d j : %d iteration : %d\n",
        // T->request_size, j, iteration);
        //(*msk_send_buffer)[j*T->num_mask_words] =
        //(*children)[i];
        (*msk_send_buffer)[j*T->num_mask_words] = (*children)[i];
        (*msk_send_buffer)[j*T->num_mask_words + 1] = iteration;
        if (MPI_SUCCESS != 
            MPI_Send ((void *)&(*msk_send_buffer).front(), j*T->num_mask_words + 2, MPI_UNSIGNED_LONG_LONG, 
                      destination, MPI_REQ_WGH_TAG, MPI_COMM_WORLD))
        {
            FERROR("MPI Send failed in processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
        //DEBUG("Successfully send %d elements to : %d\n", j*T->num_mask_words + 2, destination);
    }
    return 1;
}


static int recv_wgh(TDTree *T, vector<int> *children, vector<int> *ind_wgh)
{
    int wgh_recv_index;
    int wgh_recv_flag = 0;
    int start_pos = 0;
    int recv_count = 0;
    int iter = 0;
    int j = 0;
    int i = 0;
    int request_size = T->request_size;
    vector<int> wgh_recv_buffer(request_size + 1, 0);

    MPI_Status wgh_recv_status;
    wgh_recv_index = 0;
    MPI_Testany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_flag, &wgh_recv_status);
    if (wgh_recv_flag)
    {

//     while(!tbl_loc[wgh_recv_index])
//     {
//         T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
//         MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
//                    wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);
//         MPI_Waitany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_status);
//     }
            
        recv_count = 0;
        MPI_Get_count (&wgh_recv_status, MPI_INT, &recv_count);

        wgh_recv_buffer = (*T->wgh_recv_pool[wgh_recv_index]);
        T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 1, 0);
        MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 1, MPI_INT, 
                   wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]);


        iter = wgh_recv_buffer[0];
        start_pos = request_size*iter;
        (*T->iter_count)[wgh_recv_index]++;

        for (i = 0; i < recv_count - 1; i++)
        {
            if (start_pos + i > ind_wgh->size())
            {

                FERROR("location : %d is out of boundaries of the array: %d \n",
                      start_pos + i, ind_wgh->size());
                FERROR("source : %d iteration : %d recv_count: %d\n", wgh_recv_status.MPI_SOURCE, 
                      wgh_recv_buffer[0], recv_count);
                MPI_Finalize();
                exit (-1);
            }
            (*ind_wgh)[start_pos + i] += wgh_recv_buffer[i + 1];
        }
    }
    else
    {
        return -1;
    }
    return 1;

}

static bool recv_wgh_9 (TDTree *T, int k, vector<int> *children)
{
    int wgh_recv_index;
    int wgh_recv_flag = 0;
    int start_pos = 0;
    int recv_count = 0;
    int iter = 0;
    int j = 0;
    int i = 0;
    int request_size = T->request_size;
    int tree_node = 0;
    int received_tree_node = 0;

    vector<int> wgh_recv_buffer(request_size + 2, 0);

    MPI_Status wgh_recv_status;
    wgh_recv_index = 0;

    if (MPI_SUCCESS != 
        MPI_Testany (T->size, &T->wgh_recv_requests->front(), &wgh_recv_index, &wgh_recv_flag, &wgh_recv_status))
    {
        FERROR ("MPI Test any failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
        return false;
    }

    if (!wgh_recv_flag)
        return false;

    recv_count = 0;
    if (MPI_SUCCESS != 
        MPI_Get_count (&wgh_recv_status, MPI_INT, &recv_count))
    {
        FERROR ("MPI Get count failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
        return false;
    }

    //DEBUG("received : %d from : %d\n", recv_count, wgh_recv_status.MPI_SOURCE);
    wgh_recv_buffer = (*T->wgh_recv_pool[wgh_recv_index]);

    T->wgh_recv_pool[wgh_recv_index]->assign(request_size + 2, 0);
    if (MPI_SUCCESS != 
        MPI_Irecv ((void *)&T->wgh_recv_pool[wgh_recv_index]->front(), request_size + 2, MPI_INT, 
                   wgh_recv_index, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[wgh_recv_index]))
    {
        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
        return false;
    }

    iter = wgh_recv_buffer[0];
    tree_node = wgh_recv_buffer[1];
    start_pos = request_size*iter;

    received_tree_node = T->tree_nodes[tree_node]->adj.front();
    if (k != received_tree_node)
    {
        WARN ("expected weights for tree node : %d but received : %d\n", k, received_tree_node);
        WARN ("source: %d child node: %d iteration : %d elements : %d\n", wgh_recv_index, tree_node, iter, recv_count);
        return false;
    }

    for (i = 0; i < T->storage_nodes_size; i++)
    {
        if (T->storage_nodes[i] == wgh_recv_index)
            (*T->request_count)[i]++;
    }


    for (i = 0; i < recv_count - 2; i++)
    {
        if (start_pos + i > T->ind_wgh->size())
        {

            FERROR("location : %d is out of boundaries of the array: %d \n",
                  start_pos + i, T->ind_wgh->size());
            FERROR("source : %d iteration : %d recv_count: %d\n", wgh_recv_status.MPI_SOURCE, 
                  wgh_recv_buffer[0], recv_count);
            MPI_Finalize();
            return false;
        }
        (*T->ind_wgh)[start_pos + i] += wgh_recv_buffer[i + 2];
    }

    return true;
}

void parallel_wis_init_9 (TDTree *T, int size, int rank)
{
    int i = 0;

    T->size = size;
    T->rank = rank;
    T->request_size = 2048;
    T->msk_send_count = 0;
    T->wgh_send_count = 0;
    T->entries = 8192;

    if (!(T->size > 0))
    {
        FERROR("Environment size can not be zero\n");
        MPI_Finalize();
        exit (-1);
    }

    if (T->rank < 0)
    {
        FERROR("Ranks can not be negative\n");
        MPI_Finalize();
        exit (-1);
    }

    try 
    {

        T->start_masks = new vector<bigint_t *> (T->entries);
        T->end_masks = new vector<bigint_t *> (T->entries);
        for (i = 0; i < T->entries; i++)
        {
            (*T->start_masks)[i] = new bigint_t (T->num_mask_words);
            (*T->end_masks)[i] = new bigint_t (T->num_mask_words);
        }


        T->tbl_entries = new vector<int> (T->size,0);
        T->tbl_loc = new vector<int> (T->num_tree_nodes, -1);
        T->node_sizes = new vector<int> (T->num_tree_nodes, 0);
        T->iter_count = new vector<int> (T->size, 0);
        //T->rcount = new vector<int> (T->size, 0);
        // request_count vector keep track of # of requests each node has
        // to receive.
        T->request_count = NULL;

    
        T->requests = new vector<MPI_Request>(T->size);
        if (T->rank == 0)
        {
            for (i = 0; i < T->size; i++)
            {
                MPI_Irecv (&(T->requester_rank), 1, MPI_INT, i, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, &((*T->requests)[i]));
            }
        }

        // Allocating buffers to receive weight requests. Masks are sent
        // over the network along with the child number and iteration #.
        T->msk_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->msk_recv_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            MPI_Irecv((void *)&(*T->msk_recv_pool[i]).front(), T->request_size*T->num_mask_words + 2, 
                      MPI_UNSIGNED_LONG_LONG, i, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_recv_requests)[i]);
        }


        // Response with weights for masks comes with one additional
        // element, iteration #
        // Iteration # is used to find the correct position in original
        // vector for the weights received. 
        T->wgh_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->wgh_recv_pool.push_back(new vector<int>(T->request_size + 2, 0));
            MPI_Irecv ((void *)&T->wgh_recv_pool[i]->front(), T->request_size + 2, 
                       MPI_INT, i, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[i]);
        }

        T->msk_send_requests = new vector<MPI_Request>(T->size);
        T->wgh_send_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->msk_send_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            T->wgh_send_pool.push_back(new vector<int>(T->request_size + 2, 0));
        }
    }
    catch (std::bad_alloc e)
    {
        cerr << e.what() << endl;
        FERROR("Bad allocation: %s\n", e.what());
        MPI_Finalize();
        exit (-1);
    }
}

void parallel_wis_init_10 (TDTree *T, int size, int rank)
{
    int i = 0;

    T->size = size;
    T->rank = rank;
    T->request_size = 2048;
    T->msk_send_count = 0;
    T->wgh_send_count = 0;
    T->entries = 8192;

    if (!(T->size > 0))
    {
        FERROR("Environment size can not be zero\n");
        MPI_Finalize();
        exit (-1);
    }

    if (T->rank < 0)
    {
        FERROR("Ranks can not be negative\n");
        MPI_Finalize();
        exit (-1);
    }

    try 
    {

        T->start_masks = new vector<bigint_t *> (T->entries);
        T->end_masks = new vector<bigint_t *> (T->entries);
        for (i = 0; i < T->entries; i++)
        {
            (*T->start_masks)[i] = new bigint_t (T->num_mask_words);
            (*T->end_masks)[i] = new bigint_t (T->num_mask_words);
        }


        T->tbl_entries = new vector<int> (T->size,0);
        T->tbl_loc = new vector<int> (T->num_tree_nodes, -1);
        T->node_sizes = new vector<int> (T->num_tree_nodes, 0);
        T->iter_count = new vector<int> (T->size, 0);
        //T->rcount = new vector<int> (T->size, 0);
        // request_count vector keep track of # of requests each node has
        // to receive.
        T->request_count = NULL;

    
        T->requests = new vector<MPI_Request>(T->size);
        if (T->rank == 0)
        {
            for (i = 0; i < T->size; i++)
            {
                MPI_Irecv (&(T->requester_rank), 1, MPI_INT, i, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, &((*T->requests)[i]));
            }
        }

        // Allocating buffers to receive weight requests. Masks are sent
        // over the network along with the child number and iteration #.
        T->msk_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->msk_recv_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            MPI_Irecv((void *)&(*T->msk_recv_pool[i]).front(), T->request_size*T->num_mask_words + 2, 
                      MPI_UNSIGNED_LONG_LONG, i, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_recv_requests)[i]);
        }


        // Response with weights for masks comes with one additional
        // element, iteration #
        // Iteration # is used to find the correct position in original
        // vector for the weights received. 
        T->wgh_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->wgh_recv_pool.push_back(new vector<int>(T->request_size + 2, 0));
            MPI_Irecv ((void *)&T->wgh_recv_pool[i]->front(), T->request_size + 2, 
                       MPI_INT, i, MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[i]);
        }

        T->msk_send_requests = new vector<MPI_Request>(T->size);
        T->wgh_send_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->msk_send_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
            T->wgh_send_pool.push_back(new vector<int>(T->request_size + 2, 0));
        }
    }
    catch (std::bad_alloc e)
    {
        cerr << e.what() << endl;
        FERROR("Bad allocation: %s\n", e.what());
        MPI_Finalize();
        exit (-1);
    }
}

void parallel_wis_cleanup_9 (TDTree *T)
{
    int i;
     // joining threads
    DEBUG("join threads\n");

    // Cancel waiting requests.
    for (i = 0; i < T->size; i++)
    {
        if (T->rank == 0)
        {
            MPI_Cancel (&(*T->requests)[i]);
        }
        MPI_Cancel(&(*T->wgh_recv_requests)[i]);
        MPI_Cancel(&(*T->msk_recv_requests)[i]);
    }

    delete T->tbl_entries;
    delete T->tbl_loc;
    delete T->node_sizes;
    delete T->iter_count;
    delete T->requests;
    delete T->start_masks;
    delete T->end_masks;


    delete T->msk_recv_requests;
    for (i = 0; i < T->size; i++)
    {
        delete T->msk_recv_pool[i];
        delete T->wgh_recv_pool[i];
        
    }
    delete T->wgh_recv_requests;

    delete T->msk_send_requests;
    delete T->wgh_send_requests;

    for (i = 0; i < T->size; i++)
    {
        delete  T->msk_send_pool[i];
        delete  T->wgh_send_pool[i];
    }
}

void parallel_wis_cleanup_10 (TDTree *T)
{
    int i;
     // joining threads
    DEBUG("join threads\n");

    // Cancel waiting requests.
    for (i = 0; i < T->size; i++)
    {
        if (T->rank == 0)
        {
            MPI_Cancel (&(*T->requests)[i]);
        }
        MPI_Cancel(&(*T->wgh_recv_requests)[i]);
        MPI_Cancel(&(*T->msk_recv_requests)[i]);
    }

    delete T->tbl_entries;
    delete T->tbl_loc;
    delete T->node_sizes;
    delete T->iter_count;
    delete T->requests;
    delete T->start_masks;
    delete T->end_masks;


    delete T->msk_recv_requests;
    for (i = 0; i < T->size; i++)
    {
        delete T->msk_recv_pool[i];
        delete T->wgh_recv_pool[i];
        
    }
    delete T->wgh_recv_requests;

    delete T->msk_send_requests;
    delete T->wgh_send_requests;

    for (i = 0; i < T->size; i++)
    {
        delete  T->msk_send_pool[i];
        delete  T->wgh_send_pool[i];
    }
}
static int compute_nonnice_table_parallel_9 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
    list<int>::iterator it_li;
    int kk;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());
    int children_size = (int) children.size();
    int child_tbl_flag = 0;
    // If a child table is stored in current node, child_tbl_flag will
    // be set to 1.  Which means current node should be ready to serve
    // weight requests. 
    for (i = 0; i < children_size; i++)
    {
        if (T->rank == (*T->tbl_loc)[children[i]])
        {
            child_tbl_flag = 1;
            break;
        }
    }


    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    DEBUG("mask generation\n");
    parallel_mask_gen (T, k, w, T->start_masks, T->end_masks, T->entries);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;

    vector<bigint_t> indset_masks;
    vector<int> ind_wgh;
    unsigned long long ind_counter = 0;
    
    T->indset_masks = &indset_masks;
    T->ind_wgh = &ind_wgh;

 
	bigint_t current_mask(T->num_mask_words);
	bigint_t end_mask(T->num_mask_words);

	TDSolution temp_set(T->num_mask_words);

    double stime = MPI_Wtime();
    int request_flag = 0;
    int request_index = -1;
    int last_send_pos = 0;

    unsigned long long start_pos, end_pos, wchunk, counter, stwork;
    wchunk = T->entries/T->size;
    start_pos = (wchunk * T->rank);
    end_pos = (wchunk * (T->rank + 1));

    vector<int> requests(T->size, 0);
    vector<int> terminate(T->size, -1);

    if (children_size > 0)
        T->request_count = new vector<int> (children_size, 0);

    int iteration = 0;

    DEBUG("DP loop\n");
    counter = start_pos;
    while (1)
    {
        
        // DEBUG("counter %d end_pos %d\n", counter, end_pos);
        if (counter < end_pos)
        {
            current_mask = *(*T->start_masks)[counter];
            end_mask = *(*T->end_masks)[counter];
            counter ++;
        }

        while(current_mask < end_mask)
        {

            if ((T->rank == 0) && (counter < end_pos))
            {

                if (MPI_SUCCESS != 
                    MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE))
                {
                    FERROR ("MPI Testany failed processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }

                while (request_flag)
                {
                    if (MPI_SUCCESS != 
                        MPI_Irecv (&(T->requester_rank), 1, MPI_INT, request_index, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, 
                                   &((*T->requests)[request_index])))
                    {
                        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize ();
                    }

                    requests[request_index] = 1;
                    
                    request_flag = 0;
                    request_index = -1;
                    T->requester_rank = -1;
                    if (MPI_SUCCESS != 
                        MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE))
                    {
                        FERROR ("MPI Testany failed processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }
                }
                
                // Send work to the requester.
                for (i = 1; i < T->size && (end_pos > counter); i++)
                {
                    if (requests[i])
                    {
                        end_pos --;
                        if (MPI_SUCCESS != 
                            MPI_Send (&end_pos, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD))
                        {
                            FERROR("MPI Send failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                            MPI_Finalize();
                        }
                        requests[i] = 0;
                    }
                }
            }
            
            is_independent=true;

            for(i = 0; i < w ; i++) 
            {            
                if(current_mask.test_bit(nbr_mask_vec[i]->k))
                {
                    // current_mask has the bit corresponding to
                    // nbr_mask_vec[i] set - so and with the appropriate
                    // vector for a hit 
                    has_bit=false;
                    // Do the loop below manually rather than calling the overloaded AND
                    for(j=0; j<current_mask.S; j++)
                    {
                        if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                        {
                            // Some nbr is found in word j
                            has_bit=true;
                            break;
                        }
                    }
                    // breaks to here
                    if( has_bit )
                    {                
                        // Set represented by current_mask is not independent
                        is_independent=false;
                        int rpos=0;
                        // Find the rightmost bit in current_mask - this
                        // could be sped up by looking at the words 
                        // SPEEDUP
                        while(!(current_mask.test_bit(rpos)))
                            rpos++;

                        // Now advance current_mask by a potentially massive amount
                        // Would be faster here to figure out which word we are in and 
                        // do this w/o test_bit
                        int b=0;
                        while(current_mask.test_bit(rpos+b))
                            b++;
                        for(int a=0;a<b;a++)
                            current_mask.xor_bit(rpos+a);
                        current_mask.or_bit(rpos+b);
                        break;                   
                    }
                }
            }


            if(is_independent)
            {   
                // ind_counter keep track of how many independent sets
                // are in the indset_mask vector
                ind_counter ++;
                T->indset_masks->push_back(current_mask);
                T->ind_wgh->push_back(0);
                ++current_mask;
                T->info->total_table_entries++;

                if ((ind_counter - last_send_pos) >= T->request_size)
                {
                    DEBUG("send iterations inside loop\n");
                    if (send_masks_9 (T, k, last_send_pos, iteration, &k_bag_vec, &cvec, &children) > 0)
                    {
                        last_send_pos = ind_counter;
                        iteration ++;
                    }
                    else
                    {
                        FERROR("no free buffers to send\n");
                    }
                }
            }
        }


        // Once assigned work is finished, ask for more work
        if (!(counter < end_pos))
        {
            if (T->rank > 0)
            {
                MPI_Request maskrequest;
                MPI_Status maskstatus;
                int wflag = 0;
                int count = 0;
            
                //DEBUG("process : %d ask for work from : %d\n", T->rank, 0);
                if (MPI_SUCCESS !=
                    MPI_Send(&(T->rank), 1, MPI_INT, 0, MPI_REQ_WRK_TAG, MPI_COMM_WORLD))
                {
                    FERROR ("MPI Send failed at processor :%d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize ();
                }

                if (MPI_SUCCESS != 
                    MPI_Recv (&stwork, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_WRK_TAG, MPI_COMM_WORLD, &maskstatus))
                {
                    FERROR ("MPI Recv failed at processor :%d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize ();
                }

                count = 0;
                if (MPI_SUCCESS != 
                    MPI_Get_count (&maskstatus, MPI_UNSIGNED_LONG_LONG, &count))
                {
                    FERROR ("MPI Get count failed at processor :%d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize ();
                }

                if (!(count > 0))
                    break;

                current_mask.zeroize ();
                end_mask.zeroize ();
                current_mask = *(*T->start_masks)[stwork];
                end_mask = *(*T->end_masks)[stwork];
            }
            else
            {
                //Sending terminate signal
                DEBUG("Send termination signal\n");
                for (i = 1; i < T->size; i++)
                {
                    if (MPI_SUCCESS != 
                        MPI_Send (&end_pos, 0, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD))
                    {
                        FERROR ("MPI Send failed in processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }
                }
                break;
            }
        }
    }

    DEBUG("iterations counter: %d\n", iteration);
    DEBUG("ind_counter: %d last_send_pos: %d\n", ind_counter, last_send_pos);
    if (ind_counter - last_send_pos > 0)
    {
        DEBUG("send iterations outside loop: %d\n", iteration);
        if (send_masks_9(T, k, last_send_pos, iteration, &k_bag_vec, &cvec, &children) > 0)
        {
            last_send_pos = ind_counter;
            iteration ++;
        }
    }

    
    // Now do the DP over the children for this independent set - if this is a leaf, then
    // the body of the loop is never executed
    int ind_set_tbl_size = (int)T->indset_masks->size();
    MPI_Request wgh_s_request;
    MPI_Request wgh_r_request;
    int imaskcount = 0;
    int request_size = T->request_size;
    MPI_Status wstatus[T->size];
    request_index = 0;
    request_flag = 0;
    vector<int> *wgh_response;


    DEBUG("computing mask weights\n");
    DEBUG("mask count %d tble size %d\n", imaskcount, ind_set_tbl_size);
    if (!(ind_set_tbl_size > 0))
    {
       vector<unsigned long long> *msk_send_buffer;
       int destination = 0;
       for(i = 0; i < (int)children_size; i++)
       {
           destination = (*T->tbl_loc)[children[i]];
           msk_send_buffer =  T->msk_send_pool[destination];
           
           if (destination != T->rank)
           {
               (*msk_send_buffer)[0] = children[i];
               (*msk_send_buffer)[1] = iteration;

               if (MPI_SUCCESS != 
                   MPI_Send ((void *)&(*msk_send_buffer).front(), 2, MPI_UNSIGNED_LONG_LONG, 
                             destination, MPI_REQ_WGH_TAG, MPI_COMM_WORLD))
               {
                   FERROR("MPI Send failed in processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                   MPI_Finalize();
               }
           }
           else
           {
               (*T->request_count)[i] = iteration;
           }
        }
    }

    //DEBUG("sending out requests\n");
    DEBUG("# of nbrs %d\n", children_size);
    DEBUG("child:location: ");
    for ( i = 0; i < children_size; i++)
    {
        GEN("%d:%d ", children[i], (*T->tbl_loc)[children[i]]);
    }
    GEN("\n");


    int local_children = 0;

    // This extra messages will help to make sure that any request for
    // weights get lost during the communication. Whoever has table
    // should receive this many requests.
    //iteration_counter --;
    //if (iteration_counter > 0)
    //iteration = iteration_counter;

    DEBUG("processor : %d iterations %d \n", T->rank, iteration);
    
    int it_len = (int)ceil(((double)iteration)/BIGINT_WORD_SIZE);    
    if (it_len > 0)
    {
        T->iteration_counter = new bigint_t (it_len);
    }

    int complete_flag = 0;
    if (child_tbl_flag)
    {
         
    }

    while (1)
    {
        if (!complete_flag)
        {
            while (recv_wgh_9 (T, k, &children));

            complete_flag = 1;
            for (i = 0; i < iteration; i++)
            {
                if (!T->iteration_counter->test_bit(i))
                {
                    DEBUG("Did not receive :%d\n", i);
                    complete_flag = 0;
                    break;
                }
            }
        }

        while(recv_msk_send_wgh_9 (T));
        
        if (complete_flag)
            break;
    }

    
    DEBUG("finished recv weights\n");
    // Done with child DP
    TDSolution *temp_sol;
    int parent_weight = 0;

    // Updating weight for intersection with the parent.
    int maxvalue = 0;
    vector<int> weight = T->G->get_weight();
    for (i = 0; i < ind_set_tbl_size; i++)
    {
        temp_set.mask->zeroize();
        parent_weight=0;
        current_mask = (*T->indset_masks)[i];
    
        for(j = 0; j < w; j++)
        {                
            if(current_mask.test_bit(j))
            { 
                (*T->ind_wgh)[i] += weight[k_bag_vec[j]];
                if(parent_location_vec[k_bag_vec[j]]!=-1)
                {
                    temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
                    parent_weight += weight[k_bag_vec[j]];
                }
            }
        }
        
        // temp_mask is a mask of this set \cap parent_bag
        hash_ptr=NULL;
        
        // k is not the root - so do the usual
        HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (k == T->root_node)
        {
            if ((*T->ind_wgh)[i] > maxvalue)
            {
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
                *(added_set->mask)=*(temp_set.mask);
                *(added_set->orig_mask)=current_mask;
                
                added_set->value = (*T->ind_wgh)[i];
                maxvalue = (*T->ind_wgh)[i];
                
                HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
        else
        {
            
            if(hash_ptr)
            {
                if ((*T->ind_wgh)[i] - parent_weight > hash_ptr->value)
                {
                    hash_ptr->value = (*T->ind_wgh)[i] - parent_weight;
                }
            }
            else
            {
                // There is no entry for this mask - add it to the table
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
            *(added_set->mask)=*(temp_set.mask);
            *(added_set->orig_mask)=current_mask;
            
            added_set->value = (*T->ind_wgh)[i] - parent_weight;
            
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
    }
    
    DEBUG("Hash tables ready to save\n");
    
    // Find a location to save the table.
    int minloc = -1;
    int minent = INT_MAX;
    for (i = (T->size - 1); i >= 0; i--)
    {
        if ((*T->tbl_entries)[i] < minent)
        {
            minent = (*T->tbl_entries)[i];
            minloc = i;
        }
    }
    (*T->tbl_loc)[k] = minloc;
    DEBUG("table for node : %d will be saved in : %d\n", k, minloc);

    // Before reaching this point we have calculate weights for the
    // masks in indset_masks vector.
	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
	T->info->total_pc_table_entries+=(unsigned long long)table_size;
 
    DEBUG("size of the table : %d\n", table_size);
    DEBUG("Time for DP loop : %lf\n", MPI_Wtime() - stime);

    hash_ptr = NULL;
    i = 0;
    int masksize = T->num_mask_words + 1;
    vector<unsigned long long> maskweight_array(table_size*(masksize));

    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        for (j = 0; j < T->num_mask_words; j++)
        {
            maskweight_array[i*(masksize) + j] = hash_ptr->mask->words[j];
        }
        maskweight_array[i*(masksize) + j] = hash_ptr->value;
        i++;
    }

    vector<int> n_elements_vec(T->size, 0);
    vector<int> w_recv_counts (T->size, 0);
    vector<int> w_recv_disp (T->size, 0);

    DEBUG("In gather phase\n");
    
    if (MPI_SUCCESS != 
        MPI_Allgather(&table_size, 1, MPI_INT, &n_elements_vec.front(), 1, MPI_INT, MPI_COMM_WORLD))
    {
        FERROR ("MPI Allgather failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
    }

    int total_entries = 0;
    for (i = 0; i < T->size; i++)
    {
        total_entries += n_elements_vec[i];
        w_recv_counts[i] = n_elements_vec[i];
        if (i > 0)
        {
            w_recv_disp[i] = w_recv_disp[i - 1] + w_recv_counts[i - 1];
        }
    }

    (*T->tbl_entries)[minloc] += total_entries;
    (*T->node_sizes)[k] = total_entries;

    // Deriving datatype
    MPI_Datatype mask_type;
    if (MPI_SUCCESS != 
        MPI_Type_contiguous (T->num_mask_words + 1, MPI_UNSIGNED_LONG_LONG, &mask_type))
    {
        FERROR ("MPI Type contiguous failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
    }

    if (MPI_SUCCESS != 
        MPI_Type_commit(&mask_type))
    {
        FERROR ("MPI Type commit failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
    }


    while(recv_msk_send_wgh_9 (T));

    vector<unsigned long long> maskweights_vec(total_entries*masksize, 0);    

    DEBUG("starting Gather v\n");
    if (MPI_SUCCESS != 
        MPI_Gatherv((void *)&maskweight_array.front(), w_recv_counts[T->rank], mask_type, 
                    (void *)&maskweights_vec.front(), &w_recv_counts.front(), &w_recv_disp.front(), mask_type, minloc, MPI_COMM_WORLD))
    {
        FERROR ("MPI Gatherv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
    }

    unsigned long long max = 0;
    DEBUG("total entries : %d\n", total_entries);
    DEBUG("entering loop to update hash table \n");

    stime = MPI_Wtime ();
    // Loop to  update hash table after gathering values at other
    // processors 
    DEBUG("gather all results to %d\n", minloc);
    for (i = 0; i < total_entries && (T->rank == minloc); i++)
    {

       if (T->root_node == k)
        {
            unsigned long long cval = maskweights_vec[i*masksize + T->num_mask_words];
            if ( cval > max)
            {
                max = cval;
            }
            continue;
        }

        TDSolution *added_set;
        added_set=new TDSolution(T->num_mask_words);
        if (!added_set)
        {
            FERROR("no enough memory \n");
            MPI_Finalize();
            exit(-1);
        }

        for (j = 0; j < T->num_mask_words; j++)
        {
            added_set->mask->words[j] = maskweights_vec[i*masksize + j];
        }
        added_set->value= maskweights_vec[i*masksize + j];

        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (!hash_ptr)
        {
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
        }
        else
        {
            if (added_set->value > hash_ptr->value)
            {
                hash_ptr->value = added_set->value;
            }
        }
    }

    DEBUG("completed table size : %d\n", (int)HASH_COUNT(T->tree_nodes[k]->hash_table));
    if (T->root_node == k)
    {
        DEBUG("maximum value %d\n", max);
		// CSG adding just to check answer
        if (T->rank == minloc)
            fprintf(stderr,"Max value is %d\n",max);
    }
    DEBUG("Finish hash table update\n");
    DEBUG("Time for hash table update : %lf\n", MPI_Wtime() - stime);

    // clear some vectors to gain more memory.
    maskweights_vec.clear();
    maskweight_array.clear();
    w_recv_disp.clear();
    w_recv_counts.clear();
    T->indset_masks->clear();
    T->ind_wgh->clear();
    T->request_count->clear();

    for (i = 0; i < children_size; i++)
    {
        (*T->tbl_entries)[(*T->tbl_loc)[children[i]]] -= (*T->node_sizes)[children[i]];
        (*T->node_sizes)[children[i]] = 0;
        (*T->tbl_loc)[children[i]] = -1;
    }


	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	
    // Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}


static int compute_nonnice_table_parallel_10 (TDTree *T, int k)
{
    DEBUG("In function : %s \n", __FUNCTION__);
    list<int>::iterator it_li;
    int kk;

    int i, j;
    int w, pw;
    int coi, coi_size;
    int parent = T->tree_nodes[k]->adj.front();

    // converting parent bag and node[k] bag into vectors for easy
    // accessiblity of elements.
    vector<int> parent_bag_vec (T->tree_nodes[parent]->bag.begin(), T->tree_nodes[parent]->bag.end());
    vector<int> k_bag_vec (T->tree_nodes[k]->bag.begin(), T->tree_nodes[k]->bag.end());
    vector<int> parent_location_vec(T->G->get_num_nodes(), -1);

    w = (int)k_bag_vec.size ();
    pw = (int)parent_bag_vec.size();

    for (i = 0; i < pw; i++)
    {
        parent_location_vec[parent_bag_vec[i]] = i;
    }

    // Create a list of neighbor masks (effectively creating a symmetric adj. matrix)
    vector<int_bigint *> nbr_mask_vec(w);
    vector<Graph::Node> nodes = T->G->get_nodes();
    list<int> nbrs;

    for(i = 0; i < w; i++)
    {
        nbr_mask_vec[i]=new int_bigint(T->num_mask_words);
        nbr_mask_vec[i]->k=i;
        for(j = 0;j < w; j++)
        {
            // Consider all poss. edges bag_vec[i]-->bag_vec[j]
            nbrs = nodes[k_bag_vec[i]].get_nbrs();
            for(it_li = nbrs.begin(); it_li != nbrs.end(); ++it_li)
            {
                if((i != j) && (k_bag_vec[j] == *it_li))
                {       
                    nbr_mask_vec[i]->w->or_bit(j);
                    break;
                }
            }
        }
    }

    vector<bool *> cvec(T->tree_nodes[k]->adj.size() - 1);
    vector<int> children(T->tree_nodes[k]->adj.begin(), T->tree_nodes[k]->adj.end());
    children.erase(children.begin());
    int children_size = (int) children.size();

    DEBUG("children size: %d\n", children_size);

    for (i = 0; i < T->tree_nodes[k]->adj.size() - 1; i++)
    {
        cvec[i] = new bool[T->G->get_num_nodes()];

        for (j = 0; j < T->G->get_num_nodes(); j++)
            cvec[i][j] = false;

        for (it_li = T->tree_nodes[children[i]]->bag.begin(); it_li != T->tree_nodes[children[i]]->bag.end(); ++it_li)
            cvec[i][*it_li] = true;
    }

    // preparing to split search space into pieces. 
    DEBUG("mask generation\n");

    parallel_mask_gen (T, k, w, T->start_masks, T->end_masks, T->entries);

	bool is_independent, has_bit;
	T->tree_nodes[k]->obj_val=0;
	TDSolution *hash_ptr=NULL;

    vector<bigint_t> indset_masks;
    vector<int> ind_wgh;
    unsigned long long ind_counter = 0;
    
    T->indset_masks = &indset_masks;
    T->ind_wgh = &ind_wgh;

 
	bigint_t current_mask(T->num_mask_words);
	bigint_t end_mask(T->num_mask_words);

	TDSolution temp_set(T->num_mask_words);

    double stime = MPI_Wtime();
    int request_flag = 0;
    int request_index = -1;
    int last_send_pos = 0;

    unsigned long long start_pos, end_pos, wchunk, counter, stwork;
    int compute_node_size = T->compute_nodes.size();
    int compute_node_rank = -1;

    compute_node_rank = std::find(T->compute_nodes.begin(), T->compute_nodes.end(), T->rank) - T->compute_nodes.begin();

    wchunk = T->entries/compute_node_size;
    start_pos = (wchunk * compute_node_rank);
    end_pos = (wchunk * (compute_node_rank + 1));

    //DEBUG("compute node size : %d rank : %d trank : %d \n", compute_node_size, compute_node_rank, T->rank);

    vector<int> requests(T->size, 0);
    vector<int> terminate(T->size, -1);

    int iteration = 0;

    for (i = 0; i < T->storage_nodes_size; i++)
    {
        (*T->request_count)[i] = 0;
    }


    DEBUG("DP loop\n");
    stime = MPI_Wtime ();
    counter = start_pos;
    while (1)
    {
        
        //DEBUG("counter %d end_pos %d\n", counter, end_pos);
        if (counter < end_pos)
        {
            current_mask = *(*T->start_masks)[counter];
            end_mask = *(*T->end_masks)[counter];
            counter ++;
        }

        while(current_mask < end_mask)
        {

            if ((T->rank == T->compute_nodes[0]) && (counter < end_pos))
            {

                if (MPI_SUCCESS != 
                    MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE))
                {
                    FERROR ("MPI Testany failed processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }

                while (request_flag)
                {
                    if (MPI_SUCCESS != 
                        MPI_Irecv (&(T->requester_rank), 1, MPI_INT, request_index, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, 
                                   &((*T->requests)[request_index])))
                    {
                        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize ();
                    }

                    requests[request_index] = 1;
                    
                    request_flag = 0;
                    request_index = -1;
                    T->requester_rank = -1;
                    if (MPI_SUCCESS != 
                        MPI_Testany (T->size, &T->requests->front(), &request_index, &request_flag, MPI_STATUS_IGNORE))
                    {
                        FERROR ("MPI Testany failed processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }
                }
                
                // Send work to the requester.
                for (i = 1; i < T->size && (end_pos > counter); i++)
                {
                    if (requests[i])
                    {
                        end_pos --;
                        if (MPI_SUCCESS != 
                            MPI_Send (&end_pos, 1, MPI_UNSIGNED_LONG_LONG, i, MPI_WRK_TAG, MPI_COMM_WORLD))
                        {
                            FERROR("MPI Send failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                            MPI_Finalize();
                        }
                        requests[i] = 0;
                    }
                }
            }
            
            is_independent=true;

            for(i = 0; i < w ; i++) 
            {            
                if(current_mask.test_bit(nbr_mask_vec[i]->k))
                {
                    // current_mask has the bit corresponding to
                    // nbr_mask_vec[i] set - so and with the appropriate
                    // vector for a hit 
                    has_bit=false;
                    // Do the loop below manually rather than calling the overloaded AND
                    for(j=0; j<current_mask.S; j++)
                    {
                        if(current_mask.words[j] & nbr_mask_vec[i]->w->words[j])
                        {
                            // Some nbr is found in word j
                            has_bit=true;
                            break;
                        }
                    }
                    // breaks to here
                    if( has_bit )
                    {                
                        // Set represented by current_mask is not independent
                        is_independent=false;
                        int rpos=0;
                        // Find the rightmost bit in current_mask - this
                        // could be sped up by looking at the words 
                        // SPEEDUP
                        while(!(current_mask.test_bit(rpos)))
                            rpos++;

                        // Now advance current_mask by a potentially massive amount
                        // Would be faster here to figure out which word we are in and 
                        // do this w/o test_bit
                        int b=0;
                        while(current_mask.test_bit(rpos+b))
                            b++;
                        for(int a=0;a<b;a++)
                            current_mask.xor_bit(rpos+a);
                        current_mask.or_bit(rpos+b);
                        break;                   
                    }
                }
            }


            if(is_independent)
            {   
                // ind_counter keep track of how many independent sets
                // are in the indset_mask vector
                ind_counter ++;
                T->indset_masks->push_back(current_mask);
                T->ind_wgh->push_back(0);
                ++current_mask;
                T->info->total_table_entries++;

                if (children_size > 0)
                {
                    if ((ind_counter - last_send_pos) >= T->request_size)
                    {
                        if (send_masks_9 (T, k, last_send_pos, iteration, &k_bag_vec, &cvec, &children) > 0)
                        {
                            last_send_pos = ind_counter;
                            iteration ++;
                        }
                        while(recv_wgh_9(T, k, &children));
                    }
                }
            }
        }


        // Once assigned work is finished, ask for more work
        if (!(counter < end_pos))
        {
            if (T->rank != T->compute_nodes[0])
            {
                MPI_Request maskrequest;
                MPI_Status maskstatus;
                int wflag = 0;
                int count = 0;
            
                //DEBUG("process : %d ask for work from : %d\n", T->rank, 0);
                if (MPI_SUCCESS !=
                    MPI_Send(&(T->rank), 1, MPI_INT, T->compute_nodes[0], MPI_REQ_WRK_TAG, MPI_COMM_WORLD))
                {
                    FERROR ("MPI Send failed at processor :%d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize ();
                }

                if (MPI_SUCCESS != 
                    MPI_Recv (&stwork, 1, MPI_UNSIGNED_LONG_LONG, T->compute_nodes[0], MPI_WRK_TAG, MPI_COMM_WORLD, &maskstatus))
                {
                    FERROR ("MPI Recv failed at processor :%d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize ();
                }

                count = 0;
                if (MPI_SUCCESS != 
                    MPI_Get_count (&maskstatus, MPI_UNSIGNED_LONG_LONG, &count))
                {
                    FERROR ("MPI Get count failed at processor :%d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize ();
                }

                if (!(count > 0))
                    break;

                current_mask.zeroize ();
                end_mask.zeroize ();
                current_mask = *(*T->start_masks)[stwork];
                end_mask = *(*T->end_masks)[stwork];
            }
            else
            {
                //Sending terminate signal
                int compute_node_size = T->compute_nodes.size();
                DEBUG("Send termination signal\n");
                for (i = 1; i < compute_node_size; i++)
                {
                    if (MPI_SUCCESS != 
                        MPI_Send (&end_pos, 0, MPI_UNSIGNED_LONG_LONG, T->compute_nodes[i], MPI_WRK_TAG, MPI_COMM_WORLD))
                    {
                        FERROR ("MPI Send failed in processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }
                }
                break;
            }
        }
    }

    DEBUG("iterations counter: %d\n", iteration);
    DEBUG("ind_counter: %d last_send_pos: %d\n", ind_counter, last_send_pos);
    if (children_size > 0)
    {
        if (ind_counter - last_send_pos > 0)
        {
            //DEBUG("send iterations outside loop: %d\n", iteration);
            if (send_masks_9(T, k, last_send_pos, iteration, &k_bag_vec, &cvec, &children) > 0)
            {
                last_send_pos = ind_counter;
                iteration ++;
            }
        }
        while(recv_wgh_9(T, k, &children));
    }

    DEBUG("Time for DP Loop : %f\n", MPI_Wtime() - stime);


    // Updating expected requests from each storage nodes.
    for (i = 0; i < T->storage_nodes_size; i++)
    {
        (*T->request_expect)[i] = 0;
    }

    for (i = 0; i < children_size; i++)
    {
        for (j = 0; j < T->storage_nodes_size; j++)
        {
            if (T->storage_nodes[j] == (*T->tbl_loc)[children[i]])
            {
                if (iteration > 0)
                    (*T->request_expect)[j] += iteration;
                else
                    (*T->request_expect)[j] += 1;
            }
        }
    }

    
    // Now do the DP over the children for this independent set - if this is a leaf, then
    // the body of the loop is never executed
    int ind_set_tbl_size = (int)T->indset_masks->size();
    MPI_Request wgh_s_request;
    MPI_Request wgh_r_request;
    int imaskcount = 0;
    int request_size = T->request_size;
    MPI_Status wstatus[T->size];
    request_index = 0;
    request_flag = 0;
    vector<int> *wgh_response;


    DEBUG("computing mask weights\n");
    DEBUG("mask count %d tble size %d\n", imaskcount, ind_set_tbl_size);
    if (!(ind_set_tbl_size > 0))
    {
       vector<unsigned long long> *msk_send_buffer;
       int destination = 0;
       for(i = 0; i < (int)children_size; i++)
       {
           destination = (*T->tbl_loc)[children[i]];
           msk_send_buffer =  T->msk_send_pool[destination];
           (*msk_send_buffer)[0] = children[i];
           (*msk_send_buffer)[1] = iteration;
           
           if (MPI_SUCCESS != 
               MPI_Send ((void *)&(*msk_send_buffer).front(), 2, MPI_UNSIGNED_LONG_LONG, 
                         destination, MPI_REQ_WGH_TAG, MPI_COMM_WORLD))
           {
               FERROR("MPI Send failed in processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
               MPI_Finalize();
           }
        }
    }

    //DEBUG("sending out requests\n");
    DEBUG("# of nbrs %d\n", children_size);
    DEBUG("child:location: ");
    for ( i = 0; i < children_size; i++)
    {
        GEN("%d:%d ", children[i], (*T->tbl_loc)[children[i]]);
    }
    GEN("\n");


    int complete_flag = 1;

    DEBUG("receiving weights from children\n");
    stime = MPI_Wtime();
    while (children_size > 0)
    {
        while(recv_wgh_9(T, k, &children));
        complete_flag = 1;
        for (i = 0; i < T->storage_nodes_size; i++)
        {
            //DEBUG("expect : %d received : %d\n", (*T->request_expect)[i], (*T->request_count)[i]);
            if ((*T->request_count)[i] < (*T->request_expect)[i])
                complete_flag = 0;
        }
        
        if (complete_flag)
            break;
    }
    DEBUG("Time for weight receive: %f\n", MPI_Wtime() - stime);
    DEBUG("finished recv weights\n");
    
    // Done with child DP
    TDSolution *temp_sol;
    int parent_weight = 0;

    // Updating weight for intersection with the parent.
    int maxvalue = 0;
    vector<int> weight = T->G->get_weight();
    for (i = 0; i < ind_set_tbl_size; i++)
    {
        temp_set.mask->zeroize();
        parent_weight=0;
        current_mask = (*T->indset_masks)[i];
    
        for(j = 0; j < w; j++)
        {                
            if(current_mask.test_bit(j))
            { 
                (*T->ind_wgh)[i] += weight[k_bag_vec[j]];
                if(parent_location_vec[k_bag_vec[j]]!=-1)
                {
                    temp_set.mask->set_bit(parent_location_vec[k_bag_vec[j]]);
                    parent_weight += weight[k_bag_vec[j]];
                }
            }
        }
        
        // temp_mask is a mask of this set \cap parent_bag
        hash_ptr=NULL;
        
        // k is not the root - so do the usual
        HASH_FIND(hh,T->tree_nodes[k]->hash_table, temp_set.mask->words,
                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

        if (k == T->root_node)
        {
            if ((*T->ind_wgh)[i] > maxvalue)
            {
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
                *(added_set->mask)=*(temp_set.mask);
                *(added_set->orig_mask)=current_mask;
                
                added_set->value = (*T->ind_wgh)[i];
                maxvalue = (*T->ind_wgh)[i];
                
                HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
        else
        {
            
            if(hash_ptr)
            {
                if ((*T->ind_wgh)[i] - parent_weight > hash_ptr->value)
                {
                    hash_ptr->value = (*T->ind_wgh)[i] - parent_weight;
                }
            }
            else
            {
                // There is no entry for this mask - add it to the table
                TDSolution *added_set;
                added_set=new TDSolution(T->num_mask_words);
                
            *(added_set->mask)=*(temp_set.mask);
            *(added_set->orig_mask)=current_mask;
            
            added_set->value = (*T->ind_wgh)[i] - parent_weight;
            
            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
            }
        }
    }
    

    if (MPI_SUCCESS != 
        MPI_Send ((void *)&k, 1, MPI_INT, 0, MPI_TBL_LOC_TAG, MPI_COMM_WORLD))
    {
        FERROR ("MPI Send failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    if (MPI_SUCCESS != 
        MPI_Wait (&(*T->tbl_loc_requests)[0], MPI_STATUS_IGNORE))
    {
        FERROR ("MPI Wait failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    int tbl_location = T->tblloc;
    (*T->tbl_loc)[k] = tbl_location;

    if (MPI_SUCCESS != 
        MPI_Irecv ((void *)&T->tblloc, 1, MPI_INT, 0, MPI_TBL_LOC_TAG, MPI_COMM_WORLD, &(*T->tbl_loc_requests)[0]))
    {
        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    DEBUG("Table for node : %d will be saved at : %d\n", k, tbl_location);
    // Before reaching this point we have calculate weights for the
    // masks in indset_masks vector.
	int table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);
	T->info->total_pc_table_entries+=(unsigned long long)table_size;
 
    DEBUG("size of the table : %d\n", table_size);

    DEBUG("Preparing to send the table\n");
    stime = MPI_Wtime ();
    hash_ptr = NULL;
    i = 0;
    int masksize = T->num_mask_words + 1;

    vector<unsigned long long> *maskweight_array;
    if (k != T->root_node)
        maskweight_array = new vector<unsigned long long>((table_size + 1)*(masksize), 0);

    maxvalue = 0;
    HASH_ITER(hh, T->tree_nodes[k]->hash_table, hash_ptr, temp_sol)
    {
        if (k != T->root_node)
        {
            for (j = 0; j < T->num_mask_words; j++)
            {
                (*maskweight_array)[i*(masksize) + j] = hash_ptr->mask->words[j];
            }
            (*maskweight_array)[i*(masksize) + j] = hash_ptr->value;
            i++;
        }
        else
        {
            if (hash_ptr->value > maxvalue)
                maxvalue = hash_ptr->value;
        }
    }


    if (k != T->root_node)
    {

        // Transmit 'k' value with masks and weights for storage
        (*maskweight_array)[masksize*(table_size) + T->num_mask_words] = k;
        // Deriving datatype
        MPI_Datatype mask_type;
        if (MPI_SUCCESS != 
            MPI_Type_contiguous (T->num_mask_words + 1, MPI_UNSIGNED_LONG_LONG, &mask_type))
        {
            FERROR ("MPI Type contiguous failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize ();
        }
        
        if (MPI_SUCCESS != 
            MPI_Type_commit(&mask_type))
        {
            FERROR ("MPI Type commit failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize ();
        }
        

        DEBUG("Send : %d elements to store at : %d\n", table_size + 1, tbl_location);
        if (MPI_SUCCESS != 
            MPI_Send ((void *)&maskweight_array->front(), table_size + 1, mask_type, tbl_location, MPI_TBL_SAVE_TAG, MPI_COMM_WORLD))
        {
            FERROR ("MPI Send failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
        DEBUG("Send successful !!\n");
    }
    else
    {
        int max = 0;

        DEBUG("max value : %d\n", maxvalue);
        if (MPI_SUCCESS != 
            MPI_Reduce ((void *)&maxvalue, (void *)&max, 1, MPI_INT, MPI_MAX, 0, *T->compute_nodes_comm))
        {
            FERROR ("MPI Reduce failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
        
        if (T->rank == T->compute_nodes[0])
        {
            T->info->opt_obj = max;
            GEN("maximum value is : %d\n", max);
        }
        DEBUG("max value : %d\n", maxvalue);
    }
    DEBUG("Time to send the table : %f\n", MPI_Wtime () - stime);

    // clear some vectors to gain more memory.
    if (k != T->root_node)
        delete maskweight_array;
    T->indset_masks->clear();
    T->ind_wgh->clear();

	for(i = 0; i < w;i++)
		delete nbr_mask_vec[i];
	
    // Get rid of the cvec's
	for(i=0;i< (int)T->tree_nodes[k]->adj.size()-1;i++)
		delete [] cvec[i];

	// Return the size of the table
	return table_size;
}

void parallel_wis_head_init(TDTree *T, int size, int rank)
{
    int i = 0;
    T->size = size;
    T->rank = rank;

    if (!(T->size > 0))
    {
        FERROR("Environment size can not be zero\n");
        MPI_Finalize();
        exit (-1);
    }

    if (T->rank < 0)
    {
        FERROR("Ranks can not be negative\n");
        MPI_Finalize();
        exit (-1);
    }

    try 
    {
        T->tbl_entries = new vector<int> (T->size,0);
        T->tbl_loc = new vector<int> (T->num_tree_nodes, -1);
        T->node_sizes = new vector<int> (T->size, 0);
        T->all_table_sizes = new vector<int> (T->num_tree_nodes, 0);

        T->tbl_size = new vector<int> (T->size);

        // Response with weights for masks comes with one additional
        // element, iteration #
        // Iteration # is used to find the correct position in original
        // vector for the weights received. 
        T->tbl_loc_requests = new vector<MPI_Request>(T->size);
        T->tbl_size_requests = new vector<MPI_Request>(T->size);
        T->compute_term_requests = new vector<MPI_Request>(T->size);

        for (i = 1; i < T->size; i++)
        {
            if (MPI_SUCCESS != 
                MPI_Irecv ((void *)&(*T->tbl_size)[i], 1, MPI_INT, i, MPI_TBL_SIZE_TAG, MPI_COMM_WORLD, &(*T->tbl_size_requests)[i]))
            {
                FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }

            if (MPI_SUCCESS != 
                MPI_Irecv ((void *)&T->tblloc, 1, MPI_INT, i, MPI_TBL_LOC_TAG, MPI_COMM_WORLD, &(*T->tbl_loc_requests)[i]))
            {
                FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }

            if (MPI_SUCCESS != 
                MPI_Irecv ((void *)&T->tmp, 1,
                           MPI_INT, i, MPI_COMP_TERM_TAG, MPI_COMM_WORLD, &(*T->compute_term_requests)[i]))
            {
                FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize ();
            }
        }
    }
    catch (std::bad_alloc e)
    {
        cerr << e.what() << endl;
        FERROR("Bad allocation: %s\n", e.what());
        MPI_Finalize();
        exit (-1);
    }

}

void parallel_wis_head_cleanup (TDTree *T)
{
    int i;
    for (i = 1; i < T->size; i++)
    {
        if (MPI_SUCCESS != 
            MPI_Cancel(&(*T->tbl_size_requests)[i]))
        {
            FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }

        if (MPI_SUCCESS != 
            MPI_Cancel(&(*T->tbl_loc_requests)[i]))
        {
            FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }

        if (MPI_SUCCESS != 
            MPI_Cancel(&(*T->compute_term_requests)[i]))
        {
            FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
    }
    delete T->tbl_entries;
    delete T->tbl_loc;
    delete T->tbl_size;

    delete T->tbl_size_requests;
    delete T->tbl_loc_requests;
    delete T->compute_term_requests;
}

void parallel_wis_compute_init (TDTree *T, int size, int rank)
{
    int i = 0;

    T->size = size;
    T->rank = rank;
    T->request_size = REQUEST_SIZE;
    T->msk_send_count = 0;
    T->wgh_send_count = 0;
    T->entries = 8192;

    if (!(T->size > 0))
    {
        FERROR("Environment size can not be zero\n");
        MPI_Finalize();
        exit (-1);
    }

    if (T->rank < 0)
    {
        FERROR("Ranks can not be negative\n");
        MPI_Finalize();
        exit (-1);
    }

    try 
    {
        T->start_masks = new vector<bigint_t *> (T->entries);
        T->end_masks = new vector<bigint_t *> (T->entries);
        for (i = 0; i < T->entries; i++)
        {
            (*T->start_masks)[i] = new bigint_t (T->num_mask_words);
            (*T->end_masks)[i] = new bigint_t (T->num_mask_words);
        }

        T->tbl_loc = new vector<int> (T->num_tree_nodes, -1);
        T->request_expect = new vector<int> (T->storage_nodes_size, 0);
        T->request_count = new vector<int> (T->storage_nodes_size, 0);

        if (T->rank == T->compute_nodes[0])
        {
            // request_count vector keep track of # of requests each node has
            // to receive.
               T->requests = new vector<MPI_Request>(T->size);

            for (i = 0; i < T->size; i++)
            {
                if (MPI_SUCCESS != 
                    MPI_Irecv (&(T->requester_rank), 1, MPI_INT, i, MPI_REQ_WRK_TAG, MPI_COMM_WORLD, &((*T->requests)[i])))
                {
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }
            }
        }

        // Response with weights for masks comes with one additional
        // element, iteration #
        // Iteration # is used to find the correct position in original
        // vector for the weights received. 
        T->wgh_recv_requests = new vector<MPI_Request>(T->size);
        for (i = 0; i < T->size; i++)
        {
            T->wgh_recv_pool.push_back(new vector<int>(T->request_size + 2, 0));
            if (MPI_SUCCESS != 
                MPI_Irecv ((void *)&T->wgh_recv_pool[i]->front(), T->request_size + 2, MPI_INT, i, 
                           MPI_WGH_TAG, MPI_COMM_WORLD, &(*T->wgh_recv_requests)[i]))
            {
                FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }
        }

        for (i = 0; i < T->size; i++)
            T->msk_send_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));

        T->tbl_loc_requests = new vector<MPI_Request>(1);
        if (MPI_SUCCESS != 
            MPI_Irecv ((void *)&T->tblloc, 1, MPI_INT, 0, MPI_TBL_LOC_TAG, MPI_COMM_WORLD, &(*T->tbl_loc_requests)[0]))
        {
            FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }

    }
    catch (std::bad_alloc e)
    {
        cerr << e.what() << endl;
        FERROR("Bad allocation: %s\n", e.what());
        MPI_Finalize();
        exit (-1);
    }
}

void parallel_wis_compute_cleanup (TDTree *T)
{
    int i;
    for (i = 0; i < T->size; i++)
    {
        if (T->rank == T->compute_nodes[0])
        {
            if (MPI_SUCCESS != 
                MPI_Cancel (&(*T->requests)[i]))
            {
                
                FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }
        }

        if (MPI_SUCCESS != 
            MPI_Cancel(&(*T->wgh_recv_requests)[i]))
        {
            FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
    }

    if (MPI_SUCCESS != 
        MPI_Cancel (&(*T->tbl_loc_requests)[0]))
    {
        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    for (i = 0; i < T->size; i++)
    {
        delete  T->msk_send_pool[i];
        delete  T->wgh_recv_pool[i];
    }

    T->start_masks->clear();
    T->end_masks->clear();
    delete T->start_masks;
    delete T->end_masks;
}


void parallel_wis_storage_init (TDTree *T, int size, int rank)
{
    int i = 0;

    T->size = size;
    T->rank = rank;

    if (!(T->size > 0))
    {
        FERROR("Environment size can not be zero\n");
        MPI_Finalize();
        exit (-1);
    }

    if (T->rank < 0)
    {
        FERROR("Ranks can not be negative\n");
        MPI_Finalize();
        exit (-1);
    }

    T->request_size = REQUEST_SIZE;
    T->msk_send_count = 0;
    T->wgh_send_count = 0;
    T->entries = 8192;

    T->tbl_entries = new vector<int> (T->num_tree_nodes, 0);

    T->msk_recv_requests = new vector<MPI_Request>(T->size);
    for (i = 0; i < T->size; i++)
    {
        T->msk_recv_pool.push_back(new vector<unsigned long long>(T->request_size*T->num_mask_words + 2, 0));
        if (MPI_SUCCESS != 
            MPI_Irecv((void *)&(*T->msk_recv_pool[i]).front(), T->request_size*T->num_mask_words + 2, 
                      MPI_UNSIGNED_LONG_LONG, i, MPI_REQ_WGH_TAG, MPI_COMM_WORLD, &(*T->msk_recv_requests)[i]))
        {
            FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
    }


    // T->wgh_send_requests = new vector<MPI_Request>(T->size);
    for (i = 0; i < T->size; i++)
        T->wgh_send_pool.push_back(new vector<int>(T->request_size + 2, 0));


    T->storage_term_requests = new vector<MPI_Request>(1);
    if (MPI_SUCCESS != 
        MPI_Irecv((void *)&T->tmp, 1, MPI_INT, 0, MPI_STOR_TERM_TAG, MPI_COMM_WORLD, &(*T->storage_term_requests)[0]))
    {
        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    T->child_free_requests = new vector<MPI_Request>(1);
    if (MPI_SUCCESS != 
        MPI_Irecv((void *)&T->del_ch, 1, MPI_INT, MPI_ANY_SOURCE, MPI_CHILD_FREE_TAG, MPI_COMM_WORLD, &(*T->child_free_requests)[0]))
    {
        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

}

void parallel_wis_storage_cleanup (TDTree *T)
{
    int i;
    for (i = 0; i < T->size; i++)
    {
        if (MPI_SUCCESS != 
            MPI_Cancel (&(*T->msk_recv_requests)[i]))
        {
            FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
            MPI_Finalize();
        }
    }

    if (MPI_SUCCESS != 
        MPI_Cancel(&(*T->storage_term_requests)[0]))
    {
        FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    if (MPI_SUCCESS != 
        MPI_Cancel(&(*T->child_free_requests)[0]))
    {
        FERROR ("MPI Cancel failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize();
    }

    for (i = 0; i < T->size; i++)
    {
        delete  T->msk_recv_pool[i];
        delete  T->wgh_send_pool[i];
    }
}


void parallel_wis_head_start (TDTree *T)
{

    int compute_node_index = -1;
    int compute_node_flag = 0;
    MPI_Status compute_node_status;

    int storage_node_index = -1;
    int storage_node_flag = 0;
    MPI_Status storage_node_status;

    int storage_nodes_size;
    MPI_Request storage_term_request;

    int i;
    int node;
    int location;
    int source;
    int min = INT_MAX;

    int entries;
    int ch_last;
    
    int count_tbl_loc = 0;

    DEBUG("Starting head node\n");
    storage_nodes_size = T->storage_nodes.size();
    while (1)
    {

        //Start: Looking for table size data from storage nodes
        {
            if (MPI_SUCCESS != 
                MPI_Testany (T->size, &T->tbl_size_requests->front(), 
                             &storage_node_index, &storage_node_flag, &storage_node_status))
            {
                FERROR ("MPI Test any failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }
            
            while (storage_node_flag)
            {
                entries = (*T->tbl_size)[storage_node_index];
                //(*T->node_sizes)[storage_node_index] = entries;
                //DEBUG("node : %d entries : %d \n", storage_node_index, entries);
                
                if (MPI_SUCCESS != 
                    MPI_Irecv ((void *)&(*T->tbl_size)[storage_node_index], 1, MPI_INT, 
                               storage_node_index, MPI_TBL_SIZE_TAG, MPI_COMM_WORLD, &(*T->tbl_size_requests)[storage_node_index]))
                {
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }


                storage_node_index = -1;
                storage_node_flag = 0;

                if (MPI_SUCCESS != 
                    MPI_Testany (T->size, &T->tbl_size_requests->front(), 
                                 &storage_node_index, &storage_node_flag, &storage_node_status))
                {
                    FERROR ("MPI Test any failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }

            }
        }

        // Start: Looking for table location requests from compute nodes
        {
            if (MPI_SUCCESS != 
                MPI_Testany (T->size, &T->tbl_loc_requests->front(), &compute_node_index, 
                             &compute_node_flag, &compute_node_status))
            {
                FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }

            while (compute_node_flag)
            {
                node = T->tblloc;
                source = compute_node_index;

                if (MPI_SUCCESS != 
                    MPI_Irecv ((void *)&T->tblloc, 1, MPI_INT, compute_node_index, MPI_TBL_LOC_TAG, 
                               MPI_COMM_WORLD, &(*T->tbl_loc_requests)[compute_node_index]))
                {
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }
                //DEBUG("Got request for tblloc for : %d\n",  node);

                if ((*T->tbl_loc)[node] < 0)
                {
                    min = INT_MAX;
                    location = -1;
                    for (i = 0; i < storage_nodes_size; i++)
                    {
                        // DEBUG("min rank:%d value: %d\n", T->storage_nodes[i], (*T->node_sizes)[T->storage_nodes[i]]);
                        if ((*T->node_sizes)[T->storage_nodes[i]] < min)
                        {
                            location = T->storage_nodes[i];
                            min = (*T->node_sizes)[T->storage_nodes[i]];
                        }
                    }
                    (*T->tbl_loc)[node] = location;
                    (*T->node_sizes)[location] ++;

                    DEBUG("Save node : %d at : %d\n", node, location);
                }
                else
                {
                    location = (*T->tbl_loc)[node];
                }
            

                if (MPI_SUCCESS != 
                    MPI_Send ((void *)&location, 1, MPI_INT, source, MPI_TBL_LOC_TAG, MPI_COMM_WORLD))
                {
              
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }

                compute_node_flag = 0;
                compute_node_index = -1;

                if (MPI_SUCCESS != 
                    MPI_Testany (T->size, &T->tbl_loc_requests->front(), &compute_node_index, 
                                 &compute_node_flag, &compute_node_status))
                {
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }
            }
        }
        // End : Looking for table location requests from compute nodes
        
        
        // Start: Looking for terminate signal from compute nodes
        {            
            if (MPI_SUCCESS != 
                MPI_Testany (T->size, &T->compute_term_requests->front(), &compute_node_index, 
                             &compute_node_flag, &compute_node_status))
            {
                FERROR ("MPI Test any failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize ();
            }
    
            if (compute_node_flag)
            {
                DEBUG("Received terminate signal from compute node : %d\n", compute_node_index);
                if (MPI_SUCCESS != 
                    MPI_Irecv ((void *)&T->tmp, 1,
                               MPI_INT, compute_node_index, MPI_COMP_TERM_TAG, MPI_COMM_WORLD, &(*T->compute_term_requests)[compute_node_index]))
                {
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                }
            
                storage_nodes_size = T->storage_nodes.size ();
                for (i = 0; i < storage_nodes_size; i++)
                {
                    DEBUG("Sending terminate signal to storage node : %d\n", T->storage_nodes[i]);
                    if (MPI_SUCCESS != 
                        MPI_Send ((void *)&T->tmp, 1, MPI_INT, T->storage_nodes[i], MPI_STOR_TERM_TAG, MPI_COMM_WORLD))
                    {
                    
                        FERROR ("MPI Send failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }
                }

                break;
            }
        }
        // End : Looking for terminate signal from compute nodes
    }
}


void parallel_wis_storage_start (TDTree *T)
{

    int head_node_flag = 0;
    MPI_Status head_node_status;

    int storage_node_flag;
    MPI_Status storage_node_status;

    int compute_node_flag = 0;
    int compute_node_index = 0;
    int compupte_node_status;

    MPI_Status tbl_recv_status;

    int element_count, count;
    vector<unsigned long long> *table_buffer;

    int child_free_flag;
    MPI_Status child_free_status;

    int i, j, k;
    int k_pos;
    int table_size;
    int entries = 0;
    TDSolution *hash_ptr;

    // Deriving datatype
    int masksize = T->num_mask_words + 1;
    MPI_Datatype mask_type;
    if (MPI_SUCCESS != 
        MPI_Type_contiguous (masksize, MPI_UNSIGNED_LONG_LONG, &mask_type))
    {
        FERROR ("MPI Type contiguous failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
    }

    if (MPI_SUCCESS != 
        MPI_Type_commit(&mask_type))
    {
        FERROR ("MPI Type commit failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
        MPI_Finalize ();
    }

    DEBUG("Starting storage node\n");
    while (1)
    {
        // Start: Reply to weight requests
        {
            while(recv_msk_send_wgh_9(T));
        }
        // End: Reply to weight requests

        // Start: Receive table from compute nodes.
        {
            if (MPI_SUCCESS != 
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_TBL_SAVE_TAG, MPI_COMM_WORLD, &storage_node_flag, &storage_node_status))
            {
                FERROR ("MPI Iprobe failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }
        
            while (storage_node_flag)
            {

                if (MPI_SUCCESS != 
                    MPI_Get_count (&storage_node_status, mask_type, &element_count))
                {
                
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }
            
                if (element_count > 0)
                {

                    table_buffer = new vector<unsigned long long> (masksize * element_count);
                    if (!table_buffer)
                    {
                        FERROR("Can not allocate memory for receiving table\n");
                        MPI_Finalize ();
                    }
                
                    if (MPI_SUCCESS != 
                        MPI_Recv ((void *)&table_buffer->front(), element_count, mask_type, storage_node_status.MPI_SOURCE, 
                                  storage_node_status.MPI_TAG, MPI_COMM_WORLD, &tbl_recv_status))
                    {
                        FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }

                    k = (*table_buffer)[masksize*(element_count - 1) + T->num_mask_words];

                    for (i = 0; i < (element_count - 1); i++)
                    {

                        TDSolution *added_set;
                        added_set=new TDSolution(T->num_mask_words);
                        if (!added_set)
                        {
                            FERROR("no enough memory \n");
                            MPI_Finalize();
                            exit(-1);
                        }

                        for (j = 0; j < T->num_mask_words; j++)
                        {
                            added_set->mask->words[j] = (*table_buffer)[i*masksize + j];
                        }
                        added_set->value= (*table_buffer)[i*masksize + j];

                        HASH_FIND(hh,T->tree_nodes[k]->hash_table, added_set->mask->words,
                                  T->num_mask_words*sizeof(BIGINT_WORD), hash_ptr);

                        if (!hash_ptr)
                        {
                            HASH_ADD_KEYPTR(hh, T->tree_nodes[k]->hash_table, added_set->mask->words,
                                            T->num_mask_words*sizeof(BIGINT_WORD), added_set);
                        }
                        else
                        {
                            if (added_set->value > hash_ptr->value)
                            {
                                hash_ptr->value = added_set->value;
                            }
                        }
                    }
                    
                    delete table_buffer;

                    table_size = HASH_COUNT(T->tree_nodes[k]->hash_table);

                    if ((*T->tbl_entries)[k] < table_size)
                    {
                        entries -= (*T->tbl_entries)[k];
                        entries += table_size;
                        (*T->tbl_entries)[k] = table_size;
                    }

                    if (MPI_SUCCESS != 
                        MPI_Send ((void *)&entries, 1, MPI_INT, 0, MPI_TBL_SIZE_TAG, MPI_COMM_WORLD))
                    {
                        FERROR ("MPI Send failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                        MPI_Finalize();
                    }
                    DEBUG("Received table for node : %d from : %d\n", k, storage_node_status.MPI_SOURCE);
                }

                storage_node_flag = 0;
                if (MPI_SUCCESS != 
                    MPI_Iprobe(MPI_ANY_SOURCE, MPI_TBL_SAVE_TAG, MPI_COMM_WORLD, &storage_node_flag, &storage_node_status))
                {
                    FERROR ("MPI Iprobe failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }
            }
        }
        // End: Receive table from compute nodes


        // Start: Free children
        {
            if (MPI_SUCCESS != 
                MPI_Test (&(*T->child_free_requests)[0], &child_free_flag, &child_free_status))
            {
                FERROR ("MPI Test failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }

            if (child_free_flag)
            {
                int k;
                list<int>::iterator ii;

                child_free_flag = 0;
                DEBUG("Delete child table for : %d\n", T->del_ch);

                k = T->del_ch;
                if (MPI_SUCCESS != 
                    MPI_Irecv((void *)&T->del_ch, 1, MPI_INT, MPI_ANY_SOURCE, MPI_CHILD_FREE_TAG, MPI_COMM_WORLD, &(*T->child_free_requests)[0]))
                {
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }

                ii= T->tree_nodes[k]->adj.begin(); 
                ++ii;
                for(  ; ii!=T->tree_nodes[k]->adj.end(); ++ii)
                {
                    //DEBUG("Deleting : %d\n", (int)*ii);
                    if(T->tree_nodes[*ii]->hash_table)
                    {
                        TDSolution *current_set, *tmp;
                        HASH_ITER(hh, T->tree_nodes[*ii]->hash_table, current_set, tmp)
                        {
                            HASH_DEL(T->tree_nodes[*ii]->hash_table,current_set);
                            //print_message(0,"About to delete TDSolution\n");
                            //current_set->mask->print(stderr);
                            delete current_set;
                        }
                        T->tree_nodes[*ii]->hash_table=NULL;
                    }
                }
            }
        }
        // End: Free children


        // Start: Check for terminate signal from head node
        {
            if (MPI_SUCCESS != 
                MPI_Test (&(*T->storage_term_requests)[0], &head_node_flag, &head_node_status))
            {
                FERROR ("MPI Test failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                MPI_Finalize();
            }

            if (head_node_flag)
            {
                DEBUG ("Terminate storage node\n");
                if (MPI_SUCCESS != 
                    MPI_Irecv((void *)&T->tmp, 1, MPI_INT, 0, MPI_STOR_TERM_TAG, MPI_COMM_WORLD, &(*T->storage_term_requests)[0]))
                {
                    FERROR ("MPI Irecv failed at processor : %d %s:%d\n", T->rank, __FILE__, __LINE__);
                    MPI_Finalize();
                }
                break;
            }
        }
        // End: Check for terminate signal from head node
    }
}

#endif
