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

#ifndef _TD_SUBSET_H_
#define _TD_SUBSET_H_

#include <uthash.h>
#include <bigint.h>
using namespace std;

class TDSolution
{
private:
	/**
	* length of the masks stored in the hash table
	*/
	int mask_length;

public:
	TDSolution(int num_words);
	TDSolution();
	~TDSolution();

	// Copy constructor
	TDSolution(const TDSolution& rhs);
	// Assignment operator
	TDSolution& operator=(const TDSolution& rhs);
	// The < operator just compares based on mask
	friend bool operator<(const TDSolution& a, const TDSolution &b);
	// The mask field is used to represent the intersection of
	// a set with its parent's bag, since this is all the parent needs to know
	// when he does his DP operation.  When stored in the hash table,
	// the mask is adjusted so that it represents vertices in the parent bag
	bigint_t *mask;

	// The orig_mask field is used to represent the original mask
	// of this set..  So it must be the case that the set represented
	// by orig_mask contains the set represented by mask.  However, note that
	// orig_mask and mask cannot be directly compared position-by-position
	// since the orig_mask bits reference nodes in this tree node's bag
	// while the mask bits reference nodes in the parent's bag.
	bigint_t *orig_mask;

	// This is the value/weight attached to the subset
	// CSG changing 6/10/2011
	int value;

	// set mask length
	void set_mask_length(int len);

	// get mask length
	int get_mask_length() const;

	// handle for UThash
	UT_hash_handle hh;
};

//
//bool forget_compare(TDSolution a, TDSolution b);
//bool aux_mask_compare(TDSolution a, TDSolution b);
//bool aux_mask_ptr_compare(TDSolution *a, TDSolution *b);

#endif

