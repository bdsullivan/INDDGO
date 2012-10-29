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

TDSolution::TDSolution(int num_words)
{
	this->mask = new bigint_t(num_words);
	this->orig_mask = new bigint_t(num_words);
	this->value = 0;
	this->mask_length = num_words;
}

TDSolution::TDSolution()
{
	this->mask = NULL;
	this->orig_mask = NULL;
	this->value = 0;
	this->mask_length = 0;

}

TDSolution::~TDSolution()
{
	if (this->mask)
		delete this->mask;
	this->mask = NULL;

	if (this->orig_mask)
		delete this->orig_mask;
	this->orig_mask = NULL;

	this->value = 0;

}

// Copy constructor
TDSolution::TDSolution(const TDSolution& rhs)
{
	this->mask = NULL;
	if (rhs.mask)
	{
		this->mask = new bigint_t(rhs.mask->S);
		*this->mask = *rhs.mask;
	}

	this->orig_mask = NULL;
	if (rhs.orig_mask)
	{
		this->orig_mask = new bigint_t(rhs.orig_mask->S);
		*this->mask = *rhs.orig_mask;
	}

	this->value = 0;
	if (rhs.value)
		this->value = rhs.value;

	this->mask_length = 0;
	if (rhs.mask_length)
		this->mask_length = rhs.mask_length;
}

TDSolution& TDSolution::operator=(const TDSolution& rhs)
{
	this->mask = NULL;
	if (rhs.mask)
	{
		this->mask = new bigint_t(rhs.mask->S);
		*this->mask = *rhs.mask;
	}

	this->orig_mask = NULL;
	if (rhs.orig_mask)
	{
		this->orig_mask = new bigint_t(rhs.orig_mask->S);
		*this->mask = *rhs.orig_mask;
	}

	this->value = 0;
	if (rhs.value)
		this->value = rhs.value;

	this->mask_length = 0;
	if (rhs.mask_length)
		this->mask_length = rhs.mask_length;

	return *this;
}


// set mask length
void TDSolution::set_mask_length(int len)
{
	this->mask_length = len;
}

// get mask length
int TDSolution::get_mask_length() const
{
	return this->mask_length;
}

