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

#ifndef _TD_TABLE_H_
#define _TD_TABLE_H_

class TDSubset
{
public:
    TDSubset();
    // The < operator just compares based on mask
    friend bool operator<(const TDSubset& a, const TDSubset &b);
    // The TDSubset generally lives inside a TDTable 
    // (which lives inside a TDTreeNode).   The mask represents a subset
    // of vertices from the bag belonging to this TDTreeNode

    BIGINT mask;

    // The aux_mask is used for convenience when doing various 
    // tasks involved in the dynamic programming

    BIGINT aux_mask;
    BIGINT orig_aux_mask;
    
    // This is the value/weight attached to the subset
    double value;
    double temp_val;
    // Flag to indicate whether or not the subset is processed or not
    bool processed;
};


class TDTable
{
    friend ostream &operator<<(ostream &, const TDTable &);

public:
    TDTable();
    ~TDTable();
    list<TDSubset> sets;
};

bool forget_compare(TDSubset a, TDSubset b);
bool aux_mask_compare(TDSubset a, TDSubset b);

#endif

