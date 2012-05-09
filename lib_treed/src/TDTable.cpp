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

TDSubset::TDSubset()
{
    this->processed=false;
}

TDTable::TDTable()
{
}

TDTable::~TDTable()
{
    this->sets.clear();
}

/**
* Overloading the << operator for the TDTable class.
*/
ostream &operator<<(ostream &output, const TDTable &T)
{
    list<int>::const_iterator i;
    output<<"Table contains "<<T.sets.size()<<" subsets"<<endl;
    
    list<TDSubset>::const_iterator ss;
    int j=1;
    for(ss=T.sets.begin();ss!=T.sets.end();++ss)
    {
        output<<"\nSubset "<<j<<"(cost is "<<(*ss).value<<", mask is "<<(unsigned long)(*ss).mask<<")\n";
        j++;
    }
    
    return output;
}

bool operator <(const TDSubset &a, const TDSubset &b)
{
    return a.mask < b.mask;
    
}

bool forget_compare(TDSubset a, TDSubset b)
{
    return ((a.mask & a.aux_mask) < (b.mask & b.aux_mask));

}

bool aux_mask_compare(TDSubset a, TDSubset b)
{
    return a.aux_mask < b.aux_mask;

}
