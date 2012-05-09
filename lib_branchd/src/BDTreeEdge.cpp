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
#include "BranchDecomposition.h"

BDTreeEdge::~BDTreeEdge()
{
  this->middle_set.clear();

}

ostream &operator<<(ostream &output, const BDTreeEdge &T)
{
    list<int>::const_iterator ii;

    output<<"("<<T.start<<","<<T.end<<"):";
    output<<"\t{";
    ii=T.middle_set.begin();
    if(ii!=T.middle_set.end())
        output<<*ii;
    ++ii;
    while(ii!=T.middle_set.end())
    {
        output<<","<<*ii;
        ++ii;
    }
    output<<"}\n";

    return output;
}


