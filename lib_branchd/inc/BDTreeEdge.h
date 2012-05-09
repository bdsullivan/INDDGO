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

#ifndef _BD_TREE_EDGE_H_
#define _BD_TREE_EDGE_H_

class BDTreeEdge
{
     friend ostream &operator<<(ostream &, const BDTreeEdge &);

public:
    BDTreeEdge(){};
    ~BDTreeEdge();

    int start, end;         // the indices in the BDTreeNode array of the BDTreeNodes joined by this edge 
    list<int> middle_set;   // list of vertices from the Graph in the middle set

};

#endif

