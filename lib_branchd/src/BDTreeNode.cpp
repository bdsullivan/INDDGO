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

BDTreeNode::BDTreeNode()
{
    this->num_leaf_nbrs=0;
}
/**
* Destructor for the BDTree.
*/
BDTreeNode::~BDTreeNode()
{
  this->edges.clear();
}


/*
* Overloading the << operator for BDTreeNode.
*/
ostream &operator<<(ostream &output, const BDTreeNode &T)
{
    list<int>::const_iterator ii;

    output<<"TreeNode "<<T.id;
    if(T.graph_edge_start!=-1)
        output<<" (graph edge is "<<T.graph_edge_start<<"-"<<T.graph_edge_end<<")\n";
    else
        output<<" (interior)\n";
    //output<<"  Nbrs: ";
    //for(list<int>::const_iterator ii=T.adj.begin();ii!=T.adj.end();++ii)
    //    output<<*ii<<" ";
    //output<<endl;
    output<<"Incident tree edges:\n";
    for(ii=T.edges.begin(); ii!=T.edges.end(); ++ii)
    {
        output<<"\tEdge "<<*ii<<":"<<T.tree->edges[*ii];
    }
    output<<endl;
    return output;
}




