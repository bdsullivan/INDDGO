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

#ifndef _BD_TREE_NODE_H_
#define _BD_TREE_NODE_H_

// This is here so that the tree node can reference the tree that contains it!
// The real BDTree is defined later - will this work???
class BDTree;

class BDTreeNode
{
     friend ostream &operator<<(ostream &, const BDTreeNode &);
public:
    BDTreeNode();       // Constructor - will we need parameters?
    ~BDTreeNode();      // Destructor

    int id;                                 // Unique id for this node - corresponds to index in the BDTree
    int graph_edge_start, graph_edge_end;   // If the BDTreeNode is a leaf, then this is
                                            // the information aboupt the edge in the original graph 
                                            // that this BDTreeNode represents.
   
    // List of indices of edges in the BDTreeEdge array incident to this BDTreeNode
    // Once the tree is valid (ternary) and rooted, then the edge on the path to the root
    // node is first in this list!! The other two edges (left/right) are listed arbitrarily
    list<int> edges;   
    
    int num_leaf_nbrs;

    BDTree *tree; // pointer to the tree containing the node
};

#endif
