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

#ifndef NODE_H_
#define NODE_H_
#include <list>
#include <vector>

using namespace std;

namespace Graph {
    class Node
    {
public:

    Node();
    virtual ~Node();
    int get_label() const;
    list<int> get_nbrs() const;
    list<int> *get_nbrs_ptr();
    int get_degree() const;
    void set_label(int label);
    void set_nbr(list<int> nbr);
    void add_nbr(int i);
    void remove_nbr(int u);
    void delete_node();
    void sort_nbr();
    void reverse_nbr();
    void unique_nbr();

    /**
     * \brief get the index of the neighbor whose index is the largest below n
     */
    std::list<int>::const_reverse_iterator get_largest_neighbor_below(int n);
    Node & operator =(const Node & n);

    //Declare friend classes
    friend class Graph;
    friend class GraphUtil;
    friend class GraphProperties;
    friend class GraphDisplay;
    friend class GraphEOUtil;
    friend class GraphCreator;

private:
    list<int> nbrs;
    int label;
    };
}

#endif /* NODE_H_ */
