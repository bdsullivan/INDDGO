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

#include "Node.h"
#include "GraphException.h"
#include <algorithm>
namespace Graph {
    Node::Node(){
        this->label = -1;
    }

    int Node::get_label() const {
        return label;
    }

    list<int> Node::get_nbrs() const {
        return nbrs;
    }

    list<int> *Node::get_nbrs_ptr(){
        return &(this->nbrs);
    }

    void Node::set_label(int label){
        this->label = label;
    }

    void Node::set_nbr(list<int> nbr){
        this->nbrs = nbr;
    }

    void Node::add_nbr(int i){
        nbrs.push_back(i);
    }

    Node::~Node(){
    }

    void Node::remove_nbr(int u){
        list<int>::iterator iter;
        for(iter = nbrs.begin(); iter != nbrs.end(); ++iter){
            if(*iter == u){
                nbrs.erase(iter);
                break;
            }
        }
    }

    int Node::get_degree() const {
        return nbrs.size();
    }

    void Node::delete_node(){
        nbrs.clear();
    }

    void Node::sort_nbr(){
        nbrs.sort();
    }

    void Node::unique_nbr(){
        nbrs.unique();
    }

    /**
     * Used to save some compares when calculating triangles. **ASSUMES NEIGHBOR LIST IS SORTED**
     * \param[in] n the index below which to return the largest neighbor
     * \return x the index of the largest neighbor below n
     */
    int Node::get_largest_neighbor_below(int n){
    //FIXME: probably should be able to handle unsorted case, but how to do without sacrificing efficiency?
    //       also maybe do a binary instead of linear search
        int v;
        std::list<int>::reverse_iterator it;

        for(it = this->nbrs.rbegin(); it != this->nbrs.rend(); ++it){
            if(*it < n) {
                return *it;
            }
        }
        return -1;
    }

    Node & Node::operator=(const Node & n){
        nbrs = n.nbrs;
        label = n.label;

        return *this;
    }
}
