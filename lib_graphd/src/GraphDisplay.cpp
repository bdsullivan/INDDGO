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

#include "GraphDisplay.h"
#include <stdio.h>

namespace Graph {
    GraphDisplay::GraphDisplay(){
        // TODO Auto-generated constructor stub
    }

    GraphDisplay::~GraphDisplay(){
        // TODO Auto-generated destructor stub
    }

    /**
     * Debugging function that prints the current adjacency lists to stdout.
     */
    void GraphDisplay::show_adjlists(Graph *g){
        std::list<int>::iterator it;
        int i;

        for(i = 0; i < g->Graph::capacity; i++){
            printf("%d: ", i);
            if(g->nodes[i].label == -1){
                printf("VACANT");
            }
            else {
                it = g->nodes[i].nbrs.begin();
                while(it != g->nodes[i].nbrs.end()){
                    printf("%d ", *it);
                    it++;
                }
            }
            printf("\n");
        }
    } // show_adjlists
}
