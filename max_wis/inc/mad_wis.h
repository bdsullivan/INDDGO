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

#ifndef __MAD_WIS_H__
#define __MAD_WIS_H__
#ifdef __MADNESS__
#include <map>
using namespace std;

//void async_table_compute(TDTree *T, int node);
void async_table_compute(map<int, TDTreeNode*>& m, TDTreeNode& node);
void madness_async_table_compute (TDTree *T, int node);
int compute_nonnice_table_map(map<int, TDTreeNode*>& m, TDTreeNode *node);

#endif /* __MADNESS__ */
#endif /* __MAD_WIS_H__ */
