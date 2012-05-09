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

#ifndef _BRANCH_DEC_H
#define _BRANCH_DEC_H


#if WIN32 || _WIN32
 // Define these here if using Windows
 // Otherwise, define them in the makefile
 #define HAS_METIS   1
 #define HAS_AMD     1
 #define HAS_GOBLIN  1
#endif

#if HAS_GOBLIN
#include <goblin.h>
#endif
#include "BDTreeEdge.h"
#include "BDTreeNode.h"
#include "BDTree.h"
//#include "ArpackUtils.h" 

#endif


