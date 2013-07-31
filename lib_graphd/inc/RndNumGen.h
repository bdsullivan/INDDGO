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

#ifndef RNDNUMGEN_H_
#define RNDNUMGEN_H_

namespace Graph {
    void init_lcgrand(int stream, int seed);
    double lcgrand(int stream);
    void random_permutation(int *perm, int n);
    int rand_int(int min, int max);
}

#endif /* RNDNUMGEN_H_ */
