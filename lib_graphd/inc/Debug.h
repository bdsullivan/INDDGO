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

#ifndef DEBUG_H_
#define DEBUG_H_
#include <iostream>
#include <stdarg.h>

// Taken from http://oopweb.com/CPP/Documents/DebugCPP/Volume/techniques.html

/**
 * The higher the number, the more output you see.
 * Setting DEBUG_LEVEL=0 will print no debugging statements.
 */
#define DEBUG_LEVEL 0

void print_message(int level, const char *format, ...);
void fatal_error(const char *format, ...);
template<typename Container> void print(int level, const Container& cntn){
    #ifdef _NDEBUG
    // Don't do anything if _NDEBUG is defined since
    // we have compiled without DEBUGging enabled
    return;
    #else
    if(level > DEBUG_LEVEL){
        return;
    }

    typename Container::const_iterator i = cntn.begin();
    for(; i != cntn.end(); ++i){ //iterate through the container
        std::cerr << *i << " ";         //print each element
    }
    std::cerr << std::endl;
    #endif
} // print

;

#endif /* DEBUG_H_ */
