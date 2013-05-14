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

#include "Debug.h"
#include <stdio.h>
#include <stdlib.h>
#include "GraphException.h"
#include "Log.h"
/**
 * Print the message to stderr if level <= DEBUG_LEVEL.
 */
void print_message(int level, const char *format, ...){
    #ifdef _NDEBUG
    // Don't do anything if _NDEBUG is defined since
    // we have compiled without DEBUGging enabled
    return;
    #else

    if(level > DEBUG_LEVEL){
        return;
    }

    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    #endif /* NDEBUG */
} // print_message

/**
 * Prints a fatal error message to stderr and exits.
 */
void fatal_error(const char *format, ...){
    fprintf(stderr,"Fatal Error!\n");
    char buffer[1024];
    va_list args;

    va_start(args, format);
    vfprintf(stderr, format, args);
    sprintf(buffer, format, args);
    va_end(args);

    fprintf(stderr,"Exiting\n");

    FERROR("%s", buffer);
    const std::string desc(buffer);
    throw Graph::GraphException(desc);
}

