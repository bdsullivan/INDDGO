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

#include <stdlib.h>
#include <getopt.h>
#include <sstream>

#include "Graph.h"
#include "Debug.h"
#include "Log.h"
#include "GraphException.h"

using namespace std;

/**
 * Parse all of our options
 * \param[in] argc argument count
 * \param[in] argv arguments
 * \param[out] infile file to read as input
 * \param[out] intype type of intput file
 * \param[out] outfile file to write output to
 * \param[out] methods list of methods we want to run
 */
int parse_options(int argc, char **argv, string& infile, string& intype, string& outfile, vector<string>& methods){
    int flags, opt;
    while((opt = getopt(argc, argv, "i:t:o:m:")) != -1){
        switch(opt){
        case 'i':
            infile = optarg;
            break;
        case 't':
            intype = optarg;
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'm':
            string tempstr;
            string method;
            string token;
            tempstr = optarg;
            stringstream convert(tempstr);
            /* split the string on commas */
            while(getline(convert, token, ',')){
                methods.push_back(token);
            }
            break;
        }
    }

    return 0;
} // parse_options

int main(int argc, char **argv){
    string infile;
    string outfile;
    string intype;
    vector<string> methods;
    parse_options(argc, argv, infile, intype, outfile, methods);
    cout << "done parsing options\n";
    cout << "Input  file: " << infile << "\n";
    cout << "Input  type: " << intype << "\n";
    cout << "Output file: " << outfile << "\n";
    cout << "Methods    :";
    for(vector<string>::iterator it = methods.begin(); it != methods.end(); ++it){
        cout << " " << *it;
    }
    cout << "\n";
}

