#!/usr/bin/env python2

#This file is part of INDDGO.
#
#Copyright (C) 2012, Oak Ridge National Laboratory
#
#This product includes software produced by UT-Battelle, LLC under Contract No.
#DE-AC05-00OR22725 with the Department of Energy.
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the New BSD 3-clause software license (LICENSE).
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#LICENSE for more details.
#
#For more information please contact the INDDGO developers at:
#inddgo-info@googlegroups.com

from optparse import OptionParser
from sys import exit

class graph_properties:
    #add aditional properties here
    properties = ("edges","diameter","eigenspectrum","vertices")

    def __init__(self,graph_path):
        self.items = dict.fromkeys(self.properties)
        self.process_graph(graph_path)

    def process_graph(self, graph_path):
        raw_file = []
        with open(graph_path, 'r') as graph_file:
            for line in graph_file:
                processed_line = line.partition("#")[0].strip()
                if len(processed_line) > 0:
                    raw_file.append(processed_line)

        for line in raw_file:
            split_up_line = line.partition(" ")
            if split_up_line[0] in self.items:
                self.items[split_up_line[0]] = split_up_line[2]
            else:
                print "Tried to add unknown element, line in file was: ", line

    def to_string(self):
        print_string = []
        for element in self.properties:
            if self.items[element] is not None:
                print_string.append(str(self.items[element]))
            else:
                print_string.append("")
        return ",".join(print_string)

    def properties_string(self):
        return ",".join(self.properties)

usage_string = "usage: %prog [-h] [-o OUTPUT] input [input ...]"
parser = OptionParser(usage=usage_string)
parser.add_option("-o", "--output", action="store", type="string", dest="output", help="output file name")
(options, args) = parser.parse_args() 
if(len(args) < 1):
    parser.error("too few arguments")
    exit()

graphs = []
for graph_name in args:
    graphs.append(graph_properties(graph_name))

if options.output is None:
    print graphs[0].properties_string()
    for graph in graphs:
        print graph.to_string()
else:
    with open(options.output,"w") as out_file:
        out_file.write(graphs[0].properties_string())
        out_file.write("\n")
        for graph in graphs:
            out_file.write(graph.to_string())
            out_file.write("\n")
