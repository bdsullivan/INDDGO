#Installing INDDGO

To install INDDGO, follow the steps listed in the INSTALL file.

As described in INSTALL and thirdparty.txt, INDDGO can be compiled with 
several optional packages. The METIS software package is highly recommended
for performance, as it enables several more efficient elimination heuristics, 
as well as a much faster graph triangulation routine.

For visualization, we have provided routines which produce DOT output 
compatible with the graphviz package (http://www.graphviz.org/). These 
files can be passed directly to any one of the layout executables (such as 
neato or dot) after graphviz is installed.

Note that although it is not a prerequisite for using INDDGO, if you 
wish to use the optional functionality of solving mixed integer programming 
formulations of MWIS, you should have glpsol installed and available.

#Running INDDGO binaries

The main INDDGO binaries are placed in the "bin" subdirectory of your
source tree.

The default binaries created by INSTALL include executables to:
* generate partial k-trees (./bin/gen_pkt)
* compute elimination orderings, tree decompositions, and maximum weighted independent sets (./bin/serial_wis and/or ./bin/parallel_wis)
* provide various visualization output (./bin/td_viz)

Currently, all tree decomposition functionality is accessible through the 
weighted independent set executables - flags are provided to disable the 
problem-specific dynamic programming portion of the code (e.g. -decompose_only).

All binaries respect the -h flag for printing a comprehensive usage message detailing required input and all options. Note that not all options are supported yet in the parallel_wis executable.

When running parallel versions of INDDGO, you must set the following
environment variables:

* PTD_NUM_THREADS
* OMP_NUM_THREADS
* MAD_NUM_THREADS

NOTE: It is assumed the input graph for weighted independent set is connected. If it is not, 
these executables will find and run on the largest connected component. 
All components will be written to separate .comp files, and the one used
will be identified through a message to standard out. 

##Examples of common use cases
Two sample input graph files, as well as an example tree decomposition
file and elimination order file can be found in the sample_graphs directory.
 All examples use these files. The optimal MWIS value for each is in the 
 comments at the top of the DIMACS file. Additional information on both input and 
 output formatting for the weighted independent set binaries is in max_wis/README.

#Generating Test Data
 Create 3 partial k-trees on 1000 vertices with width 20, 60 percent of 
 full k-tree edges, randomized vertex labels, and names of the form 
 mygraph.1000.20.60_.dimacs
* ./bin/gen_pkt -t 3 -n 1000 -k 20 -p 60 -r -fn ./sample_graphs/mygraph

#Find Tree Decomposition
 Generate a tree decomposition (serial) and save elimination ordering
 and decomposition to file for future use.
* ./bin/serial_wis -f sample_graphs/1dc.64.dimacs -gavril -mind -decompose_only -w sample_graphs/1dc.64.tree -eorder sample_graphs/1dc.64.mind.eorder

#Run MWIS (serial, find set)
 Computes max weighted independent set and saves vertices to solution file 
 named <inputfile>.WIS.sol
* ./bin/serial_wis -f sample_graphs/1dc.64.dimacs -gavril -mind 

#Run MWIS (serial, find objective) 
 Computes objective (set weight) only, but reduces memory use drastically. 
./bin/serial_wis -f sample_graphs/1dc.64.dimacs -gavril -mind -no_reconstruct -del_ch

#Run MWIS (all parallel) 
 Complete run from graph to final solution. Requires PARMETIS.
* export MAD_NUM_THREADS=4
* export OMP_NUM_THREADS=4
* export PTD_NUM_THREADS=4 
* mpirun -n 4 ./bin/parallel_wis -f sample_graphs/1dc.64.dimacs -gavril -mind -pbag -parmetis

#Run MWIS (only DP parallel) 
 Reads EO from file, serial tree decomposition, parallel (MADNESS) DP.  
* export MAD_NUM_THREADS=4
* export OMP_NUM_THREADS=4
* export PTD_NUM_THREADS=4 
* mpirun -n 4 ./bin/parallel_wis -ord sample_graphs/1dc.64.mind.eorder -f sample_graphs/1dc.64.dimacs -gavril 

#Generate visualization
 Create a basic DOT file of the tree decomposition for use with graphviz where
 bags are labelled with the vertices they contain 
 (sample graphviz command: neato -Tpdf -o tree.pdf tree.dot). 
* ./bin/td_viz -f ./sample_graphs/1dc.64.dimacs -t ./sample_graphs/1dc.64.tree -e 

More information on the output formats for the weighted independent set executables is in the README in max_wis.

##Getting Help

To get help with or discuss INDDGO, join the mailing list at inddgo-info@googlegroups.com
