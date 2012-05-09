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

#include "GraphDecomposition.h"
#include "BranchDecomposition.h"
#include <goblin.h>

void usage(char *str)
{
    fprintf(stderr,"Usage: %s -f DIMACS_file <-v>\n-v will run in verbose mode",str);
    return;
}

int main(int argc, char **argv)
{
     if(argc==1 || (argc==2 && strcmp(argv[1],"-h")==0) || (argc==2 && strcmp(argv[1],"--help")==0) 
        || (argc==2 && strcmp(argv[1],"--h")==0) )
    {
        usage(argv[0]);
        exit(-1);
    }

    int i;
    list<int>::iterator ii;
    char *DIMACS_file=NULL;
    bool has_graph=false,verbose=false;
    Graph *G=NULL;
    time_t start, stop;
    int num_pushes;
    
    // Parse arguments
    for(i=0;i<argc;i++)
    {
      if(strcmp(argv[i],"-v")==0)
        verbose=true;

      if(strcmp(argv[i],"-f")==0)
        {
            DIMACS_file=argv[i+1];
            // Read in the file
            G=new Graph(DIMACS_file, true);
           
            has_graph=true;
            if(verbose)
            {
                printf("Read in graph from file: %s\n",DIMACS_file);
                cout << *G;
            }
        }
    }
    // Make sure we have a graph
    if(!has_graph)
        fatal_error("Did not load graph\n");

    if(!G->check_connected())
        fatal_error("Graph is not connected\n");

    if(!G->check_two_connected())
        fatal_error("Graph is not 2 connected!\n");

    //cout << *G;

    BDTree btree(G);
    bool ts;
    btree.create_star();

    if(verbose){
	//cout << "Initial star." << endl;
	//cout<<btree;
	btree.write_graphviz_file(false,"star.gviz",true, false);
    }

    list<int> X;

    // Controller object has some methods which are used independently from
    // graph objects. 
    goblinController *CT;
    CT = new goblinController();
    //Turn off printing of dots.
    CT->traceLevel = 0;


    // if(verbose){
    //   cout << "Running two separation function at vertex 0." << endl;
    // }
     list<int> X1, Y1;
     int new_node;
     // start = clock(); 
    // ts = btree.two_separation(0, &X1, &Y1, CT);
    // stop = clock(); 
    // printf("Checked for two separation in %f seconds.\n",
    // 	   ((double)(stop-start))/CLOCKS_PER_SEC);    
    // if(ts){
    // 	print_message(0, "Found a valid 2 separation!\n");
    // 	print_message(0, "Splitting node 0\n");
    // 	// Split the node
    // 	new_node = btree.split_node(0,&X1,EDGE_SPLIT);
    // }
    // else
    // 	print_message(0, "No 2 separation at vertex 0!\n");
    
    
    // btree.write_graphviz_file(true,"btree_firstsplit.gviz",true, false);

    //Push at both nodes affected by split
    if(btree.nodes[0].edges.size() > 3)
	{
	    num_pushes=btree.push(0);
	    printf("Found %d pushes at 0 \n",num_pushes);
	}
    // if(btree.nodes[new_node].edges.size() > 3)
    // 	{
    // 	    num_pushes=btree.push(new_node);
    // 	    printf("Found %d pushes at %d\n",num_pushes, new_node);
    // 	}

    btree.write_graphviz_file(true,"btree_firstpush.gviz",true, false);
    
    
    X1.clear(); 
    Y1.clear();

    ts = false;
    int bv = 0;
    bool found;
    int ms;
    while(ts == false)
	{
	    //we need to try to separate a vertex with at least 4 edges adjacent to it
	    found = false;
	    while(found == false && bv < btree.num_nodes)
		{
		    if(btree.nodes[bv].edges.size() >3)
			found = true; 
		    else
			bv++;
		} 
	    if(found == false)
		{
		    print_message(0, "Did not find a(nother) node of btree to split.\n");
		    break;
		}
		
		if(verbose)
		    cout << "Running two separation function at vertex " << bv << endl;
		start = clock();
		ts= btree.two_separation(bv, &X1, &Y1, CT);
		stop = clock(); 
		printf("Checked for two separation in %f seconds.\n",
		       ((double)(stop-start))/CLOCKS_PER_SEC);    
		if(ts){
		    print_message(0, "Found a valid 2 separation!\n");
		    print_message(0, "Splitting node %d\n", bv);
		    // Split the node
		    new_node = btree.split_node(bv,&X1,EDGE_SPLIT, &ms);
		    if(btree.nodes[bv].edges.size() > 3)
			{
			    num_pushes=btree.push(bv);
			    printf("Found %d pushes at %d \n",num_pushes, bv);
			}
		    if(btree.nodes[new_node].edges.size() > 3)
			{
			    num_pushes=btree.push(new_node);
			    printf("Found %d pushes at %d \n",num_pushes, new_node);
			}
		    ts = false;
		    bv = 0;
		}
		else
		    {
			print_message(0, "No 2 separation at vertex %d!\n", bv);
			bv++;
		    }
	}
    
    btree.write_graphviz_file(true,"btree_final.gviz",true, false);
	    
    //make valgrind happier - you need to move these to bottom of 
    //function when you get rid of this exit call.
    delete G;

    return 1;
    

    while(!btree.is_valid)
    {
        // Find a node with degree greater than 3 and perform a greedy split
        for(i=0;i<btree.num_nodes && !btree.is_valid;i++)
        {
            if(btree.nodes[i].edges.size()>3)
            {
                X.clear();
                btree.find_greedy_split(i,&X);
                
                print_message(0,"\nSplitting node %d\nX=",i);
                print(0,X);

                for(ii=X.begin();ii!=X.end();++ii)
                    cerr<<btree.edges[*ii];

                btree.split_node(i,&X,NODE_SPLIT, &ms);


                cerr<<btree;
                //exit(-1);
            }

        }
    }
    cout<<btree;

  
    btree.write_graphviz_file(true,"btree.gviz",true, false);

    //we need to destroy the controller differently!!
    delete CT;

}
