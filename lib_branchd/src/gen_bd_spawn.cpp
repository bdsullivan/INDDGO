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
    fprintf(stderr,"Usage: %s -f DIMACS_file -p eig_prob <-v>\n"
        "\t-p eig_prob will have p*deg(v) edges on each side of the split\n"
        "\t\t (default is p=.3333)\n"
        "\t-nopush will not use pushing\n"
        "\t-initpush will do only the initial pushing if nopush is set\n"
        "\t-no2sep will not use 2-separators\n"
        "\t-seed seed will effectively seed the RNG\n"
        "\t-root u will root the BD by breaking the edge u and adding new edges and nodes\n"
        "\t-subtree s will print out a list of the nodes in the subtree rooted at s\n"
        "\t-exhaust max_bits\n"
        "\t-v will run in verbose mode\n\n"
        ,str);
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

    time_t global_start=clock();

    int i,j,new_node, num_pushes=0, seed=0, nsplits=0, root=-1,subtree_root=-1,new_ms,max_bits=-1,
        num_spawns=0;
    // Set default p to 1/3
    double p=0.3333; 
    list<int>::iterator ii;
    char *DIMACS_file=NULL;
    bool has_graph=false,verbose=false,do_push=true,do_two_sep=true, do_init_push=false;
    Graph *G=NULL;
    time_t start,stop;
    //for printing out the graphviz representations of each move on the bd 
    bool gviz_all = false; 
    bool gviz_el = false; //edge labels on
    bool gviz_th = true;  //thicknesses off
    char gvizfile[20];
    int step = 0;

    //    time_t begin=clock();
    // Controller object has some methods which are used independently from
    // graph objects. 
    goblinController *CT;
    CT = new goblinController();
    //Turn off printing of dots.
    CT->traceLevel = 0;

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
            G->write_graphviz_file("orig.gviz");
            char metis_file[100];
            sprintf(metis_file,"%s.metis",argv[i+1]);
            G->write_METIS_file(metis_file);
            sprintf(metis_file,"%s.hmetis",argv[i+1]);
            G->write_HMETIS_file(metis_file);
            has_graph=true;
            if(verbose)
            {
                print_message(0,"Read in graph from file: %s\n",DIMACS_file);
                cout << *G;
            }
        }
        if(strcmp(argv[i],"-p")==0)
            p=atof(argv[i+1]);
        if(strcmp(argv[i],"-nopush")==0)
            do_push=false;
        if(strcmp(argv[i],"-no2sep")==0)
            do_two_sep=false;
        if(strcmp(argv[i],"-seed")==0)
            seed=atoi(argv[i+1]);
        if(strcmp(argv[i],"-root")==0)
            root=atoi(argv[i+1]);
         if(strcmp(argv[i],"-subtree")==0)
            subtree_root=atoi(argv[i+1]);
         if(strcmp(argv[i],"-exhaust")==0)
             max_bits=atoi(argv[i+1]);
    }
    // Make sure we have a graph
    if(!has_graph)
        fatal_error("Did not load graph\n");

    if(!G->check_connected())
        fatal_error("Graph is not connected\n");

    if(!G->check_two_connected())
        fatal_error("Graph is not 2 connected!\n");

    // "seed" the rng
    for(i=0;i<seed;i++)
        lcgrand(0);
    // Create the tree
    BDTree btree(G);
    // Initialize to the star configuration 
    btree.create_star();

    if(gviz_all)
    {
        sprintf(gvizfile, "bd.%d.gviz", step);
        btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
        step++;
    }

    // Find candidate pushes
    if(do_init_push)
    {
        num_pushes=btree.push(0);
        print_message(0,"Found %d initial pushes\n",num_pushes);
    }

    if(gviz_all && num_pushes > 0)
    {
        sprintf(gvizfile, "bd.%d.gviz", step);
        btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
        step++;
    }

    if(do_two_sep)
    {
        // Look for 2-separations 
        list<int> X1, Y1;
        bool ts = false;
        int bv = 0;
        bool found;
        while(!ts)
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
            if(!found)
            {
                print_message(0, "Did not find a(nother) node of btree to split.\n");
                break;
            }

            if(verbose)
                print_message(0,"Running two separation function at vertex %d\n", bv);
            X1.clear();
            Y1.clear();
            start = clock();
	    ts= btree.two_separation(bv, &X1, &Y1, CT);
            //ts = btree.bf_two_separation(bv, &X1, &Y1);
            stop = clock(); 
            print_message(0, "Checked for two separation in %f seconds.\n",
                ((double)(stop-start))/CLOCKS_PER_SEC);    
            if(ts)
            {
                print_message(0, "Found a valid 2 separation!\n");
                print_message(0, "Splitting node %d\n", bv);
                // Split the node
                new_node = btree.split_node(bv,&X1,EDGE_SPLIT,&new_ms);
                print_message(0,"After split - new edge has middle set of size %d\n",  new_ms);

                if(gviz_all)
                {
                    sprintf(gvizfile, "bd.%d.gviz", step);
                    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
                    step++;
                }

                if(do_push)
                {                    
                    //Push the vertex you split
                    num_pushes=btree.push(bv);
                    print_message(0, "Found %d pushes at %d \n",num_pushes, bv);

                    if(gviz_all && num_pushes > 0)
                    {
                        sprintf(gvizfile, "bd.%d.gviz", step);
                        btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
                        step++;
                    }


                    //Push the new vertex created 
                    num_pushes=btree.push(new_node);
                    print_message(0, "Found %d pushes at %d \n",num_pushes, new_node);

                    if(gviz_all && num_pushes > 0)
                    {
                        sprintf(gvizfile, "bd.%d.gviz", step);
                        btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
                        step++;
                    }


                }
                //Add one to our count, then reset ts and bv so we restart search for 2-seps.

                nsplits++;
                ts = false;
                bv = 0;
            }
            else
            {
                print_message(0, "No 2 separation at vertex %d!\n", bv);
                bv++;
            }
        }
        print_message(0, "Finished with two-separations. Split %d nodes.\n", nsplits);
    }

    // Now run eigenvector splitting until BD is valid    
    int split_node=0,new_node_1, new_node_2;
    nsplits=0;
    list<int> A, B, CA, CB, partition, partition2;
    vector<int> candidates(btree.num_nodes);
    while(!btree.is_valid)
    {
        // Find a BDTreeNode with at least 4 neighbors
        // fill candidates with possibilities
        // This is not smart since we really should just update
        // the candidates as we split and add nodes...
        // but it's probably in the noise anyway
        j=0;
        split_node = -1;
        for(i=0;i<btree.num_nodes;i++)
        {
            if(btree.nodes[i].edges.size()>=4)
            { 
                candidates[j]=i;
                j++;          
            }
        }
      
        print_message(1,"generating random int in 0...%d\n",j-1);
        split_node=rand_int(0,j-1);
        split_node=candidates[split_node];
        print_message(1,"split_node is %d (degree=%d)\n",split_node,btree.nodes[split_node].edges.size());
        
        A.clear(); B.clear(); CA.clear(); CB.clear();
        nsplits++;	

        // Run eigenvector heuristic - if fill_extra=true then we will 
        // we are adding "obvious" edges to A and B within this function!
        start=clock(); 
        if(btree.num_interior_nodes==1)
            btree.eigenvector_split(split_node,&A, &B, p, true);
        else
        {
            btree.eigenvector_leaf_split(split_node,&A, &B, p, true);
#if 0
            // Tried to force interior edge to be on one side of split - appears to 
            // be much worse than the leaf split but I left the code in here...
            btree.eigenvector_split(split_node,&A, &B, p, true);
            // Add the interior edge to A or B before running max flow
            int interior_edge=-1;
            if(btree.num_interior_nodes!=1)
            {
                for(ii=btree.nodes[split_node].edges.begin();ii!=btree.nodes[split_node].edges.end();++ii)
                {
                    if(btree.edges[*ii].start==split_node)
                        j=btree.edges[*ii].end;
                    else
                        j=btree.edges[*ii].start;
                    if(btree.nodes[j].graph_edge_start==-1 && btree.nodes[j].graph_edge_end==-1)
                    {
                        interior_edge=*ii;
                        break;
                    }
                }
                if(interior_edge==-1)
                    fatal_error("%s:  couldn't find interior edge??\n",__FUNCTION__);
            }
            list<int> A_ms,B_ms;
            btree.middle_set_union(&A,&A_ms);
            btree.middle_set_union(&B,&B_ms);
            A_ms.sort(); B_ms.sort();
            vector<int> ms_intersection(A_ms.size()+B_ms.size(),0);
            vector<int>::iterator ms_vv;
            btree.edges[interior_edge].middle_set.sort();
            ms_vv=set_intersection(A_ms.begin(),A_ms.end(),btree.edges[interior_edge].middle_set.begin(),
                btree.edges[interior_edge].middle_set.end(),ms_intersection.begin());
            int A_inter_size=ms_vv-ms_intersection.begin();
            ms_vv=set_intersection(B_ms.begin(),B_ms.end(),btree.edges[interior_edge].middle_set.begin(),
                btree.edges[interior_edge].middle_set.end(),ms_intersection.begin());
            int B_inter_size=ms_vv-ms_intersection.begin();
            if(A_inter_size<=B_inter_size)
                A.push_back(interior_edge);
            else
                B.push_back(interior_edge);
#endif
        }

        stop=clock();
        print_message(0,"Computed eigenvector with p=%f in %f seconds.\n",p,
            ((double)(stop-start))/CLOCKS_PER_SEC);
        print_message(0,"Eigenvector A:\n");
        print(0,A);
        print_message(0,"Eigenvector B:\n");
        print(0,B);

        // Should do this only if A and B require it 
        if(A.size()+B.size() < btree.nodes[split_node].edges.size())
        {
            // Run the max flow to get an actual splitting of the edges
            start=clock(); btree.spawn_maxflow_partition(split_node,&A,&B,&CA,&CB, CT); stop = clock(); 
            print_message(0,"Computed partition given eigenvector results in %f seconds.\n",
                ((double)(stop-start))/CLOCKS_PER_SEC);
           
            // Create  the edge set
            partition.clear();
            for(ii=A.begin();ii!=A.end();++ii)
                partition.push_back(*ii);
            for(ii=CA.begin();ii!=CA.end();++ii)
                partition.push_back(*ii);
            partition.sort();
        }
        else
        {
            // Just use A
            partition.clear();
            for(ii=A.begin();ii!=A.end();++ii)
                partition.push_back(*ii);
        }
        
        // Check the size of the exhaust BEFORE splitting the node
        int exhaust_bits=btree.nodes[split_node].edges.size() - A.size() - B.size();
        print_message(1,"Exhaust size is %d bits\n",exhaust_bits);

        int exhaust_ms_size=-1;
        time_t exh_start=0, exh_stop=0;
        list<int> C;
        // Check to see if we can/want to exhaust
        if( exhaust_bits <= max_bits)
        {
            print_message(0,"Checking exhaust\n");
            // Check this exhaustively
            exh_start=clock();
            C.clear();
            exhaust_ms_size=btree.best_partition(split_node,&A, &B, &C);
            exh_stop=clock();         
        }

        
        print_message(1,"Splitting %d with edge partition of size %d\n",split_node,partition.size());
        print(1,partition);
        if(btree.num_interior_nodes==1)
        {
            // initial star - use split node
            new_node_2=-1;
            new_node_1 = btree.split_node(split_node,&partition,EDGE_SPLIT,&new_ms);
            btree.write_graphviz_file(true,"init.gviz",true,true);
        }
        else
        {
            // Later on in the process - use spawn node
            int ms1,ms2;
            btree.spawn_nodes(split_node,&partition,EDGE_SPLIT,&ms1,&ms2, &new_node_1, &new_node_2);
            print_message(0,"After spawn - new middle sets of size %d,%d (max ms is %d)\n",ms1,ms2,btree.max_middle_set_size);
            print_message(0,"new_nodes: %d,%d\n",new_node_1, new_node_2);
            //char spawn_file[100];
            //sprintf(spawn_file,"spawn_%d.gviz",num_spawns);
            //btree.write_graphviz_file(true,spawn_file,true,true);
            num_spawns++;
        }

        // Check for validity here!!!
        if(btree.is_valid)
            break;

        // CSG - this fails because we had a spawn of a leaf that was an old node -
        // so new node1 and 2 are getting set incorrectly in spawn_node
        if(new_node_1!=-1)
        {
            // Try to push the newly created nodes
            if(btree.nodes[new_node_1].edges.size()>=4 && do_push)
            {
                num_pushes=btree.push(new_node_1);
                print_message(0,"\n\tFound %d pushes for newly introduced node %d\n",num_pushes,new_node_1);
                if(gviz_all && num_pushes > 0)
                {
                    sprintf(gvizfile, "bd.%d.gviz", step);
                    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
                    step++;
                }
            }
            // Check for validity here!!!
            if(btree.is_valid)
                break;
        }

        if(new_node_2!=-1)
        {
            if(btree.nodes[new_node_2].edges.size()>=4 && do_push)
            {
                num_pushes=btree.push(new_node_2);
                print_message(0,"\n\tFound %d pushes for newly introduced node %d\n",num_pushes,new_node_2);
                if(gviz_all && num_pushes > 0)
                {
                    sprintf(gvizfile, "bd.%d.gviz", step);
                    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
                    step++;
                }
            }
        }
        // Check for validity here!!!
        if(btree.is_valid)
            break;
    }

    if(root!=-1)
    {
        // This seems to be working when all graphviz flags are off, but something doesn't seem right
        // if last param is set to 2, probably because middle set is empty and we get 0 pen width??!
        // BDS - fixed this by making minimum penwidth 1 (i.e. pw = |mid set| unless |mid set| = 0, in which 
        // case, pw = 1. 
        btree.write_graphviz_file(false,"before_root.gviz",false,true);
        //cout<<btree;
        btree.root(root);
        btree.write_graphviz_file(false,"after_root.gviz",false,true);
        //cout<<btree;
    }

    if(verbose)
        cout<<btree;

    print_message(0,"%d splits performed\n",nsplits);
    int max_ms=0;
    for(i=0;i<btree.num_edges;i++)
        if((int)btree.edges[i].middle_set.size()>max_ms)
            max_ms=btree.edges[i].middle_set.size();
    print_message(0,"max middle set is %d\n",max_ms);

    vector<int> hist(max_ms+1,0);
    for(i=0;i<btree.num_edges;i++)
        hist[btree.edges[i].middle_set.size()]++;
    print(0,hist);

    if(gviz_all && num_pushes > 0)
    {
        sprintf(gvizfile, "bd.%d.gviz", step);
        btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
        step++;
    }

    time_t global_stop=clock();
    printf("%s %3.3f %3.3f %d %d\n",DIMACS_file,p,(double)(global_stop-global_start)/CLOCKS_PER_SEC,nsplits,max_ms);
    fflush(stdout);
    
    delete G;

    delete CT;

    return 1;
}

