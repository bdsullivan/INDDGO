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
	    "\t-bfpush will use a (slow) brute force push routine. This overrides -nopush and -initpush options.\n"
	    "\t-bfpushOUT will use a brute force push routine that only considers pushes where the pair of special edges"
	    "\t includes the new edge created in the previous split. This overrides all other push-related options.\n"
	    "\t-no2sep will not use 2-separators\n"
	    "\t-no3sep will not use 3-separators\n"	    
	    "\t-seed seed will effectively seed the RNG\n"
	    "\t-root u will root the BD by breaking the edge u and adding new edges and nodes\n"
	    "\t-subtree s will print out a list of the nodes in the subtree rooted at s\n"
	    "\t-exhaust max_bits\n"
	    "\t-brute2 will try to use brute force to find a 2 separation (very very slow)\n"
        "\t-brute3 will try to use brute force to find a 3 separation (very very slow)\n"
	    "\t-isort will try to split nodes with the fewest # of interior edges first\n"
	    "\t-dsort will try to split the nodes with highest degree first\n"
	    "\t-read_bwidth <tree_file> will read in a file produced by Bill Cook's bwidth binary\n"
        "\t\t and write tree_file.gviz graphviz file (you still have to provide DIMACs file too)\n"
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

    int i,j,new_node, num_pushes=0, seed=0, nsplits=0, root=-1,subtree_root=-1,new_ms,max_bits=-1;
    int newe;
    // Set default p to 1/3
    double p=0.3333; 
    list<int>::iterator ii;
    char *DIMACS_file=NULL, *bwidth_file=NULL;
    bool has_graph=false,verbose=false,do_push=true,do_two_sep=true, do_init_push=false, do_brute2=false, 
        do_brute3=false, do_three_sep = true, isort=false, dsort=false, bf_push = false, bf_pushOUT = false, read_bwidth=false;
    Graph *G=NULL;
    time_t start,stop;
    //for printing out the graphviz representations of each move on the bd 
    bool gviz_all = false; 
    bool gviz_el = true;  //edge labels on
    bool gviz_th = true;  //thicknesses off
    char gvizfile[20];
    int step = 0;

    time_t global_start=clock();
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
	    if(strcmp(argv[i],"-bfpush")==0)
		bf_push=true;
	    if(strcmp(argv[i],"-bf_pushOUT")==0)
		bf_pushOUT=true;
	    if(strcmp(argv[i],"-no2sep")==0)
		do_two_sep=false;
	    if(strcmp(argv[i],"-no3sep")==0)
		do_three_sep=false;
	    if(strcmp(argv[i],"-brute2")==0)
		do_brute2=true;
        if(strcmp(argv[i],"-brute3")==0)
		do_brute3=true;
	    if(strcmp(argv[i],"-seed")==0)
		seed=atoi(argv[i+1]);
	    if(strcmp(argv[i],"-root")==0)
		root=atoi(argv[i+1]);
	    if(strcmp(argv[i],"-subtree")==0)
		subtree_root=atoi(argv[i+1]);
	    if(strcmp(argv[i],"-exhaust")==0)
		max_bits=atoi(argv[i+1]);
	    if(strcmp(argv[i],"-isort")==0)
		isort=true;
	    if(strcmp(argv[i],"-dsort")==0)
		dsort=true;
	    if(strcmp(argv[i],"-read_bwidth")==0)
		{
		    read_bwidth=true;
		    bwidth_file=argv[i+1];
		}
	}
    // Make sure we have a graph
    if(!has_graph)
        fatal_error("Did not load graph\n");

    if(!G->check_connected())
        fatal_error("Graph from file %s is not connected\n",DIMACS_file);

    if(!G->check_two_connected())
        fatal_error("Graph from file %s is not 2 connected!\n",DIMACS_file);

    // "seed" the rng
    for(i=0;i<seed;i++)
        lcgrand(0);
    // Create the tree
    BDTree btree(G);

    if(read_bwidth)
    {
        btree.read_bwidth_file(bwidth_file);
        char bwidth_gviz_file[100];
        sprintf(bwidth_gviz_file,"%s.gviz",bwidth_file);
        btree.write_graphviz_file(false,bwidth_gviz_file,false,false);
        fprintf(stderr,"bwidth tree %s written in GVIZ format to %s\n",bwidth_file,bwidth_gviz_file);
        exit(-1);
    }
    // Initialize to the star configuration 
    btree.create_star();

    if(do_brute2)
	{    
	    bool has_2sep=btree.brute_force_two_sep();
	    if(has_2sep)
		print_message(0,"2 separation found for graph %s!!\n",DIMACS_file);
	}

    if(do_brute3)
	{    
	    bool has_3sep=btree.brute_force_three_sep();
	    if(has_3sep)
            print_message(0,"3 separation found for graph %s!!\n",DIMACS_file);
	}


    if(gviz_all)
	{
	    sprintf(gvizfile, "bd.%d.gviz", step);
	    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
	    step++;
	}

    // Find candidate pushes
    if(bf_push || bf_pushOUT)
	{
	    num_pushes = btree.doall_bfpushes(0, false, 0);
	    print_message(0,"BRUTE FORCE: Found %d initial pushes\n",num_pushes);
	}
    else if(do_init_push || do_push)
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
	    //until bf_two_separation is finished, we need this to keep from splitting while 
	    //lists are not populated!
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
			    new_node = btree.split_node(bv,&X1,EDGE_SPLIT,&new_ms, &newe);
			    print_message(0,"After 2sep split - new edge has middle set of size %d\n",  new_ms);

			    if(gviz_all)
				{
				    sprintf(gvizfile, "bd.%d.gviz", step);
				    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
				    step++;
				}
		
			    if(bf_push || bf_pushOUT)
				{
				    if(bf_pushOUT)
					num_pushes = btree.doall_bfpushes(bv, true, newe);
				    else
					num_pushes = btree.doall_bfpushes(bv, false, 0);
				    print_message(0, "BRUTE FORCE: Found %d pushes at %d \n",num_pushes, bv);
				    if(gviz_all && num_pushes > 0)
					{
					    sprintf(gvizfile, "bd.%d.gviz", step);
					    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
					    step++;
					}
			
				    //Push the new vertex created 
				    if(bf_pushOUT)
					num_pushes = btree.doall_bfpushes(new_node, true, newe);
				    else
					num_pushes=btree.doall_bfpushes(new_node, false, 0);
				    print_message(0, "BRUTE FORCE: Found %d pushes at %d \n",num_pushes, new_node);
			
				    if(gviz_all && num_pushes > 0)
					{
					    sprintf(gvizfile, "bd.%d.gviz", step);
					    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
					    step++;
					}
			
				}//end of brute force
			    else if(do_push)
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

    if(do_three_sep)
	{
	    if(!do_two_sep)
		{
		    print_message(0, "Cannot look for three separations with two separations disabled. Rerun without -no2sep flag.\n");
		}
	    else
		{
		    // Look for 3-separations 
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
				print_message(0,"Running three separation function at vertex %d, which has degree %d\n", bv, btree.nodes[bv].edges.size());
			    X1.clear();
			    Y1.clear();
			    start = clock();
			    ts= btree.three_separation(bv, &X1, &Y1, CT);
			    stop = clock(); 
			    print_message(0, "Checked for three separation in %f seconds.\n",
					  ((double)(stop-start))/CLOCKS_PER_SEC);    
			    if(ts)
				{
				    print_message(0, "Found a valid 3 separation!\n");
				    print_message(0, "Splitting node %d\n", bv);
				    // Split the node
				    new_node = btree.split_node(bv,&X1,EDGE_SPLIT,&new_ms, &newe);
				    print_message(0,"After 3sep split - new edge has middle set of size %d\n",  new_ms);

				    if(gviz_all)
					{
					    sprintf(gvizfile, "bd.%d.gviz", step);
					    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
					    step++;
					}


				    if(bf_push || bf_pushOUT)
					{
					    if(bf_pushOUT)
						num_pushes = btree.doall_bfpushes(bv, true, newe);
					    else
						num_pushes = btree.doall_bfpushes(bv, false, 0);
					    print_message(0, "BRUTE FORCE: Found %d pushes at %d \n",num_pushes, bv);
					    if(gviz_all && num_pushes > 0)
						{
						    sprintf(gvizfile, "bd.%d.gviz", step);
						    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
						    step++;
						}
			
					    //Push the new vertex created 
					    if(bf_pushOUT)
						num_pushes = btree.doall_bfpushes(new_node, true, newe);
					    else
						num_pushes=btree.doall_bfpushes(new_node, false, 0);
					    print_message(0, "BRUTE FORCE: Found %d pushes at %d \n",num_pushes, new_node);
			
					    if(gviz_all && num_pushes > 0)
						{
						    sprintf(gvizfile, "bd.%d.gviz", step);
						    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
						    step++;
						}
			
					}//end of brute force
				    else if(do_push)
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
				    print_message(0, "No 3 separation at vertex %d!\n", bv);
				    bv++;
				}
			}
		    print_message(0, "Finished with three-separations. Split %d nodes (total).\n", nsplits);
		}
	}

    // Now run eigenvector splitting until BD is valid    
    int split_node=0;
    nsplits=0;
    list<int> A, B, CA, CB, partition, partition2;
    // Make candidates an int_int - if isort=false, only use the first data entry
    vector<int_int> candidates(btree.num_nodes);
    vector<int> sorted_candidates(btree.num_nodes);
    
    while(!btree.is_valid)
	{
	    // Find a BDTreeNode with at least 4 neighbors
	    // fill candidates with possibilities
	    // This is not smart since we really should just update
	    // the candidates as we split and add nodes...
	    // but it's probably in the noise anyway
	    j=0;
	    split_node = -1;
	    int min_interior_edges = GD_INFINITY;
	    int max_degree= -1;
	    for(i=0;i<btree.num_nodes;i++)
		{
		    if(btree.nodes[i].edges.size()>=4)
			{
			    // Shouldn't have to do this but something not right elsewhere??
			    btree.calculate_num_leaf_nbrs(i);

			    candidates[j].p1=i;
			    if(isort)
				{
				    // Set p2 to be the # of interior edges
				    candidates[j].p2 = (int)(btree.nodes[i].edges.size())-btree.nodes[i].num_leaf_nbrs;
				    // Sanity check
				    if(candidates[j].p2<=0 && nsplits>0)
					{
					    // This shouldn't happen!
					    print_message(0,"node has %d<=0 interior edges??\n",candidates[j].p2);
					    cerr<<btree.nodes[i];
					}
				    if(candidates[j].p2<=min_interior_edges)
					min_interior_edges=candidates[j].p2;
				}
			    if(dsort)
				{
				    // Set p2 to be the degree
				    candidates[j].p2= (int)(btree.nodes[i].edges.size());
				    if(candidates[j].p2>=max_degree)
					max_degree=candidates[j].p2;
				}                
               
			    j++;          
			}
		}
	    if(isort)
		{
		    // pick a random candidate that has min_interior_edges interior edges
		    print_message(1,"min interior edges is %d\n",min_interior_edges);
		    j=0;
		    for(i=0;i<btree.num_nodes;i++)
			{         
			    if(btree.nodes[i].edges.size()>=4 && 
			       (btree.nodes[i].edges.size()-btree.nodes[i].num_leaf_nbrs==min_interior_edges))
				{ 
				    sorted_candidates[j]=i;
				    j++;          
				}
			}
		    split_node=rand_int(0,j-1);
		    split_node=sorted_candidates[split_node];
		}
	    if(dsort)
		{
		    print_message(1,"max degree is %d\n",max_degree);
		    j=0;
		    for(i=0;i<btree.num_nodes;i++)
			{         
			    if(btree.nodes[i].edges.size()>=4 && 
			       (btree.nodes[i].edges.size()==max_degree))
				{ 
				    sorted_candidates[j]=i;
				    j++;          
				}
			}
		    split_node=rand_int(0,j-1);
		    split_node=sorted_candidates[split_node];

		}
	    if(!isort && !dsort)
		{
		    // pick a random candidate
		    print_message(1,"generating random int in 0...%d\n",j-1);
		    split_node=rand_int(0,j-1);
		    split_node=candidates[split_node].p1;
		    print_message(1,"split_node is %d (degree=%d)\n",split_node,btree.nodes[split_node].edges.size());
		}
	    A.clear(); B.clear(); CA.clear(); CB.clear();
	    nsplits++;	

#if 0
	    // Try multiple splits - just get the entire sorted list implied by eigenvector
	    list<int> V;
	    btree.eigenvector_list(split_node, &V);
	    // Copy A to a vector - this should really be in the above function
	    vector<int> V_vec(V.size(),0);
	    i=0;
	    for(ii=V.begin();ii!=V.end();++ii)
		{
		    V_vec[i]=*ii;
		    i++;
		}
	    int best_ms=999,current_ms,num_each_side,dim=btree.nodes[split_node].edges.size();
	    double q,best_q;
	    for(i=1;i<=19;i++)
		{
		    // Try a bunch of q's and pick the best one
		    q=(double)i/20;
		    num_each_side = (int)ceil((double)dim * q);
		    if(num_each_side<=1)
			num_each_side=2;
		    A.clear(); B.clear();
		    for(j=0;j<num_each_side;j++)
			{
			    A.push_back(V_vec[j]);
			    B.push_back(V_vec[dim-1-j]);
			}
		    if(A.size()+B.size() < btree.nodes[split_node].edges.size())
			{
			    // Run the max flow to get an actual splitting of the edges
			    CA.clear();CB.clear();
			    current_ms=btree.split_maxflow_partition(split_node,&A,&B,&CA,&CB, CT); 
			    if(current_ms<best_ms)
				{
				    best_q=q; best_ms=current_ms;
				}
			}
		}
	    print_message(0,"Best q=%f\n",best_q);
	    // Restore the best partition found (best_q)
	    q=best_q;
	    num_each_side = (int)ceil((double)dim * q);
	    if(num_each_side<=1)
		num_each_side=2;
	    A.clear(); B.clear();
	    for(j=0;j<num_each_side;j++)
		{
		    A.push_back(V_vec[j]);
		    B.push_back(V_vec[dim-1-j]);
		}
	    if(A.size()+B.size() < btree.nodes[split_node].edges.size())
		{
		    // Run the max flow to get an actual splitting of the edges
		    CA.clear();CB.clear();
		    current_ms=btree.split_maxflow_partition(split_node,&A,&B,&CA,&CB, CT); 
		    if(current_ms<best_ms)
			{
			    best_q=q; best_ms=current_ms;
			}
		}
	    // Create  the edge set
	    partition.clear();
	    for(ii=A.begin();ii!=A.end();++ii)
		partition.push_back(*ii);
	    for(ii=CA.begin();ii!=CA.end();++ii)
		partition.push_back(*ii);
	    partition.sort();
#else
        // Run eigenvector heuristic - note that we set fill_extra=true so that
        // we are adding "obvious" edges to A and B within this function!
        start=clock();   btree.eigenvector_split(split_node,&A, &B, p, false);   stop=clock();
        print_message(1,"Computed eigenvector with p=%f in %f seconds.\n",p,
            ((double)(stop-start))/CLOCKS_PER_SEC);
        print_message(1,"Eigenvector A:\n");
        print(1,A);
        int predicted_width = -1;

	    if(verbose)
		print_message(0, "size of A: %d, size of B: %d, num edges: %d\n", A.size(), B.size(), btree.nodes[split_node].edges.size());

	    // Should do this only if A and B require it 
	    if(A.size()+B.size() < btree.nodes[split_node].edges.size())
		{
		    // Run the max flow to get an actual splitting of the edges
		    start=clock();
		    // CSG - sep. 30 - try to run many maxflows and find the best while A and B have at least 2 nodes
		    //BDS  - changed to partition function (not partition2)
		    predicted_width= btree.split_maxflow_partition(split_node,&A,&B,&CA,&CB, CT); 

		    stop = clock();
		    if(verbose)
			print_message(0,"Step %d:  Computed partition given eigenvector results in %f seconds.\n",step,

				      ((double)(stop-start))/CLOCKS_PER_SEC,step);

		    // Create  the edge set
		    partition.clear();
		    for(ii=A.begin();ii!=A.end();++ii)
			partition.push_back(*ii);
		    for(ii=CA.begin();ii!=CA.end();++ii)
			partition.push_back(*ii);
		    partition.sort();

		    if(0)
			{
			    partition2.clear();
			    for(ii=B.begin();ii!=B.end();++ii)
				partition2.push_back(*ii);
			    for(ii=CB.begin();ii!=CB.end();++ii)
				partition2.push_back(*ii);
			    partition2.sort();
			    list<int> m1, m2;
			    btree.middle_set_union(&partition, &m1);
			    btree.middle_set_union(&partition2, &m2);
			    print_message(0, "Middle set Union of A (size=%d)\n",m1.size());
			    print(0, m1); 
			    print_message(0, "Middle set Union of B(size=%d)\n",m2.size());
			    print(0, m2);
			}

		    print_message(1,"Splitting %d with partition of size %d (deg is %d)\n",
				  split_node, partition.size(), btree.nodes[split_node].edges.size());

		    CA.sort();
		    print_message(1,"MaxFlow additions to A:\n");
		    print(1,CA);

		    print_message(1,"Splitting %d with partition of size %d (deg is %d)\n",
				  split_node, partition.size(), btree.nodes[split_node].edges.size());	    

		}
	    else
		{
		    // Just use A
		    partition.clear();
		    for(ii=A.begin();ii!=A.end();++ii)
			partition.push_back(*ii);
		}
#endif

	    // Check the size of the exhaust BEFORE splitting the node
	    int exhaust_bits=btree.nodes[split_node].edges.size() - A.size() - B.size();
	    print_message(1,"Exhaust size is %d bits\n",exhaust_bits);

	    int exhaust_ms_size=-1;
	    time_t exh_start=0, exh_stop=0;
	    list<int> C;
	    // Check to see if we can/want to exhaust
	    if( exhaust_bits <= max_bits)
		{
		    print_message(1,"Checking exhaust\n");
		    // Check this exhaustively
		    exh_start=clock();
		    C.clear();
		    exhaust_ms_size=btree.best_partition(split_node,&A, &B, &C);
		    exh_stop=clock();         
		    print_message(1,"exhaust size is %d\n",exhaust_ms_size);
		}

	    // Now split the node

	    print_message(1,"Splitting %d with edge partition of size %d\n",split_node,partition.size());
	    print(1,partition);
	    int new_node= btree.split_node(split_node,&partition,EDGE_SPLIT,&new_ms, &newe);
       

	    if(verbose)
		print_message(0,"After split - new middle set of size %d\n",new_ms);

	    if(predicted_width!=-1 && predicted_width!=new_ms)
		print_message(0,"%d != %d\n",predicted_width,new_ms);

	    // Compare with the exhaust if we did it
	    if(exhaust_ms_size!=-1)
		{
		    // Print out the middle set sizes found
		    if(exhaust_ms_size != new_ms)
			{
			    // The exhaust allegedly beat max flow
			    print_message(0,"2^%d exhaust: %d,%d [%f secs]\n",exhaust_bits,
					  exhaust_ms_size,new_ms,(double)(exh_stop-exh_start)/CLOCKS_PER_SEC);
			    print_message(0,"exhaust[%d] != maxflow[%d]\n\tEXH: ",exhaust_ms_size,new_ms);
			    C.sort();
			    print(0,C);
			    print_message(0,"\tFLW: ");
			    CA.sort();
			    print(0,CA);
			}
		}            

	    print_message(1,"\t\t\t\tAfter eig.  new middle set %d\n",
			  btree.edges[btree.num_edges-1].middle_set.size());

	    if(gviz_all)
		{
		    sprintf(gvizfile, "bd.%d.gviz", step);
		    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
		    step++;
		}


	    // Check for validity here (ternary tree)!!!
	    if(btree.is_valid)
		break;

	    
	    // Now push the two modified nodes
	    if(btree.nodes[split_node].edges.size()>=4 && (bf_push || bf_pushOUT))
		{	
		    if(bf_pushOUT)
			num_pushes = btree.doall_bfpushes(split_node, true, newe);
		    else
			num_pushes = btree.doall_bfpushes(split_node, false, 0);
		    print_message(0, "BRUTE FORCE: Found %d pushes at %d \n",num_pushes, split_node);
		    if(gviz_all && num_pushes > 0)
			{
			    sprintf(gvizfile, "bd.%d.gviz", step);
			    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
			    step++;
			}
		
		}//end of brute force
	    else if(btree.nodes[split_node].edges.size()>=4 && do_push)
		{
		    num_pushes = btree.push(split_node);
		    print_message(1,"Found %d pushes for recently split node %d\n",num_pushes,split_node);
		    if(gviz_all)
			{
			    sprintf(gvizfile, "bd.%d.gviz", step);
			    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
			    step++;
			}
		}
	    
	    // Check for validity here (ternary tree)!!!
	    if(btree.is_valid)
		break;

		if(btree.nodes[new_node].edges.size()>=4 && (bf_push || bf_pushOUT))
		{
		    if(bf_pushOUT)
			num_pushes = btree.doall_bfpushes(new_node, true, newe);
		    else
			num_pushes = btree.doall_bfpushes(new_node, false, 0);
		    print_message(0, "BRUTE FORCE: Found %d pushes at %d \n",num_pushes, new_node);
		    if(gviz_all && num_pushes > 0)
			{
			    sprintf(gvizfile, "bd.%d.gviz", step);
			    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
			    step++;
			}
		}//end of brute force
	    else if(btree.nodes[new_node].edges.size()>=4 && do_push)
		{
		    num_pushes=btree.push(new_node);
		    print_message(1,"Found %d pushes for newly introduced node %d\n",num_pushes,new_node);
		    if(gviz_all && num_pushes > 0)
			{
			    sprintf(gvizfile, "bd_push.%d.gviz", step);
			    btree.write_graphviz_file(false,gvizfile,gviz_el, gviz_th);
			    step++;
			}

		}
	}
    if(root!=-1)
	{ 
	    btree.write_graphviz_file(false,"before_root.gviz",false,true);
	    //cout<<btree;
	    btree.root(root);
	    btree.write_graphviz_file(false,"after_root.gviz",false,true);
	    //cout<<btree;
	}

    //if(verbose)
    //   cout<<btree;

#if 0
    // This section was just for generating a specific plot when running
    // c:\Users\tcg\PROJECTS\SGD\gaudi\code\trunk\branch_decomposition\Release>BranchDecomposition.exe -p .
    // 33 -no2sep -root 10 -subtree 330 -f ..\data\ch130.tsp.del.100.dimacs
    // Check to see if a subtree is desired
    if(subtree_root!=-1)
	{
	    int roots[24]={400,328,267,292,263,
			   403,438,257,251,302,
			   276,452,294,405,364,
			   379,349,369,330,443,
			   338,420,291,425};
	    char *colors[6]={"red","blue","green","orange","purple","yellow"};
	    list<int> subtree;
	    for(i=0;i<24;i++)
		{
		    subtree.clear();
		    btree.find_subtree(roots[i], &subtree);
		    print_message(1,"Subtree rooted at %d:\n",roots[i]);
		    for(ii=subtree.begin();ii!=subtree.end();++ii)
			{
			    printf("%d [label=\"\",style=filled,fillcolor=%s,color=%s];\n",*ii,colors[i%6],colors[i%6]);
			}
		    printf("\n\n");
		}
	}
#endif


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

    //write a file with thick/thin lines.
    btree.write_graphviz_file(false,"final.gviz",true, false);

    G->write_graphviz_file("orig_graph.gviz");
    //write a file with thick/thin lines and edge labels
    //btree.write_graphviz_file(false,"final.gviz",true, true);

    time_t global_stop=clock();
    //max memory used during entire process.
    int max_mem = getHWmem();
    printf("%s %3.3f %3.3f %d %d %dKB\n",DIMACS_file,p,(double)(global_stop-global_start)/CLOCKS_PER_SEC,nsplits,max_ms, max_mem);

    fflush(stdout);

    delete G;

    delete CT;

    return 1;
}

