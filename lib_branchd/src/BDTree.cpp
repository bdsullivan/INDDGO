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
#define GEN_BD_PICS 0

#define ARPACK_DEBUG 0

/**
* Constructor given a graph *G.
*/
BDTree::BDTree(Graph *G)
{
    // Go ahead and allocate space for the ternary tree
    // CSG added 3 and 2 to edges[] and nodes[] sizes to account for 
    //eventual rooting of tree!
    this->edges=new BDTreeEdge[2*G->get_num_edges()-3 + 3];
    this->nodes=new BDTreeNode[2*G->get_num_edges()-2 + 2];
    for(int i=0;i<(2*G->get_num_edges()-2+2);i++)
        this->nodes[i].tree=this;

    this->num_nodes=0;
    this->num_interior_nodes=0;
    this->is_valid=false;
    this->is_rooted=false;
    this->G=G;

    this->max_middle_set_size=0;
}

/**
* Destructor for the BDTree.
*/
BDTree::~BDTree()
{
    delete [] this->edges;
    delete [] this->nodes;
}

/*
* Overloading the << operator for BDTree.
*/
ostream &operator<<(ostream &output, const BDTree &T)
{
    int i;

    //if(T.graph_file)
    //    output<<"Branch decomposition for graph from file "<<T.graph_file<<endl;
    if(T.is_valid)
        output<<"Decomposition is valid\n\n";
    else
        output<<"Decomposition is not valid!\n\n";

    output<<"Tree has "<<T.num_nodes<<" nodes and "<<T.num_interior_nodes<<" interior nodes\n";

    for(i=0;i<T.num_nodes;i++)
        output<<(T.nodes[i]);

    return output;
}

/**
* Function to write a graphviz representation (DOT language) of the branch decomposition.
* The spline parameter allows the user to specify whether or not spline curves are used 
* in edge layout. This is (strongly) not recommended for large graphs (> 100 vertices).
* If the tree is rooted, this creates digraph output. 
**/

void BDTree::write_graphviz_file(bool spline, const char *gviz_file, bool edge_labels, bool thickness)
{

    int i;
    FILE *out; 
    list<int>::iterator ii,jj;

    if( (out = fopen(gviz_file, "w")) == NULL)
        fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, gviz_file);

    if(this->is_rooted)
        fprintf(out, "Digraph BD{\n");
    else
        fprintf(out, "Graph BD{\n"); 

    if(spline == true)
        fprintf(out, "splines=true;\n");  
    // CSG added Aug. 20 - seems better on smallish (30 node) BD's
    fprintf(out,"overlap=scale\n");

    // Print out the BDTreeNodes
    int nbr;
    list<int> interior_nodes;
    list<int> leaf_nodes;
    int rs, rt;
    for(i=0;i<this->num_nodes;i++)
    {
        if(this->nodes[i].edges.size()==3)
            interior_nodes.push_back(i);
        else
            leaf_nodes.push_back(i);

        for(ii=this->nodes[i].edges.begin(); ii!=this->nodes[i].edges.end();++ii)
        {
            if(this->edges[*ii].start!=i)
                nbr=this->edges[*ii].start;
            else
                nbr=this->edges[*ii].end;

            // We have the edge i-nbr
            if(i<nbr)
            {
                // Print out the edge to the file
                if(this->is_rooted)
                {
                    rs=this->edges[*ii].start;
                    rt=this->edges[*ii].end;
                    fprintf(out,"%d -> %d", rs, rt);
                }//end if_rooted
                else
                {
                    fprintf(out,"%d -- %d",i, nbr);
                }//not rooted

                //begin the options box
                if(thickness || edge_labels)
                {
                    fprintf(out, " [");
                }

                if(thickness || GEN_BD_PICS)
                {
                    int pw=this->edges[*ii].middle_set.size();
                    if(pw == 0)//handle the edge from the root vertex
                        pw = 1;
#if GEN_BD_PICS
                    fprintf(out, "penwidth=20");
#else
                    fprintf(out, "penwidth=%d", pw);
#endif
                    //add a comma if another option
                    if(edge_labels)
                        fprintf(out, ", ");
                }

                if(edge_labels)
                {
                    // Now label the edge with the middle set
                    if(this->edges[*ii].middle_set.size()==0)
                        fprintf(out,"label=\"{}\"");
                    else
                    {
                        // CSG switching to labels
                        fprintf(out," label=\"{");                
                        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                            fprintf(out,"%d ",this->G->nodes[*jj].label);

                        fprintf(out,"}\"");
                    }
                }

                //end the options box
                if(thickness || edge_labels)
                {
                    fprintf(out, "]");
                }

                //always end the line
                fprintf(out, ";\n");

            }
        }
    }

#if GEN_BD_PICS
    for(i=0;i<this->num_nodes;i++)
        fprintf(out,"%d [label=\"\",style=filled,color=black,fillcolor=black];\n",i);
#else
    // Remove the labels on the interior nodes since they don't really mean anything
    for(ii=interior_nodes.begin();ii!=interior_nodes.end();++ii)
        fprintf(out,"%d [label=\"\"];\n",*ii);

    //add labels on the root and leaf nodes
    for(ii = leaf_nodes.begin(); ii != leaf_nodes.end(); ++ii)
    {
        // This is a leaf node - label it with the graph edge it represents
        // CSG switching this to labels
        if(this->nodes[*ii].graph_edge_start > -1)
            fprintf(out,"%d [label=\"(%d,%d)\"];\n", *ii, this->G->nodes[this->nodes[*ii].graph_edge_start].label,
            this->G->nodes[this->nodes[*ii].graph_edge_end].label);
        else // you found the root node.
            fprintf(out, "%d [label=\"root\"];\n", *ii);
    }
#endif


    fprintf(out, "}\n");

    fflush(out);  
    fclose(out); 
}


/**
* Function to create an initial "star" configuration where we have a single central node,
* and m leafs, one for each edge in the graph being decomposed.
*/
void BDTree::create_star()
{
    int i,j=0,k=0;
    list<int>::iterator ii, jj;

    // Loop through the edges in the graph by traversing the node's adjacency lists, adding
    // edges u-v if u<v
    // Use j for edge counting G, k for node counting in the BDTree

    // First create the central BDTreeNode
    this->nodes[0].id=0;
    // Only leaf nodes in the BDTree correspond to edges in the graph
    this->nodes[0].graph_edge_start=-1;
    this->nodes[0].graph_edge_end=-1;

    // We will add the edges incident to the central node as we construct the star
    k++;
    for(i=0;i<this->G->get_capacity();i++)
    {
        if(this->G->nodes[i].label!=-1) // Should be GD_UNDEFINED?
        {
            for(ii=this->G->nodes[i].nbrs.begin();ii!=this->G->nodes[i].nbrs.end();++ii)
            {
                if(i<*ii)
                {
                    // We have a new edge in G of the right form - set the indices of the BDTreeNodes
                    this->edges[j].start=0;
                    this->edges[j].end=k;
                    // Compute the middle set for this edge- I believe that if i and *ii both have degree>=2, then
                    // the middle set is {i,*ii}. Otherwise
                    if(this->G->get_degree(i)>=2 && this->G->get_degree(*ii)>=2)
                    {
                        this->edges[j].middle_set.push_back(i);
                        this->edges[j].middle_set.push_back(*ii);
                    }
                    else
                    {
                        if(this->G->get_degree(i)==1)
                            // *ii has only one neighbor, namely i
                            this->edges[j].middle_set.push_back(*ii);
                        if(this->G->get_degree(*ii)==1)
                            // i has only one neighbor, namely *ii
                            this->edges[j].middle_set.push_back(i);
                    }

                    // Now add the BDTreeNode
                    this->nodes[k].id=k;

                    // BDTreeEdge in position j is adjacent to BDTreeNode k
                    this->nodes[k].edges.push_back(j);
                    // All edges are adjacent to the central BDTreeNode
                    this->nodes[0].edges.push_back(j);
                    // Need to add information about mapping b/w leaf BDTreeNodes and edges in orig. graph 
                    // NOTE: Not using labels here!!
                    this->nodes[k].graph_edge_start=i;
                    this->nodes[k].graph_edge_end=*ii;
                    // record/increment the leaf nbrs
                    this->nodes[k].num_leaf_nbrs=0;
                    this->nodes[0].num_leaf_nbrs++;

                    j++;
                    k++;
                }
            }
        }
    }

    // Sanity check
    if(j!=this->G->get_num_edges())
        fatal_error("%s:  Didn't find the right # of edges! %d != %d\n",__FUNCTION__,j,G->get_num_edges());

    // Set the # of nodes and edges in the tree
    this->num_interior_nodes=1;
    this->num_nodes=k;
    this->num_edges=j;
    this->max_middle_set_size=2;

    return;

}

/**
* Splits the tree node u into two nodes.  A should be a list of at least 2 tree node indices
* that are currently adjacent to u.  After the split, u will still be adjacent to these nodes,
* and letting B denote the set of other neighbors of A, they will now be adjacent to a new node v.
* The new BDTreeEdge u-v is added to the tree.  The value of type should be either NODE_SPLIT or
* EDGE_SPLIT which determines how to interpret the contents of the list A (either node indices
* or edge indices).  Sets the value of new_ms_size to be the size of the middle set of the
* new edge created from the splitting. Returns the new node index.
*/
int BDTree::split_node(int u, list<int> *A, int type, int *new_ms_size)
{
    int i;
    return this->split_node(u, A, type, new_ms_size, &i); 
}

/**
* Splits the tree node u into two nodes.  A should be a list of at least 2 tree node indices
* that are currently adjacent to u.  After the split, u will still be adjacent to these nodes,
* and letting B denote the set of other neighbors of A, they will now be adjacent to a new node v.
* The new BDTreeEdge u-v is added to the tree.  The value of type should be either NODE_SPLIT or
* EDGE_SPLIT which determines how to interpret the contents of the list A (either node indices
* or edge indices).  Sets the value of new_ms_size to be the size of the middle set of the
* new edge created from the splitting. Sets new_e_index to be the index of the newly created tree 
* edge in the edges[] array. 
* Returns the new node index.
*/
int BDTree::split_node(int u, list<int> *A, int type, int *new_ms_size, int *new_e_index)
{
    // Make sure that u currently has 4 or more neighbors in the BDTree
    if(this->nodes[u].edges.size()<=3)
    {
        cerr<<this->nodes[u];
        fatal_error("%s:  trying to split node %d that has degree %d<=3!\n",__FUNCTION__,u,this->nodes[u].edges.size());
    }

    // After the split, u will be adjacent (in the BDTree) to the nodes whose indices are contained in the partition, and
    // the new node will be adjacent to all of the other (old) neighbors of u
    // Old:  A - u - B
    // New:  A - u - v -B
    // where A and b are sets of nodes and v is the new BDTreeNode

    vector<int>::iterator vv;
    list<int>::iterator ii,jj;
    list<BDTreeEdge>::iterator ee;
    list<int> complement, B, new_u_edges;
    int nbr, new_node_index;

    vector<bool> A_vec(this->num_nodes,false);
    vector<bool> A_union_vec(this->G->get_num_nodes(),false);
    vector<bool> B_union_vec(this->G->get_num_nodes(),false);

    list<int> AA;
    // if type is EDGE_SPLIT, then A is a list of edge indices.
    // just populate an alternate list with the indices of the adjacent
    // nodes - then we can use the same code.
    if(type==EDGE_SPLIT)
    {
        print_message(10,"Translating edges - %d nodes, %d edges in BDtree\n",this->num_nodes,this->num_edges);
        // Loop over all the edge indices in A and add the node index to a new list AA
        //print_message(0,"Middle sets of split edges\n");
        for(ii=A->begin();ii!=A->end();++ii)
        {
            //print_message(0,"\t");
            print(1,this->edges[*ii].middle_set);

            print_message(1,"Translating edge %d (%d-%d) into nodes\n",*ii,this->edges[*ii].start,
                this->edges[*ii].end);
            if(this->edges[*ii].start==u)
                AA.push_back(this->edges[*ii].end);
            else
                if(this->edges[*ii].end==u)
                    AA.push_back(this->edges[*ii].start);
                else
                {
                    // This is an error - neither the edge start or end hits the split node u
                    cerr<<*this;
                    fatal_error("%s:  u=%d not found in edge %d (%d-%d)\n",__FUNCTION__,u,
                        *ii,this->edges[*ii].start,this->edges[*ii].end);           
                }
        }
    }
    // AA should now be a node list - now use the pointer C to point to the list we actually 
    // want to use (either the original list A, or the new list just created AA)
    list<int> *C;
    if(type==NODE_SPLIT)
        C=A;
    else
    {
        C=&AA;
        print_message(10,"%s:  Translated initial edge list into node list\n",__FUNCTION__);
        print(10,*A);
        print_message(10,"Nodes:\n");
        print(10,AA);
    }

    //print_message(0,"Splitting %d: sizes are %d,%d\n",u,C->size(),this->nodes[u].edges.size()-C->size());

    // Make sure that u is not in the list C!
    for(ii=C->begin();ii!=C->end();++ii)
        if((*ii)==u)
        {
            print_message(0,"%s:  Split(%d,A) called where %d is in the list of tree nodes A!\n",__FUNCTION__,u,u);
            print(0,*C);
            exit(-1);
        }

        // Create an incidence vector for the node list C
        for(ii=C->begin();ii!=C->end();++ii)
            A_vec[*ii]=true;

        new_node_index=this->num_nodes;
        this->nodes[new_node_index].id=new_node_index;
        // The new node is interior so the graph edge info should be undefined
        this->nodes[new_node_index].graph_edge_start=-1;
        this->nodes[new_node_index].graph_edge_end=-1;

        // Create the set B
        // Loop over all BDTreeEdges incident to u
        for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
        {       
            // *ii is the index in the edges array of the edge incident to u
            // nbr is the non-u endpoint of the edge - check to see if it is in A by looking at A_vec
            nbr=-1;
            if(this->edges[*ii].start != u )
            {
                nbr=this->edges[*ii].start;
                //this->edges[*ii].=new_node_index;
            }
            else
            {
                nbr=this->edges[*ii].end;
                //this->edges[*ii].end=new_node_index;
            }
            // Sanity check
            if(nbr==-1)
                fatal_error("%s:  Did not find node %d in an incident edge! (%d-%d)",
                __FUNCTION__,u,this->edges[*ii].start,this->edges[*ii].end);

            // The middle set of these updated edge should stay the same
            if(!A_vec[nbr])
            {
                // This is an edge in the "B half" - these edges need to be updated as they are no longer
                // incident to u
                B.push_back(nbr);
                for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                    B_union_vec[*jj]=true;

                this->nodes[new_node_index].edges.push_back(*ii);
                print_message(10,"Updating former edge %d: (%d,%d)\n",*ii,this->edges[*ii].start,this->edges[*ii].end);
                this->edges[*ii].start=nbr;//min(nbr,u);
                this->edges[*ii].end=new_node_index;//max(nbr,u);
                // This edge is no longer incident to u - surely there is a better way to remove this...
#if 0
                // This does not work...
                this->nodes[u].edges.erase(ii);
                //this->nodes[u].edges.remove(*ii);
#endif

            }
            else
            {
                // This is an edge in the "A half" - these edges should remain unchanged
                for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                    if(!A_union_vec[*jj])
                        A_union_vec[*jj]=true;

                new_u_edges.push_back(*ii);
            }
        }


        // Update u's list of incident edges - this at least prevents us from having to call remove, but I'm still
        // not elated with this approach...
#if 1
        this->nodes[u].edges.clear();
        this->nodes[u].edges=new_u_edges;
#endif

        // Now add a new BDTreeEdge to the tree
        int new_edge_index=this->num_edges;
        this->edges[new_edge_index].start = u;
        this->edges[new_edge_index].end = new_node_index;
        this->nodes[u].edges.push_back(new_edge_index);
        this->nodes[new_node_index].edges.push_back(new_edge_index);
        this->num_edges++;
	
	*new_e_index = new_edge_index;

        // Compute the middle set of the new edge
        // The middle set of the new edge should be A_union \cap B_union
        //set_intersection(A_union.begin(),A_union.end(),B_union.begin(),B_union.end(),this->edges[new_edge_index].middle_set.begin());
        // CSG - adding this Sep.27 since the split_node middle set appears to occassionally differ from
        // the middle set calculated in maxflow
        print_message(1,"LEFT: ");
        for(int i=0;i<this->G->get_num_nodes();i++)
            if(A_union_vec[i])
                print_message(1,"%d ",i);
        print_message(1,"\nRIGHT: ");
        for(int i=0;i<this->G->get_num_nodes();i++)
            if(B_union_vec[i])
                print_message(1,"%d ",i);
        print_message(1,"\n");
        for(int i=0;i<this->G->get_num_nodes();i++)
            if(A_union_vec[i] && B_union_vec[i])
                this->edges[new_edge_index].middle_set.push_back(i);

        print_message(1,"New middle set\n");
        print(1,this->edges[new_edge_index].middle_set);

        // Set the value of new_ms_size
        *new_ms_size = (int)this->edges[new_edge_index].middle_set.size();

        // Update max middle set size
        if((int)this->edges[new_edge_index].middle_set.size() > this->max_middle_set_size)
        {
            print_message(1,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n%s:  max size increasing from %d -> %d\n",__FUNCTION__,
                this->max_middle_set_size , this->edges[new_edge_index].middle_set.size());
            this->max_middle_set_size = this->edges[new_edge_index].middle_set.size();
        }

#if 0
        // Debugging
        int left_size=0, right_size=0;
        for(int i=0;i<this->G->get_num_nodes();i++)
        {
            if(A_union_vec[i])
                left_size++;
            if(B_union_vec[i])
                right_size++;
        }

        print_message(0,"Splitting %d\n\tleft_size=%d\n\tright_size=%d\n",u,left_size,right_size);
        print_message(0,"\tNew edge has middle set size %d\n\t",this->edges[new_edge_index].middle_set.size());
        print(0,this->edges[new_edge_index].middle_set);

#endif
        // Update internal node counters
        this->num_nodes++;
        this->num_interior_nodes++;

        if(this->num_interior_nodes == this->G->get_num_edges()-2)
            // The decomposition is now valid
            this->is_valid=true;

        // Just (re-)calculate the # of leaf nbrs of u
        this->calculate_num_leaf_nbrs(u);
        // Calculate the # of leaf nbrs of the new node
        this->calculate_num_leaf_nbrs(new_node_index);
        return new_node_index;
}

/**
* Takes the edges adjacent to node u and creates two new nodes, v and w
* that are both adjacent to u.  The list A is a list of neighboring LEAF nodes or edges
*leading to LEAF nodes that will now be adjacent to v.  The remaining LEAF nodes or edges
not in A will be adjacent to v.  It is assumed that u is incident to exactly 
one edge that does not lead to a leaf.  The value of type should be either
NODE_SPLIT or EDGE_SPLIT depending on what A represents.  Returns the index of the
new node added to the tree.  Sets the value of new_ms_size1 and new_ms_size2 to be
the size of the middle set of the newly added edges u-v and u-w.
*/
int BDTree::spawn_nodes(int u, list<int> *A, int type, int *new_ms_size1, int *new_ms_size2,
    int *n1, int *n2)

{
    // Make sure that u currently has 4 or more neighbors in the BDTree
    if(this->nodes[u].edges.size() < 4)
    {
        cerr<<this->nodes[u];
        fatal_error("%s:  trying to split node %d that has degree %d < 4!\n",__FUNCTION__,u,this->nodes[u].edges.size());
    }

    // After the split, u will be adjacent (in the BDTree) to its old non-leaf neighbor (call it t),
    // and will also be adjacent to two new nodes, v and w which are the two new nodes spawned by
    // the paritition. v will be adjacent to edges/nodes in A, and w will be adjacent to edges/nodes
    // not in A.
    //
    // Old:  (INTERIOR) t- u - A
    //                     |
    //                     B
    //
    // New:  (INTERIOR) t - u - v - A
    //                      |
    //                      w - B


    vector<int>::iterator vv;
    list<int>::iterator ii,jj;
    list<BDTreeEdge>::iterator ee;
    list<int> complement, B, new_u_edges;
    int nbr, new_v_index, new_w_index;

    vector<bool> A_vec(this->num_nodes,false);
    vector<bool> A_union_vec(this->G->get_num_nodes(),false);
    vector<bool> B_union_vec(this->G->get_num_nodes(),false);
    vector<bool> interior_union_vec(this->G->get_num_nodes(),false);

    list<int> AA;
    // if type is EDGE_SPLIT, then A is a list of edge indices.
    // just populate an alternate list with the indices of the adjacent
    // nodes - then we can use the same code.
    if(type==EDGE_SPLIT)
    {
        print_message(10,"Translating edges - %d nodes, %d edges in BDtree\n",this->num_nodes,this->num_edges);
        // Loop over all the edge indices in A and add the node index to a new list AA
        //print_message(0,"Middle sets of split edges\n");
        for(ii=A->begin();ii!=A->end();++ii)
        {
            print_message(1,"Translating edge %d (%d-%d) into nodes\n",*ii,this->edges[*ii].start,
                this->edges[*ii].end);
            if(this->edges[*ii].start==u)
                AA.push_back(this->edges[*ii].end);
            else
                if(this->edges[*ii].end==u)
                    AA.push_back(this->edges[*ii].start);
                else
                {
                    // This is an error - neither the edge start or end hits the split node u
                    cerr<<*this;
                    fatal_error("%s:  u=%d not found in edge %d (%d-%d)\n",__FUNCTION__,u,
                        *ii,this->edges[*ii].start,this->edges[*ii].end);           
                }
        }
    }
    // AA should now be a node list - now use the pointer C to point to the list we actually 
    // want to use (either the original list A, or the new list just created AA)
    list<int> *C;
    if(type==NODE_SPLIT)
        C=A;
    else
    {
        C=&AA;
        print_message(10,"%s:  Translated initial edge list into node list\n",__FUNCTION__);
        print(10,*A);
        print_message(10,"Nodes:\n");
        print(10,AA);
    }

    // Make sure that u is not in the list C and that all nodes in C are leaf nodes
    for(ii=C->begin();ii!=C->end();++ii)
    {
        // Create an incidence vector for the node list C
        A_vec[*ii]=true;
        if((*ii)==u)
        {
            print_message(0,"%s:  called where %d is in the list of tree nodes A!\n",__FUNCTION__,u,u);
            print(0,*C);
            exit(-1);
        }
        if(this->nodes[*ii].graph_edge_start==-1 || this->nodes[*ii].graph_edge_end==-1)
        {
            print_message(0,"%s:  list of tree nodes includes non-leaf node %d. Removing it from list!\n",__FUNCTION__,*ii);
            C->remove(*ii);            
        }
    }


    // Classify all the BDTreeEdges incident to u - all but one will be involved in the spawning
    list<int> v_edges, w_edges, v_nbrs, w_nbrs;
    int num_interior_nbrs=0, interior_edge, interior_nbr;
    for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
    {       
        // *ii is the index in the edges array of the edge incident to u
        // nbr is the non-u endpoint of the edge - check to see if it is in A by looking at A_vec
        nbr=-1;
        if(this->edges[*ii].start != u )
            nbr=this->edges[*ii].start;
        else
            nbr=this->edges[*ii].end;

        // Sanity check
        if(nbr==-1)
            fatal_error("%s:  Did not find node %d in an incident edge! (%d-%d)",
            __FUNCTION__,u,this->edges[*ii].start,this->edges[*ii].end);

        // Now make sure that the nbr is a leaf!
        if(this->nodes[nbr].graph_edge_start!=-1)
        {
            if(!A_vec[nbr])
            {
                for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                    B_union_vec[*jj]=true;
                w_nbrs.push_back(nbr);
                w_edges.push_back(*ii);
            }
            else
            {

                // This is an edge in the "A half" - these edges are adjacent to v
                for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                    A_union_vec[*jj]=true;
                v_nbrs.push_back(nbr);
                v_edges.push_back(*ii);                
            }
        }
        else
        {
            // nbr is an interior node - this means that *ii represents the edge that "points" to the interior
            // of the tree
            interior_nbr=nbr;
            num_interior_nbrs++;
            interior_edge=*ii;
        }
    }
    // Make sure we encountered only one interior neighbor
    if(num_interior_nbrs!=1)
        fatal_error("%s:  Encountered %d!=1 interior neighbors of node %d\n",num_interior_nbrs,u);

    // Now update the BDTree with the information - need to consider when w_edges or v_edges contains only a single edge
    this->nodes[u].edges.clear();
    this->nodes[u].edges.push_back(interior_edge);

    // Create vector for middle set of interior edge 
#if 1
    for(ii=this->edges[interior_edge].middle_set.begin(); 
        ii!=this->edges[interior_edge].middle_set.end(); ++ii)
        interior_union_vec[*ii]=true;
#else
    // This is slow - used this to check computation of middle set
    list<int> interior_ms;
    this->middle_set_union(&(this->nodes[interior_nbr].edges),&interior_ms);
    for(ii=interior_ms.begin(); ii!=interior_ms.end(); ++ii)
        interior_union_vec[*ii]=true;
#endif

    int interior_ms_size=this->edges[interior_edge].middle_set.size();

    if(v_edges.size() > 1)
    {
        // Add the new node v
        new_v_index=this->num_nodes;
        this->num_nodes++;
        this->num_interior_nodes++;
        this->nodes[new_v_index].id=new_v_index;

        this->nodes[new_v_index].graph_edge_start=-1;
        this->nodes[new_v_index].graph_edge_end=-1;
        jj=v_nbrs.begin();
        for(ii=v_edges.begin();ii!=v_edges.end();++ii)
        {
            this->nodes[new_v_index].edges.push_back(*ii);
            this->edges[*ii].start=*jj;
            this->edges[*ii].end=new_v_index;
            ++jj;
        }

        // Add u-v edge
        int new_uv_edge=this->num_edges;
        this->edges[new_uv_edge].start = u;
        this->edges[new_uv_edge].end = new_v_index;
        this->nodes[u].edges.push_back(new_uv_edge);
        this->nodes[new_v_index].edges.push_back(new_uv_edge);
        this->num_edges++;

        // Compute new middle set
        for(int i=0;i<this->G->get_num_nodes();i++)
        {
            if(A_union_vec[i] && (interior_union_vec[i] || B_union_vec[i]) )
                this->edges[new_uv_edge].middle_set.push_back(i);
        }
        print_message(0,"New uv middle set (size=%d, interior ms=%d)\n",(int)this->edges[new_uv_edge].middle_set.size(),
            interior_ms_size);
        print(0,this->edges[new_uv_edge].middle_set);

        *new_ms_size1 = (int)this->edges[new_uv_edge].middle_set.size();

    }
    else
    {
        // there is only one node/edge on this side of the split - don't need to add a new node
        ii=v_edges.begin();
        this->nodes[u].edges.push_back(*ii);
        *new_ms_size1 = (int)this->edges[*ii].middle_set.size();

        // set new_v_index=-1 so that we set n1 correctly upon returning!!
        new_v_index=-1;
    }

    if(w_edges.size()>1)
    {
        // Add the new node w
        new_w_index=this->num_nodes;
        this->num_nodes++;
        this->num_interior_nodes++;
        this->nodes[new_w_index].id=new_w_index;

        this->nodes[new_w_index].graph_edge_start=-1;
        this->nodes[new_w_index].graph_edge_end=-1;
        jj=w_nbrs.begin();
        for(ii=w_edges.begin();ii!=w_edges.end();++ii)
        {
            this->nodes[new_w_index].edges.push_back(*ii);
            this->edges[*ii].start=*jj;
            this->edges[*ii].end=new_w_index;
            ++jj;
        }

        // u-w edge
        int new_uw_edge=this->num_edges;
        this->edges[new_uw_edge].start = u;
        this->edges[new_uw_edge].end = new_w_index;
        this->nodes[u].edges.push_back(new_uw_edge);
        this->nodes[new_w_index].edges.push_back(new_uw_edge);
        this->num_edges++;

        // Compute new middle set
        for(int i=0;i<this->G->get_num_nodes();i++)
        {
            if(B_union_vec[i] && (interior_union_vec[i] || A_union_vec[i]) )
                this->edges[new_uw_edge].middle_set.push_back(i);
        }
        print_message(0,"New uw middle set (size=%d, interior ms=%d)\n",(int)this->edges[new_uw_edge].middle_set.size(),
            interior_ms_size);
        print(0,this->edges[new_uw_edge].middle_set);
        *new_ms_size2 = (int)this->edges[new_uw_edge].middle_set.size();
    }
    else
    {
        // there is only one node/edge on this side of the split - don't need to add a new node
        ii=w_edges.begin();
        this->nodes[u].edges.push_back(*ii);
        *new_ms_size2 = (int)this->edges[*ii].middle_set.size();

        // set new_v_index=-1 so that we set n1 correctly upon returning!!
        new_w_index=-1;
    }


    // Update max middle set size
    if(*new_ms_size1 > this->max_middle_set_size)
    {
        print_message(0,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n%s:  max size increasing from %d -> %d\n",__FUNCTION__,
            this->max_middle_set_size , *new_ms_size1);
        this->max_middle_set_size = *new_ms_size1;
    }

    if(*new_ms_size2 > this->max_middle_set_size)
    {
        print_message(0,"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n%s:  max size increasing from %d -> %d\n",__FUNCTION__,
            this->max_middle_set_size , *new_ms_size2);
        this->max_middle_set_size = *new_ms_size2;
    }

#if 0
    // Debugging
    int left_size=0, right_size=0;
    for(int i=0;i<this->G->get_num_nodes();i++)
    {
        if(A_union_vec[i])
            left_size++;
        if(B_union_vec[i])
            right_size++;
    }

    print_message(0,"Splitting %d\n\tleft_size=%d\n\tright_size=%d\n",u,left_size,right_size);
    print_message(0,"\tNew edge has middle set size %d\n\t",this->edges[new_edge_index].middle_set.size());
    print(0,this->edges[new_edge_index].middle_set);

#endif


    if(this->num_interior_nodes == this->G->get_num_edges()-2)
        // The decomposition is now valid
        this->is_valid=true;

    *n1=new_v_index;
    *n2=new_w_index;

    // Calculate # of leaf nbrs
    // New:  (INTERIOR) t - u - v - A
    //                      |
    //                      w - B

    this->nodes[u].num_leaf_nbrs=0;
    this->calculate_num_leaf_nbrs(new_v_index);
    this->calculate_num_leaf_nbrs(new_w_index);

    return 1;
}

/**
* Returns the size of the middle set if the split were to be made. Does not 
* actually modify BDTree.
*/
int BDTree::split_size(int u, list<int> *A, int type)
{
    // Make sure that u currently has 4 or more neighbors in the BDTree
    if(this->nodes[u].edges.size()<=3)
    {
        cerr<<this->nodes[u];
        fatal_error("%s:  trying to split node %d that has degree %d<=3!\n",__FUNCTION__,u,this->nodes[u].edges.size());
    }

    // After the split, u will be adjacent (in the BDTree) to the nodes whose indices are contained in the partition, and
    // the new node will be adjacent to all of the other (old) neighbors of u
    // Old:  A - u - B
    // New:  A - u - v -B
    // where A and b are sets of nodes and v is the new BDTreeNode

    vector<int>::iterator vv;
    list<int>::iterator ii,jj;
    list<BDTreeEdge>::iterator ee;
    int nbr;

    vector<bool> A_vec(this->num_nodes,false);
    vector<bool> A_union_vec(this->G->get_num_nodes(),false);
    vector<bool> B_union_vec(this->G->get_num_nodes(),false);

    list<int> AA;
    // if type is EDGE_SPLIT, then A is a list of edge indices.
    // just populate an alternate list with the indices of the adjacent
    // nodes - then we can use the same code.
    if(type==EDGE_SPLIT)
    {
        // Loop over all the edge indices in A and add the node index to a new list AA
        //print_message(0,"Middle sets of split edges\n");
        for(ii=A->begin();ii!=A->end();++ii)
        {
            //print_message(0,"\t");
            print(1,this->edges[*ii].middle_set);

            print_message(10,"Translating edge %d (%d-%d) into nodes\n",*ii,this->edges[*ii].start,
                this->edges[*ii].end);
            if(this->edges[*ii].start==u)
                AA.push_back(this->edges[*ii].end);
            else
                if(this->edges[*ii].end==u)
                    AA.push_back(this->edges[*ii].start);
                else
                {
                    // This is an error - neither the edge start or end hits the split node u
                    cerr<<*this;
                    fatal_error("%s:  u=%d not found in edge %d (%d-%d)\n",__FUNCTION__,u,
                        *ii,this->edges[*ii].start,this->edges[*ii].end);           
                }
        }
    }
    // AA should now be a node list - now use the pointer C to point to the list we actually 
    // want to use (either the original list A, or the new list just created AA)
    list<int> *C;
    if(type==NODE_SPLIT)
        C=A;
    else
        C=&AA;

    //print_message(0,"Splitting %d: sizes are %d,%d\n",u,C->size(),this->nodes[u].edges.size()-C->size());

    // Make sure that u is not in the list C!
    for(ii=C->begin();ii!=C->end();++ii)
        if((*ii)==u)
            fatal_error("%s:  Split(u,A) called where u is in the list of tree nodes A!\n",__FUNCTION__);

    // Create an incidence vector for the node list C
    for(ii=C->begin();ii!=C->end();++ii)
        A_vec[*ii]=true;


    // Loop over all BDTreeEdges incident to u
    for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
    {       
        // *ii is the index in the edges array of the edge incident to u
        // nbr is the non-u endpoint of the edge - check to see if it is in A by looking at A_vec
        nbr=-1;
        if(this->edges[*ii].start != u )
            nbr=this->edges[*ii].start;
        else
            nbr=this->edges[*ii].end;

        // Sanity check
        if(nbr==-1)
            fatal_error("%s:  Did not find node %d in an incident edge! (%d-%d)",
            __FUNCTION__,u,this->edges[*ii].start,this->edges[*ii].end);

        // The middle set of these updated edge should stay the same
        if(!A_vec[nbr])
        {
            // This is an edge in the "B half" - these edges need to be updated as they are no longer
            // incident to u
            for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                B_union_vec[*jj]=true;       

        }
        else
        {
            // This is an edge in the "A half" - these edges should remain unchanged
            for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                if(!A_union_vec[*jj])
                    A_union_vec[*jj]=true;
        }
    }

    // Compute the middle set of the new edge
    list<int> new_middle_set;
    for(int i=0;i<this->G->get_num_nodes();i++)
        if(A_union_vec[i] && B_union_vec[i])
            new_middle_set.push_back(i);

    return new_middle_set.size();
}


/**
* A greedy algorithm that finds the edge sets X and Y such that
* all edges are adjacent to u, and such that the middle set of the
* new edge introduced after splitting u in this manner has minimal
* cardinality.  Returns the cardinality of the middle set of the edge
* to be introduced if this greedy split were made.
*/
int BDTree::find_greedy_split(int u, list<int> *X)
{
    // First make sure u has degree at least 4
    if(this->nodes[u].edges.size()<4)
        fatal_error("%s:  Node u=%d has degree %d.  Can only split nodes with degree 4 or more\n",
        __FUNCTION__,u,this->nodes[u].edges.size());

    // If u has degree d, then we will have to exhaust over 2^d possible lists X, so make sure
    // this is a reasonable attempt
    if(this->nodes[u].edges.size()>26)
        fatal_error("%s:  Node u=%d has degree %d - this is too much work!\n",__FUNCTION__,u,
        this->nodes[u].edges.size());

    int i=0,j,deg=this->nodes[u].edges.size();

    vector<int> edge_vec(deg);
    list<int>::iterator ii;

    for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
    {
        // Sanity check
        if(this->edges[*ii].start!=u && this->edges[*ii].end!=u)
        {
            fatal_error("%s: when creating edge_vec, edge %d does not have u=%d as an endpoint!\n"
                "Endpoints are %d-%d\n",__FUNCTION__,*ii,u,this->edges[*ii].start,this->edges[*ii].end);

        }
        edge_vec[i]=*ii;
        i++;
    }

    vector<int> X_vec(this->G->get_num_nodes(),0);
    vector<int> Y_vec(this->G->get_num_nodes(),0);

    int order, best_order=1<<30, best_i=-1;

    for(i=1;i<(1<<deg);i++)
    {
        X_vec.clear();
        Y_vec.clear();
        int num_bits=0;
        for(j=0;j<deg;j++)
        {
            if(i & (1<<j))
            {
                // bit j is set in i, so the jth edge incident to u is in X
                // we need to include it's middle set in X_vec
                for(ii=this->edges[edge_vec[j]].middle_set.begin();ii!=this->edges[edge_vec[j]].middle_set.end();++ii)
                    X_vec[*ii]=1;
                num_bits++;
            }
            else
            {
                // bit j is not set in i, so the jth edge incident to u is in Y
                // we need to include it's middle set in Y_vec
                for(ii=this->edges[edge_vec[j]].middle_set.begin();ii!=this->edges[edge_vec[j]].middle_set.end();++ii)
                    Y_vec[*ii]=1;
            }
        }

        if(num_bits>=2)
        {
            // Now use X_vec to Y_vec to compute the order of the separation
            order=0;
            for(j=0;j<this->G->get_num_nodes();j++)
            {
                if(X_vec[j]==1 && Y_vec[j]==1)
                    order++;
            }
            if(order<best_order)
            {
                best_order=order;
                best_i=i;
            }
        }
    }
    print_message(0,"best_i=%d\n",best_i);

    // We can now reconstruct the EDGE set X - but to split, we actually want the
    // BDTreeNode indices.  This is what we put in the output list X.
    X->clear();
    for(j=0;j<deg;j++)
    {
        if(best_i & (1<<j))
        {
            // bit j is set in i, so the jth edge incident to u is in X
            // Get the other endpoint of edge[edge_vec[j]]
            // sanity check
            if(this->edges[edge_vec[j]].start!=u && this->edges[edge_vec[j]].end!=u)
                fatal_error("%s: edge_vec[j=%d] does not have u=%d as an endpoint!\n"
                "Endpoints are %d-%d\n",__FUNCTION__,j,u,this->edges[edge_vec[j]].start,this->edges[edge_vec[j]].end);

            print_message(0,"Adding edge_vec[%d] endpoints to X\n",j);
            cerr<<this->edges[edge_vec[j]];
            // Grab the non-u part of the edge
            if(this->edges[edge_vec[j]].start==u)
                X->push_back(this->edges[edge_vec[j]].end);
            else
                X->push_back(this->edges[edge_vec[j]].start);


            //X->push_back(edge_vec[j]);
        }
    }
    X->sort();
    return best_order;
}


//find & run all bfpushes arising if you do a push at u.
//Modified Nov 2 to be more efficient at checking nodes which may have a push.
//if newe_only is set to true, we only consider pushes involving the newly created edge at that vertex, and 
// u_edge should be set to the desired edge for the initial vertex u.
//Returns the number of splits (pushes) enacted.
int BDTree::doall_bfpushes(int u, bool newe_only, int u_edge)
{
    int num_pushes = 0;
    int currv;
    int curre;
    int newv, newe;
    list<int> queue;
    list<int> Equeue;
    bool found;
    queue.push_back(u); 
    if(newe_only)
	Equeue.push_back(u_edge);
    
    while(queue.size() > 0)
	{
	    //push vertex off end of queue
	    currv = queue.back(); 
	    queue.pop_back();
	    found = false;
	    if(newe_only)
		{
		    curre = Equeue.back(); 
		    Equeue.pop_back();
		}
	    if(this->nodes[currv].edges.size() > 3)
		{
		    if(newe_only)			    
			found = this->bf_push(currv, curre, &newv, &newe);
		    else
			found = this->bf_push(currv, &newv, &newe);
		    
		    if(found)
			{
			    num_pushes++;
			    if(this->nodes[currv].edges.size() >3)
				{
				    queue.push_back(currv); 
				    if(newe_only)
					Equeue.push_back(newe);
				}
			    if(this->nodes[newv].edges.size() >3)
				{
				    queue.push_back(newv);
				    if(newe_only)
					Equeue.push_back(newe);
				}
			}
		}
	}
    return num_pushes;
}//end of if doall_bfpushes



//Oct 29 - BDS
// brute force version of pushing - this just looks for a single separation
// at u which satisfies the pushing inequality and makes that split. 
// User is responsible for looping over tree repeatedly looking for more.
// Let Ti be a partial branch-decomposition of G and
// let u be a node of Ti with degree at least four. If D is the set of edges of
// Ti incident with u,  for each e ∈ D, let Me be the
// middle of the separation of G arising from e.
// Then if e1,e2 ∈ D are distinct edges and 
// Me1 ∪ Me2 ∩ (\cup Me : e ∈ D\{e1,e2}) ≤ max Me1 Me2
// Then we split with X = e1,e2  and Y =D\X.
// Returns true if a push was made, returns false otherwise.
bool BDTree::bf_push(int u, int *newv, int *newe)
{
    int e1e2size, i;
    list<int>::iterator ii, jj, kk, ll;
    vector<bool> Aunionvec(this->G->get_num_nodes(), false);
    vector<bool> e1vec(this->G->get_num_nodes(), false);
    vector<bool> e2vec(this->G->get_num_nodes(), false);
    *newv = -1; 
    *newe = -1;
    
    //consider each pair of edges occurring at u, looking for one that allows a safe split
    for(ii=this->nodes[u].edges.begin(); ii!=this->nodes[u].edges.end();++ii)
	{
	    //reset e1vec to all false
	    fill(e1vec.begin(), e1vec.end(), false);
	    //set the vector for e1 = *ii.
	    for(jj = this->edges[*ii].middle_set.begin(); jj!= this->edges[*ii].middle_set.end(); ++jj)
		{
		    e1vec[*jj] = true;
		}
	    jj=ii;
	    ++jj;
	    while(jj != this->nodes[u].edges.end())
		{
		    print_message(1, "Considering e1 = %d, e2 = %d at vertex %d\n", *ii, *jj, u);
		    //reset e2vec to all false
		    fill(e2vec.begin(), e2vec.end(), false);
		    //find e2's middle set vector
		    for(kk = this->edges[*jj].middle_set.begin(); kk!= this->edges[*jj].middle_set.end(); ++kk)
			{
			    e2vec[*kk] = true;
			}   
		    // now we need to know elements in e1 OR e2 and also in 
		    // the union of all the edges other than e1, e2.
		    //compute the union:
		    //reset Aunionvec to all false
		    fill(Aunionvec.begin(), Aunionvec.end(), false);
		    for(kk=this->nodes[u].edges.begin();kk!=this->nodes[u].edges.end();++kk)
			{
			    //not e1 or e2
			    if(*kk != *ii && *kk != *jj)
				{
				    for(ll = this->edges[*kk].middle_set.begin(); ll!= this->edges[*kk].middle_set.end(); ++ll)
					{
					    Aunionvec[*ll] = true;
					}
				} //if kk is not e1, e2 
				}//end of loop over edges at u
		    //and now we can find the size of the intersection
		    e1e2size = 0;
		    for(i = 0; i < G->get_num_nodes(); i++)
			{
			    if(Aunionvec[i] && (e1vec[i] || e2vec[i]))
				e1e2size++;
			}//loop over vertices of G
		    // now we check the inequality
		    if((e1e2size < (int) this->edges[*ii].middle_set.size()) || (e1e2size < (int) this->edges[*jj].middle_set.size()))
			{
			    list<int> A; 
			    A.push_back(*ii); 
			    A.push_back(*jj); 
			    //we need to split and return out of function.
			    //set the values of newv and newe as needed.
			    *newv = this->split_node(u, &A, EDGE_SPLIT, &i, newe);
			    print_message(0,"BFpush split node %d with edges %d and %d and got a new middle set of size %d\n", u, *ii, *jj, i);
			    return true;
			}//end of if splitting
			++jj;
			}//end of for jj loop
		}  //end of for ii loop
	    return false;
	}//end of function


//Nov 1 - BDS
// specialized version of pushing only looking for separations where
// one of the two 'chosen' edges was just introduced by a split. 
// this edge is specified by its index in the edge array, e1. This is otherwise
// a brute force version of pushing (it looks at all possible other choices for e2)
// and makes the split if one is found.
// User is responsible for looping over tree repeatedly looking for more.
// Returns true if a push was made, returns false otherwise.
bool BDTree::bf_push(int u, int e1, int *newv, int *newe)
{
    int e1e2size, i;
    list<int>::iterator ii, jj, kk, ll;
    vector<bool> Aunionvec(this->G->get_num_nodes(), false);
    vector<bool> e1vec(this->G->get_num_nodes(), false);
    vector<bool> e2vec(this->G->get_num_nodes(), false);
    *newv = -1; 
    *newe = -1;
    
    //set the vector for e1 
    for(jj = this->edges[e1].middle_set.begin(); jj!= this->edges[e1].middle_set.end(); ++jj)
	{
	    e1vec[*jj] = true;
	}
    
    //consider each pair of edges (e1, e2) occurring at u, looking for one that allows a safe split
    for(jj=this->nodes[u].edges.begin(); jj!=this->nodes[u].edges.end();++jj)
	{
	    if(*jj != e1)
		{
		    print_message(1, "Considering e1 = %d, e2 = %d at vertex %d\n", e1, *jj, u);
		    //reset e2vec to all false
		    fill(e2vec.begin(), e2vec.end(), false);
		    //find e2's middle set vector
		    for(kk = this->edges[*jj].middle_set.begin(); kk!= this->edges[*jj].middle_set.end(); ++kk)
			{
			    e2vec[*kk] = true;
			}   
		    // now we need to know elements in e1 OR e2 and also in 
		    // the union of all the edges other than e1, e2.
		    //compute the union:
		    //reset Aunionvec to all false
		    fill(Aunionvec.begin(), Aunionvec.end(), false);
		    for(kk=this->nodes[u].edges.begin();kk!=this->nodes[u].edges.end();++kk)
			{
			    //not e1 or e2
			    if(*kk != *ii && *kk != *jj)
				{
				    for(ll = this->edges[*kk].middle_set.begin(); ll!= this->edges[*kk].middle_set.end(); ++ll)
					{
					    Aunionvec[*ll] = true;
					}
				} //if kk is not e1, e2 
			}//end of loop over edges at u
		    //and now we can find the size of the intersection
		    e1e2size = 0;
		    for(i = 0; i < G->get_num_nodes(); i++)
			{
			    if(Aunionvec[i] && (e1vec[i] || e2vec[i]))
				e1e2size++;
			}//loop over vertices of G
		    // now we check the inequality
		    if((e1e2size < (int) this->edges[*ii].middle_set.size()) || (e1e2size < (int) this->edges[*jj].middle_set.size()))
			{
			    list<int> A; 
			    A.push_back(*ii); 
			    A.push_back(*jj); 
			    //we need to split and return out of function.
			    //set newv and newe as appropriate.
			    *newv = this->split_node(u, &A, EDGE_SPLIT, &i, newe);
			    print_message(0,"BFpush split node %d with edges %d and %d and got a new middle set of size %d\n", u, *ii, *jj, i);
			    return true;
			}//end of if splitting
		    ++jj;
		}//end if
	}//end of for jj loop
    return false;
}//end of function
    
    
int BDTree::push(int u)
{
    int i,j;
    vector<bool> middle_set_vec(this->G->get_num_nodes()+1,false);
    list<int>::iterator ii,jj,kk,ll;

    int num_in_union=0;
    vector<int> positions(this->G->get_num_nodes()+1,-1);
    vector<int> inv_positions(this->G->get_num_nodes()+1,-1);
    vector<int> edge_perm(this->nodes[u].edges.size());
    i=0;
    for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
    {
        edge_perm[i]=*ii;
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!middle_set_vec[*jj])
                num_in_union++;
            middle_set_vec[*jj]=true;
        }
        i++;
    }
    // middle_set_vec[j]=true for every j that is in the middle set of some edge adjacent to u
    // positions[i] is the position of i in a 0-based list of the true entries in middle_set_vec
    j=0;
    for(i=0;i<this->G->get_num_nodes()+1;i++)
    {
        if(middle_set_vec[i])
        {
            positions[i]=j;
            inv_positions[j]=i;
            j++;
        }
    }

    vector<list<int> > cols(num_in_union);
    vector<list<int> > rows(this->nodes[u].edges.size());

    // Create a sparse matrix with lists
    // A 1 in position i,j indicates that inv_positions[j] is in the middle set of edge_perm[i]
    // We want to find columns in this matrix that sum to 2. The edges corresponding to these rows
    // are the e1,e2 pairs for the push.
    i=0;
    for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
    {
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            // The row corresponding to the ith edge contains *jj
            cols[positions[*jj]].push_back(i);
            rows[i].push_back(positions[*jj]);
        }
        rows[i].sort();
        i++;
    }

    // Sort the column lists
    for(i=0;i<num_in_union;i++)
        cols[i].sort();

    list<int_int> pushes;
    int_int push;
    list<int> eunion;
    list<int> middle_set;
    list<list<int> > new_middle_sets;
    int num_pushes=0;
    vector<bool> used_cols(num_in_union,false);
    vector<bool> pushed(this->num_nodes+1,false);
    while(true)
    {
        // Look for a column with exactly two entries
        for(i=0;i<num_in_union;i++)
            if(cols[i].size()==2 && !used_cols[i] && !pushed[cols[i].front()] && !pushed[cols[i].back()])
            {
                used_cols[i]=true;
                break;
            }

            if(i==num_in_union)
                // No columns found with 2 entries ==> no more pushes possible
                break;

            // We found a push in column i
            push.p1=cols[i].front();
            push.p2=cols[i].back();
            pushed[push.p1]=true; 
            pushed[push.p2]=true;

            // Now compute the middle set of this push
            // eunion is the union of the middle sets of e1 and e2
            eunion.clear();
            for(ii=rows[push.p1].begin();ii!=rows[push.p1].end();++ii)
                eunion.push_back(*ii);
            for(ii=rows[push.p2].begin();ii!=rows[push.p2].end();++ii)
                eunion.push_back(*ii);
            eunion.sort(); 
            eunion.unique();

            // Now we want to find the intersection of eunion with the union of the remaining middle sets
            // eunion is a list of column indices, so go through each column referenced in eunion and
            // look for something other than e1 or e2 - if something is found, then this column is added to
            // the middle set
            middle_set.clear();
            for(ii=eunion.begin(); ii!=eunion.end(); ++ii)
            {
                for(jj=cols[*ii].begin();jj!=cols[*ii].end();++jj)
                {
                    if(*jj!=push.p1 && *jj!=push.p2)
                    {
                        // This column (middle set entry) belongs to the middle set of p1 or p2 AND something else,
                        // so it is in the intersection - add it to the middle set
                        middle_set.push_back(*ii);
                        break;
                    }
                }
                // breaks to here
            }
            // Now remember the two edges being pushed, and the middle set of the new edge
            pushes.push_back(push);
            new_middle_sets.push_back(middle_set);
            print_message(1,"u=%d: found push (%d,%d)\n",u,push.p1,push.p2);
            print(1,middle_set);

            // Now clear out the rows for push.p1 and push.p2 - then overwrite one of the rows with the new middle set!!
            rows[push.p1].clear();
            rows[push.p2].clear();
            rows[push.p1]=middle_set;
    }

    print_message(1,"\n\n\n%d pushes found. Starting to split\n\n\n",num_pushes);

    // Now modify the tree with the pushes
    list<int> push_list;
    int splitting_node=u;
    list<int> new_middle_set;
    list<list<int> >::iterator mm=new_middle_sets.begin();
    for(list<int_int>::iterator tt=pushes.begin();tt!=pushes.end();++tt)
    {
        print_message(1,"  tt  is (%d,%d)\n",tt->p1,tt->p2);
        print_message(1,"e[tt] is (%d,%d)\n",edge_perm[tt->p1],edge_perm[tt->p2]);
        push_list.push_back(edge_perm[(*tt).p1]);
        push_list.push_back(edge_perm[(*tt).p2]);
        print_message(1,"splitting_node=%d; edges=(%d,%d)\n",splitting_node,edge_perm[(*tt).p1],edge_perm[(*tt).p2]);
        print_message(1,"\tclaimed middle set: ");
        new_middle_set.clear();
        for(ii=mm->begin();ii!=mm->end();++ii)
        {
            new_middle_set.push_back(inv_positions[*ii]);
            print_message(1,"%d ",inv_positions[*ii]);
        }
        print_message(1,"\n");


        if(!this->verify_push(u,edge_perm[(*tt).p1],edge_perm[(*tt).p2]))
            // The push satisfies the test inequality but is not actually valid!!
            print_message(1,"%s:  Invalid push %d, %d, %d\n",__FUNCTION__,u,edge_perm[(*tt).p1],edge_perm[(*tt).p2]);
        else
        {
            print_message(1,"%s:  valid push. Old middle set sizes: %d, %d; New middle set size: %d\n",
                __FUNCTION__,this->edges[edge_perm[(*tt).p1]].middle_set.size(),  
                this->edges[edge_perm[(*tt).p2]].middle_set.size(),
                new_middle_set.size());
            // Manually split the node
            int new_node_index=this->num_nodes;
            int new_edge_index=this->num_edges;
            // increment counters
            this->num_edges++;
            this->num_interior_nodes++;
            this->num_nodes++;
            // The new edge joins splitting_node - new_node_index
            this->edges[new_edge_index].start=splitting_node;
            this->edges[new_edge_index].end=new_node_index;
            this->edges[new_edge_index].middle_set=new_middle_set;
            print_message(1,"\t\t\t\tAfter push  new middle set %d\n",
                new_middle_set.size());
            this->nodes[new_node_index].edges.push_back(edge_perm[(*tt).p1]);
            this->nodes[new_node_index].edges.push_back(edge_perm[(*tt).p2]);
            this->nodes[new_node_index].edges.push_back(new_edge_index);
            this->nodes[new_node_index].id=new_node_index;
            this->nodes[new_node_index].edges.sort();
            // CSG, Aug. 22 - this new node is interior
            this->nodes[new_node_index].graph_edge_start=-1;
            this->nodes[new_node_index].graph_edge_end=-1;

            // Now remove the pushed edges from splitting_node
            this->nodes[splitting_node].edges.remove(edge_perm[(*tt).p1]);
            this->nodes[splitting_node].edges.remove(edge_perm[(*tt).p2]);
            // The splitting_node is now adjacent to the newly created edge
            this->nodes[splitting_node].edges.push_back(new_edge_index);
            print_message(1,"Splitting node %d now has %d nbrs\n",splitting_node,this->nodes[splitting_node].edges.size());
            // Update the endpoints of these two edges!!
            // edge 1
            if(this->edges[edge_perm[(*tt).p1]].end==splitting_node)
                this->edges[edge_perm[(*tt).p1]].end=new_node_index;
            else
                this->edges[edge_perm[(*tt).p1]].start=new_node_index;
            // edge 2
            if(this->edges[edge_perm[(*tt).p2]].end==splitting_node)
                this->edges[edge_perm[(*tt).p2]].end=new_node_index;
            else
                this->edges[edge_perm[(*tt).p2]].start=new_node_index;

            // Update max middle set size
            if((int)this->edges[new_edge_index].middle_set.size() > this->max_middle_set_size)
            {
                print_message(0,"%s:  max size increasing from %d -> %d\n", __FUNCTION__,
                    this->max_middle_set_size , this->edges[new_edge_index].middle_set.size());
                this->max_middle_set_size = this->edges[new_edge_index].middle_set.size();
            }

            // update # leaf nbrs for splitting_node and new_node_index
            this->calculate_num_leaf_nbrs(splitting_node);
            this->calculate_num_leaf_nbrs(new_node_index);         


            // We performed a push;
            num_pushes++;

            print_message(1,"\t\t%d interior nodes (need %d for validity)\n",
                this->num_interior_nodes,this->G->get_num_edges()-2);

            // Check to see if the decomposition is now valid
            if(this->num_interior_nodes == this->G->get_num_edges()-2)
            {
                // The decomposition is now valid
                this->is_valid=true;
                // Stop pushing!!!
                break;
            }

            // Check the degree of the node we are splitting
            if(this->nodes[splitting_node].edges.size()<=3)
            {
                print_message(1,"Splitting node has degree %d!\n",this->nodes[splitting_node].edges.size());
                break;
            }
        }
        // Advance to the next middle set
        ++mm;
    }

    // update u's leaf count
    this->calculate_num_leaf_nbrs(u);

    print_message(1,"%d pushes found\n",num_pushes);
    return num_pushes;   
}

bool BDTree::verify_push(int u, int e1, int e2)
{
    if(this->nodes[u].edges.size()<=3)
    {
        print_message(0,"Trying to verify push for node %d which has degree %d<=3!\n",
            u,this->nodes[u].edges.size());
        return false;
    }
    // First try to find a w in M_e1 \cap M_e2 that is not in any other edges incident to u
    vector<int> vec1(this->G->get_num_nodes()+1,0);
    vector<int> vec2(this->G->get_num_nodes()+1,0);
    int i;

    list<int>::iterator ii,jj;

    for(ii=this->edges[e1].middle_set.begin();ii!=this->edges[e1].middle_set.end();++ii)
        vec1[*ii]++;
    for(ii=this->edges[e2].middle_set.begin();ii!=this->edges[e2].middle_set.end();++ii)
        vec1[*ii]++;
    // vec1[i] = 2 for all i\in (M_e1 \cap M_e2)

    // Now fill up vec2
    for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
    {
        if(*ii != e1 && *ii !=e2)
            for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                vec2[*jj]++;
    }

    // Now try to find w where vec1[w]=2 and vec2[w]=0
    bool found_w=false;
    for(i=0;i<=this->G->get_num_nodes();i++)
    {
        if(vec1[i] == 2 && vec2[i]==0)
        {
            print_message(1,"Found w=%d for edges %d,%d\n",i,e1,e2);
            found_w=true;
        }
    }
    if(!found_w)
    {
        // This should possibly be a fatal error?
        print_message(0,"%s:  Didn't find w!\n",__FUNCTION__);
        print(0,vec1);
        print(0,vec2);
        return false;
    }

    // Now check the inequality
    int max_mid_e1_e2, mid_others=0;
    for(i=0;i<=this->G->get_num_nodes();i++)
        if(vec1[i]>0 && vec2[i]>0) // i is in the middle set of (e1 or e2) AND some other edge
            mid_others++;
    if(this->edges[e1].middle_set.size()>=this->edges[e2].middle_set.size())
        max_mid_e1_e2=this->edges[e1].middle_set.size();
    else
        max_mid_e1_e2=this->edges[e2].middle_set.size();
    print_message(1,"%s:  new middle_set size computed to be %d; max(m(e1),m(e2))=%d\n",
        __FUNCTION__,mid_others,max_mid_e1_e2);
    if(mid_others <= max_mid_e1_e2)
        return true;
    else
    {
        print_message(1,"pushing inequality violated %d not <= %d!\n",mid_others,
            max_mid_e1_e2);
        return false;
    }


    return true;

}

/**
* Creates the middle set graph for v. 
* Optionally populates the CRS data structures xadj and adjncy of the Graph.  
* Note that the adjncy vector (and all information about the returned Graph)
* will contain entries between 0 and num_in_union-1 where 
* num_in_union is the total # of vertices from the original graph contained 
* in the middle sets of all edges incident to v.  To get back to the indexing 
* in the original graph and vice versa, the 
* perm array maps from {0,1,...,num_in_union-1} to {0,1,...|V|-1} and iperm
* goes the other direction.
* Note: the neighbor lists will be sorted, and the degree array will be accurate!
*/
Graph *BDTree::construct_middle_set_graph(int v, vector<int> *perm, vector<int> *iperm, bool fill_CRS)
{
    int i,j;
    vector<bool> middle_set_vec(this->G->get_num_nodes()+1,false);
    list<int>::iterator ii,jj,kk,ll;

    int num_in_union=0, max_label=-1;
    list<int> labels;
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!middle_set_vec[*jj])
            {
                num_in_union++;
                labels.push_back(*jj);
                if(*jj>max_label)
                    max_label=*jj;
            }
            middle_set_vec[*jj]=true;
        }
    }

    // Sort the labels
    labels.sort();
    // Create a permutation vector to go from labels->index
    iperm->resize(max_label+1,-1);
    perm->resize(labels.size()+1,-1);

    j=0;
    for(ii=labels.begin();ii!=labels.end();++ii)
    {
        iperm->at(*ii)=j;
        perm->at(j)=*ii;
        j++;
    }

    // H will have num_in_union nodes
    Graph *H=new Graph(num_in_union);


    // Set the labels
    for(i=0;i<num_in_union;i++)
        H->nodes[i].label=perm->at(i);

    // Now go through the middle sets again to add the edges
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        // Assumes middle sets are sorted! There will be duplicate edges at the end of this
        // which will be removed.
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            kk=jj;
            ++kk;
            for( ; kk!=this->edges[*ii].middle_set.end(); ++kk)
            {
                // *jj and *kk are both in the middle set of edge *ii which is adjacent to v
                // Add an edge in H between *jj and *kk

                H->nodes[iperm->at(*jj)].nbrs.push_back(iperm->at(*kk));
                H->nodes[iperm->at(*kk)].nbrs.push_back(iperm->at(*jj));
            }          
        }
    }

    // We will have duplicate edges in this graph
    int num_edges=0;
    for(i=0;i<num_in_union;i++)
    {
        H->nodes[i].nbrs.sort();
        H->nodes[i].nbrs.unique();
        num_edges += H->nodes[i].nbrs.size();
    }

    //BDS - number of edges was being double counted (Aug 10)
    num_edges = num_edges/2;

    // public/private issue
    H->set_num_edges(num_edges);
    H->set_num_nodes(num_in_union);
    H->set_simple(1);
    H->canonical=true;

    //set the degrees - do we need this?
    H->recompute_degrees();

    //H->set_symmetric();

    // Now populate CRS if needed
    if(fill_CRS) 
        H->populate_CRS();

    return H;
}

/**
* Helper function to check for a two separation in the branch decomposition
* at a specified vertex v. If one is found, X and Y are populated 
* with edge indices so that |M(X,Y)| = 2.
* returns true if one is found, false otherwise.
*/
bool BDTree::two_separation(int v, list<int> *X, list<int> *Y, goblinController *CT)
{
    list<int>::iterator ii, jj;
    bool ts_debug = false;
    vector<int> edge_index;
    vector<int> HtoG; 
    vector<int> GtoH;

    //We first build an auxiliary graph H, which is symmetric with sorted adj lists
    //and a populated degree array.
    //We fill in the CRS representation and use  HtoG, and GtoH to 
    //encode the information about how the middle sets map into the graph H. 
    //We assume edges of H will be represented as (i,j) with 
    // i < j when mapped into the aux. graph F.
    Graph *H = this->construct_middle_set_graph(v,&HtoG, &GtoH, true);

    if(ts_debug)
        cerr << *H;

    int u = -1, i = 0, j = 0;
    int  k;
    //we'll pick a vertex in H with degree 3 or more.  
    while(u == -1 && i < H->get_num_nodes())
    {
        if(H->get_degree(i) >= 3) // need to make sure equality is okay here
            u = i;
        i++;
    }
    if(u == -1)
    {
        print_message(0, "%s: Warning: The middle set graph at %d had no vertices of degree 3. Proceeding as if there is no 2-separation at this node.\n", __FUNCTION__, v);
        return false;
    }

    if(ts_debug)
        print_message(0, "Selected u = %d\n", u);

    //we'll use the first three neighbors to form our edges e,f,g.
    int e = H->xadj[u]; 
    int f = e+1; 
    int g = f+1;

    print_message(0, "e = (%d, %d), f = (%d, %d), g = (%d, %d)\n", HtoG[H->adjncy[e]], HtoG[u], HtoG[H->adjncy[f]], HtoG[u], HtoG[H->adjncy[g]], HtoG[u]);

    if((u != H->get_num_nodes()-1 && g >= H->xadj[u+1]) || g > (int) H->adjncy.size()-1)
        fatal_error("%s: We chose a vertex in H which should have had degree at least 3, but there were not 3 neighbors in the CRS adjcny array!\n", __FUNCTION__);

    //We'll now create a goblin graph on which to run max flow. 

    int F_num_vertices = 2 + 2*H->get_num_nodes();
    int source = F_num_vertices-1;
    int sink = F_num_vertices-2;
    TNode n = F_num_vertices;

    //we'll construct the goblin graph directly. 
    sparseDiGraph *F; 
    //I wanted to use a single graph object here and remove the 
    //few edges that change between iterations, but I can't find a function for
    //removing an arc in goblin (cancelArc and deleteArcs() seem promising in 
    //SparseRepresentation, but that's a protected part of sparseDigraph!
    //note hte only edges that differ between the different f's are the ones 
    //from the source & to the sink.

    int a1, a2, b1; int testnum = 1;
    //partitioning with f,g forming source in A and e the sink in B 
    sparseDiGraph Ffg_e(n, *CT);
    //partitioning with e,g forming source in A and f the sink in B 
    sparseDiGraph Feg_f(n, *CT);

    bool found = false; 
    TCap flowValue;
    TNode *dist;    

    // We will need two nodes for every non source/sink node
    // Every such node i will have a copy with index (i + shift)
    int shift=H->get_num_nodes();

    //we will first try two max flow problems looking for an appropriate 
    //partition which splits e and f. 
    while(found != true)
    {
        if(testnum == 1)
        {
            F = &Ffg_e; 
            a1 = f; 
            a2 = g; 
            b1= e; 
        }
        else if(testnum == 2)
        {
            F = &Feg_f; 
            a1 = e; 
            a2 = g; 
            b1 = f; 
        }
        else
        {
            //we need to look for a partition with e and f on the same side.
            break;
        }

        //insert the dummy arcs to do node-disjoint paths.
        for(i = 0; i < H->get_num_nodes(); i++)
        {
            F->InsertArc(i, i+shift);
        }

        //insert the arcs to/from the source & sink.
        F->InsertArc(source, u);
        F->InsertArc(source, H->adjncy[a1]);
        F->InsertArc(source, H->adjncy[a2]);
        F->InsertArc(u+shift, sink);
        F->InsertArc(H->adjncy[b1]+shift, sink);

        if(ts_debug)
        {
            print_message(1, "Adding edges from  %d (source) to %d, %d, %d\n", source, u, H->adjncy[a1], H->adjncy[a2]);
            print_message(1, "Adding edges from  %d,%d to %d (sink)\n", u, H->adjncy[b1], sink);
        }

        //we now need to add the edges of H.
        // The neighbors of node i are adjncy[xadj[i]],adjncy[xadj[i]+1],...adjncy[xadj[i+1]-1]	
        //only add edges (i,j) with j > i, so we don't have to look at i = H->get_num_nodes()-1.
        for(i = 0; i < H->get_num_nodes()-1; i++)
        {
            for(j = H->xadj[i]; j < H->xadj[i+1]; j++)
            {
                if(H->adjncy[j] > i)
                {
                    F->InsertArc(i+shift, H->adjncy[j]);
                    F->InsertArc(H->adjncy[j]+shift, i);
                }
            }	
        }		

        // set the capacities on the edges to be between 0 and 1
        static_cast<sparseRepresentation*>((*F).Representation())->SetCUCap(1);
        static_cast<sparseRepresentation*>((*F).Representation())->SetCLCap(0);
        static_cast<sparseRepresentation*>((*F).Representation())->SetCDemand(0);

        //run the max flow.
        flowValue = F->MaxFlow(source, sink);

        if(ts_debug)
        {
            cout << "a1 = " << a1 << " a2 =  " << a2 << " b1 = " << b1 << endl;
            cout << "Max Flow " << flowValue << endl;
            //tried using print_message, but %d and %f don't print a TCap correctly.
            //also, printing to stderr failed.
        }

        //This is totally a cheat - Max Flow is supposed to populate the distance labels, 
        //but it fails to do so. Discovered from shortestpath code in Goblin 
        //that NodeColours holds the same information.
        dist = F->GetNodeColours();

        //if this works out we need to set found to true.

        vector<int> S_vec(H->get_num_nodes(), false);
        vector<int> T_vec(H->get_num_nodes(), false);
        list<int> middle_set;
        bool left, right;
        //find the sets X and Y
        //first we find the partition of the nodes of H.

        //We need to set the nodes in the source & sink
        S_vec[u] = true; 
        S_vec[H->adjncy[a1]] = true; 
        S_vec[H->adjncy[a2]] = true; 
        T_vec[u] = true;
        T_vec[H->adjncy[b1]] = true;

        for(i = 0; i < shift; i++)
        {
            if( dist[i]!=NoNode && dist[shift+i]==NoNode) //dist[i]==NoNode && dist[shift+i]!=NoNode )
            {
                // These are the nodes in the new middle set that are in both S and T
                // since the actual node and its dummy copy are on different sides!!!!!
                middle_set.push_back(i);
                S_vec[i]=true; 
                T_vec[i]=true;
            }
            else
            {
                // The node lives in either S or T
                if(dist[i]==NoNode)
                    T_vec[i]=true;
                else
                    S_vec[i]=true;
            }
        }


        // Classify the BDTreeEdges adjacent to our node.
        for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
        {
            left = right = true;
            print(1,this->edges[*ii].middle_set);
            // Look to see if edge[*ii] has a middle set containing nodes on either side of the cut
            for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
            {
                if(!S_vec[GtoH[*jj]])
                    left=false;
                if(!T_vec[GtoH[*jj]]) 
                    right=false;
            }
            if(left && !right)
            {
                X->push_back(*ii); 
                print_message(1, "Left!\n");
            }
            else if(right)
            { 
                Y->push_back(*ii);
                print_message(1, "Right!\n");
            }
            else
            {
                // The edge has a middle set who straddles the cut
                //this invalidates our 2-separation.			
                fatal_error("Not really a two separation!\n");
            }
        }

        //the question is now whether there is a vertex of H adjacent only
        //to edges which are not reachable - this is equivalent to there being 
        //a vertex in the union of the middle sets of the edges in Y which is 
        //not also in the union of the middle sets of edges in X.

        //Try calculating this based on the middle set graph partition alone.
        for(i = 0; i < H->get_num_nodes(); i++)
        {
            if(S_vec[i] == false && T_vec[i]== true)
            {
                print_message(0,"Vertex %d is %d in S and %d in T\n", i, S_vec[i], T_vec[i]);
                found = true;
                break;
            }
        }
        //otherwise, we have to clear X,Y and move on.
        if(found == false)
        {
            X->clear();
            Y->clear();
            testnum++;
        }    
    }//end of while loop

    if(!found)
    {
        print_message(0, "Looking for a 2-sep with e,f together\n");
        //We need to look for a separation where we don't care where
        // g ends up, but we have e,f on the same side, order 2, and A minimal.
        // we will let every other edge of H take a turn at forming the sink. 
        int sink_edges = 0;
        while(found == false && sink_edges < H->get_num_edges()){ 
            print_message(1, "Searching with %d of %d as sink\n", sink_edges, H->get_num_edges());
            F = new sparseDiGraph(n, *CT);
            //insert the dummy arcs to do node-disjoint paths.
            for(k = 0; k < H->get_num_nodes(); k++)
            {
                F->InsertArc(k, k+shift);
            }

            //insert the arcs from the source to the 3 vertices defined by e and f
            F->InsertArc(source, u);
            F->InsertArc(source, H->adjncy[e]);
            F->InsertArc(source, H->adjncy[f]);

            //Oct 19 - patched, but now it can run up to E-2 max flows, 
            //where E is the number of edges in H (could be up to M^2/2 
            //if M is the number of vertices in a middle set of some 
            //edge that touches v.

            //As of now, the only way I see to fix this is test 
            //e,f versus p as a sink where p loops over all 
            //the other edges in the graph. 
            //We could possibly also add failed p's to the source side, 
            //but this is not implemented now.

            //we now need to add the edges of H and our sink edges.
            // The neighbors of node i are adjncy[xadj[i]],adjncy[xadj[i]+1],...adjncy[xadj[i+1]-1]	
            //only add edges (i,j) with j > i, so we don't have to look at i = H->get_num_nodes()-1.
            //we will also insert arcs from the current specified edge 
            //(given by our index sink_edges) to the sink
            int edge_num = 0;
            for(i = 0; i < H->get_num_nodes()-1; i++)
            {
                for(j = H->xadj[i]; j < H->xadj[i+1]; j++)
                {
                    if(H->adjncy[j] > i)
                    {
                        F->InsertArc(i+shift, H->adjncy[j]);
                        F->InsertArc(H->adjncy[j]+shift, i);
                        if(edge_num == sink_edges){
                            F->InsertArc(i+shift, sink);
                            F->InsertArc(H->adjncy[j]+shift, sink);
                        }
                        edge_num++;
                    }
                }	
            }		

            // set the capacities on the edges to be between 0 and 1
            static_cast<sparseRepresentation*>((*F).Representation())->SetCUCap(1);
            static_cast<sparseRepresentation*>((*F).Representation())->SetCLCap(0);
            static_cast<sparseRepresentation*>((*F).Representation())->SetCDemand(0);

            //run the max flow.
            flowValue = F->MaxFlow(source, sink);
            print_message(1, "Flow value %d\n", (int)flowValue);

            if(flowValue == 2) // we want a 2-separation. 
            {
                //This is totally a cheat - Max Flow is supposed to populate the distance labels, 
                //but it fails to do so. Discovered from shortestpath code in Goblin 
                //that NodeColours holds the same information.
                dist = F->GetNodeColours();

                //if this works out we need to set found to true.

                vector<int> S_vec(H->get_num_nodes(), false);
                vector<int> T_vec(H->get_num_nodes(), false);
                list<int> middle_set;
                bool left, right;

                //We need to set the nodes in the source - we don't know what's in the sink.
                S_vec[u] = true; 
                S_vec[H->adjncy[e]] = true; 
                S_vec[H->adjncy[f]] = true; 


                //find the sets X and Y
                //first we find the partition of the nodes of H.
                for(i = 0; i < shift; i++)
                {
                    if( dist[i]!=NoNode && dist[shift+i]==NoNode )
                    {
                        // These are the nodes in the new middle set that are in both S and T
                        // since the actual node and its dummy copy are on different sides!!!!!
                        middle_set.push_back(i);
                        S_vec[i]=true; 
                        T_vec[i]=true;
                    }
                    else
                    {
                        // The node lives in either S or T
                        if(dist[i]==NoNode)
                            T_vec[i]=true;
                        else
                            S_vec[i]=true;
                    }
                }

                // Classify the BDTreeEdges adjacent to our node.
                for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
                {
                    //xe should be the edge index we're looking at.
                    left = right = true;
                    print(1,this->edges[*ii].middle_set);
                    // Look to see if edge[*ii] has a middle set containing nodes on either side of the cut
                    for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
                    {
                        if(!S_vec[GtoH[*jj]])
                            left=false;
                        if(!T_vec[GtoH[*jj]]) 
                            right=false;
                    }
                    if(left && !right)
                    {
                        X->push_back(*ii); 
                        print_message(1, "Left!\n");
                    }
                    else if(right)
                    { 
                        Y->push_back(*ii);
                        print_message(1, "Right!\n");
                    }
                    else
                    {
                        // The edge has a middle set who straddles the cut
                        //this invalidates our 2-separation.
                        fatal_error("Not really a two separation!\n");
                    }
                }

                //the question is now whether there is a vertex of H adjacent only
                //to edges which are not reachable - this is equivalent to there being 
                //a vertex in the union of the middle sets of the edges in Y which is not
                //also in the union of the middle sets of edges in X.

                //Try calculating this based on the middle set graph partition alone.
                for(i = 0; i < H->get_num_nodes(); i++)
                {
                    if(S_vec[i] == false && T_vec[i]== true)
                    {
                        print_message(0,"Vertex %d is %d in S and %d in T\n", i, S_vec[i], T_vec[i]);
                        found = true; 
                        break;
                    }
                }

                if(found == false)
                {
                    X->clear();
                    Y->clear();
                }    
            }//flow value = 2 
            else //found a 3-flow 
                found = false;
            delete F;
            sink_edges++;
        }//end of while loop
    }//if !found
    H->free_CRS();
    delete H;
    return found;
}//end of two_separation



/**
* Helper function to check for a three separation in the branch decomposition
* at a specified vertex v. If one is found, X and Y are populated 
* with edge indices so that |M(X,Y)| = 3.
* returns true if one is found, false otherwise.
* NOT FINISHED!
*/
bool BDTree::three_separation(int v, list<int> *X, list<int> *Y, goblinController *CT)
{
    list<int>::iterator ii, jj;
    bool ts_debug = false;
    vector<int> edge_index;
    vector<int> HtoG; 
    vector<int> GtoH;
    bool found = false;

    //We first build an auxiliary graph H, which is symmetric with sorted adj lists
    //and a populated degree array. 
    // HtoG, and GtoH encode the information about how the 
    //branch decomposition edges and middle sets map into the graph H. 
    // Important that degree array is accurate and neighbors are sorted.
    Graph *H = this->construct_middle_set_graph(v, &HtoG, &GtoH, false);

    // PROBLEMS:
    //I think when 4 nodes occur, we can get a 'trivial' 3-separation that will continue
    //trying to split indefinitely. If H has 3 nodes, we run into trouble because it's not 3-connected. 
    //I'm not sure why we can get H with 3 nodes, but this is happening sometimes right now.
    if(H->get_num_nodes() <= 4)
    {
	print_message(0, "H is too small!\n");
        delete H;
        return false;
    }

    if(ts_debug)
        cerr << *H;

    int u = -1;
    int i,k;
    int n1,n2,n3;

    //check if there's a vertex of degree 3 with two adjacent neighbors
    list<int> deg3nodes;
    for(i = 0; i < H->get_num_nodes(); i++)
    {
        if(H->get_degree(i) == 3)
            deg3nodes.push_back(i);	    
    }
    for(ii = deg3nodes.begin(); ii != deg3nodes.end(); ++ii)
    {
        print_message(1, "looking at node %d\n", *ii);
        //get the neighbors of our node.
        jj = H->nodes[*ii].nbrs.begin();
        n1 = *jj; 
        ++jj;
        n2 = *jj;
        ++jj; 
        n3 = *jj;
        print_message(1, "neighbors %d %d %d\n", n1, n2, n3);
        //look for n2 and n3 in neighbors of n1
        jj = H->nodes[n1].nbrs.begin(); 
        while(jj != H->nodes[n1].nbrs.end() && *jj < n2)
            ++jj; 
        if(jj == H->nodes[n1].nbrs.end())
            break;
        else if(*jj == n2)
        {
            found = true; 
            u = HtoG[*ii];
            break;
        }
        else
        {
            while(jj != H->nodes[n1].nbrs.end() && *jj < n3)
                ++jj; 
            if(*jj == n3)
            {
                found = true; 
                u = HtoG[*ii];
                break;
            }
        }	    
    }
    if(found)
    {
        print_message(1, "Selected u = %d\n", u);
        bool cl;
        //We want to set X to be all edges where u is in the middle set, and Y 
        //to be the rest of the edges incident with v
        for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
        {
            cl = false;
            //search the middle sets for u
            jj=this->edges[*ii].middle_set.begin();
            while (jj!=this->edges[*ii].middle_set.end() && *jj <= u)
            {
                if(*jj == u)
                {
                    X->push_back(*ii); 
                    cl = true;
                }
                ++jj;
            }
            if(cl == false)
                Y->push_back(*ii);
        }
        print_message(0, "found a type 1 3-separation!\n"); 
        print(1, *X); 
        print(1, *Y);
        delete H; 
        return found;
    }
    else
    {
     	//if not, choose a vertex u arbitrarily - x in Cook/Seymour paper (CSTSP)
	int x = rand_int(0, H->get_num_nodes()-1);

	print_message(1, "Chose %d as an arbitrary vertex in H\n", x);
	
        if(ts_debug)
            print_message(0, "Selected x = %d\n", x);

	//we'll use two distinct neighbors of x
	jj = H->nodes[x].nbrs.begin(); 
	int w = *jj;
	++jj;
	int y = *jj;

	//Find another neighbor of y which is not equal to x or w, call it z.
	int z =-1; //z in CSTSP
	jj = H->nodes[y].nbrs.begin(); 
	while(jj != H->nodes[y].nbrs.end() && z == -1)
	    {
		if(*jj != x && *jj != w)
		    z = *jj; 
		else 
		    ++jj;
	    }
	
	print_message(1, "x = %d, w = %d, y = %d, z = %d\n", x, w, y, z);
	print_message(1, "deg x = %d, deg w = %d, deg y = %d, deg z = %d\n", H->get_degree(x), H->get_degree(w), H->get_degree(y), H->get_degree(z));
		      
	//We need to find a 3 separation (A,B) so left L and right R both meet {w,x,y,z} and left and right both have 
	//size at least two. We use M to denote the middle set (that is L \cup M are endpoints of A, R \cup M are endpoints of B, and 
	// (L,R,M) is a partition of the vertex set into three pieces.
	// Recall that we have w -- x -- y -- z (-- are edges in H).
	//There are three cases (note that A & B are interchangeable, so we always consider when something is in A).
	// 1: x is in L. One can see that w and y must also occur in (L \cup M). If y were in L, then z would be in L \cup M, which 
	//    contradicts our assuptions. Thus y is in M and z must be in R for both L and R to meet {w,x,y,z}. 
	// NEW: we set the source to be x and the sink to be z.
	// 2: y is in L. We may assume x is in M (see case 1). Then z is in L \cup M and w must be in R. 
	// NEW: we set the source to be y and the sink to be w.
	// 3: w is in L. We may assume x, y are in M (cases 1 & 2). Then z must lie in R to meet our criteria. 
	// NEW: we set the source to be w and the sink to be z.
	
	// In all cases, if we find a 3-separation we need to check that both L and R have size at least 2.
	
	//We'll now create a goblin graph on which to run max flow. 	
	int F_num_vertices = 2*(H->get_num_nodes());
	TNode n = F_num_vertices;
	    
	//we'll construct the goblin graph directly. 
	sparseDiGraph *F; 
	TCap flowValue;
	TNode *dist;
	int source;
	int sink;
	
	int cases = 1;	
	while(found == false && cases < 4){
	    F = new sparseDiGraph(n, *CT);
	    
	    print_message(1, "In case %d\n", cases);

	    // We will need two nodes for every non source/sink node
	    // Every such node i will have a copy with index (i + shift)
	    int shift=H->get_num_nodes(); 
	    //source and sink are nodes of H and we will ignore their shifted versions.
		
	    //We need to look for a separation where right and left meet 
	    // {u,e,f,g} and both right and left have size at least 2.
	    // Here we fill in the appropriate choices for source and sink.

	    if(cases == 1) //NEW: source = x, sink = z
		{
		    source = x; 
		    sink = z; 
		}
	    else if(cases == 2) //NEW: source = y, sink = w 
		{
		    source = y; 
		    sink = w; 	    
		}
	    else // NEW: source = w, sink = z
		{
		    source = w; 
		    sink = z; 	    
		}
	    
	    
	    //insert the dummy arcs to do node-disjoint paths for all non source/sink nodes.
	    for(k = 0; k < H->get_num_nodes(); k++)
		{
		    if(k != source && k != sink)
			F->InsertArc(k, k+shift);
		}
	    
	    //we now need to add the edges of H.
	    //only add edges (i,j) with j > i, so we don't have to look at i = H->get_num_nodes()-1.
	    //if one endpoint is the source or sink we don't add the appropriate shifted version.
	    for(i = 0; i < H->get_num_nodes(); i++)
		{
		    //handle source & sink separately
		    if(i == sink)
			{
			    print_message(1, "handling sink\n");
			    for(jj = H->nodes[i].nbrs.begin(); jj != H->nodes[i].nbrs.end(); ++jj)
				if(*jj != sink)
				    F->InsertArc(*jj+shift, i);  
			}
		    else if(i == source)
			{
			    print_message(1, "handling source\n");
			    for(jj = H->nodes[i].nbrs.begin(); jj != H->nodes[i].nbrs.end(); ++jj)
				F->InsertArc(i, *jj);  
			}
		    else
			{
			    print_message(1, "handling vertex %d\n", i);
			    for(jj = H->nodes[i].nbrs.begin(); jj != H->nodes[i].nbrs.end(); ++jj)
				{
				    //don't deal with edges to the source or sink
				    if(*jj > i && *jj != source && *jj != sink)
					{
					    F->InsertArc(i+shift, *jj);
					    F->InsertArc(*jj+shift, i);
					}
				}	
			}
		}		
	    
	    // set the capacities on the edges to be between 0 and 1
	    static_cast<sparseRepresentation*>((*F).Representation())->SetCUCap(1);
	    static_cast<sparseRepresentation*>((*F).Representation())->SetCLCap(0);
	    static_cast<sparseRepresentation*>((*F).Representation())->SetCDemand(0);
	    
	    //run the max flow.
	    flowValue = F->MaxFlow(source, sink);
	    print_message(1, "Flow value = %d\n", (int)flowValue);
	    
	    if(flowValue == 3) // we want a 3-separation. 
		{
		    found = true; // we will set this to false if the separation isn't valid.    
		    //This is totally a cheat - Max Flow is supposed to populate the distance labels, 
		    //but it fails to do so. Discovered from shortestpath code in Goblin 
		    //that NodeColours holds the same information.
		    dist = F->GetNodeColours();
		    
		    //if this works out we need to set found to true.
		    vector<int> S_vec(H->get_num_nodes(), false);
		    vector<int> T_vec(H->get_num_nodes(), false);
		    //list<int> middle_set;
		    bool left, right;
		    
		    //find the sets X and Y
		    
		    //first we find the partition of the nodes of H.
		    for(i = 0; i < shift; i++)
			{
			    //source & sink don't play by these rules.
			    if(i != source && i != sink)
				{ 
				    
				    if( dist[i]!=NoNode && dist[shift+i]==NoNode )
					{
					    // These are the nodes in the new middle set that are in both S and T
					    // since the actual node and its dummy copy are on different sides!!!!!
					    //middle_set.push_back(i);
					    S_vec[i]=true; 
					    T_vec[i]=true;
					}
				    else
					{
					    // The node lives in either S or T
					    if(dist[i]==NoNode)
						T_vec[i]=true;
					    else
						S_vec[i]=true;
					}
				}
			}

		    //print_message(0, "Middle set: ");
		    //print(0, middle_set);
		    
		    //We need to set the nodes in the source and sink in S_vec and T_vec.
		    if(cases == 1) //source x, sink z
			{
			    S_vec[x]=true;
			    T_vec[z]=true;
			}
		    else if(cases == 2) //source = y, sink = w
			{
			    S_vec[y]=true;
			    T_vec[w]=true;
			}
		    else // case 3 - source = w, sink = z
			{
			    S_vec[ w]=true;
			    T_vec[z]=true;
			}		
		    
		    // Classify the BDTreeEdges adjacent to our node.
		    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
			{
			    left = right = true;
			    print(1,this->edges[*ii].middle_set);
			    // Look to see if edge[*ii] has a middle set containing nodes on either side of the cut
			    for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
				{
				    if(!S_vec[GtoH[*jj]])
					left=false;
				    if(!T_vec[GtoH[*jj]]) 
					right=false;
				}
			    if(left)
				{
				    X->push_back(*ii); 
				    print_message(1, "Left!\n");
				}
			    else if(right)
				{ 
				    Y->push_back(*ii);
				    print_message(1, "Right!\n");
				}
			    else
				{
				    // The edge has a middle set who straddles the cut
				    //this invalidates our 3-separation.
				    fatal_error("Not really a three separation!\n");
				}
			}
		    
		    // Checks for a 'valid' 3-separation:
		    // need: left and right both have size at least two (there are at least 2 
		    // vertices in H adjacent only to edges that are reachable/not reachable). 
		    // also, must ensure appropriate members of {w,x,y,z} ended up in L and R
		    
		    //print_message(0, "Checking for L,R >=2 \n");
		    //We need to find at least two vertices with S_vec = true, T_vec = false, and two with the opposite.
		    //Try calculating this based on the middle set graph partition alone.
		    int numL = 0, numR = 0;
		    i = 0; 
		    while(i < H->get_num_nodes() && (numL < 2 || numR < 2))
			{
			    if(S_vec[i] == false && T_vec[i]== true)
				{
				    print_message(1,"Vertex %d is %d in S and %d in T\n", i, S_vec[i], T_vec[i]);
				    numL++;
				}
			    else if(S_vec[i] == true && T_vec[i]== false)
				{
				    print_message(1,"Vertex %d is %d in S and %d in T\n", i, S_vec[i], T_vec[i]);
				    numR++;
				}
			    i++;
			}//end while loop
		    
		    if(numL < 2 || numR < 2)
			{
			    found = false;
			    print_message(1, "Failed to have 2 vertices on L or R\n");
			}    
		    
		    if(found == false)
			{
			    X->clear();
			    Y->clear();
			}
		    else
			print_message(0, "Found a type 2 3-separation\n");
		}//if (flowvalue == 3)
	    else 
		found = false;
    	    delete F;
	    cases = cases+1;
	}//while cases < 4 & found == false
    }//end of else
    
    //will need one more case to deal with A \subseteq {wx,xy,yz}.
    
    delete H;
    return found;    
}//end of three_separation




/**
* Take the middle set union of all non-interior BDTreeEdges adjacent to BDTreeNode v
* and create H, the subgraph of G induced by these nodes.  Assumes that v has exactly one
* "interior" BDTreeEdge.  Adds a source and sink node to
* H where the source is adjacent to all nodes in middle set union(A) and sink is adjacent to middle set
* union(B).  Add dummy nodes so that all nodes have effective capacity 1 and solve max flow.
* At this point, you know the nodes in the middle set union on each side of the split, and this
* allows you to partition the BDTreeEdges into two halves.  CA and CB are filled as before.  Note that
* the interior edge will not be added to CA or CB and the value of the max flow is exactly the value
* of the middle set of an edge connecting the two newly spawned nodes.  However, when we calculate the
* actual middle set of the edges joining v with the two newly created nodes, it will possibly be larger
* than the value of max flow since we have to account for the interior edge which is ignored in the 
* max flow calculation...
*/
void BDTree::spawn_maxflow_partition(int v, list<int> *A, list<int> *B, list<int>*CA, list<int> *CB, 
    goblinController *CT)
{
    int source,sink,i,j;
    list<int> S, T;
    list<int>::iterator ii,jj; 

    // Find the middle set union of each edge partition - A and B are lists of edges in the 
    // BDTreeEdge Array
    this->middle_set_union(A, &S);
    this->middle_set_union(B, &T);

    // Find the interior_edge - this is only relevant if we are spawning two new nodes from v
    // and is meaning less if we are actually splitting v
    int interior_edge=-1;
    if(this->num_interior_nodes!=1)
    {
        for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
        {
            if(this->edges[*ii].start==v)
                j=this->edges[*ii].end;
            else
                j=this->edges[*ii].start;
            if(this->nodes[j].graph_edge_start==-1 && this->nodes[j].graph_edge_end==-1)
            {
                interior_edge=*ii;
                break;
            }
        }
        if(interior_edge==-1)
            fatal_error("%s:  couldn't find interior edge??\n",__FUNCTION__);
    }

    print_message(0, "S(ms_union(A)):\n");
    print(0, S);
    print_message(0, "T(ms_union(B)):\n");
    print(0, T);

    // Find the # of nodes in the middle set of all edges adjacent to v - but not the interior edge
    list<int> leaf_middle_set=this->nodes[v].edges;
    if(interior_edge!=-1)
        leaf_middle_set.remove(interior_edge);
    // C must contain A U B
    list<int> C;
    this->middle_set_union(&(leaf_middle_set),&C);
    print_message(0, "C:\n");
    print(0, C);

    // The maxflow graph is the subgraph induced by C with new edges added between source and 
    // members of S, and between sink and members of T
    int maxflow_n=C.size();

    // Create vectors to do the labeling
    // perm[i] is the index of the i-th node in the subgraph
    // iperm[k] is the index of node k from the original graph in the subgraph 
    vector<int> perm(maxflow_n,-1);
    vector<int> iperm(this->G->get_num_nodes(),-1); 
    //the middle set union can only contain nodes of G, why +1?
    C.sort();
    i=0;
    for(ii=C.begin();ii!=C.end();++ii)
    {
        perm[i]=*ii;
        iperm[*ii]=i;
        i++;
    }

    // We will need two nodes for every non source/sink node
    // Every such node i will have a copy with index (i + shift)
    int shift=maxflow_n;

    // Now account for doubling the # of nodes
    maxflow_n= 2*maxflow_n;

    //last two nodes will be source and sink
    source = maxflow_n;
    sink = maxflow_n+1;
    maxflow_n+=2;  

    // Goblin graph creation
    TNode n = maxflow_n; 
    print_message(0,"Creating GOBLIN graph with %d nodes, source=%d, sink=%d\n",
        maxflow_n,source,sink);
    sparseDiGraph H(n, *CT);

    int maxflow_m=0;
    vector<bool> S_vec(this->G->get_num_nodes(),false);
    vector<bool> T_vec(this->G->get_num_nodes(),false);
    // First, create edges at the source and sink. 
    for(ii = S.begin(); ii != S.end();++ii)
    {
        print_message(1,"Adding edge %d-%d\n",source, iperm[*ii]);
        H.InsertArc(source, iperm[*ii]);
        maxflow_m++;
        // added sep.21
        S_vec[*ii]=true;
    }

    for(ii = T.begin(); ii != T.end();++ii)
    {
        print_message(1,"Adding edge %d-%d\n",iperm[*ii],sink);  
        H.InsertArc(iperm[*ii]+shift,sink);  
        maxflow_m++;
        // added sep.21
        T_vec[*ii]=true;
    }

    // Done w/ source/sink edges

    // Now add the edges in the induced subgraph of G to the Goblin graph
    // Create a bool vector to make it easy to see if a vertex should be in H
    vector<bool> in_H(this->G->get_num_nodes(),false);
    for(ii=C.begin();ii!=C.end();++ii)
        in_H[*ii]=true;

    for(i=0;i<this->G->get_num_nodes();i++)
    {
        if(in_H[i])
        {
            // Add dummy edge
            H.InsertArc(iperm[i],iperm[i]+shift);
            maxflow_m++;
            for(ii=this->G->nodes[i].nbrs.begin();ii!=this->G->nodes[i].nbrs.end();++ii)
            {
                if( in_H[*ii] && (iperm[i] < iperm[*ii]) )
                {
                    // Add the edge i-*ii to Goblin graph
                    print_message(1,"Adding edge %d-%d\n",iperm[i],iperm[*ii]);  
                    H.InsertArc(iperm[i]+shift,iperm[*ii]);
                    H.InsertArc(iperm[*ii]+shift,iperm[i]);
                    maxflow_m+=2;
                }
            }
        }
    }

    print_message(0,"GOBLIN graph has %d edges\n",maxflow_m);

#if 0
    H.ExportToAscii("test.txt", 0); 
#endif
    // set the capacities on the edges to be between 0 and 1
    static_cast<sparseRepresentation*>(H.Representation()) -> SetCUCap(1);
    static_cast<sparseRepresentation*>(H.Representation()) -> SetCLCap(0);
    static_cast<sparseRepresentation*>(H.Representation()) -> SetCDemand(0);

    // Run the max flow 
    TCap flowValue = H.MaxFlow(source,sink);
    print_message(0,"Max Flow is %f\n",flowValue);

    //This is totally a cheat - Max Flow is supposed to populate the distance labels, but 
    //it fails to do so. Discovered from shortestpath code in Goblin that NodeColours holds
    //the same information.
    TNode* dist = H.GetNodeColours();

    list<int> middle_set;

    for(i = 0; i < shift; i++)
    {

        print_message(10,"Checking perm[%d]=%d: (%d,%d)\n",i,perm[i],dist[i],dist[shift+i]);
        if( (dist[i]==NoNode && dist[shift+i]!=NoNode) ||
            (dist[i]!=NoNode && dist[shift+i]==NoNode) )
        {
            // These are the nodes in the new middle set that are in both S and T
            // since the actual node and its dummy copy are on different sides!!!!!
            //print_message(0,"Adding middle set node %d\n",i);
            //middle_set.push_back(perm[i]);
            S_vec[perm[i]]=T_vec[perm[i]]=true;
        }
        else
        {
            // The node lives in either S or T
            if(dist[i]==NoNode)
                T_vec[perm[i]]=true;
            else
                S_vec[perm[i]]=true;
        }
    }

    // Compute the middle set here - it should ALWAYS be the same as max flow
    for(i = 0; i < shift; i++)
    {
        if(S_vec[perm[i]] && T_vec[perm[i]])
        {
            middle_set.push_back(perm[i]);
        }
    }

    print_message(0,"New middle set has size %d\n",middle_set.size());
    print(0,middle_set);

    // Check for equality!
    if((int)flowValue != (int)middle_set.size())
        fatal_error("%s:  flow value %f != ms size %d\n",(int)flowValue, (int)middle_set.size());

    int num_s=0, num_t=0;
    print_message(0, "left: ");
    for(i=0;i<this->G->get_num_nodes();i++)
    {
        // changed to perm[i], sep. 21
        if(S_vec[i])
        {
            print_message(0,"%d ",i);
            num_s++;
        }
    }
    print_message(0,"\n");
    print_message(0, "right: ");
    for(i=0;i<this->G->get_num_nodes();i++)
    {
        // changed to perm[i], sep. 21
        if(T_vec[i])
        {
            print_message(0,"%d ",i);
            num_t++;
        }
    }
    print_message(0,"\n");
    print_message(1,"S: %d elements\nT: %d elements\n",num_s,num_t);

    // Now we can classify all the edges - just fill CA and CB with them
    //these could be smaller.
    vector<bool> in_A(this->num_edges+1,false);
    vector<bool> in_B(this->num_edges+1,false);

    for(ii=A->begin();ii!=A->end();++ii)
        in_A[*ii]=true;
    for(ii=B->begin();ii!=B->end();++ii)
        in_B[*ii]=true;

    // Classify the BDTreeEdges adjacent to v
    bool left,right;
    CA->clear(); CB->clear();
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        if(!in_A[*ii] && !in_B[*ii])
        {
            // This edge is not yet classified
            print(1,this->edges[*ii].middle_set);

            // A non-interior edge's middle set should always exist entirely in S or entirely in T.
            left=right=true;
            // Look to see if edge[*ii] has a middle set containing nodes on either side of the cut
            for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
            {
                if(!S_vec[*jj])
                    left=false;
                if(!T_vec[*jj]) 
                    right=false;
            }
            if(left)
                CA->push_back(*ii); 
            else if(right) 
                CB->push_back(*ii);
            else
            {
                if(*ii != interior_edge)
                {
                    print_message(0,"Non-interior edge %d straddles cut!!!\n",*ii);
                    cerr<<this->edges[*ii];

                    fatal_error("%s:  Non-interior edge %d straddling cut with %d interior nodes in BDTree??\n",
                        __FUNCTION__,*ii,this->num_interior_nodes);
                }
            }
        }
    }

}

/**
* NEW: Oct 12, 2010
* Take the middle set union of all BDTreeEdges adjacent to BDTreeNode v. These will form the vertices of 
* a new graph H, to which we add a source S and a sink T.
* The source is adjacent to nodes in middle set union(A) and sink is adjacent to middle set
* union(B).  
* To get the remaining edges in H, we take the middle set of each unclassified BDTreeEdges (in vector U), and 
* create a clique in H. [Note: This naturally creates duplicate edges, which we may ignore?]
* We add dummy nodes so that all nodes have effective capacity 1 and solve max flow.
* At this point, you have a partition of the edges of H, which should allow 
* you to partition the BDTreeEdges into two halves.  CA and CB are filled as before.
* Returns the predicted size of the new middle set edge.
*/
int BDTree::split_maxflow_partition(int v, list<int> *A, list<int> *B, list<int>*CA, list<int> *CB, 
    goblinController *CT)
{
    int source,sink,i;
    list<int> S, T;
    list<int>::iterator ii,jj; 

    // Find the middle set union of each edge partition - A and B are lists of edges in the 
    // BDTreeEdge Array
    print_message(1,"A:\n");
    print(1,*A);
    print_message(1,"B:\n");
    print(1,*B);

    this->middle_set_union(A, &S);
    this->middle_set_union(B, &T);
    print_message(1, "MS_union(A):\n");
    print(1, S);
    print_message(1, "MS_union(B):\n");
    print(1, T);

    //find the edges adjacent to v which still need classifying. Say there are E of them stored in a vector U.
    //create a vector which is true at the indices of edges in A, B (in BDTreeEdge).
    //added 3 to potential number of edges in case tree is rooted.
    vector<bool> class_edges(2*((this->G)->get_num_edges())-3 +3, false); 
    for(ii = A->begin(); ii != A->end(); ++ii)
    {
        class_edges[*ii] = true;
    }
    for(ii = B->begin(); ii != B->end(); ++ii)
    {
        class_edges[*ii] = true;
    }

    // loop over the edges adjacent to v, check if they're in A\cup B using class_edges. If not, put them into U.
    // So U is a vector of edges that need to be assigned to A or B
    int E = (int)(this->nodes[v].edges.size() - B->size() - A->size()); 
    vector<int> U(E, 0);
    i = 0;
    for(ii = this->nodes[v].edges.begin(); ii != this->nodes[v].edges.end(); ++ii)
    {
        if(class_edges[*ii] != true)
        {
            U[i] = *ii;
            i++;
        }
    }

    // Find the # of nodes in the middle set of all edges adjacent to v
    // C must contain A U B
    list<int> C;
    this->middle_set_union(&(this->nodes[v].edges),&C);
    print_message(1, "C:\n");
    print(1, C);

    // The maxflow graph is a graph with a clique on the members of the middle sets of each unclassified edge around
    //the tree node, plus edges from the source to members of S, and the sink to members of T
    int maxflow_n=C.size();

    // Create vectors to do the labeling
    vector<int> perm(maxflow_n,-1);
    vector<int> iperm(this->G->get_num_nodes(),-1); 

    C.sort();
    i=0;
    for(ii=C.begin();ii!=C.end();++ii)
    {
        perm[i]=*ii;
        iperm[*ii]=i;
        i++;
    }

    // We will need two nodes for every non source/sink node
    // Every such node i will have a copy with index (i + shift)
    int shift=maxflow_n;

    // Now account for doubling the # of nodes
    maxflow_n= 2*maxflow_n;

    //last two nodes will be source and sink
    source = maxflow_n;
    sink = maxflow_n+1;
    maxflow_n+=2;  

    // Goblin graph creation
    TNode n = maxflow_n; 
    print_message(1,"Creating GOBLIN graph with %d nodes, source=%d, sink=%d\n",
        maxflow_n,source,sink);
    sparseDiGraph H(n, *CT);

    int maxflow_m=0;
    vector<bool> S_vec(this->G->get_num_nodes(),false);
    vector<bool> T_vec(this->G->get_num_nodes(),false);
    // First, create edges at the source and sink. 
    for(ii = S.begin(); ii != S.end();++ii)
    {
        print_message(1,"Adding edge %d-%d\n",source, iperm[*ii]);
        H.InsertArc(source, iperm[*ii]);
        maxflow_m++;
        S_vec[*ii]=true;
    }

    for(ii = T.begin(); ii != T.end();++ii)
    {
        print_message(1,"Adding edge %d-%d\n",iperm[*ii],sink);  
        H.InsertArc(iperm[*ii]+shift,sink);  
        maxflow_m++;
        T_vec[*ii]=true;
    }

    // Done w/ source/sink edges

    // Add dummy edge between each vertex and its copy
    for(i = 0; i < (int)C.size(); i++)
    {
        H.InsertArc(i, i+shift);
        maxflow_m++;
    }

    // Now add the edges corresponding to the cliques (middle sets of edges in U). 
    int xe, xu, xv;
    //loop over the edges in U
    for(i = 0; i < E; i++)
    {
        xe = U[i];
        //make a vector of the vertices in the middle set
        vector<int> M(this->edges[xe].middle_set.size(), 0);
        copy(this->edges[xe].middle_set.begin(), this->edges[xe].middle_set.end(), M.begin());
        //Loop over pairs of vertices in the middle set, adding edges to Goblin graph using iperm to look up index.
        for(xu = 0; xu < (int)M.size(); xu++)
        {
            for(xv = xu+1; xv < (int)M.size(); xv++)
            {
                print_message(1,"Adding edge %d-%d\n",iperm[M[xu]],iperm[M[xv]]);  
                H.InsertArc(iperm[M[xu]]+shift,iperm[M[xv]]);
                H.InsertArc(iperm[M[xv]]+shift,iperm[M[xu]]);
                maxflow_m+=2;
            }
        }
    }

    print_message(1,"GOBLIN graph has %d edges\n",maxflow_m);

#if 0
    H.ExportToAscii("test.txt", 0); 
#endif

    // set the capacities on the edges to be between 0 and 1
    static_cast<sparseRepresentation*>(H.Representation()) -> SetCUCap(1);
    static_cast<sparseRepresentation*>(H.Representation()) -> SetCLCap(0);
    static_cast<sparseRepresentation*>(H.Representation()) -> SetCDemand(0);

    // Run the max flow 
    TCap flowValue = H.MaxFlow(source,sink);
    print_message(1,"Max Flow is %f\n",flowValue);

    //This is totally a cheat - Max Flow is supposed to populate the distance labels, but 
    //it fails to do so. Discovered from shortestpath code in Goblin that NodeColours holds
    //the same information.
    TNode* dist = H.GetNodeColours();

    list<int> middle_set;

    //reset S & T vectors
    //    for(i = 0; i < this->G->get_num_nodes(); i++)
    //	{
    //	    S_vec[i]  = false; 
    //	    T_vec[i] = false;
    //	}

    for(i = 0; i < shift; i++)
    {
        //print_message(1,"%d %d\n",i,dist[i]);

        if( //(dist[i]==NoNode && dist[shift+i]!=NoNode) ||
            (dist[i]!=NoNode && dist[shift+i]==NoNode) )
        {
            // These are the nodes in the new middle set that are in both S and T
            // since the actual node and its dummy copy are on different sides!!!!!
            print_message(1,"Middle set node perm[%d]=%d (dist[i]=%u, dist[shift+i]=%u)\n",i,perm[i],
                dist[i], dist[shift+i]);
            middle_set.push_back(perm[i]);
            S_vec[perm[i]]=true; 
            T_vec[perm[i]]=true;
        }
        else
        {
            // The node lives in either S or T
            if(dist[i]==NoNode)
                T_vec[perm[i]]=true;
            else
                S_vec[perm[i]]=true;
        }
    }

    // Note - actual middle set could change if we have edges whose middle sets straddle cut
    print_message(1,"New middle set has size %d\n",middle_set.size());
    print(1,middle_set);

    // Classify the BDTreeEdges in U.

    // Print out left and right before classifying edges
    print_message(1, "left/right before edge classification\nleft: ");
    for(i=0;i<this->G->get_num_nodes();i++)
    {
        if(S_vec[i])
        {
            print_message(1,"%d ",i);
        }
    }
    print_message(1,"\n");
    print_message(1, "right: ");
    for(i=0;i<this->G->get_num_nodes();i++)
    {
        if(T_vec[i])
        {
            print_message(1,"%d ",i);
        }
    }
    print_message(1,"\n");
    bool left,right;
    CA->clear(); CB->clear();


    print_message(1, "A has %d edges, B has %d edges\n", A->size(), B->size());
    print_message(1, "Classifying %d edges\n", E);

    for(i = 0; i < E; i++)
    {
        xe = U[i];
        print(1,this->edges[xe].middle_set);

        // We would like for an edge's middle set to always exist entirely in S or entirely in T and believe this 
        //method guarantees this. 
        left=right=true;
        // Look to see if edge[*ii] has a middle set containing nodes on either side of the cut
        for(jj=this->edges[xe].middle_set.begin();jj!=this->edges[xe].middle_set.end();++jj)
        {
            if(!S_vec[*jj])
                left=false;
            if(!T_vec[*jj]) 
                right=false;
        }
        if(left)
        {
            CA->push_back(xe); 
            print_message(1, "Left!\n");
        }
        else if(right)
        { 
            CB->push_back(xe);
            print_message(1, "Right!\n");
        }
        else
        {
            // The edge has a middle set who straddles the cut
            print_message(0,"Edge %d straddles cut\n",xe);
            print(0,this->edges[xe].middle_set);
            int num_in_A=0, num_in_B=0;
            // Need to decide where to put this edge as its middle set contains nodes on both sides
            for(jj=this->edges[xe].middle_set.begin();jj!=this->edges[xe].middle_set.end();++jj)
            {
                if(S_vec[*jj] && !T_vec[*jj] )
                {
                    // The middle set contains a node in S that is not in T
                    num_in_A++;
                }
                if(T_vec[*jj] && !S_vec[*jj] )
                {
                    // The middle set contains a node in T that is not in S
                    num_in_B++;
                }                   
            }

            if(num_in_A>=num_in_B)
            {
                print_message(0,"Adding straddling edge to LEFT\n");
                CA->push_back(xe);
                // Now update S_vec so that we can possibly make better choices with other edges whose
                // middle sets straddle the cut
                for(jj=this->edges[xe].middle_set.begin();jj!=this->edges[xe].middle_set.end();++jj)
                    S_vec[*jj]=true;
            }
            else
            {
                print_message(0,"Adding straddling edge to RIGHT\n");
                CB->push_back(xe);
                // Now update T_vec so that we can possibly make better choices with other edges whose
                // middle sets straddle the cut
                for(jj=this->edges[xe].middle_set.begin();jj!=this->edges[xe].middle_set.end();++jj)
                    T_vec[*jj]=true;
            }
        }
    }

    fflush(stderr);
    // Debugging - calc final middle set here
    int num_s=0, num_t=0;
    print_message(1, "left: ");
    for(i=0;i<this->G->get_num_nodes();i++)
    {
        if(S_vec[i])
        {
            print_message(1,"%d ",i);
            num_s++;
        }
    }
    print_message(1,"\n");
    print_message(1, "right: ");
    for(i=0;i<this->G->get_num_nodes();i++)
    {
        if(T_vec[i])
        {
            print_message(1,"%d ",i);
            num_t++;
        }
    }
    print_message(1,"\n");
    int final_ms_size=0;
    print_message(1, "Final middle set: ");
    for(i=0;i<this->G->get_num_nodes();i++)
        if(S_vec[i] && T_vec[i])
        {
            final_ms_size++;
            print_message(1,"%d ",i);
        }
        print_message(1,"\n%s:  final middle set size=%d\n",__FUNCTION__,final_ms_size);

        return final_ms_size;
}


void BDTree::eigenvector_list(int v, list<int> *A)
{
    int i,j,k,dim;
    time_t start,stop;
    list<int>::iterator ii,jj;

    dim=(int)this->nodes[v].edges.size();

    // Allocate the matrix
    double **F=new double*[dim];
    F[0]=new double[dim*dim];
    for(i=1;i<dim;i++)
        F[i]=F[i-1]+dim;

    // Map the edge indices to 0,1,...,dim-1
    // Need the max index and the number in the union
    int num_in_union=0;
    int max_index=-1;
    int max_vertex=-1;
    vector<int> midset_vec(this->G->get_num_nodes()+1,0);
    list<int> S;
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        if(*ii>max_index)
            max_index=*ii;
        // Now go through the middle sets to find the union
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!midset_vec[*jj])
            {
                midset_vec[*jj]=1;
                num_in_union++;
                S.push_back(*jj);

                if(*jj>max_vertex)
                    max_vertex=*jj;
            }
        }
    }
    S.sort();

    // debugging
    print_message(1,"Found S\n");
    print_message(1,"S: ");
    print(1,S);

    // Create permutations to go back and forth between edge and node labels
    vector<int> edgeperm(max_index+1,-1);
    vector<int> iedgeperm(dim,-1);
    i=0;
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        print_message(1,"Edge %02d: ",i);
        print(1,this->edges[*ii].middle_set);
        edgeperm[*ii]=i;
        iedgeperm[i]=*ii;
        i++;
    }
    vector<int> nodeperm(S.size(),0);
    vector<int> inodeperm(max_vertex+1,0);
    i=0;
    for(ii=S.begin();ii!=S.end();++ii)
    {
        nodeperm[i]=*ii;
        inodeperm[*ii]=i;
        i++;
    }

    print_message(1,"Found permutations\n");

    // Allocate for a binary incidence matrix and zeroize it
    int **mid_mat;
    mid_mat=new int*[this->nodes[v].edges.size()];
    mid_mat[0]=new int[this->nodes[v].edges.size()*num_in_union];
    for(i=1;i<(int)this->nodes[v].edges.size();i++)
        mid_mat[i]=mid_mat[i-1]+num_in_union;
    memset(mid_mat[0],0,this->nodes[v].edges.size()*num_in_union*sizeof(int));

    // New - make this an array of lists
    list<int> *mid_lists=new list<int>[this->nodes[v].edges.size()];

    // Now create the rows - one for each edge
    start=clock();
    i=j=0;
    vector<int> Nv(num_in_union,0);
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        // Row i in this matrix corresponds to edge iperm[i]
        // Each row contains num_in_union binary entries - if column j is filled, then
        // that means that 
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            mid_mat[i][inodeperm[*jj]]=1;
            mid_lists[i].push_back(inodeperm[*jj]);
            // Row i contains vertex *jj
            Nv[inodeperm[*jj]]++;
        }
        mid_lists[i].sort();
        i++;
    }
    stop=clock();

    print_message(1,"Created midmat matrix in %f secs\nCreating F (num_in_union is %d)\n", 
        (double)(stop-start)/CLOCKS_PER_SEC, num_in_union);

    for(i=0;i<num_in_union;i++)
        print_message(1,"Nv %d[%d]: %d %f\n",i,nodeperm[i],Nv[i],1/(double)(Nv[i]-1));

#if 0
    // Dump the matrix
    ii=this->nodes[v].edges.begin();
    for(i=0;i<this->nodes[v].edges.size();i++)
    {
        for(j=0;j<num_in_union;j++)
            printf("%d ",mid_mat[i][j]);
        printf("| ");
        print(0,this->edges[*ii].middle_set);
        ++ii;
    }
    printf("\n");
    for(j=0;j<num_in_union;j++)
        printf("%d: %d\n",j,Nv[j]);
#endif

    // Now we can create F - it is symmetric
    // The problem here is that mid_mat is extremely sparse
    // We are better off storing this information as
    // the postions of nonzeros in each row and then just computing
    // intersections
    start=clock();
    int nnz=0;
    vector<int> I(num_in_union,0);
    vector<int>::iterator vv,ww;
    double val;

    for(i=0;i<dim;i++)
    {
        for(j=i;j<dim;j++)
        {
            F[i][j]=0;

            if(i==j)
            {
                F[i][i]=this->edges[iedgeperm[i]].middle_set.size();
                if(F[i][i]!=0)
                    nnz++;
                /*F[i][i]=0;
                for(k=0;k<num_in_union;k++)
                F[i][i]+=(double)(mid_mat[i][k]);
                print_message(1,"%d==%f?\n",this->edges[iedgeperm[i]].middle_set.size(),F[i][i]);
                if(F[i][i]!=0)
                nnz++;*/
            }
            else
            {
                // Take row i and row j in the mid_mat matrix and look for columns where
                // both rows have a 1 - the column position represents a v\in S with
                // {i,j}\in N_v
#if 0
                int intersection_size=0;
                for(k=0;k<num_in_union;k++)
                {
                    if(mid_mat[i][k]==1 && mid_mat[j][k]==1)
                    {
                        val=(-1/((double)(Nv[k])-1));
                        print_message(1,"pos=%d; Adding %f to F[%d,%d]\n",k,val,i,j);
                        F[i][j]+=val;
                        intersection_size++;
                    }
                }
#else

                // Use set intersection instead
                //I.clear();
                ww=set_intersection(mid_lists[i].begin(),mid_lists[i].end(),
                    mid_lists[j].begin(),mid_lists[j].end(),I.begin());
                // Length of intersection is ww-I.begin()

                for(int m=0;m!=ww-I.begin();m++)
                {
                    print_message(1,"I[m]=I[%d]=%d;Nv[I[m]]=%d\n",m,I[m],Nv[I[m]]);
                    val=-1/((double)Nv[I[m]]-1);
                    F[i][j]+=val;
                }
#endif
                /* for(int x=0;x<intersection_size;x++)
                print_message(0,"%d ",I[x]);

                if(intersection_size>=1)
                print_message(0,"\n%d==%d?\n",intersection_size,vv-I.begin());*/
                F[j][i]=F[i][j];
                if(F[i][j]!=0)
                    nnz+=2;
            }            
        }
    }

    // Done with mid_lists
    delete [] mid_lists;


#if 0
    printf("\n\n{");
    for(i=0;i<dim;i++)
    {
        printf("{%3.3f",F[i][0]);
        for(j=1;j<dim;j++)
        {
            printf(",%3.3f",F[i][j]); 
        }
        printf("},\n");
    }
    printf("}\n\n\n");
#endif

    double **FF= new double*[nnz];
    for(i=0;i<nnz;i++)
        FF[i]=new double[3];
    k=0;
    for(i=0;i<dim;i++)
    {
        for(j=0;j<dim;j++)
        {
            if(F[i][j]!=0)
            {
                FF[k][0]=i;
                FF[k][1]=j;
                FF[k][2]=F[i][j];
                k++;
            }
        }
    }
    stop=clock();

    print_message(1,"Created FF matrix in %f seconds\n",(double)(stop-start)/CLOCKS_PER_SEC);

    // We have the matrix - now invoke ARPACK
    double *Evals, **Evecs;
    int numeigs=2;
    Evals=new double[numeigs];
    Evecs=new double*[dim];
    for(i=0;i<dim;i++)
        Evecs[i]=new double[numeigs];
    //Evecs[0]=new double[dim];
    //Evecs[1]=new double[dim];
    print_message(1,"nnz=%d (dim^2=%d)\nInvoking ARPACK...\n",nnz,dim*dim);

    // call ARPACK
    start=clock();
    dsaupd(dim,FF,nnz,numeigs,Evals,Evecs);
    stop=clock();
    print_message(1, "ARPACK required %f seconds for %d x %d matrix\n",(double)(stop-start)/CLOCKS_PER_SEC,dim,dim);
    print_message(1,"Evals:\n");
    for(i=0;i<numeigs;i++)
        print_message(1,"%f ",Evals[i]);

    print_message(1,"Evec_2:\n");
    for(i=0;i<dim;i++)
        print_message(1,"%f ",Evecs[i][1]);

    // Copy the eigenvector so we can keep track of the edges
    vector<int_double> xx(dim);
    for(i=0;i<dim;i++)
    {
        // set xx[i].i to be the other endpoint of edge iedgeperm[i]
        //if(this->edges[iedgeperm[i]].start==v)
        //    xx[i].i = this->edges[iedgeperm[i]].end;
        //else
        //    xx[i].i = this->edges[iedgeperm[i]].start;

        //set the first part to be the edge index 
        xx[i].i = iedgeperm[i];

        xx[i].d=Evecs[i][1];
    }
    sort(xx.begin(),xx.end(),int_d_sort);
    // Now send the indices to A
    A->clear();
    for(i=0;i<dim;i++)
        A->push_back(xx[i].i);

    for(int i = 0; i < nnz; i++)
        delete[] FF[i];
    delete[] FF; 

    delete[] F[0];
    delete[] F; 
    delete[] Evals;
    for(int i = 0; i < dim; i++)
        delete[] Evecs[i];
    delete[] Evecs;
    delete[] mid_mat[0];
    delete[] mid_mat;
}

/**
* Splits the edges adjacent to BDTreeNode v into two parts, A and B, each containing
* about p times the total degree(v) edges.  Uses the Cook-Seymour eigenvector heuristic.
* If fill_extra=true, then after getting the eigenvector split, all unassigned edges whose
* middle sets are completely contained in the middle set union(A) or middle_set_union(B) are added.
* Hicks seems to use p to denote the proportion of nodes in the middle set union of the edges
* on either side of the split, and not the # of edges.  
*/
void BDTree::eigenvector_split(int v, list<int> *A, list<int> *B, double p, bool fill_extra)
{
    int i,j,k,dim;
    time_t start,stop;
    list<int>::iterator ii,jj;

    dim=(int)this->nodes[v].edges.size();

    // Allocate the matrix
    double **F=new double*[dim];
    F[0]=new double[dim*dim];
    for(i=1;i<dim;i++)
        F[i]=F[i-1]+dim;

    // Map the edge indices to 0,1,...,dim-1
    // Need the max index and the number in the union
    int num_in_union=0;
    int max_index=-1;
    int max_vertex=-1;
    vector<int> midset_vec(this->G->get_num_nodes()+1,0);
    list<int> S;
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        if(*ii>max_index)
            max_index=*ii;
        // Now go through the middle sets to find the union
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!midset_vec[*jj])
            {
                midset_vec[*jj]=1;
                num_in_union++;
                S.push_back(*jj);

                if(*jj>max_vertex)
                    max_vertex=*jj;
            }
        }
    }
    S.sort();

    // debugging
    print_message(1,"Found S\n");
    print_message(1,"S: ");
    print(1,S);

    // Create permutations to go back and forth between edge and node labels
    vector<int> edgeperm(max_index+1,-1);
    vector<int> iedgeperm(dim,-1);
    i=0;
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        print_message(1,"Edge %02d: ",i);
        print(1,this->edges[*ii].middle_set);
        edgeperm[*ii]=i;
        iedgeperm[i]=*ii;
        i++;
    }
    vector<int> nodeperm(S.size(),0);
    vector<int> inodeperm(max_vertex+1,0);
    i=0;
    for(ii=S.begin();ii!=S.end();++ii)
    {
        nodeperm[i]=*ii;
        inodeperm[*ii]=i;
        i++;
    }

    print_message(1,"Found permutations\n");

    // Allocate for a binary incidence matrix and zeroize it
    int **mid_mat;
    mid_mat=new int*[this->nodes[v].edges.size()];
    mid_mat[0]=new int[this->nodes[v].edges.size()*num_in_union];
    for(i=1;i<(int)this->nodes[v].edges.size();i++)
        mid_mat[i]=mid_mat[i-1]+num_in_union;
    memset(mid_mat[0],0,this->nodes[v].edges.size()*num_in_union*sizeof(int));

    // New - make this an array of lists
    list<int> *mid_lists=new list<int>[this->nodes[v].edges.size()];

    // Now create the rows - one for each edge
    start=clock();
    i=j=0;
    vector<int> Nv(num_in_union,0);
    for(ii=this->nodes[v].edges.begin();ii!=this->nodes[v].edges.end();++ii)
    {
        // Row i in this matrix corresponds to edge iperm[i]
        // Each row contains num_in_union binary entries - if column j is filled, then
        // that means that 
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            mid_mat[i][inodeperm[*jj]]=1;
            mid_lists[i].push_back(inodeperm[*jj]);
            // Row i contains vertex *jj
            Nv[inodeperm[*jj]]++;
        }
        mid_lists[i].sort();
        i++;
    }
    stop=clock();

    print_message(1,"Created midmat matrix in %f secs\nCreating F (num_in_union is %d)\n", 
        (double)(stop-start)/CLOCKS_PER_SEC, num_in_union);

    for(i=0;i<num_in_union;i++)
        print_message(1,"Nv %d[%d]: %d %f\n",i,nodeperm[i],Nv[i],1/(double)(Nv[i]-1));

#if 0
    // Dump the matrix
    ii=this->nodes[v].edges.begin();
    for(i=0;i<this->nodes[v].edges.size();i++)
    {
        for(j=0;j<num_in_union;j++)
            printf("%d ",mid_mat[i][j]);
        printf("| ");
        print(0,this->edges[*ii].middle_set);
        ++ii;
    }
    printf("\n");
    for(j=0;j<num_in_union;j++)
        printf("%d: %d\n",j,Nv[j]);
#endif

    // Now we can create F - it is symmetric
    // The problem here is that mid_mat is extremely sparse
    // We are better off storing this information as
    // the postions of nonzeros in each row and then just computing
    // intersections
    start=clock();
    int nnz=0;
    vector<int> I(num_in_union,0);
    vector<int>::iterator vv,ww;
    double val;

    for(i=0;i<dim;i++)
    {
        for(j=i;j<dim;j++)
        {
            F[i][j]=0;

            if(i==j)
            {
                F[i][i]=this->edges[iedgeperm[i]].middle_set.size();
                if(F[i][i]!=0)
                    nnz++;
                /*F[i][i]=0;
                for(k=0;k<num_in_union;k++)
                F[i][i]+=(double)(mid_mat[i][k]);
                print_message(1,"%d==%f?\n",this->edges[iedgeperm[i]].middle_set.size(),F[i][i]);
                if(F[i][i]!=0)
                nnz++;*/
            }
            else
            {
                // Take row i and row j in the mid_mat matrix and look for columns where
                // both rows have a 1 - the column position represents a v\in S with
                // {i,j}\in N_v
#if 0
                int intersection_size=0;
                for(k=0;k<num_in_union;k++)
                {
                    if(mid_mat[i][k]==1 && mid_mat[j][k]==1)
                    {
                        val=(-1/((double)(Nv[k])-1));
                        print_message(1,"pos=%d; Adding %f to F[%d,%d]\n",k,val,i,j);
                        F[i][j]+=val;
                        intersection_size++;
                    }
                }
#else

                // Use set intersection instead
                //I.clear();
                ww=set_intersection(mid_lists[i].begin(),mid_lists[i].end(),
                    mid_lists[j].begin(),mid_lists[j].end(),I.begin());
                // Length of intersection is ww-I.begin()

                for(int m=0;m!=ww-I.begin();m++)
                {
                    print_message(1,"I[m]=I[%d]=%d;Nv[I[m]]=%d\n",m,I[m],Nv[I[m]]);
                    val=-1/((double)Nv[I[m]]-1);
                    F[i][j]+=val;
                }
#endif
                /* for(int x=0;x<intersection_size;x++)
                print_message(0,"%d ",I[x]);

                if(intersection_size>=1)
                print_message(0,"\n%d==%d?\n",intersection_size,vv-I.begin());*/
                F[j][i]=F[i][j];
                if(F[i][j]!=0)
                    nnz+=2;
            }            
        }
    }

    // Done with mid_lists
    delete [] mid_lists;


#if ARPACK_DEBUG
    // Show the matrix in Mathematica friendly format
    printf("\n\n{");
    for(i=0;i<dim;i++)
    {
        printf("{%3.3f",F[i][0]);
        for(j=1;j<dim;j++)
        {
            printf(",%3.3f",F[i][j]); 
        }
        printf("},\n");
    }
    printf("}\n\n\n");
#endif

    double **FF= new double*[nnz];
    for(i=0;i<nnz;i++)
        FF[i]=new double[3];
    k=0;
    for(i=0;i<dim;i++)
    {
        for(j=0;j<dim;j++)
        {
            if(F[i][j]!=0)
            {
                FF[k][0]=i;
                FF[k][1]=j;
                FF[k][2]=F[i][j];
                k++;
            }
        }
    }
    stop=clock();

    print_message(1,"Created FF matrix in %f seconds\n",(double)(stop-start)/CLOCKS_PER_SEC);

    // We have the matrix - now invoke ARPACK
    double *Evals, **Evecs;
    int numeigs=2;
    Evals=new double[numeigs];
    Evecs=new double*[dim];
    for(i=0;i<dim;i++)
        Evecs[i]=new double[numeigs];
    //Evecs[0]=new double[dim];
    //Evecs[1]=new double[dim];
    print_message(1,"nnz=%d (dim^2=%d)\nInvoking ARPACK...\n",nnz,dim*dim);

    // call ARPACK
    start=clock();
    dsaupd(dim,FF,nnz,numeigs,Evals,Evecs);
    stop=clock();
    print_message(1, "ARPACK required %f seconds for %d x %d matrix\n",(double)(stop-start)/CLOCKS_PER_SEC,dim,dim);

#if ARPACK_DEBUG
    printf("Evals:\n");
    for(i=0;i<numeigs;i++)
        printf("%f ",Evals[i]);

    printf("Evec_2:\n");
    for(i=0;i<dim;i++)
        printf("%f ",Evecs[i][1]);
#endif
    // Copy the eigenvector so we can keep track of the edges
    vector<int_double> xx(dim);
    for(i=0;i<dim;i++)
    {
        // set xx[i].i to be the other endpoint of edge iedgeperm[i]
        //if(this->edges[iedgeperm[i]].start==v)
        //    xx[i].i = this->edges[iedgeperm[i]].end;
        //else
        //    xx[i].i = this->edges[iedgeperm[i]].start;

        //set the first part to be the edge index 
        xx[i].i = iedgeperm[i];

        xx[i].d=Evecs[i][1];
    }
    sort(xx.begin(),xx.end(),int_d_sort);
    //debugging
    for(i=0;i<dim;i++)
        print_message(1,"%d %f\n",xx[i].i,xx[i].d);

    // Finally ready to fill A and B
    int partition_size=(int)ceil((double)dim * p);// was floor - changed for hicks comparison
    if(partition_size<=1)
        partition_size=2;
    vector<bool> assigned(this->num_edges);
    vector<int> ms_intersection(this->G->get_num_nodes(),-1);
    list<int> ms_union;
    A->clear(); B->clear();
    for(i=0;i<partition_size;i++)
    {
#if 1
        // A gets the smallest, B gets the largest
        A->push_back(xx[i].i);
        assigned[xx[i].i]=true;
        B->push_back(xx[dim-1-i].i);
        assigned[xx[dim-1-i].i]=true;
#else
        // This is used to try and force disjoint middle set unions
        // CSG - inserting check here to make sure we have disjoint middle set unions!!
        ms_union.clear();
        this->edges[xx[i].i].middle_set.sort();
        this->middle_set_union(B,&ms_union);
        vv=set_intersection( this->edges[xx[i].i].middle_set.begin(), this->edges[xx[i].i].middle_set.end(),
            ms_union.begin(),ms_union.end(),ms_intersection.begin());
        if(vv-ms_intersection.begin() == 0)
        {
            A->push_back(xx[i].i);
            assigned[xx[i].i]=true;
        }

        ms_union.clear();
        this->edges[xx[dim-1-i].i].middle_set.sort();
        this->middle_set_union(A,&ms_union);
        vv=set_intersection( this->edges[xx[dim-1-i].i].middle_set.begin(), this->edges[xx[dim-1-i].i].middle_set.end(),
            ms_union.begin(),ms_union.end(),ms_intersection.begin());
        if(vv-ms_intersection.begin() == 0)
        {
            B->push_back(xx[dim-1-i].i);
            assigned[xx[dim-1-i].i]=true;        
        }
#endif
    }
    print_message(1,"%s:  size(A)=%d; size(B)=%d\n",__FUNCTION__,A->size(), B->size());

    print_message(1,"A MS Union before fill\n");
    list<int> ms_check;
    this->middle_set_union(A,&ms_check);
    print(1,ms_check);
    print_message(1,"B MS Union before fill\n");
    ms_check.clear();
    this->middle_set_union(B,&ms_check);
    print(1,ms_check);

    /*
    for(ii=A->begin();ii!=A->end();++ii)
    {
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].start].graph_edge_start,this->nodes[this->edges[*ii].start].graph_edge_end);
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].end].graph_edge_start,this->nodes[this->edges[*ii].end].graph_edge_end);
    }
    print_message(0,"B before fill\n");
    for(ii=B->begin();ii!=B->end();++ii)
    {
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].start].graph_edge_start,this->nodes[this->edges[*ii].start].graph_edge_end);
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].end].graph_edge_start,this->nodes[this->edges[*ii].end].graph_edge_end);
    }*/

    if(fill_extra)
    {
        // Assign all unclassified edges to a side of the split if they don't affect the 
        // middle set union of that side - this functionality used to be in maxflow but we
        // moved it here to assist in the exhaust checking        
        list<int> A_ms, B_ms;
        this->middle_set_union(A,&A_ms);
        this->middle_set_union(B,&B_ms);

        vector<bool> A_verts((this->G)->get_num_nodes(), false);
        vector<bool> B_verts((this->G)->get_num_nodes(), false);
        for(ii = A_ms.begin(); ii != A_ms.end(); ++ii)
            A_verts[*ii] = true;
        for(ii = B_ms.begin(); ii != B_ms.end(); ++ii)
            B_verts[*ii] = true;
        bool subset_A, subset_B; 
        list<int>::iterator jj;
        int num_additions=0;
        for(ii = this->nodes[v].edges.begin(); ii != this->nodes[v].edges.end(); ++ii)
        {
            if(!assigned[*ii])
            {
                subset_A = true; 
                subset_B = true;
                for(jj = this->edges[*ii].middle_set.begin(); jj != this->edges[*ii].middle_set.end(); ++jj)
                {
                    if(A_verts[*jj] == false)
                        subset_A = false; 

                    if(B_verts[*jj] == false)
                        subset_B = false;

                    if(!subset_A && !subset_B)
                        // No point to keep looking!
                        break;
                }
                if(subset_A && subset_B)
                {
                    // Here is an edge that could go either way! Pick a random side
                    print_message(1,"Edge %d can go either way!\n",*ii);
                    print(1,this->edges[*ii].middle_set);
                    int side=rand_int(0,1);
                    print_message(1,"side=%d\n",side);
                    if(side==0)
                    {
                        // Send to A
                        A->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }
                    else
                    {
                        B->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }

                }
                else
                {
                    if(subset_A)
                    {
                        A->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }
                    else if(subset_B)
                    {
                        B->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }
                }
            }
        }
        print_message(1,"%d additions made\n",num_additions);
    }


#if 0
    // Debugging - will print out properties of middle set implied by the eig. results
    vector<bool> A_vec(this->G->get_num_nodes()+1,false);
    vector<bool> B_vec(this->G->get_num_nodes()+1,false);
    int left_size=0, right_size=0, middle_size=0;
    for(ii=A->begin();ii!=A->end();++ii)
    {
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!A_vec[*jj])
            {
                A_vec[*jj]=true;
                left_size++;
            }
        }
    }
    for(ii=B->begin();ii!=B->end();++ii)
    {
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!B_vec[*jj])
            {
                B_vec[*jj]=true;
                right_size++;
            }
        }
    }

    for(i=0;i<this->G->get_num_nodes()+1;i++)
    {
        if(A_vec[i] && B_vec[i])
            middle_size++;
    }
    print_message(0,"After eigenvector_split(%d) - left_size=%d; right_size=%d; middle_size=%d\n",
        v,left_size, right_size, middle_size);
    print_message(0,"Need to add %d more edges to split\n",
        this->nodes[v].edges.size()-A->size()-B->size());
    print_message(0,"LEFT\n");
    j=0;
    for(i=0;i<this->G->get_num_nodes()+1;i++)
    {
        if(A_vec[i])
        {
            print_message(0,"%04d ",i);
            j++;
            if(j%10==0)
                print_message(0,"\n");
        }
    }
    print_message(0,"\nRIGHT\n");
    j=0;
    for(i=0;i<this->G->get_num_nodes()+1;i++)
    {
        if(B_vec[i])
        {
            print_message(0,"%04d ",i);
            j++;
            if(j%10==0)
                print_message(0,"\n");
        }
    }
    print_message(0,"\n");

#endif

    for(int i = 0; i < nnz; i++)
        delete[] FF[i];
    delete[] FF; 

    delete[] F[0];
    delete[] F; 
    delete[] Evals;
    for(int i = 0; i < dim; i++)
        delete[] Evecs[i];
    delete[] Evecs;
    delete[] mid_mat[0];
    delete[] mid_mat;
}

/**
* Splits the leaf edges adjacent to BDTreeNode v into two parts, A and B, each containing
* about p times the total degree(v) edges.  Uses the Cook-Seymour eigenvector heuristic.
* If fill_extra=true, then after getting the eigenvector split, all unassigned edges whose
* middle sets are contained in A or B are added.
* Hicks seems to use p to denote the proportion of nodes in the middle set union of the edges
* on either side of the split, and not the # of edges.  
*/
void BDTree::eigenvector_leaf_split(int v, list<int> *A, list<int> *B, double p, bool fill_extra)
{
    int i,j,k,dim;
    time_t start,stop;
    list<int>::iterator ii,jj;

    list<int> edge_list=this->nodes[v].edges;
    // Find the interior edge and remove it
    k=-1;
    for(ii=edge_list.begin();ii!=edge_list.end();++ii)
    {
        if(this->edges[*ii].start==v)
            j=this->edges[*ii].end;
        else
            j=this->edges[*ii].start;
        if(this->nodes[j].graph_edge_start==-1 && this->nodes[j].graph_edge_end==-1)
        {
            k=*ii;
            break;
        }
    }
    if(k==-1)
    {
        char err_file[100];
        sprintf(err_file,"no_interior.gviz");
        this->write_graphviz_file(false,err_file,true,true);
        fatal_error("%s:  Couldn't find interior edge for node %d (%d nbrs)?? Gviz written to no_interior.gviz\n",
            __FUNCTION__,v,this->nodes[v].edges.size());
    }
    // Remove k from the edge_list
    edge_list.remove(k);

    dim=(int)edge_list.size();

    // Allocate the matrix
    double **F=new double*[dim];
    F[0]=new double[dim*dim];
    for(i=1;i<dim;i++)
        F[i]=F[i-1]+dim;

    // Map the edge indices to 0,1,...,dim-1
    // Need the max index and the number in the union
    int num_in_union=0;
    int max_index=-1;
    int max_vertex=-1;
    vector<int> midset_vec(this->G->get_num_nodes()+1,0);
    list<int> S;
    for(ii=edge_list.begin();ii!=edge_list.end();++ii)
    {
        if(*ii>max_index)
            max_index=*ii;
        // Now go through the middle sets to find the union
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!midset_vec[*jj])
            {
                midset_vec[*jj]=1;
                num_in_union++;
                S.push_back(*jj);

                if(*jj>max_vertex)
                    max_vertex=*jj;
            }
        }
    }
    S.sort();

    // debugging
    print_message(1,"Found S\n");
    print_message(1,"S: ");
    print(1,S);

    // Create permutations to go back and forth between edge and node labels
    vector<int> edgeperm(max_index+1,-1);
    vector<int> iedgeperm(dim,-1);
    i=0;
    for(ii=edge_list.begin();ii!=edge_list.end();++ii)
    {
        print_message(1,"Edge %02d: ",i);
        print(1,this->edges[*ii].middle_set);
        edgeperm[*ii]=i;
        iedgeperm[i]=*ii;
        i++;
    }
    vector<int> nodeperm(S.size(),0);
    vector<int> inodeperm(max_vertex+1,0);
    i=0;
    for(ii=S.begin();ii!=S.end();++ii)
    {
        nodeperm[i]=*ii;
        inodeperm[*ii]=i;
        i++;
    }

    print_message(1,"Found permutations\n");

    // Allocate for a binary incidence matrix and zeroize it
    int **mid_mat;
    mid_mat=new int*[edge_list.size()];
    mid_mat[0]=new int[edge_list.size()*num_in_union];
    for(i=1;i<(int)edge_list.size();i++)
        mid_mat[i]=mid_mat[i-1]+num_in_union;
    memset(mid_mat[0],0,edge_list.size()*num_in_union*sizeof(int));

    // New - make this an array of lists
    list<int> *mid_lists=new list<int>[edge_list.size()];

    // Now create the rows - one for each edge
    start=clock();
    i=j=0;
    vector<int> Nv(num_in_union,0);
    for(ii=edge_list.begin();ii!=edge_list.end();++ii)
    {
        // Row i in this matrix corresponds to edge iperm[i]
        // Each row contains num_in_union binary entries - if column j is filled, then
        // that means that 
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            mid_mat[i][inodeperm[*jj]]=1;
            mid_lists[i].push_back(inodeperm[*jj]);
            // Row i contains vertex *jj
            Nv[inodeperm[*jj]]++;
        }
        mid_lists[i].sort();
        i++;
    }
    stop=clock();

    print_message(1,"Created midmat matrix in %f secs\nCreating F (num_in_union is %d)\n", 
        (double)(stop-start)/CLOCKS_PER_SEC, num_in_union);

    for(i=0;i<num_in_union;i++)
        print_message(1,"Nv %d[%d]: %d %f\n",i,nodeperm[i],Nv[i],1/(double)(Nv[i]-1));

#if 0
    // Dump the matrix
    ii=edge_list.begin();
    for(i=0;i<edge_list.size();i++)
    {
        for(j=0;j<num_in_union;j++)
            printf("%d ",mid_mat[i][j]);
        printf("| ");
        print(0,this->edges[*ii].middle_set);
        ++ii;
    }
    printf("\n");
    for(j=0;j<num_in_union;j++)
        printf("%d: %d\n",j,Nv[j]);
#endif

    // Now we can create F - it is symmetric
    // The problem here is that mid_mat is extremely sparse
    // We are better off storing this information as
    // the postions of nonzeros in each row and then just computing
    // intersections
    start=clock();
    int nnz=0;
    vector<int> I(num_in_union,0);
    vector<int>::iterator vv,ww;
    double val;

    for(i=0;i<dim;i++)
    {
        for(j=i;j<dim;j++)
        {
            F[i][j]=0;

            if(i==j)
            {
                F[i][i]=this->edges[iedgeperm[i]].middle_set.size();
                if(F[i][i]!=0)
                    nnz++;
                /*F[i][i]=0;
                for(k=0;k<num_in_union;k++)
                F[i][i]+=(double)(mid_mat[i][k]);
                print_message(1,"%d==%f?\n",this->edges[iedgeperm[i]].middle_set.size(),F[i][i]);
                if(F[i][i]!=0)
                nnz++;*/
            }
            else
            {
                // Take row i and row j in the mid_mat matrix and look for columns where
                // both rows have a 1 - the column position represents a v\in S with
                // {i,j}\in N_v
#if 0
                int intersection_size=0;
                for(k=0;k<num_in_union;k++)
                {
                    if(mid_mat[i][k]==1 && mid_mat[j][k]==1)
                    {
                        val=(-1/((double)(Nv[k])-1));
                        print_message(1,"pos=%d; Adding %f to F[%d,%d]\n",k,val,i,j);
                        F[i][j]+=val;
                        intersection_size++;
                    }
                }
#else

                // Use set intersection instead
                //I.clear();
                ww=set_intersection(mid_lists[i].begin(),mid_lists[i].end(),
                    mid_lists[j].begin(),mid_lists[j].end(),I.begin());
                // Length of intersection is ww-I.begin()

                for(int m=0;m!=ww-I.begin();m++)
                {
                    print_message(1,"I[m]=I[%d]=%d;Nv[I[m]]=%d\n",m,I[m],Nv[I[m]]);
                    val=-1/((double)Nv[I[m]]-1);
                    F[i][j]+=val;
                }
#endif
                /* for(int x=0;x<intersection_size;x++)
                print_message(0,"%d ",I[x]);

                if(intersection_size>=1)
                print_message(0,"\n%d==%d?\n",intersection_size,vv-I.begin());*/
                F[j][i]=F[i][j];
                if(F[i][j]!=0)
                    nnz+=2;
            }            
        }
    }

    // Done with mid_lists
    delete [] mid_lists;


#if 0
    printf("\n\n{");
    for(i=0;i<dim;i++)
    {
        printf("{%3.3f",F[i][0]);
        for(j=1;j<dim;j++)
        {
            printf(",%3.3f",F[i][j]); 
        }
        printf("},\n");
    }
    printf("}\n\n\n");
#endif

    double **FF= new double*[nnz];
    for(i=0;i<nnz;i++)
        FF[i]=new double[3];
    k=0;
    for(i=0;i<dim;i++)
    {
        for(j=0;j<dim;j++)
        {
            if(F[i][j]!=0)
            {
                FF[k][0]=i;
                FF[k][1]=j;
                FF[k][2]=F[i][j];
                k++;
            }
        }
    }
    stop=clock();

    print_message(1,"Created FF matrix in %f seconds\n",(double)(stop-start)/CLOCKS_PER_SEC);

    // We have the matrix - now invoke ARPACK
    double *Evals, **Evecs;
    int numeigs=2;
    Evals=new double[numeigs];
    Evecs=new double*[dim];
    for(i=0;i<dim;i++)
        Evecs[i]=new double[numeigs];
    //Evecs[0]=new double[dim];
    //Evecs[1]=new double[dim];
    print_message(1,"nnz=%d (dim^2=%d)\nInvoking ARPACK...\n",nnz,dim*dim);

    // call ARPACK
    start=clock();
    dsaupd(dim,FF,nnz,numeigs,Evals,Evecs);
    stop=clock();
    print_message(1, "ARPACK required %f seconds for %d x %d matrix\n",(double)(stop-start)/CLOCKS_PER_SEC,dim,dim);
    print_message(1,"Evals:\n");
    for(i=0;i<numeigs;i++)
        print_message(1,"%f ",Evals[i]);

    print_message(1,"Evec_2:\n");
    for(i=0;i<dim;i++)
        print_message(1,"%f ",Evecs[i][1]);

    // Copy the eigenvector so we can keep track of the edges
    vector<int_double> xx(dim);
    for(i=0;i<dim;i++)
    {
        // set xx[i].i to be the other endpoint of edge iedgeperm[i]
        //if(this->edges[iedgeperm[i]].start==v)
        //    xx[i].i = this->edges[iedgeperm[i]].end;
        //else
        //    xx[i].i = this->edges[iedgeperm[i]].start;

        //set the first part to be the edge index 
        xx[i].i = iedgeperm[i];

        xx[i].d=Evecs[i][1];
    }
    sort(xx.begin(),xx.end(),int_d_sort);
    //debugging
    for(i=0;i<dim;i++)
        print_message(1,"%d %f\n",xx[i].i,xx[i].d);

    // Finally ready to fill A and B
    int partition_size=(int)ceil((double)dim * p);// was floor - changed for hicks comparison
    if(partition_size<=1)
        partition_size=2;
    vector<bool> assigned(this->num_edges);
    for(i=0;i<partition_size;i++)
    {
        // A gets the smallest, B gets the largest
        A->push_back(xx[i].i);
        assigned[xx[i].i]=true;
        B->push_back(xx[dim-1-i].i);
        assigned[xx[dim-1-i].i]=true;        
    }

    // CSG, Sep. 20 - look for interior edge in A and B partition
    int nbr;
    for(ii=A->begin();ii!=A->end();++ii)
    {
        if(this->edges[*ii].end==v)
            nbr=this->edges[*ii].start;
        else
            nbr=this->edges[*ii].end;

        // Now look to see if this is an interior edge
        if(this->nodes[nbr].graph_edge_start==-1 && this->nodes[nbr].graph_edge_end==-1)
            print_message(0,"Interior edge %d on A side of split!!\n",*ii);
    }
    for(ii=B->begin();ii!=B->end();++ii)
    {
        if(this->edges[*ii].end==v)
            nbr=this->edges[*ii].start;
        else
            nbr=this->edges[*ii].end;

        // Now look to see if this is an interior edge
        if(this->nodes[nbr].graph_edge_start==-1 && this->nodes[nbr].graph_edge_end==-1)
            print_message(0,"Interior edge %d on B side of split!!\n",*ii);
    }

    print_message(1,"A MS Union before fill\n");
    list<int> ms_check;
    this->middle_set_union(A,&ms_check);
    print(1,ms_check);
    print_message(1,"B MS Union before fill\n");
    ms_check.clear();
    this->middle_set_union(B,&ms_check);
    print(1,ms_check);

    /*
    for(ii=A->begin();ii!=A->end();++ii)
    {
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].start].graph_edge_start,this->nodes[this->edges[*ii].start].graph_edge_end);
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].end].graph_edge_start,this->nodes[this->edges[*ii].end].graph_edge_end);
    }
    print_message(0,"B before fill\n");
    for(ii=B->begin();ii!=B->end();++ii)
    {
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].start].graph_edge_start,this->nodes[this->edges[*ii].start].graph_edge_end);
    print_message(0,"%d (%d-%d)\n",*ii,this->nodes[this->edges[*ii].end].graph_edge_start,this->nodes[this->edges[*ii].end].graph_edge_end);
    }*/

    if(fill_extra)
    {
        // Assign all unclassified edges to a side of the split if they don't affect the 
        // middle set union of that side - this functionality used to be in maxflow but we
        // moved it here to assist in the exhaust checking        
        list<int> A_ms, B_ms;
        this->middle_set_union(A,&A_ms);
        this->middle_set_union(B,&B_ms);

        vector<bool> A_verts((this->G)->get_num_nodes(), false);
        vector<bool> B_verts((this->G)->get_num_nodes(), false);
        for(ii = A_ms.begin(); ii != A_ms.end(); ++ii)
            A_verts[*ii] = true;
        for(ii = B_ms.begin(); ii != B_ms.end(); ++ii)
            B_verts[*ii] = true;
        bool subset_A, subset_B; 
        list<int>::iterator jj;
        int num_additions=0;
        for(ii = edge_list.begin(); ii != edge_list.end(); ++ii)
        {
            if(!assigned[*ii])
            {
                subset_A = true; 
                subset_B = true;
                for(jj = this->edges[*ii].middle_set.begin(); jj != this->edges[*ii].middle_set.end(); ++jj)
                {
                    if(A_verts[*jj] == false)
                        subset_A = false; 

                    if(B_verts[*jj] == false)
                        subset_B = false;

                    if(!subset_A && !subset_B)
                        // No point to keep looking!
                        break;
                }
                if(subset_A && subset_B)
                {
                    // Here is an edge that could go either way! Pick a random side
                    print_message(1,"Edge %d can go either way!\n",*ii);
                    print(1,this->edges[*ii].middle_set);
                    int side=rand_int(0,1);
                    print_message(1,"side=%d\n",side);
                    if(side==0)
                    {
                        // Send to A
                        A->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }
                    else
                    {
                        B->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }

                }
                else
                {
                    if(subset_A)
                    {
                        A->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }
                    else if(subset_B)
                    {
                        B->push_back(*ii);
                        assigned[*ii] = true;
                        num_additions++;
                    }
                }
            }
        }
        print_message(1,"%d additions made\n",num_additions);
    }


#if 0
    // Debugging - will print out properties of middle set implied by the eig. results
    vector<bool> A_vec(this->G->get_num_nodes()+1,false);
    vector<bool> B_vec(this->G->get_num_nodes()+1,false);
    int left_size=0, right_size=0, middle_size=0;
    for(ii=A->begin();ii!=A->end();++ii)
    {
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!A_vec[*jj])
            {
                A_vec[*jj]=true;
                left_size++;
            }
        }
    }
    for(ii=B->begin();ii!=B->end();++ii)
    {
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!B_vec[*jj])
            {
                B_vec[*jj]=true;
                right_size++;
            }
        }
    }

    for(i=0;i<this->G->get_num_nodes()+1;i++)
    {
        if(A_vec[i] && B_vec[i])
            middle_size++;
    }
    print_message(0,"After eigenvector_split(%d) - left_size=%d; right_size=%d; middle_size=%d\n",
        v,left_size, right_size, middle_size);
    print_message(0,"Need to add %d more edges to split\n",
        edge_list.size()-A->size()-B->size());
    print_message(0,"LEFT\n");
    j=0;
    for(i=0;i<this->G->get_num_nodes()+1;i++)
    {
        if(A_vec[i])
        {
            print_message(0,"%04d ",i);
            j++;
            if(j%10==0)
                print_message(0,"\n");
        }
    }
    print_message(0,"\nRIGHT\n");
    j=0;
    for(i=0;i<this->G->get_num_nodes()+1;i++)
    {
        if(B_vec[i])
        {
            print_message(0,"%04d ",i);
            j++;
            if(j%10==0)
                print_message(0,"\n");
        }
    }
    print_message(0,"\n");

#endif

    for(int i = 0; i < nnz; i++)
        delete[] FF[i];
    delete[] FF; 

    delete[] F[0];
    delete[] F; 
    delete[] Evals;
    for(int i = 0; i < dim; i++)
        delete[] Evecs[i];
    delete[] Evecs;
    delete[] mid_mat[0];
    delete[] mid_mat;
}

/**
* specialized middle_set_union allows you to record what new vertices are in the union
* that you have not seen before (where prior viewings are recorded by W)
* E is a list of indices into the BDTreeEdge array.
* V will be filled with a sorted list of the vertices that are
* in the middle set of at least one edge e\in E.
* If W != NULL, W[v] will be set to true for all v in V. 
* newv counts the number of these found vertices for which W[v] 
* was false at the start of this function.
*/
int BDTree::middle_set_union(list<int> *E, list<int> *V, vector<bool>* W, int *newv)
{
    bool fillw = true;
    if(W == NULL)
        fillw = false;

    int found = 0;
    list<int>::iterator ii,jj;
    vector<bool> midset_vec(this->G->get_num_nodes()+1,false);
    V->clear();
    for(ii=E->begin();ii!=E->end();++ii)
    {
        // Now go through the middle sets to find the union
        for(jj=this->edges[*ii].middle_set.begin();jj!=this->edges[*ii].middle_set.end();++jj)
        {
            if(!midset_vec[*jj])
            {
                midset_vec[*jj]=true;
                V->push_back(*jj);
                if(fillw)
                {
                    //set The found vector to true for this vertex *jj
                    if((*W)[*jj] == false)
                    {
                        (*W)[*jj] = true;
                        found++;
                    }
                }
            }
        }
    }
    V->sort();

    if(fillw)
        *newv = found;

    // We shouldn't have duplicates since we use the midset_vec
    // V->unique();

    return V->size();
}

/**
* E is a list of indices into the BDTreeEdge array.
* V will be filled with a sorted list of the vertices that are
* in the middle set of at least one edge e\in E.
*/
int BDTree::middle_set_union(list<int> *E, list<int> *V)
{
    return middle_set_union(E, V, NULL, NULL);

}


/*
* To make the notion of inwards precise, we
* root the branch-decomposition by selecting an arbitrary
* edge a-b of the tree and adding new tree nodes
* r and s and new edges a-s, s-b, and s-r and by
* removing the edge a-b. Let T denote the tree we
* obtain in this way, and call r the root node of T . No
* edge of G is assigned to node r, so the separations
* corresponding to a-s and s-b are the same as the
* separation for the old edge a-b; the separation corresponding
* to s-r is E(G).
*         r
* a-b ==> |
*       a-s-b 
* Furthermore, we set the edge endpoints throughout the tree so that 
* the 'end' is nearer the root, effectively forcing all edges to point 
* towards the root. 
* Returns the index of the new root node r.
*/
int BDTree::root(int edge_index)
{
    // Make sure the tree is valid
    if(!this->is_valid)
        fatal_error("%s:  Cannot root tree unless it is ternary/valid\n",__FUNCTION__);

    // Make sure the edge_index is meaningful
    if(edge_index > this->num_edges || edge_index<0)
        fatal_error("%s:  Bad edge_index %d\n",__FUNCTION__,edge_index);

    //BDS - added Aug 26
    this->is_rooted = true;

    // Remember old endpoints of this edge (a-b) and the middle set
    int a=this->edges[edge_index].start;
    int b=this->edges[edge_index].end;

    // Now remove this edge from nodes a and b
    this->nodes[a].edges.remove(edge_index);
    this->nodes[b].edges.remove(edge_index);

    // Create new nodes r and s
    int r=this->num_nodes+1;
    int s=this->num_nodes;
    this->nodes[r].id=r;
    this->nodes[s].id=s;
    this->num_nodes+=2;
    // s is an interior node, so increment that?
    this->num_interior_nodes++;

    // The former edge[edge_index] is now a-s (start/end are arbitrary)
    this->edges[edge_index].start=a;
    this->edges[edge_index].end=s;
    this->nodes[a].edges.push_back(edge_index);
    this->nodes[s].edges.push_back(edge_index);
    // middle set is ok

    // Add new edge s-b
    this->edges[this->num_edges].start=s;
    this->edges[this->num_edges].end=b;
    this->edges[this->num_edges].middle_set.clear();
    // s-b as same middle set as former a-b (also same as a-s which now lives where a-b used to...)
    this->edges[this->num_edges].middle_set=this->edges[edge_index].middle_set;
    this->nodes[s].edges.push_back(this->num_edges);
    this->nodes[b].edges.push_back(this->num_edges);

    // Add new edge s-r
    this->edges[this->num_edges+1].start=s;
    this->edges[this->num_edges+1].end=r;
    this->nodes[s].edges.push_back(this->num_edges+1);
    this->nodes[r].edges.push_back(this->num_edges+1);
    // empty middle set since the root doesn't represent anything
    this->edges[this->num_edges+1].middle_set.clear();

    // Added two new edges
    this->num_edges+=2;

    // None of the new nodes, r or s correspond to an edge in the original graph
    this->nodes[r].graph_edge_start=-2; // just a way to mark it in gviz
    this->nodes[r].graph_edge_end=-2;
    this->nodes[s].graph_edge_start=-1;
    this->nodes[s].graph_edge_end=-1;

    this->root_node=r;

    // Now traverse the tree and re-order the edges as necessary...
    vector<bool> visited(this->num_nodes,false);
    int next_node, nbr, num_visited=0;
    list<int> stack;
    stack.push_back(r);

    while(num_visited<this->num_nodes)
    {
        print(1,stack);
        next_node=stack.front();
        print_message(1,"next node=%d\n",next_node);
        stack.pop_front();
        if(!visited[next_node])
        {
            visited[next_node]=true;
            num_visited++;

            // Process this node
            for(list<int>::iterator ii=this->nodes[next_node].edges.begin();ii!=this->nodes[next_node].edges.end();++ii)
            {
                // Set the edge endpoints so that the end is nearer the root, effectively forcing all edges to point 
                // towards the root
                // Find the nbr
                if(this->edges[*ii].end==next_node)
                    nbr=this->edges[*ii].start;
                else
                    nbr=this->edges[*ii].end;

                // If we haven't visited the nbr yet, then it is farther away, and next_node should be the end of the edge
                if(!visited[nbr])
                {
                    this->edges[*ii].start=nbr;
                    this->edges[*ii].end=next_node;
                    stack.push_back(nbr);
                }
            }
        }
    }

    return r;
}

/**
* Assuming the tree is rooted, this fills the subtree list with the indices of the
* nodes beneath u in the tree.
*/
void BDTree::find_subtree(int t, list<int> *subtree)
{
    if(!this->is_rooted)
        fatal_error("%s: requires rooted tree\n");

    subtree->clear();
    // If u is a leaf node, then the subtree is just u
    if(this->nodes[t].edges.size()==1)
    {
        subtree->push_back(t);
        return;
    }

    int u;
    list<int> stack;
    list<int>::iterator ii;
    stack.push_back(t);
    while(!stack.empty())
    {
        // Grab the first thing off the stack
        u=stack.front();
        stack.pop_front();
        // u is in the subtree
        subtree->push_back(u);
        // If u is a leaf, then nothing else to do. Otherwise, there are two edges that end
        // at u, say v->u and w->u - need to add them to the stack
        for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
        {
            if(this->edges[*ii].end==u)
                stack.push_back(this->edges[*ii].start);
        }
    }
    return;
}


/**
* Assumes that you want to split u with A and B being edge lists on opposite sides
* of the split.  Letting t=|nbrs(u)|-|A|-|B|, the function exhausts all 2^t possible 
* ways of adding the remaining edges into the split.  Returns the size of the lowest
* middle set found, and fills the list C with the edge indices of the corresponding
* set.
*/
int BDTree::best_partition(int u, list<int> *A, list<int> *B, list<int> *C)
{

    list<int>::iterator ii;
    list<int> A_union_B;
    A_union_B=*A;
    for(ii=B->begin();ii!=B->end();++ii)
        A_union_B.push_back(*ii);
    A_union_B.sort();
    A->sort();
    B->sort();
    this->nodes[u].edges.sort();

    // Find the missing edges not in A nor B
    vector<int> missing_edges(this->nodes[u].edges.size(),-1);
    vector<int>::iterator vv;
    vv=set_difference(this->nodes[u].edges.begin(),this->nodes[u].edges.end(),
        A_union_B.begin(),A_union_B.end(),missing_edges.begin());

    int num_missing=vv-missing_edges.begin();
    print_message(1,"%d edges must be distributed between A and B\n",num_missing);

    int i,j, k, min_middle_set=99999,best_i=-1;
    for(i=0;i<(1<<num_missing);i++)
    {
        C->clear();
        *C=*A;
        for(j=0;j<num_missing;j++)
        {
            // Check for bit j
            if(i & (1<<j))
                // Add missing_edges[j] to the partition
                C->push_back(missing_edges[j]);
        }
        if( (k=this->split_size(u,C,EDGE_SPLIT))<min_middle_set)
        {
            min_middle_set=k;
            best_i=i;
        }
    }
    // Now reconstruct the best partition since we know i
    C->clear();
    for(j=0;j<num_missing;j++)
    {
        // Check for bit j
        if(best_i & (1<<j))
            // Add missing_edges[j] to C
            C->push_back(missing_edges[j]);
    }

    return min_middle_set;
}


/**
* Removes all pairs of vertices u,v from the graph and returns true
* if the removal of such a pair results in 2 or more components, and false 
* otherwise.
*/
bool BDTree::brute_force_two_sep()
{
    int u,v,w,u_label,v_label,num_remaining=this->G->get_num_nodes()-2;
    bool has_2sep=false;
    list<int> members;
    list<int> u_nbrs, v_nbrs;
    list<int>::iterator ii;

    for(u=0;u<this->G->get_num_nodes();u++)
    {
        print_message(1,"%s:  %d of %d\n",__FUNCTION__,u,this->G->get_num_nodes());
        u_label=this->G->nodes[u].label;
        u_nbrs=this->G->nodes[u].nbrs;
        for(ii=u_nbrs.begin();ii!=u_nbrs.end();++ii)
            this->G->remove_edge(u,*ii);
        this->G->nodes[u].label=GD_UNDEFINED;
        for(v=u+1;v<this->G->get_num_nodes();v++)
        {
            v_label=this->G->nodes[v].label;
            v_nbrs=this->G->nodes[v].nbrs;
            for(ii=v_nbrs.begin();ii!=v_nbrs.end();++ii)
                this->G->remove_edge(v,*ii);
            this->G->nodes[v].label=GD_UNDEFINED;

            w=0;
            while(w==u || w==v || this->G->nodes[w].label==GD_UNDEFINED)
                w++;

            members.clear();
            int comp_size=this->G->find_component(w,&members);
            if(comp_size!=num_remaining)
            {
                print_message(0,"Removing (%d,%d) disconnects the graph (%d!=%d)--used w=%d\n",
                    u,v,comp_size,num_remaining,w);
                has_2sep=true;
                members.sort();
                for(list<int>::iterator ii = members.begin(); ii != members.end(); ii++) 
                    printf("%d ", *ii);
                printf("\n"); 
                fflush(stdout);
            }

            this->G->nodes[v].label=v_label;
            for(ii=v_nbrs.begin();ii!=v_nbrs.end();++ii)
                this->G->add_edge(v,*ii);
        }
        this->G->nodes[u].label=u_label;
        for(ii=u_nbrs.begin();ii!=u_nbrs.end();++ii)
            this->G->add_edge(u,*ii);

    }

    return has_2sep;

}

bool BDTree::brute_force_three_sep()
{
    int u,v,w,x,u_label,v_label,w_label, num_remaining=this->G->get_num_nodes()-3;
    bool has_3sep=false;
    list<int> members;
    list<int> u_nbrs, v_nbrs, w_nbrs;
    list<int>::iterator ii;

    for(u=0;u<this->G->get_num_nodes();u++)
    {
        print_message(0,"%s:  %d of %d\n",__FUNCTION__,u,this->G->get_num_nodes());
        u_label=this->G->nodes[u].label;
        u_nbrs=this->G->nodes[u].nbrs;
        for(ii=u_nbrs.begin();ii!=u_nbrs.end();++ii)
            this->G->remove_edge(u,*ii);
        this->G->nodes[u].label=GD_UNDEFINED;
        for(v=u+1;v<this->G->get_num_nodes();v++)
        {
            v_label=this->G->nodes[v].label;
            v_nbrs=this->G->nodes[v].nbrs;
            for(ii=v_nbrs.begin();ii!=v_nbrs.end();++ii)
                this->G->remove_edge(v,*ii);
            this->G->nodes[v].label=GD_UNDEFINED;
            for(w=v+1;w<this->G->get_num_nodes();w++)
            {
                w_label=this->G->nodes[w].label;
                w_nbrs=this->G->nodes[w].nbrs;
                for(ii=w_nbrs.begin();ii!=w_nbrs.end();++ii)
                    this->G->remove_edge(w,*ii);
                this->G->nodes[w].label=GD_UNDEFINED;

                x=0;
                while(x==u || x==v || x==w || this->G->nodes[x].label==GD_UNDEFINED)
                    x++;

                members.clear();
                int comp_size=this->G->find_component(x,&members);
                if(comp_size!=num_remaining)
                {
                    print_message(0,"Removing (%d,%d,%d) disconnects the graph (%d!=%d)--used w=%d\n",
                        u,v,w,comp_size,num_remaining,x);
                    has_3sep=true;
                    members.sort();
                    for(list<int>::iterator ii = members.begin(); ii != members.end(); ii++) 
                        printf("%d ", *ii);
                    printf("\n"); 
                    fflush(stdout);
                }

                this->G->nodes[w].label=w_label;
                for(ii=w_nbrs.begin();ii!=w_nbrs.end();++ii)
                    this->G->add_edge(w,*ii);
            }

            this->G->nodes[v].label=v_label;
            for(ii=v_nbrs.begin();ii!=v_nbrs.end();++ii)
                this->G->add_edge(v,*ii);
        }
        this->G->nodes[u].label=u_label;
        for(ii=u_nbrs.begin();ii!=u_nbrs.end();++ii)
            this->G->add_edge(u,*ii);

    }

    return has_3sep;

}


/**
* Manually counts the # of neighboring BDTreeNodes of u that are leaves.
*/
int BDTree::calculate_num_leaf_nbrs(int u)
{
    int nbr;
    list<int>::iterator ii;

    this->nodes[u].num_leaf_nbrs=0;
    for(ii=this->nodes[u].edges.begin();ii!=this->nodes[u].edges.end();++ii)
    {
        nbr=-1;
        if(this->edges[*ii].end==u)
            nbr=this->edges[*ii].start;
        if(this->edges[*ii].start==u)
            nbr=this->edges[*ii].end;
            if(nbr==-1)
                print_message(0,"\n\n\nFERROR! Neighboring edge does not begin or end at u=%d??\n",u);

            // Check to see if nodes[nbr] represents a leaf
        if(this->nodes[nbr].graph_edge_start!=GD_UNDEFINED && this->nodes[nbr].graph_edge_end!=GD_UNDEFINED)
            this->nodes[u].num_leaf_nbrs++;
    }    

    return this->nodes[u].num_leaf_nbrs;

}

bool BDTree::read_bwidth_file(const char *infile)
{
    FILE *in;
    if( (in=fopen(infile,"r"))==NULL)
        fatal_error("Can't open %s for reading\n",infile);

    int i,j,k,m,n;
    //BDS - Nov 2 - compiler was complaining about your fscanf for this.
    float weight;

    // Read in nnodes and nedges from file
    fscanf(in,"%d %d",&n,&m);

    // Read in edge list
    vector<int> edge_list(2*m,0);
    for(i=0;i<m;i++)
    {
        // read in i j weight triples - not using weight right now
        fscanf(in,"%d %d %f\n",&(edge_list[2*i]), &(edge_list[2*i+1]), &weight);
    }

    fscanf(in,"%d",&j);
    // j is # nodes in tree
    if(j!=2*m-2)
        fatal_error("Not the right # of nodes in BDTree?\n");

    // Next line is -1 -1 - represents root node that is interior
    i=0;
    this->nodes[i].graph_edge_end=-1; 
    this->nodes[i].graph_edge_start=-1;
    fscanf(in,"%d %d",&j,&k);
    if(j!=-1 || k!=-1)
        fatal_error("Didn't expect this line\n");
    
    i=1;
    while(!feof(in))
    {
        // Read in details of tree node i
        fscanf(in,"%d %d\n",&j,&k);

        // We have a BDTreeEdge between BDTreeNodes i and j
        this->edges[i-1].start=i; this->edges[i-1].end=j;
        // Both i and j are adjacent to this edge
        this->nodes[i].edges.push_back(i-1);
        this->nodes[j].edges.push_back(i-1);

        // If k>=0, tree node i is a leaf and represents
        // the k'th edge in the edge_list
        if(k>=0)
        {
            this->nodes[i].graph_edge_start=edge_list[2*k];
            this->nodes[i].graph_edge_end=edge_list[2*k+1];

            // We know the middle set of the edge we just added
            this->edges[i-1].middle_set.clear();
            this->edges[i-1].middle_set.push_back(edge_list[2*k]);
            this->edges[i-1].middle_set.push_back(edge_list[2*k+1]);

        }
        else
        {
            // else it is interior and we don't know anything yet
            this->edges[i-1].middle_set.clear();
            this->nodes[i].graph_edge_start = -1;
            this->nodes[i].graph_edge_end = -1;
        }

        i++;
    }
    fclose(in);
    // Make sure the tree makes sense!
    if(i!=2*m-2)
    {
        print_message(0,"%s:  Didn't read the expected # of nodes?? %d != %d\n",__FUNCTION__,
            i,2*m-2);
        return false;
    }

    this->num_edges=2*m-3;
    this->num_nodes=2*m-2;
    this->num_interior_nodes=m-2;

    this->is_valid=true;
    this->is_rooted=true;
    return true;
}
