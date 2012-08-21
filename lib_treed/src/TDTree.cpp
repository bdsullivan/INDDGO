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

#include "WeightedMutableGraph.h"
#include "Debug.h"
#include "TreeDecomposition.h"
#include "GraphDecomposition.h"
#include "Util.h"
#include "RndNumGen.h"
#include <math.h>
#include <algorithm>

/**
 * The constructor for the TDTree class.
 */
TDTree::TDTree(Graph::WeightedMutableGraph *H)
{
	this->G=H;
	// This was new list<int>(this->G->get_capacity()); ??
	this->node_locations = new list<int>[this->G->get_capacity()];
	string ss;
	ss=H->get_input_file();
	if(!ss.empty())
	{
		this->graph_file=new char[100];
		sprintf(this->graph_file,"%s", H->get_input_file().c_str());
	}
	else
		this->graph_file = NULL;

	this->num_tree_nodes=0;
	this->num_tree_nodes_processed=0;
	this->num_leafs=0;
	this->nice=false;
	this->rooted=false;
	this->num_mask_words=0;
	this->width =  0;
	return;
}

/**
 * Default constructor for TDTree.
 */
TDTree::TDTree()
{
	this->graph_file=NULL;
	this->node_locations=NULL;
	this->num_tree_nodes=0;
	this->num_tree_nodes_processed=0;
	this->num_leafs=0;
	this->nice=false;
	this->rooted=false;
	this->num_mask_words=0;
	this->width = 0;
	return;
}

/**
 * Destructor for TDTree class.
 */
TDTree::~TDTree()
{
	int i;
	int k=(int)this->tree_nodes.size();

	for(i=0;i<k;i++)
	{
		if(this->tree_nodes[i])
		{
			delete this->tree_nodes[i];
			this->tree_nodes[i]=NULL;
		}        
	}

	if(this->graph_file)
	{
		delete [] this->graph_file;
		this->graph_file=NULL;
	}

	if(this->node_locations)
	{
		delete [] this->node_locations;
		this->node_locations=NULL;
	}

	return;
}

/** 
 * The copy constructor to do something like TDTree T=S;  This is a "deep copy".
 */
TDTree::TDTree(const TDTree& rhs)
{

	// Points to the same graph
	this->G=rhs.G;
	// Set these fields to be the same
	this->nice=rhs.nice;
	this->num_tree_nodes=rhs.num_tree_nodes;
	this->num_leafs=rhs.num_leafs;
	this->rooted=rhs.rooted;
	this->root_node=rhs.root_node;
	int cap=this->G->get_capacity();
	this->node_locations = new list<int>[cap];
	for(int i=0;i<cap;i++)
	{
		for(list<int>::const_iterator ii=rhs.node_locations[i].begin();ii!=rhs.node_locations[i].end();++ii)
			this->node_locations[i].push_back(*ii);
	}
	int size=rhs.tree_nodes.size();

	this->tree_nodes.resize(size);
	// BDS - March 28, 2011 - corrected to allow for gaps in tree node array
	for(int i=0;i<size;i++)
	{
		if(rhs.tree_nodes[i] != NULL)
		{
			this->tree_nodes[i]=new TDTreeNode();
			*(this->tree_nodes[i]) = *(rhs.tree_nodes[i]);
		}
		else
			this->tree_nodes[i] = NULL;
	}

	if(rhs.graph_file)
	{
		this->graph_file=new char[100];
		strcpy(this->graph_file,(const char *)rhs.graph_file);
	}
	else
		this->graph_file=NULL;

}

/**
 * Assignment operator for TDTree.  Use is something like 
 * TDTree new_tree; 
 * new_tree=old_tree;
 * This will make new_tree into a deep_copy of old_tree.
 */
TDTree& TDTree::operator=(const TDTree& rhs)
{
	if (this != &rhs)
	{
		// Points to the same graph
		this->G=rhs.G;
		// Set these fields to be the same
		this->nice=rhs.nice;
		this->num_tree_nodes=rhs.num_tree_nodes;
		this->num_leafs=rhs.num_leafs;
		this->rooted=rhs.rooted;
		this->root_node=rhs.root_node;
		int cap=this->G->get_capacity();
		this->node_locations = new list<int>[cap];
		for(int i=0;i<cap;i++)
		{
			for(list<int>::const_iterator ii=rhs.node_locations[i].begin();ii!=rhs.node_locations[i].end();++ii)
				this->node_locations[i].push_back(*ii);
		}
		int size=rhs.tree_nodes.size();
		this->tree_nodes.resize(size);
		//Feb 24 2011 allowing for NULL in array of tree_nodes
		for(int i=0;i<size;i++)
		{
			if(rhs.tree_nodes[i] != NULL)
			{
				this->tree_nodes[i]=new TDTreeNode();
				*(this->tree_nodes[i]) = *(rhs.tree_nodes[i]);
			}
			else
				this->tree_nodes[i] = NULL;
		}

		if(rhs.graph_file)
		{
			this->graph_file=new char[100];
			strcpy(this->graph_file,(const char *)rhs.graph_file);
		}
		else
			this->graph_file=NULL;

	}
	return *this;

}
/*
 * Overloading the << operator for TDTree.
 */
ostream &operator<<(ostream &output, const TDTree &T)
{
	int i;

	if(T.graph_file)
		output<<"Tree decomposition for graph from file "<<T.graph_file<<endl;
	output<<"Tree has "<<T.num_tree_nodes<<" nodes\n";
	for(i=0;i<(int)T.tree_nodes.size();i++)
	{
		if(T.tree_nodes[i])
			output<<*(T.tree_nodes[i]);
	}

	return output;
}


/**
 * Merges any parent(P) /child(C) TDTreeNode pair 
 * which have identical bags (P->bag has same members as C->bag) so that P now has 
 * the children of C added to its children list and C is removed from the decomposition.
 * Deletes nodes that are merged and leaves null pointers in tree_nodes array. 
 */
void TDTree::remove_duplicate_bags()
{
	print_message(1, "In remove_duplicate_bags()\n");

	int p;
	list<int>::iterator cit, git, dit; 
	TDTreeNode *parent, *child; 
	vector<int>::iterator vit; 
	vector<int> diff(2*(this->width));
	list<int> dead_nodes;

	for(p = 0; p < this->num_tree_nodes; p++)
	{
		//check to make sure we haven't already merged this "parent"
		if(find (dead_nodes.begin(), dead_nodes.end(), p) ==  dead_nodes.end())
		{
			print_message(1,"Parent %d\n", p); 
			parent = this->tree_nodes[p];
			cit = parent->adj.begin(); 
			cit++;//skip the parent
			while( cit != parent->adj.end())
			{
				child = this->tree_nodes[*cit];
				//necessary but not sufficient condition for equality
				if(child->bag.size() == parent->bag.size())
				{
					print_message(1,"Child %d has same bag size\n", *cit);
					vit = set_symmetric_difference(parent->bag.begin(), 
                                                   parent->bag.end(), 
                                                   child->bag.begin(), 
                                                   child->bag.end(),
                                                   diff.begin());
					//actual equality
					if((vit - diff.begin()) == 0)
					{
						print_message(1,"Child %d is a match!\n", child->id);
						dead_nodes.push_back(child->id);
						//add the child's children to the parent's adjacency lis
						git = child->adj.begin(); 
						git++; //skip the parent
						while(git != child->adj.end())
						{			  
							print_message(1, "Updating node %d to have parent %d\n", *git, p);
							parent->adj.push_back(*git);
							(this->tree_nodes[*git])->adj.pop_front(); 
							(this->tree_nodes[*git])->adj.push_front(p); 
							++git;
						}
						print_message(1,"Removing child\n");
						//remove child from parent's adjlist
						cit = parent->adj.erase(cit);
					}//end of merging
					else
						cit++;
				}//end of basic size test
				else
					cit++;
			}//end of loop over children
		}//test for dead parents
	}//end of topdown walk

	print_message(1, "Removing merged nodes from tree_nodes vector.\n");
	//remove the merged nodes and adjust the num_tree_nodes
	this->num_tree_nodes -= dead_nodes.size();
	for( dit = dead_nodes.begin(); dit != dead_nodes.end(); ++dit)
	{	 
		print_message(1, "Removing node %d.\n", *dit);
		delete this->tree_nodes[*dit]; 
		this->tree_nodes[*dit] = NULL;
	}
}

/**
 * Working from the top down, merges any parent(P) /child(C) TDTreeNode pair 
 * where one is a subset of the other(P->bag \subset C->bag or P->bag \superset C->bag).
 * Note this also removes duplicate bags.
 * Deletes nodes that are merged and leaves null pointers in tree_nodes array. 
 */
void TDTree::remove_subset_bags()
{
	print_message(1, "In remove_subset_bags()\n");

	int p;
	list<int>::iterator cit, git, dit; 
	TDTreeNode *parent, *child; 
	vector<int>::iterator vit; 
	vector<int> diff(2*(this->width));
	list<int> dead_nodes;

	//get a pre-order walk of the tree
	vector<int> topdown(this->num_tree_nodes);
	this->pre_order_walk(&topdown);
	int pind;

	for(pind = 0; pind < (int) topdown.size(); pind++)
	{
		//get the pth entry in the topdown walk
		p = topdown[pind];
		//check to make sure we haven't already merged this "parent"
		if(find (dead_nodes.begin(), dead_nodes.end(), p) ==  dead_nodes.end())
		{
			print_message(1,"Parent %d\n", p); 
			parent = this->tree_nodes[p];
			cit = parent->adj.begin(); 
			cit++;//skip the parent
			while( cit != parent->adj.end())
			{
				child = this->tree_nodes[*cit];

				//child is smaller
				if(!(child->bag.size() > parent->bag.size()))
				{
					print_message(1,"Child %d has equal or smaller bag size\n", *cit);
					//containment
					if(includes(parent->bag.begin(), parent->bag.end(), child->bag.begin(), child->bag.end()))
					{
						print_message(1,"Child %d is a subset!\n", child->id);
						dead_nodes.push_back(child->id);
						//add the child's children to the parent's adjacency lis
						git = child->adj.begin(); 
						git++; //skip the parent
						while(git != child->adj.end())
						{			  
							print_message(1, "Updating node %d to have parent %d\n", *git, p);
							parent->adj.push_back(*git);
							(this->tree_nodes[*git])->adj.pop_front(); 
							(this->tree_nodes[*git])->adj.push_front(p); 
							++git;
						}
						print_message(1,"Removing child\n");
						//remove child from parent's adjlist
						cit = parent->adj.erase(cit);
					}//end of merging
					else
						cit++;
				}//end of child smaller
				else //if( child->bag.size() > parent->bag.size())
				{
					print_message(1,"Child %d has larger bag size\n", *cit);
					//containment
					if(includes(child->bag.begin(), child->bag.end(), parent->bag.begin(), parent->bag.end()))
					{
						print_message(1,"Child %d is a superset!\n", child->id);
						dead_nodes.push_back(child->id);
						//add the child's children to the parent's adjacency list
						git = child->adj.begin(); 
						git++; //skip the parent
						while(git != child->adj.end())
						{			  
							print_message(1, "Updating node %d to have parent %d\n", *git, p);
							parent->adj.push_back(*git);
							(this->tree_nodes[*git])->adj.pop_front(); 
							(this->tree_nodes[*git])->adj.push_front(p); 
							++git;
						}
						print_message(1,"Removing child\n");
						//remove child from parent's adjlist
						cit = parent->adj.erase(cit);

						//update the parent's bag
						parent->bag.clear();
						for(list<int>::iterator myit = child->bag.begin(); myit != child->bag.end(); myit++)
						{
							parent->bag.push_back(*myit);
						}	
						cit++;
					}//end of merging
					else
						cit++;

				}
			}//end of loop over children
		}//test for dead parents
	}//end of topdown walk

	print_message(1, "Removing merged nodes from tree_nodes vector.\n");
	//remove the merged nodes and adjust num_tree_nodes
	this->num_tree_nodes -= dead_nodes.size();
	for( dit = dead_nodes.begin(); dit != dead_nodes.end(); ++dit)
	{	 
		print_message(1, "Removing node %d.\n", *dit);
		delete this->tree_nodes[*dit]; 
		this->tree_nodes[*dit] = NULL;
	}
}



/**
 * Checks to ensure the current TDTree is a valid rooted decomposition, 
 * then transforms it into a nice decomposition with the same width 
 * and num_tree_nodes + Sum_{nodes v with >=2 children} [2(deg(v)-1)] + 
 * Sum_{nodes v with a single child w} [|B(v)\B(w)| + |B(w)\B(v)|]  
 * tree nodes. The root remains the same.
 */
void TDTree::make_nice()
{
	//  Assumes the TD is valid, so don't check it here
	this->compute_width();// this allows us to be efficient in creating vectors to hold bags

	// Now pre-process to remove any duplicate parent-child bags, as these could 
	// cause resulting tree to fail nice criterion
	print_message(1, "Before removing duplicate bags, number of nodes is %d\n", this->num_tree_nodes);
	this->remove_duplicate_bags();
	print_message(1, "After removing duplicate bags, number of nodes is %d\n", this->num_tree_nodes);

	// First we need to allocate enough space for all the new nodes. 
	// To do this, we will need to loop over the current nodes and compute how many nodes we need.
	int i;
	TDTreeNode *curr_node, *par_node;
	int new_num_nodes = this->tree_nodes.size();
	// in case we removed duplicate bags and have nulls in array
	int old_num_nodes = this->num_tree_nodes;
	int max_tnode_id = new_num_nodes-1;
	vector<int> v(this->G->get_num_nodes());
	vector<int>::iterator vit;
	int_int logans;
	for(i = 0; i < max_tnode_id+1; i++)// was new_num_nodes 
	{
		print_message(1,"i=%d of %d\n",i,max_tnode_id);
		curr_node = this->tree_nodes[i];
		if(curr_node != NULL)
		{
			if(curr_node->adj.front() != curr_node->id) // not root
			{
				// we add the size of the symm diff with the parent - 1
				par_node = this->tree_nodes[*(curr_node->adj.begin())]; 
				vit = set_symmetric_difference(curr_node->bag.begin(), 
                                               curr_node->bag.end(), 
                                               par_node->bag.begin(), 
                                               par_node->bag.end(), 
                                               v.begin()); 
				new_num_nodes += (int)(vit - v.begin())-1 ; 
				print_message(1, "Adding %d for node %d (chain)\n", (int)(vit-v.begin())-1, curr_node->id);
			}
			if(curr_node->adj.size() >2 ) //  parent, 2+ children - joins
			{
				// number of nodes in binary tree with adj.size()-1 leaves, not counting root (it is current node)
				// new_num_nodes += 2* (int)pow(2.0, (int)ceil(log((double)curr_node->adj.size()-1)/log(2.0))) -2 + (curr_node->adj.size()-1);
				logans = mylog2(curr_node->adj.size() -1);
				if(logans.p2 > 1)
					logans.p1++; 
				new_num_nodes += (1<<(logans.p1+1)) -2 + (curr_node->adj.size()-1);
				print_message(1, "Adding %d == %d  for node %d (join)\n", 
                              2*(int)pow(2.0, (int)ceil(log((double)curr_node->adj.size()-1)/log(2.0))) -2 + (curr_node->adj.size()-1), 
                              (1<<(logans.p1+1)) -2 + (curr_node->adj.size()-1),
                              curr_node->id);
			}
		}
	}    // Done with accounting

	print_message(1, "Max number of nodes in new TD: %d\n", new_num_nodes);

	// get a pre-order walk on the original tree.
	vector<int> topdown(old_num_nodes);
	this->pre_order_walk(&topdown);
	print_message(1, "Original number of nodes: %d\n" , old_num_nodes);
	print(1, topdown);

	// Now we resize the tree_node array. The nodes will be allocated as they're created.
	print_message(1, "current size of tree_nodes: %d\n", this->tree_nodes.size());
	this->tree_nodes.resize(new_num_nodes,0);
	print_message(1, "Resized tree_nodes to %d\n", new_num_nodes );
	int tmp;
	// Create appropriate join nodes for those with multiple children, working with order from pre-order walk
	for(i = 0; i < old_num_nodes; i++)
	{
		if((this->tree_nodes[topdown[i]])->adj.size() > 2)
		{
			tmp = this->replace_by_join_tree(topdown[i], max_tnode_id);  
			this->num_tree_nodes += (tmp  - max_tnode_id);
			max_tnode_id = tmp;
		}
	}

	// get a pre-order walk on the new tree.
	old_num_nodes = this->num_tree_nodes;
	topdown.resize(old_num_nodes);
	this->pre_order_walk(&topdown);
	print_message(1, "Post-joins number of nodes: %d\n", old_num_nodes);
	print(1, topdown);

	// Create chains of nodes to address symmetric differences between nodes and their single children.
	for(i = 0; i < old_num_nodes; i++)
	{
		if((this->tree_nodes[topdown[i]])->adj.size() == 2)// one parent, one child
		{		    
			// replace with a chain
			tmp = this->expand_symdiff(topdown[i], max_tnode_id);
			// update the number of nodes in the tree
			this->num_tree_nodes += (tmp  - max_tnode_id);
			max_tnode_id = tmp;
			print_message(1, "number of nodes is now %d\n", this->num_tree_nodes);
		}	    
	}

	this->num_tree_nodes=0;
	int size=this->tree_nodes.size();
	for(i=0;i<size;i++)
	{
		if(this->tree_nodes[i] && this->tree_nodes[i]->adj.size()==1)
			this->num_leafs++;
		if(this->tree_nodes[i])
			this->num_tree_nodes++;
	}
	print_message(1,"num_tree_nodes after recalculating=%d\n",this->num_tree_nodes);

	this->nice = true;
}

/*
 * Replaces the node at P = tree_nodes[old_node] (which has one child C) by a chain of nodes, 
 * each of which has one child, 
 * differing from its parent only by 'forgetting' or 'introducing' a single vertex 
 * (following the rules of a nice decomposition). 
 * The original child becomes the bottom node of this new chain, which introduces
 * exactly |B(P)/B(C)| + |B(C)\B(P)| - 1 new nodes.
 * max_id should be the highest node id currently in the tree decomposition.
 * Returns the max id currently in use by a node after the chain has been added.
 */
int TDTree::expand_symdiff(int old_node, int max_id)
{
	vector<int> pdiff(this->width);
	vector<int> cdiff(this->width);
	vector<int>::iterator pit, cit,vit;
	TDTreeNode *parent, *child; 
	parent = this->tree_nodes[old_node]; 
	child = this->tree_nodes[*(++(parent->adj.begin()))];
	//find the vertices in parent & not child
	pit = set_difference(parent->bag.begin(), 
                         parent->bag.end(), 
                         child->bag.begin(), 
                         child->bag.end(), 
                         pdiff.begin()); 
	int numP = (int)(pit-pdiff.begin());

	//find the vertices in child and not parent
	cit = set_difference(child->bag.begin(), 
                         child->bag.end(),
                         parent->bag.begin(), 
                         parent->bag.end(), 
                         cdiff.begin()); 
	int numC = (int)(cit-cdiff.begin());

	int num_added = 0;
	int curr_id = max_id+1;
	TDTreeNode *new_node;

	print_message(1, "P = %d, C = %d, %d in P-C, %d in C-P\n", parent->id, child->id, numP, numC);

	if(numP+numC==1)
		return max_id; //we don't need to make a change here.

	parent->adj.remove(child->id); // remove the original child from the parent. 

	//create new nodes starting at the parent, deleting one vertex in pit from the bag at each step
	vit = pdiff.begin();
	while(num_added < numP)
	{
		//Feb 24 allocating only the number we need.
		this->tree_nodes[curr_id] = new TDTreeNode(); 
		this->tree_nodes[curr_id]->id=curr_id;
		new_node = this->tree_nodes[curr_id]; 
		new_node->bag = parent->bag;  //copy the bag of parent
		new_node->bag.remove(*vit); // remove an element that was in original parent and not child
		new_node->adj.push_front(parent->id); //set the parent
		parent->adj.push_back(curr_id); //set as a child of previous node
		parent = new_node; //advance the parent pointer
		curr_id++; 
		num_added++;
		++vit; 
	}

	print_message(1, "Added %d nodes to account for those in parent, not child.\n", num_added);

	//now we create new nodes that add the vertices which were in the original child
	vit = cdiff.begin();
	num_added = 0;
	while(num_added < numC-1)
	{
		//Feb 24 - allocating as we need nodes
		this->tree_nodes[curr_id] = new TDTreeNode(); 
		this->tree_nodes[curr_id]->id=curr_id;
		new_node = this->tree_nodes[curr_id]; 
		new_node->bag = parent->bag;  //copy the bag of parent
		new_node->bag.push_back(*vit); // add an element that was in original child and not parent
		new_node->bag.sort();//bags must be sorted
		new_node->adj.push_front(parent->id); //set the parent
		parent->adj.push_back(curr_id); //set as a child of previous node
		parent = new_node; //advance the parent pointer
		curr_id++; 
		num_added++;
		++vit; 
	}
	print_message(1, "Added %d nodes to account for those in child, not parent.\n", num_added);

	//the last node we added now needs to be the parent of the original child node.
	parent->adj.push_back(child->id);
	child->adj.pop_front(); //remove original parent reference
	child->adj.push_front(parent->id); //update with new parent

	print_message(1, "Replacing node %d with a chain of %d nodes. Max id now %d\n", old_node, numC+numP-1, curr_id-1);

	return curr_id-1;	
}


/*
 * Replaces the node at tree_nodes[old_node] by a binary tree with deg(old_node) leaves, each of which
 * has one of old_node's children as its only child.
 * Currently this creates a binary tree that is as balanced as possible, but this can be modified in future
 * versions of the function. 
 * Returns the max id currently in use by a node after the join tree was created & introduced.
 * This should be max_id + 2*deg(old_node)-1.
 */
int TDTree::replace_by_join_tree(int old_node, int max_id)
{
	TDTreeNode *curr_node = this->tree_nodes[old_node];
	TDTreeNode *new_node1, *new_node2;
	int curr_max_id = max_id;
	int num_leaves = curr_node->adj.size() -1; //The number of children 
	int num_levels;
	int max;
	int_int log_ans = mylog2(num_leaves);
	num_levels = log_ans.p1 + 1; 
	if(log_ans.p2 > 1)
		num_levels++; 

	int i;
	int level = 0;
	// Create the binary tree on nodes oldnode along with new nodes max_id + 1 through 
	// max_id + 2*num_leaves-1.

	//the top node is oldnode - we need to remember its children to reattach them at the bottom.
	list<int> old_children = curr_node->adj; //copy the list
	int old_parent = old_children.front();
	old_children.pop_front(); //remove the old parent
	//clear the old node's children 
	curr_node->adj.clear();
	curr_node->adj.push_front(old_parent);
	level++; //old_node is level 0

	print_message(1, "Modified old node as parent\n");

	//create two new children - we always have to do at least this first level.
	this->tree_nodes[curr_max_id+1] = new TDTreeNode(); 
	this->tree_nodes[curr_max_id+1]->id=curr_max_id+1;
	this->tree_nodes[curr_max_id+2] = new TDTreeNode(); 
	this->tree_nodes[curr_max_id+2]->id=curr_max_id+2;
	new_node1 = this->tree_nodes[curr_max_id+1];
	new_node2 = this->tree_nodes[curr_max_id+2];
	curr_max_id +=2;
	new_node1->bag = curr_node->bag; 
	new_node2->bag = curr_node->bag; 
	new_node1->adj.push_front(old_node);//add the parent 
	new_node2->adj.push_front(old_node);//add the parent 
	curr_node->adj.push_back(new_node1->id); //add as child to the old node
	curr_node->adj.push_back(new_node2->id); //add as child to the old node
	level++; // first two children were on level 1

	//we want to give two children to all the nodes in the previous level, until we reach the last level with the 
	// 'internal leaves'.
	while(level < num_levels-1)
	{
		print_message(1, "Starting level %d and curr_max_id = %d\n", level, curr_max_id);
		//i loops over the ids of the nodes in level above
		max = curr_max_id;
		for(i = curr_max_id - (1 << (level-1))  + 1; i < max+1; i++)
		{	    
			print_message(1, "creating children for node %d\n", i);
			curr_node = this->tree_nodes[i];
			//create a pair of children for every node
			this->tree_nodes[curr_max_id+1] = new TDTreeNode(); 
			this->tree_nodes[curr_max_id+1]->id=curr_max_id+1;
			this->tree_nodes[curr_max_id+2] = new TDTreeNode(); 
			this->tree_nodes[curr_max_id+2]->id=curr_max_id+2; 
			new_node1 = this->tree_nodes[curr_max_id+1];
			new_node2 = this->tree_nodes[curr_max_id+2];
			curr_max_id+=2;
			new_node1->bag = curr_node->bag; 
			new_node2->bag = curr_node->bag; 
			new_node1->adj.push_front(i);//add the parent 
			new_node2->adj.push_front(i);//add the parent 
			curr_node->adj.push_back(new_node1->id); //add as child to the previous level node
			curr_node->adj.push_back(new_node2->id); //add as child to the previous level node
		}	   
		level++;
	}

	//now handle creating the leaves - level == num_levels-1 (we counted from zero)
	//i loops over the ids of the nodes in level above    
	//we simultanously attach the old children to the leaves
	print_message(1, "Starting level %d and curr_max_id = %d\n", level, curr_max_id);
	i = curr_max_id - (1 << (level-1)) + 1;
	list<int>::iterator iit = old_children.begin();
	int leaf_count = 0;
	TDTreeNode *old_child1, *old_child2;
	while(leaf_count < num_leaves)
	{
		if(leaf_count < num_leaves-1)
		{	    
			curr_node = this->tree_nodes[i];
			//create a pair of children for every node 
			this->tree_nodes[curr_max_id+1] = new TDTreeNode(); 
			this->tree_nodes[curr_max_id+1]->id=curr_max_id+1;
			this->tree_nodes[curr_max_id+2] = new TDTreeNode(); 
			this->tree_nodes[curr_max_id+2]->id=curr_max_id+2;
			new_node1 = this->tree_nodes[curr_max_id+1];
			new_node2 = this->tree_nodes[curr_max_id+2];
			curr_max_id+=2;
			new_node1->bag = curr_node->bag; 
			new_node2->bag = curr_node->bag; 
			new_node1->adj.push_front(i);//add the parent 
			new_node2->adj.push_front(i);//add the parent 
			curr_node->adj.push_back(new_node1->id); //add as child to the previous level node
			curr_node->adj.push_back(new_node2->id); //add as child to the previous level node
			old_child1 = this->tree_nodes[*iit];
			iit++; 
			old_child2 = this->tree_nodes[*iit];
			iit++; 
			old_child1->adj.pop_front(); //remove the old parent
			old_child2->adj.pop_front(); //remove the old parent
			old_child1->adj.push_front(new_node1->id); //add the new parent
			old_child2->adj.push_front(new_node2->id); //add the new parent
			new_node1->adj.push_back(old_child1->id); //add them as children
			new_node2->adj.push_back(old_child2->id);
			leaf_count+=2;
		}
		else
		{
			curr_node = this->tree_nodes[i];
			//create a pair of children, but only one has an original child. 
			this->tree_nodes[curr_max_id+1] = new TDTreeNode(); 
			this->tree_nodes[curr_max_id+1]->id=curr_max_id+1;
			this->tree_nodes[curr_max_id+2] = new TDTreeNode(); 
			this->tree_nodes[curr_max_id+2]->id=curr_max_id+2;
			new_node1 = this->tree_nodes[curr_max_id+1];
			new_node2 = this->tree_nodes[curr_max_id+2];
			curr_max_id+=2;
			new_node1->bag = curr_node->bag; 
			new_node2->bag = curr_node->bag; 
			new_node1->adj.push_front(i);//add the parent 
			new_node2->adj.push_front(i);//add the parent 
			curr_node->adj.push_back(new_node1->id); //add as child to the previous level node
			curr_node->adj.push_back(new_node2->id); //add as child to the previous level node
			old_child1 = this->tree_nodes[*iit];
			iit++; 
			old_child1->adj.pop_front(); //remove the old parent
			old_child1->adj.push_front(new_node1->id); //add the new parent
			new_node1->adj.push_back(old_child1->id); //add them as children
			leaf_count+=2;
		}
		i++;
	}

	print_message(1, "Replacing node %d with a tree with %d leaves and %d levels. Max id now %d\n", old_node, num_leaves, num_levels, curr_max_id);

	return curr_max_id;

}


/** Constructs a tree decomposition of this->G using
 * Bodlaender-Koster's algorithm 2 from Treewidth Computations I, Upper Bounds.   
 * The graph G is assumed to have 
 * a single connected component since tree decompositions
 * of disconnected components can be joined together in
 * an arbitrary way.
 */
void TDTree::construct_BK(vector<int> *elim_order)
{
	// First make sure this->G has one connected component
	if(this->G->get_num_components()!=1)
	{
		fatal_error("%s:  Tree decomposition routines operate only on connected graphs."
                    "  This graph has %d connected components\n",__FUNCTION__,
                    this->G->get_num_components());
	}

	// We change the graph G, so we should use a copy. Note that this is potentially
	// expensive for large graphs
	Graph::WeightedMutableGraph H=*(this->G);

	int i, v, minpos, num_fwd_neighbors;
	list<int> neighbors;
	list<int>::iterator jj;
	int num_nodes=H.get_num_nodes();
	// We know beforehand that we will have a single TDTreeNode for every vertex in G
	this->tree_nodes.resize(num_nodes,NULL);
	this->num_tree_nodes=num_nodes;

	// Allocate room for the TDTreeNodes
	for(i=0;i<num_nodes;i++)
	{
		TDTreeNode *next_tree_node;
		next_tree_node=new TDTreeNode;
		next_tree_node->id=i;
		// tree_nodes[i] is a pointer
		this->tree_nodes[i]=next_tree_node;
	}

	// Assume that elim_order is num_nodes long 
	// We should not encounter any positions in G's node[] array
	// that are undefined as long as the ordering is valid.
	Graph::GraphEOUtil eoutil;
	Graph::Node *n1;
	for(i=0;i<num_nodes;i++)
	{
		v=elim_order->at(i);
		n1=H.get_node(v);
		if(n1->get_label()==GD_UNDEFINED)
			fatal_error("%s:  Trying to eliminate GD_UNDEFINED node %d\n",__FUNCTION__,v);

		// Find the forward neighbors of v wrt the elim_order
		num_fwd_neighbors=eoutil.find_forward_neighbors(&H, v,elim_order,i+1,&neighbors,&minpos);

		print_message(1,"%s:  i=%d of %d; minpos=%d\n",__FUNCTION__,i,num_nodes,minpos);

		this->tree_nodes[v]->bag.push_front(v);
		for(jj=neighbors.begin();jj!=neighbors.end();++jj)
			this->tree_nodes[v]->bag.push_front(*jj);
		// Sort the bag
		this->tree_nodes[v]->bag.sort();

		if(neighbors.size()>0)
		{
			// Check for loop?
			if(v==elim_order->at(minpos))
			{
				fatal_error("%s:  Trying to create loop in TDTree (i=%d, v=%d, minpos=%d)\n",__FUNCTION__,
                            i,v,elim_order->at(minpos));
			}
			// Now add the tree_edge between tree_node[v] and tree_node[minpos]
			print_message(1,"Adding 1-based tree_edge %d-%d\n",v+1,elim_order->at(minpos)+1);
			this->tree_nodes[v]->adj.push_back(elim_order->at(minpos));
			this->tree_nodes[elim_order->at(minpos)]->adj.push_back(v);
		}

		// Now run the "forward elimination" routine on v
		if(neighbors.size()>0)
		{
			print_message(1,"Eliminating vertex %d with %d fwd neighbors\n",
                          v,num_fwd_neighbors);
			H.eliminate_vertex(v,&neighbors,false);
		}
		else
			print_message(1,"i=%d of %d:  Vertex %d has no neighbors in BK - continuing...\n",
                          i,H.get_capacity(),v);
	}    
	
	for(i=0;i<num_tree_nodes;i++)
	{
		if(this->tree_nodes[i])
		{
			if(this->tree_nodes[i]->adj.size()==1)
				this->num_leafs++;
		}
	}
	
	return;
}


/** Constructs a tree decomposition of this->G using
 * Gavril's algorithm.   The graph G is assumed to have 
 * a single connected component since tree decompositions
 * of disconnected components can be joined together in
 * an arbitrary way.
 */
void TDTree::construct_gavril(vector<int> *elim_order)
{  
	// First make sure this->G has one connected component
	if(this->G->get_num_components()!=1)
	{
		fatal_error("%s:  Tree decomposition routines operate only on connected graphs."
                    "  This graph has %d connected components\n",__FUNCTION__,
                    this->G->get_num_components());
	}

	int i, k, m, v, minpos, num_fwd_neighbors;
	list<int> neighbors;
	list<int>::iterator ii;
	vector<int>::iterator vi;
	vector<int> W;
	int num_nodes=this->G->get_num_nodes();
	vector<int> t(num_nodes,-1);

	// k counts the # of nodes in the TDTree
	k=0;

	// For Gavril, we will have at most n vertices in the TDTree
	this->tree_nodes.resize(num_nodes,NULL);

	t[elim_order->at(num_nodes-1)]=k;
	// Create the first tree_node
	TDTreeNode *next_tree_node;
	next_tree_node=new TDTreeNode;
	next_tree_node->id=k;

	// Add the newly created node to the tree
	this->tree_nodes[k]=next_tree_node;
	this->tree_nodes[k]->bag.push_front(elim_order->at(num_nodes-1));

	int ctr=0;
	
	Graph::GraphEOUtil eoutil;

	for(i=num_nodes-2; i >= 0 ; i--)
	{
		ctr++;
		if(ctr%1000==0)
        {
            DEBUG("%d tree nodes done\n",ctr);
			//fprintf(stderr,"%d tree nodes done\n",ctr);p
        }

		neighbors.clear();
		v=elim_order->at(i);

		num_fwd_neighbors=eoutil.find_forward_neighbors(G, v,elim_order,i+1,&neighbors,&minpos);
		if(num_fwd_neighbors==0)
			fatal_error("%s:  No forward neighbors found for vertex %d\n"
                        "This is currently a fatal error - perhaps the graph has\n"
                        "more than one connected component that contains two or more\n"
                        "vertices???\n",__FUNCTION__,v);

		m=elim_order->at(minpos);        
		neighbors.sort();

		this->tree_nodes[t[m]]->bag.sort();

		// CSG - Below is why we have a vector of treenodes
		// This is a comparison of two lists
		if(neighbors==(this->tree_nodes[t[m]]->bag))
		{
			print_message(1,"Adding to existing bag\n");
			this->tree_nodes[t[m]]->bag.push_back(elim_order->at(i));
			t[elim_order->at(i)]=t[m];
		}
		else
		{
			print_message(1,"Creating new bag\n");
			k++;
			neighbors.push_back(elim_order->at(i));
			t[elim_order->at(i)]=k;
			// Now add a new TDTreeNode to the tree
			next_tree_node=new TDTreeNode;
			next_tree_node->id=k;

			// Now add the edge k,t[m]
			next_tree_node->adj.push_front(t[m]);
			this->tree_nodes[k]=next_tree_node;
			this->tree_nodes[t[m]]->adj.push_back(k);
			for(ii=neighbors.begin();ii!=neighbors.end();++ii)
				this->tree_nodes[k]->bag.push_back(*ii);

			// Keep the bag sorted
			this->tree_nodes[t[m]]->bag.sort();
		}
	}

	this->num_tree_nodes=k+1; // This was k!
	// CSG changing loop below upper index from: (int)this->tree_nodes.size()
	for(i=0;i<this->num_tree_nodes;i++)
	{
		if(this->tree_nodes[i] && this->tree_nodes[i]->adj.size()==1)
			this->num_leafs++;
	}

	return;
}


/** Constructs a tree decomposition of this->G from supernodal elimination tree.   
 * Uses cholmod, a part of SuiteSparse by Tim Davis. Exits with errors if not installed.
 * IMPORTANT: The graph does not need to be triangulated for this! Saves a huge amount of time and space, 
 * so please update calling routines to avoid triangulation when creating trees this way.
 * The graph G is assumed to have  a single connected component since tree decompositions
 * of disconnected components can be joined together in an arbitrary way.
 */
void TDTree::construct_superetree(vector<int> *elim_order)
{  
  if(!HAS_SUITESPARSE)
	{
		fatal_error("%s:  Superetree construction only available with SuiteSparse installed and configured in make.inc"
			    , __FUNCTION__);
	}
  
  // First make sure this->G has one connected component
	if(this->G->get_num_components()!=1)
	{
		fatal_error("%s:  Tree decomposition routines operate only on connected graphs."
                    "  This graph has %d connected components\n",__FUNCTION__,
                    this->G->get_num_components());
	}

	cholmod_common Common ;
	SuperTree *S ;
	int j; 

	//input graph in suitesparse format.
	Long n = (Long)G->get_num_nodes(); 
	       
	Graph::GraphUtil util; 
	util.populate_CRS(G);
	// Nodes are numbered 0,1,...,n-1
	// The neighbors of node i are adjncy[xadj[i]],adjncy[xadj[i]+1],...adjncy[xadj[i+1]-1]
	vector<int> xadj = G->get_xadj(); 
	vector<int> adjncy = G->get_adjncy();
	util.free_CRS(G);
	
	Long *Ap =new Long[xadj.size()];
	for(int i= 0; i < xadj.size(); i++){
	    Ap[i] = xadj.at(i);
	}
	
	Long *Ai =new Long[adjncy.size()];
	for(int i= 0; i < adjncy.size(); i++){
	  Ai[i] = adjncy.at(i);
	}

	Long *perm =new Long[G->get_num_nodes()];
	for(int i= 0; i < G->get_num_nodes(); i++){
	  perm[i] = elim_order->at(i);
	}
	
	cholmod_l_start (&Common) ;             // creates Common for CHOLMOD
	//set up the methods for superwrap.
	Common.nmethods = 1 ;
	Common.method [0].ordering = CHOLMOD_GIVEN ;

	
	S = SuperWrap (n, Ap, Ai, perm, &Common) ;  // allocates and creates S
	

	//Populate the tree decomposition from S
	this->tree_nodes.resize(S->ns,NULL);	
	for(int k = 0; k < S->ns; k++)
	  { 
	  // Create the k^th tree_node
	  TDTreeNode *next_tree_node;
	  next_tree_node=new TDTreeNode;
	  next_tree_node->id=k;	
	  this->tree_nodes[k]=next_tree_node;
	  //set up the bag of the tree_node, sorted
	  for(int l = S->Sp[k]; l < S->Sp[k+1]; l++)
	    {
	      this->tree_nodes[k]->bag.push_back(S->Si[l]);
	    }
	  this->tree_nodes[k]->bag.sort();
	}
	

	//add all the edges
	for(int i = 0; i < S->ns; i++){
	  j =  S->Sparent[i];
	  if(j == -1){
	    //this is the root, it is its own parent
	    (this->tree_nodes[i])->adj.push_front(i);
	  }
	  else
	    {
	      (this->tree_nodes[i])->adj.push_front(j);
	      (this->tree_nodes[j])->adj.push_back(i);
	    }
	}

	SuperFree (&S, &Common) ;                   // frees S
	cholmod_l_finish (&Common) ;                   // frees contents of Common
	delete[] perm;
}




/**
 * Fills up the node_locations[] array so that node_locations[i] is a list of
 * the indices of the tree nodes whose bags contain connected node i.
 */
void TDTree::record_node_locations()
{
	int i,k;
	list<int>::iterator j;
	int size=(int)this->tree_nodes.size();

	for(i=0;i<size;i++)
	{
		if(this->tree_nodes[i] != NULL)
		{
			for(j=this->tree_nodes[i]->bag.begin();
				j!=this->tree_nodes[i]->bag.end();++j)
			{
				k=*j;
				// Node k in the graph G is in the bag for tree node i
				this->node_locations[k].push_front(i);
			}
		}
	}
	return;
}

/**
 * Recursive function for depth first search for the path between two nodes in the tree.
 * The pred[] array should already be allocated and be of size tree_nodes.size().
 */
void TDTree::find_path(int start, int end, int *pred)
{
	print_message(1, "Entering find_path %d %d\n", start, end);

	int v;
	list<int>::iterator i1;

	for(i1=this->tree_nodes[start]->adj.begin();i1!=this->tree_nodes[start]->adj.end();++i1)
	{
		v = *i1;
		if(v==end)
		{
			// We reached the end node
			pred[end]=start;
			break;
		}
		if(pred[v]==-1)
		{
			// Haven't visited this node before - search deeper from here
			pred[v]=start;
			this->find_path(v,end,pred);
		}

	}
	// breaks to here
	print_message(1, "Leaving find_path %d %d\n", start, end);
	return;
}

/**
 * Function to find the path between two vertices in the tree.  If a path
 * is found, then the path list is populated with the path as
 * i, j, ..., m where i!=start and m!=end.  In other words, we only include the
 * "interior" of the path.  
 */
void TDTree::find_path(int start, int end, list<int> *path)
{
	// Create a pred[] array to trace the path
	int j,k,*pred;
	int size=this->tree_nodes.size();
	pred=new int[size];
	memset(pred,-1,size*sizeof(int));

	print_message(100,"\nBefore:\n");
	for(k=0;k<size;k++)
		print_message(100,"pred[%d]=%d\n",k,pred[k]);

	find_path(start,end,pred);

	print_message(100,"\nAfter:\n");
	for(k=0;k<size;k++)
		print_message(100,"pred[%d]=%d\n",k,pred[k]);

	// The pred array should contain the information about the path from
	// start to end
	// Now fill the path list - start and end are NOT included in the path
	k=0;
	j=pred[end];
	while(j!=start)
	{
		k++;
		path->push_front(j);
		j=pred[j];

		// Sanity check
		if(k > size)
			fatal_error("%s:  Found too many nodes in find_path from %d to %d\n",
                        __FUNCTION__,start,end);
	}

	// clean up
	delete [] pred;
	return;
}

/**
 * Determines whether or not the decomposition is valid by explicitly checking
 * the three conditions. Note that this is potentially very expensive!
 */
bool TDTree::verify()
{
	////////////////
	// Condition 1 //
	////////////////
	// First make sure every connected node in the graph is in some bag
	// Create a binary vector in_tree where in_tree[i]=0 until we find connected node
	// i in some bag
	int cap=this->G->get_capacity();
	print_message(1,"%s:  G->capacity=%d\n",__FUNCTION__,cap);
	int *in_tree;
	in_tree=new int[cap];
	memset(in_tree,0,cap*sizeof(int));

	int i, j;
	list<int>::iterator i1,i2,i3;
	int size=(int)this->tree_nodes.size();
	for(i=0;i<size;i++)
	{
		if(this->tree_nodes[i]!= NULL )
		{
			// Sort the bags - we need this for Criteria 3 anyway
			this->tree_nodes[i]->bag.sort();
			for(i1=this->tree_nodes[i]->bag.begin();
				i1!=this->tree_nodes[i]->bag.end();++i1)
			{
				in_tree[*i1]++;
				print_message(1, "Found vertex %d in tree node %d\n",*i1,i);
				print(1,(this->tree_nodes[i]->bag));
			}
		}
	}

	// Now make sure we don't have in_tree[i]==0 anywhere
	bool found_all=true;
	Graph::Node *n1;
	list<int> *nbr_list;
	for(i=0;i<cap;i++)
	{
		n1=this->G->get_node(i);
		nbr_list=n1->get_nbrs_ptr();
		if(in_tree[i]==0 && nbr_list->size()>0)
		{
			// If vertex i is disconnected from the graph, this "might" be ok for now...
			print_message(0,"Graph vertex %d not found in any tree node bags!\n",i);
			found_all=false;
		}
	}
	// Clean up
	delete [] in_tree;
	// Return false if invalid
	if(!found_all)
		return false;

	print_message(1,"Condition 1 is satisfied\n");

	////////////////
	// Condition 2 //
	////////////////
	// Now make sure every edge is in some bag
	// For each connected node in the graph, record the indices of the bags
	// that contain the node
	this->record_node_locations();
	// indicator variables for the edges
	bool found_i1=false, found_i=false;

	bool found_all_edges=true;

	for(i=0;i<cap;i++)
	{
		n1=this->G->get_node(i);
		nbr_list=n1->get_nbrs_ptr();
		// Check every edge i-*i1 where i<*i1
		for(i1=nbr_list->begin();i1!=nbr_list->end();++i1)
		{
			print_message(1,"Checking for edge %d-%d (%d edges)\n",
                          i,*i1,nbr_list->size());
			if(*i1>i)
			{
				// We are searching the bags for the edge i-*i1
				if(this->node_locations[*i1].size()<this->node_locations[i].size())
					j=*i1;
				else
					j=i;
				// j is either *i1 or i, whichever has fewer bags to try and search fewer bags

				for(i2=this->node_locations[j].begin();i2!=this->node_locations[j].end();++i2)
				{
					// *i2 will represent the index of the tree node
					// Search the bag of tree node *i2 for both *i1 and i
					// Reset the indicator variables
					found_i1=false; found_i=false;
					for(i3=this->tree_nodes[*i2]->bag.begin();
						i3!=this->tree_nodes[*i2]->bag.end();++i3)
					{
						if(*i3==*i1) found_i1=true;
						if(*i3==i) found_i=true;
						if(found_i1 && found_i)
							// We found the edge i-*i1 in the bag of tree node *i2
							break;

					}
					if(found_i1 && found_i)
						// We found the edge i-*i1 in the bag of tree node *i2
						break;
				}
				// We should have found_i1=found_i=true once we get here
				if(!found_i1 || !found_i)
				{
					print_message(0,"Did not find edge %d-%d in any bags!\n",i,*i1);
					found_all_edges=false;
					// Don't return yet since we want to see if we are missing more than one edge
					//return false;
				}  
			}
		}
	}

	if(found_all_edges)
		print_message(1,"Condition 2 is satisfied\n");

	////////////////
	// Condition 3 //
	////////////////
	// For every pair of tree nodes i and j, if bag_i contains v and bag_j contains v, then
	// every node on the path in the tree between tree nodes i and j must contain v
	// The check here is naive - surely this can be made more efficient by finding a set of longest
	// paths in the graph so that for every i,j w/ i<j we have a path that contains both
	// i and j

	list<int> path;
	vector<int> I(cap);
	vector<int> J(cap);
	vector<int>::iterator v, v2;
	int tree_size=(int)this->tree_nodes.size();
	for(i=0;i<tree_size;i++)
	{
		if(this->tree_nodes[i])
		{
			for(j=i+1;j<tree_size;j++)
			{
				if(this->tree_nodes[j])
				{
					// Consider tree nodes i and j
					// First, find the nodes in the intersection
					print_message(1,"Finding intersection of bags %d(%d nodes) and %d(%d nodes)\n",
                                  i,this->tree_nodes[i]->bag.size(),
                                  j,this->tree_nodes[j]->bag.size());

					v=set_intersection(
						this->tree_nodes[i]->bag.begin(),this->tree_nodes[i]->bag.end(),
						this->tree_nodes[j]->bag.begin(),this->tree_nodes[j]->bag.end(),
						I.begin());
					print_message(1,"Intersection has %d elements\n",(int)(v-I.begin()));
					// Only bother to find the path if the intersection is non-empty 
					// v is an iterator to the end of the intersection.
					if(I.begin()!=v)
					{
						// Intersection is non-empty - find the path
						path.clear();
						this->find_path(i,j,&path);
						// Now make sure that the bag of each tree node in the path contains I
						for(i2=path.begin();i2!=path.end();++i2)
						{
							// Check tree node *i2
							v2=set_intersection(I.begin(),v,
                                                this->tree_nodes[*i2]->bag.begin(),
                                                this->tree_nodes[*i2]->bag.end(),
                                                J.begin());
							// We must have I==J - sufficient to check lengths
							if((v2-J.begin()) != (v-I.begin()))
							{
								print_message(1,"Intersection set of size %d!=%d not found for tree node %d "
                                              "on the path from %d to %d\n",(v2-J.begin()),(v-I.begin()),
                                              *i2,i,j);
								return false;
							}
						}
					}
				}
			}
		}
	}

	print_message(1, "Condition 3 is satisfied\n");

	if(!found_all_edges)
		return false;

	// Everything passed so it must be valid!
	return true;
}

/**
 * Wrapper to write DIMACS file without considering preordering.
 */
void TDTree::write_DIMACS_file(const char *DIMACS_file)
{
    this->write_DIMACS_file(DIMACS_file, false);
}

/**
 * Writes the tree decomposition to the given file in
 * our interpretation of the DIMACS format for tree decomposition.
 * If preordered = true, finds a preorder walk and relabels tree nodes accordingly.
 * Assumes tree is rooted. Also writes edges as 'e p c' where p is the parent of c in 
 * the rooted tree.
 */
void TDTree::write_DIMACS_file(const char *DIMACS_file, bool preordered)
{
	FILE *out;
	if( (out=fopen(DIMACS_file,"w"))==NULL)
		fatal_error("%s:  Couldn't open %s for writing\n",__FUNCTION__,DIMACS_file);

	int i, j, m=0;
	list<int>::iterator L;

	// Count the # of edges in the tree decomposition

	int size=(int)this->tree_nodes.size();
	for(i=0;i<size;i++)
	{
		if(this->tree_nodes[i] != NULL)
			m+=this->tree_nodes[i]->adj.size();
	}
	m=m/2;// We double counted all the edges!

	// Print out the "problem" information
	fprintf(out,"p treed %d %d\n",this->num_tree_nodes,m);
	fflush(out);

	vector<int> walk(this->num_tree_nodes,GD_UNDEFINED);
	//if writing pre-order, find the walk and store it.
	if(preordered)
		this->pre_order_walk(&walk);
	else
	{
	    for(int i = 0; i < this->num_tree_nodes; ++i)
            walk[i] = i;
	}

	vector<int> perm(this->num_tree_nodes,GD_UNDEFINED);
	for(int i = 0; i < this->num_tree_nodes; ++i)
	    perm[walk[i]] = i;


	
	// Print out the edges so that 'e p c' means p is the parent of c	
	for(i=0;i<size;i++)
	{   
		if(this->tree_nodes[i] != NULL)
		{
		    L = this->tree_nodes[i]->adj.begin(); 
		    ++L; //skip the parent
		    for(L;L!=this->tree_nodes[i]->adj.end();++L)
			{
			    j=*L;
			    fprintf(out,"e %d %d\n",perm[i]+1,perm[j]+1);
			    fflush(out);
			}
		}
	}


	// Print out the bags
	for(i=0;i<size;i++)
	{   
		if(this->tree_nodes[i] != NULL)
		{
		    fprintf(out,"B %d %u ",perm[i]+1,(int)this->tree_nodes[i]->bag.size());
			
			L=this->tree_nodes[i]->bag.begin();
			j=*L;
			Graph::Node *n1;
			if(this->G->get_capacity()==this->G->get_num_nodes())
				j++;    // Just change from 1-based to 0-based
			else
			{
				n1=this->G->get_node(j);
				j=n1->get_label();
				
			}
			fprintf(out,"%d",j);
			++L;
			for(;L!=this->tree_nodes[i]->bag.end();++L)
			{
				// Get the actual node label if necessary (we have a pointer to G
				j=*L;
				if(this->G->get_capacity()==this->G->get_num_nodes())
					j++;    // Just change from 1-based to 0-based
				else
				{
					n1=this->G->get_node(j);
					j=n1->get_label();
					
				}

				fprintf(out," %d",j);
				// NOTE!! The bags store the DIMACS labels (1,2,..,num_total_nodes) so if there
				// are disconnected nodes, then the union of all the bags
				// will NOT be equal to (1,2,...,num_total_nodes)

			}
			fprintf(out,"\n");
		}//if not null
	}//for loop to print bags

	fflush(out);
	fclose(out);

	return;
}

/**
 * Populates the TDTree data structures by reading the provided file.  Returns
 * the size of the largest bag minus one (aka the width).
 */
int TDTree::read_DIMACS_file(const char *DIMACS_file)
{
	int width=-1;
	FILE *in;
	if( (in=fopen(DIMACS_file,"r"))==NULL)
		fatal_error("%s:  Can't open %s for reading\n",__FUNCTION__,DIMACS_file);

	int retval, i, j, k, m, n, id, num, start, end;
	char line[200], format[200];

	// CSG - adding a second pass through the file to address re-numbering issues
	// Just look at the edges and create perm and iperm
	list<int> bag_labels;
	list<int> edge_labels;
	while(!feof(in))
	{
		retval=fscanf(in,"%2c",line);
		if(feof(in))
			break;

		if(retval!=1)
			fatal_error("%s:  fscanf read %d char's (expected to read 1!)\n",__FUNCTION__,
                        retval);

		switch(line[0])
		{
            case 'p':
                // This is the problem line.  
                retval=fscanf(in,"%s %d %d",format,&n,&m);
                if(strcmp(format,"treed")!=0 && strcmp(format,"TreeD")!=0 &&
                   strcmp(format,"TREED")!=0 && strcmp(format,"Treed")!=0)
                    fatal_error("%s:  Unknown problem type of %s in TreeD::read_DIMACS_file\n",
                                __FUNCTION__);


            case 'f':
                break;

            case 'c':
                break;

            case 'B':
                // "Bag" line - of the form B id k n_1 n_2 ... n_k
                // In other words, there are k nodes in the bag
                retval=fscanf(in,"%d %d",&id,&num);
                bag_labels.push_back(id);
                break;

            case 'e':
                // Edge line - of the form e start end
                retval=fscanf(in,"%d %d",&start,&end);
                if(retval!=2)
                    fatal_error("%s:  DIMACS read error - didn't understand edge line\n",
                                __FUNCTION__);
                edge_labels.push_back(start);
                edge_labels.push_back(end);
                break;

            default:
                // We encountered some character that we don't understand
                fatal_error("%s:  DIMACS read error - didn't understand character %c to start a line\n",
                            __FUNCTION__,line[0]);
		}// end switch

		// Advance to next line if there is format at the end of it
		while (!feof(in) && (retval=getc(in)) != '\n') ;  
	}// end while

	// Make sure edge_labels and bag_labels are identical
	edge_labels.sort();
	edge_labels.unique();
	bag_labels.sort();
	if(edge_labels!=bag_labels)
		fatal_error("%s:  labels are not the same\n",__FUNCTION__);
	else
		print_message(1,"labels are identical\n");
	// Find max. label
	int max_label=-1;
	list<int>::iterator ii;
	for(ii=bag_labels.begin();ii!=bag_labels.end();++ii)
		if(*ii>max_label)
			max_label=*ii;

	// We need to re-label everything
	vector<int> perm(max_label+1,-1);
	int ctr=0;
	for(ii=bag_labels.begin();ii!=bag_labels.end();++ii)
	{
		perm[*ii]=ctr;
		ctr++;
	}

	rewind(in);

	j=0;
	// Use j to count the tree nodes read in
	while(!feof(in))
	{
		retval=fscanf(in,"%2c",line);
		if(feof(in))
			break;

		if(retval!=1)
			fatal_error("%s:  fscanf read %d char's (expected to read 1!)\n",__FUNCTION__,
                        retval);

		switch(line[0])
		{
            case 'p':
                // This is the problem line.  
                retval=fscanf(in,"%s %d %d",format,&n,&m);
                if(strcmp(format,"treed")!=0 && strcmp(format,"TreeD")!=0 &&
                   strcmp(format,"TREED")!=0 && strcmp(format,"Treed")!=0)
                    fatal_error("%s:  Unknown problem type of %s in TreeD::read_DIMACS_file\n",
                                __FUNCTION__);

                print_message(1,"Reading tree decomposition with %d nodes in tree and %d edges\n",
                              n, m);

                // Allocate n TDTreeNodes
                if((int)this->tree_nodes.size()<n)
                    this->tree_nodes.resize(n,NULL);
                for(i=0;i<n;i++)
                {
                    TDTreeNode *next_tree_node;
                    next_tree_node=new TDTreeNode;
                    next_tree_node->id=i;
                    this->tree_nodes[i]=next_tree_node;
                }
                this->num_tree_nodes=n;
                print_message(1,"Allocated for tree to read in\n");

                break;
            case 'f':
                // This is the "file line" - the next entry is the name of the file
                // containing the original graph for this decomposition
                retval=fscanf(in,"%s",this->graph_file);
                print_message(1,"Scanned in file name %s\n",this->graph_file);
                //strcpy(this->graph_file,format);
                if(retval==0)
                    fatal_error("%s:  NULL filename found in TreeD::read_DIMACS_file??\n",
                                __FUNCTION__);
                print_message(100,"DIMACS: read graph_file=%s\n",this->graph_file);

                break;

            case 'c':
                // Comment line - skip and move on
                break;

            case 'B':
                // "Bag" line - of the form B id k n_1 n_2 ... n_k
                // In other words, there are k nodes in the bag
                retval=fscanf(in,"%d %d",&id,&num);
                if(num>width)
                    width=num;
                // We read in id - stick this in position perm[id];
                print_message(1,"Read in bag id %d - translated to %d\n",id,perm[id]);
                id=perm[id];

                print_message(1,"Reading in bag for TreeNode %d (j=%d)\n\n",id,j);

                for(i=0;i<num;i++)
                {
                    fscanf(in,"%d",&k);

                    k--; // The file information is 1-based!
                    if(k<0)
                        fatal_error("%s:  Read in negative vertex number of %d?\n",__FUNCTION__,k);
                    this->tree_nodes[id]->bag.push_back(k);
                    this->tree_nodes[id]->bag_vec.push_back(k);
                }

                print_message(1,"Tree Node %d has %d vertice in its bag\n",j,
                              this->tree_nodes[id]->bag.size());

                j++;
                break;

            case 'e':
                // Edge line - of the form e start end
                retval=fscanf(in,"%d %d",&start,&end);
                if(retval!=2)
                    fatal_error("%s:  DIMACS read error - didn't understand edge line\n",
                                __FUNCTION__);
                // We need to add this information to our adjacency lists
                //start--;end--;
                start=perm[start];
                end=perm[end];
                print_message(1,"Read in tree edge %d-%d\n",start,end);
                this->tree_nodes[start]->adj.push_back(end);
                this->tree_nodes[end]->adj.push_back(start);
                // We thankfully don't have to worry about disconnected nodes since
                // we have a tree by assumption
                break;

            default:
                // We encountered some character that we don't understand
                fatal_error("%s:  DIMACS read error - didn't understand character %c to start a line\n",
                            __FUNCTION__,line[0]);
		}// end switch

		// Advance to next line if there is format at the end of it
		while (!feof(in) && (retval=getc(in)) != '\n') ;  
	}// end while

	fclose(in);

	// return the width
	width--;
	return width;

}

/**
 * Function to write a graphviz representation (DOT language) of the tree decomposition to the specified file. 
 * The spline parameter allows the user to specify whether or not spline curves are used in edge layout. 
 * This is (strongly) not recommended for large graphs (> 100 vertices).
 * The final parameter takes one of the options GV_BAG_LABELS, GV_SCALE_BAGS, GV_TREE_ONLY
 * GV_BAG_LABELS lists the entire bag of vertices inside each node's label. 
 * GV_SCALE_BAGS scales the bags to indicate bag size, and labels them with the treenode ID. Note this may 
 * create dimacs files that neato cannot handle on very large graphs! 
 * GV_TREE_ONLY writes a file that just draws the underlying tree.
 * GV_COLORS creates a colored version of the decomposition showing the underlying tree with nodes colored
 * red if they have size k+1, blue if size k, and gradations of purple for smaller nodes, getting lighter as 
 * you get smaller.
 **/
void TDTree::write_graphviz_file(bool spline, const char *GVIZ_file, int style)
{
	int i, j; 
	FILE *out; 
	Graph::WeightedMutableGraph *G = this->G;
	list<int>::iterator L;
	char rgb[8];
	int k=0, rcolor; 

	if( (out = fopen(GVIZ_file, "w")) == NULL)
		fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, GVIZ_file);

	fprintf(out, "Graph TD{\n"); 
	/* You can choose one of the two options below if you want to eliminate node overlap 
     * (the scale option makes a larger, more symmetric drawing than the false option). 
     * The problem is that the scale factors for these graphs are enormous: output files are in the range of 3MB (png) for a 50 node tree, 
     * and PDF rendering fails due to canvas size restrictions. 
     * Current config allows node overlap, which also makes diagrams difficult to read. We should think on how to fix this.
     */
	
	if(style != GV_BAG_LABELS)
		fprintf(out, "overlap=false;\n");
	if(spline == true)
		fprintf(out, "splines=true;\n");  
	if(style == GV_COLORS)
		k = this->compute_width()+1; //stuff to calculate tree-width or max bag size or whatever;
	print_message(1, "k : %d\n", k);

	// Print out the bags
	// Feb 24 2011 - changed back to tree_nodes.size() and check for NULL.
	int size=this->tree_nodes.size();
	Graph::Node *n1;
	for(i=0;i<size;i++)
	{
		if(this->tree_nodes[i] != NULL)
		{
			if(style == GV_BAG_LABELS)
			{
				fprintf(out,"%d [label=\"",i+1);
				print_message(10, "Bag %d:\n", i+1);
				// Sort the bags before outputting
				this->tree_nodes[i]->bag.sort();
				for(L=this->tree_nodes[i]->bag.begin();
					L!=this->tree_nodes[i]->bag.end();++L )
				{
					// Get the actual node label if necessary (we have a pointer to G)
					j=*L;
					print_message(10, "Looking for %d\n", j);
					n1=this->G->get_node(j);
					j=n1->get_label();
					//j=G->nodes[j].label;
					print_message(10, "Found %d\n", j);
					fprintf(out,"%d ",j);
					// NOTE!! The bags store the DIMACS labels (1,2,..,num_total_nodes) so if there
					// are disconnected nodes, then the union of all the bags
					// will NOT be equal to (1,2,...,num_total_nodes)
				}
				fprintf(out,"\"];\n");
			}
			else if(style == GV_SCALE_BAGS)
			{
				// you could set shape to circle or something else here.
				fprintf(out,"%d [width=%f, shape=ellipse];\n",i+1, 
                        (double)(this->tree_nodes[i]->bag.size())/sqrt((double)G->get_capacity()));
			}
			else if(style == GV_TREE_ONLY)
			{
				fprintf(out,"%d;\n",i+1);
			}
			else if(style == GV_COLORS)
			{
				//bags of max size are red
				if((int)(this->tree_nodes[i]->bag.size()) == k)
					sprintf(rgb, "#FF0000");
				//bags with max size - 1 are blue
				else if((int)(this->tree_nodes[i]->bag.size()) == k-1)
					sprintf(rgb, "#0000FF");
				// all smaller bags are gradations of purple, getting lighter as they get smaller
				else
				{
					rcolor = 255- (255/(k-2))*(this->tree_nodes[i]->bag.size());
					sprintf(rgb, "#EE%02xEE", rcolor);
				}
				fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\"];\n",i+1, rgb);
			}
			else
				fatal_error("%s : style specified was not GV_TREE_ONLY, GV_BAG_LABELS, or GV_SCALE_BAGS.\n", __FUNCTION__);
		}
	}

	// Print out the edges
	//Feb 24 2011 - updated to allow for NULL entries in tree_nodes array
	for(i=0;i<size;i++)
	{   
		if(this->tree_nodes[i] != NULL)
		{
			for(L=this->tree_nodes[i]->adj.begin();L!=this->tree_nodes[i]->adj.end();++L)
			{
				j=*L;
				if(i < j)
					// Make this 1-based in true DIMACS spirit
					fprintf(out,"%d -- %d;\n",i+1,j+1);
			}
		}
	}

	fprintf(out, "}\n");

	fflush(out);  
	fclose(out); 
}




/*
  This function writes a graphviz representation (DOT language) of the tree decomposition with node labels colored by the color_vector.
  The color_vector is a vector<double> with a length equal to the number of nodes in the graph with range min_color-max_color.
  The colormap used is defined in get_rgb_value which can be found in Util.cpp.

  If min_color-max_color are not provided, then the first function will calculate them and call the second function.
  Style should be one of GV_COLORS, GV_BAG_LABELS.
 */


void TDTree::write_scored_graphviz_file(bool spline, const char *GVIZ_file, const vector<double>& color_vector, int style, const bool log_range)
{
    size_t length = color_vector.size();
    
    if (!length > 0)
    {
        FERROR("Empty color vector, Abort !!\n");
        return;
    }

    int max_color = color_vector[0];
    int min_color = max_color;

    //loop to find range
    for(size_t i=0;i<length;++i)
    {
        if(color_vector[i] > max_color)
            max_color = color_vector[i];
      
        if(color_vector[i] < min_color)
            min_color = color_vector[i];
    }

  this->write_scored_graphviz_file(spline, GVIZ_file, color_vector, max_color, min_color, log_range, style);

}

void TDTree::write_scored_graphviz_file(bool spline, const char *GVIZ_file, const vector<double>& color_vector, const double & max_color, const double & min_color, int style, const bool log_range)
{
    int i, j; 
    FILE *out; 
    Graph::WeightedMutableGraph *G = this->G;
    list<int>::iterator L;
    char rgb[8];
    int k=0, rcolor; 

    if( (out = fopen(GVIZ_file, "w")) == NULL)
        fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, GVIZ_file);

    fprintf(out, "Graph TD{\n"); 
    /* You can choose one of the two options below if you want to eliminate node overlap 
     * (the scale option makes a larger, more symmetric drawing than the false option). 
     * The problem is that the scale factors for these graphs are enormous: output files are in the range of 3MB (png) for a 50 node tree, 
     * and PDF rendering fails due to canvas size restrictions. 
     * Current config allows node overlap, which also makes diagrams difficult to read. We should think on how to fix this.
     */
	    
  //  fprintf(out, "overlap=true;\n");
  if(spline==true)
    fprintf(out, "splines=true;\n");
  else
    fprintf(out, "splines=false;\n");
    print_message(1, "k : %d\n", k);

    double new_max;
    double new_min;

    //scales and then takes log of color vector, for unevenly distributed color_vectors
    if(log_range)
    {
        new_max = 10;
        new_min = 1;
    }
    else
    {
        new_max = max_color;
        new_min = min_color;
    }

    // Print out the bags
    // Feb 24 2011 - changed back to tree_nodes.size() and check for NULL.
    int size=this->tree_nodes.size();
    Graph::Node *n1;
  
    for(i=0;i<size;i++)
    {
      if(style == GV_BAG_LABELS)
	{
	  double val;
	  if(this->tree_nodes[i] != NULL){
	    fprintf(out,"%d [label=<",i+1);
	    print_message(10, "Bag %d:\n", i+1);
	    
	    //	cout<<i<<"\n";
	    this->tree_nodes[i]->bag.sort();
	    //cout<<i<<" two\n";
	    int count = 0;
	    for(L=this->tree_nodes[i]->bag.begin();
		L!=this->tree_nodes[i]->bag.end();++L )
	      {
		count++;
		//  cout<<i<<" "<<*L<<"\n";
		// Get the actual node label if necessary (we have a pointer to G)
		j=*L;
		print_message(10, "Looking for %d\n", j);
		n1=this->G->get_node(j);
		j=n1->get_label();
		//j=G->nodes[j].label;
		print_message(10, "Found %d\n", j);
		if(log_range)
		  {
		    val = (new_max-new_min)*(color_vector[j-1]-min_color)/(max_color-min_color) + new_min;
		    //	cout<<new_max<<" "<<new_min<<" "<<color_vector[j-1]<<" "<<min_color<<" "<<max_color<<" "<<val<<" \n";
		    val = log(val);
		    
		  }
		else
		  val = color_vector[j-1];

		if(max_color==min_color)
		  val = .5;
		else
		  {
		    if(log_range)
		      val = (val-log(new_min))/(log(new_max)-log(new_min));
		    else
		      val = (val-new_min)/(new_max-new_min);
		    
		  }
		//get hex rgb value to print for font color
		get_rgb_value(rgb,val);
		
		//Must use graphviz html output for labels
		fprintf(out,"<FONT COLOR=\"%s\">",rgb);
		fprintf(out,"%d ",j);
		
		//Need to insert line breaks every so often to prevent really long bags from occurring
		if(count==25)	 
		  {
		    fprintf(out,"<BR/>");
		    count = 0;
		  }
		fprintf(out,"</FONT>\n");
		
		
		// NOTE!! The bags store the DIMACS labels (1,2,..,num_total_nodes) so if there
		// are disconnected nodes, then the union of all the bags
		// will NOT be equal to (1,2,...,num_total_nodes)
	      }
	    fprintf(out,">];\n");
	  }      
	}//end BAG_LABELS

      else if(style == GV_COLORS)
	{
	  if(this->tree_nodes[i] != NULL)
	    {
	      double stat = get_statistics(this->tree_nodes[i]->bag_vec, color_vector, 1); //stat is hard coded right now.			  
	      //	      cout << this->tree_nodes[i]->bag_vec[0] << " " << this->tree_nodes[i]->bag_vec[1] << "\n";
	      //cout << stat << "\n";
	      get_rgb_value(rgb, (stat-min_color)/(max_color-min_color));
	      //cout << (stat-min_color)/(max_color-min_color) << "\n";

	      fprintf(out,"%d [label=\"\",width=%f, shape=circle, style=\"filled\", fillcolor=\"%s\"];\n", i+1,
		      (double)(5*this->tree_nodes[i]->bag.size())/this->width, rgb);	  
	    } 
	}

      else
	{
	  fatal_error("%s : Style specified was not GV_BAG_LABELS or GV_COLOR"); 
	}

    }


    // Print out the edges
    //Feb 24 2011 - updated to allow for NULL entries in tree_nodes array
    for(i=0;i<size;i++)
    {   
        if(this->tree_nodes[i] != NULL)
        {
            for(L=this->tree_nodes[i]->adj.begin();L!=this->tree_nodes[i]->adj.end();++L)
            {
                j=*L;
                if(i < j)
                    // Make this 1-based in true DIMACS spirit
                    fprintf(out,"%d -- %d;\n",i+1,j+1);
            }
        }
    }

    fprintf(out, "}\n");

    fflush(out);  
    fclose(out); 
}





/*
 * Highlights the subtree of the tree decomposition which contains the vertex v. Writes a
 * gviz file to the specified output file. 
 * The final parameter takes one of the options GV_BAG_LABELS, GV_SCALE_BAGS, GV_TREE_ONLY
 * GV_BAG_LABELS lists the entire bag of vertices inside each node's label. 
 * GV_SCALE_BAGS scales the bags to indicate bag size, and labels them with the treenode ID. Note this may 
 * create dimacs files that neato cannot handle on very large graphs! 
 * GV_TREE_ONLY creates small circles for all nodes in the tree decomposition.
 */

void TDTree::highlight_subtree_gviz(int v, const char *GVIZ_file, int style)
{
  
	int i, j; 
	FILE *out; 
	Graph::WeightedMutableGraph *G = this->G;
	list<int>::iterator L;
	char rgb[8];
	//You can modify this to make the uniform circle bags bigger or smaller.
	double fixed_width = .5;
	//for now we'll make the highlights a light blue                
	sprintf(rgb, "#00ced1");



	//this should be a valid test for whether or not the locations have been reocrded.
	//it does not guarantee they are up to date, necessarily.
	if(this->node_locations[0].size() == 0)
	{
		this->record_node_locations();
	}
	int size=(int) this->tree_nodes.size();
	vector<bool> containsv(size,false);
	//use v-1 because indices are 0-n-1, while dimacs node labels are 1-n
	for(L = (this->node_locations[v-1]).begin(); L != (this->node_locations[v-1]).end(); ++L) 
	{
		//	    printf("Node %d is in bag %d \n", v, *L);
		containsv[*L] = true;
	}

	if( (out = fopen(GVIZ_file, "w")) == NULL)
		fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, GVIZ_file);

	fprintf(out, "Graph TD{\n"); 
	//fprintf(out, "overlap=false;\n");
	Graph::Node *n1;
	for(i=0;i<size;i++)
	{
		if(this->tree_nodes[i] != NULL)
		{
			if(style == GV_BAG_LABELS)
			{
				if(containsv[i])
					fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\", label=\"",i+1, rgb);
				else
					fprintf(out,"%d [label=\"",i+1);

				print_message(10, "Bag %d:\n", i+1);
				// Sort the bags before outputting
				this->tree_nodes[i]->bag.sort();
				for(L=this->tree_nodes[i]->bag.begin(); L!=this->tree_nodes[i]->bag.end();++L )
				{
					// Get the actual node label if necessary (we have a pointer to G)
					j=*L;
					print_message(10, "Looking for %d\n", j);
					n1=this->G->get_node(j);
					j=n1->get_label();
					//j=G->nodes[j].label;
					print_message(10, "Found %d\n", j);
					fprintf(out,"%d ",j);
					// NOTE!! The bags store the DIMACS labels (1,2,..,num_total_nodes) so if there
					// are disconnected nodes, then the union of all the bags
					// will NOT be equal to (1,2,...,num_total_nodes)
				}
				fprintf(out,"\"];\n");
			}
			else if(style == GV_SCALE_BAGS)
			{
				if(containsv[i])
					fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\", width=%f, shape=ellipse];\n",i+1, rgb,
                            (double)(this->tree_nodes[i]->bag.size())/sqrt((double)G->get_capacity()));
				else
					fprintf(out,"%d [width=%f, shape=ellipse];\n",i+1, 
                            (double)(this->tree_nodes[i]->bag.size())/sqrt((double)G->get_capacity()));
			}
			else if(style == GV_TREE_ONLY)
			{
				if(containsv[i])
					fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\", width=%f, shape=ellipse];\n",i+1, rgb, fixed_width);
				else
					fprintf(out,"%d [width=%f, shape=ellipse];\n",i+1, fixed_width);


			}
			else
				fatal_error("%s: style parameter was not one of GV_BAG_LABELS, GV_TREE_ONLY, or GV_SCALE_BAGS.\n", __FUNCTION__);

		}
	}	   

	// Print out the edges
	for(i=0;i<size;i++)
	{   
		if(this->tree_nodes[i] != NULL)
		{
			for(L=this->tree_nodes[i]->adj.begin();L!=this->tree_nodes[i]->adj.end();++L)
			{
				j=*L;
				// Make this 1-based in true DIMACS spirit
				if(i < j)
					fprintf(out,"%d -- %d;\n",i+1,j+1);
			}
		}
	}

	fprintf(out, "}\n");

	fflush(out);  
	fclose(out); 
}

/*
  This function writes a graphviz representation (DOT language) of the tree decomposition with node labels colored by the color_vector.
  The color_vector is a vector<double> with a length equal to the number of nodes in the graph with range min_color-max_color.
  The colormap used is defined in get_rgb_value which can be found in Util.cpp.

  If min_color-max_color are not provided, then the first function will calculate them and call the second function.
*/

void TDTree::highlight_subtree_scored_graphviz_file(int v, const char *GVIZ_file, const vector<double>& color_vector)
{
    size_t length = color_vector.size();

    int max_color = color_vector[0];
    int min_color = max_color;

    //loop to find range
    for(size_t i=0;i<length;++i)
    {
        if(color_vector[i] > max_color)
            max_color = color_vector[i];
      
        if(color_vector[i] < min_color)
            min_color = color_vector[i];
    }

    this->highlight_subtree_scored_graphviz_file(v, GVIZ_file, color_vector, max_color, min_color);
}


/*
 * Highlights the subtree of the tree decomposition which contains the vertex v. * Writes a gviz file to the specified output file. 
 * Bag display as the same as scored_graphviz with GV_BAG_LABELS. 
 */
void TDTree::highlight_subtree_scored_graphviz_file(int v, const char *GVIZ_file, const vector<double>& color_vector, const double max_color, const double min_color)
{
    //Blair - still need to add a check whether v is a valid node of G.

	int i, j; 
	FILE *out; 
	Graph::WeightedMutableGraph *G = this->G;
	list<int>::iterator L;
	char node_rgb[8];
	//for now we'll make the highlights a purple to make heat map labels readable
	sprintf(node_rgb, "#f4bbff");
	char rgb[8];
	int k=0, rcolor; 

	//this should be a valid test for whether or not the 
	//locations have been recorded. It does not guarantee they are 
	//up to date, necessarily.
	if(this->node_locations[0].size() == 0)
	{
		this->record_node_locations();
	}

	int size=(int) this->tree_nodes.size();
	vector<bool> containsv(size,false);
	//use v-1 because indices are 0-n-1, while dimacs node labels are 1-n
	for(L = (this->node_locations[v-1]).begin(); L != (this->node_locations[v-1]).end(); ++L)
	{
        //printf("Node %d is in bag %d \n", v, *L);
        containsv[*L] = true;
	}

	if( (out = fopen(GVIZ_file, "w")) == NULL)
		fatal_error("%s:  Error opening file %s for writing graphviz output\n",__FUNCTION__, GVIZ_file);

	fprintf(out, "Graph TD{\n"); 
	
	/* You can choose one of the two options below if you want to eliminate node overlap 
	 * (the scale option makes a larger, more symmetric drawing than the false option). 
	 * The problem is that the scale factors for these graphs are enormous: output files are in the range of 3MB (png) for a 50 node tree, 
	 * and PDF rendering fails due to canvas size restrictions. 
	 * Current config allows node overlap, which also makes diagrams difficult to read. We should think on how to fix this.
	 */
	
	fprintf(out, "overlap=true;\n");
	if(0)//currently disabled - Blair Mar 21, 2012.
        fprintf(out, "splines=true;\n");
	else
        fprintf(out, "splines=false;\n");
	
	print_message(1, "k : %d\n", k);

    Graph::Node *n1;
    double val;
    for(i=0;i<size;i++)
    {
        // Sort the bags before outputting

        if(this->tree_nodes[i] != NULL){

            if(containsv[i])
                fprintf(out,"%d [style=\"filled\", fillcolor=\"%s\", label=<",i+1, node_rgb);
            else
                fprintf(out,"%d [label=<",i+1);
	
            this->tree_nodes[i]->bag.sort();
	
            int count = 0;
            for(L=this->tree_nodes[i]->bag.begin();
                L!=this->tree_nodes[i]->bag.end();++L )
            {
                count++;
	
                // Get the actual node label if necessary (we have a pointer to G)
                j=*L;
                print_message(10, "Looking for %d\n", j);
                n1=this->G->get_node(j);
                j=n1->get_label();
                //j=G->nodes[j].label;
                print_message(10, "Found %d\n", j);
		  
	 
                if(max_color==min_color)
                    val = .5;
                else
                    val = (color_vector[j-1]-min_color)/(max_color-min_color);
		  
	 
                //get hex rgb value to print for font color
                get_rgb_value(rgb,val);
	       
                //Must use graphviz html output for labels
                fprintf(out,"<FONT COLOR=\"%s\">",rgb);
                fprintf(out,"%d ",j);

                //Need to insert line breaks every so often to prevent really long bags from occurring
                if(count==25)	 
                {
                    fprintf(out,"<BR/>");
                    count = 0;
                }
                fprintf(out,"</FONT>\n");
	    
	      
                // NOTE!! The bags store the DIMACS labels (1,2,..,num_total_nodes) so if there
                // are disconnected nodes, then the union of all the bags
                // will NOT be equal to (1,2,...,num_total_nodes)
            }
            fprintf(out,">];\n");
        }
    }

  
    // Print out the edges
    //Feb 24 2011 - updated to allow for NULL entries in tree_nodes array
    for(i=0;i<size;i++)
    {   
        if(this->tree_nodes[i] != NULL)
        {
            for(L=this->tree_nodes[i]->adj.begin();L!=this->tree_nodes[i]->adj.end();++L)
            {
                j=*L;
                if(i < j)
                    // Make this 1-based in true DIMACS spirit
                    fprintf(out,"%d -- %d;\n",i+1,j+1);
            }
        }
    }

    fprintf(out, "}\n");

    fflush(out);  
    fclose(out); 


}






/**
 * Returns the width of the decomposition, defined to be one less than the size
 * of the largest bag. If the value of TDTree::width is nonzero, then this value is
 * returned to avoid recomputing.  Otherwise, loop over the bags and set the value 
 * of width
 */
int TDTree::compute_width()
{
	if(this->width>0)
		return this->width;
	int max_bag_size=-GD_INFINITY;

	int size=(int)this->tree_nodes.size();
	for(int i=0;i<size;i++)
	{
		if(this->tree_nodes[i] != NULL)
		{
			if((int)this->tree_nodes[i]->bag.size()>max_bag_size)
				max_bag_size=this->tree_nodes[i]->bag.size();
		}
	}

	this->width=max_bag_size-1;
	return max_bag_size-1;
}


/**
 * Removes the edge u-v from the tree. Returns false if the edge was not in the
 * tree already, true otherwise.
 */
bool TDTree::remove_tree_edge(int u, int v)
{
	int s1, s2;
	s1=this->tree_nodes[u]->adj.size();
	this->tree_nodes[u]->adj.remove(v);
	s2=this->tree_nodes[u]->adj.size();
	if(s1==s2)
	{
		print_message(0,"%s:  Tree edge u-v not found? (v not in u's adj list)\n",
                      __FUNCTION__);
		return false;
	}

	s1=this->tree_nodes[v]->adj.size();
	this->tree_nodes[v]->adj.remove(u);
	s2=this->tree_nodes[v]->adj.size();
	if(s1==s2)
	{
		print_message(0,"%s:  Tree edge u-v not found? (u not in v's adj list)\n",
                      __FUNCTION__);
		return false;
	}


	return true;
}

/**
 * Adds the tree edge u-v.  Returns false if u-v is already in the tree.
 */
bool TDTree::add_tree_edge(int u, int v)
{
	list<int>::iterator ii;

	for(ii=this->tree_nodes[u]->adj.begin();ii!=this->tree_nodes[u]->adj.end();++ii)
	{
		if(*ii==v)
		{
			print_message(0,"%s:  Tree edge u-v already exists (v in u's adj list)\n",
                          __FUNCTION__);
			return false;
		}
	}

	for(ii=this->tree_nodes[v]->adj.begin();ii!=this->tree_nodes[v]->adj.end();++ii)
	{
		if(*ii==u)
		{
			print_message(0,"%s:  Tree edge u-v already exists (u in v's adj list)\n",
                          __FUNCTION__);
			return false;
		}
	}

	// The edge doesn't already exist --> add it
	this->tree_nodes[v]->adj.push_back(u);
	this->tree_nodes[u]->adj.push_back(v);

	return true;
}

/** Constructs a nice tree decomposition of this->G using
 * a modified version of Kloks' algorithm.  The graph G is assumed to have 
 * a single connected component since tree decompositions
 * of disconnected components can be joined together in
 * an arbitrary way. Furthermore, G is assumed to be triangulated, 
 * and the given elimination ordering perfect for G (zero-fill). 
 * The algorithm simultaneously constructs a k-tree H which contains G. 
 * The maximum clique in G must have size at most k+1.
 * The boolean indicates whether to push new node creation farther 
 * down the tree when possible in the one-child case.  The tree is rooted
 * with root node 0 and we adopt the convention that a tree node's parent
 * is the first node in the adjacency list, followed by the children.
 * Since this is a nice TD, we have at most 2 children.
 */
void TDTree::construct_knice(vector<int> *elim_order, int k, bool descend_one)
{  
	// First make sure this->G has one connected component
	if(this->G->get_num_components()!=1)
	{
		fatal_error("%s:  Tree decomposition routines operate only on connected graphs."
                    "  This graph has %d connected components\n",__FUNCTION__,
                    this->G->get_num_components());
	}

	int i, j; // counters
	int x, w; // vertices
	int minpos, s;
	int tmp;
	// neighbor lists
	list<int> Cw, Cx; 
	// iterators
	list<int>::iterator itx, itw, adjit;
	vector<int>::iterator vv;

	// number of nodes
	int n = (this->G)->get_num_nodes();

	// We're not going to modify G but we do need a graph H to hold our k-tree. Initialize it as a copy of G.
	Graph::WeightedMutableGraph H=*(this->G);

	// This keeps track of a treenode so that the bag associated with it contains C[i] \cup {i}
	// initialized to -1 (an invalid node ID)
	vector<int> bnode(n, -1);

	// t counts the # of nodes in the TDTree
	int t=0;

	// For Kloks, we will have at most 4n vertices in the TDTree
	this->tree_nodes.resize(4*n,NULL);

	// Set the tree node for all k+1 vertices
	for(i = 0; i < k+1; i++)
	{
		bnode[elim_order->at(n-1-i)]=t;
		print_message(1,"initial bag for bnode[%d]\n", elim_order->at(n-1-i));
	}

	// Create the tree_node
	TDTreeNode *next_tree_node;
	next_tree_node=new TDTreeNode;
	next_tree_node->id=t;


	// Make the last  k+1 nodes a clique in H,and set them to be in the first bag of the TD.

	//Create the bag
	for(vv=elim_order->end()-(k+1);vv!=elim_order->end();++vv)
		next_tree_node->bag.push_back(*vv);
	next_tree_node->bag.sort();
	Graph::GraphProperties properties;

	properties.make_clique(&H, &(next_tree_node->bag));

	// It is its own parent, since this is the root node
	next_tree_node->adj.push_front(next_tree_node->id);

	// Add the newly created node to the tree
	this->tree_nodes[t]=next_tree_node;
	t++; 

	print_message(1, "Initialization complete\n");
	//array of bools to indicate the neighbors of x - this is supposedly about 8 times more space efficient than chars.
	vector<bool> xnbr(n, false);
	int numnew;
	int exclude, found;
	TDTreeNode  *currtnode,*child1, *child2; 
	TDTreeNode *newnode1, *newnode2;
	int s1, s2;
	bool moved;
	int parentid; 
	Graph::GraphEOUtil eoutil;

	print_message(1, "Position %d is the lowest index for being in the last (k+1) clique\n", n-k-1);
	for(i=n-(k+2); i >= 0 ; i--)
	{
		//Clean out data structures
		for(j = 0; j < n; j++)
			xnbr[j] = false;
		Cx.clear();
		Cw.clear();

		//Find the set of neighbors of x in G which occur later in the elim. ordering. 
		//The follower will be given by minpos.
		x= elim_order->at(i);
		s = eoutil.find_forward_neighbors(this->G, x, elim_order, i+1, &Cx, &minpos);
		w = elim_order->at(minpos);
		print_message(1, "x: %d s: %d w: %d k: %d minpos: %d\n", x,s,w,k,minpos);

		if(s <= k)
		{
			if(minpos > n-(k+2)) //We need to set Cw to be the last k+1 clique
			{ 
				for(j = n-(k+1); j < n; j++)
				{
					if(j != minpos)
						Cw.push_back(elim_order->at(j));
				}
			}
			else  //Find the set of neighbors of w in H which are later in the order. 
			{
				tmp = eoutil.find_forward_neighbors(&H, w, elim_order, minpos+1, &Cw, &tmp);
				if(tmp != k)
					fatal_error("%s: Error - already processed vertex does not have k forward neighbors in H\n", __FUNCTION__);
			}

			print_message(1, "Found the forward neighbors\n");        

			//Add edges to H from x to (k-s) elements of Cw
			//Transform Cx into a 0/1 incidence vector.
			itx = Cx.begin();
			while(itx != Cx.end())
			{
				xnbr[*itx] = true;
				itx++;
			}

			//iterator for Cw
			itw = Cw.begin();
			numnew = 0;
			exclude = Graph::rand_int(1, k-(s-1));
			//choose the element of Cw to exclude
			while(itw != Cw.end() && s < k)
			{
				if(xnbr[*itw] == false)
				{
					numnew++;
					if(numnew != exclude)
					{
						xnbr[*itw] = true; //this way xnbr will be exactly the set of forward neighbors of x
						H.add_edge(x, *itw);
						s++;
					}
				}
				++itw;
			}
		}//end of loop that adds edges to make H a k-tree

		if(s != k)
			fatal_error("You have failed!\n");

		//Find the starting tree node
		currtnode = this->tree_nodes[bnode[w]]; 
		print_message(1, "Starting the search at bnode[%d] = %d\n", w, bnode[w]);
		print(1, (currtnode->bag));

		while((currtnode->adj).size() > 1)
		{
			while((currtnode->adj).size() == 3) // this means it has two children
			{
				print_message(1, "Tree node %d has two children\n", currtnode->id);
				adjit = (currtnode->adj).begin();
				//find the two children
				adjit++;
				child1 = this->tree_nodes[*adjit]; 
				adjit++; 
				child2 = this->tree_nodes[*adjit]; 
				s1 = (child1->adj).size(); 
				s2 = (child2->adj).size(); 
				parentid = currtnode->id;
				// choose the one with less children
				if(s1 <= s2)
					currtnode = child1; 
				else 
					currtnode = child2;

				// ADDED CHECK HERE
				if(bnode[w] == parentid)
				{
					bnode[w] = currtnode->id; //move bnode so we dont' have to do this again if w is a follower.
					print_message(1, "We changed bnode for w = %d\n", w);
				}
			}
			if((currtnode->adj).size() == 2) // single child
			{
				print_message(1, "Tree node %d has one child\n", currtnode->id);
				adjit = (currtnode->adj).begin();
				adjit++; 
				child1 = this->tree_nodes[*adjit]; 
				moved = false;

				if(descend_one) // check the child's bag to see if we can move downwards and try again.
				{
					found = 0;
					itx = child1->bag.begin(); 
					while(itx != child1->bag.end())
					{
						if(xnbr[*itx] == true)
							found++;
						++itx;
					}
					if(found == k) // change currtnode and return to top of loop
					{
						print_message(1, "We descended!\n");
						currtnode = child1;
						moved = true;
					}
				}
				if(moved == false)
				{
					// create two new children for currtnode
					newnode1 =new TDTreeNode;
					newnode2 = new TDTreeNode;
					newnode1->id=t;
					newnode2->id=t+1;
					// set the bags to be the same as the current node
					newnode1->bag = currtnode->bag;
					newnode2->bag = currtnode->bag; 

					// set the parent to currtnode - give newnode1 child1 and newnode1 no children
					newnode1->adj.push_front(child1->id);
					newnode1->adj.push_front(currtnode->id);
					newnode2->adj.push_front(currtnode->id);

					// change the parent of child1
					child1->adj.pop_front();
					child1->adj.push_front(newnode1->id);

					// remove child1 and add newnode1 and newnode2 as children of currtnode
					(currtnode->adj).erase(adjit);
					(currtnode->adj).push_back(newnode1->id);
					(currtnode->adj).push_back(newnode2->id);

					// Add the newly created nodes to the tree
					this->tree_nodes[t]=newnode1;
					this->tree_nodes[t+1] = newnode2; 

					//increment the node counter by 2
					t = t+2;            

					// this appears to be the problem - why do we do this?
					// CHANGED
					if(bnode[w] == currtnode->id){
						print_message(1, "We're changing bnode (newnode2) for w = %d from %d to %d \n", w, bnode[w], newnode2->id);
						// printf("from \n");
						print(1, (this->tree_nodes[bnode[w]]->bag));
						// printf("to \n");
						print(1, (newnode2->bag));
						bnode[w] = newnode2->id; 
					}
					// move the current node to the newnode2 leaf
					currtnode = newnode2;


				}//end of modification steps
			}//end of 2 case
		}//end of loop for degree 2/3

		//currtnode is a leaf! 
		print_message(1, "Tree node %d is a leaf\n", currtnode->id);

		//currnode contains an extra vertex, which we must find and remove in a new child bag.
		if((int)currtnode->bag.size() == k+1) 
		{
			// create a new child
			newnode1 =new TDTreeNode;

			// create a bag with these neighbors plus x
			for(itx=currtnode->bag.begin();itx!=currtnode->bag.end();++itx)
				newnode1->bag.push_back(*itx);

			adjit=newnode1->bag.begin();
			// remove the vertex that was not a neighbor of x
			while(adjit != newnode1->bag.end())
			{
				if(xnbr[*adjit] == false)
				{
					print_message(1, "Removed non-neighbor %d from the bag\n", *adjit);
					newnode1->bag.erase(adjit);
					break;
				}
				++adjit;
			}


			newnode1->id=t;

			// set the parent to currtnode and set it as a child of that node
			(newnode1->adj).push_front(currtnode->id);
			(currtnode->adj).push_back(newnode1->id);

			// Add the newly created node to the tree
			this->tree_nodes[t]=newnode1;

			//increment the node and bag counters 
			t = t+1;            
			//numbags++;

			currtnode = newnode1;
		}

		// the currtnode is a leaf with bag exactly the forward neighbors of x
		print_message(1, "Tree node %d is a perfect leaf!\n", currtnode->id);

		// create a new child        
		// create a bag with these neighbors plus x
		newnode1 =new TDTreeNode;
		newnode1->id=t;
		for(itx=currtnode->bag.begin();itx!=currtnode->bag.end();++itx)
			newnode1->bag.push_back(*itx);
		newnode1->bag.push_back(x);
		newnode1->bag.sort();

		// set the parent to currtnode and set it as a child of that node
		(newnode1->adj).push_front(currtnode->id);
		(currtnode->adj).push_back(newnode1->id);

		// Add the newly created node to the tree
		this->tree_nodes[t]=newnode1;

		// increment the node counter
		t++;
		bnode[x] = newnode1->id; //set the bnode of x

	}//end of loop over simplicial vertices
	int size=(int)this->tree_nodes.size();
	for(i=0;i<size;i++)
	{
		if(this->tree_nodes[i] && this->tree_nodes[i]->adj.size()==1)
			this->num_leafs++;
	}

	this->num_tree_nodes = t;
	this->root_node = 0;
	this->rooted=true;
	this->nice=true;

	return;
}

/**
 * Manually verifies if the decomposition is nice.  Returns true and sets nice=true 
 * if the decomposition is ``nice'' and false otherwise.
 * Modified Feb 24, 2011 to fix problems not detecting anomalies and allow for "gaps" in tree_nodes vector.
 */
bool TDTree::is_nice()
{
	// First make sure it is valid!
	if(!this->verify())
	{
		fatal_error("TD is not valid!\n");
		return false;
	}

	int i=0;
	list<int>::iterator ii,jj,kk;
	vector<TDTreeNode *>::iterator tt;
	TDTreeNode *curr, *child1, *child2;

	for(tt = this->tree_nodes.begin(); tt != this->tree_nodes.end(); ++tt)
	{   
		if(*tt != NULL && (*tt)->id != GD_UNDEFINED)//check for a missing node in the vector
		{
			i++;
			curr = *tt;
			print_message(1,"Testing tree node %d of %d\n",i,this->num_tree_nodes);

			// Make sure that all tree nodes have at most 3 neighbors
			if(curr->adj.size()>3)
			{
				print_message(0,"Tree Node %d has more %d>3 neighbors\n",i,curr->adj.size());
				this->nice=false;
				return false;
			}

			// If tree node has degree 3, then the last two in the list are the children
			// Make sure bags are all the same
			else if( curr->adj.size()==3)
			{
				ii=curr->adj.begin();
				++ii;
				child1 = this->tree_nodes[*ii];
				++ii;
				child2 = this->tree_nodes[*ii];
				// Make sure bags are ==
				if( (curr->bag) != (child1->bag))
				{
					print_message(0,"Bag not equal to first child's bag!\n");
					this->nice=false;
					return false;
				}
				if( (curr->bag) != (child2->bag))
				{
					print_message(0,"Bag not equal to second child's bag!\n");
					this->nice=false;
					return false;
				}
			}

			// One child case
			else if( curr->adj.size()==2)
			{
				ii=curr->adj.begin();
				++ii;
				child1 = this->tree_nodes[*ii];


				// Make sure bags have proper relationship
				if( curr->bag.size()==child1->bag.size()+1 )
				{
					// The parent's bag is larger - make sure the child's is contained in it
					// The bags should be sorted!

					if(!includes(curr->bag.begin(),curr->bag.end(),
                                 child1->bag.begin(),child1->bag.end()))
					{
						print_message(0,"Parent larger:  Do not appear to have proper bag containment!\n");
						print(0,curr->bag);
						print(0,(child1->bag));
						this->nice=false;
						return false;
					}

				}
				else if(curr->bag.size() == child1->bag.size() -1)
				{
					if(!includes(child1->bag.begin(),child1->bag.end(),
                                 curr->bag.begin(),curr->bag.end()))
					{
						print_message(0,"Child larger:  Do not appear to have proper bag containment!\n");
						print(0,curr->bag);
						print(0,(child1->bag));
						this->nice=false;
						return false;
					}
				}
				else
				{
					print_message(0,"Bags have improper size relationship!\n");
					this->nice=false;
					return false;
				}

			}//one child case

		}//non-null node
	}

	// We appeared to pass the tests
	this->nice=true;
	return true;
}

/**
 * Finds & sets a root node to minimize the maximum number of tables
 * required to be stored at any given time (assuming no reconstruction is required).
 * Reference: "Memory Requirements for Table Computations" by Aspvall et al.
 * algorithm needs to be one of the following: 
 * TD_ROOT_ASPVALL - original tabreq() function assumes child tables are deleted sequentially due to DFS
 * TD_ROOT_ALLCHILDREN - tries to account for all child tables needing to exist simultaneously to create parent.
 */
int TDTree::root_mintables(int algorithm)
{

	int i = 1;
	int v,w, j;
	list<int>::iterator it;

	//make a copy since we're going to remove nodes as we go
	TDTree T(*this);

	int max_stacks = T.num_tree_nodes+1;    

	//array of stacks where S[i] has leaf vertices of remaining tree with subtrees requiring i tables
	vector< list<int> > S(max_stacks);

	//vector to hold larg and next.larg values (stored in p1 and p2, resp)
	vector<int_int> tab_req(T.tree_nodes.size());
	vector<int> scores(T.tree_nodes.size()); // keeps track of the total table requirements

	int num_processed = 0;

	//cout << "processing tree nodes\n";

	//push all degree one vertices of T on S[1]
	//Simultaneously initialize p1 and p2 to 0 for all tree nodes
	for(int j = 0; j < (int) T.tree_nodes.size(); j++)
	{
		if(T.tree_nodes[j] != NULL)
		{
			{
				tab_req[j].p1 = 0; 
				tab_req[j].p2 = 0;
				scores[j] = 0;
			}
			if(T.tree_nodes[j]->is_leaf())
			{
				S[1].push_front(j);
				scores[j] = 1;
				//cout << "leaf " << j << "\n";
			}
		}	
	}


	cout << "Done with leaves. Num tree nodes: " << T.num_tree_nodes << "\n";

	//While T still has more than one vertex left
	while(num_processed < T.num_tree_nodes - 1 )
	{
		i = 1; 

		//find the smallest i for which S[i] is nonempty
		while(S[i].size() == 0)
			i++;

		//Pop a vertex w from S[i] and let v be the (only) neighbor of w which is still in T
		w = S[i].front();
		S[i].pop_front();
		v = T.tree_nodes[w]->adj.front(); 

		//remove w from v's adjacency list
		T.tree_nodes[v]->adj.remove(w);
		
		//remove w from the tree
		//free(T.tree_nodes[w]);
		delete T.tree_nodes[w];
		T.tree_nodes[w] = NULL; 	   
		num_processed++;


		//Update the table requirement values of v. 
		switch(algorithm)
		{
            case TD_ROOT_ASPVALL: 
                tab_req[v].p2 = max(min(i, tab_req[v].p1), tab_req[v].p2); 
                tab_req[v].p1 = max(i, tab_req[v].p1); 
                break; 
            case TD_ROOT_ALLCHILDREN: 
                tab_req[v].p2 = max(tab_req[v].p1 + i, tab_req[v].p2); 
                tab_req[v].p1++; 
                break; 
		}    

		//if v has only one neighbor left, we push it onto a stack & decrease num nodes in T.
		if(T.tree_nodes[v]->adj.size() == 1)
		{
			switch(algorithm)
			{
                case TD_ROOT_ASPVALL: 
                    j = max(tab_req[v].p1, 1+tab_req[v].p2); 
                    j = max(j, 2);
                    break; 
                case TD_ROOT_ALLCHILDREN: 
                    j = tab_req[v].p1 + tab_req[v].p2;
                    j = max(j, 2);
                    break; 
			}    

			//cout << "Trying to add to stack " << j << "\n";
			S[j].push_front(v); 
			scores[v] = j;
			//cout << "Added " << v << " to stack " << j << "\n";
		}

	}

	//there should only be one left. - note that this calls size() every time through the loop!
	for(j = 0; j < (int)T.tree_nodes.size(); j++)
	{
		if(T.tree_nodes[j] != NULL)
		{
			//		    cout << "root " << j << " has tab_req = " << max(tab_req[j].p1, 1+tab_req[j].p2);/* << " which should be at most " << ( mylog2((4*(T.num_tree_nodes + 1)) / 3) ).p1 << endl;*/
			switch(algorithm)
			{
                case TD_ROOT_ASPVALL: 
                    this->root(j);
                    break;
                case TD_ROOT_ALLCHILDREN: 
                    list<int> myS;
                    list<int>::iterator ii;
                    int t;

                    this->rooted=true;
                    this->root_node=j;

                    //sort the other neighbors in decreasing order based on tab_req
                    score_sort(&((this->tree_nodes[j])->adj), &scores);

                    // Set the root's parent to itself (leave the other existing neighbors intact)
                    this->tree_nodes[j]->adj.push_front(j);

                    myS.push_front(j);
                    while(!myS.empty())
                    {
                        t=myS.front();
                        //			cout << "processing node " << t << endl;
                        // Reorder the adj. lists of all t's children, setting t to be the new parent
                        ii=this->tree_nodes[t]->adj.begin();
                        // Advance past the parent
                        ++ii;

                        while(ii!=this->tree_nodes[t]->adj.end())
                        {
                            //sort the neighbors in decreasing order based on tab_req
                            score_sort(&(this->tree_nodes[*ii]->adj), &scores);

                            if(*(this->tree_nodes[*ii]->adj.begin())!=t)
                            {
                                // Remove t from the list of adjacent nodes and put it at the front
                                // since t is the parent
                                int orig_size=(int)this->tree_nodes[*ii]->adj.size();
                                this->tree_nodes[*ii]->adj.remove(t);
                                if((int)this->tree_nodes[*ii]->adj.size()!=orig_size -1)
                                    print_message(0, "Didn't find alleged parent in existing neighbors??\n");
                                this->tree_nodes[*ii]->adj.push_front(t);
                            }

                            // Add *ii to the end of the list
                            myS.push_back(*ii);
                            ++ii; 
                        }
                        myS.pop_front();
                    }			    
			}//end of switch statement
			return j;
		}
	}

	//error.
	return -1;
}



/**
 * Roots the tree at node k and reorders all adjacency lists in the tree so that
 * parent is first followed by children.  The root's parent is itself.
 */
void TDTree::root(int k)
{
	list<int> S;
	list<int>::iterator ii;
	int t;


	this->rooted=true;
	this->root_node=k;
	// Set the root's parent to itself (leave the other existing neighbors intact)
	this->tree_nodes[k]->adj.push_front(k);

	S.push_front(k);
	while(!S.empty())
	{
		t=S.front();
		// Reorder the adj. lists of all t's children, setting t to be the new parent
		ii=this->tree_nodes[t]->adj.begin();
		// Advance past the parent
		++ii;

		while(ii!=this->tree_nodes[t]->adj.end())
		{
			if(*(this->tree_nodes[*ii]->adj.begin())!=t)
			{
				// Could do size check here to make sure we find t
				// Remove t from the list of adjacent nodes and put it at the front
				// since t is the parent
				int orig_size=(int)this->tree_nodes[*ii]->adj.size();
				this->tree_nodes[*ii]->adj.remove(t);
				if((int)this->tree_nodes[*ii]->adj.size()!=orig_size -1)
					print_message(0, "Didn't find alleged parent in existing neighbors??\n");

				this->tree_nodes[*ii]->adj.push_front(t);
			}
			// Add *ii to the end of the list
			S.push_back(*ii);
			++ii; // CSG, June 10 - this was missing?
		}
		S.pop_front();
	}


	return;

}

/**
 * Computes a post order walk of the tree so that if we process the elements in the 
 * order given by the walk[], then we know that we will never encounter a descendant
 * of a previously processed node.
 * Important: This is a walk that results from a DFS starting at the root, as required for Aspvall.
 * This is a Breadth First Search.
 */
void TDTree::post_order_walk(vector<int> *walk)
{
	if((int)walk->capacity()< this->num_tree_nodes)
	{
		print_message(0,"%s:  walk array is too small!\n",__FUNCTION__);
		return;
	}

	list<int> S;
	int t = 0;
	int i=this->num_tree_nodes-1;
	int size=(int)this->tree_nodes.size();
	char *processed=new char[size];
	memset(processed,0,(int)size*sizeof(char));
	list<int>::iterator ii;

	S.push_front(this->root_node);
	while(!S.empty())
	{
		t=S.front();
		print_message(1,"S.front is %d; size is %d\n",t,S.size());
		if(processed[t]==1)
		{
			print_message(1,"Setting walk[%d]=%d\n",i,t);
			walk->at(i)=t;
			i--;
			S.pop_front();
		}
		else
		{
			processed[t]=1;
			// Now get the children of t
			print_message(1,"%d has %d neighbors\n",t,this->tree_nodes[t]->adj.size());

			ii=this->tree_nodes[t]->adj.begin();
			++ii;
			while(ii!=this->tree_nodes[t]->adj.end())
			{
				S.push_back(*ii);
				++ii;
			}
		}
	}

	delete [] processed;

}

void TDTree::post_order_walk_DFS(vector<int> *walk)
{
	if((int)walk->capacity()< this->num_tree_nodes)
	{
		print_message(0,"%s:  walk array is too small!\n",__FUNCTION__);
		return;
	}

	list<int> S;
	int t;
	int i=this->num_tree_nodes-1;
	int size=(int)this->tree_nodes.size();
	char *processed=new char[size];
	memset(processed,0,(int)size*sizeof(char));
	list<int>::iterator ii;
	S.push_front(this->root_node);
	while(!S.empty())
	{
		t=S.front();
		walk->at(i)=t;
		i--;
		processed[t]=1;
		S.pop_front();
		ii=this->tree_nodes[t]->adj.begin();
		++ii;
		while(ii!=this->tree_nodes[t]->adj.end())
		{
			if(!processed[*ii])
				S.push_front(*ii);
			++ii;
		}
	}

	delete [] processed;
}


/**
 * Computes a pre order walk of the tree so that if we process the elements in the 
 * order given by the walk[], then we know that we will never encounter an ancestor
 * of a previously processed node.
 */
void TDTree::pre_order_walk(vector<int> *walk)
{
	if((int)walk->capacity()< this->num_tree_nodes)
	{
		print_message(0,"%s:  walk array is too small!\n",__FUNCTION__);
		return;
	}

	list<int> S;
	int t;
	int i=0;
	//Feb 24 2011 - allowing size() != num_tree_nodes
	int size=(int)this->tree_nodes.size();
	char *processed=new char[size+1];
	memset(processed,0,(int)(size+1)*sizeof(char));
	list<int>::iterator ii;

	S.push_front(this->root_node);
	while(!S.empty())
	{
		t=S.front();
		print_message(1,"S.front is %d; size is %d\n",t,S.size());
		if(processed[t]==1)
		{
			print_message(1,"Setting walk[%d]=%d\n",i,t);
			walk->at(i)=t;
			i++;
			S.pop_front();
		}
		else
		{
			processed[t]=1;
			// Now get the children of t
			print_message(1,"%d has %d neighbors\n",t,this->tree_nodes[t]->adj.size());

			ii=this->tree_nodes[t]->adj.begin();
			++ii;
			while(ii!=this->tree_nodes[t]->adj.end())
			{
				S.push_back(*ii);
				++ii;
			}
		}
	}

	delete [] processed;
}


/**
 * Returns the type of node if the decomposition is nice. Possible return values are
 * TD_LEAF_NODE, TD_INTRODUCE_NODE, TD_FORGET_NODE, TD_JOIN_NODE, TD_NONLEAF_NODE.
 */
int TDTree::get_node_type(int k)
{
	if(this->tree_nodes[k]->adj.size()==1 && this->nice==true)
		return TD_NICE_LEAF_NODE;

	if(this->tree_nodes[k]->adj.size()==1 && this->nice==false)
		return TD_NONNICE_LEAF_NODE;

	if(!this->nice)
		// Not nice - must be non-leaf node if we get here
		return TD_NONLEAF_NODE;       

	// Must be nice and not leaf
	if(this->tree_nodes[k]->adj.size()==3)
		return TD_JOIN_NODE;

	// We must have deg 2 - 
	// if parent is larger, introduce
	// if child is larger, forget

	int child=this->tree_nodes[k]->adj.back();
	if(this->tree_nodes[k]->bag.size() > this->tree_nodes[child]->bag.size())
		return TD_INTRODUCE_NODE;

	return TD_FORGET_NODE;
}


/**
 * Computes the table of node tnode in the tree, using the user-supplied
 * table_function.  This is a function that takes in a TDTree and an integer k
 * and then computes a problem-specific table for this node.  This function should
 * return the size of the table for this node. Has a third parameter that allows
 * some options on the computation.
 */
int TDTree::compute_table(int(*table_function)(TDTree *T, int k),int tnode)
{
	return table_function(this,tnode);

}


/**
 * Fills the TDTree S with the subtree induced by the provided indices. Currently
 * does not verify that the resulting TDTree is, in fact, a tree.  All edge information
 * between the tree nodes is translated to use the new labels.  The only neighbors included
 * are those in the indices list!
 */
bool TDTree::create_subtree(list<int> *indices, TDTree *S)
{
	int i;
	list<int>::iterator ii,jj;

	S->num_tree_nodes=indices->size();
	S->G=this->G;

	S->tree_nodes.resize(S->num_tree_nodes);
	for(i=0;i<S->num_tree_nodes;i++)
		S->tree_nodes[i]=new TDTreeNode();

	// Create an array to do the renumbering, etc.
	vector<int> perm(indices->size(),-1);
	vector<int> iperm(this->num_tree_nodes+1,-1);
	i=0;
	for(ii=indices->begin();ii!=indices->end();++ii)
	{
		perm[i]=*ii;
		iperm[*ii]=i;
		i++;
	}

	for(i=0;i<S->num_tree_nodes;i++)
	{
		// Copy it, but then we will have to manually adjust some data
		*(S->tree_nodes[i]) = *(this->tree_nodes[perm[i]]);
		S->tree_nodes[i]->id=i;

		// The bag is ok as it represents vertices from the orig. graph

		// Need to fix the adj. list
		S->tree_nodes[i]->adj.clear();
		print_message(1,"Adjusting adj list for new node %d (old index %d)\n",i,perm[i]);
		cerr<<*(this->tree_nodes[perm[i]]);
		for(jj=this->tree_nodes[perm[i]]->adj.begin();jj!=this->tree_nodes[perm[i]]->adj.end();++jj)
		{
			print_message(1,"Orig. tree node %d -> %d; neighbor is %d->%d\n",i,perm[i],*jj,iperm[*jj]);
			if(iperm[*jj]!=-1)
			{
				// This will make sure we store only the edges to other nodes in the subtree
				// Also, if we are copying the root, don't preserve that property
				if(iperm[*jj]!=i)
					S->tree_nodes[i]->adj.push_back(iperm[*jj]);
			}
		}
	}

	return true;

}


void TDTree::clear()
{
	this->G=NULL;
	this->nice=false;
	this->node_locations->clear();
	this->num_tree_nodes=0;
	this->rooted=false;
	// This should clear out all the TDTreeNode information
	this->tree_nodes.clear();
	this->width=0;

}

/*
 * For rooted tree decompositions only (fatal error if called on unrooted).
 * Populates the vector B (which will be sized to equal the number of bags in the TD)
 * with B[i] = the number of vertices in bag i which are also in the bags of at least 2 children of i. 
 * If i is a leaf, B[i] is the size of the bag (aka the width).
 */
int TDTree::find_child_boundaries(vector<int> *B)
{

	int i; 
	list<int> *currbag; 
	TDTreeNode *currnode;
	list<int>::iterator it, jt; 

	if(!this->rooted)
		fatal_error("%s: tried to compute child boundaries on non-rooted TD.\n", __FUNCTION__);

	if(B->size() != this->tree_nodes.size()){
		B->resize((int)this->tree_nodes.size(), 0);
	}

	//do a post order walk. this guarantees adj vectors for the children have been computed.
	vector<int> walk(this->num_tree_nodes, 0);
	this->post_order_walk(&walk);    

	//compute the adjacency vectors for each bags while simulaneously summing over its children
	vector< vector<int> > adjvecs;
	int ccount;
	adjvecs.resize((int)this->tree_nodes.size());
	for(i = 0; i < this->num_tree_nodes; i++)
	{
		adjvecs[walk[i]].resize((this->G)->get_capacity(), 0);
		currnode = this->tree_nodes[walk[i]];
		currbag = &(currnode->bag);
		if(currnode->adj.size() > 1)
		{
			for(it = currbag->begin(); it != currbag->end(); ++it)
			{
				ccount = 0;
				(adjvecs[walk[i]])[*it] = 1;
				jt = (currnode->adj).begin(); 
				++jt; //skip the parent, we only want children
				while(jt != (currnode->adj).end()) 
				{
					if((adjvecs[*jt])[*it] == 1)
						ccount++;
					++jt;
				}
				if(ccount > 1) //found a member of this bag which is in >=2 children.
					(*B)[walk[i]]++;
			}
		}
		else
		{
			//set the adjacency list for this one
			for(it = currbag->begin(); it != currbag->end(); ++it)
				(adjvecs[walk[i]])[*it] = 1;
			//set B to be the size of the bag.
			(*B)[walk[i]] = currbag->size();
		}
	}
	return 0;    
}


void TDTree::free_children(int k)
{
	list<int>::iterator ii=this->tree_nodes[k]->adj.begin(); 
	// Advance past the parent
	++ii;
	// Delete the children
	for( ; ii!=this->tree_nodes[k]->adj.end(); ++ii)
	{
		// Don't make it NULL - just clear the hash table!
		if(this->tree_nodes[*ii]->hash_table)
		{
			TDSolution *current_set, *tmp;
			HASH_ITER(hh, this->tree_nodes[*ii]->hash_table, current_set, tmp)
			{
				HASH_DEL(this->tree_nodes[*ii]->hash_table,current_set);
				delete current_set;
			}
			this->tree_nodes[*ii]->hash_table=NULL;
		}
	}
	return;

}

void TDTree::free_table(int k)
{

	// Don't make it NULL - just clear the hash table!
	if(this->tree_nodes[k]->hash_table)
	{
		TDSolution *current_set, *tmp;
		HASH_ITER(hh, this->tree_nodes[k]->hash_table, current_set, tmp)
		{
			HASH_DEL(this->tree_nodes[k]->hash_table,current_set);
			delete current_set;
		}
		this->tree_nodes[k]->hash_table=NULL;
	}

	return;

}

void TDTree::fill_bag_vecs()
{
	int i,k;
	list<int>::iterator ii;
	int size=(int)this->tree_nodes.size();
	for(k=0;k<size;k++)
	{
		if(this->tree_nodes[k])
		{
			this->tree_nodes[k]->bag_vec.resize(this->tree_nodes[k]->bag.size());
			i=0;
			for(ii=this->tree_nodes[k]->bag.begin();ii!=this->tree_nodes[k]->bag.end();++ii)
			{
				this->tree_nodes[k]->bag_vec[i]=*ii;
				i++;
			}
		}
	}
}

/*
 * Goes through the tree and removes all nodes from the list T from the bags
 * of all tree nodes 
 */
void TDTree::refine_tree(list<int> *T)
{
	vector<int> I(this->width+1,-1);
	vector<int>::iterator uu,vv;
	int i,max_bag_size=-1;
	for(i=0;i<this->num_tree_nodes;i++)
	{
		// Fill I with all vertices in i's bag that are not in T - we have to keep these
		// around as they might be in the opt. solution
		vv=set_difference(this->tree_nodes[i]->bag.begin(),this->tree_nodes[i]->bag.end(),
                          T->begin(),T->end(),I.begin());
		this->tree_nodes[i]->bag.clear();
		this->tree_nodes[i]->bag_vec.clear();
		for(uu=I.begin();uu!=vv;++uu)
		{
			this->tree_nodes[i]->bag.push_back(*uu);
			//this->tree_nodes[i]->bag_vec.push_back(*uu);
		}
		this->tree_nodes[i]->bag.sort();
		this->tree_nodes[i]->bag_vec.insert(this->tree_nodes[i]->bag_vec.end(), this->tree_nodes[i]->bag.begin(), this->tree_nodes[i]->bag.end());

		if((int)this->tree_nodes[i]->bag.size()>max_bag_size)
			max_bag_size=(int)this->tree_nodes[i]->bag.size();
	}
	// Reset tree width - not currently doing this since it messes up the WIS output
	// Not sure what the right thing to do is...
	// this->width=max_bag_size-1;
	return;

}

/*
 * Fills the intersection_sizes[] array with the size of the intersection of each
 * node in the tree with its parent.  Assigns size of 0 to root node.
 * Assumes that intersection_sizes[] array is big enough to hold everything.
 */
void TDTree::compute_parent_child_intersections(vector<int> *intersection_sizes)
{
	int i;
	vector<int> I(this->width+1);
	vector<int>::iterator uu;

	int size=(int)this->tree_nodes.size();
	for(i=0;i<size;i++)
	{
		// CSG adding July 5
		if(this->tree_nodes[i])
		{
			int current_node=i;
			int parent_node=this->tree_nodes[i]->adj.front();

			// if current_node=parent_node, then this is the root - assign 0
			if(current_node==parent_node)
				intersection_sizes->at(current_node)=0;
			else
			{
				uu=set_intersection(this->tree_nodes[current_node]->bag.begin(),
                                    this->tree_nodes[current_node]->bag.end(),
                                    this->tree_nodes[parent_node]->bag.begin(),
                                    this->tree_nodes[parent_node]->bag.end(),
                                    I.begin());
				intersection_sizes->at(current_node)=uu-I.begin();
			}
		}
	}
}


/**
 * Fills vertices list with all vertices contained in the bag of tree node
 * v or any node beneath v in the tree.
 */
void TDTree::compute_subtree(int v, list<int> *subtree)
{
	int j;
	list<int>::iterator ii,jj;
	list<int> S;
	// Do a BFS to get a list of all tree nodes below v in the tree
	subtree->clear();
	subtree->push_back(v);
	ii=this->tree_nodes[v]->adj.begin();
	++ii;
	for( ;ii!=this->tree_nodes[v]->adj.end();++ii)
	{
		S.push_back(*ii);
		subtree->push_back(*ii);
	}

	while(!S.empty())
	{
		// While we still have nodes on the stack, keep searching
		// Get the first entry in S, remove it and then add his children
		// to the list
		j=S.front();
		S.pop_front();
		ii=this->tree_nodes[j]->adj.begin();
		++ii;
		// This loop doesn't do anything if j is a leaf node
		for( ;ii!=this->tree_nodes[j]->adj.end();++ii)
		{
			S.push_back(*ii);
			subtree->push_back(*ii);
		}
	}

	return;
}

void TDTree::compute_subtree_vertices(int v, list<int> *vertex_list)
{
	// Compute the subtree
	list<int> subtree;
	this->compute_subtree(v,&subtree);

	// Now get the union of all the bags in the subtree
	vertex_list->clear();
	for(list<int>::iterator ii=subtree.begin();ii!=subtree.end();++ii)
	{
		for(list<int>::iterator jj=this->tree_nodes[*ii]->bag.begin();
			jj!=this->tree_nodes[*ii]->bag.end();++jj)
			vertex_list->push_back(*jj);
	}

	vertex_list->sort();
	vertex_list->unique();
}

/**
 * Goes through tree from leaves up, removing vertices from unnecessary "higher" bags
 * once they and all their neighbors have been encountered. This should help trim bag
 * sizes on decompositions arising from triangulation where vertices are needed for the 
 * TD to be valid for the chordal graph, but not the original. After sweep, algorithm 
 * removes all duplicate and subset bags.
 */
void TDTree::refine()
{
	//first, remove duplicates and subsets
	this->remove_subset_bags(); 

	list<int> extraneous, subtree;
	vector<bool> in_subtree(this->G->get_capacity());
	vector<bool> subnode(this->tree_nodes.size());

	vector<int> walk(this->num_tree_nodes);
	this->post_order_walk(&walk); 

	list<int> st_vertices;
	list<int>::iterator eit, vit, bagit, nbrit; 
	int curr, cind; 
	bool extra;
	int s;
	Graph::Node *n1;
	list<int> *nbr_list;
	//no need to do this at the root (thus the -1)
	for(cind = 0; cind < this->num_tree_nodes-1; cind++)
	{
		curr = walk[cind];    

		this->compute_subtree_vertices(curr,&st_vertices);
		this->compute_subtree(curr,&subtree);
		print_message(1,"Subtree vertices for node %d:\n", curr);
		print(1,st_vertices);

		//fill in an incidence vector for subtree to make neighbor lookups easy
		in_subtree.assign(this->G->get_capacity(), false); 
		for(vit = st_vertices.begin(); vit != st_vertices.end(); ++vit){
			in_subtree[*vit] = true;
		}

		subnode.assign(this->tree_nodes.size(), false); 
		for(vit = subtree.begin(); vit != subtree.end(); ++vit){
			subnode[*vit] = true;
		}

		extraneous.clear();

		//work through the vertices in the current bag. 
		bagit=this->tree_nodes[curr]->bag.begin();
		while(bagit!=this->tree_nodes[curr]->bag.end())
		{
			//If all neighbors are in subtree, add to extraneous.
			extra = true;
			n1=this->G->get_node(*bagit);
			nbr_list=n1->get_nbrs_ptr();
			for(nbrit = nbr_list->begin();
				nbrit != nbr_list->end(); 
				++nbrit)
			{
				if(!(in_subtree[*nbrit]))
					extra = false;
			}
			if(extra){
				print_message(1, "Marking %d as extra.\n", *bagit);
				extraneous.push_back(*bagit);
			}
			++bagit;
		}
		//remove the extraneous vertices found in this bag from remainder of tree

		for(int i = 0; i < (int)this->tree_nodes.size(); i++)
		{
			if(this->tree_nodes[i] != NULL && subnode[i] != true)
			{
				for(eit = extraneous.begin(); eit != extraneous.end(); ++eit)
				{
					if(i == 0) print_message(1, "%d is extraneous at %d\n", *eit, curr);
					if(1)
						s = this->tree_nodes[i]->bag.size(); 
					this->tree_nodes[i]->bag.remove(*eit); 	 
					if(1) 
						if((int)this->tree_nodes[i]->bag.size()< s)
							print_message(0, "Removing %d from %d\n", *eit, i);
				}
			}//end if node not in subtree at curr
		}//end of loop over all tree nodes to remove extraneous
	}//end of processing post-order walk

	this->remove_subset_bags(); 
}
