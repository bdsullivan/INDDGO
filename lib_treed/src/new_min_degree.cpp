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

/**
 * Eliminates the vertex and returns the minimum degree of those nodes whose
 * degree changes when adding the fill edges.  Sets affected[u]=true for all vertices
 * u whose degree changes, leaves affected[u] untouched otherwise.  The affected[] array is
 * assumed to have size of at least this->capacity. Appends the
 * corrected degrees of all affected nodes to the updates list. This list is sorted by degree before returning
 */
int Graph::eliminate_vertex(int v, list<int_int> *updates, bool *affected, bool *ordered)
{
     if(!this->canonical)
         print_message(0,"%s:  Graph must be in canonical form\n",__FUNCTION__);

     if(v<0||v>=this->capacity)
         fatal_error("%s:  eliminate_vertex() called with vertex %d and there are "
         "%d possible nodes\n",__FUNCTION__,v,this->capacity);

     if(this->degree[v]==0)
         // Nothing to do
         return TD_INFINITY;

     int i,j;
     int_int old_nn, nn;
     list<int_int>::iterator aa;
     // nn.p2 is index, nn.p1 is degree

     // v is affected!
     affected[v]=true;

     // First, remove all elements of Nbrs(v) from the updates list since we are going
     // to be modifying these entries
     print_message(1,"Update list begins with %d elements\n",updates->size());
     this->fill_adj_vec(v);

     if(this->adj_vec[v])
         fatal_error("Self edge %d-%d in elim. vertex!\n",v,v);

     aa=updates->begin();
     while(aa!=updates->end())
     {
         print_message(1, "Checking (%d,%d)\n",aa->p2,aa->p1);
         if(this->adj_vec[aa->p2] || aa->p2==v || ordered[aa->p2])
         {
             print_message(1,"Erasing (%d,%d)\n",aa->p2,aa->p1);
             // this vertex will be updated below - remove it here!
             aa=updates->erase(aa);
         }
         else
         {
             if(aa!=updates->end())
                 ++aa;
         }
     }
     print_message(1,"Update list now has %d elements\n",updates->size());

     // Handle special case where v has degree 1
     if(this->degree[v]==1)
     {
         print_message(1,"Degree[%d]=1\n",v);
         j=this->nodes[v].nbrs.front();
         // j is the one neighbor of v

         int orig_size=this->nodes[j].nbrs.size();
         this->nodes[j].nbrs.remove(v);
         if(this->nodes[j].nbrs.size() != orig_size-1)
             fatal_error("Didn't succeed in removing vertex v=%d from %d's nlist\n",v,j);

         this->nodes[v].nbrs.pop_front();
         this->num_edges--;
         this->degree[v]=0;
         this->degree[j]--;
         nn.p2=j;
         nn.p1=this->degree[j];
         updates->push_back(nn);
         affected[j]=true;

         updates->sort();
         return this->degree[j];
     }

     int min_deg=TD_INFINITY;
     list<int>::iterator kk;
     vector<int>::iterator ii,vv;
     vector<int> nadj_nodes(this->capacity);

     for(kk=this->nodes[v].nbrs.begin();kk!=this->nodes[v].nbrs.end();++kk)
     {   
         if(ordered[*kk])
             fatal_error("Encountered ordered node %d in elim vertex\n",*kk);

         // Compute the list of nodes that are in v's nbr list but not in *kk's 
         this->nodes[v].nbrs.sort();
         this->nodes[*kk].nbrs.sort();
         print_message(1,"Calling set_diff(%d[%d],%d[%d])\n",v,this->nodes[v].nbrs.size(),
             *kk,this->nodes[*kk].nbrs.size());
         vv=set_difference(this->nodes[v].nbrs.begin(),this->nodes[v].nbrs.end(),
             this->nodes[*kk].nbrs.begin(),this->nodes[*kk].nbrs.end(),nadj_nodes.begin());
         print_message(1,"set_diff(%d,%d) contains %d elements\n",v,*kk,vv-nadj_nodes.begin());

         if(nadj_nodes.size()>0)
         {
             for(ii=nadj_nodes.begin();ii!=vv;++ii)
             {
                 //if(*kk==*ii)
                 //    print_message(0,"WARNING!! Adding self loop %d-%d in elim vertex\n",*kk,*ii);

                 if(*kk!=*ii)
                 {
                     // Don't add self edges!
                     this->nodes[*kk].nbrs.push_back(*ii);
                     this->nodes[*ii].nbrs.push_back(*kk);
                     this->degree[*kk]++;
                     this->degree[*ii]++;
                     this->num_edges++;
                 }
             }
             // Need to make sure nbr lists are always sorted
             this->nodes[*kk].nbrs.sort();
             // this->nodes[*kk].nbrs.unique(); // shouldn't be necessary...
         }
         // Now delete the v-*kk edges!
         int orig_size=this->nodes[*kk].nbrs.size();
         this->nodes[*kk].nbrs.remove(v);
         if(this->nodes[*kk].nbrs.size() != orig_size-1)
             fatal_error("Didn't succeed in removing vertex v=%d from %d's nlist\n",v,*kk);
     }

     // Now populate the list with the new degrees
     for(kk=this->nodes[v].nbrs.begin();kk!=this->nodes[v].nbrs.end();++kk)
     {
         // Insert check for v here!
         if(*kk==v)
             fatal_error("Encountered kk=v when adding to updates list!!\n");
         // Node *kk is affected
         affected[*kk]=true;
         nn.p1=this->degree[*kk];
         nn.p2=*kk;
         updates->push_back(nn);
         if(this->degree[*kk]<min_deg)
             min_deg=this->degree[*kk];
     }
     // Sort updates by degree
     updates->sort();

     // v now has no neighbors
     this->nodes[v].nbrs.clear();

     // Return the lowest degree of the affected vertices  
     return min_deg;
}


/**
 * Finds an elimination ordering using the min_degree heuristic. Experimental code
 * This is currently slower than other implementation, but I think there is some room 
 * for improvement.
 */
int Graph::find_min_degree_ordering(vector<int> *ordering, int start_v)
{
    if(!this->canonical)
        print_message(0,"%s:  Graph must be in canonical form\n",__FUNCTION__);

    int i, j, m, min_deg, min_pos=-1;
  
    // Make sure ordering[] is big enough for the graph
    if((int)ordering->size()<this->capacity)
        fatal_error("%s:  The ordering[] array is not big enough: %d<%d\n",
        __FUNCTION__,ordering->size(),this->capacity);

    bool *affected=new bool[this->capacity];
    bool *ordered=new bool[this->capacity];
    list<int_int> *min_list, *updates;
    min_list = new list<int_int>;
    updates=new list<int_int>;
    int_int jj,kk,nn;
    list<int_int>::iterator aa,bb;
    for(i=0;i<this->capacity;i++)
        ordered[i]=false;

    // Make a copy of this graph 
    Graph G=*this;

    // Make sure the degrees are accurate
    G.recompute_degrees();

    // Begin with the provided start_v ("best" ordering probably comes
    // by selecting start_v to be a vertex of minimum degree)
    i=0;
    print_message(1,"Setting ordering[0]=%d\n",start_v);
    ordering->at(i)=start_v;
    ordered[start_v]=true;
    G.eliminate_vertex(start_v, min_list, affected, ordered);
    G.degree[start_v]=TD_INFINITY;
    int num_elimed;

    // Create the initial sorted min_list
    min_list->clear();
    for(m=0;m<G.capacity;m++)
    {
        affected[m]=false;
        // We could do this in random order here as well...
        // Make sure location m is valid
        if(G.nodes[m].label!=UNDEFINED && ordered[m]==false)
        {
            nn.p1=G.degree[m];
            nn.p2=m;
            min_list->push_back(nn);
        }
    }
    min_list->sort();
    updates->clear();

    print_message(MIN_DEG_DEBUG,"Initial min_list has %d entries\n",min_list->size());

    i=1;
    while(i < this->num_nodes)
    {
        print_message(MIN_DEG_DEBUG,"i=%d of %d\n",i,G.num_nodes);
        while(!min_list->empty())
        {
            // Grab the front of the list
            jj=min_list->front();
            if(updates->size()==0)
            {
                // This is just to handle the very beginning when updates is empty
                kk.p1=TD_INFINITY;
                kk.p2=UNDEFINED;
            }
            else
                kk=updates->front();

            print_message(MIN_DEG_DEBUG,"i=%d: min_list=(%d,%d); updates=(%d,%d)\n",i,jj.p2,jj.p1,kk.p2,kk.p1);

            while(jj.p1 <= kk.p1 && !min_list->empty())
            {
                print_message(MIN_DEG_DEBUG,"i=%d: Ordering %d,%d,%d from min_list (updates is %d,%d)\n",i,jj.p2,jj.p1,affected[jj.p2],
                    kk.p2,kk.p1);
                // Sanity check
                if(ordered[jj.p2])
                    fatal_error("Node %d is already ordered in min_list!!\n",jj.p2);
                
                // Next min element is in min_list - order it and remove all the affected nodes from
                // the front of min_list
                // Eliminate the vertex and add it to the ordering
                G.eliminate_vertex(jj.p2, updates, affected, ordered);
                ordering->at(i)=jj.p2;
                ordered[jj.p2]=true;
                i++;

                // Now remove all the affected nodes from the front of min_list
                print_message(MIN_DEG_DEBUG,"min_list has %d entries before removals\n",min_list->size());
                aa=min_list->begin();
                while(affected[aa->p2]==true || ordered[aa->p2])
                {
                    min_list->pop_front();                    
                    if(min_list->empty())
                        break;
                    aa=min_list->begin();
                }                
                print_message(MIN_DEG_DEBUG,"min_list has %d entries after removals\n\n",min_list->size());

                if(min_list->empty())
                    break;
                // Update the value of jj and kk
                jj=min_list->front();
                kk=updates->front();

            }
            
            while(kk.p1 <= jj.p1 && !min_list->empty())
            {
                print_message(MIN_DEG_DEBUG,"i=%d: Ordering %d,%d,%d from updates (min_list is %d,%d)\n",i,kk.p2,kk.p1,
                     affected[kk.p2],jj.p2,jj.p1);

                  // Sanity check
                if(ordered[kk.p2])
                    fatal_error("Node %d is already ordered in updates!!\n",kk.p2);
                
                // Next min degree element is in updates list - order it and remove all the affected nodes from
                // the front of min_list
                G.eliminate_vertex(kk.p2, updates, affected, ordered);
                ordering->at(i)=kk.p2;
                ordered[kk.p2]=true;
                i++;

                /// CSGCSGCSG
                // Need to remove kk from the updates list - otherwise we might encounter it again!!

               // Now remove all the affected nodes from the front of min_list
                print_message(MIN_DEG_DEBUG,"min_list has %d entries before removals\n", min_list->size());
                aa=min_list->begin();
                while(affected[aa->p2]==true || ordered[aa->p2])
                {
                    min_list->pop_front();                    
                    if(min_list->empty())
                        break;
                    aa=min_list->begin();
                }                
                print_message(MIN_DEG_DEBUG,"min_list has %d entries after removals\n\n", min_list->size());

                // Remove any ordered nodes from the beginning of updates list
                aa=updates->begin();
                while(ordered[aa->p2])
                {
                    updates->pop_front();
                    aa=updates->begin();
                }
                
                if(min_list->empty())
                    break;
                // Update the value of jj and kk
                jj=min_list->front();
                kk=updates->front();
            }
        }
        // min_list is now empty - so every node has been hit since we last updated
        // Switch updates and min_list;
#if 1
        list<int_int> *temp;
        temp=min_list;
        min_list=updates;
        updates=temp;
#else

        // Copy the updates list to min_lists, set affected to false, clear updates and continue
        print_message(MIN_DEG_DEBUG,"min_list is empty! restarting...\n");
        for(aa=updates->begin();aa!=updates->end();++aa)
            if(ordered[aa->p2]==false)
                min_list->push_back(*aa);
#endif

        // CSG - think about this more - do we need to consider only those nodes in the
        // new min_list?
        updates->clear();
        for(m=0;m<this->capacity;m++)
        {
            if(ordered[m]==false)
                affected[m]=false;
            else
                affected[m]=true;
        }

    }

    print_message(MIN_DEG_DEBUG,"Ordering done with i=%d\n",i);
    print(MIN_DEG_DEBUG,*ordering);

    delete [] affected;
    delete [] ordered;

    return 1;
}
