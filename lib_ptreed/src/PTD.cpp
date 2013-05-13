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

#include "PTD.h" 
#ifdef VTRACE 
#include "vt_user.h" 
#endif 
#include <sys/types.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace madness;

PTD::PTD(World &world, Graph::VertexWeightedGraph g) :
    WorldObject<PTD>(world), G(g), H(g), xtime(MPI_Wtime()) 
{
    r = world.rank();
    s = world.size();
}

PTD::~PTD() {}

int PTD::create_bag(int v, int i, list<string> &lstring) 
{
#ifdef VTRACE
    VT_TRACER("create_bag");
#endif
    list<int> neighbors;
    int minpos = INT_MAX;
    ostringstream os;
    
    //Graph::GraphEOUtil eo;
    //eo.find_forward_neighbors(&G, v, &ordering, i, &neighbors,
    //&minpos);

    // We can not use Graph::GraphEOUtil::find_forward_neighbors since
    // it is not thread safe and modifying G's adj_vec. In order to
    //use this function one have to make copies of "G" which is
    //expensive for larger graph. Hence we are using new
    //implementation for finding forward neighbors.

        Graph::Node *n = G.get_node(v);
        list<int> *nl = n->get_nbrs_ptr();
        list<int>::iterator itn = nl->begin();
        int j;

        for (; itn != nl->end(); ++itn)
        {
            j = *itn;
            if (position[j] > i)
            {
                neighbors.push_back(j);
                if (position[j] < minpos)
                    minpos = position[j];
            }
        }
        
    neighbors.push_back(v);
    neighbors.sort();
 
    list<int>::iterator it = neighbors.begin();
    os << "e " << ordering[minpos] + 1 << " " << v + 1 << endl;
    os << "B " << v + 1 << " " << neighbors.size() << " ";
    for ( ; it != neighbors.end(); ++it)
        os << *it + 1 << " ";

    os << endl;
    lstring.push_back(os.str());
    
    return 1;
}


int PTD::owner(int i) 
{
	return (i % size());
}

vector<int>* PTD::get_elim_order() 
{
	return &ordering;
}

void PTD::set_elim_order(vector<int> &order)
{
    ordering = order;
}

int PTD::compute_elim_order(int order_type, bool t, const char *s) 
{
#ifdef VTRACE
    VT_TRACER("comp_eorder");
#endif

    this->ordering.resize(G.get_num_nodes());
    this->position.resize(G.get_num_nodes());

    if (t)
        read_ordering_file(s, &ordering);
    else
        eoutil.find_elimination_ordering(&G, &this->ordering, order_type, false);
    
    int os = ordering.size();

    for (int i = 0; i < os; i++)
        position[ordering[i]] = i;

	return 1;
}

int PTD::triangulate() 
{
#ifdef VTRACE
    VT_TRACER("trangulate");
#endif
	int width = 0;

#if HAS_METIS
	width = eoutil.METIS_triangulate(&G, &(this->ordering));
#else
	width = eoutil.triangulate(&G, &(this->ordering));
#endif

    this->set_width(width);
	return width;
}


int PTD::get_width() const {
	return this->width;
}

void PTD::set_width(int w)
{
	this->width = w;
}

//Function that executed by different threads. 
static void *thread_process_bag(void *a) 
{
    TPTD *t = reinterpret_cast<TPTD*>(a);
    PTD *p = t->get_ptd();
    vector<int> range = t->get_range_vector();
    vector<int> *ordering = p->get_elim_order();
    list<string> lstring;
    ostringstream s;

    for (int i = range[0]; i < range[1]; i++)
        p->create_bag(ordering->at(i), i, lstring);

    ScopedMutex<MutexFair>m(p->tlock);
    while (!lstring.empty() && (p->fw.good()))
    {
        p->fw << lstring.front();
        lstring.pop_front();
    }
    p->fw.flush();
    
}

void PTD::process_bag() 
{
    int os = ordering.size();
    os --;

    int ws = size();
    int wr = rank();
    int xr = os / ws;
    int xm = os % ws;
    int startpos = 0;

    int nthreads = 4;
    nthreads = atoi(getenv("PTD_NUM_THREADS"));
    if (!(nthreads > 0))
        fatal_error("error in PTD_NUM_THREADS env variable: %d\n", nthreads);

    ostringstream s;
    s << this->get_file_name();
    s << "." << wr;
    fw.open(s.str().c_str(), fstream::out);

    if (wr < xm)
    {
        xr ++;
        startpos = wr*xr;
    }
    else
    {
        startpos = wr*xr + xm;
    }
    
    //DEBUG("rank: %d xr: %d xm: %d startpos: %d endpos: %d \n", wr, xr, xm, startpos, startpos + xr);
    // Now divide work for threads.
    int tr = xr / nthreads;
    int tm = xr % nthreads;
    pthread_t pid[nthreads];
    int id[nthreads];
    vector<TPTD *> tptd(nthreads);
    int tstart = 0;

    // creating threads based on PTD_NUM_THREADDS
    for (int i = 0; i < nthreads; i++)
    {
        pthread_t t;
        vector<int> range;
        int itr = tr;
        if (i < tm)
        {
            itr++;
            tstart = startpos + i*itr;
        }
        else
        {
            tstart = startpos + i*itr + tm;
        }

        range.push_back(tstart);
        range.push_back(tstart + itr);

        id[i] = i;
        tptd[i] = new TPTD(this, range, id[i]);

        int rc = pthread_create(&pid[i], NULL, &thread_process_bag, tptd[i]);
        if (rc)
        {
            FERROR("error creating threads: %d\n", rc);
        }
    }


    // Join all threads before returning from this function
    for (int i = 0; i < nthreads; i++)
    {
        int rc = pthread_join(pid[i], NULL);
        if (rc)
        {
            FERROR("error while joining pthreads: %d\n", rc);
        }
         
        // delete all TPTD classes which created before pthread_create
        delete tptd[i];
    }
    fw.close();
}

int PTD::size() 
{
    return s;
}

int PTD::rank() 
{
    return r;
}

string PTD::get_random_name(char *x)
{
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";

    int len = 8;
    string s;
    ostringstream os;
    int seed = get_int(x);
    srand(seed);

    string fff(x);
    size_t found = fff.rfind("/");
    if (found != string::npos)
    {
        found ++;
        string::iterator it = fff.begin();
        fff.erase(it, it + found);
    }
    os << fff;
    os << "." << size() << ".";

    for (int i = 0; i < len; ++i)
        s.push_back(alphanum[rand() % (sizeof(alphanum) - 1)]);

    os << s;

    return os.str();
}

void PTD::set_file_name(string s)
{
    this->tmp_filename = s;
}

string PTD::get_file_name() const
{
    return this->tmp_filename;
}

int PTD::write_elim_order(char *filename)
{

    int os = (this->ordering).size();
    FILE *ef = fopen(filename, "w");
    if (!ef)
      {
        FERROR("can not open %s for writing\n", filename);
        return(-1);
      }
    for (int i = 0; i < os; i++)
      {
        fprintf (ef, "%d\n", (this->ordering)[i] + 1);
      }
    fclose (ef);
}


void PTD::merge_files(char *x)
{
    ostringstream os;
    ostringstream bs;

    string s(x);
    int nthreads = atoi(getenv("PTD_NUM_THREADS"));
    int si = size();
    string tmp = this->get_file_name();

    //BDS - this should respect the file name given with the -w option! 
    os << s << "." << size() << ".tree";
    
    bs << tmp << ".bag";

    this->set_tree_name (os.str());
    //DEBUG("tree: %s\n", os.str().c_str());
    ofstream outfile(os.str().c_str());
    ofstream bagfile(bs.str().c_str());
    string line;

    outfile << "p treed " << G.get_num_nodes() << " "<< G.get_num_nodes() - 1 << endl;
    for (int i = 0; i < si; i++)
    {
        ostringstream tf;
        tf << tmp << "." << i;
        ifstream infile(tf.str().c_str());
        if (!infile)
        {
            cerr << "Error reading file" << tf.str().c_str() << endl;
            fatal_error("Error reading file: %s\n", tf.str().c_str());
        }
        
        while (infile.good())
        {
            getline(infile, line);
            if (line[0] == 'e')
                outfile << line << endl;
            else if (line[0] == 'B')
                bagfile << line << endl;
        }
        infile.close();
        remove(tf.str().c_str());
    }

    ifstream inbagfile(bs.str().c_str());

    outfile << "B " << ordering.back() + 1 << " " << 1 << " " << ordering.back() + 1 << endl;
    while (inbagfile.good())
    {
        getline(inbagfile, line);
        if (line.size() > 1)
            outfile << line << endl;
    }
    inbagfile.close();
    remove(bs.str().c_str());
    outfile.close();

    // Refine generated tree decomposition
    refine_td(os.str());
}

void PTD::refine_td(string x)
{
    //DEBUG("refine tree: %s\n", x.c_str());
    TDTree *T = new TDTree();
    T->width = T->read_DIMACS_file(x.c_str());
    T->root(ordering.back());
    //DEBUG("width: %d\n", T->width);
    int os = ordering.size();
    for (int i = os - 2; i >=0; i--)
    {
        TDTreeNode *n = T->tree_nodes[ordering[i]];
        if (n->bag.size() > 1)
        {
            TDTreeNode *p = T->tree_nodes[n->adj.front()];
            n->bag.remove(ordering[i]);
            if (p->bag == n->bag)
            {
                //DEBUG("p: %d collapse into n: %d\n", p->id, n->id);
                list<int>::iterator it = p->adj.begin();
                n->adj.remove(p->id);
                it = p->adj.begin();
                for (; it != p->adj.end(); ++it)
                {
                    if (*it != n->id && *it != p->id)
                    {
                        TDTreeNode *t = T->tree_nodes[*it];
                        // Check for already deleted nodes.
                        if (!t->adj.empty())
                        {
                            // add current node as parent or neighbor of t
                            if(t->adj.front() == p->id){
                                t->adj.push_front(n->id);
                            } else {
                                t->adj.push_back(n->id);
                            }
                            // remove old parent node from t
                            t->adj.remove(p->id);
                            // add t as a new non-parent neighbor of n
                            n->adj.push_back(t->id);
                        }
                    }
                }

                // update parent of n
                if(p->id == p->adj.front())
                {  // p was the old tree root
                        n->adj.push_front(n->id);
                } 
                else 
                { // p was not the old tree root
                    n->adj.remove(p->adj.front());
                    n->adj.push_front(p->adj.front());
                }
                // Deleting parent node after merging bags.
                p->adj.clear();
                p->bag.clear();
            }
            n->bag.push_back(ordering[i]);
            n->bag.sort();
        }
    }
    
    int m = 0;
    vector<int> perm;
    for (int i = 0; i < os; i++)
    {
        TDTreeNode *n = T->tree_nodes[i];
        if (n->adj.size() > 0)
        {
            perm.push_back(m);
            m++;
        }
    }

    DEBUG("m: %d perm: %d\n", m, perm.size());
    for (int i = 0, j = 0; i < os; i++)
    {
        TDTreeNode *n = T->tree_nodes[i];
        if (n->adj.size() > 0)
            n->id = j++;
    }

    srand(time(NULL));
    random_shuffle(perm.begin(), perm.end());

    ofstream outfile(x.c_str());

	outfile << "p treed " << m << " " << m - 1 << endl;
    for (int i = os - 1; i >= 0; i--)
    {
        TDTreeNode *n = T->tree_nodes[ordering[i]];
        if (n->adj.size() > 0)
        {
            list<int>::iterator it = n->adj.begin();
            for (; it != n->adj.end(); ++it)
            {
                if (perm[n->id] < perm[T->tree_nodes[*it]->id])
                    outfile << "e " << perm[n->id] + 1 << " " << perm[T->tree_nodes[*it]->id] + 1 << endl; 
            }
        }
    }

    for (int i = 0; i < os; i++)
    {
        TDTreeNode *n = T->tree_nodes[i];
        if (n->adj.size() > 0)
        {
            list<int>::iterator it = n->bag.begin();
            outfile << "B " << perm[n->id] + 1 << " " << (int)T->tree_nodes[i]->bag.size() << " ";
            for (; it != n->bag.end(); ++it)
            {
                outfile << *it + 1 << " ";
            }
            outfile << endl;
        }
    }
    //DEBUG("finish refinement\n");
}

int PTD::get_int(char *x)
{
    if (!x)
        return 0;

    int len = strlen(x);
    int value = 0;
    for (int i = 0; i < len; i++)
        value += x[i];

    return value;
}

void PTD::create_tree_decomposition(char *filename, TDTree *&T)
{
    ostringstream os;
    string s(filename);
    os << s << "." << size() << ".tree";

    T = new TDTree(&G);
    this->set_file_name(this->get_random_name(filename).c_str());
    this->process_bag();
    world.gop.fence();
    if (this->rank() == 0)
        this->merge_files(filename);
    world.gop.fence();
    //DEBUG("read: %s\n", os.str().c_str());
    T->width = T->read_DIMACS_file(os.str().c_str());

	// Determine how big the masks need to be based on the tree width
	// CSG changing July 6!
	// Width is 1 less than the largest bag! Then we need an extra word
	// to handle a shift of BIGINT_WORD_SIZE!!
	T->num_mask_words = (int) ceil(
			((double) T->width + 2) / BIGINT_WORD_SIZE);
    // Root the tree
    T->root(0);

	// Reset (*T)'s graph is the original, non-triangulated graph!
	T->G = &H;
	DEBUG("width=%d; num_words=%d\n", T->width, T->num_mask_words);
}

void PTD::set_tree_name (string t)
{
    this->tree_name = t;
}

string PTD::get_tree_name() const
{
    return tree_name;
}

int PTD::parmetis_elim_order(int order_type)
{
#ifdef HAS_PARMETIS
    this->ordering.resize(G.get_num_nodes(), -1);
    this->position.resize(G.get_num_nodes(), -1);
    world.gop.fence();
    if (order_type)
        eoutil.parmetis_elimination_ordering(&G, this->ordering, order_type, false, MPI_COMM_WORLD);    
    else
        eoutil.parmetis_elimination_ordering(&G, this->ordering, -1, false, MPI_COMM_WORLD);    
    int os = ordering.size();
    for (int i = 0; i < os; i++)
        position[this->ordering[i]] = i;
#endif
	return 1;
}











