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

#include "ORD.h" 
#include "metislib.h" 
#include "parmetis.h"

#ifdef VTRACE 
#include "vt_user.h" 
#endif 
#include "GraphCreatorFile.h"

ORD::ORD(MPI_Comm c, char *f) : comm(c) {
    filename.assign(f);
    cflag = true;
}

ORD::ORD(char *f) {
    filename.assign(f);
    cflag = true;
}

ORD::~ORD(){}

void ORD::buildGraph() {
    Graph::GraphCreatorFile creator(this->filename);
    Graph::GraphUtil util;
    Graph::GraphProperties properties;
    
    int wr;
    int size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &wr);

    G = creator.create_weighted_mutable_graph();
    properties.make_canonical(G);
    char compname[512];
    if (!properties.is_connected(G))
    {
        Graph::GraphReaderWriterFactory factory;
        Graph::VertexWeightedGraph *H;

        // Since we are going to proceed with processing largest
        // component , set a mark as not original and write largest
        // component into a file.
        this->set_original(false);

        vector<list<int> *> members;
        int comp;
        comp = util.find_all_components(G, &members);

        DEBUG("found %d components\n", comp);
        int i = 0;
        int gsize = 0;
        int index = 0;
        for (int j = 0; j < comp; j++)
        {
            i = members[j]->size();
            if (i > gsize)
            {
                gsize = i;
                index = j;
                sprintf(compname, "%s_%d_comp", filename.c_str(), gsize);
                DEBUG("new larger comp : %d\n", gsize);
            }
        }

        this->set_filename(compname);
        H = creator.create_component(G, members[index], true);
        delete G;
        G = H;

        // Change the label into 1-base DIMACS format
        int j = 0;
        int gs = G->get_num_nodes();
        Graph::Node *n;
        for (int i = 0 ; i < gs; i++)
        {
            n = G->get_node(i);
            n->set_label(i+1);
        }

    }
}

int ORD::find_elim_ordering() {
    int ws;
    int wr;

    char eoname[512];
    char eoname_other[512];

    // Get size and rank from the communicator
    MPI_Comm_size(comm, &ws);
    MPI_Comm_rank(comm, &wr);

    double xtime = MPI_Wtime();
    sprintf(eoname, "%s.order.%d", this->filename.c_str(), ws);
    sprintf(eoname_other, "%s.order_other.%d", this->filename.c_str(), ws);

    DEBUG("size: %d, rank %d \n", ws, wr);
    int n = G->get_num_nodes();
    int x = n/ws;
    int xm = n%ws;
    int i = 0;
    DEBUG("n: %d x: %d xm: %d \n", n, x, xm);

    vector<int> xadj;
    vector<int> adjncy;

    vector<int> vtxdist(ws + 1, 0);
    vector<int> sizes(2*ws,0);
    vector<int> ordering(x+1, 0);
    vector<int> recvcnt(ws, 0);
    vector<int> displ(ws, 0);

    int numflag = 0;




    int options[10];

    options[0] = 0;
    vtxdist[0] = 0;
    for (i = 1; i <= ws; i++)
    {
        vtxdist[i] = vtxdist[i - 1] + x;
        if (i <= xm)
            vtxdist[i]++;
    }

    // prepareing displacement and receive counts to use with MPI_Gatherv
    for (i = 0; i < ws; i++)
    {
        recvcnt[i] = x;
        if (i < xm)
            recvcnt[i] ++;

        if (i > 0)
            displ[i] += displ[i - 1] + recvcnt[i - 1];
    }

    DEBUG("range: %d, %d\n", vtxdist[wr], vtxdist[wr + 1]);
    int j = 0;
    xadj.push_back(0);
    for (i = vtxdist[wr]; i < vtxdist[wr + 1]; i++)
    {
        Graph::Node *no = G->get_node(i);
        list<int> *l = no->get_nbrs_ptr();
        list<int>::iterator it = l->begin();

        for (; it != l->end(); ++it)
        {
            adjncy.push_back(*it);
            j++;
        }
        xadj.push_back(j);
    }

    if (METIS_OK != ParMETIS_V3_NodeND(&vtxdist.front(), &xadj.front(), &adjncy.front(), &numflag, options, &ordering.front(), &sizes.front(), &comm))
    {
        FERROR("error occured while processing parmetis, aborting\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    DEBUG("output from ParMETIS\n");
    double parmet_time = MPI_Wtime() - xtime;

    vector<int> recvbuf;
    n = G->get_num_nodes();
    if (wr == 0)
    {
        recvbuf = vector<int>(n, 0);
    }

    if (MPI_SUCCESS !=
        MPI_Gatherv((void *)&ordering.front(), recvcnt[wr], MPI_INT,
                    (void *)&recvbuf.front(), &recvcnt.front(), &displ.front(), MPI_INT,
                    0, comm))
    {
        FERROR("MPI error occured at Gatherv, Abort!\n");
        MPI_Abort(comm, -1);
    }

    vector<int> eo(n, 0);
    if (wr == 0)
    {
        for (int i = 0; i < n; i++)
        {
            eo[recvbuf[i]] = i;
        }

        FILE *f = fopen(eoname_other, "w");
        for (int i = 0; i < n; i++)
            fprintf(f, "%d\n", eo[i] + 1);
        fclose(f);
        DEBUG("ParMetis NodeND elimination ordering is in : %s\n", eoname_other);
    }

    ordering.clear();
    ordering.resize(recvcnt[wr], 0);

    if (MPI_SUCCESS !=
        MPI_Scatterv ((void *)&eo.front(), &recvcnt.front(), &displ.front(), MPI_INT,
                      (void *)&ordering.front(), recvcnt[wr], MPI_INT,
                      0, comm))
    {
        FERROR("MPI error occured at Scatterv, Abort! \n");
        MPI_Abort(comm, -1);
    }

    DEBUG("Scatterv completed\n");

    Graph::GraphCreatorFile gf;
    Graph::VertexWeightedGraph *wg;
    Graph::GraphEOUtil eoutil;
    Graph::GraphProperties prop;
    list<int>members(ordering.begin(), ordering.end());

    wg = gf.create_component(G, &members, false);
    prop.make_canonical(wg);

    vector<int> ord(recvcnt[wr], 0);
    vector<int> ordsend(recvcnt[wr, 0]);
    double xxtime = MPI_Wtime();
    eoutil.find_elimination_ordering(wg, &ord, GD_AMD, false);
    DEBUG("eo time : %f\n", MPI_Wtime() - xxtime);
    
    int sz = recvcnt[wr];

    for (int i = 0; i < sz; i++)
        ordsend[i] = wg->get_node(ord[i])->get_label();


    recvbuf.assign(n, -1);
    if (MPI_SUCCESS !=
        MPI_Gatherv((void *)&ordsend.front(), recvcnt[wr], MPI_INT, (void *)&recvbuf.front(), &recvcnt.front(), &displ.front(), MPI_INT, 0, comm))
    {
        FERROR("MPI error occured at Gatherv, Abort!\n");
        MPI_Abort(comm, -1);
    }

    double p_amd_time = MPI_Wtime() - xtime;
    if (wr == 0)
    {
        FILE *f = fopen(eoname, "w");
        for (int i = 0; i < n && wr == 0; i++)
            fprintf(f, "%d\n", recvbuf[i]);
        fclose(f);
    } 
    DEBUG("ordering is written into %s\n", eoname);
    DEBUG("%f,%f\n", parmet_time, p_amd_time);

    return 0;
}

void ORD::outputMetisFormat() {
    int n = G->get_num_nodes();
    cout << n << " " << G->get_num_edges() << endl;
    for (int i = 0; i < n; i++)
    {
        Graph::Node *no = G->get_node(i);
        list<int> *l = no->get_nbrs_ptr();
        list<int>::iterator it = l->begin();

        for (; it != l->end(); ++it)
        {
            cout << (*it)+1 << " ";
        }        
        cout << endl;
    }
}

void ORD::set_filename(char *fn) {
    this->filename.assign(fn);
}

const char *ORD::get_filename() const {
    return this->filename.c_str();
}


Graph::VertexWeightedGraph *ORD::get_graph() {
    return this->G;
}

bool ORD::is_original() {
    return this->cflag;
}

void ORD::set_original(bool x) {
    this->cflag = x;
}
