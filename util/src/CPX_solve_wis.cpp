#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ilcplex/cplex.h>
#include <glpk.h>

#include "GraphDecomposition.h"

void check_CPX_status(char *func, int status){
    if(status){
        fprintf (stderr, "CPX nonzero status after %s: %d\n",func,status);
        exit(-1);
    }
};

void write_ind_set_model(const char *DIMACS_file, const char *model_file,
                         Graph::VertexWeightedGraph *G){
    int i, j;

    FILE *out = fopen(model_file, "w");
    if(!out){
        fatal_error("Can't open %s for writing\n", model_file);
    }

    // Write out the model
    fprintf(out,
            "# Automatically generated model file for max independent set\n"
            "# Assumes 1-based node names with no holes!!\n"
            "# Original graph file is %s\n"
            "# Model section\n\n"
            "param num_nodes ;\n"
            "set Nodes := 1..num_nodes ;\n"
            "set E within Nodes cross Nodes ;\n"
            "param w{Nodes};\n"
            "var x{i in Nodes} binary;\n\n"
            // Change to minimization for archaic MPS reasons
            "minimize Cost: sum{i in Nodes} -x[i]*w[i];\n\n"
            "subject to Independent{i in Nodes, j in Nodes: (i,j) in E}:\n"
            "x[i]+x[j]<=1;\n\n"
            "solve;\n\n", DIMACS_file);
    /*
       fprintf(
            out,
            "printf \" \\nMaximum Independent Set: %%d\\n\", sum{i in Nodes} x[i]*w[i] > \"%s\";\n",
            sol_file);
       fprintf(out, "for { i in Nodes }\n"
            "{\n"
            "   printf{0..0: x[i]!=0} \"%%d (%%3.1f)\\n\", i, w[i] >> \"%s\";\n",
            sol_file);
       fprintf(out, "}\n\n");*/
    fflush(out);
    fprintf(stderr,"Writing data section\n");

    fprintf(out, "# Data Section\n\ndata;\n\n");
    fflush(out);
    // Write out the problem data
    fprintf(out, "param num_nodes:=%d;\n", G->get_num_nodes());
    fflush(out);
    fprintf(out, "param w:=\n");
    fflush(out);
    vector<Graph::Node> nodes = G->get_nodes();
    vector<int> weight = G->get_weight();

    for(i = 0; i < G->get_num_nodes(); i++){
        fprintf(out, "%d %d\n",nodes[i].get_label(), weight[i]);
    }
    fprintf(out, ";\n");
    fprintf(out, "set E:=\n");
    Graph::GraphProperties properties;

    vector<int> adj_vec;
    for(i = 0; i < G->get_num_nodes(); i++){
        // This assumes 1-based labeling!
        int current_key = properties.fill_adj_vec(G, i);
        //int current_key=G->fill_adj_vec(i);
        adj_vec = G->get_adj_vec();
        for(j = i + 1; j < G->get_num_nodes(); j++){
            if(adj_vec[j] == current_key){
                fprintf(out, "(%d,%d)\n", i + 1, j + 1);
            }
        }
        fprintf(out, "\n");
    }
    fprintf(out, ";\nend;\n");
    fclose(out);
};

void read_DIMACS(const char *DIMACS_file, Graph::VertexWeightedGraph *G){
    int i;
    Graph::GraphCreatorFile creator;
    creator.set_file_name(DIMACS_file);
    creator.set_graph_type("DIMACS");

    G = creator.create_weighted_mutable_graph();
    if(!G){
        fprintf(stderr,"Could not create graph for DIMACS file %s\n",DIMACS_file);
        exit(-1);
    }

    // Add weights of 1 to all vertices if we didn't load any from the file
    bool has_weights = false;
    vector<int> weight = (G)->get_weight();
    for(i = 0; i < (G)->get_capacity(); i++){
        if(weight[i] != 0){
            has_weights = true;
            break;
        }
    }

    if(!has_weights){
        for(i = 0; i < (G)->get_capacity(); i++){
            weight[i] = 1;
        }
    }

    fprintf(stderr,"Graph loaded with %d nodes, %d edges\n",G->get_num_nodes(),
            G->get_num_edges());
};

void create_mps_file(char *gmpl_file, char *mps_file, bool verbose){
    // Read in the model
    glp_prob *lp;
    glp_tran *tran;
    if(!verbose){
        glp_term_out(GLP_OFF);
    }
    lp = glp_create_prob();
    tran = glp_mpl_alloc_wksp();
    fprintf(stderr,"reading model\n");
    int ret = glp_mpl_read_model(tran, gmpl_file, 0);
    if(ret){
        fprintf(stderr, "Error translating model\n");
        exit(-1);
    }
    ret = glp_mpl_generate(tran, NULL);
    if(ret){
        fprintf(stderr, "Error generating model\n");
        exit(-1);
    }

    glp_mpl_build_prob(tran, lp);
    ret = glp_write_mps(lp, GLP_MPS_FILE, NULL, mps_file);
    if(ret){
        fprintf(stderr, "Error on writing MPS file\n");
        exit(-1);
    }
    glp_mpl_free_wksp(tran);
    glp_delete_prob(lp);
    fprintf(stderr,"GLPK created MPS file %s\n",mps_file);
} // create_mps_file

void usage(char *str){
    fprintf(stderr,"%s <DIMACS_file> <max_seconds> [-n num_threads]\n"
                   "\t Uses CPLEX with num_threads to solve WIS with max_seconds time limit\n",str);
    exit(-1);
}

int main(int   argc,char *argv[]){
    if(argc < 3){
        usage(argv[0]);
    }
    char *DIMACS_file = argv[1];
    // set time limit
    int max_secs = atoi(argv[2]);
    int num_threads = -1;
    if(argc >= 4){
        for(int i = 3; i < argc; i++){
            if(strcmp(argv[i],"-n") == 0){
                num_threads = atoi(argv[i + 1]);
            }
        }
    }

    // Load the graph from DIMACS_file
    Graph::VertexWeightedGraph *G;
    Graph::GraphCreatorFile creator;
    creator.set_file_name(DIMACS_file);
    creator.set_graph_type("DIMACS");
    G = creator.create_weighted_mutable_graph();

    char mps_file[100], mod_file[100];

    int status;
    fprintf(stderr,"G has %d nodes, %d edges\n",G->get_num_nodes(),
            G->get_num_edges());

    // Write a GMPL model
    sprintf(mod_file,"%s.mod",DIMACS_file);
    write_ind_set_model(DIMACS_file, mod_file, G);
    // Write mps file
    sprintf(mps_file,"%s.mps",DIMACS_file);
    create_mps_file(mod_file,mps_file,false);

    CPXENVptr env = NULL;
    env = CPXopenCPLEX (&status);
    check_CPX_status("CPXopenCPLEX",status);

    CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);

    clock_t start = clock(),stop;
    CPXLPptr mip;
    mip = CPXcreateprob(env,&status,mps_file);
    status = CPXreadcopyprob (env, mip, mps_file, NULL);
    check_CPX_status("CPXreadcopyprob",status);
    CPXsetdblparam(env, CPX_PARAM_TILIM, max_secs);

    if(num_threads != -1){
        // Set # of threads
        CPXsetintparam(env,CPX_PARAM_THREADS,num_threads);
    }

    int num_cols = CPXgetnumcols(env,mip);
    int num_rows = CPXgetnumrows(env,mip);

    int mip_status = CPXmipopt(env,mip);
    double obj;
    CPXgetobjval(env,mip,&obj);
    int objval = (int)obj;
    stop = clock();
    fprintf(stdout,"%s %d %d %d %d %d %3.3f\n",DIMACS_file,num_cols, num_rows,mip_status,
            -objval,getHWmem(),(double)(stop - start) / CLOCKS_PER_SEC);
    fflush(stdout);
    exit(-1);
} // main

