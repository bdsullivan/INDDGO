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
#include <fstream>
#include <limits>
#include <string>
#include <sstream>
#include <stdio.h>

/*BDS - added to get memory highwater mark*/
int parseLine(char *line){
    int i = strlen(line);
    while(*line < '0' || *line > '9'){
        line++;
    }
    line[i - 3] = '\0';
    i = atoi(line);
    return i;
}

int getHWmem(){
    #if WIN32
    return -1;
    #endif
    //Note: this value is in KB!
    FILE *file = fopen("/proc/self/status", "r");
    if(file == NULL){
        fprintf(stderr, "Error opening processor status file!\n");
        return -1;
    }
    else {
        int result = -1;
        char line[128];

        while(fgets(line, 128, file) != NULL){
            if(strncmp(line, "VmHWM:", 6) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
} // getHWmem

int getMemInfo(const char *x, int n){
    FILE *file = fopen("/proc/self/status", "r");
    if(file == NULL){
        fprintf(stderr, "Error opening processor status file!\n");
        return -1;
    }
    else {
        int result = -1;
        char line[128];
        while(fgets(line, 128, file) != NULL){
            if(strncmp(line, x, n) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
} // getMemInfo

//Reads in a color file with # for commented lines, #min #max optional range arguments, and a single space seperated line of
//floats (or ints), that will be placed in color_vector.  The max_color and min_color arguments are populated if the
//return value is true
bool read_color_file(const char input_file[], double & max_color, double & min_color, vector<double> & color_vector){
    ifstream k_core_input;
    k_core_input.open(input_file);

    if(!k_core_input.is_open()){
        std::cerr << "Error opening color file \n" << input_file;
        exit(EXIT_FAILURE);
    }
    bool max_flag, min_flag;
    max_flag = false;
    min_flag = false;
    max_color = 0;
    min_color = 0;

    while(k_core_input.good()){
        string line;
        getline(k_core_input,line);

        cerr << "Got color line: " << line << endl;
        if(line.substr(0,1) == "#"){
            if(strcmp(line.substr(1,3).c_str(),"max") == 0){
                size_t pos1 = line.find(" ");
                size_t pos2 = line.find_last_not_of("\t\f\v\n\r ");
                cout << "p1max: " << pos1 << " p2max: " << pos2 << "\n";

                max_color = atof(line.substr(pos1 + 1,pos2 + 1).c_str());
                max_flag = true;
            }

            if(strcmp(line.substr(1,3).c_str(),"min") == 0){
                size_t pos1 = line.find(" ");
                size_t pos2 = line.find_last_not_of("\t\f\v\n\r ");

                min_color = atof(line.substr(pos1 + 1,pos2 + 1).c_str());
                min_flag = true;
            }
        }
        else if(line.length() != 0){
            size_t pos1 = 0;
            size_t pos2 = 0;
            size_t end = line.find_last_not_of("\t\f\v\n\r ");
            //	int count = 0;
            while(pos1 < end){
                //count++;
                pos2 = line.find_first_of("\t\f\v\n\r ",pos1 + 1);

                color_vector.push_back(atof(line.substr(pos1,pos2).c_str()));

                pos1 = pos2;
            }
        }
    }

    k_core_input.close();

    return(min_flag && max_flag);
} // read_color_file

//Input a vector of node indices, a vector of node scores (values associated with each node), and a flag
//that indicates the statistic you desire to be calculated
//FLAG should be one of GD_XXX defined in Util.h

double get_statistics(const vector<int> & nodes, const vector<double> & scores, const int FLAG){
    double statistic = 0.0;
    const int size = nodes.size();

    if(FLAG == GD_STAT_MEAN){
        //compute average
        for(int i = 0; i < size; ++i){
            statistic += scores[nodes[i]];
        }

        statistic = statistic / double(size);
    }
    else if(FLAG == GD_STAT_STD){
        //compute standard deviation
        double mean = 0;
        for(int i = 0; i < size; ++i){
            mean += scores[nodes[i]];
        }
        mean = mean / double(size);

        double square_sum = 0;
        for(int i = 0; i < size; ++i){
            square_sum += pow((double)(scores[nodes[i]] - mean),2);
        }
        statistic = square_sum / double(size);

        statistic = sqrt(statistic);
    }
    else if(FLAG == GD_STAT_MED){
        //compute median - currently using builtin vector sort
        //Need to select sub-vector associated with nodes given in function call
        vector<double> sorted(size);
        for(int i = 0; i < size; ++i){
            sorted[i] = scores[nodes[i]];
        }

        sort(sorted.begin(), sorted.end());

        int middle = size / 2;
        if(size % 2 == 1){
            statistic = sorted[middle];
        }
        else {
            //Adjust downwards since we index from 0, not 1
            statistic = (sorted[middle] + sorted[middle - 1]) / 2;
        }
    }
    else if(FLAG == GD_STAT_COUNT){
        //Count number of non-negative scored vertices in bag.
        for(int i = 0; i < size; ++i){
            if(!(scores[nodes[i]] < 0)){
                statistic++;
            }
        }
    }
    else {
        fatal_error("Statistic must be one of GD_STAT_XXX from Util.h\n");
    }
    return statistic;
} // get_statistics

//Input a value between 0-1, get out an rgb heatmap value
void get_rgb_value(char rgb[], double val){
    int red, green, blue;
    double temp_color_value;

    //Red
    temp_color_value = 0;
    if(val > 0.875){
        temp_color_value = -4.0 * val + 4.5;
    }
    else if(val > .625){
        temp_color_value = 1.0;
    }
    else if(val > .375){
        temp_color_value = 4.0 * val - 1.5;
    }

    red = (int)(255 * temp_color_value);

    //green
    temp_color_value = 0;
    if(val > 0.875){
        temp_color_value = 0;
    }
    else if(val > .625){
        temp_color_value = -4.0 * val + 3.5;
    }
    else if(val > .375){
        temp_color_value = 1.0;
    }
    else if(val > .125){
        temp_color_value = 4.0 * val - 0.5;
    }

    green = (int)(255 * temp_color_value);

    //blue
    temp_color_value = 0;
    if(val > 0.625){
        temp_color_value = 0;
    }
    else if(val > .375){
        temp_color_value = -4.0 * val + 2.5;
    }
    else if(val > .125){
        temp_color_value = 1.0;
    }
    else if(val >= 0.0){
        temp_color_value = 4.0 * val + .5;
    }

    blue = (int)(255 * temp_color_value);

    sprintf(rgb,"#%02x%02x%02x",red,green,blue);
} // get_rgb_value

bool decreasing(int a, int b){
    // Use this to sort int's in decreasing order
    if(a < b){
        return false;
    }
    return true;
}

//checks if the double in i is < the double in j.
bool int_d_sort(int_double i, int_double j){
    return (i.d < j.d);
};

//checks if the int in i is < the int in j.
bool int_d_isort(int_double a, int_double b){
    return (a.i < b.i);
};

/*
 * Function to compare two int_int's (p1 first, then p2)
 */
bool operator <(const int_int& a, const int_int &b){
    if(a.p1 != b.p1){
        return a.p1 < b.p1;
    }
    else {
        return a.p2 <= b.p2;
    }
}

/**
 * Function to compare two int_int's for equality.
 */
bool operator ==(const int_int& a, const int_int &b){
    if((a.p1 == b.p1) && (a.p2 == b.p2) ){
        return true;
    }
    else {
        return false;
    }
}

// Assuming that scores[i] is a non-negative integer for all i in list, reorders list in ascending score order.
void score_sort(list<int> *mylist, vector<int> *scores){
    list<int_int> aux_list;
    list<int_int>::iterator ait;

    while(!(mylist->empty())){
        int_int curr;
        curr.p2 = mylist->front();
        mylist->pop_front();
        curr.p1 = (*scores)[curr.p2];
        aux_list.push_front(curr);
    }

    aux_list.sort();

    /*for(ait = aux_list.begin(); ait!=aux_list.end(); ++ait)
     * cout << *ait;
     * cout << "\n";
     */

    int_int curr;
    while(!aux_list.empty()){
        curr = aux_list.front();
        mylist->push_back(curr.p2);
        aux_list.pop_front();
    }

    return;
} // score_sort

/**
 * Function to extract the values of n and m from a DIMACS file so that
 * we can allocate correctly sized data structures.
 */
void get_DIMACS_dimensions(const char *DIMACS_file, int *n, int *m){
    char line[100], format[100], x;
    int retval;
    FILE *in;

    if( (in = fopen(DIMACS_file,"r")) == NULL){
        fatal_error("%s:  Error opening %s for reading\n",__FUNCTION__,DIMACS_file);
    }

    while(!feof(in)){
        retval = fscanf(in,"%2c",line);
        if(feof(in)){
            break;
        }

        if(retval != 1){
            fatal_error("%s:  fscanf read %d char's (expected to read 1!)\n",__FUNCTION__,retval);
        }

        switch(line[0])
        {
        case 'p':
            // This is the "problem line" - get n and m here
            // Make sure we don't already know n and m!
            retval = fscanf(in,"%s %d %d",format, n, m);

            print_message(1,"DIMACS: read n=%d, m=%d\n", *n, *m);

            // Simple error checking
            if((n <= 0) || (m <= 0) || (retval != 3) ){
                fatal_error("%s:  DIMACS read error - didn't understand problem line!\n"
                            "p %s %d %d\n",__FUNCTION__,format,*n,*m);
            }

            if((strncmp(format,"e",1) != 0) && (strncmp(format,"E",1) != 0) &&
               (strncmp(format,"edge",4) != 0) && (strncmp(format,"EDGE",4) != 0) ){
                fatal_error("%s:  DIMACS read error - problem line - FORMAT must be one of\n"
                            " e, E, edge, or EDGE\n",__FUNCTION__);
            }

            fclose(in);
            return;
            break;

        case 'c':
            // Comment line - skip and move on
            break;

        case 'n':
            // Node line - skip
            break;

        case 'e':
            // Edge line - skip
            break;
        case 'x':
            // in TSP files - skip
            break;

        default:
            // We encountered some character that we don't understand
            fatal_error("%s:  DIMACS read error - didn't understand character %c to start a line\n",
                        __FUNCTION__,line[0]);
            break;
        }        // end switch
                 // Advance to next line if there is format at the end of it
        while(!feof(in) && (x = getc(in)) != '\n'){
            ;
        }
    }    // end while

    // If we make it here, then we never found the dimension!
    fatal_error("%s:  Unable to find dimensions in DIMACS file %s\n",__FUNCTION__,DIMACS_file);
} // get_DIMACS_dimensions

ostream &operator<<(ostream &output, const VertexLabel &VL){
    output << "VertexLabel id=" << VL.id << "(" << VL.entries.size() << " entries)\n";
    list<int>::const_iterator i;
    for(i = VL.entries.begin(); i != VL.entries.end(); ++i){
        output << *i << " ";
    }
    output << endl;

    return output;
}

ostream &operator<<(ostream &output, const int_int &i){
    output << "(" << i.p1 << "," << i.p2 << ")";
    return output;
}

ostream &operator<<(ostream &output, const int_double &i){
    output << "(" << i.i << "," << i.d << ")";
    return output;
}

bool operator <(const VertexLabel& a, const VertexLabel &b){
    // Assumes that the entries lists are sorted in decreasing order!!
    list<int>::const_iterator ai,bi;

    ai = a.entries.begin();
    bi = b.entries.begin();
    while(ai != a.entries.end() && bi != b.entries.end()){
        if( (*ai) < (*bi)){
            return true;
        }
        else if( (*ai) > (*bi)){
            return false;
        }
        // entries are the same - move to next entry in the list
        ++ai;
        ++bi;
    }
    // We've reacheed the end of at least one of the lists
    if(a.entries.size() < b.entries.size()){
        // b is longer, so consider a<b
        return true;
    }
    else {
        return false;
    }
} // <

bool empty_intersection(list<int> *S1, list<int> *S2){
    // Assumes lists are sorted in increasing order!
    list<int>::iterator i1 = S1->begin();
    list<int>::iterator i2 = S2->begin();

    while(i1 != S1->end() && i2 != S2->end()){
        // Check for equality
        if(*i1 == *i2){
            return false;
        }
        else if(*i1 < *i2){
            i1++;
        }
        else {
            i2++;
        }
    }
    // Empty intersection if we get here
    return true;
} // empty_intersection

int read_separator_file(const char *sep_file, list<int> *set, int max_val){
    FILE *in;
    if((in = fopen(sep_file,"r")) == NULL){
        fatal_error("Can't open %s to read separator\n",sep_file);
    }

    int x,i;

    i = 0;
    while(!feof(in)){
        fscanf(in,"%d ",&x);

        if(x > max_val){
            fatal_error("Separator file %s contained label %d, which is larger than %d, the number of vertices in the current graph.\n",
                        sep_file, x, max_val);
        }
        if(x < 1){
            fatal_error("Separator file %s contained label %d, which is less than 1, the minimum vertex label in G.",
                        sep_file, x, max_val);
        }
        // Assuming 1-based to 0-based - no messing with labels here for now...
        set->push_back(x - 1);
        i++;
    }
    fclose(in);

    // return the # of elements read in
    return i;
} // read_separator_file

int read_ordering_file(const char *ordering_file, vector<int> *ordering){
    FILE *in;
    if((in = fopen(ordering_file,"r")) == NULL){
        fatal_error("Can't open %s to read ordering\n",ordering_file);
    }

    int x,i;

    i = 0;
    while(!feof(in)){
        fscanf(in,"%d ",&x);
        // Assuming 1-based to 0-based - no messing with labels here for now...
        ordering->at(i) = x - 1;
        i++;
    }
    fclose(in);

    // return the # of elements read in
    return i;
} // read_ordering_file

int read_SCOTCH_ordering_file(const char *ordering_file, vector<int> *ordering){
    FILE *in;
    if((in = fopen(ordering_file,"r")) == NULL){
        fatal_error("Can't open %s to read ordering\n",ordering_file);
    }

    int x,n,junk, i;

    fscanf(in,"%d\n",&n);
    vector<int> inv_order(n);

    i = 0;
    while(!feof(in)){
        fscanf(in,"%d %d\n",&junk, &x);
        // Assuming 1-based to 0-based - no messing with labels here for now...
        inv_order[i] = x - 1;
        i++;
    }
    fclose(in);
    for(i = 0; i < n; i++){
        ordering->at(inv_order[i]) = i;
    }

    // return the # of elements read in
    return i;
} // read_SCOTCH_ordering_file

/**
 * Forces the labels in a DIMACS file to be between 1 and n where
 * n is the # of nodes on the first line of the input DIMACS_file.  Does
 * the relabeling and writes to new_file in DIMACS format.
 */
void normalize_DIMACS_file(const char *DIMACS_file, const char *new_file){
    char line[100], format[100], x;
    int i,j,m,n,retval, max_label = -1;
    FILE *in;
    bool has_zero = false;

    if( (in = fopen(DIMACS_file,"r")) == NULL){
        fatal_error("%s:  Error opening %s for reading\n",__FUNCTION__,DIMACS_file);
    }

    n = 0;
    m = 0;
    while(!feof(in)){
        retval = fscanf(in,"%2c",line);
        if(feof(in)){
            break;
        }

        if(retval != 1){
            fatal_error("%s:  fscanf read %d char's (expected to read 1!)\n",__FUNCTION__,retval);
        }

        int i,j;
        switch(line[0])
        {
        case 'p':
            // This is the "problem line" - get n and m here
            retval = fscanf(in,"%s %d %d",format, &n, &m);

            // Simple error checking
            if((n <= 0) || (m <= 0) || (retval != 3) ){
                fatal_error("%s:  DIMACS read error - didn't understand problem line in file %s!\n"
                            "p %s %d %d\n",__FUNCTION__,DIMACS_file,format,n,m);
            }

            if((strncmp(format,"e",1) != 0) && (strncmp(format,"E",1) != 0) &&
               (strncmp(format,"edge",4) != 0) && (strncmp(format,"EDGE",4) != 0) ){
                fatal_error("%s:  DIMACS read error - problem line - FORMAT must be one of\n"
                            " e, E, edge, or EDGE\n",__FUNCTION__);
            }
            break;

        case 'c':
            // Comment line - skip and move on
            break;

        case 'n':
            // Node line - skip
            break;

        case 'e':
            // Edge line - of the form e start end
            retval = fscanf(in,"%d %d",&i,&j);
            // Simple error checking - make sure we know n and m already!
            if((n == 0) || (m == 0) ){
                fatal_error("%s:  DIMACS read error - edge line found before problem line!\n",__FUNCTION__);
            }

            if(retval != 2){
                fatal_error("%s:  DIMACS read error - didn't understand edge line\n",__FUNCTION__);
            }

            // The line appears to be valid - determine if i or j is a new max_label
            if(i > max_label){
                max_label = i;
            }
            if(j > max_label){
                max_label = j;
            }
            if((i == 0) || (j == 0) ){
                has_zero = true;
            }
            break;
        case 'x':
            // skip
            break;

        default:
            // We encountered some character that we don't understand
            fatal_error("%s:  DIMACS read error - didn't understand character %c to start a line\n",
                        __FUNCTION__,line[0]);
            break;
        }        // end switch
                 // Advance to next line if there is format at the end of it
        while(!feof(in) && (x = getc(in)) != '\n'){
            ;
        }
    }    // end while

    // Now rewind the file and go through it again
    // Create the t[] array
    vector<bool> t(max_label + 1,false);
    rewind(in);

    while(!feof(in)){
        retval = fscanf(in,"%2c",line);
        if(feof(in)){
            break;
        }

        if(retval != 1){
            fatal_error("%s:  fscanf read %d char's (expected to read 1!)\n",__FUNCTION__,retval);
        }

        int i,j;
        switch(line[0])
        {
        case 'p':
            // This is the "problem line" - we already know n and m

            break;

        case 'c':
            // Comment line - skip and move on
            break;

        case 'n':
            // Node line - skip
            break;

        case 'e':
            // Edge line - of the form e start end
            retval = fscanf(in,"%d %d",&i,&j);
            // Set t[i] and t[j] to true for these 2 labels
            t[i] = true;
            t[j] = true;

            break;
        case 'x':
            // skip
            break;

        default:
            // We encountered some character that we don't understand
            fatal_error("%s:  DIMACS read error - didn't understand character %c to start a line\n",
                        __FUNCTION__,line[0]);
            break;
        }        // end switch
                 // Advance to next line if there is format at the end of it
        while(!feof(in) && (x = getc(in)) != '\n'){
            ;
        }
    }    // end while

    // Now go through t and create the perm[] array to translate the indices from
    // 1...max_label to 1...n
    // or 0, .. max label if has_zero flag is set.
    vector<int> perm(max_label + 1,-1);
    j = 1;
    if(has_zero){
        perm[0] = 1;
        j++;
    }
    for(i = 1; i <= max_label; i++){
        if(t[i]){
            perm[i] = j;
            j++;
        }
    }

    // Now go through the file one last time and write out perm[i] perm[j] for the
    // edge i, j
    rewind(in);
    FILE *out = fopen(new_file,"w");
    if(!out){
        fatal_error("%s:  Couldn't open %s for writing\n",__FUNCTION__,new_file);
    }

    fprintf(stderr,"Writing to %s\n",new_file);
    // Write a comment
    fprintf(out,"c Reformatted version of original file %s\n",DIMACS_file);
    fprintf(out,"c Original file had max_label of %d\n",max_label);
    fprintf(out,"p edge %d %d\n",n,m);
    fflush(out);

    while(!feof(in)){
        retval = fscanf(in,"%2c",line);
        if(feof(in)){
            break;
        }

        if(retval != 1){
            fatal_error("%s:  fscanf read %d char's (expected to read 1!)\n",__FUNCTION__,retval);
        }

        int i,j;
        max_label = -1;
        switch(line[0])
        {
        case 'p':
            // This is the "problem line" - we already know n and m
            break;

        case 'c':
            // Comment line - skip and move on
            break;

        case 'n':
            // Node line - skip
            break;

        case 'e':
            // Edge line - of the form e start end
            retval = fscanf(in,"%d %d",&i,&j);
            // Write out the edge perm[i] perm[j]
            fprintf(out,"e %d %d\n",perm[i],perm[j]);
            break;
        case 'x':
            // skip
            break;

        default:
            // We encountered some character that we don't understand
            fatal_error("%s:  DIMACS read error - didn't understand character %c to start a line\n",
                        __FUNCTION__,line[0]);
            break;
        }        // end switch
                 // Advance to next line if there is format at the end of it
        while(!feof(in) && (x = getc(in)) != '\n'){
            ;
        }
    }    // end while

    fclose(out);
    fclose(in);
    return;
} // normalize_DIMACS_file

/* Returns a string naming the elimination ordering routine corresponding to the given integer.
 * If input is not one of the GD_ constants defined in GraphDecomposition.h, this returns "Unknown".
 */
const char *eo_name(int eo_id){
    int index = eo_id;
    if((index > GD_NUM_ELIM_HEURISTICS) || (index < 0) ){
        return EO_NAMES[GD_NUM_ELIM_HEURISTICS];        //this is set to unknown
    }
    else {
        return EO_NAMES[index];
    }
}

/** Returns an int_int where p1 is the position of the high bit
 * (floor of log 2 of x), and p2 is the # of bits set in x.
 * If p2==1, then x is a power of two. If p2>1, then x is not a power of 2.
 */
int_int mylog2(int x){
    int_int ans;
    ans.p1 = 0;
    ans.p2 = 0;
    // The # of bits in an int is its size (in bytes)*8
    for(int i = 0; i < (int)(sizeof(int)) << 3; i++){
        if(x & (1 << i)){
            // bit i is set
            ans.p1 = i;
            ans.p2++;
        }
    }
    return ans;
}

std::string str_to_up(string s){
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

void split(const std::string& s, char sep, vector<int>& v){
    stringstream ss(s);      //< a stringstream to do our parsing for us
    int k = -1;    //< just some variable
    string temp;
    while(std::getline(ss, temp, sep)){
        stringstream convert(temp);
        convert >> k;
        if(k != -1){
            v.push_back(k);
        }
        k = -1;
    }
}     // split

/**
 * Write a degree distribution out to file.
 * \param[in] filename filename to write output to
 * \param[in] dist a vector<int>, indexed on degree
 */
void write_degree_distribution(string filename, const vector<int> &dist){
    ofstream outfile;

    outfile.open(filename.c_str());

    if(!outfile.is_open()){
        cerr << "Error opening " << filename << "for writing\n";
    }
    else {
        int i;
        for(i = 0; i < dist.size(); i++){
            outfile << i << " " <<  dist[i] << "\n";
        }
    }

    outfile.close();
} // write_degree_distribution

/**
 * Write an eccentricity distribution out to file.
 * \param[in] filename filename to write output to
 * \param[in] dist a vector<int>, indexed on degree
 */
void write_eccentricity_distribution(string filename, const vector<double> &dist){
    ofstream outfile;

    outfile.open(filename.c_str());

    if(!outfile.is_open()){
        cerr << "Error opening " << filename << "for writing\n";
    }
    else {
        int i;
        for(i = 0; i < dist.size(); i++){
            outfile << i << " " <<  dist[i] << "\n";
        }
    }

    outfile.close();
} // write_eccentricity_distribution

/**
 * Write a expansion out to file.
 * \param[in] filename filename to write output to
 * \param[in] dist a vector<int>, indexed on degree
 */
void write_expansion(string filename, const vector<double> &expansion){
    ofstream outfile;

    outfile.open(filename.c_str());

    if(!outfile.is_open()){
        cerr << "Error opening " << filename << "for writing\n";
    }
    else {
        int i;
        for(i = 0; i < expansion.size(); i++){
            outfile << i << " " <<  expansion[i] << "\n";
        }
    }

    outfile.close();
} // write_expansion

/**
 * Write a kcore list out to file.
 * \param[in] filename filename to write output to
 * \param[in] kcores a vector<int>, indexed on vertex number
 */
void write_kcores(string filename, const vector<int> &kcores){
    ofstream outfile;

    outfile.open(filename.c_str());

    if(!outfile.is_open()){
        cerr << "Error opening " << filename << "for writing\n";
    }
    else {
        int i;
        for(i = 0; i < kcores.size(); i++){
            outfile << i << " " << kcores[i] << "\n";
        }
    }

    outfile.close();
} // write_k_cores

/**
 * Write an eccentricity list out to file.
 * \param[in] filename filename to write output to
 * \param[in] ecc a vector<int>, indexed on vertex number
 */
void write_eccentricity(string filename, const vector<int> &ecc){
    ofstream outfile;

    outfile.open(filename.c_str());

    if(!outfile.is_open()){
        cerr << "Error opening " << filename << "for writing\n";
    }
    else {
        int i;
        for(i = 0; i < ecc.size(); i++){
            outfile << i << " " << ecc[i] << "\n";
        }
    }

    outfile.close();
} // write_k_cores

/**
 * Write a betweenness centrality vecotr out to file.
 *
 * \param[in] filename filename to write output to
 * \param[in] bc a vector<uint64_t>, indexed on vertex number
 */
void write_betweenness(string filename, const vector<double> &bc){
    ofstream outfile;

    outfile.open(filename.c_str());

    if(!outfile.is_open()){
        cerr << "Error opening " << filename << "for writing\n";
    }
    else {
        outfile.precision(10);
        int i;
        for(i = 0; i < bc.size(); i++){
            outfile << fixed << i << " " << bc[i] << "\n";
        }
    }

    outfile.close();
} // write_betweenness

/**
 * Write a delta hyperbolicity vector out to file.
 *
 * \param[in] filename filename to write output to
 * \param[in] bc a vector<uint64_t>, indexed on vertex number
 */
void write_delta_hyperbolicity(string filename, const vector< vector<double> > &delta){
    ofstream outfile;

    outfile.open(filename.c_str());

    if(!outfile.is_open()){
        cerr << "Error opening " << filename << "for writing\n";
    }
    else {
        outfile.precision(10);
        for(int idx = 0; idx < delta.size(); idx++){
            for(int jdx = 0; jdx < delta[idx].size(); jdx++){
                outfile << delta[idx][jdx] << " ";
            }
            outfile << endl;
        }
    }

    outfile.close();
} // write_betweenness

/**
 * Write an all pairs shortest path matrix  out to file.
 *
 * \param[in] filename filename to write output to
 * \param[in] apsp a vector< vector<int> >, indexed on vertex number, assumed to be |V|x|V|
 */
void write_apsp_matrix(string filename, vector< vector<int> > &apsp){
    int i;
    int n = apsp.size();  //assuming it's sized correctly
    cerr << "Writing apsp matrix of size: " << n << endl;

    FILE *outfile;
    outfile = fopen(filename.c_str(), "w+");

    if(!outfile){
        cerr << "Error opening " << filename << "for writing\n";
    }

    fwrite((void *) &n, 1, sizeof(int), outfile);

    //#pragma omp parallel for default(none) share(n, apsp)
    int ret;
    for(i = 0; i < n; i++){
        int *arr = &apsp[i].front();
        fwrite(arr, n, sizeof(int), outfile);
    }

    fclose(outfile);
} // write_apsp_matrix

/**
 * Read an all pairs shortest path matrixfrom file.
 *
 * \param[in] filename filename to read
 * \param[out] apsp a vector< vector<int> >, indexed on vertex number, assumed to be |V|x|V|
 */
void read_apsp_matrix(string filename, vector< vector<int> > &apsp){
    int i;
    int n;

    FILE *infile;
    infile = fopen(filename.c_str(), "r");

    if(!infile){
        cerr << "Error opening " << filename << "for writing\n";
    }

    fread((void *) &n, 1, sizeof(int), infile);
    cerr << "Reading apsp matrix of size: " << n << endl;

    apsp.resize(n);

    //#pragma omp parallel for default(none) share(n, apsp)
    int ret;
    for(i = 0; i < n; i++){
        apsp[i].resize(n);
        int *arr = &apsp[i].front();
        fread(arr, n, sizeof(int), infile);
    }

    fclose(infile);
} // write_apsp_matrix

