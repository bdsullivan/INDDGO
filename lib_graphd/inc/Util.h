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

#ifndef UTIL_H_
#define UTIL_H_

#define GD_STAT_MEAN 1 //arithmetic mean
#define GD_STAT_STD 2 //standard deviation
#define GD_STAT_MED 3 // median
#define GD_STAT_COUNT 4 //number of vertices with non-negative score

class int_double
{
friend std::ostream &operator<<(std::ostream &, const int_double &);
public:
int i;
double d;
};

/*
 * An int_int class useful for sorting a list and
 * keeping track of the resulting permutation.
 */
class int_int
{
friend std::ostream &operator<<(std::ostream &, const int_int &);
public:
int p1;
int p2;

friend bool operator<(const int_int& a, const int_int &b);
friend bool operator==(const int_int& a, const int_int &b);
};

void get_DIMACS_dimensions(const char *DIMACS_file, int *n, int *m);
class VertexLabel
{
friend std::ostream &operator<<(std::ostream &, const VertexLabel &);

public:
int id;
// The entries list is kept in decreasing order
std::list<int> entries;

// Use this to determine lexicographic ordering
friend bool operator<(const VertexLabel& a, const VertexLabel &b);
};

bool decreasing(int a, int b);
bool int_d_sort(int_double i, int_double j);
//checks if the int in i is < the int in j.
bool int_d_isort(int_double i, int_double j);
void score_sort(std::list<int> *mylist, std::vector<int> *scores);

bool empty_intersection(std::list<int> *S1, std::list<int> *S2);

int read_ordering_file(const char *ordering_file, std::vector<int> *ordering);
int read_SCOTCH_ordering_file(const char *ordering_file, std::vector<int> *ordering);
int read_separator_file(const char *sep_file, std::list<int> *set, int max_val);
int parseLine(char *line);
int getHWmem();

//Aaron's Functions: Read in a color vector from file, returns whether range was provided
//produce a linear scaled RGB value.
bool read_color_file(const char input_file[], double & max_color, double & min_color, vector<double> & color_vector);
void get_rgb_value(char rgb[], double val);
double get_statistics(const vector<int> & nodes, const vector<double> & scores, const int FLAG);

void normalize_DIMACS_file(const char *DIMACS_file,const char *new_file);

const char *eo_name(int eo_id);

// String manipulation functions
std::string str_to_up(std::string s);
void split(const std::string& s, char sep, vector<int>& v);      // used in metisgraph reader

/**
 * \brief write a degree distribution out to a file
 */
void write_degree_distribution(string filename, const vector<int> &dist);

/**
 * \brief write an eccentricity distribution out to a file
 */
void write_eccentricity_distribution(string filename, const vector<double> &dist);

/**
 * \brief write a k_core list out to a file
 */
void write_kcores(string filename, const vector<int> &kcores);

/**
 * \brief write an eccentricity list out to a file
 */
void write_eccentricity(string filename, const vector<int> &ecc);

/**
 * \brief write a betweenness centrality list out to a file
 */
void write_betweenness(string filename, const vector<double> &bc);

/**
 * \brief write a delta hyperbolicity distribution out to a file
 */
void write_delta_hyperbolicity(string filename, const vector< vector<double> > &delta);

/**
 * \brief write an expansion list out to a file
 */
void write_expansion(string filename, const vector<double> &expansion);

/**
 * \brief Write an all pairs shortest path matrix  out to file.
 */
void write_apsp_matrix(string filename, vector< vector<int> > &apsp);

/**
 * \brief Reads an all pairs shortest path matrix from file
 */
void read_apsp_matrix(string filename, vector< vector<int> > &apsp);

int_int mylog2(int x);

#endif /* UTIL_H_ */
