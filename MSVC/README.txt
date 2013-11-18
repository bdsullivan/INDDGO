All the projects in this directory should build without issue
except for the CPX_solve_wis project as it requires both
GLPK (open source) and CPLEX (not open source).  If you have
both GLPK and CPLEX, then it should be straightforward to
add the appropriate paths and libs to the project properties
by just editing what is in there now as that is specific
to a particular machine with the paths, etc.

Here is a compile line for CPX_solve_wis :
 g++ CPX_solve_wis.cpp -O2 -g  -I../../lib_graphd/inc
-I/opt/ibm/ILOG/CPLEX_Studio125/cplex/include/
-L/opt/ibm/ILOG/CPLEX_Studio125/cplex/lib/x86-64_sles10_4.1/static_pic/
-L../../lib/ -lgraphd -lpthread -lcplex -lglpk
