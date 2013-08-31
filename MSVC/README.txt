All the projects in this directory should build without issue
except for the CPX_solve_wis project as it requires both
GLPK (open source) and CPLEX (not open source).  If you have
both GLPK and CPLEX, then it should be straightforward to
add the appropriate paths and libs to the project properties
by just editing what is in there now as that is specific
to a particular machine with the paths, etc.