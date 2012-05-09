CC=g++
CCC=gcc
CFLAGS= -Wall -O3 -DHAS_GOBLIN -g 
VALGRIND_FLAGS= -O0 -g
PFLAGS= -O3 -pg

# Set this to 0 if you do not have METIS
HAS_METIS =1
ifeq ($(HAS_METIS),1)

METIS_TEST_SRC = ./src/metis_test.cpp
METIS_TEST_EXE = ./bin/metis_test
VAL_METIS_TEST_EXE = ./bin/val_metis_test
PROF_METIS_TEST_EXE = ./bin/prof_metis_test
METIS_LIBDIR= -L../../metis-4.0/ 
METIS_INCDIR= ../../metis-4.0/Lib 
METIS_LIB= -lmetis
else

METIS_TEST_SRC =
METIS_TEST_EXE =
VAL_METIS_TEST_EXE = 
PROF_METIS_TEST_EXE = 
METIS_LIBDIR=
METIS_INCDIR=
METIS_LIB=

endif

HAS_AMD=1
ifeq ($(HAS_AMD),1)

AMD_INCDIR= -I../../AMD/Include -I../../UFconfig
AMD_LIBDIR= -L../../AMD/Lib
AMD_LIB= -lamd

else

AMD_INCDIR=
AMD_LIBDIR=
AMD_LIB=

endif

ifeq ($(HAS_METIS),1)
CFLAGS := $(CFLAGS) -DHAS_METIS
endif
ifeq ($(HAS_AMD),1)
CFLAGS := $(CFLAGS) -DHAS_AMD
endif

INC_DIR = -I../graph_decomposition/inc -I./inc -I$(GOBLIN_DIR)/include -I$(METIS_INCDIR) $(AMD_INCDIR)
GOBLIN_DIR = ../../goblin.2.8b30
BDTEST_EXE = bin/bdtest
BDTEST_SRC = src/bdtest.cpp
GEN_BD_SPLIT_EXE = bin/gen_bd_split
GEN_BD_SPLIT_SRC = src/gen_bd_split.cpp
GEN_BD_SPAWN_EXE = bin/gen_bd_spawn
GEN_BD_SPAWN_SRC = src/gen_bd_spawn.cpp


DATA_DIR=../sample_graphs/

VALGRIND_COMMAND = valgrind --leak-check=yes --show-reachable=yes -v

DOCS_DIR = ./doc/*
DOXYFILE = ./Doxyfile
ARPACK_LIBDIR= ../../ARPACK

HAS_GOBLIN= 1

LIBS = -lm -lgd $(METIS_LIB) $(AMD_LIB)
LIB_DIR = -L../graph_decomposition/lib/ -L./lib $(METIS_LIBDIR) $(AMD_LIBDIR)

SRCS=./src/BDTree.cpp ./src/BDTreeEdge.cpp ./src/BDTreeNode.cpp #./src/ArpackUtils.cpp  
OBJS=$(SRCS:.cpp=.o) 

all: #$(OBJS)
	$(error	This code is deprecated, and incompatible with the current lib_graphd. It was written using a previous version of the graph decomposition library. Most functionality should be re-attainable, but new GraphEOUtil classes, etc need to be incorporated into the code) 
	gcc ./src/pathfinder.c ./src/paths.c -I./inc -o bin/pathfinder
	$(CC) $(CFLAGS) $(INC_DIR) $(BDTEST_SRC) $(OBJS) -L$(ARPACK_LIBDIR) -L$(GOBLIN_DIR) $(LIB_DIR) $(LIBS) -lgoblin -larpack_x86_64 -lgfortran   -o $(BDTEST_EXE)
	$(CC) $(CFLAGS) $(INC_DIR) $(GEN_BD_SPLIT_SRC) $(OBJS) -L$(ARPACK_LIBDIR) -L$(GOBLIN_DIR) $(LIB_DIR) $(LIBS) -lgoblin -larpack_x86_64 -lgfortran  -o $(GEN_BD_SPLIT_EXE)
	$(CC) $(CFLAGS) $(INC_DIR) $(GEN_BD_SPAWN_SRC) $(OBJS) -L$(ARPACK_LIBDIR) -L$(GOBLIN_DIR) $(LIB_DIR) $(LIBS) -lgoblin -larpack_x86_64 -lgfortran  -o $(GEN_BD_SPAWN_EXE)


.cpp.o  :
	$(CC) -c $(CFLAGS) $(INC_DIR) $< -o $@
.c.o	:	
	$(CCC) -c $(CFLAGS) $(INC_DIR) $< -o $@

docs: $(DOXYFILE)
	doxygen $(DOXYFILE)

clean: 
	-rm -rf $(OBJS)
	-rm -rf $(BDTEST_EXE)
	-rm -rf $(GEN_BD_SPLIT_EXE)
	-rm -rf $(GEN_BD_SPAWN_EXE)
	-rm -rf $(DOCS_DIR)

