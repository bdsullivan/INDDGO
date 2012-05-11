include ../make.inc

# Optimal WIS for this one is 79 
TEST_GRAPH = $(SRC_DIR)/sample_graphs/pkt.200.70.75_0.dimacs

INC_DIR=  -I./inc -I$(GRAPH)/inc -I$(TREE)/inc/  $(METIS_INCLUDES) $(AMD_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR) -I$(PTREE)/inc/

LIBS = -lm $(INDDGO_LIB) -lgraphd -ltreed $(AMD_LIB) $(ARPACK_LIB) $(MPILIBS) $(MADNESS_LIB) $(GMP_LIB) $(METIS_LIB)

WEIGHTED_IND_SET_SRC = ./src/weighted_ind_set.cpp ./src/main.cpp
WEIGHTED_IND_SET_EXE = $(INDDGO_BIN)/serial_wis

MAD_WIS_SRC = ./src/MadnessTableProcessor.cpp ./src/mad_main.cpp ./src/weighted_ind_set.cpp 
MAD_WIS_EXE = $(INDDGO_BIN)/parallel_wis

#other executables compiled into local bin directory; not for prime-time use
PROF_WEIGHTED_IND_SET_EXE = ./bin/prof_weighted_ind_set

PARALLEL_WIS_SRC = ./src/parallel_wis.cpp ./src/parallel_main.cpp ./src/weighted_ind_set.cpp
PARALLEL_WIS_EXE = ./bin/parallel_wis

PARMETIS_EO_SRC = ./src/ORD.cpp ./src/parmetis_main.cpp ./src/weighted_ind_set.cpp
PARMETIS_EO_EXE = ./bin/parmetiseo

PGAVRIL_SRC =  ./src/pgavril_main.cpp ./src/weighted_ind_set.cpp
PGAVRIL_EXE = ./bin/pgavril

SHAREDSRC = ./src/DPinfo.cpp 
SHAREDOBJ = ./src/DPinfo.o $(INDDGO_LIB_DIR)/libgraphd.a $(PARMETISOBJ)

SRCS = $(MAD_WIS_SRC) $(PARALLEL_WIS_SRC) $(PROF_WEIGHTED_IND_SET_SRC) $(WEIGHTED_IND_SET_SRC) $(SHAREDSRC)
OBJS=$(SRCS:.cpp=.o)
EXES =  $(WEIGHTED_IND_SET_EXE) $(PROF_WEIGHTED_IND_SET_EXE) $(PARALLEL_WIS_EXE) $(PARMETIS_EO_EXE) $(MAD_WIS_EXE)


ALL= $(WEIGHTED_IND_SET_EXE)

#Optional partial parallel binaries
#ifeq ($(IS_PARALLEL),1) 
#ALL:= $(ALL) pgavril parallel_wis
#endif

ifeq ($(HAS_MADNESS),1) 
ALL:= $(ALL) $(MAD_WIS_EXE)
endif

#Optional parallel elimination ordering binary
#ifeq ($(HAS_PARMETIS),1) 
#ALL:= $(ALL) parmetiseo
#endif


all: $(ALL)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR) $< -o $@

# Binaries

$(WEIGHTED_IND_SET_EXE): $(SHAREDOBJ) 
	$(CXX) $(CFLAGS) $(WEIGHTED_IND_SET_SRC) $(SHAREDOBJ) $(INC_DIR)  $(LIBS) -o $(WEIGHTED_IND_SET_EXE)

$(MAD_WIS_EXE): $(SHAREDOBJ) $(MAD_WIS_SRC)
	$(CXX) $(CFLAGS) $(MAD_WIS_SRC) $(SHAREDOBJ) $(INC_DIR) $(LIBS) -lptreed -o $(MAD_WIS_EXE)

parmetiseo:
	$(CXX) $(CFLAGS) $(PARMETIS_EO_SRC) $(SHAREDOBJ) $(INC_DIR)  $(LIBS) -o $(PARMETIS_EO_EXE)

pgavril: $(SHAREDOBJ)
	$(CXX) $(CFLAGS) $(PGAVRIL_SRC) $(SHAREDOBJ) $(INC_DIR)  $(LIBS) -lptreed -o $(PGAVRIL_EXE)

parallel_wis: $(SHAREDOBJ)
	$(CXX) $(CFLAGS) $(PARALLEL_WIS_SRC) $(SHAREDOBJ) $(INC_DIR) $(LIBS) -o $(PARALLEL_WIS_EXE)
clean:
	-rm -rf $(OBJS)
	-rm -rf $(EXES)

# make test will run a simple test using basic binary
test: weighted_ind_set
	$(WEIGHTED_IND_SET_EXE) -gavril  -mind -f $(TEST_GRAPH)  -v

# make check will run WIS on TEST_GRAPH using all kinds of options and will make
# sure we get the same answer
check: weighted_ind_set
	$(WEIGHTED_IND_SET_EXE) -gavril  -mind -f $(TEST_GRAPH) > gav_mind.out
	$(WEIGHTED_IND_SET_EXE) -gavril  -mind -f $(TEST_GRAPH) -make_nice > gav_mind_make_nice.out
	$(WEIGHTED_IND_SET_EXE) -gavril  -mind -f $(TEST_GRAPH) -pc > gav_mind_pc.out
	$(WEIGHTED_IND_SET_EXE) -gavril  -mind -f $(TEST_GRAPH) -make_nice -pc > gav_mind_make_nice_pc.out
	$(WEIGHTED_IND_SET_EXE) -gavril -mind -f $(TEST_GRAPH) -nice -pc > nice_mind.out
	$(WEIGHTED_IND_SET_EXE) -gavril  -mind -f $(TEST_GRAPH) -make_nice -pc > gav_mind_make_nice_pc.out
	$(WEIGHTED_IND_SET_EXE) -gavril -mind -f $(TEST_GRAPH) -nice -pc -nonniceDP > nice_mind_nonniceDP.out

