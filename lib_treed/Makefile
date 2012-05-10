include ../make.inc

DOCS_DIR = ./doc/*
DOXYFILE = ./Doxyfile

INC_DIR= -I./inc -I$(GRAPH)/inc  $(METIS_INCLUDES) $(AMD_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR)

LIBS = -lm -lgraphd $(METIS_LIB) $(AMD_LIB) $(ARPACK_LIB) $(MPILIBS) $(MADNESS_LIB) $(GMP_LIB) $(PARMETIS_LIB)

LIB_DIR = -L./lib -L$(GRAPH)/lib/ -L$(METIS_LIB_DIR) $(AMD_LIBDIR) $(ARPACK_LIBDIR) -L$(GRAPH)/lib $(MADNESS_LIBDIR) $(GMP_LIBDIR) -L$(PARMETIS_LIBDIR) $(INDDGO_LIB)


TD_LIB_NAME = libtreed.a

#optional local libraries for profiling and valgrind
TD_PROF_LIB = ./lib/libproftd.a
TD_VALGRIND_LIB = ./lib/libvaltd.a

SRCS=  ./src/TDTree.cpp ./src/TDTreeNode.cpp ./src/TDSolution.cpp ./src/TDDynamicProgramming.cpp ./src/TDMadTreeNode.cpp
OBJS=$(SRCS:.cpp=.o)

all: $(INDDGO_LIB)/$(TD_LIB_NAME)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

# The "normal" TD library
$(INDDGO_LIB)/$(TD_LIB_NAME): $(OBJS) $(PARMETISOBJ)
	$(AR) $(ARFLAGS) $@ $(OBJS) $(PARMETISOBJ)
	ranlib $@

locallib: $(INDDGO_LIB)/$(TD_LIB_NAME)
	cp $(INDDGO_LIB)/$(TD_LIB_NAME) ./lib/$(TD_LIB_NAME)

# The profiled library
td_profiled: $(OBJS) $(PARMETISOBJ)
	$(CXX) -c $(PFLAGS) $(INC_DIR) $(LIB_DIR) $(SRCS)
	$(AR) $(ARFLAGS) $(TD_PROF_LIB) $(OBJS) $(PARMETISOBJ)
	ranlib $(TD_PROF_LIB)

# The valgrind-prepared library
td_valgrind: $(OBJS) $(PARMETISOBJ)
	$(CXX) -c $(VALGRIND_FLAGS) $(INC_DIR) $(LIB_DIR) $(SRCS)
	$(AR) $(ARFLAGS) $(TD_VALGRIND_LIB) $(OBJS) $(PARMETISOBJ)
	ranlib $(TD_VALGRIND_LIB)

clean:
	-rm -rf $(OBJS)
	-rm -rf $(INDDGO_LIB)/$(TD_LIB_NAME) ./lib/$(TD_LIB_NAME)
	-rm -rf $(TD_VALGRIND_LIB)
	-rm -rf $(TD_PROF_LIB)
	-rm -rf $(DOCS_DIR)

docs: $(DOXYFILE)
	doxygen $(DOXYFILE)



