include ../make.inc

#DOCS_DIR = ./doc/*
DOX_DIR = ./doc/
DOXYFILE = ./Doxyfile

INC_DIR= -I./inc -I$(GRAPH)/inc  $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR)

LIBS = -L./lib -lm -L$(GRAPH)/lib -lgraphd $(INDDGO_LIB) $(METIS_LIB) $(SSPARSE_LIB) $(ARPACK_LIB) $(MPILIBS) $(MADNESS_LIB) $(GMP_LIB) $(PARMETIS_LIB)

TD_LIB_NAME = libtreed.a

#optional local libraries for profiling and valgrind
TD_PROF_LIB = ./lib/libproftd.a
TD_VALGRIND_LIB = ./lib/libvaltd.a

SRCS=  ./src/TDTree.cpp ./src/TDTreeNode.cpp ./src/TDSolution.cpp ./src/TDDynamicProgramming.cpp ./src/TDMadTreeNode.cpp ./src/TDUtil.cpp
OBJS=$(SRCS:.cpp=.o)

ifeq ($(HAS_SUITESPARSE),1)
OBJS = $(SRCS:.cpp=.o) ./src/SuperWrap.o
endif

all: $(INDDGO_LIB_DIR)/$(TD_LIB_NAME)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

.c.o:
	$(CC) $(CFLAGS) -c $(INC_DIR)  $< -o $@

# The "normal" TD library
$(INDDGO_LIB_DIR)/$(TD_LIB_NAME): $(OBJS) $(PARMETISOBJ)
	$(AR) $(ARFLAGS) $@ $(OBJS) $(PARMETISOBJ)
	ranlib $@

locallib: $(INDDGO_LIB_DIR)/$(TD_LIB_NAME)
	cp $(INDDGO_LIB_DIR)/$(TD_LIB_NAME) ./lib/$(TD_LIB_NAME)

# The profiled library
td_profiled: $(OBJS) $(PARMETISOBJ)
	$(CXX) -c $(PFLAGS) $(INC_DIR) $(LIBS) $(SRCS)
	$(AR) $(ARFLAGS) $(TD_PROF_LIB) $(OBJS) $(PARMETISOBJ)
	ranlib $(TD_PROF_LIB)

# The valgrind-prepared library
td_valgrind: $(OBJS) $(PARMETISOBJ)
	$(CXX) -c $(VALGRIND_FLAGS) $(INC_DIR) $(LIBS) $(SRCS)
	$(AR) $(ARFLAGS) $(TD_VALGRIND_LIB) $(OBJS) $(PARMETISOBJ)
	ranlib $(TD_VALGRIND_LIB)

clean:
	-rm -rf $(OBJS)
	-rm -rf $(INDDGO_LIB_DIR)/$(TD_LIB_NAME) ./lib/$(TD_LIB_NAME)
	-rm -rf $(TD_VALGRIND_LIB)
	-rm -rf $(TD_PROF_LIB)

clean_doc:
	-rm -rf $(DOX_DIR)

doc: $(DOXYFILE)
	doxygen $(DOXYFILE)



