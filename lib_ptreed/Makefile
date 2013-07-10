include ../make.inc

DOX_DIR = ./doc/
DOXYFILE = ./Doxyfile

INC_DIR= -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR) $(SLEPC_INCDIR) $(PETSC_INCDIR)

LIBS = -lm $(INDDGO_LIB) -ltreed -lgraphd $(METIS_LIB) $(SSPARSE_LIB) $(ARPACK_LIB) $(MPILIBS) $(MADNESS_LIB) $(GMP_LIB) $(PARMETIS_LIB) $(SLEPC_LIB) $(PETSC_LIB)

PTD_LIB_NAME = libptreed.a

SRCS = ./src/PTD.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(INDDGO_LIB_DIR)/$(PTD_LIB_NAME)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(INDDGO_LIB_DIR)/$(PTD_LIB_NAME): $(OBJS) 
	$(AR) $(ARFLAGS) $@ $(OBJS) 
	ranlib $@
clean: 
	-rm -f $(INDDGO_LIB_DIR)/$(PTD_LIB_NAME) $(OBJS)

clean_doc:
	-rm -rf $(DOX_DIR)

doc: $(DOXYFILE)
	doxygen $(DOXYFILE)
