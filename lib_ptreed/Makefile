include ../make.inc

INC_DIR= -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(AMD_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR)

LIBS = -lm -ltreed -lgraphd $(METIS_LIB) $(AMD_LIB) $(ARPACK_LIB) $(MPILIBS) $(MADNESS_LIB) $(GMP_LIB) $(PARMETIS_LIB)

LIB_DIR = -L$(GRAPH)/lib/ -L$(METIS_LIB_DIR) $(AMD_LIBDIR) $(ARPACK_LIBDIR) -L$(GRAPH)/lib $(MADNESS_LIBDIR) $(GMP_LIBDIR) -L$(PARMETIS_LIBDIR) $(INDDGO_LIB)

PTD_LIB_NAME = libptreed.a

SRCS = ./src/PTD.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(INDDGO_LIB)/$(PTD_LIB_NAME)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(INDDGO_LIB)/$(PTD_LIB_NAME): $(OBJS) 
	$(AR) $(ARFLAGS) $@ $(OBJS) 
	ranlib $@
clean: 
	-rm -f $(INDDGO_LIB)/$(PTD_LIB_NAME) $(OBJS)
