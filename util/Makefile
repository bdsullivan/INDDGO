include ../make.inc

INC_DIR = -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(AMD_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR)

LIBS = -lm -lgraphd -ltreed $(METIS_LIB) $(AMD_LIB) $(ARPACK_LIB) $(GMP_LIB) 

LIB_DIR = -L$(INDDGO_LIB) -L$(METIS_LIB_DIR) $(AMD_LIBDIR) $(ARPACK_LIBDIR) $(GMP_LIBDIR) 

PKT_SRCS =  ./src/gen_pkt.cpp 
PKT_OBJ =  $(PKT_SRCS:.cpp=.o)
PKT_EXE = $(INDDGO_BIN)/gen_pkt

all: $(PKT_EXE)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(PKT_EXE): $(PKT_OBJ) 
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIB_DIR) $(LIBS) -o $@ 

clean: 
	-rm -rf $(PKT_OBJ) $(PKT_EXE)
