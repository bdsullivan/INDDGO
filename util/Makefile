include ../make.inc

INC_DIR = -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR)

LIBS = -lm $(INDDGO_LIB) -lgraphd -ltreed $(METIS_LIB) $(SSPARSE_LIB) $(ARPACK_LIB) $(GMP_LIB) 

PKT_SRCS =  ./src/gen_pkt.cpp 
PKT_OBJ =  $(PKT_SRCS:.cpp=.o)
PKT_EXE = $(INDDGO_BIN)/gen_pkt

all: $(PKT_EXE)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(PKT_EXE): $(PKT_OBJ) 
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

clean: 
	-rm -rf $(PKT_OBJ) $(PKT_EXE)
