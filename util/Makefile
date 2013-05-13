include ../make.inc

INC_DIR = -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR)

LIBS = -lm $(INDDGO_LIB) -ltreed -lgraphd $(METIS_LIB) $(SSPARSE_LIB) $(ARPACK_LIB) $(GMP_LIB) 

PKT_SRCS =  ./src/gen_pkt.cpp 
PKT_OBJ =  $(PKT_SRCS:.cpp=.o)
PKT_EXE = $(INDDGO_BIN)/gen_pkt

STATS_SRCS =  ./src/td_stats.cpp 
STATS_OBJ =  $(STATS_SRCS:.cpp=.o)
STATS_EXE = $(INDDGO_BIN)/td_stats

E2D_SRCS = ./src/edge2dimacs.cpp
E2D_OBJ = $(E2D_SRCS:.cpp=.o)
E2D_EXE = $(INDDGO_BIN)/edge2dimacs

M2D_SRCS = ./src/metis2dimacs.cpp
M2D_OBJ = $(M2D_SRCS:.cpp=.o)
M2D_EXE = $(INDDGO_BIN)/metis2dimacs

all: $(PKT_EXE) $(STATS_EXE) $(E2D_EXE) $(M2D_EXE)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(PKT_EXE): $(PKT_OBJ) 
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(STATS_EXE): $(STATS_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(E2D_EXE): $(E2D_OBJ)  ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(M2D_EXE): $(M2D_OBJ)  ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

clean: 
	-rm -rf $(PKT_OBJ) $(PKT_EXE) $(STATS_OBJ) $(STATS_EXE) $(E2D_OBJ) $(E2D_EXE) $(M2D_EXE) $(M2D_OBJ)
