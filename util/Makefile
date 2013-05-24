include ../make.inc
DOX_DIR = ./doc/
DOXYFILE = ./Doxyfile

INC_DIR = -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR)

LIBS = -lm $(INDDGO_LIB) -ltreed -lgraphd $(METIS_LIB) $(SSPARSE_LIB) $(ARPACK_LIB) $(GMP_LIB) 

PKT_SRCS =  ./src/gen_pkt.cpp 
PKT_OBJ =  $(PKT_SRCS:.cpp=.o)
PKT_EXE = $(INDDGO_BIN)/gen_pkt

TD_STATS_SRCS =  ./src/td_stats.cpp 
TD_STATS_OBJ =  $(TD_STATS_SRCS:.cpp=.o)
TD_STATS_EXE = $(INDDGO_BIN)/td_stats

G_STATS_SRCS =  ./src/graph_stats.cpp 
G_STATS_OBJ =  $(G_STATS_SRCS:.cpp=.o)
G_STATS_EXE = $(INDDGO_BIN)/graph_stats

E2D_SRCS = ./src/edge2dimacs.cpp
E2D_OBJ = $(E2D_SRCS:.cpp=.o)
E2D_EXE = $(INDDGO_BIN)/edge2dimacs

E2A_SRCS = ./src/edge2adjlist.cpp
E2A_OBJ = $(E2A_SRCS:.cpp=.o)
E2A_EXE = $(INDDGO_BIN)/edge2adjlist

M2D_SRCS = ./src/metis2dimacs.cpp
M2D_OBJ = $(M2D_SRCS:.cpp=.o)
M2D_EXE = $(INDDGO_BIN)/metis2dimacs

TRI_SRCS =  ./src/count-tri.cpp 
TRI_OBJ =  $(TRI_SRCS:.cpp=.o)
TRI_EXE = $(INDDGO_BIN)/count-tri

SRT_SRCS =  ./src/shortest-paths.cpp 
SRT_OBJ =  $(SRT_SRCS:.cpp=.o)
SRT_EXE = $(INDDGO_BIN)/shortest-paths

all: $(PKT_EXE) $(TD_STATS_EXE) $(E2D_EXE) $(E2A_EXE) $(M2D_EXE) $(G_STATS_EXE) $(TRI_EXE) $(SRT_EXE)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(PKT_EXE): $(PKT_OBJ) 
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(TD_STATS_EXE): $(TD_STATS_OBJ) 
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(G_STATS_EXE): $(G_STATS_OBJ) 
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(E2D_EXE): $(E2D_OBJ)  ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(E2A_EXE): $(E2A_OBJ)  ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(M2D_EXE): $(M2D_OBJ)  ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(TRI_EXE): $(TRI_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(SRT_EXE): $(SRT_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

clean: 
	-rm -rf $(PKT_OBJ) $(PKT_EXE) $(TD_STATS_OBJ) $(TD_STATS_EXE) $(E2D_OBJ) $(E2D_EXE) $(M2D_EXE) $(M2D_OBJ) $(G_STATS_EXE) $(G_STATS_OBJ) $(TRI_OBJ) $(TRI_EXE) $(SRT_OBJ) $(SRT_EXE) $(E2A_EXE) $(E2A_OBJ)

doc: $(DOXYFILE)
	doxygen $(DOXYFILE)

clean_doc:
	rm -rf $(DOX_DIR)
