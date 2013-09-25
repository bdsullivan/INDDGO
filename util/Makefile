include ../make.inc
DOX_DIR = ./doc/
DOXYFILE = ./Doxyfile

INC_DIR = -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(UTHASH_INCDIR) $(MADNESS_INCDIR) $(GMP_INCDIR) $(PARMETIS_INCDIR) $(SLEPC_INCDIR) $(PETSC_INCDIR)

LIBS = -lm $(INDDGO_LIB) -ltreed -lgraphd $(SSPARSE_LIB) $(ARPACK_LIB) $(GMP_LIB) $(PARMETIS_LIB) $(METIS_LIB) $(SLEPC_LIB) $(PETSC_LIB)

PKT_SRCS =  ./src/gen_pkt.cpp 
PKT_OBJ =  src/orbtimer.o $(PKT_SRCS:.cpp=.o)
PKT_EXE = $(INDDGO_BIN)/gen_pkt

TD_STATS_SRCS =  ./src/td_stats.cpp 
TD_STATS_OBJ =  $(TD_STATS_SRCS:.cpp=.o)
TD_STATS_EXE = $(INDDGO_BIN)/td_stats

G_STATS_SRCS =  ./src/graph_stats.cpp 
G_STATS_OBJ =  src/orbtimer.o $(G_STATS_SRCS:.cpp=.o)
G_STATS_EXE = $(INDDGO_BIN)/graph_stats

GCONVERT_SRCS = ./src/gconvert.cpp
GCONVERT_OBJ = src/orbtimer.o $(GCONVERT_SRCS:.cpp=.o)
GCONVERT_EXE = $(INDDGO_BIN)/gconvert

TRI_SRCS =  ./src/count-tri.cpp 
TRI_OBJ =  src/orbtimer.o $(TRI_SRCS:.cpp=.o)
TRI_EXE = $(INDDGO_BIN)/count-tri

SRT_SRCS =  ./src/shortest-paths.cpp 
SRT_OBJ =  $(SRT_SRCS:.cpp=.o)
SRT_EXE = $(INDDGO_BIN)/shortest-paths

ECC_SRCS = ./src/v_eccentricity.cpp
ECC_OBJ = $(ECC_SRCS:.cpp=.o)
ECC_EXE = $(INDDGO_BIN)/v_eccentricity

EXP_SRCS = ./src/expansion.cpp
EXP_OBJ = $(EXP_SRCS:.cpp=.o)
EXP_EXE = $(INDDGO_BIN)/expansion

all: $(PKT_EXE) $(TD_STATS_EXE) $(GCONVERT_EXE) $(E2A_EXE) $(M2D_EXE) $(G_STATS_EXE) $(TRI_EXE) $(SRT_EXE) $(ECC_EXE) $(EXP_EXE)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(PKT_EXE): $(PKT_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(TD_STATS_EXE): $(TD_STATS_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(G_STATS_EXE): $(G_STATS_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(GCONVERT_EXE): $(GCONVERT_OBJ)  ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(TRI_EXE): $(TRI_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(SRT_EXE): $(SRT_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@ 

$(ECC_EXE): $(ECC_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@

$(EXP_EXE): $(EXP_OBJ) ../lib/libgraphd.a
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $@

clean: 
	-rm -rf $(PKT_OBJ) $(PKT_EXE) $(TD_STATS_OBJ) $(TD_STATS_EXE) $(GCONVERT_OBJ) $(GCONVERT_EXE) $(M2D_EXE) $(M2D_OBJ) $(G_STATS_EXE) $(G_STATS_OBJ) $(TRI_OBJ) $(TRI_EXE) $(SRT_OBJ) $(SRT_EXE) $(ECC_OBJ) $(ECC_EXE) $(EXP_OBJ) $(EXP_EXE) $(E2A_EXE) $(E2A_OBJ)

doc: $(DOXYFILE)
	doxygen $(DOXYFILE)

clean_doc:
	rm -rf $(DOX_DIR)
