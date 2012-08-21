include ../make.inc

VIZ_SRC = ./src/viz_main.cpp ./src/viz.cpp
VIZ_EXE = $(INDDGO_BIN)/td_viz

INC_DIR=  -I./inc -I$(GRAPH)/inc -I$(TREE)/inc $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(UTHASH_INCDIR) $(GMP_INCDIR) $(MADNESS_INCDIR) $(PARMETIS_INCDIR)

LIBS = -lm $(INDDGO_LIB) -lgraphd -ltreed $(METIS_LIB) $(SSPARSE_LIB) $(ARPACK_LIB) $(GMP_LIB) 


VIZ_OBJ=$(VIZ_SRC:.cpp=.o)

all: $(VIZ_EXE)

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $(INC_DIR)  $< -o $@

$(VIZ_EXE): $(VIZ_OBJ)
	$(CXX) $(CPPFLAGS) $(CFLAGS) $^ $(INC_DIR) $(LIBS) $(LDFLAGS) -o $(VIZ_EXE)

clean: 
	-rm -rf $(VIZ_OBJ)
	-rm -rf $(VIZ_EXE)

