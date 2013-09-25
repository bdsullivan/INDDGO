include ../../make.inc

METIS_PROG_DIR=$(METIS_SRC_DIR)/programs

ifeq ($(HAS_METIS),1)
CFILES = $(METIS_PROG_DIR)/smbfactor.c
endif

LIBNAME = libgraphd.a

GVIZ=../bin/gviz
# include directories
INCLUDES = -I../inc/ -I/usr/local/include $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(PARMETIS_INCDIR) $(PETSC_INCDIR) $(SLEPC_INCDIR)
# libraries
LIBS = -L/usr/local/lib $(INDDGO_LIB) $(SSPARSE_LIB) -lm $(METIS_LIB) $(PARMETIS_LIB) $(SLEPC_LIB) $(PETSC_LIB)
# source
SRC := Log.cpp GraphReader.cpp Node.cpp \
	GraphCreator.cpp GraphCreatorFile.cpp Graph.cpp VertexWeightedGraph.cpp GraphProperties.cpp \
	GraphUtil.cpp Debug.cpp GraphEOUtil.cpp GraphDisplay.cpp RndNumGen.cpp Util.cpp GraphWriter.cpp \
	GraphReader.cpp

# source files.
OBJ = $(SRC:.cpp=.o) 
POBJ =  
#compile smbfactor only if using metis.
ifeq ($(HAS_METIS),1)
OBJ = $(SRC:.cpp=.o)  smbfactor.o
endif

ifeq ($(HAS_PARMETIS),1)
OBJ = 
POBJ = $(SRC:.cpp=.o)  smbfactor.o $(PARMETISOBJ)
endif

.SUFFIXES: .cpp .c
.PHONY: all libgraphd

libgraphd: $(INDDGO_LIB_DIR)/$(LIBNAME)

all: libgraphd  $(GVIZ)

#This is a problem if they're not using gcc 
%.o: $(METIS_PROG_DIR)/%.c
	$(CC) $(INCLUDES) $(CFLAGS) -U__GNUC__ -c $< 

.cpp.o:
	$(CXX) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) -c $< 	

$(INDDGO_LIB_DIR)/$(LIBNAME): $(OBJ) $(POBJ)
	$(AR) -rv $(INDDGO_LIB_DIR)/$(LIBNAME) $^

#builds the library to the lib directory inside lib_graphd
locallib: $(INDDGO_LIB_DIR)/$(LIBNAME)
	cp $(INDDGO_LIB_DIR)/$(LIBNAME) ../lib/$(LIBNAME)

$(GVIZ): gviz.o 
	$(CXX) -I../inc $(CFLAGS) $< -o $@ $(INCLUDES) -lgraphd $(LIBS) $(LDFLAGS)

clean:
	rm -f *.o \#* PI* ${EXECUTABLES} core *~ $(OBJ) $(OUT) $(INDDGO_LIB_DIR)/$(LIBNAME) ../lib/$(LIBNAME) $(GVIZ)

# makefile tutorial: http://www.cs.duke.edu/~ola/courses/programming/Makefiles/node1.html










