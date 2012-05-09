include ../../make.inc

METIS_PROG_DIR=$(METIS_SRC_DIR)/programs

ifeq ($(HAS_METIS),1)
CFILES = $(METIS_PROG_DIR)/smbfactor.c 
endif


LIBNAME = libgraphd.a

GVIZ=../bin/gviz
# include directories
INCLUDES = -I../inc/ -I/usr/local/include $(METIS_INCLUDES) -I$(AMD)/Include -I$(AMD)/../UFconfig $(PARMETIS_INCDIR) 
# library paths
LIBDIR= -L$(INDDGO_LIB) -L/usr/local/lib -L$(METIS_LIB_DIR) -L$(AMD)/Lib  -L$(PARMETIS_LIBDIR)
LIBS = $(LIBDIR) $(AMD_LIB) -lm $(METIS_LIB) $(PARMETIS_LIB)
# source
SRC := Log.cpp GraphReader.cpp DIMACSGraphReader.cpp AdjMatrixGraphReader.cpp Node.cpp GraphReaderWriterFactory.cpp \
	GraphCreator.cpp GraphCreatorFile.cpp Graph.cpp WeightedGraph.cpp MutableGraph.cpp WeightedMutableGraph.cpp GraphProperties.cpp\
	GraphUtil.cpp Debug.cpp GraphEOUtil.cpp GraphDisplay.cpp RndNumGen.cpp Util.cpp AdjMatrixGraphWriter.cpp GraphWriter.cpp\
	DIMACSGraphWriter.cpp GraphVizGraphWriter.cpp MetisGraphWriter.cpp MetisGraphReader.cpp
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

all: $(INDDGO_LIB)/$(LIBNAME)  $(GVIZ)

%.o: $(METIS_PROG_DIR)/%.c
	gcc $(INCLUDES) $(CFLAGS) -U__GNUC__ -c $< 

.cpp.o:
	$(CC) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) -c $< 	

$(INDDGO_LIB)/$(LIBNAME): $(OBJ) $(POBJ)
	ar -rv $(INDDGO_LIB)/$(LIBNAME) $^

#builds the library to the lib directory inside lib_graphd
locallib: $(INDDGO_LIB)/$(LIBNAME)
	cp $(INDDGO_LIB)/$(LIBNAME) ../lib/$(LIBNAME)


$(GVIZ): gviz.cpp
	$(CC) -I../inc $(CFLAGS) $< -o $@ $(INCLUDES) $(LIBS) -L$(INDDGO_LIB) -lgraphd

clean:
	rm -f *.o \#* PI* ${EXECUTABLES} core *~ $(OBJ) $(OUT) $(INDDGO_LIB)/$(LIBNAME) ../lib/$(LIBNAME) $(GVIZ)

# makefile tutorial: http://www.cs.duke.edu/~ola/courses/programming/Makefiles/node1.html










