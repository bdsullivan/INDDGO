#
# This file is part of INDDGO.
#
# Copyright (C) 2012, Oak Ridge National Laboratory 
#
# This product includes software produced by UT-Battelle, LLC 
# under Contract No.  DE-AC05-00OR22725 with the Department of Energy. 
#
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the New BSD 3-clause software license (see LICENSE file). 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# LICENSE for more details.
#
#  For more information please contact the INDDGO developers at: 
#  inddgo-info@googlegroups.com
#
include ../../make.inc


#-----Macros---------------------------------
# set up compiler and options
#CXX = g++
#-----------------------------------------------------------------------------
# To build the test library

CXXFLAGS = $(INCLUDES) -g

INCLUDES = -I../inc -I. -I.. $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(PARMETIS_INCDIR)
LIBDIR= -L. -L.. $(METIS_LIB_DIR) $(SSPARSE_LIB) -L../lib $(INDDGO_LIB) $(MPILIB_DIR)
LIB=  $(MPILIBS) -lpthread -lgraphd -lgtest -lgtest_main $(METIS_LIB) 
LFLAGS = -g -O0

#----changes are below this line----------------------------------------------

#EXECUTABLES = $(LIBGTEST) testDimacsReader testGraphReader testGraphCreatorFile 
#EXECUTABLES = bin/testDimacsReader bin/testGraphReader bin/testGraphCreatorFile \
bin/testGraph bin/testWeightedGraph  \
bin/testVertexWeightedGraph bin/testAdjMatrixGraphWriter bin/testAdjMatrixGraphReader \
bin/testDIMACSGraphWriter bin/testGraphVizGraphWriter bin/testGraphProperty bin/testGraphUtil\
bin/testMetisGraphWriter bin/testMetisGraphReader

EXECUTABLES = bin/testGraphReader bin/testGraphCreatorFile \
bin/testGraph bin/testWeightedGraph  \
bin/testVertexWeightedGraph bin/testAdjMatrixGraphWriter bin/testAdjMatrixGraphReader \
bin/testDIMACSGraphWriter bin/testGraphVizGraphWriter bin/testGraphProperty bin/testGraphUtil\
bin/testMetisGraphWriter bin/testMetisGraphReader

ifeq ($(HAS_PARMETIS), 1)
EXECUTABLES = $(EXECUTABLES) testParMetis testGraphEOUtil
endif

SRC=DIMACSReaderTest.cpp GraphReaderTest.cpp GraphCreatorFileTest.cpp \
GraphTest.cpp WeightedGraphTest.cpp MutableGraphTest.cpp\
VertexWeightedGraphTest.cpp AdjMatrixGraphWriterTest.cpp AdjMatrixGraphReaderTest.cpp\
DIMACSGraphWriterTest.cpp GraphVizGraphWriterTest.cpp GraphPropertyTest.cpp\
GraphUtilTest.cpp GraphEOUtilTest.cpp MetisGraphWriterTest.cpp MetisGraphReaderTest.cpp\
ParmetisTest.cpp

ifeq ($(HAS_GTEST), 1)

all: bin ${EXECUTABLES} $(INDDGO_LIB_DIR)/libgraphd.a

bin:
	mkdir -p bin

.cpp.o:
	$(CXX) $(INCLUDES) $(CFLAGS) -O0 -c $<

bin/testDimacsReader: DIMACSReaderTest.o 
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testGraphReader: GraphReaderTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testGraphCreatorFile: GraphCreatorFileTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testGraph: GraphTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testWeightedGraph: WeightedGraphTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testMutableGraph: MutableGraphTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testVertexWeightedGraph: VertexWeightedGraphTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testAdjMatrixGraphWriter: AdjMatrixGraphWriterTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testAdjMatrixGraphReader: AdjMatrixGraphReaderTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testDIMACSGraphWriter: DIMACSGraphWriterTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testGraphVizGraphWriter: GraphVizGraphWriterTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testGraphProperty: GraphPropertyTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testGraphUtil: GraphUtilTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testMetisGraphWriter: MetisGraphWriterTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

bin/testMetisGraphReader: MetisGraphReaderTest.o
	${CXX} ${LFLAGS} -o $@   $< ${INCLUDES} $(LIBDIR) $(LIB)

ifeq ($(HAS_PARMETIS), 1)
bin/testGraphEOUtil: GraphEOUtilTest.o
	${CXX} ${LFLAGS} -o $@   $< $(PARMETISOBJ) ${INCLUDES} $(LIBDIR) $(LIB)

bin/testParMetis: ParmetisTest.o 
	${CXX} ${LFLAGS} -o $@   $^ $(PARMETISOBJ) ${INCLUDES} $(LIBDIR) $(LIB)
endif

#---- no change below this line--------------------------------------------------

test:
	sh run_test.sh

endif

clean:
	rm -f *.o \#* PI* ${EXECUTABLES} core *~ bin/*
