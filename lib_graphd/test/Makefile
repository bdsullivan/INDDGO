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

#INCLUDES = -I ../inc -I. -I $(GTEST_DIR)/include -I$(GTEST_DIR) -I.. $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(PARMETIS_INCDIR)
INCLUDES = -I../inc -I. -I.. $(METIS_INCLUDES) $(SSPARSE_INCDIR) $(PARMETIS_INCDIR)
LIBDIR= -L. -L.. $(METIS_LIB_DIR) $(SSPARSE_LIB) -L../lib $(INDDGO_LIB)
#LIB=  $(MPILIBS) -lpthread -lgraphd -lgtest $(METIS_LIB) 
#LIB=  $(MPILIBS) -lpthread -lgtest $(METIS_LIB) 
LIB=  $(MPILIBS) -lpthread -lgraphd -lgtest -lgtest_main $(METIS_LIB) 
LFLAGS = $(CFLAGS)
#LIBGRAPHDECOMP=$(GRAPH)/lib/libgraphd.a
#LIBGRAPHDECOMP=$(SRC_DIR)/lib/libgraphd.a

#----changes are below this line----------------------------------------------

EXECUTABLES = $(LIBGTEST) testDimacsReader testGraphReader testGraphCreatorFile \
testGraphReaderWriterFactory testGraph testWeightedGraph testMutableGraph \
testWeightedMutableGraph testAdjMatrixGraphWriter testAdjMatrixGraphReader \
testDIMACSGraphWriter testGraphVizGraphWriter testGraphProperty testGraphUtil\
testGraphEOUtil testMetisGraphWriter testMetisGraphReader testParMetis


SRC=DIMACSReaderTest.cpp GraphReaderTest.cpp GraphCreatorFileTest.cpp \
GraphReaderWriterFactoryTest.cpp GraphTest.cpp WeightedGraphTest.cpp MutableGraphTest.cpp\
WeightedMutableGraphTest.cpp AdjMatrixGraphWriterTest.cpp AdjMatrixGraphReaderTest.cpp\
DIMACSGraphWriterTest.cpp GraphVizGraphWriterTest.cpp GraphPropertyTest.cpp\
GraphUtilTest.cpp GraphEOUtilTest.cpp MetisGraphWriterTest.cpp MetisGraphReaderTest.cpp\
ParmetisTest.cpp

ifeq ($(HAS_GTEST), 1)

all : ${EXECUTABLES}

.cpp.o:
	$(CXX) $(INCLUDES) $(CFLAGS) -c $<

testDimacsReader: DIMACSReaderTest.o 
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphReader: GraphReaderTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphCreatorFile: GraphCreatorFileTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphReaderWriterFactory: GraphReaderWriterFactoryTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraph: GraphTest.o
	#${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)
	${CXX} ${LFLAGS} -o bin/$@   $< ${INCLUDES} $(LIBDIR) $(LIB)
testWeightedGraph: WeightedGraphTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testMutableGraph: MutableGraphTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testWeightedMutableGraph: WeightedMutableGraphTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testAdjMatrixGraphWriter: AdjMatrixGraphWriterTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testAdjMatrixGraphReader: AdjMatrixGraphReaderTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testDIMACSGraphWriter: DIMACSGraphWriterTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphVizGraphWriter: GraphVizGraphWriterTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphProperty: GraphPropertyTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphUtil: GraphUtilTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphEOUtil: GraphEOUtilTest.o
	${CXX} ${LFLAGS} -o bin/$@ $< $(LIBGRAPHDECOMP) $(PARMETISOBJ)  ${INCLUDES} $(LIB) $(LIBDIR)

testMetisGraphWriter: MetisGraphWriterTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testMetisGraphReader: MetisGraphReaderTest.o
	${CXX} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testParMetis: ParmetisTest.o 
	${CXX} ${LFLAGS} -o bin/$@ $^ $(LIBGRAPHDECOMP) $(PARMETISOBJ) ${INCLUDES} $(LIB) $(LIBDIR)
#---- no change below this line--------------------------------------------------

#gtest-all.o: ${GTEST_DIR}/src/gtest-all.cc
#	$(CXX) -c $< $(CXXFLAGS)	
#
#gtest_main.o: ${GTEST_DIR}/src/gtest_main.cc
#	$(CXX) -c $< $(CXXFLAGS)	
#
#$(LIBGTEST): gtest-all.o gtest_main.o
#	ar -rv ../lib/$(LIBGTEST) $^

test:
	sh run_test.sh

endif

clean:
	rm -f *.o \#* PI* ${EXECUTABLES} core *~ bin/*
