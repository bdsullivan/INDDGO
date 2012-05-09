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
#CC = g++
#-----------------------------------------------------------------------------
# To build the test library

CXXFLAGS = $(INCLUDES) -g

INCLUDES = -I ../inc -I. -I $(GTEST_DIR)/include -I$(GTEST_DIR) -I.. $(METIS_INCLUDES) -I$(AMD)/Include -I$(AMD)/../UFconfig $(PARMETIS_INCDIR)
LIBDIR= -L. -L.. $(METIS_LIB_DIR) -L$(AMD)/Lib -L../lib -L $(INDDGO_LIB)
LIB=  $(MPILIBS) -lpthread -lgraphd -lgtest $(METIS_LIB) $(AMD_LIB)
LFLAGS = $(CFLAGS)
LIBGRAPHDECOMP=$(GRAPH)/lib/libgraphd.a

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
	$(CC) $(INCLUDES) $(CFLAGS) -c $<

testDimacsReader: DIMACSReaderTest.o 
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphReader: GraphReaderTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphCreatorFile: GraphCreatorFileTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphReaderWriterFactory: GraphReaderWriterFactoryTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraph: GraphTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)
testWeightedGraph: WeightedGraphTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testMutableGraph: MutableGraphTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testWeightedMutableGraph: WeightedMutableGraphTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testAdjMatrixGraphWriter: AdjMatrixGraphWriterTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testAdjMatrixGraphReader: AdjMatrixGraphReaderTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testDIMACSGraphWriter: DIMACSGraphWriterTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphVizGraphWriter: GraphVizGraphWriterTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphProperty: GraphPropertyTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphUtil: GraphUtilTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testGraphEOUtil: GraphEOUtilTest.o
	${CC} ${LFLAGS} -o bin/$@ $< $(LIBGRAPHDECOMP) $(PARMETISOBJ)  ${INCLUDES} $(LIB) $(LIBDIR)

testMetisGraphWriter: MetisGraphWriterTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testMetisGraphReader: MetisGraphReaderTest.o
	${CC} ${LFLAGS} -o bin/$@   $(LIBGRAPHDECOMP) $< ${INCLUDES} $(LIB) $(LIBDIR)

testParMetis: ParmetisTest.o 
	${CC} ${LFLAGS} -o bin/$@ $^ $(LIBGRAPHDECOMP) $(PARMETISOBJ) ${INCLUDES} $(LIB) $(LIBDIR)
#---- no change below this line--------------------------------------------------

gtest-all.o: ${GTEST_DIR}/src/gtest-all.cc
	$(CC) -c $< $(CXXFLAGS)	

gtest_main.o: ${GTEST_DIR}/src/gtest_main.cc
	$(CC) -c $< $(CXXFLAGS)	

$(LIBGTEST): gtest-all.o gtest_main.o
	ar -rv ../lib/$(LIBGTEST) $^

test:
	sh run_test.sh

endif

clean:
	rm -f *.o \#* PI* ${EXECUTABLES} core *~ bin/*
