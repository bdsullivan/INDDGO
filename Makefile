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


include make.inc

MADLIB = ./madness/deploy/lib/libMADworld.a

.PHONY : deps all wis viz util cleandeps cleanwis cleanviz cleanmisc clean cleangraph cleantree testgraph testtree cleanutil
#all: wis viz util doc 
all: wis viz util
deps: $(MADLIB)
libs: graph tree ptree
test: testgraph testtree

#----------------------------------------------------------------
#  Targets for documentation
#
doc:
	test -d $(GRAPH)
	@($(CD) "$(GRAPH)";\
	$(MAKE) doc;\
	$(CD) ..;)
	test -d $(TREE)
	@($(CD) "$(TREE)";\
	$(MAKE) doc;\
	$(CD) ..;)
ifeq ($(HAS_MADNESS), 1)
	test -d $(PTREE)
	@($(CD) "$(PTREE)";\
	$(MAKE) doc;\
	$(CD) ..;)
endif


#----------------------------------------------------------------
#  Targets for decomposition libraries
#

graph:
	test -d $(GRAPH)
	@($(CD) "$(GRAPH)/src";\
	$(MAKE) ;\
	$(CD) ..;)

tree: deps
	test -d $(TREE)
	@($(CD) "$(TREE)";\
	$(MAKE) ;\
	$(CD) ..;)

ptree: 
ifeq ($(HAS_MADNESS), 1)
	test -d $(PTREE)
	@($(CD) "$(PTREE)";\
	$(MAKE) ;\
	$(CD) ..;)
endif

valtree:
	test -d $(TREE)
	@($(CD) "$(TREE)";\
	$(MAKE) td_valgrind ;\
	$(CD) ..;)


#----------------------------------------------------------------
#  Targets for weighted independent set
#
wis: libs deps
	test -d $(WIS)
	@($(CD) "$(WIS)";\
	$(MAKE) ;\
	$(CD) ..;)


#----------------------------------------------------------------
#  Targets for visualization
#

viz: graph tree 
	test -d $(VIZ)
	@($(CD) "$(VIZ)";\
	$(MAKE) ;\
	$(CD) ..;)

#----------------------------------------------------------------
#  Targets for utilities
#
util: graph tree
	test -d $(UTIL)
	@($(CD) "$(UTIL)";\
	$(MAKE) ;\
	$(CD) ..;)



#----------------------------------------------------------------
#  Target for installing INDDGO-distributed MADNESS runtime. 
#

$(MADLIB):
ifeq ($(HAS_MADNESS), 1)
	test -d $(MADNESS)
	@($(CD) "$(MADNESS)";\
	$(AUTOGEN);\
	$(CONFIGURE) --prefix=$(MADNESS_INSTALL_DIR);\
	$(MAKE) -j 4 libraries; $(MAKE) install;\
	$(CD) ..;)
endif


#----------------------------------------------------------------
#  Targets for test frameworks
#

testgraph:
ifeq ($(HAS_GTEST), 1)
	test -d $(GRAPH)/test
	@($(CD) "$(GRAPH)/test";\
	$(MAKE) ;\
	$(RUN_TEST) ;\
	$(CD) ..;)
else
	@(echo "Testing not available without HAS_GTEST enabled.")
endif

testtree:
ifeq ($(HAS_GTEST), 1)
	test -d $(TREE)/test
	@($(CD) "$(TREE)/test";\
	$(MAKE) ;\
	$(RUN_TEST) ;\
	$(CD) ..;)
else
	@(echo "Testing not available without HAS_GTEST enabled.")
endif


#----------------------------------------------------------------
#  Targets for cleanup
#

# Default cleanup for INDDGO-distributed source and libraries only.
clean:
	test -d $(SRC_DIR)
	@( \
	for f in $(GRAPH)/src $(GRAPH)/test $(PTREE) $(TREE) $(WIS) $(VIZ) $(UTIL); \
	do \
		$(CD) "$$f"; \
		$(MAKE) clean;\
		$(CD) ..; \
	done );
	find . -name '*~' -print0 | xargs -0 rm -f

cleandeps: cleanmad


# Targets for cleaning INDDGO-distributed subdirectories
cleanmad:
	test -d $(SRC_DIR)
	@($(CD) "$(MADNESS)";\
	$(MAKE) distclean;\
	$(CD) ..;)
	rm -f $(MADLIB)

cleangraph: 
	test -d $(SRC_DIR)
	@($(CD) "$(GRAPH)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleantree:
	test -d $(SRC_DIR)
	@($(CD) "$(TREE)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleanwis: 
	test -d $(SRC_DIR)
	@($(CD) "$(WIS)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleanviz: 
	test -d $(SRC_DIR)
	@($(CD) "$(VIZ)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleanutil: 
	test -d $(SRC_DIR)
	@($(CD) "$(UTIL)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

clean_doc:
	test -d $(SRC_DIR)
	@( \
	for f in $(GRAPH) $(PTREE) $(TREE) ; \
	do \
		$(CD) "$$f"; \
		$(MAKE) clean_doc;\
		$(CD) ..; \
	done );
