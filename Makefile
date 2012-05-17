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
all: wis viz util 
deps: $(MADLIB)
libs: graph tree ptree
test: testgraph testtree

#----------------------------------------------------------------
#  Targets for decomposition libraries
#

graph:
	@($(CD) "$(GRAPH)/src";\
	$(MAKE) ;\
	$(CD) ..;)

tree: deps
	@($(CD) "$(TREE)";\
	$(MAKE) ;\
	$(CD) ..;)

ptree: 
ifeq ($(HAS_MADNESS), 1)
	@($(CD) "$(PTREE)";\
	$(MAKE) ;\
	$(CD) ..;)
endif

valtree:
	@($(CD) "$(TREE)";\
	$(MAKE) td_valgrind ;\
	$(CD) ..;)


#----------------------------------------------------------------
#  Targets for weighted independent set
#
wis: libs deps
	@($(CD) "$(WIS)";\
	$(MAKE) ;\
	$(CD) ..;)


#----------------------------------------------------------------
#  Targets for visualization
#

viz: graph tree 
	@($(CD) "$(VIZ)";\
	$(MAKE) ;\
	$(CD) ..;)

#----------------------------------------------------------------
#  Targets for utilities
#
util: graph 
	@($(CD) "$(UTIL)";\
	$(MAKE) ;\
	$(CD) ..;)



#----------------------------------------------------------------
#  Target for installing INDDGO-distributed MADNESS runtime. 
#

$(MADLIB):
ifeq ($(HAS_MADNESS), 1)
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
	@($(CD) "$(GRAPH)/test";\
	$(MAKE) ;\
	$(RUN_TEST) ;\
	$(CD) ..;)
else
	@(echo "Testing not available without HAS_GTEST enabled.")
endif

testtree:
ifeq ($(HAS_GTEST), 1)
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
	@($(CD) "$(MADNESS)";\
	$(MAKE) distclean;\
	$(CD) ..;)
	rm -f $(MADLIB)

cleangraph: 
	@($(CD) "$(GRAPH)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleantree:
	@($(CD) "$(TREE)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleanwis: 
	@($(CD) "$(WIS)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleanviz: 
	@($(CD) "$(VIZ)";\
	  $(MAKE) clean;\
	  $(CD) ..;)

cleanutil: 
	@($(CD) "$(UTIL)";\
	  $(MAKE) clean;\
	  $(CD) ..;)






