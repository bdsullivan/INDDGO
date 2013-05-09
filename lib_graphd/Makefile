include ../make.inc

DOX_DIR = ./doc/
DOXYFILE = ./Doxyfile

.PHONY: test

all:
	cd src;	make all;

test:
	cd test; 
	make all;
	./run_test.sh;

doc: $(DOXYFILE)
	doxygen $(DOXYFILE)

clean:
	cd src;	make clean;
	cd ..;
	cd test;make clean;

clean_doc:
	rm -rf $(DOX_DIR)

# makefile tutorial: http://www.cs.duke.edu/~ola/courses/programming/Makefiles/node1.html
