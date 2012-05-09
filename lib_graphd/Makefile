include ../make.inc

all:
	cd src;	make all;

test:
	cd test; 
	make all;
	./run_test.sh;

clean:
	cd src;	make clean;
	cd ..;
	cd test;make clean;

# makefile tutorial: http://www.cs.duke.edu/~ola/courses/programming/Makefiles/node1.html
