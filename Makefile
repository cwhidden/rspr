CC=g++
CC64=x86_64-w64-mingw32-g++
CFLAGS=-O2
C64FLAGS=$(CFLAGS)
BOOST_GRAPH=-lboost_graph-mt
all:
	$(CC) $(CFLAGS) -o rspr rspr.cpp
.PHONY: test
.PHONY: debug
.PHONY: hyb
test:
	./rspr.exe -fpt <test_trees/trees2.txt;
	./rspr.exe -fpt <test_trees/trees3.txt;
	./rspr.exe -fpt <test_trees/trees4.txt;
	./rspr.exe -fpt <test_trees/trees5.txt;
	./rspr.exe -fpt <test_trees/trees6.txt;
bb-test:
	./rspr.exe -bb <test_trees/trees1.txt;
	./rspr.exe -bb <test_trees/trees2.txt;
	./rspr.exe -bb <test_trees/trees3.txt;
	./rspr.exe -bb <test_trees/trees4.txt;
	./rspr.exe -bb <test_trees/trees5.txt;
	./rspr.exe -bb <test_trees/trees6.txt;
debug:
	#make "CFLAGS= -g -pg -fprofile-arcs -ftest-coverage" all
	CFLAGS=" -g -pg -fprofile-arcs -ftest-coverage $(CFLAGS)"
hyb:
	$(CC64) $(C64FLAGS) $(BOOST_GRAPH) -o hyb hyb.cpp
w32:
	$(CC) $(CFLAGS) -o rspr rspr.cpp
w64:
	$(CC64) $(C64FLAGS) -o rspr rspr.cpp
