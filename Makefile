CC=g++
CC64=CC
CFLAGS=-O3
OMPFLAGS=-fopenmp
C64FLAGS=$(CFLAGS)
BOOST_GRAPH=-lboost_graph-mt
BOOST_ANY=-L/lib/libboost*
LFLAGS=#$(BOOST_GRAPH) $(BOOST_ANY)
DEBUGFLAGS=-g -O0 -std=c++0x
PROFILEFLAGS=-pg
OBJS=rspr spr_supertree
all: $(OBJS)

rspr: rspr.cpp *.h
	$(CC) $(CFLAGS) -o rspr rspr.cpp
spr_supertree: spr_supertree.cpp *.h
	$(CC) $(CFLAGS) -o spr_supertree spr_supertree.cpp


.PHONY: test
.PHONY: debug
.PHONY: hyb
.PHONY: profile

test: test.cpp *.h
	$(CC) $(lFLAGS) $(DEBUGFLAGS) -o test test.cpp

bb-test:
	./rspr.exe -bb <test_trees/trees1.txt;
	./rspr.exe -bb <test_trees/trees2.txt;
	./rspr.exe -bb <test_trees/trees3.txt;
	./rspr.exe -bb <test_trees/trees4.txt;
	./rspr.exe -bb <test_trees/trees5.txt;
	./rspr.exe -bb <test_trees/trees6.txt;
debug:
	$(CC) $(LFLAGS) $(DEBUGFLAGS) -o rspr rspr.cpp
	$(CC) $(LFLAGS) $(DEBUGFLAGS) -o spr_supertree spr_supertree.cpp
profile:
	$(CC) $(LFLAGS) $(DEBUGFLAGS) $(PROFILEFLAGS) -o rspr rspr.cpp
	$(CC) $(LFLAGS) $(DEBUGFLAGS) $(PROFILEFLAGS) -o spr_supertree spr_supertree.cpp
hyb:
	$(CC64) $(LFLAGS) $(C64FLAGS) -o hyb hyb.cpp
w32:
	$(CC) $(LFLAGS) $(CFLAGS) -o rspr rspr.cpp
w64:
	$(CC64) $(LFLAGS) $(C64FLAGS) -o rspr rspr.cpp
	$(CC64) $(LFLAGS) $(C64FLAGS) -o spr_supertree spr_supertree.cpp
omp:
	$(CC) $(CFLAGS) $(OMPFLAGS) -o rspr-omp rspr.cpp
	$(CC) $(CFLAGS) $(OMPFLAGS) -o spr_supertree-omp spr_supertree.cpp
omp-debug:
	$(CC) $(LFLAGS) $(DEBUGFLAGS) $(OMPFLAGS) -o rspr-omp rspr.cpp
	$(CC) $(LFLAGS) $(DEBUGFLAGS) $(OMPFLAGS) -o spr_supertree-omp spr_supertree.cpp
