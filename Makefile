CC=g++
CC64=x86_64-w64-mingw32-g++
CFLAGS=-O2
C64FLAGS=$(CFLAGS)
BOOST_GRAPH=-lboost_graph-mt
BOOST_ANY=
LFLAGS=$(BOOST_GRAPH) $(BOOST_ANY)
DEBUGFLAGS=-g -pg
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
	$(CC) $(LFLAGS) $(DEBUGFLAGS) -o rspr rspr.cpp
hyb:
	$(CC64) $(LFLAGS) $(C64FLAGS) -o hyb hyb.cpp
w32:
	$(CC) $(LFLAGS) $(CFLAGS) -o rspr rspr.cpp
w64:
	$(CC64) $(LFLAGS) $(C64FLAGS) -o rspr rspr.cpp
