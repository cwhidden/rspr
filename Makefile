CC=g++
CFLAGS=-O2
all:
	$(CC) $(CFLAGS) -o rspr rspr.cpp
.PHONY: test
.PHONY: debug
test:
	./rspr.exe -fpt <test_trees/trees1.txt;
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
	make "CFLAGS= -g -pg -fprofile-arcs -ftest-coverage" all
