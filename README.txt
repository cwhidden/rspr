################################################################################
rspr

################################################################################

Usage: rspr [OPTIONS]
Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of rooted binary trees from STDIN in newick format. By default, computes a
3-approximation of the rSPR distance. Supports arbitrary labels.

Copyright 2009-2010 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 22, 2010
Version 1.01

This file is part of rspr.

rspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rspr.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
ALGORITHM
*******************************************************************************

These options control what algorithm is used

-fpt        Calculate the exact rSPR distance with an FPT algorithm

-bb         Calculate the exact rSPR distance with a branch-and-bound
            FPT algorithm. This is the default option.

-approx		Calculate just a linear -time 3-approximation of the rSPR distance

*******************************************************************************
OPTIMIZATIONS
*******************************************************************************

These options control the use of optimized branching. All optimizations are
enabled by default. Specifying any subset of -cob, -cab, and -sc will use
just that subset of optimizations.

-allopt		Use -cob -cab -sc. This is the default option

-noopt		Use 3-way branching for all FPT algorithms

-cob		Use "cut one b" improved branching

-cab		Use "cut all b" improved branching

-sc			Use "separate components" improved branching

*******************************************************************************
UNROOTED COMPARISON OPTIONS
*******************************************************************************

-unrooted   Compare the first input tree to each other input tree.
            Output the best found distance and agreement forest
-unrooted_min_approx    Compare the first input tree to each other input tree.
                        Run the exact algorithms on the pair with the
                        minimum approximate rspr distance

*******************************************************************************
OTHER OPTIONS
*******************************************************************************
-cc         Calculate a potentially better approximation with a quadratic time
            algorithm

-q          Quiet; Do not output the input trees or approximation
*******************************************************************************

Example:
$ ./rspr.exe -fpt <test_trees/trees2.txt
T1: ((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))))
T2: (((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)))

F1: (((1,2),(3,4)),(7,8)) 14 13 5 12 11 6 9 10 (15,16)
F2: ((7,8),((1,2),(3,4))) 14 5 13 12 6 11 9 (15,16) 10
approx drSPR=9

3 4
F1: ((((1,2),(3,4)),(7,8)),((10,(11,12)),(13,(15,16)))) 14 6 9 5
F2: (((7,8),((1,2),(3,4))),(((11,12),10),(13,(15,16)))) 14 6 9 5
exact drSPR=4

################################################################################

CONTACT INFORMATION

Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR

################################################################################

FILES

Forest.h				Forest data structure
gen_rooted_trees.pl		Generate all rootings of an unrooted binary tree
gpl.txt					The GPL license
Makefile				Makefile
Node.h					Node data structure
README.txt				This README
rspr.cpp				Calculate rSPR distances between pairs of trees
test_trees/				Folder of test tree pairs

################################################################################

INSTALLATION

rSPR is a command-line program written in C++. To use it, simply
compile rspr.cpp and execute the resulting program. On systems with
the g++ compiler and make program, the included make file will
compile rspr; simply run `make'.

################################################################################

INPUT

rSPR requires pairs of Newick format binary rooted trees with arbitrary labels
as input.  A sample Newick tree is shown below:

((1,2),(3,4),(5,6));

rSPR can also compare a rooted reference tree to an unrooted test tree.
First use gen_rooted_trees.pl to generate all rootings of the unrooted
test tree. Then use the -unrooted or -unrooted_min_approx options and
input the test tree and the set of rootings. rSPR will find the best
rooting of the test tree with the -unrooted option and guess the best 
rooting based on the approximation algorithm with the
-unrooted_min_approx option.

################################################################################

OUTPUT

rspr writes to standard output.

A sample command line and output are shown below:

/////////////////////

$ ./rspr < test_trees/trees2.txt
T1: ((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))))
T2: ((((3,4),(8,(2,((11,12),1)))),((15,16),(7,(6,5)))),(14,((10,13),9)))

F1: ((3,4),(5,6)) 13 14 10 (11,12) 9 1 8 7 2 (15,16)
F2: ((3,4),(6,5)) 13 10 14 (11,12) 1 9 8 2 7 (15,16)
approx drSPR=12

4
F1: ((((1,2),(3,4)),((5,6),7)),((9,10),14)) 13 (11,12) 8 (15,16)
F2: ((((3,4),(2,1)),(7,(6,5))),(14,(10,9))) 13 (11,12) 8 (15,16)
exact BB drSPR=4

/////////////////////

The first set of lines show the input trees. The second set of lines are the
approximate agreement forests and the corresponding approximate rSPR distance.
The third set of lines are the maximum agreement forests and the corresponding
exact rSPR distance. When calculating exact distances, the distance
currently being considered is printed on the first line of this section.

Each component of an agreement forest corresponds to an rSPR operation. 
The set of rSPR operations required to turn one tree into the other can
be found by applying rSPR operations that move these components to their
correct place in the other tree.

An agreement forest may contain p (rho) as a component. This represents
the root of the trees and indicates that an extra rSPR operation is
required to correctly root the tree.

################################################################################

EFFICIENCY

The 3-approximation algorithm runs in O(n) time, where n is the number of
leaves in the trees.

The unoptimized FPT and branch-and-bound algorithms run in O(3^k n) time, where
k is the rSPR distance and n is the number of leaves in the trees. The
branch-and-bound algorithm should be significantly faster in practice.

Using all 3 of the -cob -cab and -sc optimizations improves the running times of
the algorithms to O(2.42^k n) time. This provides a significant improvement in
practice and is provably correct, thus this is the default.

When using the -unrooted option, the exact algorithms run in O(3^k n^2) time.

NOTE: This is an exponential algorithm that exactly solves an NP-hard problem.
Thus the algorithms may not finish in a reasonable amount of time for large
rSPR distances (> 20 without optimizations and > 50 with optimizations).

################################################################################

REFERENCES

For more information on the algorithms see:

Whidden, C., Beiko, R.G., Zeh, N. Fast FPT Algorithms for Computing
Rooted Agreement Forests: Theory and Experiments (Extended Abstract)
Accepted to SEA 2010.

Whidden, C., Zeh, N. A Unifying View on Approximation and FPT of
Agreement Forests. In: WABI 2009. LNCS, vol. 5724, pp. 390.401.
Springer-Verlag (2009). Available at
http://www.springerlink.com/content/n56q2846v645p655/

Whidden, C. A Unifying View on Approximation and FPT of Agreement Forests.
Masters Thesis. Dalhousie University, Canada. 2009. Available at
www.cs.dal.ca/~whidden

################################################################################

CITING rSPR

If you use rSPR in your research, please cite:
Whidden, C., Beiko, R.G., Zeh, N. Fast FPT Algorithms for Computing
Rooted Agreement Forests: Theory and Experiments (Extended Abstract)
Accepted to SEA 2010.

Whidden, C., Zeh, N. A Unifying View on Approximation and FPT of
Agreement Forests. In: WABI 2009. LNCS, vol. 5724, pp. 390.401.
Springer-Verlag (2009).

################################################################################
