################################################################################
rspr

################################################################################

Usage: rspr [OPTIONS]
Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of rooted binary trees from STDIN in newick format. Supports arbitrary labels.
The second tree may be multifurcating. 

Can also compare the first input tree to each other tree with -total or
compute a pairwise distance matrix with -pairwise.

Copyright 2009-2021 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
February 12, 2021
Version 1.3.1

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
            FPT algorithm. This is enabled by default.

-approx     Calculate just a linear-time 3-approximation of the rSPR distance

-split_approx
-split_approx x  Calculate the exact rSPR distance if it is k or less and
                 otherwise use the exponential-time approximation

-cluster_test   Use the cluster reduction to speed up the exact algorithm.
                This is enabled by default.

-total          Find the total SPR distance from the first input tree to
                the rest of the list of trees. Uses the other algorithm
                options as specified (including unrooted options).

*******************************************************************************
OPTIMIZATIONS
*******************************************************************************

These options control the use of optimized branching. All optimizations are
enabled by default. Specifying any subset of -cob, -cab, and -sc will use
just that subset of optimizations.

-allopt    Use -cob -cab -sc and a new set of improvements. This is the
           default
option

-noopt     Use 3-way branching for all FPT algorithms

-cob       Use "cut one b" improved branching

-cab       Use "cut all b" improved branching

-sc        Use "separate components" improved branching

*******************************************************************************
MULTIFURCATING COMPARISON OPTIONS
*******************************************************************************

-support x     Collapse bipartitions with less than x support
-length x      Collapse bipartitions with branch lengths less than or
                equal to x

*******************************************************************************
UNROOTED COMPARISON OPTIONS
*******************************************************************************

-unrooted   Compare the first input tree to each other input tree.
            Output the best found distance and agreement forest.
            This option can be used with gen_rooted_trees.pl to provide
            the rootings.
            Note that this option is a bit unintuitive to maintain
            compatibility with previous versions of rSPR.
            If -total or -pairwise analysis is used then there is no need
            to specify rootings.

-unrooted_min_approx    Compare the first input tree to each other input tree.
                        Run the exact algorithms on the pair with the
                        minimum approximate rspr distance

-simple_unrooted        Root the gene trees using
                        a bipartition balanced accuracy measure
                        (fast but potentially less accurate). Only
                        used with -total.

*******************************************************************************
PAIRWISE COMPARISON OPTIONS
*******************************************************************************

-pairwise
-pairwise a b
-pairwise a b c d        Compare each input tree to each other tree and output
                         the resulting SPR distance matrix. If -unrooted is
                         enabled this will compute the "best rooting" SPR
                         distance by testing each rooting of the trees. The
                         optional arguments a b c d compute only rows a-b and/or
                         columns c-d of the matrix.

-no-symmetric-pairwise   By default, -pairwise will ignore the symmetric lower
                         left triangle of the matrix. With this option the
                         lower triangle is filled in.

-pairwise_max x          Use with -pairwise to only compute distances at most x.
                         Larger values are output as -1. Very efficient for
                         small distances (e.g. 1-10).

*******************************************************************************
OTHER OPTIONS
*******************************************************************************
-cc         Calculate a potentially better approximation with a quadratic time
            algorithm

-q          Quiet; Do not output the input trees or approximation
*******************************************************************************

Example:
$ ./rspr < test_trees/trees2.txt
T1: ((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))))
T2: (((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)))

T1: ((((0,1),(2,3)),((4,5),(6,7))),(((8,9),(10,11)),((12,13),(14,15))))
T2: (((6,7),((0,(1,(13,4))),(2,3))),(((10,(5,11)),9),((12,(14,15)),8)))
approx F1: ((0,(2,3)),((9,(10,11)),(12,(14,15)))) 13 5 8 4 (6,7) 1
approx F2: ((0,(2,3)),(((10,11),9),(12,(14,15)))) 13 5 8 4 1 (6,7)
approx drSPR=6

C1_1: ((((0,1),(2,3)),((4,5),(6,7))),(((8,9),(10,11)),((12,13),(14,15))))
C1_2: (((6,7),((0,(1,(13,4))),(2,3))),(((10,(5,11)),9),((12,(14,15)),8)))
cluster approx drSPR=4

4
F1_1: ((((0,1),(2,3)),(6,7)),((9,(10,11)),(12,(14,15)))) 5 4 13 8
F1_2: (((6,7),((0,1),(2,3))),(((10,11),9),(12,(14,15)))) 13 5 4 8
cluster exact drSPR=4

F1: ((((1,2),(3,4)),(7,8)),((10,(11,12)),(13,(15,16)))) 6 5 14 9
F2: (((7,8),((1,2),(3,4))),(((11,12),10),(13,(15,16)))) 14 6 5 9
total exact drSPR=4

################################################################################

CONTACT INFORMATION

Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR

################################################################################

FILES


ClusterForest.h   Cluster Decomposition
Forest.h          Forest data structure
gen_rooted_trees.pl    Generate all rootings of an unrooted binary tree
gpl.txt           The GPL license
LCA.h             Compute LCAs of tree leaves
Makefile          Makefile
Node.h            Node data structure
README.txt        This README
rspr.h            Library to calculate rSPR distances between pairs of trees
rspr.cpp          Calculate rSPR distances between pairs or sets of trees
test_trees/       Folder of test tree pairs
SiblingPair.h     Sibling Pair class

################################################################################

INSTALLATION

rSPR is a command-line program written in C++. To use it, simply
compile rspr.cpp and execute the resulting program. On systems with
the g++ compiler and make program, the included make file will
compile rspr; simply run `make'.

################################################################################

INPUT

rSPR requires pairs of Newick format trees with arbitrary labels
as input. The first tree must be binary and rooted. The second tree
may be multifurcating and rooted. A sample Newick tree is shown below:

((1,2),(3,4),(5,6));

rSPR can also compare a rooted reference tree to an unrooted test tree.
First use gen_rooted_trees.pl to generate all rootings of the unrooted
test tree. Then use the -unrooted or -unrooted_min_approx options and
input the test tree and the set of rootings. rSPR will find the best
rooting of the test tree with the -unrooted option and guess the best 
rooting based on the approximation algorithm with the
-unrooted_min_approx option. Alternatively, the -total option with
the -unrooted or -unrooted_min_approx options will provide just the
distance. The -total option with -simple_unrooted will use a faster
biparition based measure to approximate the optimal rooting.

The -support x option can be used to collapse poorly supported branches
of the second tree.

With the -pairwise option, rSPR will compare each pair of input trees
and output the results as a distance matrix. To save time, only the
upper right triangle is output as the lower left triangle is symmetric.
Use the included fill_matrix program to fill in missing values or the
-no-symmetric-pairwise option to explicitly compute these values.
Optional arguments to -pairwise can be used to compute subsets of the
matrix (e.g. for partitioning computation over multiple processes).
The -pairwise_max x option can be used to quickly find trees with
SPR distance at most x when x is small (e.g. 1-10).


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

OUTPUT WITH CLUSTERING

/////////////////////

$ ./rspr < test_trees/cluster_test 
T1: (((x,((b1,b3),b2)),y),(f,(a,c)))
T2: (((x,y),f),((a,((b1,b2),b3)),c))

F1: (((0,((1,2),3)),4),(5,(6,7))) 
F2: (((0,4),5),((6,((1,3),2)),7)) 
approx drSPR=9


CLUSTERS
C1_1: ((1,2),3) 
C1_2: ((1,3),2) 
cluster approx drSPR=3

1 
F1_1: (1,2) 3 
F1_2: (1,2) 3 
cluster exact drSPR=1

C2_1: (((0,(1,2)),4),(5,(6,7))) 
C2_2: (((0,4),5),((6,(1,2)),7)) 
cluster approx drSPR=6

2 
F2_1: (5,(6,7)) (1,2) (0,4) 
F2_2: (5,(6,7)) (1,2) (0,4) 
cluster exact drSPR=2

F1: (f,(a,c)) b2 (b1,b3) (x,y) 
F2: (f,(a,c)) b2 (b1,b3) (x,y) 
total exact drSPR=3

/////////////////////

When clustering is enabled (as it is by default), each solved
cluster is displayed along with its approximate and exact distance in
an intermediate representation with labels mapped from 0-(N-1) where
N is the number of labels. The final agreement forest and distance
are output last.

################################################################################

OUTPUT WITH PAIRWISE

$ cat test_trees/big_test* | ./rspr -pairwise
0,46,0,46
,0,46,50
,,0,46
,,,0

$ cat test_trees/big_test* | ./rspr -pairwise | ./fill_matrix
0,46,0,46
46,0,46,50
0,46,0,46
46,50,46,0

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

In addition, this version contains new improvements that
give a bound of O(2^k n). This provides another significant improvement
and is provably correct so these options are also enabled by default.

For much larger trees, the -split_approx option will compute an
exponential time approximation of the distance that is exact for
small distances and generally within a few percent of the optimal 
distance otherwise.

When using the -unrooted option, the exact algorithms run in O(2^k n^2) time.

The cluster reduction improves the running time of the
algorithm to O(2^k n) time where k is the largest rSPR distance of
any cluster (as opposed to the full rSPR distance). This provides a large
speedup when the trees are clusterable.

With the -pairwise option on m rooted trees, the program takes O(m^2
2^k n) time, where k is the largest SPR distance computed. With -unrooted
this becomes O(m^2 2^k n^3). The -pairwise_max x option limits k to x, 
but does not use clustering and is slow for large distances.
-

NOTE: This is an exponential algorithm that exactly solves an NP-hard problem.
Thus the algorithms may not finish in a reasonable amount of time for large
rSPR distances (> 20 without optimizations and > 70 with optimizations).

################################################################################

REFERENCES

For more information on the algorithms see:

Whidden, C., Zeh, N., Beiko, R.G.  Fixed-Parameter and Approximation
Algorithms for Maximum Agreement Forests of Multifurcating Trees.
(Submitted). 2013. Preprint available at
http://arxiv.org/abs/1305.0512

Whidden, C., Zeh, N., Beiko, R.G.  Supertrees based on the subtree
prune-and-regraft distance. Syst. Biol. 63 (4): 566-581. 2014.
doi:10.1093/sysbio/syu023.

Whidden, C., Beiko, R.G., Zeh, N. Fixed-Parameter Algorithms for Maximum
Agreement Forests. SIAM Journal on Computing 42.4 (2013), pp. 1431-1466.
Available at http://epubs.siam.org/doi/abs/10.1137/110845045

Whidden, C. Efficient Computation of Maximum Agreement Forests and their
Applications. PhD Thesis. Dalhousie University, Canada. 2013. Available at
www.cs.dal.ca/~whidden

Whidden, C., Beiko, R.G., Zeh, N. Fast FPT Algorithms for Computing
Rooted Agreement Forests: Theory and Experiments. Experimental Algorithms.
Ed. by P. Festa. Vol. 6049. Lecture Notes in Computer Science. Springer
Berlin Heidelberg, 2010, pp. 141-153. Available at
http://link.springer.com/chapter/10.1007/978-3-642-13193-6_13

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

Whidden, C., Beiko, R.G., Zeh, N.  Computing the SPR Distance of Binary
Rooted Trees in O(2^k n) Time. (In Preparation). 2013.

Whidden, C., Beiko, R.G. Zeh, N.  Fixed-Parameter and Approximation
Algorithms for Maximum Agreement Forests of Multifurcating Trees.
(Submitted). 2013.

Whidden, C., Zeh, N., Beiko, R.G.  Supertrees based on the subtree
prune-and-regraft distance. Syst. Biol. 63 (4): 566-581. 2014.
doi:10.1093/sysbio/syu023.

Whidden, C., Beiko, R.G., Zeh, N. Fixed-Parameter Algorithms for Maximum
Agreement Forests. SIAM Journal on Computing 42.4 (2013), pp. 1431-1466.
Available at http://epubs.siam.org/doi/abs/10.1137/110845045

Whidden, C., Beiko, R.G., Zeh, N. Fast FPT Algorithms for Computing
Rooted Agreement Forests: Theory and Experiments. Experimental Algorithms.
Ed. by P. Festa. Vol. 6049. Lecture Notes in Computer Science. Springer
Berlin Heidelberg, 2010, pp. 141-153. Available at
http://link.springer.com/chapter/10.1007/978-3-642-13193-6_13

Whidden, C., Zeh, N. A Unifying View on Approximation and FPT of
Agreement Forests. In: WABI 2009. LNCS, vol. 5724, pp. 390.401.
Springer-Verlag (2009).

################################################################################

