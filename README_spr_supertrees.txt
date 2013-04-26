################################################################################
spr_supertree

################################################################################

Usage: spr_supertree [OPTIONS]
       spr_supertree-omp [OPTIONS]
Calculate supertrees that minimize the SPR distance from the input
trees. By default calculates a rooted SPR supertree from a list
of rooted binary trees from STDIN in newick format.

Copyright 2011-2013 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/SPR_Supertrees
April 26, 2013
Version 1.2.0

This file is part of spr_supertrees.

spr_supertrees is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

spr_supertrees is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with spr_supertrees.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
ALGORITHM
*******************************************************************************

These options control what algorithm is used to determine the SPR distance
from the supertree to the input trees

-fpt        Calculate the exact rSPR distance with an FPT algorithm

-bb         Calculate the exact rSPR distance with a branch-and-bound
            FPT algorithm. This is the default option.

-approx		Calculate just a linear -time 3-approximation of the rSPR distance

-max k	Calculate the exact rSPR distance if it is k or less and
				otherwise use the 3-approximation

*******************************************************************************
UNROOTED COMPARISON OPTIONS
*******************************************************************************

-unrooted   Compare the supertree to each rooting of the input trees.
            Use the best found distance.

-unrooted_min_approx    Compare the supertree to each rooting of the
												input trees.
                        Run the exact algorithm on the rooting with the
                        minimum approximate rspr distance

*******************************************************************************
OTHER OPTIONS
*******************************************************************************
-cc         Calculate a potentially better approximation with a quadratic time
            algorithm

-valid_trees    Output the set of trees that appear valid

-multi_trees    Output the set of multifurcating or invalid trees

################################################################################

CONTACT INFORMATION

Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/SPR_Supertrees

################################################################################

FILES

ClusterForest.h   Cluster Decomposition
Forest.h			  	Forest data structure
gen_rooted_trees.pl		Generate all rootings of an unrooted binary tree
gpl.txt				  	The GPL license
LCA.h             Compute LCAs of tree leaves
Makefile			  	Makefile
Node.h				  	Node data structure
README.txt				This README
rspr.h			    	Calculate rSPR distances between pairs of trees
spr_supertrees    Compute supertrees that minimize spr distance

################################################################################

INSTALLATION

SPR Supertrees is a command-line program written in C++. To use it, simply
compile spr_supertree.cpp and execute the resulting program. On systems
with the g++ compiler and make program, the included make file will
compile spr_supertrees; simply run `make'.

SPR Supertrees can also use multiple cores on SMP machines using OpenMP.
Compile with the -fopenmp flag or run `make omp'.

################################################################################

INPUT

SPR Supertrees requires a list of Newick format binary rooted trees with arbitrary labels
as input.  A sample Newick tree is shown below:

((1,2),(3,4),(5,6));

SPR Supertrees can also compare a rooted reference tree to an unrooted test
tree. Use the -unrooted or -unrooted_min_approx options and
input a list of rooted or unrooted Newick format binary trees.
rSPR will find the best rooting of each input tree with respect to the
current supertree using the -unrooted option and guess the best 
rooting based on the approximation algorithm with the
-unrooted_min_approx option.

################################################################################

OUTPUT

rspr writes to standard output.

A sample command line and output are shown below:

/////////////////////

$ ./spr_supertree -i 1 < test_trees/trees2.txt 
NUM_ITERATIONS=1
skipped 0 lines with no opening bracket 
skipped 0 multifurcating or invalid trees
skipped 0 trees with less than 4 leaves
2 gene trees remaining

Initial Supertree:  ((15,14),(13,12))
Adding leaf 11  (5/16)
(((16,15),(14,13)),12)
Adding leaf 10  (6/16)
(((16,15),(14,13)),(12,11))
Adding leaf 9   (7/16)
(((16,15),(14,13)),((12,11),10))
Adding leaf 8   (8/16)
((((16,15),(14,13)),9),((12,11),10))
Adding leaf 7   (9/16)
(((((16,15),(14,13)),9),((12,11),10)),8)
Adding leaf 6   (10/16)
(((((16,15),(14,13)),9),((12,11),10)),(8,7))
Adding leaf 5   (11/16)
(((((16,15),(14,13)),9),((12,11),10)),((8,7),6))
Adding leaf 4   (12/16)
(((((16,15),(14,13)),9),((12,11),10)),((8,7),(6,5)))
Adding leaf 3   (13/16)
(((((16,15),(14,13)),9),((12,11),10)),((8,7),((6,4),5)))
Adding leaf 2   (14/16)
(((((16,15),(14,13)),9),((12,11),10)),((8,7),((6,(4,3)),5)))
Adding leaf 1   (15/16)
(((((16,15),(14,13)),9),((12,11),10)),((8,7),((6,(4,3)),(5,2))))
Adding leaf 0   (16/16)
(((((16,15),(14,13)),9),((12,11),10)),((8,7),((6,(4,3)),((5,2),1))))

Initial Supertree:
(((((16,15),(14,13)),9),((12,11),10)),((8,7),((6,(4,3)),((5,2),1))))
Total Distance: 5
Current Supertree:
(((((16,15),(14,13)),9),((12,11),10)),((6,(8,7)),((4,3),((5,2),1))))
Total Distance: 4
Final Supertree:
(((((16,15),(14,13)),9),((12,11),10)),((6,(8,7)),((4,3),((5,2),1))))
Final Distance: 4

/////////////////////

The first set of lines indicate the options chosen, the number of invalid
trees and the number of valid trees.
The program then builds a supertree greedily by placing the most
frequent taxa first. Finally, the program applies 25 iterations of
global SPR rearrangements (or a user-specified number using the -i option
as shown here ) and outputs the best tree and distance found at the end of
each iteration.


################################################################################

EFFICIENCY

The 3-approximation algorithm runs in O(n) time, where n is the number of
leaves in the trees.

the exact algorithms run in  O(2.42^k n) time, where $k$ is the computed
SPR distance.

When using the -unrooted option, the exact algorithms run in O(2.42^k n^2) time.

When using the -max x option, the exact algorithms will run up to
a distance of x and then the approximation is used. This provides
a running time of O(n + 2.42^x n) or O(n + 2.42^x n^2) for rooted
trees and allows for a trade-off between space and efficiency.

Since there are O(n^2) possible SPR rearrangements, the total running
time is O(i * n^2 * X), where i is the number of iterations and X
is the running time of the chosen SPR computation method.
NOTE: The exact algorithms are exponential algorithms that exactly solve an NP-hard problem.
Thus the algorithms may not finish in a reasonable amount of time for large
rSPR distances (>50) .

################################################################################

REFERENCES

For more information on the algorithms see:

Whidden, C., Zeh, N., Beiko, R. G. Subtree Prune-and-Regraft Supertrees.

Whidden, C., Beiko, R. G., Zeh, N. Rooted Agreement Forests: Theory and
Experiments (Extended Abstract). Accepted to SEA 2010.

Whidden, C., Zeh, N. A Unifying View on Approximation and FPT of
Agreement Forests. In: WABI 2009. LNCS, vol. 5724, pp. 390.401.
Springer-Verlag (2009). Available at
http://www.springerlink.com/content/n56q2846v645p655/

Whidden, C. A Unifying View on Approximation and FPT of Agreement Forests.
Masters Thesis. Dalhousie University, Canada. 2009. Available at
www.cs.dal.ca/~whidden

################################################################################

CITING SPR Supertrees

If you use SPR Supertrees in your research, please cite:
Whidden, C., Zeh, N., Beiko R. G. Subtree Prune-and-Regraft Supertrees.

################################################################################
