################################################################################
spr_supertree

################################################################################

Usage: spr_supertree [OPTIONS]
       spr_supertree-omp [OPTIONS]
Calculate supertrees that minimize the SPR distance from the input
trees. By default calculates a rooted SPR supertree from a list
of rooted binary trees from STDIN in newick format. An initial
tree is built by greedily adding taxa in decreasing order of
ocurrence. The tree is then improved by SPR rearrangements.
Additional options allow for unrooted and/or multifurcating input trees.

Copyright 2013-2014 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/SPR_Supertrees
March 3, 2014
Version 1.2.1

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
from the supertree to the input trees. By default -bb is used.

-fpt        Calculate the exact rSPR distance with an FPT algorithm

-bb         Calculate the exact rSPR distance with a branch-and-bound
            FPT algorithm. This is the default option.

-approx     Calculate just a linear -time 3-approximation of the rSPR distance

-max k      Calculate the exact rSPR distance if it is k or less and
            otherwise use the 3-approximation

-split_approx
-split_approx x  Calculate the exact rSPR distance if it is k or less and
                 otherwise use the exponential-time approximation

*******************************************************************************
OPTIMIZATIONS
*******************************************************************************

These options control the use of optimized branching. All optimizations are
enabled by default. Specifying any subset of -cob, -cab, and -sc will use
just that subset of optimizations. See the README for more information.

-allopt   Use -cob -cab -sc and a new set of optimizations. This is the default
          option

-noopt    Use 3-way branching for all FPT algorithms

-cob      Use "cut one b" improved branching

-cab      Use "cut all b" improved branching

-sc       Use "separate components" improved branching

-bipartition_cluster x  Do not consider supertree rearrangements that violate
                        biparitions supported by x% of gene trees containing
                        at least two members from each side of the bipartition.
                        Enabled by default with x=0.5

*******************************************************************************
MULTIFURCATING COMPARISON OPTIONS
*******************************************************************************

-allow_multi   Allow multifurcating gene trees

-lgt_multi   Allow multifurcating input tree for LGT analysis

-support x     Collapse bipartitions with less than x support

*******************************************************************************
MULTIFURCATING NONBINARY NODE MOVE OPTIONS
*******************************************************************************

-lgt_move_parent    Moves the parent node of the nonbinary nodes to the destination
                    (Note: More taxa are transferred)

-lgt_move_individual_node    Moves each individual child node of the nonbinary node
                             to the individual child of nonbinary destination node
                             (Note: Overcounts the number of transfers)

-lgt_maintain_list    Moves the parent node of the nonbinary nodes to the destination
                      while maintaining the list of children of both source and
                      destination nonbinary node

*******************************************************************************
UNROOTED COMPARISON OPTIONS
*******************************************************************************

-unrooted   Compare the supertree to each rooting of the input trees.
            Use the best found distance

-unrooted_min_approx    Compare the supertree to each rooting of the
                        input trees.
                        Run the exact algorithm on the rooting with the
                        minimum approximate rspr distance

-simple_unrooted        Root the gene trees at each iteration using
                        a bipartition balanced accuracy measure
                        (fast but potentially less accurate)
                        Reports an unrooted SPR distance comparison
                        at the end of each iteration for comparable
                        iteration scores

-simple_unrooted x      Root the gene trees at the first x iterations

-simple_unrooted_fast   The same as -simple_unrooted but does not use
                        an unrooted comparison at the end of each
                        iteration

-outgroup FILE          Root the gene trees with the outgroup taxa
                        listed in FILE, one per line. Trees with a
                        polyphyletic outgroup are considered invalid.

-reroot                 Reroot the super tree at each iteration using
                        a bipartition balanced accuracy measure

-rspr_reroot            Root trees using the SPR distance instead
                        of the bipartition balanced accuracy



*******************************************************************************
SEARCH STRATEGY OPTIONS
*******************************************************************************

-i x    Run for x iterations of the global rearrangement search

-r x    Only consider transfers of length x in the global rearrangement
        search. Default is infinite (All SPRs). For NNI search use
        -r 1

-include_only <file>  Build the supertree only from taxa included in
                      <file>, one per line

-initial_tree <file>  Begin the search with the tree in <file>

-num_leaves x         Build the supertree from the x taxa that are found
                      in the largest number of trees

-random_insert_order  Insert taxa in random order when building the
                      greedy addition tree. The default order is
                      descending occurence

-rf_ties              Break SPR distance ties with the RF distance

*******************************************************************************
LGT ANALYSIS
*******************************************************************************

-lgt_analysis          Conduct an LGT analysis with the initial user-specified
                       or greedy addition tree

-lgt_evaluate          Print inferred transfers for each tree with the initial
                       user-specified or greedy addition tree

-lgt_csv               Output the LGT analysis seperated by commas rather than
                       spaces.

-lgt_groups FILE       Specify a set of groups (e.g. genus or class) to analyze
                       with -lgt_analysis. The group FILE contains a set of
                       groups consisting of a group name on one line, group
                       members one per line, and a blank line to seperate each
                       group.
                       
*******************************************************************************
OTHER OPTIONS
*******************************************************************************
-time                  Print iteration and total CPU time used at each
                       iteration

-cc                    Calculate a potentially better approximation with a
                       quadratic time algorithm

-valid_trees           Output the set of trees that appear valid
-valid_trees_rooted    Output the set of trees that appear valid after applying
                       any rooting options.

-multi_trees           Output the set of multifurcating or invalid trees

################################################################################

CONTACT INFORMATION

Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/SPR_Supertrees

################################################################################

FILES

ClusterForest.h   Cluster Decomposition
ClusterInstance.h Cluster Decomposition
Forest.h          Forest data structure
gen_rooted_trees.pl Generate all rootings of an unrooted binary tree
gpl.txt           The GPL license
LCA.h             Compute LCAs of tree leaves
lgt.h             LGT Analysis
Makefile          Makefile
Node.h            Node data structure
README.txt        This README
rspr.h            Calculate rSPR distances between pairs of trees
SiblingPair.h     Sibling pair data structure
spr_supertree.cpp Main file
spr_supertree     Compute supertrees that minimize spr distance
UndoMachine.h     Structure to record and undo tree alterations

################################################################################

INSTALLATION

SPR Supertrees is a command-line program written in C++. To use it, simply
compile spr_supertree.cpp and execute the resulting program. On systems
with the g++ compiler and make program, the included make file will
compile spr_supertree; simply run `make'.

SPR Supertrees can also use multiple cores on SMP machines through OpenMP.
Compile with the -fopenmp flag or run `make omp'. The multicore executable
will be called spr_supertree-omp

################################################################################

INPUT

SPR Supertrees requires a list of Newick format trees with arbitrary labels
as input.  A sample Newick tree is shown below:

((1,2),(3,4),(5,6));

By default the trees must be rooted and binary.
If you wish to allow multifurcating input trees use the -allow_multi
option. Bipartitions with less than x% support can be collapsed with
-support x.

SPR Supertrees can also construct a rooted tree from unrooted gene
trees. Use the -unrooted, -unrooted_min_approx, -simple_unrooted, or
-simple_unrooted_fast options rSPR will find the best rooting of each input
tree with respect to the current supertree using the -unrooted option, guess
the best rooting based on the approximation algorithm with the
-unrooted_min_approx option, and guess the best rooting based on
a bipartition balanced accuracy measure with the -simple_unrooted or
-simple_unrooted_fast options. These are much faster but may be less accurate.

The -outgroup <FILE> option roots gene trees based on a list of outgroup taxa.
This option ignores gene trees with a polyphyletic outgroup or no outgroup
members. To root these trees, one can construct a supertree from just the trees
where the outgroup is monophyletic and then root the remainder of the trees
with the -simple_unrooted 1 option.

With the -lgt_analysis option, the program conducts
an LGT analysis of an initial or greedy addition supertree. The gene trees
should be rooted, either as input or using -simple_unrooted_fast. This analysis
considers a single minimal reconciliation scenario between the supertree and
each gene tree. The output is a series of matrices (comma-seperated with the
-lgt_csv option) showing the number of inferred SPR moves, transfers, and
transfers ignoring direction between groups of taxa or to "mixed" portions of
the tree. The -lgt_groups <FILE> option is required and specifys taxonomic
groups or individual taxa.

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
each iteration. To build larger trees the -r x option will limit the
SPR rearrangements to transfers of length at most x. For example,
-r 1 uses only NNI rearrangements.


################################################################################

EFFICIENCY

The 3-approximation algorithm runs in O(n) time, where n is the number of
leaves in the trees.

the exact algorithms run in  O(2.42^k n) time, where $k$ is the computed
SPR distance. Using a set of new optimizations we conjecture that
the running time has been improved to O(2^k n) time.

When using the -unrooted option, the exact algorithms run in O(2.42^k
n^2) time. (conjectured O(2^k n^2)). The -simple_unrooted option
has the same worst case performance as the regular exact algorithms.

When using the -max x option, the exact algorithms will run up to
a distance of x and then the approximation is used. This provides
a running time of O(n + 2^x n) or O(n + 2^x n^2) for rooted
trees and allows for a trade-off between space and efficiency.
The -split_approx x option works similarly but is both much more
accurate and slower. -split_approx is recommended over -max.

Since there are O(n^2) possible SPR rearrangements, the total running
time is O(i * n^2 * X), where i is the number of iterations and X
is the running time of the chosen SPR computation method.
NOTE: The exact algorithms are exponential algorithms that exactly solve an
NP-hard problem.  Thus the algorithms may not finish in a reasonable amount
of time for very large rSPR distances without the -split_approx or -max
options. For very large supertrees, it may also be necessary to
limit the scope of the search with the -r option.

The -bipartition_cluster x option ignores SPR rearrangments that violate any
bipartition that agrees with x% of the gene trees that contain at least two
taxa from each side of the bipartition. This is enabled by default with x=0.5
and grealy accelerates tree searches at the expense of some searching power.
This option can be disabled with -bipartition_cluster 1, requiring total
agreement.

################################################################################

REFERENCES

For more information on the algorithms see:

Whidden, C., Zeh, N., Beiko, R.G.  Fixed-Parameter and Approximation
Algorithms for Maximum Agreement Forests of Multifurcating Trees.
(In Preparation). 2013. Preprint available at
http://arxiv.org/abs/1305.0512

Whidden, C., Beiko, R.G., Zeh, N. Fixed-Parameter Algorithms for Maximum
Agreement Forests. SIAM Journal on Computing 42.4 (2013), pp. 1431-1466.
Available at http://epubs.siam.org/doi/abs/10.1137/110845045

Whidden, C., Zeh, N., Beiko, R.G.  Supertrees based on the subtree
prune-and-regraft distance. Syst. Biol. 63 (4): 566-581. 2014.
doi:10.1093/sysbio/syu023.

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

CITING SPR Supertrees

If you use SPR Supertrees in your research, please cite:
Whidden, C., Zeh, N., Beiko, R.G.  Supertrees based on the subtree
prune-and-regraft distance.  Syst. Biol. 63 (4): 566-581. 2014.
doi:10.1093/sysbio/syu023.

################################################################################
