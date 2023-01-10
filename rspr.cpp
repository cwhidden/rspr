/*******************************************************************************
rspr.cpp

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
-min_length x      Collapse bipartitions with branch lengths less than or
                equal to x
-multifurcating Calculate the exact rSPR distance with a branch-and-bound FPT 
                algorithm. Both trees may be multifurcating
-multi_4_branch Calculate the exact rSPR distance with a branch-and-bound FPT
                algorithm. Both trees may be multifurcating. Use 4-way branching

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

*******************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <time.h>
#include <list>
#include "rspr.h"

#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "UndoMachine.h"
#include "lgt.h"

using namespace std;


// options to pick default
bool DEFAULT_ALGORITHM=true;
bool DEFAULT_OPTIMIZATIONS=true;


bool MULTI_APPROX = false;
bool FPT = false;
bool RF = false;
bool QUIET = false;
bool UNROOTED = false;
bool ALL_UNROOTED = false;
bool SIMPLE_UNROOTED = false;
bool SIMPLE_UNROOTED_RSPR = false;
bool LCA_TEST = false;
bool CLUSTER_TEST = false;
bool TOTAL = false;
bool PAIRWISE = false;
bool PAIRWISE_SYMMETRIC = true;
int PAIRWISE_START = 0;
int PAIRWISE_END = INT_MAX;
int PAIRWISE_COL_START = 0;
int PAIRWISE_COL_END = INT_MAX;
bool PAIRWISE_MAX = false;
int PAIRWISE_MAX_SPR = INT_MAX;
bool APPROX = false;
bool LOWER_BOUND = false;
bool REDUCE_ONLY = false;
bool PRINT_ROOTED_TREES = false;
bool SHOW_MOVES = false;
bool SEQUENCE = false;
bool DEBUG_REVERSE = false;
bool RANDOM_SPR = false;
int RANDOM_SPR_COUNT = 0;
int MULTI_TEST = 0;

string USAGE =
"rspr, version 1.3.1\n"
"\n"
"usage: rspr [OPTIONS]\n"
"Calculate approximate and exact Subtree Prune and Regraft (rSPR)\n"
"distances and the associated maximum agreement forests (MAFs) between pairs\n"
"of rooted binary trees from STDIN in newick format.\n"
"Supports arbitrary labels. See the README for more information.\n"
"\n"
"Copyright 2009-2021 Chris Whidden\n"
"whidden@cs.dal.ca\n"
"github.com/cwhidden/rspr\n"
"February 12, 2021\n"
"Version 1.3.1\n"
"\n"
"This program comes with ABSOLUTELY NO WARRANTY.\n"
"This is free software, and you are welcome to redistribute it\n"
"under certain conditions; See the README for details.\n"
"\n"
"*******************************************************************************\n"
"ALGORITHM\n"
"*******************************************************************************\n"
"\n"
"These options control what algorithm is used\n"
"\n"
"-fpt        Calculate the exact rSPR distance with an FPT algorithm\n"
"\n"
"-bb         Calculate the exact rSPR distance with a branch-and-bound\n"
"            FPT algorithm. This is enabled by default.\n"
"\n"
"-approx     Calculate just a linear-time 3-approximation of the rSPR distance\n"
"\n"
"-split_approx\n"
"-split_approx x  Calculate the exact rSPR distance if it is k or less and\n"
"                 otherwise use the exponential-time approximation\n"
"\n"
"-cluster_test   Use the cluster reduction to speed up the exact algorithm.\n"
"                This is enabled by default.\n"
"\n"
"-total          Find the total SPR distance from the first input tree to\n"
"                the rest of the list of trees. Uses the other algorithm\n"
"                options as specified (including unrooted options).\n"
"\n"
"*******************************************************************************\n"
"OPTIMIZATIONS\n"
"*******************************************************************************\n"
"\n"
"These options control the use of optimized branching. All optimizations are\n"
"enabled by default. Specifying any subset of -cob, -cab, and -sc will use\n"
"just that subset of optimizations.\n"
"\n"
"-allopt    Use -cob -cab -sc and a new set of improvements. This is the\n"
"           default\n"
"option\n"
"\n"
"-noopt     Use 3-way branching for all FPT algorithms\n"
"\n"
"-cob       Use \"cut one b\" improved branching\n"
"\n"
"-cab       Use \"cut all b\" improved branching\n"
"\n"
"-sc        Use \"separate components\" improved branching\n"
"\n"
"*******************************************************************************\n"
"MULTIFURCATING COMPARISON OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-support x     Collapse bipartitions with less than x support\n"
"-min_length x      Collapse bipartitions with branch lengths less than or\n"
"                equal to x\n"
"-multifurcating Calculate the exact rSPR distance with a branch-and-bound FPT\n"
"                algorithm. Both trees may be multifurcating\n"
"-multi_4_branch Calculate the exact rSPR distance with a branch-and-bound FPT\n"
"                algorithm. Both trees may be multifurcating. Use 4-way branching\n"
"\n"
"*******************************************************************************\n"
"UNROOTED COMPARISON OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-unrooted   Compare the first input tree to each other input tree.\n"
"            Output the best found distance and agreement forest.\n"
"            This option can be used with gen_rooted_trees.pl to provide\n"
"            the rootings.\n"
"            Note that this option is a bit unintuitive to maintain\n"
"            compatibility with previous versions of rSPR.\n"
"            If -total or -pairwise analysis is used then there is no need\n"
"            to specify rootings.\n"
"\n"
"-unrooted_min_approx    Compare the first input tree to each other input tree.\n"
"                        Run the exact algorithms on the pair with the\n"
"                        minimum approximate rspr distance\n"
"\n"
"-simple_unrooted        Root the gene trees using\n"
"                        a bipartition balanced accuracy measure\n"
"                        (fast but potentially less accurate). Only\n"
"                        used with -total.\n"
"\n"
"*******************************************************************************\n"
"PAIRWISE COMPARISON OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-pairwise\n"
"-pairwise a b\n"
"-pairwise a b c d        Compare each input tree to each other tree and output\n"
"                         the resulting SPR distance matrix. If -unrooted is\n"
"                         enabled this will compute the \"best rooting\" SPR\n"
"                         distance by testing each rooting of the trees. The\n"
"                         optional arguments a b c d compute only rows a-b and/or\n"
"                         columns c-d of the matrix.\n"
"\n"
"-no-symmetric-pairwise   By default, -pairwise will ignore the symmetric lower\n"
"                         left triangle of the matrix. With this option the\n"
"                         lower triangle is filled in.\n"
"\n"
"-pairwise_max x          Use with -pairwise to only compute distances at most x.\n"
"                         Larger values are output as -1. Very efficient for\n"
"                         small distances (e.g. 1-10).\n"
"\n"
"*******************************************************************************\n"
"OTHER OPTIONS\n"
"*******************************************************************************\n"
"-cc         Calculate a potentially better approximation with a quadratic time\n"
"            algorithm\n"
"\n"
"-q          Quiet; Do not output the input trees or approximation\n"
"*******************************************************************************\n"
"\n"
"Example:\n"
"$ ./rspr < test_trees/trees2.txt\n"
"T1: ((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))))\n"
"T2: (((7,8),((1,(2,(14,5))),(3,4))),(((11,(6,12)),10),((13,(15,16)),9)))\n"
"\n"
"T1: ((((0,1),(2,3)),((4,5),(6,7))),(((8,9),(10,11)),((12,13),(14,15))))\n"
"T2: (((6,7),((0,(1,(13,4))),(2,3))),(((10,(5,11)),9),((12,(14,15)),8)))\n"
"approx F1: ((0,(2,3)),((9,(10,11)),(12,(14,15)))) 13 5 8 4 (6,7) 1\n"
"approx F2: ((0,(2,3)),(((10,11),9),(12,(14,15)))) 13 5 8 4 1 (6,7)\n"
"approx drSPR=6\n"
"\n"
"C1_1: ((((0,1),(2,3)),((4,5),(6,7))),(((8,9),(10,11)),((12,13),(14,15))))\n"
"C1_2: (((6,7),((0,(1,(13,4))),(2,3))),(((10,(5,11)),9),((12,(14,15)),8)))\n"
"cluster approx drSPR=4\n"
"\n"
"4\n"
"F1_1: ((((0,1),(2,3)),(6,7)),((9,(10,11)),(12,(14,15)))) 5 4 13 8\n"
"F1_2: (((6,7),((0,1),(2,3))),(((10,11),9),(12,(14,15)))) 13 5 4 8\n"
"cluster exact drSPR=4\n"
"\n"
"F1: ((((1,2),(3,4)),(7,8)),((10,(11,12)),(13,(15,16)))) 6 5 14 9\n"
"F2: (((7,8),((1,2),(3,4))),(((11,12),10),(13,(15,16)))) 14 6 5 9\n"
"total exact drSPR=4\n"

"*******************************************************************************/\n";

int main(int argc, char *argv[]) {
	int max_args = argc-1;
	while (argc > 1) {
		char *arg = argv[--argc];
		if (strcmp(arg, "-rf") == 0) {
			RF = true;
		}
		if (strcmp(arg, "-fpt") == 0) {
			FPT = true;
			DEFAULT_ALGORITHM=false;
		}
		else if (strcmp(arg, "-bb") == 0) {
			BB = true;
			DEFAULT_ALGORITHM=false;
		}
		else if (strcmp(arg, "-approx") == 0) {
			DEFAULT_ALGORITHM=false;
			APPROX=true;
		}
		else if (strcmp(arg, "-lower_bound") == 0) {
			DEFAULT_ALGORITHM=false;
			APPROX=true;
			LOWER_BOUND=true;
		}
		else if (strcmp(arg, "-fast_approx") == 0) {
			APPROX_CUT_ONE_B = true;
			APPROX_CUT_TWO_B = true;
//			APPROX_CUT_TWO_B_ROOT = true;
			APPROX_REVERSE_CUT_ONE_B = true;
//			APPROX_EDGE_PROTECTION = true;
		}
		else if (strcmp(arg, "-multi_approx") == 0) {
		        DEFAULT_ALGORITHM = false;
			BB = true;
			MULTI_APPROX = true;
			DEFAULT_OPTIMIZATIONS = false;
		}
		else if (strcmp(arg, "-multifurcating") == 0) {
		  //DEFAULT_ALGORITHM = false;
		  //	BB = true;
			MULTIFURCATING = true;
			//DEFAULT_OPTIMIZATIONS = false;
		}
		else if (strcmp(arg, "-multi_4_branch") == 0) {
		  //DEFAULT_ALGORITHM = false;
		  //	BB = true;
			MULTIFURCATING = true;
			//DEFAULT_OPTIMIZATIONS = false;
			MULT_4_BRANCH = true;
		}
		else if (strcmp(arg, "-debug_reverse_trees") == 0) {
		        DEBUG_REVERSE = true;
		}

		else if (strcmp(arg, "-a_cob") == 0) {
			APPROX_CUT_ONE_B = true;
		}
		else if (strcmp(arg, "-a_c2b") == 0) {
			APPROX_CUT_TWO_B = true;
		}
		else if (strcmp(arg, "-a_c2br") == 0) {
			APPROX_CUT_TWO_B_ROOT = true;
		}
		else if (strcmp(arg, "-a_rcob") == 0) {
			APPROX_REVERSE_CUT_ONE_B = true;
		}
		else if (strcmp(arg, "-a_rcob2") == 0) {
			APPROX_REVERSE_CUT_ONE_B_2 = true;
		}
		else if (strcmp(arg, "-a_protection") == 0) {
			APPROX_EDGE_PROTECTION = true;
		}
		else if (strcmp(arg, "-q") == 0)
			QUIET = true;
		else if (strcmp(arg, "-cc") == 0)
			APPROX_CHECK_COMPONENT = true;
		else if (strcmp(arg, "-unrooted") == 0)
			UNROOTED = true;
		else if (strcmp(arg, "-all_unrooted") == 0) {
			ALL_UNROOTED = true;
			UNROOTED = true;
		}
		else if (strcmp(arg, "-simple_unrooted") == 0)
			SIMPLE_UNROOTED = true;
		else if (strcmp(arg, "-simple_unrooted_rspr") == 0) {
			SIMPLE_UNROOTED = true;
			SIMPLE_UNROOTED_RSPR = true;
		}
		else if (strcmp(arg, "-simple_unrooted_leaf") == 0) {
			SIMPLE_UNROOTED = true;
			SIMPLE_UNROOTED_LEAF = 1;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					SIMPLE_UNROOTED_LEAF = atoi(arg2);
			}
		}
		else if (strcmp(arg, "-unrooted_min_approx") == 0)
			UNROOTED_MIN_APPROX = true;
		else if (strcmp(arg, "-noopt") == 0) {
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-cut_one_b") == 0 ||
				strcmp(arg, "-cob") == 0) {
			CUT_ONE_B = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-reverse_cut_one_b") == 0 ||
				strcmp(arg, "-rcob") == 0) {
			REVERSE_CUT_ONE_B = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-reverse_cut_one_b_2") == 0 ||
				strcmp(arg, "-rcob2") == 0) {
			REVERSE_CUT_ONE_B_2 = true;
		}
		else if (strcmp(arg, "-reverse_cut_one_b_3") == 0 ||
				strcmp(arg, "-rcob3") == 0) {
			REVERSE_CUT_ONE_B_3 = true;
		}
		else if (strcmp(arg, "-cut_two_b") == 0 ||
				strcmp(arg, "-c2b") == 0) {
			CUT_TWO_B = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-cut_two_b_root") == 0 ||
				strcmp(arg, "-c2br") == 0) {
			CUT_TWO_B_ROOT = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-cut_all_b") == 0 ||
				strcmp(arg, "-cab") == 0) {
			CUT_ALL_B = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-cut_ac_separate_components") == 0 ||
				strcmp(arg, "-sc") == 0) {
			CUT_AC_SEPARATE_COMPONENTS = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-cut_one_ab") == 0) {
			CUT_ONE_AB = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-h") == 0) {
			cout << USAGE;
			return 0;
		}
		else if (strcmp(arg, "-lca") == 0) {
			LCA_TEST = true;
		}
		else if (strcmp(arg, "-find_rate") == 0) {
			FIND_RATE = true;
		}
		else if (strcmp(arg, "-reduce") == 0) {
			REDUCE_ONLY = true;
		}
		else if (strcmp(arg, "-print_rooted_trees") == 0) {
			PRINT_ROOTED_TREES = true;
		}

/*		else if (strcmp(arg, "-cluster") == 0) {
			CLUSTER_REDUCTION = true;
			PREFER_RHO = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					MAX_CLUSTERS = atoi(arg2);
				cout << "MAX_CLUSTERS=" << MAX_CLUSTERS << endl;
			}
		}
*/

		else if (strcmp(arg, "-cluster_test") == 0) {
			CLUSTER_TEST = true;
			PREFER_RHO = true;
		}
		else if (strcmp(arg, "-prefer_rho") == 0) {
			PREFER_RHO = true;
		}
/*
		else if (strcmp(arg, "-memoize") == 0) {
			MEMOIZE = true;
		}
*/
		else if (strcmp(arg, "-all_mafs") == 0) {
			ALL_MAFS= true;
		}
		else if (strcmp(arg, "-total") == 0) {
			TOTAL= true;
			//PREFER_RHO = true;
		}
		else if (strcmp(arg, "-pairwise") == 0) {
			PAIRWISE= true;
			QUIET=true;
			//PREFER_RHO = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-') {
					PAIRWISE_START = atoi(arg2);
					if (max_args > argc+1) {
						char *arg2 = argv[argc+2];
						if (arg2[0] != '-') {
							PAIRWISE_END = atoi(arg2);
							if (max_args > argc+2) {
								char *arg2 = argv[argc+3];
								if (arg2[0] != '-') {
									PAIRWISE_COL_START = atoi(arg2);
									if (max_args > argc+3) {
										char *arg2 = argv[argc+4];
										if (arg2[0] != '-') {
											PAIRWISE_COL_END = atoi(arg2);
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else if (strcmp(arg, "-no-symmetric-pairwise") == 0) {
			PAIRWISE_SYMMETRIC=false;
		}
		else if (strcmp(arg, "-pairwise_max") == 0) {
			PAIRWISE=true;
			PAIRWISE_MAX=true;
			QUIET=true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-') {
					PAIRWISE_MAX_SPR = atoi(arg2);
				}
			}
		}
		else if (strcmp(arg, "-cluster_tune") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-') {
					CLUSTER_TUNE = atoi(arg2);
				}
			}
		}
		else if (strcmp(arg, "-v") == 0) {
			VERBOSE=true;
		}
		else if (strcmp(arg, "-max") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-') {
					MAX_SPR = atoi(arg2);
					CLUSTER_MAX_SPR = MAX_SPR;
				}
				if (!QUIET)
					cout << "MAX_SPR=" << MAX_SPR << endl;
			}
		}
		else if (strcmp(arg, "-cmax") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					CLUSTER_MAX_SPR = atoi(arg2);
				cout << "CLUSTER_MAX_SPR=" << CLUSTER_MAX_SPR << endl;
			}
		}
		else if (strcmp(arg, "-multi_test") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					MULTI_TEST = atoi(arg2);
				cout << "MULTI_TEST=" << MULTI_TEST << endl;
			}
		}
		else if (strcmp(arg, "-protect_edges") == 0) {
			EDGE_PROTECTION = true;
			cout << "EDGE_PROTECTION=" << EDGE_PROTECTION << endl;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-protect_edges_two_b") == 0) {
			EDGE_PROTECTION_TWO_B = true;
//			cout << "EDGE_PROTECTION=" << EDGE_PROTECTION << endl;
//			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-check_merge_depth") == 0) {
			CHECK_MERGE_DEPTH = true;
			cout << "CHECK_MERGE_DEPTH=" << CHECK_MERGE_DEPTH << endl;
//			DEFAULT_OPTIMIZATIONS=false;
		}
//		else if (strcmp(arg, "-check_fewer_pairs") == 0) {
//			check_all_pairs = false;
//			DEFAULT_OPTIMIZATIONS=false;
//		}
		else if (strcmp(arg, "-allow_abort") == 0) {
			ABORT_AT_FIRST_SOLUTION = true;
//			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-preorder_sib_pairs") == 0) {
			PREORDER_SIBLING_PAIRS = true;
			NEAR_PREORDER_SIBLING_PAIRS = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-near_preorder_sib_pairs") == 0) {
			NEAR_PREORDER_SIBLING_PAIRS = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-leaf_reduction") == 0) {
			LEAF_REDUCTION = true;
		}
		else if (strcmp(arg, "-leaf_reduction2") == 0) {
			LEAF_REDUCTION2 = true;
		}
		else if (strcmp(arg, "-split_approx") == 0) {
			SPLIT_APPROX = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					SPLIT_APPROX_THRESHOLD = atoi(arg2);
			}
			cout << "SPLIT_APPROX_THRESHOLD=" << SPLIT_APPROX_THRESHOLD
					<< endl;
		}
		
		//nocheckin
		else if (strcmp(arg, "-random_spr") == 0) {
			RANDOM_SPR = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					RANDOM_SPR_COUNT = atoi(arg2);
			}
		}

		else if (strcmp(arg, "-support") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					REQUIRED_SUPPORT = atof(arg2);
				//cout << "REQUIRED_SUPPORT=" << REQUIRED_SUPPORT
				//<< endl;
			}
		}
		else if (strcmp(arg, "-min_length") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					MIN_LENGTH = atof(arg2);
				cout << "MIN_LENGTH=" << MIN_LENGTH
						<< endl;
			}
		}
		else if (strcmp(arg, "-prefer_nonbranching") == 0) {
			PREFER_NONBRANCHING = true;
		}
		else if (strcmp(arg, "-deepest") == 0) {
			DEEPEST_ORDER = true;
		}
		else if (strcmp(arg, "-deepest_protected") == 0) {
			DEEPEST_PROTECTED_ORDER = true;
			DEEPEST_ORDER = true;
		}
		else if (strcmp(arg, "-count_losses") == 0) {
			COUNT_LOSSES = true;
		}
		else if (strcmp(arg, "-cut_lost") == 0) {
			CUT_LOST = true;
		}
		else if (strcmp(arg, "-multi_cluster") == 0) {
			MULTI_CLUSTER = true;
		}
		else if (strcmp(arg, "-show_moves") == 0) {
			SHOW_MOVES = true;
		}
		else if (strcmp(arg, "-sequence") == 0) {
			SEQUENCE = true;
		}
		else if (strcmp(arg, "--help") == 0 || strcmp(arg, "-help") == 0) {
			cout << USAGE;
			return 0;
		}
			
	}
	if (DEFAULT_OPTIMIZATIONS) {
		CUT_ALL_B=true;
		CUT_ONE_B = true;
		REVERSE_CUT_ONE_B = true;
//		REVERSE_CUT_ONE_B_2 = true;
		REVERSE_CUT_ONE_B_3 = true;
		CUT_TWO_B = true;
//		CUT_TWO_B_ROOT = true;
		CUT_AC_SEPARATE_COMPONENTS = true;
		EDGE_PROTECTION = true;
		EDGE_PROTECTION_TWO_B = true;
//		CHECK_MERGE_DEPTH = true;
//		if (ALL_MAFS == false)
//			ABORT_AT_FIRST_SOLUTION = true;
//		PREORDER_SIBLING_PAIRS = true;
		NEAR_PREORDER_SIBLING_PAIRS = true;
		LEAF_REDUCTION = true;
		LEAF_REDUCTION2 = true;
		PREFER_NONBRANCHING = true;

		APPROX_CUT_ONE_B = true;
		APPROX_CUT_TWO_B = true;
//		APPROX_CUT_TWO_B_ROOT = true;
		APPROX_REVERSE_CUT_ONE_B = true;
/* BUGGY: we aren't guaranteed that the protected edges mean
	 anything because we may cut off the only things that can merge with
	 them. It might make sense to cut a protected edge because it should
	 have already merged by then.
*/

//		APPROX_EDGE_PROTECTION = true;
			DEEPEST_PROTECTED_ORDER = true;
			DEEPEST_ORDER = true;
		if (CLUSTER_TUNE == -1) {
		    if (MULTIFURCATING)
		        CLUSTER_TUNE = 30;
		    else
			CLUSTER_TUNE = 30;
		}
	}
	PREORDER_SIBLING_PAIRS = true;
	if (DEFAULT_ALGORITHM) {
		BB=true;
		CLUSTER_TEST = true;
		PREFER_RHO = true;
	}

	// arguments that imply quiet
	if (PAIRWISE) {
		QUIET=true;
	}

	// Label maps to allow string labels
	map<string, int> label_map= map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();

	// set random seed
	srand(unsigned(time(0)));

	// Normal operation
	if (!UNROOTED && !UNROOTED_MIN_APPROX && !TOTAL && !PAIRWISE && !SEQUENCE) {
		string T1_line = "";
		string T2_line = "";
		while (getline(cin, T1_line) && getline(cin, T2_line)) {
		  Node *T1;
		  Node *T2;
		        if (DEBUG_REVERSE) {
			        T1 = build_tree(T2_line);
				T2 = build_tree(T1_line);
		        }
		        else {
			        T1 = build_tree(T1_line);
				T2 = build_tree(T2_line);
		        }

			//nocheckin
			if (RANDOM_SPR) {
			  Forest F1 = Forest(T1);
			  Forest F2 = Forest(T2);
			  sync_twins(&F1, &F2);
			  randomize_tree_with_spr(&F1, &F2, RANDOM_SPR_COUNT);
			  F1.get_component(0)->print_subtree_hlpr();
			  cout << ";" << endl;
			  F2.get_component(0)->print_subtree_hlpr();
			  cout << ";";
			  continue;
			}
			
			if (MULTI_TEST > 0) {
				vector<Node *> interior = T2->find_interior();
				vector<Node *> remove = random_select(interior, MULTI_TEST);
				vector<Node *>::iterator i;
				for(i = remove.begin(); i != remove.end(); i++) {
					(*i)->contract_node();
				}
			}
			// TODO: should we sync here to prune out additional leaves?
			if (REDUCE_ONLY) {

				T1->preorder_number();
				T1->edge_preorder_interval();
				T2->preorder_number();
				T2->edge_preorder_interval();

			  
				Forest F1 = Forest(T1);
				Forest F2 = Forest(T2);
				
				F1.print_components();
				F2.print_components();

				sync_twins(&F1, &F2);
				
				if (MULTIFURCATING) {
				  reduction_leaf_mult(&F1, &F2);
				}
				else {
				  reduction_leaf(&F1, &F2);
				}
				
				F1.print_components();
				F2.print_components();
				continue;
			}
			if (!QUIET) {
				cout << "T1: ";
				T1->print_subtree();
				cout << "T2: ";
				T2->print_subtree();
				cout << endl;
			}
			if (LCA_TEST) {
				LCA lca_query = LCA(T1);
				cout << endl;
				lca_query.debug();
				cout << endl;
				vector<Node *> leaves = T1->find_leaves();
				for(vector<Node *>::iterator i = leaves.begin(); i != leaves.end(); i++) {
					for(vector<Node *>::iterator j = i; j != leaves.end(); j++) {
						if (j==i)
							continue;
						Node *lca = lca_query.get_lca(*i, *j);
						(*i)->print_subtree_hlpr();
						cout << "\t";
						(*j)->print_subtree_hlpr();
						cout << "\t";
						lca->print_subtree();
					}
				}
				T1->delete_tree();
				T2->delete_tree();
				return 0;
		}


			T1->labels_to_numbers(&label_map, &reverse_label_map);
			T2->labels_to_numbers(&label_map, &reverse_label_map);


			if (SHOW_MOVES) {
				//cout << T1->str_subtree() << endl;
				//cout << T2->str_subtree() << endl;
				show_moves(T1, T2, &label_map, &reverse_label_map);
				T1->delete_tree();
				T2->delete_tree();
				return 0;
			}

	
//			ClusterForest F1 = ClusterForest(T1);
//			ClusterForest F2 = ClusterForest(T2);
//			ClusterForest F3 = ClusterForest(T1);
//			ClusterForest F4 = ClusterForest(T2);			
			
			ClusterForest F1 = ClusterForest(T1);
			ClusterForest F2 = ClusterForest(T2);
			Forest F3 = Forest(T1);
			Forest F4 = Forest(T2);

			if (MULTIFURCATING) {
			  T1->preorder_number();
			  T2->preorder_number();
			}
			
			
			if (RF) {
				int rf_d = rf_distance(T1, T2);
				cout << "RF Distance=" << rf_d << endl;
				T1->delete_tree();
				T2->delete_tree();
				continue;
			}
			if (CLUSTER_TEST) {
				T1->preorder_number();
				T1->edge_preorder_interval();
				T2->preorder_number();
				T2->edge_preorder_interval();
				
				
				int exact_k = rSPR_branch_and_bound_simple_clustering(T1,T2,true, &label_map, &reverse_label_map);
				//int exact_k = rSPR_branch_and_bound_simple_clustering(&F3,&F4,true, &label_map, &reverse_label_map);

				T1->delete_tree();
				T2->delete_tree();
				continue;
			}
	
			// APPROX ALGORITHM
			int approx_spr;
			int min_spr;
			if (MULTI_APPROX || MULTIFURCATING)
			{
			  F1.components[0]->preorder_number(0);
			  F2.components[0]->preorder_number(0);

			    approx_spr = rSPR_worse_3_mult_approx(&F1, &F2);
				min_spr = approx_spr / 3;				
			}
			else
			{
			    approx_spr = rSPR_worse_3_approx(&F1, &F2);
				min_spr = approx_spr / 3;
			}
			if (!(QUIET && (BB || FPT))) {
				F1.numbers_to_labels(&reverse_label_map);
				F2.numbers_to_labels(&reverse_label_map);
				cout << "approx F1: ";
				F1.print_components();
				cout << "approx F2: ";
				F2.print_components();
				// what the AF shows
				cout << "approx drSPR=" << F2.num_components()-1 << endl;
				if (LOWER_BOUND)
					cout << "lower bound drSPR=" << min_spr << endl;
				/* what we use to get the lower bound: 3 * the number of cutting rounds in
					 the approx algorithm
				*/
				//cout << "approx drSPR=" << approx_spr << endl;
				cout << "\n";
				if (MULTI_APPROX) {
				  T1->delete_tree();
				  T2->delete_tree();
				  continue;
				}
			}
	
			int k = min_spr;

			/*
			// FPT ALGORITHM
			int exact_spr = -1;
			int k = min_spr;
			if (FPT) {
				for(k = min_spr; k <= MAX_SPR; k++) {
					cout << k << " ";
					cout.flush();
					Forest F1 = Forest(F3);
					Forest F2 = Forest(F4);
					exact_spr = rSPR_FPT(&F1, &F2, k);
					if (exact_spr >= 0) {
						F1.numbers_to_labels(&reverse_label_map);
						F2.numbers_to_labels(&reverse_label_map);
						cout << endl;
						cout << "F1: ";
						F1.print_components();
						cout << "F2: ";
						F2.print_components();
						cout << "exact drSPR=" << exact_spr << endl;
						T1->delete_tree();
						T2->delete_tree();
						continue;
					}
				}
				if (exact_spr == -1)
					cout << "exact drSPR=?  " << "k=" << k << " too large" << endl;
				cout << "\n";
			}
			*/
		
			if (BB || FPT) {
				// BRANCH AND BOUND FPT ALGORITHM
				Forest F1 = Forest(F3);
				Forest F2 = Forest(F4);
				int exact_spr;
				if (MULTIFURCATING) {
				   F1.components[0]->preorder_number(0);
				   F2.components[0]->preorder_number(0);
			 
				  exact_spr = rSPR_branch_and_bound_mult_range(&F1, &F2, min_spr);
				}
				else {
				        exact_spr = rSPR_branch_and_bound(&F1, &F2);
				}
				if (exact_spr >= 0) {
					cout << "F1: ";
					F1.print_components();
					cout << "F2: ";
					F2.print_components();
					F1.numbers_to_labels(&reverse_label_map);
					F2.numbers_to_labels(&reverse_label_map);
					cout << endl;
					cout << "F1: ";
					F1.print_components();
					cout << "F2: ";
					F2.print_components();
					if (FPT)
						cout << "exact drSPR=" << exact_spr << endl;
					else
						cout << "exact BB drSPR=" << exact_spr << endl;
					T1->delete_tree();
					T2->delete_tree();
					continue;
				}
				if (exact_spr == -1) {
						if (FPT) {
							cout << "exact drSPR=?  " << "k=" << k << " too large"
								<< endl;
						}
						else {
							cout << "exact BB drSPR=?  " << "k=" << k << " too large"
								<< endl;
						}
				}
				cout << "\n";
			}
	
			// cleanup
			T1->delete_tree();
			T2->delete_tree();
		}
	}
	// Comparison between a rooted tree and all rootings of an unrooted tree
	else if (!TOTAL && !PAIRWISE && !SEQUENCE && (UNROOTED || UNROOTED_MIN_APPROX)) {
		string line = "";
		vector<Forest> trees = vector<Forest>();
		if (!getline(cin, line))
			return 0;
		Node *T1 = build_tree(line);
		if (!QUIET) {
			cout << "T1: ";
			T1->print_subtree();
		}
		T1->labels_to_numbers(&label_map, &reverse_label_map);
		Forest F1 = Forest(T1);
		while (getline(cin, line)) {
			Node *T2 = build_tree(line);
			if (!QUIET) {
				cout << "T2: ";
				T2->print_subtree();
			}
			T2->labels_to_numbers(&label_map, &reverse_label_map);
			trees.push_back(Forest(T2));
		}
		cout << endl;

		if (trees.size() == 0)
			return 0;


		// APPROX ALGORITHM
		int min_spr = INT_MAX;
		int min_i = 0;
		vector<int> approx_spr = vector<int>(trees.size());
		for (int i = 0; i < trees.size(); i++) {
			Forest F3 = Forest(F1);
			Forest F4 = Forest(trees[i]);
			approx_spr[i] = rSPR_worse_3_approx(&F3, &F4);
			if (approx_spr[i] < min_spr) {
				min_spr = approx_spr[i];
				min_i = i;
			}
			if (!(QUIET && (BB || FPT)) && !UNROOTED_MIN_APPROX) {
				F3.numbers_to_labels(&reverse_label_map);
				F4.numbers_to_labels(&reverse_label_map);
				cout << "F1: ";
				F3.print_components();
				cout << "F2: ";
				F4.print_components();
				cout << "approx drSPR=" << approx_spr[i] << endl;
				cout << "\n";
			}
		}
		// Choose a rooting with minimum approximate distance
		if (UNROOTED_MIN_APPROX) {
			Forest min_approx_forest = trees[min_i];
			trees.clear();
			trees.push_back(min_approx_forest);
				F1.numbers_to_labels(&reverse_label_map);
				min_approx_forest.numbers_to_labels(&reverse_label_map);
				cout << "F1: ";
				F1.print_components();
				cout << "F2: ";
				min_approx_forest.print_components();
				F1.labels_to_numbers(&label_map, &reverse_label_map);
				min_approx_forest.labels_to_numbers(&label_map, &reverse_label_map);
		}

		cout << "min approx drSPR=" << min_spr << endl;
		cout << "\n";

		min_spr /= 3;

		int k, exact_spr;
		if (FPT || BB) {
			// BRANCH AND BOUND FPT ALGORITHM
			for(k = min_spr; k <=MAX_SPR;  k++) {
				cout << k << " ";
				cout.flush();
				for (int i = 0; i < trees.size(); i++) {
					Forest F3 = Forest(F1);
					Forest F4 = Forest(trees[i]);
					exact_spr = rSPR_branch_and_bound(&F3, &F4, k);
					if (exact_spr >= 0) {
						//sync_twins(&F1, &trees[i]);
						F1.numbers_to_labels(&reverse_label_map);
						trees[i].numbers_to_labels(&reverse_label_map);
						F3.numbers_to_labels(&reverse_label_map);
						F4.numbers_to_labels(&reverse_label_map);
						cout << endl;
						cout << "T1: ";
						F1.print_components();
						cout << "T2: ";
						trees[i].print_components();
						cout << endl;
						cout << "F1: ";
						F3.print_components();
						cout << "F2: ";
						F4.print_components();
						cout << "exact BB drSPR=" << exact_spr << endl;
						F1.labels_to_numbers(&label_map, &reverse_label_map);
						trees[i].labels_to_numbers(&label_map, &reverse_label_map);
						break;
					}
				}
				if (exact_spr >= 0) {
					break;
				}
			}
			if (exact_spr == -1)
				cout << "exact BB drSPR=?  " << "k=" << k << " too large" << endl;
			cout << "\n";
		}

	}
	else if (TOTAL) {
		string line = "";
		vector<Node *> trees = vector<Node *>();
		vector<string> names = vector<string>();
		if (!getline(cin, line))
			return 0;
		if (UNROOTED || SIMPLE_UNROOTED)
			line = root(line);
		Node *T1 = build_tree(line);
		if (!QUIET) {
			cout << "T1: ";
			T1->print_subtree();
		}
		T1->labels_to_numbers(&label_map, &reverse_label_map);
		while (getline(cin, line)) {
			size_t loc = line.find_first_of("(");
			if (loc != string::npos) {
				string name = "";
				if (loc != 0) {
					name = line.substr(0,loc);
					line.erase(0,loc);
				}
				if (UNROOTED || SIMPLE_UNROOTED)
					line = root(line);
				Node *T2 = build_tree(line);
				if (!QUIET) {
					cout << "T2: ";
					T2->print_subtree();
				}
				T2->labels_to_numbers(&label_map, &reverse_label_map);
				if (UNROOTED)
					T2->preorder_number();
				names.push_back(name);
				trees.push_back(T2);
			}
		}
		if (!QUIET) {
			cout << endl;
		}

		if (SIMPLE_UNROOTED) {
			// reroot the gene trees based on the balanced accuracy of splits
			T1->preorder_number();
			int end = trees.size();
			#pragma omp parallel for
			for(int i = 0; i < end; i++) {
				trees[i]->preorder_number();
				Node *new_root;
				if (SIMPLE_UNROOTED_RSPR)
					new_root = find_best_root_rspr(T1, trees[i]);
				else
					new_root = find_best_root(T1, trees[i]);
				if (new_root != NULL)
					trees[i]->reroot(new_root);
					trees[i]->set_depth(0);
					trees[i]->fix_depths();
					trees[i]->preorder_number();
				if (PRINT_ROOTED_TREES) {
					trees[i]->numbers_to_labels(&reverse_label_map);
					cout << "T" <<  i+2 << ": " << trees[i]->str_subtree() << endl;
					trees[i]->labels_to_numbers(&label_map, &reverse_label_map);
				}
			}
		}
		vector<Node *> rootings;
		if (ALL_UNROOTED) {
			rootings = T1->find_descendants();
		}
		else {
			rootings = vector<Node *>();
			rootings.push_back(T1->lchild());
		}
		int best_distance = INT_MAX;
		for(int i = 0; i < rootings.size(); i++) {
			if (rootings[i] != T1)
				T1->reroot(rootings[i]);
			T1->set_depth(0);
			T1->fix_depths();
			T1->preorder_number();

			int distance;

			if (VERBOSE) {
				T1->numbers_to_labels(&reverse_label_map);
				cout << "T1: " <<  T1->str_subtree() << endl;
				T1->labels_to_numbers(&label_map, &reverse_label_map);
			}

			if (RF) {
				if (UNROOTED) {
					distance = rf_total_distance_unrooted(T1, trees);
				}
				else {
					distance = rf_total_distance(T1, trees);
				}
					cout << "total RF distance=" << distance << endl;
			}
			else if (APPROX) {
				if (UNROOTED) {
					distance = rSPR_total_approx_distance_unrooted(T1,trees);
				}
				else 
					distance = rSPR_total_approx_distance(T1,trees);
				cout << "total approx distance= " << distance << endl;
			}
			else {
				if (UNROOTED)
					distance = rSPR_total_distance_unrooted(T1,trees);
				else
					distance = rSPR_total_distance(T1,trees);
				cout << "total distance= " << distance << endl;
			}
			if (ALL_UNROOTED && distance < best_distance)
				best_distance = distance;

		}
		if (ALL_UNROOTED)
			cout << "best distance=" << best_distance << endl;
		T1->delete_tree();
		for(vector<Node *>::iterator T2 = trees.begin(); T2 != trees.end(); T2++)
			(*T2)->delete_tree();
	}
	else if (PAIRWISE) {
		string line = "";
		vector<Node *> trees = vector<Node *>();
		vector<string> names = vector<string>();
//		if (!getline(cin, line))
//			return 0;
//		Node *T1 = build_tree(line);
//		if (!QUIET) {
//			cout << "T1: ";
//			T1->print_subtree();
//		}
//		T1->labels_to_numbers(&label_map, &reverse_label_map);
		while (getline(cin, line)) {
			size_t loc = line.find_first_of("(");
			if (loc != string::npos) {
				string name = "";
				if (loc != 0) {
					name = line.substr(0,loc);
					line.erase(0,loc);
				}
				if (UNROOTED || SIMPLE_UNROOTED)
					line = root(line);
				Node *T2 = build_tree(line);
//				if (!QUIET) {
//					cout << "T2: ";
//					T2->print_subtree();
//				}
				T2->labels_to_numbers(&label_map, &reverse_label_map);
				if (UNROOTED)
					T2->preorder_number();
				names.push_back(name);
				trees.push_back(T2);
			}
		}
//		if (!QUIET) {
//			cout << endl;
//		}
		int start_i = PAIRWISE_START;
		if (start_i < 0) {
			start_i = 0;
		}
		int end_i = PAIRWISE_END;
		if (end_i > trees.size()) {
			end_i = trees.size();
		}
		int start_j = PAIRWISE_COL_START;
		if (start_j < 0) {
			start_j = 0;
		}
		int end_j = PAIRWISE_COL_END;
		if (end_j > trees.size()) {
			end_j = trees.size();
		}

		for(int i = start_i; i < end_i; i++) {
			int j = start_j;
			if (PAIRWISE_SYMMETRIC) {
				for(j = start_j; j < i && j < end_j; j++) {
					cout << ",";
				}
			}
			if (RF) {
				if (UNROOTED) {
					rf_pairwise_distance_unrooted(trees[i], trees, j, end_j);
				}
				else {
					rf_pairwise_distance(trees[i], trees, j, end_j);
				}
			}
			else {
				if (UNROOTED) {
					if (PAIRWISE_MAX) {
						rSPR_pairwise_distance_unrooted(trees[i], trees, PAIRWISE_MAX_SPR, j, end_j);
					}
					else {
						rSPR_pairwise_distance_unrooted(trees[i], trees, j, end_j, APPROX);
					}
				}
				else {
					if (PAIRWISE_MAX) {
						rSPR_pairwise_distance(trees[i], trees, PAIRWISE_MAX_SPR, j, end_j);
					}
					else {
						rSPR_pairwise_distance(trees[i], trees, j, end_j, APPROX);
					}
				}
			}
		}
		// cleanup
		for(vector<Node *>::iterator T2 = trees.begin(); T2 != trees.end(); T2++)
			(*T2)->delete_tree();
	}
	else if (SEQUENCE) {

		Node *T1 = NULL;
		Node *T2 = NULL;
		Node *T1_best_root;
		Node *T2_best_root;
		vector<Node *> T1_rootings = vector<Node *>();
		vector<Node *> T2_rootings = vector<Node *>();
		string line = "";
		string T1_name = "";
		string T2_name = "";
		int n = 0;
		stringstream ss_n;

		// get first tree
		if (!getline(cin, line))
			return 0;
		n++;
		T1_name = strip_newick_name(line);
		strip_trailing_whitespace(T1_name);
		if (T1_name == "") {
			stringstream ss;
			ss << n;
			T1_name = "T" + ss.str();
		}
		if (UNROOTED || SIMPLE_UNROOTED)
			line = root(line);
		T1 = build_tree(line);
		T1->preorder_number();
		T1->edge_preorder_interval();
		if (!QUIET) {
			cout << T1_name << ": ";
			T1->print_subtree();
			cout << endl;
		}
		T1->labels_to_numbers(&label_map, &reverse_label_map);

		// get each next tree
		while (getline(cin, line)) {
			// get next tree
			T2_name = strip_newick_name(line);
			n++;
			strip_trailing_whitespace(T2_name);
			if (T2_name == "") {
				stringstream ss;
				ss << n;
				T2_name = "T" + ss.str();
			}
			if (UNROOTED || SIMPLE_UNROOTED)
				line = root(line);
			T2 = build_tree(line);
			T2->labels_to_numbers(&label_map, &reverse_label_map);
			T2->preorder_number();
			T2->edge_preorder_interval();

			// list of current tree rootings (ALL_UNROOTED)
			T1_best_root = T1;
			if (ALL_UNROOTED) {
				T1_rootings = T1->find_descendants();
			}
			else {
				T1_rootings.clear();
				T1_rootings.push_back(T1);
			}

			// list of next_tree rootings (UNROOTED)
			T2_best_root = T2;
			if (UNROOTED) {
				T2_rootings = T2->find_descendants();
			}
			else {
				T2_rootings.clear();
				T2_rootings.push_back(T2);
			}

			Node *T1_original_root = T1->lchild();
			Node *T2_original_root = T2->lchild();

			int best_distance = INT_MAX;

			// compare each T1 rooting to T2 rooting
			for(int i = 0; i < T1_rootings.size(); i++) {
				if (T1_rootings[i] != T1 && T1_rootings[i] != T1->lchild() && T1_rootings[i] != T1->rchild()) {
					T1->reroot_clean(T1_rootings[i]);
//					string str = T1->str_subtree();
//					T1->delete_tree();
//					T1 = build_tree(str);
//					T1->preorder_number();
//					T1->edge_preorder_interval();
				}
				for(int j = 0; j < T2_rootings.size(); j++) {
					if (T2_rootings[j] != T2 && T2_rootings[j] != T2->lchild() && T2_rootings[j] != T2->rchild()) {
//						cout << T2->str_subtree() << endl;
//						cout << T2_rootings[j]->str_subtree() << endl;
						T2->reroot_clean(T2_rootings[j]);
//						string str = T2->str_subtree();
//						T2->delete_tree();
//						T2 = build_tree(str);
//						T2->preorder_number();
//						T2->edge_preorder_interval();
						}
					int k = rSPR_branch_and_bound_simple_clustering(T1, T2);
					if (k < best_distance) {
						best_distance = k;
						T1_best_root = T1_rootings[i];
						T2_best_root = T2_rootings[j];
						T1->numbers_to_labels(&reverse_label_map);
						T2->numbers_to_labels(&reverse_label_map);
//						cout << i << "," << j << endl;
//						cout << T1->str_subtree() << endl;
//						cout << T2->str_subtree() << endl;
						T1->labels_to_numbers(&label_map, &reverse_label_map);
						T2->labels_to_numbers(&label_map, &reverse_label_map);
//						cout << "d=" << k << endl;
//						cout << "new best" << endl;
//						cout << endl;
					}

					T2->reroot_clean(T2_original_root);
				}
				T1->reroot_clean(T1_original_root);
			}
			if (T1_best_root != T1 && T1_best_root->parent() != T1) {
				T1->reroot_clean(T1_best_root);
			}
			// hack TODO fix
			string str = T1->str_subtree();
			T1->delete_tree();
			T1 = build_tree(str);
			T1->preorder_number();
			T1->edge_preorder_interval();

			if (T2_best_root != T2 && T2_best_root->parent() != T2) {
				T2->reroot_clean(T2_best_root);
			}
			// hack TODO fix
			str = T2->str_subtree();
			T2->delete_tree();
			T2 = build_tree(str);
			T2->preorder_number();
			T2->edge_preorder_interval();

			// show_moves on best rooting pair
			if (SHOW_MOVES) {
				show_moves(T1, T2, &label_map, &reverse_label_map);
			}

			// print tree
			if (!QUIET) {
				cout << endl;
				T2->numbers_to_labels(&reverse_label_map);
				cout << T2_name << ": ";
				T2->print_subtree();
				T2->labels_to_numbers(&label_map, &reverse_label_map);
				cout << endl;
			}

			// new tree becomes old tree
			T1->delete_tree();
			T1 = T2;
			T1_name = T2_name;
			T2 = NULL;
			T2_name = "";
			T1_rootings.clear();
			T1_rootings.clear();
		}


		// cleanup
		if (T1 != NULL) {
			T1->delete_tree();
		}

	}
	return 0;
}

