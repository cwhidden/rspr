/*******************************************************************************
spr_supertree.cpp

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

-support x     Collapse bipartitions with less than x support

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

*******************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <list>
#include <time.h>
#include "rspr.h"

#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "UndoMachine.h"
#include "lgt.h"
#include "sparse_counts.h"
#include "node_glom.h"

using namespace std;

//#define DEBUG_ONE_TREE true


// options to pick default
bool DEFAULT_ALGORITHM=true;
bool DEFAULT_OPTIMIZATIONS=true;
bool DEFAULT_SEARCH_OPTIMIZATIONS=true;


bool FPT = false;
bool QUIET = false;
bool UNROOTED = false;
bool SIMPLE_UNROOTED = false;
bool SIMPLE_UNROOTED_FAST = false;
int SIMPLE_UNROOTED_NUM = INT_MAX;
bool REROOT = false;
bool REROOT_INITIAL = false;
bool APPROX_ROOTING = false;
bool EXACT_ROOTING = false;
bool RANDOM_ROOTING = false;
bool RANDOM_INSERT_ORDER = false;
bool APPROX = false;
bool TIMING = false;
int NUM_ITERATIONS = 25;
bool SMALL_TREES = false;
bool CONVERT_LIST = false;
bool INVALID_TREES = false;
bool VALID_TREES = false;
bool VALID_TREES_ROOTED = false;
bool LGT_ANALYSIS = false;
bool LGT_EVALUATION = false;
bool FIND_MAX_DEGREE = false;
bool MULTI_TREES = false;
int NUM_LEAVES=-1;
int APPROX_SIBLINGS = 0;
int C_SOURCE = -1;
int NUM_SOURCE = -1;
bool RF_TIES = false;
double SUPPORT_THRESHOLD = 0.5;
bool TABOO_SEARCH = false;
bool ONE_TREE_AT_A_TIME = false;
bool NODE_GLOM_CONSTRUCTION = false;
bool USE_PRECOMPUTED_DISTANCES = false;

/*variables Joel added*/
int R_DISTANCE;
bool R_VARIABLE = false;
bool R_LIMIT = false;
bool R_RAND = false;
bool R_CONTROL = false;
bool R_BIAS = false;
float R_PROB;
bool S_LIMIT = false;
int S_NUM = 10;
bool S_STATS = false;
bool D_STATS = false;
bool RANDOM_TREE = false;
bool GREEDY = false;
bool GREEDY_REFINED = false;

string USAGE =
"spr_supertrees, version 1.2.1\n"
"\n"
"usage: spr_supertrees [OPTIONS]\n"
"Calculate binary rooted supertrees that minimize the Subtree Prune and\n"
"Regraft (SPR) distance to a set of rooted or unrooted binary or\n"
"multifurcating trees from STDIN in newick format. Supports arbitrary\n"
"leaf labels. See the README for more information.\n"
"\n"
"Copyright 2011-14 Chris Whidden\n"
"whidden@cs.dal.ca\n"
"http://kiwi.cs.dal.ca/Software/SPR_Supertrees\n"
"March 3, 2014\n"
"Version 1.2.1\n"
"\n"
"This program comes with ABSOLUTELY NO WARRANTY.\n"
"This is free software, and you are welcome to redistribute it\n"
"under certain conditions; See the README for details.\n"
"\n"
"*******************************************************************************\n"
"ALGORITHM\n"
"*******************************************************************************\n"
"\n"
"These options control what algorithm is used to determine the SPR distance\n"
"from the supertree to the input trees. By default -bb is used.\n"
"\n"
"-fpt        Calculate the exact rSPR distance with an FPT algorithm\n"
"\n"
"-bb         Calculate the exact rSPR distance with a branch-and-bound\n"
"            FPT algorithm. This is the default option.\n"
"\n"
"-approx     Calculate just a linear -time 3-approximation of the rSPR distance\n"
"\n"
"-max k      Calculate the exact rSPR distance if it is k or less and\n"
"            otherwise use the 3-approximation\n"
"\n"
"-split_approx\n"
"-split_approx x  Calculate the exact rSPR distance if it is k or less and\n"
"                 otherwise use the exponential-time approximation\n"
"\n"
"*******************************************************************************\n"
"OPTIMIZATIONS\n"
"*******************************************************************************\n"
"\n"
"These options control the use of optimized branching. All optimizations are\n"
"enabled by default. Specifying any subset of -cob, -cab, and -sc will use\n"
"just that subset of optimizations. See the README for more information.\n"
"\n"
"-allopt   Use -cob -cab -sc and a new set of optimizations. This is the default\n"
"          option\n"
"\n"
"-noopt    Use 3-way branching for all FPT algorithms\n"
"\n"
"-cob      Use \"cut one b\" improved branching\n"
"\n"
"-cab      Use \"cut all b\" improved branching\n"
"\n"
"-sc       Use \"separate components\" improved branching\n"
"\n"
"-bipartition_cluster x  Do not consider supertree rearrangements that violate\n"
"                        biparitions supported by x% of gene trees containing\n"
"                        at least two members from each side of the bipartition.\n"
"                        Enabled by default with x=0.5\n"
"\n"
"*******************************************************************************\n"
"MULTIFURCATING COMPARISON OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-allow_multi   Allow multifurcating gene trees\n"
"\n"
"-support x     Collapse bipartitions with less than x support\n"
"\n"
"*******************************************************************************\n"
"UNROOTED COMPARISON OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-unrooted   Compare the supertree to each rooting of the input trees.\n"
"            Use the best found distance\n"
"\n"
"-unrooted_min_approx    Compare the supertree to each rooting of the\n"
"                        input trees.\n"
"                        Run the exact algorithm on the rooting with the\n"
"                        minimum approximate rspr distance\n"
"\n"
"-simple_unrooted        Root the gene trees at each iteration using\n"
"                        a bipartition balanced accuracy measure\n"
"                        (fast but potentially less accurate)\n"
"                        Reports an unrooted SPR distance comparison\n"
"                        at the end of each iteration for comparable\n"
"                        iteration scores\n"
"\n"
"-simple_unrooted x      Root the gene trees at the first x iterations\n"
"\n"
"-simple_unrooted_fast   The same as -simple_unrooted but does not use\n"
"                        an unrooted comparison at the end of each\n"
"                        iteration\n"
"\n"
"-outgroup FILE          Root the gene trees with the outgroup taxa\n"
"                        listed in FILE, one per line. Trees with a\n"
"                        polyphyletic outgroup are considered invalid.\n"
"\n"
"-reroot                 Reroot the super tree at each iteration using\n"
"                        a bipartition balanced accuracy measure\n"
"\n"
"-rspr_reroot            Root trees using the SPR distance instead\n"
"                        of the bipartition balanced accuracy\n"
"\n"
"\n"
"\n"
"*******************************************************************************\n"
"SEARCH STRATEGY OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-i x    Run for x iterations of the global rearrangement search\n"
"\n"
"-r x    Only consider transfers of length x in the global rearrangement\n"
"        search. Default is infinite (All SPRs). For NNI search use\n"
"        -r 1\n"
"\n"
"-include_only <file>  Build the supertree only from taxa included in\n"
"                      <file>, one per line\n"
"\n"
"-initial_tree <file>  Begin the search with the tree in <file>\n"
"\n"
"-num_leaves x         Build the supertree from the x taxa that are found\n"
"                      in the largest number of trees\n"
"\n"
"-random_insert_order  Insert taxa in random order when building the\n"
"                      greedy addition tree. The default order is\n"
"                      descending occurence\n"
"\n"
"-rf_ties              Break SPR distance ties with the RF distance\n"
"\n"
"*******************************************************************************\n"
"LGT ANALYSIS\n"
"*******************************************************************************\n"
"\n"
"-lgt_analysis          Conduct an LGT analysis with the initial user-specified\n"
"                       or greedy addition tree\n"
"\n"
"-lgt_csv               Output the LGT analysis seperated by commas rather than\n"
"                       spaces.\n"
"\n"
"-lgt_groups FILE       Specify a set of groups (e.g. genus or class) to analyze\n"
"                       with -lgt_analysis. The group FILE contains a set of\n"
"                       groups consisting of a group name on one line, group\n"
"                       members one per line, and a blank line to seperate each\n"
"                       group.\n"
"                       \n"
"*******************************************************************************\n"
"OTHER OPTIONS\n"
"*******************************************************************************\n"
"-time                  Print iteration and total CPU time used at each\n"
"                       iteration\n"
"\n"
"-cc                    Calculate a potentially better approximation with a\n"
"                       quadratic time algorithm\n"
"\n"
"-valid_trees           Output the set of trees that appear valid\n"
"-valid_trees_rooted    Output the set of trees that appear valid after applying\n"
"                       any rooting options.\n"
"\n"
"-multi_trees           Output the set of multifurcating or invalid trees\n";

Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees,
		int label);
Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees,
		vector<Node *> *best_siblings, int label);
void find_best_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &min_tie_distance,
		int &num_ties, Node **best_sibling);
vector<Node *> *find_best_siblings(Node *super_tree, vector<Node *> &gene_trees, int label, int num_siblings);
void find_best_siblings_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &min_tie_distance,
		int &num_ties, multimap<int, Node*> *best_siblings, int num_siblings);
void test_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &min_tie_distance,
		int &num_ties, Node **best_sibling);
void find_best_spr(Node *super_tree, vector<Node *> &gene_trees,
		Node *&best_spr_move, Node *&best_sibling);
void find_best_spr_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties);
void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties);
void get_support(Node *super_tree, vector<Node *> *gene_trees);
void get_support(Node *n, Node *super_tree, vector<Node *> *gene_trees);
void get_transfer_support(Node *super_tree, vector<Node *> *gene_trees);
void get_transfer_support(Node *n, Node *super_tree, vector<Node *> *gene_trees);
void get_bipartition_support(Node *super_tree, vector<Node *> *gene_trees,
		enum RELAXATION relaxed);
bool supported_spr(Node *source, Node *target);
bool pair_comparator (pair<int, int> a, pair<int, int> b);

/*Prototypes of Joel's functions*/
void find_best_spr_r(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, int r);
void find_best_spr_r(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, int r, vector<int> *original_scores);
void find_best_spr_r_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, int r);
void find_best_spr_r_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &min_tie_distance,
		int &num_ties, int r,
		vector<int> *original_scores);
void find_best_spr_r_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, int r, int origin);
void find_best_spr_r_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &min_tie_distance,
		int &num_ties, int r, int origin, int offset);
void find_best_spr(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, vector <pair <int, pair<int, int> > > &stats);
void find_best_spr_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <int, pair<int, int> > > &stats);
void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <int, pair<int, int> > > &stats);

void find_best_spr(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, vector<pair <pair<Node*,Node*>, int> > &approx_moves);
void find_best_spr_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <pair<Node*,Node*>, int> > &approx_moves);
void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <pair<Node*,Node*>, int> > &approx_moves);
void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <pair<Node*,Node*>, int> > &approx_moves, int offset);
bool sort_approx_moves(const pair<pair<Node*,Node*>, int> &a, const pair<pair<Node*,Node*>, int> &b);
int find_r();
int find_r(double probability);

void find_best_distance(Node * n, Node * super_tree, vector<Node *> &gene_trees, vector< pair<Node*, int> > &scores, int &num_zeros, int &best_distance);
//void find_best_distance_helper(Node *n, Forest *f, vector<Node *> &gene_trees, vector< pair<Node *, int> > &scores);
vector<pair<Node *, int> > get_best_scores(vector<pair<Node *, int> > &scores);
void find_best_spr_r_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, int r, int origin, vector<pair <pair<Node*,Node*>, int> > &approx_moves);
bool is_taboo(list<Node> taboo_trees, Node *super_tree);

	map<string, int> label_map;
	map<int, string> reverse_label_map;

bool is_pow_2(int n) {
	return (n) && !(n & (n - 1));
}

bool BIPARTITION_CLUSTER = false;

list<Node> taboo_trees = list<Node>();

int main(int argc, char *argv[]) {

	// ignore multifurcating trees by default
	IGNORE_MULTI = true;
	string INCLUDE_ONLY = "";
	bool OUTGROUP_ROOT = false;
	string OUTGROUP = "";
	string INITIAL_SUPER_TREE = "";
	string LGT_GROUPS = "";
	bool INITIAL_SUPER_TREE_UNROOTED = false;
	bool FIND_SUPPORT = false;
	bool FIND_BIPARTITION_SUPPORT = false;
	enum RELAXATION RELAXED_BIPARTITION_SUPPORT = NEGATIVE_RELAXED;
	bool FIND_CLADE_TRANSFERS = false;
	bool LGT_CSV = false;

	int max_args = argc-1;
	while (argc > 1) {
		char *arg = argv[--argc];
		if (strcmp(arg, "-fpt") == 0) {
			FPT = true;
			DEFAULT_ALGORITHM=false;
		}
		else if (strcmp(arg, "-bb") == 0) {
			BB = true;
			DEFAULT_ALGORITHM=false;
		}
		else if (strcmp(arg, "-approx") == 0) {
			APPROX=true;
		}
		else if (strcmp(arg, "-q") == 0)
			QUIET = true;
		else if (strcmp(arg, "-cc") == 0)
			APPROX_CHECK_COMPONENT = true;
		else if (strcmp(arg, "-unrooted") == 0)
			UNROOTED = true;
		else if (strcmp(arg, "-simple_unrooted") == 0) {
			SIMPLE_UNROOTED = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					SIMPLE_UNROOTED_NUM = atoi(arg2);
				cout << "SIMPLE_UNROOTED_NUM=" << SIMPLE_UNROOTED_NUM << endl;
			}
		}
		else if (strcmp(arg, "-simple_unrooted_rspr") == 0) {
			SIMPLE_UNROOTED = true;
			EXACT_ROOTING = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					SIMPLE_UNROOTED_NUM = atoi(arg2);
				cout << "SIMPLE_UNROOTED_NUM=" << SIMPLE_UNROOTED_NUM << endl;
			}
		}
		else if (strcmp(arg, "-simple_unrooted_random") == 0) {
			SIMPLE_UNROOTED = true;
			RANDOM_ROOTING = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					SIMPLE_UNROOTED_NUM = atoi(arg2);
				cout << "SIMPLE_UNROOTED_NUM=" << SIMPLE_UNROOTED_NUM << endl;
			}
		}
		else if (strcmp(arg, "-simple_unrooted_fast") == 0) {
			SIMPLE_UNROOTED=true;
			SIMPLE_UNROOTED_FAST=true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					SIMPLE_UNROOTED_NUM = atoi(arg2);
				cout << "SIMPLE_UNROOTED_NUM=" << SIMPLE_UNROOTED_NUM << endl;
			}
		}
		else if (strcmp(arg, "-random_insert_order") == 0) {
			RANDOM_INSERT_ORDER = true;
		}
		else if (strcmp(arg, "-reroot") == 0) {
			REROOT = true;
			REROOT_INITIAL = true;
		}
		else if (strcmp(arg, "-reroot_initial") == 0) {
			REROOT_INITIAL = true;
		}
		else if (strcmp(arg, "-rspr_root") == 0) {
			EXACT_ROOTING = true;
		}
		else if (strcmp(arg, "-random_root") == 0) {
			RANDOM_ROOTING = true;
		}
		else if (strcmp(arg, "-unrooted_min_approx") == 0) {
			UNROOTED = true;
			UNROOTED_MIN_APPROX = true;
		}
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
		else if (strcmp(arg, "-prefer_rho") == 0) {
			PREFER_RHO = true;
		}
/*
		else if (strcmp(arg, "-memoize") == 0) {
			MEMOIZE = true;
		}
		else if (strcmp(arg, "-all_mafs") == 0) {
			ALL_MAFS= true;
		}
*/
		else if (strcmp(arg, "-i") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					NUM_ITERATIONS = atoi(arg2);
				cout << "NUM_ITERATIONS=" << NUM_ITERATIONS << endl;
			}
		}
/*Joel: limit SPR radius*/
		else if (strcmp(arg, "-r") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] == 'r'){
					R_LIMIT = true;
					R_RAND = true;
					R_DISTANCE = find_r();
					cout << "RANDOM R_DISTANCE=" << R_DISTANCE << endl;
				}
				else if (arg2[0] != '-'){
					R_LIMIT = true;
					R_DISTANCE = atoi(arg2);
					cout << "R_DISTANCE=" << R_DISTANCE << endl;
				}

			}
		}
		else if (strcmp(arg, "-r_variable") == 0) {
			R_VARIABLE = true;
			DEFAULT_SEARCH_OPTIMIZATIONS=false;
		}

		else if (strcmp(arg, "-p") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (!R_RAND)
					R_PROB = 0.500;
				else if (arg2[0] != '-'){
					R_BIAS = true;
					R_PROB = atof(arg2);
					cout << "PROBABILITY=" << R_PROB << endl;
				}
			}	
			DEFAULT_SEARCH_OPTIMIZATIONS=false;
		}
/*Joel: refined greedy search*/
		else if(strcmp(arg, "-r_greedy")==0){
			GREEDY_REFINED = true;
			if(max_args > argc){
				char *arg2 = argv[argc+1];
				if(arg2[0] != '-')
					S_NUM = atoi(arg2);
				else
					S_NUM = 5;
			}
			if(max_args > argc){
				char *arg2 = argv[argc+2];
				if(arg2[0] != '-')
					R_DISTANCE = atoi(arg2);
				else
					R_DISTANCE = 1;
			}
			DEFAULT_SEARCH_OPTIMIZATIONS=false;
		}
/*Joel: limit starting points*/
		else if (strcmp(arg, "-num_start") == 0){
			S_LIMIT = true;
			if(max_args > argc) {
				char *arg2 = argv[argc+1];
				if(arg2[0] != '-'){
					S_NUM = atoi(arg2);
				}
			}
			cout << "STARTING POINTS TO CHECK: " << S_NUM << endl;
			DEFAULT_SEARCH_OPTIMIZATIONS=false;
		}
/*Joel: limit SPR radius options*/
		else if (strcmp(arg, "-stats") == 0){
			S_STATS = true;
		}
		else if (strcmp(arg, "-control") == 0) {
			R_CONTROL = true;
			DEFAULT_SEARCH_OPTIMIZATIONS=false;
		}
/*Joel: stats on approx distance vs total distance*/
		else if (strcmp(arg, "-d_stats") == 0)
			D_STATS = true;
/*Joel: greedy search*/
		else if(strcmp(arg, "-greedy")==0) {
			GREEDY = true;
			DEFAULT_SEARCH_OPTIMIZATIONS=false;
		}
/*Joel: random starting tree*/
		else if (strcmp(arg, "-rand_tree") == 0)
			RANDOM_TREE = true;			
		else if (strcmp(arg, "-max") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					MAX_SPR = atoi(arg2);
				cout << "MAX_SPR=" << MAX_SPR << endl;
			}
		}
		else if (strcmp(arg, "-cluster_max") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-') {
					CLUSTER_MAX_SPR = atoi(arg2);
				cout << "CLUSTER_MAX_SPR=" << CLUSTER_MAX_SPR << endl;
				}
			}
		}
		else if (strcmp(arg, "-num_leaves") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					NUM_LEAVES = atoi(arg2);
				cout << "NUM_LEAVES=" << NUM_LEAVES << endl;
			}
		}
		else if (strcmp(arg, "-time") == 0) {
			TIMING= true;
		}
		else if (strcmp(arg, "-clamp") == 0) {
			CLAMP= true;
		}
		else if (strcmp(arg, "-small_trees") == 0) {
			SMALL_TREES=true;
		}
		else if (strcmp(arg, "-convert_list") == 0) {
			CONVERT_LIST=true;
		}
		else if (strcmp(arg, "-valid_trees") == 0) {
			VALID_TREES=true;
		}
		else if (strcmp(arg, "-invalid_trees") == 0) {
			VALID_TREES=true;
			INVALID_TREES=true;
		}
		else if (strcmp(arg, "-valid_trees_rooted") == 0) {
			VALID_TREES_ROOTED=true;
		}
		else if (strcmp(arg, "-lgt_analysis") == 0) {
			LGT_ANALYSIS=true;
		}
		else if (strcmp(arg, "-lgt_evaluation") == 0) {
			LGT_EVALUATION=true;
		}
		else if (strcmp(arg, "-lgt_csv") == 0) {
			LGT_CSV=true;
		}
		else if (strcmp(arg, "-max_degree") == 0) {
			FIND_MAX_DEGREE=true;
		}
		else if (strcmp(arg, "-multi_trees") == 0) {
			MULTI_TREES=true;
		}
		else if (strcmp(arg, "-rf_ties") == 0) {
			RF_TIES=true;
		}
		else if (strcmp(arg, "-allow_multi") == 0) {
			IGNORE_MULTI = false;
		}
		else if (strcmp(arg, "-protect_edges") == 0) {
			EDGE_PROTECTION = true;
			DEFAULT_OPTIMIZATIONS=false;
		}
		else if (strcmp(arg, "-allow_abort") == 0) {
			ABORT_AT_FIRST_SOLUTION = true;
			DEFAULT_OPTIMIZATIONS=false;
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
				cout << "SPLIT_APPROX_THRESHOLD=" << SPLIT_APPROX_THRESHOLD
						<< endl;
			}
		}
		else if (strcmp(arg, "-support") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					REQUIRED_SUPPORT = atof(arg2);
				cout << "REQUIRED_SUPPORT=" << REQUIRED_SUPPORT
						<< endl;
			}
		}
		else if (strcmp(arg, "-approx_siblings") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					APPROX_SIBLINGS = atoi(arg2);
				cout << "APPROX_SIBLINGS=" << APPROX_SIBLINGS
						<< endl;
			}
		}
		else if (strcmp(arg, "-include_only") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					INCLUDE_ONLY = string(arg2);
				cout << "INCLUDE_ONLY=" << INCLUDE_ONLY
						<< endl;
			}
		}
		else if (strcmp(arg, "-outgroup") == 0) {
			OUTGROUP_ROOT = true;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					OUTGROUP = string(arg2);
				cout << "OUTGROUP=" << OUTGROUP
						<< endl;
			}
		}
		else if (strcmp(arg, "-initial_tree") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					INITIAL_SUPER_TREE = string(arg2);
				cout << "INITIAL_TREE=" << INITIAL_SUPER_TREE
						<< endl;
			}
		}
		else if (strcmp(arg, "-lgt_groups") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					LGT_GROUPS = string(arg2);
				cout << "LGT_GROUPS=" << LGT_GROUPS 
						<< endl;
			}
		}
		else if (strcmp(arg, "-initial_tree_unrooted") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					INITIAL_SUPER_TREE = string(arg2);
				cout << "INITIAL_TREE=" << INITIAL_SUPER_TREE
						<< endl;
				INITIAL_SUPER_TREE_UNROOTED=true;
				REROOT_INITIAL=true;
			}
		}
		else if (strcmp(arg, "-count_losses") == 0) {
			COUNT_LOSSES = true;
		}
		else if (strcmp(arg, "-cut_lost") == 0) {
			CUT_LOST = true;
		}
		else if (strcmp(arg, "-find_support") == 0) {
			FIND_SUPPORT = true;
		}
		else if (strcmp(arg, "-find_bipartition_support") == 0) {
			FIND_BIPARTITION_SUPPORT = true;
		}
		else if (strcmp(arg, "-relaxed_bipartition_support") == 0) {
			RELAXED_BIPARTITION_SUPPORT = ALL_RELAXED;
		}
		else if (strcmp(arg, "-strict_bipartition_support") == 0) {
			RELAXED_BIPARTITION_SUPPORT = STRICT;
		}
		else if (strcmp(arg, "-bipartition_cluster") == 0) {
			BIPARTITION_CLUSTER = true;
			cout << "BIPARTITION_CLUSTER=true" << endl;
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					SUPPORT_THRESHOLD = atof(arg2);
			}
			cout << "CLUSTER THRESHOLD=" << SUPPORT_THRESHOLD
					<< endl;
		}
		else if (strcmp(arg, "-find_clade_transfers") == 0) {
			FIND_CLADE_TRANSFERS = true;
		}
		else if (strcmp(arg, "-taboo_search") == 0) {
			TABOO_SEARCH = true;
		}
		else if (strcmp(arg, "-one_tree_at_a_time") == 0
				|| (strcmp(arg, "-one_tree") == 0) ) {
			ONE_TREE_AT_A_TIME = true;
		}
		else if (strcmp(arg, "-node_glom_construction") == 0
				|| (strcmp(arg, "-node_glom") == 0) ) {
			NODE_GLOM_CONSTRUCTION = true;
		}
		else if (strcmp(arg, "-precompute") == 0 ) {
			USE_PRECOMPUTED_DISTANCES = true;
		}
		else if (strcmp(arg, "-cluster_tune") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-') {
					CLUSTER_TUNE = atoi(arg2);
				}
			}
		}
		else if (strcmp(arg, "--help") == 0) {
			cout << USAGE;
			return 0;
		}
			
	}
	if (DEFAULT_OPTIMIZATIONS) {
		CUT_ALL_B=true;
		CUT_ONE_B = true;
		CUT_TWO_B = true;
//		CUT_TWO_B_ROOT = true;
		REVERSE_CUT_ONE_B = true;
		REVERSE_CUT_ONE_B_3 = true;
		CUT_AC_SEPARATE_COMPONENTS = true;
		EDGE_PROTECTION = true;
		EDGE_PROTECTION_TWO_B = true;
		PREFER_NONBRANCHING = true;

//		CHECK_MERGE_DEPTH = true;
//		ABORT_AT_FIRST_SOLUTION = true;
		PREORDER_SIBLING_PAIRS = true;
		NEAR_PREORDER_SIBLING_PAIRS = true;
		LEAF_REDUCTION = true;
		LEAF_REDUCTION2 = true;

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
			CLUSTER_TUNE = 30;
		}
	}
	if (DEFAULT_SEARCH_OPTIMIZATIONS) {
		if (BIPARTITION_CLUSTER == false) {
			BIPARTITION_CLUSTER = true;
		}
		if (R_LIMIT == false) {
			R_LIMIT = true;
			R_DISTANCE = INT_MAX - 1000;
		}
	}
	PREORDER_SIBLING_PAIRS = true;
	if (DEFAULT_ALGORITHM) {
		BB=true;
		PREFER_RHO = true;
	}



	// initialize random number generator
	srand((unsigned(time(0))));

	// Label maps to allow string labels
	vector<int> label_counts = vector<int>();
	label_map= map<string, int>();
	reverse_label_map = map<int, string>();

	string T_line = "";
	vector<Node *> gene_trees = vector<Node *>();
//	multimap<int, pair<Node*, string> > gene_tree_map
//		= multimap<int, pair<Node*, string> >();
	vector<string> gene_tree_names = vector<string>();
	set<string, StringCompare> include_only;
	set<string, StringCompare> outgroup;
	if (INCLUDE_ONLY != "") {
		include_only.insert("");
		ifstream include_only_file;
		include_only_file.open(INCLUDE_ONLY.c_str());
		if (include_only_file.is_open()) {
			string line;
			while(include_only_file.good()) {
				getline(include_only_file, line);
//				cout << "\"" << line << "\"" << endl;
				include_only.insert(line);
			}
//			cout << endl;
			include_only_file.close();
			for(set<string, StringCompare>::iterator i = include_only.begin();
					i != include_only.end(); i++) {
//				cout << "\"" << *i << "\"" << endl;
			}
		}
		else
			INCLUDE_ONLY = "";
	}
	if (OUTGROUP != "") {
		outgroup.insert("");
		ifstream outgroup_file;
		outgroup_file.open(OUTGROUP.c_str());
		if (outgroup_file.is_open()) {
			string line;
			while(outgroup_file.good()) {
				getline(outgroup_file, line);
//				cout << "\"" << line << "\"" << endl;
				outgroup.insert(line);
			}
//			cout << endl;
			outgroup_file.close();
			for(set<string, StringCompare>::iterator i = outgroup.begin();
					i != outgroup.end(); i++) {
//				cout << "\"" << *i << "\"" << endl;
			}
		}
		else
			OUTGROUP = "";
	}
	int skipped_multifurcating = 0;
	int skipped_small = 0;
	int skipped_no_bracket = 0;
	int skipped_star = 0;
	int skipped_no_outgroup = 0;
	while (getline(cin, T_line)) {
		string name = "";
		size_t loc = T_line.find_first_of("(");
		if (loc != string::npos) {
			if (loc != 0) {
				name = T_line.substr(0,loc);
				T_line.erase(0,loc);
			}
			// TODO: should we be doing this with OUTGROUP?
			if ((UNROOTED || SIMPLE_UNROOTED || OUTGROUP_ROOT)
					&& IGNORE_MULTI) {
				T_line = root(T_line);
			}
			Node *T;
			if (INCLUDE_ONLY != "")
				T = build_tree(T_line, &include_only);
			else
				T = build_tree(T_line);
			if (UNROOTED || SIMPLE_UNROOTED || OUTGROUP_ROOT) {
				T->fixroot();
			}
//			cout << T->str_subtree() << endl;
			// TODO: check that this works
			/*
			if ((UNROOTED || SIMPLE_UNROOTED) &&
					T->get_children().size() > 2) {
//					cout << T->str_subtree() << endl;
					Node *new_root = T->get_children().back();
					new_root->cut_parent();
					if (new_root->is_leaf()) {
						Node *real_root = new Node();
						real_root->add_child(new_root);
						new_root = real_root;
					}

					while(new_root->get_children().size() > 1) {
						Node *new_child = new_root->get_children().back();
						new_child->cut_parent();
						T->add_child(new_child);
					}
					new_root->add_child(T);
					T = new_root;
//					cout << T->str_subtree() << endl;
//					cout << endl;
			}
			*/

			if (T->is_leaf() && T->str() == "p") {
				if (MULTI_TREES)
					cout << name << T_line << endl;

				skipped_multifurcating++;
				continue;
			}
			int T_size = T->size();
			if (!INVALID_TREES && ((T_size <= 4)
					|| T_size == 5 && !SMALL_TREES)) {
				skipped_small++;
				continue;
			}
			if (!IGNORE_MULTI && !INVALID_TREES) {
				int T_depth = T->max_depth();
				if (T_depth <= 1 ||
						((UNROOTED || SIMPLE_UNROOTED) && T_depth <= 2)) {
					skipped_star++;
					continue;
				//cout << T->str_subtree() << endl;
				//cout << T_depth << endl;
				}
			}
			if (UNROOTED || SIMPLE_UNROOTED || OUTGROUP_ROOT)
				T->preorder_number();

			if (OUTGROUP != "") {
				// TODO: write a function that will root the tree based on the
				// outgroup
				// TODO: do we simply skip the tree if the outgroup is seperate?
				set<string>::iterator i;
				if (!outgroup_root(T, outgroup)) {
					skipped_no_outgroup++;
					continue;
				}
				else {
					T->preorder_number();
				}

				
			}

			gene_tree_names.push_back(name);
			gene_trees.push_back(T);
			//gene_tree_map.insert(make_pair(T->size(), make_pair(T, name)));
		}
		else
			skipped_no_bracket++;
	}

	cout << "skipped " << skipped_no_bracket << " lines with no opening bracket " << endl;
	if (IGNORE_MULTI)
		cout << "skipped " << skipped_multifurcating << " multifurcating or invalid trees" << endl;
	else
	cout << "skipped " << skipped_multifurcating << " invalid trees" << endl;
	if (SMALL_TREES) {
		cout << "skipped " << skipped_small << " trees with less than 3 leaves" << endl;
	}
	else {
		cout << "skipped " << skipped_small << " trees with less than 4 leaves" << endl;
		if (!IGNORE_MULTI)
		cout << "skipped " << skipped_star << " star trees" << endl;
	}
	if (OUTGROUP_ROOT)
		cout << "skipped " << skipped_no_outgroup << " trees with no outgroup or an unresolved outgroup" << endl;

/*	int end = gene_tree_map.size();
	for(int i = 0; i < end; i++) {
		multimap<int,  pair<Node *, string> >::iterator p =
			gene_tree_map.begin();
		gene_trees.push_back(p->second.first);
		gene_tree_names.push_back(p->second.second);
		gene_tree_map.erase(p);
	}
*/


	cout << gene_trees.size() << " gene trees remaining" << endl;

	if (gene_trees.size() <= 0)
		exit(0);

	for(int i = 0; i < gene_trees.size(); i++) {
		if (VALID_TREES) {
			cout << gene_tree_names[i];
			cout << gene_trees[i]->str_subtree() << endl;
		}
		if (FIND_MAX_DEGREE) {
			cout << gene_tree_names[i];
			cout << gene_trees[i]->max_degree() << endl;
		}

		gene_trees[i]->labels_to_numbers(&label_map, &reverse_label_map);
//		cout << gene_tree_names[i] << gene_trees[i]->str_subtree() << endl;

		gene_trees[i]->count_numbered_labels(&label_counts);
	}
	if (VALID_TREES || MULTI_TREES)
		exit(0);


	// iterate over the taxa by number of occurences
	multimap<int, int> labels = multimap<int, int>();
	for(int i = 0; i < label_counts.size(); i++) {
		int key = label_counts[i];
		if (RANDOM_INSERT_ORDER)
			key = rand();
		labels.insert(make_pair(key,i));
	}

	if (CONVERT_LIST) {
		for(multimap<int, int>::reverse_iterator i = labels.rbegin(); i != labels.rend(); i++) {
			int num = i->second;
			string name = reverse_label_map.find(num)->second;
			cout << num << ","
				<< name <<  endl;
		}
		exit(0);
	}

	int best_distance;
	int best_rooted_distance = INT_MAX;
	int best_tie_distance = 0;
	Node *super_tree;
	multimap<int, int>::reverse_iterator label = labels.rbegin();

	double time;
	double current_time;
	if (TIMING)
		time = clock()/(double)CLOCKS_PER_SEC;

	if (NODE_GLOM_CONSTRUCTION) {

		// copy the gene trees
		vector<Node *> gene_trees_copy = vector<Node *>(gene_trees.size());
		for (int i = 0; i < gene_trees_copy.size(); i++) {
			gene_trees_copy[i] = new Node(*gene_trees[i]);
		}


		// vector of partially joined trees
		vector<Node *> super_forest = vector<Node *>(label_counts.size());
		for (int i = 0; i < label_counts.size(); i++) {
			super_forest[i] = new Node(itos(i));
		}

		// partition neighbour counts as a sparse matrix (vector of maps)
		// TODO: make this a double and downweight multiple counts from same
		// tree
		SparseCounts<double> neighbour_counts =
				SparseCounts<double>(label_counts.size(), label_counts.size());
		vector<int> total_leaf_counts = vector<int>(label_counts.size(),0);

		// TODO: do we need a reverse mapping from old numbers to new?

		// glom the components
		//for(int i = 0; i < label_counts.size(); i++) {
		// }
		for(int i = 0; i < label_counts.size()-1; i++) {
			// count neighbours
			cout << "Iteration " << i+1 << " / " << label_counts.size()-1 << endl;
			neighbour_counts.clear();
			for (int j = 0; j < gene_trees.size(); j++ ) {
				vector<int> leaf_counts = vector<int>(label_counts.size(),0);
				count_leaves(gene_trees[j], &leaf_counts);
/*				bool tree_agrees = true;
				for(int k = 0; k < leaf_counts.size(); k++) {
					if (leaf_counts[k] > 1)
						tree_agrees = false;
				}
				if (tree_agrees)
*/
					count_neighbours(gene_trees[j], &neighbour_counts, &leaf_counts);
				for(int k = 0; k < leaf_counts.size(); k++) {
					total_leaf_counts[k] += leaf_counts[k];
				}
//				count_neighbours(gene_trees[j], &neighbour_counts, NULL);
			}
/*
			neighbour_counts.sparse_labelled_print(&reverse_label_map);
			cout << endl;
*/

			// scale the scores by the log of the component size
			vector<double> component_scale = vector<double>(super_forest.size());
			for(int j = 0; j < component_scale.size(); j++) {
				if (super_forest[j] != NULL)
					component_scale[j] = super_forest[j]->max_depth();
//					component_scale[j] = log(super_forest[j]->max_depth() + 1);
				else
					component_scale[j] = 1;
			}


			// scale the scores by the number of common genes
			vector<vector<int> > component_trees = find_component_trees(&gene_trees, &super_forest, label_counts.size());
			// find the most common pair
			vector<pair<int, int> > mcp_vector =
//				neighbour_counts.find_most_common_pairs_scaled(&component_trees);
				neighbour_counts.find_most_common_pairs_scaled(&component_scale, &component_trees, &total_leaf_counts);
/*			cout << "MCP: ";

			for(int j = 0; j < mcp_vector.size(); j++) {
			}
*/
			// break ties randomly for now
			pair<int, int> mcp = mcp_vector[rand() % mcp_vector.size()];
/*
			cout << "\t" << reverse_label_map.find(mcp.first)->second << ", " << reverse_label_map.find(mcp.second)->second << endl;
			cout << "\t" << mcp.first << ", " << mcp.second << endl;
*/

			// join the most common pair in the super_forest
/*
			for (int j = 0; j < super_forest.size(); j++) {
				if (super_forest[j] == NULL) {
					cout << "*";
				}
				else {
					cout << super_forest[j]->str_subtree();
				}
				cout << "   ";
			}
			cout << endl;
*/
			glom_super_forest(&super_forest, mcp.first, mcp.second);

			for (int j = 0; j < super_forest.size(); j++) {
				if (super_forest[j] == NULL) {
					cout << "*";
				}
				else {
					super_forest[j]->numbers_to_labels(&reverse_label_map);
					cout << super_forest[j]->str_subtree();
					super_forest[j]->labels_to_numbers(&label_map, &reverse_label_map);
				}
				cout << "   ";
			}
			cout << endl;


			for (int j = 0; j < gene_trees.size(); j++) {
/*				cout << "Gene Tree" << j << endl;
				cout << gene_trees[j]->str_subtree() << endl;
				cout << endl;
*/
				glom_gene_tree(gene_trees[j], mcp.first, mcp.second);
/*				cout << gene_trees[j]->str_subtree() << endl;
				cout << endl;
*/
				// join the most common pair
				// give any isolated partition members of either the
				//	lower number
			}
		}

		
/*
TODO: 
				SparseCounts functions
				process chosen pair in each gene tree
				super_forest cleanup (should just be one tree)
				Do we want to double the storage but allow quick looks at
					one neighbourhood?
				Store the most common pair in the sparse_counts structure
					and find it as we go along?
						* what about tiebreaking? random is fine but what if
							we want neighbourhood lookups?
							* we could store the largest X values and only consider
								those
				Keep a node to tree mapping and only look at trees that contain
				mcp.second
					* add those trees to the mapping for mcp.first
*/

		// cleanup
		super_tree = new Node("0");
		bool found = false;
		for(int i = 0; i < super_forest.size(); i++) {
			if (super_forest[i] != NULL) {
				if (found) {
					cout << "UH-OH! Some components not glommed together" << endl;
				}
				delete super_tree;
				super_tree = super_forest[i];
				super_forest[i] = NULL;
				found = true;
			}
		}
		for(int i = 0; i < gene_trees.size(); i++) {
			gene_trees[i]->delete_tree();
		}
		gene_trees = gene_trees_copy;
	}
	else {

		int leaf_num = 1;
	
	//	for(auto label = labels.rbegin(); label != labels.rend(); label++) {
	//		cout << label->second ;
	//		cout << " (" << reverse_label_map.find(label->second)->second << ") ";
	//		cout << ": " << label->first << endl;
	//	}
		if (INITIAL_SUPER_TREE != "") {
			ifstream super_tree_file;
			super_tree_file.open(INITIAL_SUPER_TREE.c_str());
			if (super_tree_file.is_open()) {
				string line;
				if(super_tree_file.good()) {
					getline(super_tree_file, line);
					if (INITIAL_SUPER_TREE_UNROOTED)
						line = root(line);
					if (INCLUDE_ONLY != "")
						super_tree = build_tree(line, &include_only);
					else
						super_tree = build_tree(line);
					super_tree->preorder_number();
					super_tree->labels_to_numbers(&label_map, &reverse_label_map);
				}
				super_tree_file.close();
			}
			else
				INITIAL_SUPER_TREE = "";
		}
		if (INITIAL_SUPER_TREE == "") {
	
			// 4 most common leaves
			vector<int> l = vector<int>(4);
			for(int i = 0; i < 4; i++) {
				l[i] = label->second;
				label++;
				//cout << l[i] << endl;
			}
		
		//	cout << endl;
		//	cout << "Starting Trees:" << endl;
		
			// create all trees on the 4 leaves
			vector<Node *> ST = vector<Node *>();
			int st_size = 0;
			stringstream ss;
			ss << "((" << l[0] << "," << l[1] << "),(" << l[2] << "," << l[3] <<"))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "((" << l[0] << "," << l[2] << "),(" << l[1] << "," << l[3] <<"))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "((" << l[0] << "," << l[3] << "),(" << l[1] << "," << l[2] <<"))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
		
			// testing other rooted trees
			ss << "(" << l[0] << "," << "(" << l[3] << ",(" << l[1] << "," << l[2] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[0] << "," << "(" << l[1] << ",(" << l[3] << "," << l[2] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[0] << "," << "(" << l[2] << ",(" << l[1] << "," << l[3] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			// testing other rooted trees
			ss << "(" << l[1] << "," << "(" << l[3] << ",(" << l[0] << "," << l[2] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[1] << "," << "(" << l[0] << ",(" << l[3] << "," << l[2] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[1] << "," << "(" << l[2] << ",(" << l[0] << "," << l[3] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			// testing other rooted trees
			ss << "(" << l[2] << "," << "(" << l[3] << ",(" << l[1] << "," << l[0] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[2] << "," << "(" << l[1] << ",(" << l[3] << "," << l[0] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[2] << "," << "(" << l[0] << ",(" << l[1] << "," << l[3] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			// testing other rooted trees
			ss << "(" << l[3] << "," << "(" << l[0] << ",(" << l[1] << "," << l[2] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[3] << "," << "(" << l[1] << ",(" << l[0] << "," << l[2] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
			ss << "(" << l[3] << "," << "(" << l[2] << ",(" << l[1] << "," << l[0] <<")))";
			ST.push_back(build_tree(ss.str()));
			ss.str("");
		
				vector<Node *> current_gene_trees = vector<Node *>();
				for(int i = 0; i < gene_trees.size(); i++) {
					if (gene_trees[i]->contains_leaf(l[0]))
						current_gene_trees.push_back(gene_trees[i]);
					else if (gene_trees[i]->contains_leaf(l[1]))
						current_gene_trees.push_back(gene_trees[i]);
					else if (gene_trees[i]->contains_leaf(l[2]))
						current_gene_trees.push_back(gene_trees[i]);
					else if (gene_trees[i]->contains_leaf(l[3]))
						current_gene_trees.push_back(gene_trees[i]);
				}
		
		
			// choose the best starting tree
	/*Joel: Random starting tree*/
			if (RANDOM_TREE)
				best_distance = 0;
			else {
				if (APPROX)
					if (UNROOTED)
						best_distance = rSPR_total_approx_distance_unrooted(ST[0], current_gene_trees); 
					else
						best_distance = rSPR_total_approx_distance(ST[0], current_gene_trees); 
				else 
					if (UNROOTED)
						best_distance = rSPR_total_distance_unrooted(ST[0], current_gene_trees); 
					else
						best_distance = rSPR_total_distance(ST[0], current_gene_trees); 
			}
			int best_tree = 0;
		//	cout << best_distance <<  ": ";
		//	cout << ST[0]->str_subtree() << endl;
			for(int j = 1; j < ST.size(); j++) {
				//cout << ST[j]->str_subtree() << endl;
				int distance;
				if (APPROX)
					if (UNROOTED)
						distance = rSPR_total_approx_distance_unrooted(ST[j], current_gene_trees);
					else
						distance = rSPR_total_approx_distance(ST[j], current_gene_trees);
				else
					if (UNROOTED || SIMPLE_UNROOTED)
						distance = rSPR_total_distance_unrooted(ST[j], current_gene_trees);
					else
						distance = rSPR_total_distance(ST[j], current_gene_trees);
		//		cout << distance <<  ": ";
		//		cout << ST[j]->str_subtree() << endl;
				if (distance < best_distance) {
					best_distance = distance;
					best_tree = j;
				}
			}
			super_tree = ST[best_tree];
		
			for(int i = 0; i < ST.size(); i++) {
				if (i != best_tree)
					ST[i]->delete_tree();
			}
			current_gene_trees.clear();
			leaf_num = 5;
		}
	
		cout << endl;
		cout << "Initial Supertree:  " << super_tree->str_subtree() << endl;
		if (TIMING)
			time = clock()/(double)CLOCKS_PER_SEC;
	
		// GREEDY ADDITION
		int x = 0;
		for(; label != labels.rend() &&
				NUM_LEAVES < 0 || leaf_num <= NUM_LEAVES; label++) {
			if (INITIAL_SUPER_TREE != "" &&
					super_tree->contains_leaf(label->second)) {
				leaf_num++;
				continue;
			}
			cout << "Adding leaf " << label->second;
			cout << "\t("<< leaf_num++ << "/" <<  labels.size() << ")";
			if (TIMING) {
				current_time = time;
				time = clock()/(double)CLOCKS_PER_SEC;
				current_time = time - current_time;
				cout << "\t" << current_time << "\t" << time;
			}
			cout << endl;
			vector<Node *> current_gene_trees = vector<Node *>();
			for(int i = 0; i < gene_trees.size(); i++) {
				if (gene_trees[i]->contains_leaf(label->second))
					current_gene_trees.push_back(gene_trees[i]);
			}
			cout << "gene_trees: " << gene_trees.size() << endl;
			cout << "current_gene_trees: " << current_gene_trees.size() << endl;
			if (SIMPLE_UNROOTED) {
				if (REROOT) {
					cout << "rerooting super_tree" << endl;
					// reroot the supertree based on the balanced accuracy of splits
					vector<Node *> descendants =
						super_tree->find_interior();
						//super_tree->find_descendants();
					Node *best_root = super_tree->lchild();
					double best_root_avg_acc = 0;
					if (APPROX_ROOTING || EXACT_ROOTING || RANDOM_ROOTING)
						best_root_avg_acc = -INT_MAX;
					int num_ties = 2;
					for(int j = 0; j < descendants.size(); j++) {
						double root_avg_acc = 0;
						int count = 0;
						super_tree->reroot(descendants[j]);
						super_tree->set_depth(0);
						super_tree->fix_depths();
						super_tree->preorder_number();
						if (APPROX_ROOTING) {
							root_avg_acc = -rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
							//root_avg_acc = -rSPR_total_distance_unrooted(super_tree, current_gene_trees);
						}
						else if (EXACT_ROOTING) {
							root_avg_acc = -rSPR_total_distance(super_tree, gene_trees);
						}
						else if (RANDOM_ROOTING) {
							root_avg_acc = 0;
						}
						else {
							int end = current_gene_trees.size();
							#pragma omp parallel for schedule(static) reduction(+: root_avg_acc) reduction(+: count)
							for(int i = 0; i < end; i++) {
								double acc;
								current_gene_trees[i]->preorder_number();
								acc = find_best_root_acc(super_tree, current_gene_trees[i]);
			//					cout <<  j << "\t" << i << "\t" << acc << endl;
								if (acc > -1) {
									root_avg_acc += acc;
									//root_avg_acc += (acc - root_avg_acc) / count;
									count++;
								}
							}
							if (count > 0)
								root_avg_acc /= count;
							int lsize = super_tree->lchild()->size_using_prenum();
							int rsize = super_tree->rchild()->size_using_prenum();
							int size = (lsize < rsize) ? lsize : rsize;
							root_avg_acc *= mylog2(size);
						}
						#pragma omp critical
						{
							if (root_avg_acc > best_root_avg_acc) {
								best_root = descendants[j];
								best_root_avg_acc = root_avg_acc;
							}
							else if (root_avg_acc == best_root_avg_acc) {
								int r = rand();
								if (r < RAND_MAX/ num_ties) {
									best_root = descendants[j];
									best_root_avg_acc = root_avg_acc;
									num_ties = 2;
								}
							}
						}
					}
					super_tree->reroot(best_root);
					super_tree->set_depth(0);
					super_tree->fix_depths();
					super_tree->preorder_number();
				}
	
				cout << "rerooting gene trees" << endl;
				// reroot the gene trees based on the balanced accuracy of splits
				super_tree->preorder_number();
				int end = current_gene_trees.size();
				#pragma omp parallel for
				for(int i = 0; i < end; i++) {
					current_gene_trees[i]->preorder_number();
					Node *new_root;
					if (EXACT_ROOTING)
						new_root = find_best_root_rspr(super_tree, current_gene_trees[i]);
					else if (RANDOM_ROOTING)
						new_root = find_random_root(super_tree, current_gene_trees[i]);
					else
						new_root = find_best_root(super_tree, current_gene_trees[i]);
					if (new_root != NULL) {
						current_gene_trees[i]->reroot(new_root);
						current_gene_trees[i]->set_depth(0);
						current_gene_trees[i]->fix_depths();
						current_gene_trees[i]->preorder_number();
					}
				}
			}
			Node *best_sibling;
			if (APPROX_SIBLINGS > 0) {
				cout << "finding approx best siblings" << endl;
				vector<Node *> *best_siblings = find_best_siblings(super_tree,
						current_gene_trees, label->second, APPROX_SIBLINGS);
				cout << "finding best sibling from " << best_siblings->size() << endl;
				best_sibling = find_best_sibling(super_tree,
						current_gene_trees, best_siblings, label->second);
				delete best_siblings;
			}
			else {
				cout << "finding best sibling" << endl;
				best_sibling = find_best_sibling(super_tree,
						current_gene_trees, label->second);
			}
	
	
			Node *node = best_sibling->expand_parent_edge(best_sibling);
	
			node->add_child(new Node(itos(label->second)));
			super_tree->numbers_to_labels(&reverse_label_map);
			cout << super_tree->str_subtree() << endl;
			super_tree->labels_to_numbers(&label_map, &reverse_label_map);
		}

	}

	cout << endl;
	super_tree->numbers_to_labels(&reverse_label_map);
	cout << "Initial Supertree: " <<  super_tree->str_subtree() << endl;
	super_tree->labels_to_numbers(&label_map, &reverse_label_map);
	Node *best_supertree = new Node(*super_tree);

	if (!LGT_ANALYSIS && !LGT_EVALUATION) {

			if (UNROOTED_MIN_APPROX)
				APPROX_ROOTING=true;
			if (REROOT_INITIAL) {
				cout << "rerooting super_tree" << endl;
				// reroot the supertree based on the balanced accuracy of splits
				vector<Node *> descendants =
					super_tree->find_interior();
					//super_tree->find_descendants();
				Node *best_root = super_tree->lchild();
				double best_root_avg_acc = 0;
				if (APPROX_ROOTING || EXACT_ROOTING || RANDOM_ROOTING)
					best_root_avg_acc = -INT_MAX;
				int num_ties = 2;
				int end = descendants.size();
				for(int j = 0; j < end; j++) {
					if (j > 0)
						cout << "\r \r";
					cout << j << "/" << end << flush;
					double root_avg_acc = 0;
					int count = 0;
					super_tree->reroot(descendants[j]);
					super_tree->set_depth(0);
					super_tree->fix_depths();
					super_tree->preorder_number();
					if (APPROX_ROOTING) {
						root_avg_acc = -rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
						//root_avg_acc = -rSPR_total_distance_unrooted(super_tree, current_gene_trees);
					}
					else if (EXACT_ROOTING) {
						root_avg_acc = -rSPR_total_distance(super_tree, gene_trees);
					}
					else if (RANDOM_ROOTING) {
						root_avg_acc = 0;
					}
					else {
						int end = gene_trees.size();
						#pragma omp parallel for schedule(static) reduction(+: root_avg_acc) reduction(+: count)
						for(int i = 0; i < end; i++) {
							gene_trees[i]->preorder_number();
							double acc;
							acc = 
								find_best_root_acc(super_tree, gene_trees[i]);
		//					cout <<  j << "\t" << i << "\t" << acc << endl;
							if (acc > -1) {
								root_avg_acc += acc;
								//root_avg_acc += (acc - root_avg_acc) / count;
								count++;
							}
						}
						if (count > 0)
							root_avg_acc /= count;
						int lsize = super_tree->lchild()->size_using_prenum();
						int rsize = super_tree->rchild()->size_using_prenum();
						int size = (lsize < rsize) ? lsize : rsize;
						root_avg_acc *= mylog2(size);
					}
					#pragma omp critical
					{
						if (root_avg_acc > best_root_avg_acc) {
							best_root = descendants[j];
							best_root_avg_acc = root_avg_acc;
						}
						else if (root_avg_acc == best_root_avg_acc) {
							int r = rand();
							if (r < RAND_MAX/ num_ties) {
								best_root = descendants[j];
								best_root_avg_acc = root_avg_acc;
								num_ties = 2;
							}
						}
					}
				}
				super_tree->reroot(best_root);
				super_tree->set_depth(0);
				super_tree->fix_depths();
				super_tree->preorder_number();
				super_tree->numbers_to_labels(&reverse_label_map);
				cout << "Rerooted Supertree: " <<  super_tree->str_subtree() << endl;
				super_tree->labels_to_numbers(&label_map, &reverse_label_map);
			}
	if (SIMPLE_UNROOTED) {
		cout << "rerooting gene trees" << endl;
		// reroot the gene trees based on the balanced accuracy of splits
		super_tree->preorder_number();
		int end = gene_trees.size();
		#pragma omp parallel for
		for(int i = 0; i < end; i++) {
			gene_trees[i]->preorder_number();
			Node *new_root;
			if (EXACT_ROOTING)
				new_root = find_best_root_rspr(super_tree, gene_trees[i]);
			else if (RANDOM_ROOTING)
				new_root = find_random_root(super_tree, gene_trees[i]);
			else
				new_root = find_best_root(super_tree, gene_trees[i]);
				//find_best_root(super_tree, gene_trees[i]);
			if (new_root != NULL) {
				gene_trees[i]->reroot(new_root);
				gene_trees[i]->set_depth(0);
				gene_trees[i]->fix_depths();
				gene_trees[i]->preorder_number();
			}
		}
	}
	super_tree->set_depth(0);
	super_tree->fix_depths();
	super_tree->preorder_number();

	if (APPROX)
		if (UNROOTED)
			best_distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
		else
			best_distance = rSPR_total_approx_distance(super_tree, gene_trees);
	else
		if (UNROOTED || RANDOM_ROOTING || (SIMPLE_UNROOTED && !SIMPLE_UNROOTED_FAST) )
		//if (UNROOTED || RANDOM_ROOTING)
			best_distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
		else
			best_distance = rSPR_total_distance(super_tree, gene_trees);
	cout << "Total Distance: " << best_distance << endl;
	if (RF_TIES) {
		int current_rf_distance = rf_total_distance(super_tree, gene_trees);
		cout << "Total RF Distance: " << current_rf_distance << endl;
		best_tie_distance = current_rf_distance;
	}
		if (TIMING) {
			current_time = time;
			time = clock()/(double)CLOCKS_PER_SEC;
			current_time = time - current_time;
			cout << "\t" << current_time << "\t" << time << endl;
		}

	}
		bool cleanup = false;
	//super_tree->numbers_to_labels(&reverse_label_map);
		if (FIND_SUPPORT) {
			get_support(super_tree, &gene_trees);
			super_tree->numbers_to_labels(&reverse_label_map);
			cout << super_tree->str_support_subtree() << endl;
			cleanup = true;
		}
		else if (FIND_BIPARTITION_SUPPORT) {
			cout << "finding bipartition support" << endl;
			super_tree->preorder_number();
			super_tree->edge_preorder_interval();
			for(int i = 0; i < gene_trees.size(); i++) {
				gene_trees[i]->preorder_number();
				gene_trees[i]->edge_preorder_interval();
			}
			get_bipartition_support(super_tree, &gene_trees,
					RELAXED_BIPARTITION_SUPPORT);
			super_tree->normalize_support();
			super_tree->numbers_to_labels(&reverse_label_map);
			cout << super_tree->str_support_subtree(true) << endl;
			cleanup = true;
		}
		else if (FIND_CLADE_TRANSFERS) {
			get_transfer_support(super_tree, &gene_trees);
			super_tree->numbers_to_labels(&reverse_label_map);
			cout << super_tree->str_support_subtree() << endl;
			cleanup = true;
		}

		if (VALID_TREES_ROOTED) {
			cout << "Rerooted Gene Trees: " <<  endl;
			for(int i = 0; i < gene_trees.size(); i++) {
				cout << gene_tree_names[i];
				gene_trees[i]->numbers_to_labels(&reverse_label_map);
				cout << gene_trees[i]->str_subtree() << endl;
				gene_trees[i]->labels_to_numbers(&label_map, &reverse_label_map);
			}
			cleanup = true;
		}

		if (LGT_ANALYSIS) {
			cout << "LGT Analysis" << endl;
			super_tree->preorder_number();
			super_tree->edge_preorder_interval();
			int num_nodes = super_tree->size();
			vector<vector<int> > transfer_counts =
				vector<vector<int> >(num_nodes, vector<int>(num_nodes, 0));
			for(int i = 0; i < gene_trees.size(); i++) {
				gene_trees[i]->preorder_number();
				gene_trees[i]->edge_preorder_interval();
			}
			add_transfers(&transfer_counts, super_tree, &gene_trees);
#ifdef DEBUG_LGT
			for(int i = 0; i < num_nodes; i++) {
				for(int j = 0; j < num_nodes; j++) {
					if (j > 0)
						cout << " ";
					cout << transfer_counts[i][j];
				}
				cout << endl;
			}
#endif

			if (LGT_GROUPS != "") {
				ifstream lgt_group_file;
				lgt_group_file.open(LGT_GROUPS.c_str());
				vector<int> pre_to_group = vector<int>(num_nodes, 0);
				vector<string> group_names = vector<string>();
				map<string, int> name_to_pre = map<string, int>();
				super_tree->numbers_to_labels(&reverse_label_map);
				super_tree->build_name_to_pre_map(&name_to_pre);
				super_tree->labels_to_numbers(&label_map, &reverse_label_map);
				group_names.push_back("Mixed");
				if (lgt_group_file.is_open()) {
					string line;
					bool new_group = true;
					int group_num = 0;
					while(lgt_group_file.good()) {
						getline(lgt_group_file, line);
						// finish group
						if (line == "")
							new_group = true;
						// new group
						else if (new_group) {
							group_names.push_back(line);
							group_num++;
							new_group = false;
						}
						// add to group
						else {
							map<string, int>::iterator i = name_to_pre.find(line);
							if (i != name_to_pre.end()) {
								int pre = i->second;
#ifdef DEBUG_LGT
								cout << pre << ": " << line << endl;
								cout << "group " << group_num << endl;
#endif
								pre_to_group[pre] = group_num;
							}
						}
					}
#ifdef DEBUG_LGT
					cout << endl;

					for(int i = 0; i < num_nodes; i++) {
						cout << pre_to_group[i] << ": " << super_tree->find_by_prenum(i)->str_subtree() << endl;
					}
#endif

					// add LCAs to groups
					add_lcas_to_groups(&pre_to_group, super_tree);

					// TODO: backfill polytomies
					// * look at children and grandchildren of a 0 group node
					// * in pre order
					// * match parent group (or this is the root)

#ifdef DEBUG_LGT
					for(int i = 0; i < num_nodes; i++) {
						cout << pre_to_group[i] << ": " << super_tree->find_by_prenum(i)->str_subtree() << endl;
					}
#endif


					cout << "INFERRED SPRS" << endl;
					if (LGT_CSV)
						cout << ",";
					for(int i = 0; i < group_names.size(); i++) {
						cout << group_names[i];
						if (LGT_CSV && (i + 1 < group_names.size())) {
							cout << ",";
						}
						else
							cout << endl;
					}

					int num_group_nodes = group_names.size();
					vector<vector<int> > group_transfer_counts =
							vector<vector<int> >(num_group_nodes,
							vector<int>(num_group_nodes, 0));

					for(int i = 0; i < num_nodes; i++) {
						for(int j = 0; j < num_nodes; j++) {
							group_transfer_counts[pre_to_group[i]][pre_to_group[j]]
								+= transfer_counts[i][j];
						}
					}
					if (!LGT_CSV)
						cout << endl;
					for(int i = 0; i < num_group_nodes; i++) {
						if (LGT_CSV)
							cout << group_names[i] << ",";
						for(int j = 0; j < num_group_nodes; j++) {
							cout << group_transfer_counts[i][j];
							if (j + 1 < group_names.size()) {
								if (LGT_CSV)
									cout << ",";
								else
									cout << " ";
							}
						}
						cout << endl;
					}
					cout << endl;
					// NOTE: these are inferred SPRs
					//       transfers are transposed
					cout << "INFERRED TRANSFERS" << endl;
					if (LGT_CSV)
						cout << ",";
					for(int i = 0; i < group_names.size(); i++) {
						cout << group_names[i];
						if (LGT_CSV && (i + 1 < group_names.size())) {
							cout << ",";
						}
						else cout << endl; }
					if (!LGT_CSV)
						cout << endl;
					for(int j = 0; j < num_group_nodes; j++) {
						if (LGT_CSV)
							cout << group_names[j] << ",";
						for(int i = 0; i < num_group_nodes; i++) {
							cout << group_transfer_counts[i][j];
							if (i + 1 < group_names.size()) {
								if (LGT_CSV)
									cout << ",";
								else
									cout << " ";
							}
						}
						cout << endl;
					}
					cout << endl;

					cout << "NONDIRECTIONAL TRANSFERS" << endl;
					if (LGT_CSV)
						cout << ",";
					for(int i = 0; i < group_names.size(); i++) {
						cout << group_names[i];
						if (LGT_CSV && (i + 1 < group_names.size())) {
							cout << ",";
						}
						else cout << endl; }
					if (!LGT_CSV)
						cout << endl;
					for(int j = 0; j < num_group_nodes; j++) {
						if (LGT_CSV)
							cout << group_names[j] << ",";
						for(int i = 0; i < num_group_nodes; i++) {
							cout << group_transfer_counts[i][j]
								+ (i == j ?  0 : group_transfer_counts[j][i]);
							if (i + 1 < group_names.size()) {
								if (LGT_CSV)
									cout << ",";
								else
									cout << " ";
							}
						}
						cout << endl;
					}
					cout << endl;
					lgt_group_file.close();
				}
			else
				LGT_GROUPS = "";
			}
			cleanup = true;
		}
		else if (LGT_EVALUATION) {
			cout << "LGT Evaluation" << endl;
			super_tree->preorder_number();
			super_tree->edge_preorder_interval();
			for(int i = 0; i < gene_trees.size(); i++) {
				gene_trees[i]->preorder_number();
				gene_trees[i]->edge_preorder_interval();
			}
			print_transfers(super_tree, &gene_trees, &gene_tree_names,
					&reverse_label_map);
			cleanup = true;
		}

		if (cleanup) {
			for(int i = 0; i < gene_trees.size(); i++) {
				gene_trees[i]->delete_tree();
			}
			super_tree->delete_tree();
			best_supertree->delete_tree();
			return 0;
		}


/*Joel: added these vectors*/
	vector< pair<Node*, int> > best_scores;
	vector< pair<Node* ,int> > scores;
	vector<pair <int, pair <int, int> > > stats;
	vector<pair <pair<Node*,Node*>, int> > approx_moves;
/*end*/
	int num_zeros=0;
	int edges_cut = 0;
	int current_distance = 0;
	if (NUM_ITERATIONS < 0)
		NUM_ITERATIONS=labels.size(); 

	// SUPERTREE IMPROVEMENT STEP
	for(int i = 0; i < NUM_ITERATIONS; i++) {
/*		if (TABOO_SEARCH) {
			cout << "TABOO " << taboo_trees.size() << endl;
			list<Node>::iterator t;
			for(t = taboo_trees.begin(); t != taboo_trees.end(); t++) {
				cout << (*t).str_subtree() << endl;
			}
			cout << endl;
		}
*/

			if (REROOT) {
				cout << "rerooting super_tree" << endl;
				// reroot the supertree based on the balanced accuracy of splits
				vector<Node *> descendants =
					super_tree->find_interior();
					//super_tree->find_descendants();
				Node *best_root = super_tree->lchild();
				double best_root_avg_acc = 0;
				if (APPROX_ROOTING || EXACT_ROOTING || RANDOM_ROOTING)
					best_root_avg_acc = -INT_MAX;
				int num_ties = 2;
				int end = descendants.size();
				for(int j = 0; j < end; j++) {
					if (j > 0)
						cout << "\r \r";
					cout << j << "/" << end - 1 << flush;
					double root_avg_acc = 0;
					int count = 0;
					super_tree->reroot(descendants[j]);
					super_tree->set_depth(0);
					super_tree->fix_depths();
					super_tree->preorder_number();
					if (APPROX_ROOTING) {
						root_avg_acc = -rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
						//root_avg_acc = -rSPR_total_distance_unrooted(super_tree, current_gene_trees);
					}
					else if (EXACT_ROOTING) {
						root_avg_acc = -rSPR_total_distance(super_tree, gene_trees);
					}
					else if (RANDOM_ROOTING) {
						root_avg_acc = 0;
					}
					else {
						int end = gene_trees.size();
						#pragma omp parallel for schedule(static) reduction(+: root_avg_acc) reduction(+: count)
						for(int i = 0; i < end; i++) {
							gene_trees[i]->preorder_number();
							double acc;
							acc = 
								find_best_root_acc(super_tree, gene_trees[i]);
		//					cout <<  j << "\t" << i << "\t" << acc << endl;
							if (acc > -1) {
								root_avg_acc += acc;
								//root_avg_acc += (acc - root_avg_acc) / count;
								count++;
							}
						}
						if (count > 0)
							root_avg_acc /= count;
						int lsize = super_tree->lchild()->size_using_prenum();
						int rsize = super_tree->rchild()->size_using_prenum();
						int size = (lsize < rsize) ? lsize : rsize;
						root_avg_acc *= mylog2(size);
					}
					#pragma omp critical
					{
						if (root_avg_acc > best_root_avg_acc) {
							best_root = descendants[j];
							best_root_avg_acc = root_avg_acc;
						}
						else if (root_avg_acc == best_root_avg_acc) {
							int r = rand();
							if (r < RAND_MAX/ num_ties) {
								best_root = descendants[j];
								best_root_avg_acc = root_avg_acc;
								num_ties = 2;
							}
						}
					}
				}
				cout << "\r \r";
				super_tree->reroot(best_root);
				super_tree->set_depth(0);
				super_tree->fix_depths();
				super_tree->preorder_number();
				super_tree->numbers_to_labels(&reverse_label_map);
				cout << "Rerooted Supertree: " <<  super_tree->str_subtree() << endl;
				super_tree->labels_to_numbers(&label_map, &reverse_label_map);
				if (TABOO_SEARCH && !is_taboo(taboo_trees, super_tree))
						taboo_trees.push_back(Node(*super_tree));
			}
		if (SIMPLE_UNROOTED && SIMPLE_UNROOTED_NUM > 0) {
			SIMPLE_UNROOTED_NUM--;
			cout << "rerooting gene trees" << endl;
			// reroot the gene trees based on the balanced accuracy of splits
			super_tree->preorder_number();
			int end = gene_trees.size();
			#pragma omp parallel for
			for(int i = 0; i < end; i++) {
				gene_trees[i]->preorder_number();
				Node *new_root;
				if (EXACT_ROOTING)
					new_root = find_best_root_rspr(super_tree, gene_trees[i]);
				else if (RANDOM_ROOTING)
					new_root = find_random_root(super_tree, gene_trees[i]);
				else
					new_root = find_best_root(super_tree, gene_trees[i]);
				if (new_root != NULL) {
					gene_trees[i]->reroot(new_root);
					gene_trees[i]->set_depth(0);
					gene_trees[i]->fix_depths();
					gene_trees[i]->preorder_number();
				}
			}
		}
		if (BIPARTITION_CLUSTER) {
			cout << "finding bipartition support" << endl;
			super_tree->preorder_number();
			super_tree->edge_preorder_interval();
			for(int i = 0; i < gene_trees.size(); i++) {
				gene_trees[i]->preorder_number();
				gene_trees[i]->edge_preorder_interval();
			}
			get_bipartition_support(super_tree, &gene_trees,
					RELAXED_BIPARTITION_SUPPORT);
			super_tree->normalize_support();
//			super_tree->unprotect_subtree();
//			super_tree->protect_supported_edges();
		}
		Node *best_subtree_root;
		Node *best_sibling;

		// SUPERTREE IMPROVEMENT METHOD

		/* one-tree-at-a-time method
		 * in: rooted super_tree, gene_trees, 
		 * opt: bipartition support, 
		 * do: compare each gene_tree to the supertree,
		 * 		propose an spr move
		 * 		accept with some criteria
		 * 			total SPR distance, sampled distance
		 * 			straight improvement or simulated annealing?
		 */
		if (ONE_TREE_AT_A_TIME) {
			// initial distance
			int min_distance;
			int min_tie_distance;
			int num_ties = 2;
			super_tree->set_depth(0);
			super_tree->fix_depths();
			super_tree->preorder_number();
			super_tree->edge_preorder_interval();
			vector<int> original_scores = vector<int>(gene_trees.size());
			vector<int> original_scores_temp = vector<int>(gene_trees.size());

			if (APPROX)
				if (UNROOTED)
					min_distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
				else
					min_distance = rSPR_total_approx_distance(super_tree, gene_trees);
			else
				if (UNROOTED || RANDOM_ROOTING || (SIMPLE_UNROOTED && !SIMPLE_UNROOTED_FAST) )
				//if (UNROOTED || RANDOM_ROOTING)
					min_distance = rSPR_total_distance_unrooted(super_tree, gene_trees, INT_MAX, &original_scores);
				else
					min_distance = rSPR_total_distance(super_tree, gene_trees, &original_scores);
			cout << "Total Distance: " << min_distance << endl;
			if (RF_TIES) {
				int current_rf_distance = rf_total_distance(super_tree, gene_trees);
				cout << "Total RF Distance: " << current_rf_distance << endl;
				min_tie_distance = current_rf_distance;
			}

			/*
			vector<pair<int, int> > sorted_gene_trees = vector<pair<int, int> >(gene_trees.size());
			for(int j = 0; j < gene_trees.size(); j++) {
				sorted_gene_trees[j] = make_pair(gene_trees[j]->size(), j);
			}
			sort(sorted_gene_trees.begin(), sorted_gene_trees.end(), pair_comparator);

			int end = (int)sqrt(gene_trees.size());

			vector<int> gene_tree_subset = vector<int>(end);
			for(int j = 0; j < end; j++) {
				gene_tree_subset[j] = sorted_gene_trees[j].second;
//				cout << sorted_gene_trees[j].first << ":" << gene_tree_subset[j] << ", ";
			}
//			cout << endl;

*/



			// for each tree
			int end = gene_trees.size();
			for(int k = 0; k < end; k++) {
//				int j = gene_tree_subset[k];
				int j = k;
//			for(int j = 0; j < gene_trees.size(); j++) {
// }
//			if (gene_trees[j]->size() < 10)
//				continue;
				if (k != 0)
					cout << "\r \r" << flush;
				cout << i+1 << "/" << NUM_ITERATIONS << "\t" << k+1 << "/" << end << flush;
				// compute an MAF
				Forest *MAF1 = NULL;
				Forest *MAF2 = NULL;
				Forest F1 = Forest(super_tree);
				Forest F2 = Forest(gene_trees[j]);
				if (!sync_twins(&F1,&F2)) {
					continue;
				}
				int distance = rSPR_branch_and_bound_simple_clustering(F1.get_component(0), F2.get_component(0), &MAF1, &MAF2);
				if (distance > 0) {
					expand_contracted_nodes(MAF1);
					expand_contracted_nodes(MAF2);
					#ifdef DEBUG_ONE_TREE
						cout << j << ": " << distance << endl;
						cout << "\tT1: "; F1.print_components();
						cout << "\tT2: "; F2.print_components();
						cout << "\tF1: "; MAF1->print_components_with_edge_pre_interval();
						cout << "\tF2: "; MAF2->print_components_with_edge_pre_interval();
					#endif
					sync_af_twins(MAF1, MAF2);
	
					// propose transfers for each component
					vector<vector<Node *> > transfers =
						vector<vector<Node *> >();
					for(int k = 0; k < MAF2->num_components(); k++) {
						if (k == 0 && !MAF2->contains_rho()) {
							continue;
						}
						Node *F1_source, *F1_target;
						if (!map_transfer(MAF2->get_component(k), &F1, MAF2,
								&F1_source, &F1_target)) {
							continue;
						}
						// check transfer validity
						if (F1_target->get_preorder_number() >=
								F1_source->get_edge_pre_start()
								&& F1_target->get_preorder_number() <=
								F1_source->get_edge_pre_end() ) {
							continue;
						}
					if (!supported_spr(super_tree->find_by_prenum(F1_source->get_preorder_number()),super_tree->find_by_prenum(F1_target->get_preorder_number()))) {
						continue;
					}

						
						// add transfer to list
						vector<Node *> transfer = vector<Node *>(2);
						transfer[0] = F1_source;
						transfer[1] = F1_target;
						transfers.push_back(transfer);

					}
					// any valid tranfers?
					if (transfers.empty()) {
						delete MAF1;
						delete MAF2;
						continue;
					}
					// pick a random transfer
					int chosen_num = rand() % transfers.size();
					vector<Node *> chosen_transfer = transfers[chosen_num];
					Node *F1_source = super_tree->find_by_prenum(chosen_transfer[0]->get_preorder_number());
					Node *F1_target = super_tree->find_by_prenum(chosen_transfer[1]->get_preorder_number());

//					if (!supported_spr(F1_target, F1_source)) {
//						continue;
//					}

					#ifdef DEBUG_ONE_TREE
						 cout << endl;
							cout << "SPR Move:" << endl;
							super_tree->numbers_to_labels(&reverse_label_map);
							cout << "Previous Super Tree: "
							<< super_tree->str_subtree() << endl;
							cout << "Subtree: " << F1_source->str_subtree() << endl;
							cout << "New Sibling: " << F1_target->str_subtree() << endl;
							super_tree->labels_to_numbers(&label_map, &reverse_label_map);	
					#endif
					
					//		cout << endl;
					//		cout << "Previous Super Tree: "
					//		<< super_tree->str_support_subtree(true) << endl;
					
					// apply the transfer
					Node *old_super_tree = new Node(*super_tree);
					int which_sibling = 0;
					Node *undo = F1_source->spr(F1_target, which_sibling);
					super_tree->set_depth(0);
					super_tree->fix_depths();
					super_tree->preorder_number();
					super_tree->edge_preorder_interval();

					//		cout << "Proposed Super Tree: "
					//		<< super_tree->str_support_subtree(true) << endl;
//					/*
					#ifdef DEBUG_ONE_TREE
							super_tree->numbers_to_labels(&reverse_label_map);
							cout << "Proposed Super Tree: " << super_tree->str_subtree() << endl;
							super_tree->labels_to_numbers(&label_map, &reverse_label_map);
					#endif
					// test the new supertree
//					*/	
					int distance;
					if (APPROX) {
						if (UNROOTED)
							distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
						else
							distance = rSPR_total_approx_distance(super_tree, gene_trees);
					}
					else {
						if (UNROOTED)
							distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
						else {
							if (USE_PRECOMPUTED_DISTANCES) {
								distance = rSPR_total_distance_precomputed(super_tree, gene_trees, &original_scores, &original_scores_temp, old_super_tree);
							}
							else {
								distance = rSPR_total_distance(super_tree, gene_trees);
							}
						}
					}
					old_super_tree->delete_tree();
							cout << "\t" << distance << "\t" << min_distance << flush;


					//		cout << "After SPR Distance Super Tree: "
					//		<< super_tree->str_support_subtree(true) << endl;
					bool good_move = false;
					if (distance < min_distance) {
						if (!TABOO_SEARCH || !is_taboo(taboo_trees, super_tree)) {
							min_distance = distance;
							if (RF_TIES)
								min_tie_distance = rf_total_distance(super_tree, gene_trees);
							num_ties = 2;
							good_move = true;
						}
					}
					else if (distance == min_distance) {
						bool check_tie = true;
						bool tie = false;
						if (RF_TIES) {
							int rf_distance = rf_total_distance(super_tree, gene_trees);
							if (rf_distance < min_tie_distance) {
								if (!TABOO_SEARCH || !is_taboo(taboo_trees, super_tree)) {
									min_tie_distance = rf_distance;
									num_ties = 2;
									good_move = true;
								}
								check_tie = false;
							}
							else if (rf_distance == min_tie_distance)
								tie = true;
							else
								check_tie = false;
						}
						if (check_tie) {
							int r = rand();
							if (r < RAND_MAX/num_ties) {
								if (!TABOO_SEARCH || !is_taboo(taboo_trees, super_tree)) {
									good_move = true;
								}
							}
							num_ties++;
						}
					}
					// rollback if worse (or simulated annealing?)
					if (good_move) {
						get_bipartition_support(super_tree, &gene_trees,
						RELAXED_BIPARTITION_SUPPORT);
						super_tree->normalize_support();
						cout << endl;
						super_tree->numbers_to_labels(&reverse_label_map);
						cout << super_tree->str_subtree() << endl;
						super_tree->labels_to_numbers(&label_map, &reverse_label_map);
						for(int i = 0; i < gene_trees.size(); i++) {
							original_scores[i] = original_scores_temp[i];
						}
					}
					else {
						// restore the previous tree
						F1_source->spr(undo, which_sibling);
						super_tree->set_depth(0);
						super_tree->fix_depths();
						super_tree->preorder_number();
						super_tree->edge_preorder_interval();
						//		cout << "Reverted Super Tree: "
						//		<< super_tree->str_support_subtree(true) << endl;
					}
			
			
					//		cout << "Reverted Super Tree: "
					//		<< super_tree->str_subtree() << endl;
				}
				if (MAF1 != NULL) {
					delete MAF1;
				}
				if (MAF2 != NULL) {
					delete MAF2;
				}
			}

			cout << endl;

		}
/*Joel: Limit starting point*/
		else if(!R_LIMIT && S_LIMIT){
			scores = vector< pair<Node *,int> >();
			if(current_distance !=0)
				find_best_distance(super_tree, super_tree, gene_trees, scores, num_zeros, current_distance);
			else
				find_best_distance(super_tree, super_tree, gene_trees, scores, num_zeros, best_distance);
			edges_cut = scores.size();
			best_scores = get_best_scores(scores);
//		for(vector<pair<Node *, int> >::const_iterator i = best_scores.begin(); i != best_scores.end(); i++){
//			cout << "NODE:\n";			
//			cout << i->first->str_subtree() << endl;
//			cout << "DISTAMCE:\n";
//			cout << i->second << endl;
//		}
//		scores.clear();
//		best_scores.clear();	
		}
/*Combining limiting starting point and SPR radius*/
		else if (R_LIMIT && S_LIMIT){
			int r = find_r();
			int min_distance = INT_MAX;
			int num_ties = 0;
			for(vector<pair<Node *, int> >::const_iterator i = best_scores.begin(); i != best_scores.end(); i++){
				Node *source = i->first;
				if(!R_CONTROL && R_RAND)
					r = find_r(R_PROB);
				if(source->parent()->lchild() == source)
					find_best_spr_r_helper(source, source->parent(), super_tree, gene_trees, best_subtree_root, best_sibling, min_distance, num_ties, r+1, 1);
				else if(source->parent()->rchild() == source)
					find_best_spr_r_helper(source, source->parent(), super_tree, gene_trees, best_subtree_root, best_sibling, min_distance, num_ties, r+1, 2);
				else
					find_best_spr_r_helper(source, source->parent(), super_tree, gene_trees, best_subtree_root, best_sibling, min_distance, num_ties, r+1, 3);
			}

		}
/*Limiting SPR Radius*/		
		else if(R_LIMIT && !S_LIMIT){
			int r = find_r();
			vector<int> *original_scores = NULL;
			if (!APPROX) {
				original_scores = new vector<int>(gene_trees.size(), 0);
				int distance;
				if (UNROOTED_MIN_APPROX) {
					distance = rSPR_total_distance_unrooted(super_tree, gene_trees, INT_MAX, original_scores);
				}
				else {
					distance = rSPR_total_distance(super_tree, gene_trees, original_scores);
				}
				if (SIMPLE_UNROOTED) {
					cout << "Rerooted Distance: " << distance << endl;
						best_rooted_distance = distance;
//					if (distance < best_rooted_distance) {
						//best_supertree->delete_tree();
						//best_supertree = new Node(*super_tree);
//					}
				}
			}
			NUM_SOURCE = super_tree->size();
			C_SOURCE = 1;
			find_best_spr_r(super_tree, gene_trees, best_subtree_root, best_sibling,r, original_scores);
			if ((!UNROOTED || UNROOTED_MIN_APPROX) && !APPROX) {
				delete original_scores;
				cout << endl;
			}
		}
/*Limiting Starting point*/
		else if(S_LIMIT && !R_LIMIT){
			int min_distance = INT_MAX;
			int num_ties = 0;
			for(vector<pair<Node *, int> >::const_iterator i = best_scores.begin(); i!= best_scores.end(); i++){
				Node *source = i->first;
				find_best_spr_helper(source, super_tree, gene_trees, best_subtree_root, best_sibling, min_distance, num_ties);
			}

		}
/*Greedy search*/
		else if(GREEDY){
			approx_moves = vector<pair <pair<Node*,Node*>, int> >();
			find_best_spr(super_tree, gene_trees, best_subtree_root, best_sibling, approx_moves);
			std::sort(approx_moves.begin(), approx_moves.end(), sort_approx_moves);
			bool check = false;
			int current_distance;
			for(int i = 0; i < approx_moves.size() && !check ; i++){
				int which_sibling = 0;
				Node *undo = approx_moves[i].first.first->spr(approx_moves[i].first.second, which_sibling);
				super_tree->set_depth(0);
				super_tree->fix_depths();
				super_tree->preorder_number();
	

				if (APPROX)
					if (UNROOTED)
						current_distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
					else
						current_distance = rSPR_total_approx_distance(super_tree, gene_trees);
				else
					if (UNROOTED)
						current_distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
					else
						current_distance = rSPR_total_distance(super_tree, gene_trees);
				if(current_distance < best_distance)
					check = true;

				else {
					approx_moves[i].first.first->spr(undo, which_sibling);
					super_tree->set_depth(0);
					super_tree->fix_depths();
					super_tree->preorder_number();
				}
			}
			super_tree->numbers_to_labels(&reverse_label_map);
			cout << "Current Supertree: " << super_tree->str_subtree() << endl;
			super_tree->labels_to_numbers(&label_map, &reverse_label_map);
			if (current_distance < best_distance) {
				best_supertree->delete_tree();
				best_supertree = new Node(*super_tree);
				best_distance = current_distance;
			}
			cout << "Total Distance: " << current_distance << endl;
			approx_moves.clear();
		}
/*Refined Greedy Search*/
		else if(GREEDY_REFINED){
			scores = vector< pair<Node *,int> >();
			if(current_distance !=0){
				find_best_distance(super_tree, super_tree, gene_trees, scores, num_zeros, current_distance);
			}
			else{
				find_best_distance(super_tree, super_tree, gene_trees, scores, num_zeros, best_distance);
			}
			edges_cut = scores.size();
			best_scores = get_best_scores(scores);
			int r = find_r();
			int min_distance = INT_MAX;
			int num_ties = 0;
			approx_moves = vector<pair <pair<Node*,Node*>, int> >();
			for(vector<pair<Node *, int> >::const_iterator i = best_scores.begin(); i != best_scores.end(); i++){
				Node *source = i->first;
				if(!R_CONTROL && R_RAND)
					r = find_r(R_PROB);
				if(source->parent()->lchild() == source)
					find_best_spr_r_helper(source, source->parent(), super_tree, gene_trees, best_subtree_root, best_sibling, min_distance, num_ties, r+1, 1, approx_moves);
				else if(source->parent()->rchild() == source)
					find_best_spr_r_helper(source, source->parent(), super_tree, gene_trees, best_subtree_root, best_sibling, min_distance, num_ties, r+1, 2, approx_moves);
				else
					find_best_spr_r_helper(source, source->parent(), super_tree, gene_trees, best_subtree_root, best_sibling, min_distance, num_ties, r+1, 3, approx_moves);
			}
			std::sort(approx_moves.begin(), approx_moves.end(), sort_approx_moves);
			bool check = false;
			int current_distance;
			for(int i = 0; i < approx_moves.size() && !check ; i++){
				int which_sibling = 0;
				Node *undo = approx_moves[i].first.first->spr(approx_moves[i].first.second, which_sibling);
				super_tree->set_depth(0);
				super_tree->fix_depths();
				super_tree->preorder_number();
	

				if (APPROX)
					if (UNROOTED)
						current_distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
					else
						current_distance = rSPR_total_approx_distance(super_tree, gene_trees);
				else
					if (UNROOTED)
						current_distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
					else
						current_distance = rSPR_total_distance(super_tree, gene_trees);
				if(current_distance < best_distance)
					check = true;

				else {
					approx_moves[i].first.first->spr(undo, which_sibling);
					super_tree->set_depth(0);
					super_tree->fix_depths();
					super_tree->preorder_number();
				}
			}
			super_tree->numbers_to_labels(&reverse_label_map);
			cout << "Current Supertree: " << super_tree->str_subtree() << endl;
			super_tree->labels_to_numbers(&label_map, &reverse_label_map);
			if (current_distance < best_distance) {
				best_supertree->delete_tree();
				best_supertree = new Node(*super_tree);
				best_distance = current_distance;
			}
			cout << "Total Distance: " << current_distance << endl;
			approx_moves.clear();
		}
/*Approx Distance vs Total Distance stats*/
		else if(D_STATS){
			stats = vector<pair <int, pair <int, int> > > ();
			find_best_spr(super_tree, gene_trees, best_subtree_root, best_sibling, stats);
		}
/*Default strategy*/	
		else {
			find_best_spr(super_tree, gene_trees, best_subtree_root, best_sibling);
		}
		if(!GREEDY && !GREEDY_REFINED) {
			if (!ONE_TREE_AT_A_TIME){
				best_subtree_root->spr(best_sibling);
				super_tree->set_depth(0);
				super_tree->fix_depths();
				super_tree->preorder_number();
			}
			if (TABOO_SEARCH && !is_taboo(taboo_trees, super_tree))
				taboo_trees.push_back(Node(*super_tree));
			super_tree->numbers_to_labels(&reverse_label_map);
			cout << "Current Supertree: " <<  super_tree->str_subtree() << endl;
			super_tree->labels_to_numbers(&label_map, &reverse_label_map);
//			super_tree->labels_to_numbers(&label_map, &reverse_label_map);
			int current_distance;
			if (APPROX)
				if (UNROOTED)
					current_distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
				else
					current_distance = rSPR_total_approx_distance(super_tree, gene_trees);
			else
				if (UNROOTED || RANDOM_ROOTING || (SIMPLE_UNROOTED && !SIMPLE_UNROOTED_FAST))
					current_distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
				else
					current_distance = rSPR_total_distance(super_tree, gene_trees);

			int current_tie_distance = best_tie_distance;
			if (RF_TIES)
				current_tie_distance = rf_total_distance(super_tree, gene_trees);
		
			if (current_distance < best_distance ||
					(RF_TIES && current_distance == best_distance
					 && current_tie_distance < best_tie_distance)) {
				best_supertree->delete_tree();
				best_supertree = new Node(*super_tree);
				best_distance = current_distance;
				best_tie_distance = current_tie_distance;
				if (R_VARIABLE && R_LIMIT && R_DISTANCE > 1) {
					R_DISTANCE--;
					cout << "R_DISTANCE=" << R_DISTANCE << endl;
				}
			}
			else if (R_VARIABLE && R_LIMIT) {
				if (SIMPLE_UNROOTED) {
					int current_rooted_distance =
						rSPR_total_distance(super_tree, gene_trees);
					if (current_rooted_distance < best_rooted_distance) {
						R_DISTANCE--;
						cout << "R_DISTANCE=" << R_DISTANCE << endl;
					}
					else {
						super_tree->delete_tree();
						super_tree = new Node(*best_supertree);
						R_DISTANCE++;
						cout << "R_DISTANCE=" << R_DISTANCE << endl;
					}
				}
			}
			cout << "Total Distance: " << current_distance << endl;
			if (RF_TIES)
				cout << "Total RF Distance: " << current_tie_distance << endl;
		}
/*Approx Distance vs Total Distance stats*/
		if(D_STATS){
			cout << "STATISTICS: " << endl;
			for(vector<pair<int, pair<int,int> > >::const_iterator i = stats.begin(); i != stats.end(); i++){
				cout << i->first << ":" << i->second.first << ":" << i->second.second << endl;
			}
			stats.clear();
		}
		if (TIMING) {
			current_time = time;
			time = clock()/(double)CLOCKS_PER_SEC;
			current_time = time - current_time;
			cout << "\t" << current_time << "\t" << time << endl;
		}

		scores.clear();
		best_scores.clear();
	}
/*Stats of starting location*/
	if(S_STATS){
		cout << "RATIO IS: " << (double)num_zeros << " over " << (double)(NUM_ITERATIONS*edges_cut) << endl;
	}

	super_tree->delete_tree();
	super_tree=best_supertree;
	super_tree->numbers_to_labels(&reverse_label_map);
	cout << "Final Supertree: " <<  super_tree->str_subtree() << endl;
	cout << "Final Distance: " << best_distance << endl;
	if (RF_TIES) {
		cout << "Final RF Distance: " << best_tie_distance << endl;
	}


	// cleanup
	for(int i = 0; i < gene_trees.size(); i++) {
		gene_trees[i]->delete_tree();
	}
	super_tree->delete_tree();

	return 0;

}


Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees, int label) {
	Node *best_sibling;
	Node *new_leaf = new Node(itos(label));
	int min_distance = INT_MAX;
	int min_tie_distance = INT_MAX;
	int num_ties = 0;
	find_best_sibling_helper(super_tree, new_leaf, super_tree, gene_trees,
			min_distance, min_tie_distance, num_ties, &best_sibling);
	delete new_leaf;
	return best_sibling;
}

void find_best_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &min_tie_distance,
		int &num_ties,
		Node **best_sibling) {

	if (n->lchild() != NULL) {
		find_best_sibling_helper(n->lchild(), new_leaf, super_tree,
				gene_trees, min_distance, min_tie_distance, num_ties, best_sibling);
	}
	if (n->rchild() != NULL) {
		find_best_sibling_helper(n->rchild(), new_leaf, super_tree,
				gene_trees, min_distance, min_tie_distance, num_ties, best_sibling);
	}
	test_sibling_helper(n, new_leaf, super_tree,
			gene_trees, min_distance, min_tie_distance, num_ties, best_sibling);
}

Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees,
		vector<Node *> *best_siblings, int label) {
	Node *best_sibling;
	Node *new_leaf = new Node(itos(label));
	int min_distance = INT_MAX;
	int min_tie_distance = INT_MAX;
	int num_ties = 0;
	for(int i = 0; i < best_siblings->size(); i++) {
		test_sibling_helper((*best_siblings)[i], new_leaf, super_tree, gene_trees,
				min_distance, min_tie_distance, num_ties, &best_sibling);
	}
	delete new_leaf;
	return best_sibling;
}

void test_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &min_tie_distance,
		int &num_ties, Node **best_sibling) {

	// TODO: modify this to keep a single intermediate node rather
	//than recreating and destroying it
//	cout << "Previous Super Tree" << super_tree->str_subtree() << endl;
	int status = -1;
	if (n->parent() != NULL)
		if (n->parent()->lchild() == n)
			status = 1;
		else
			status = 2;
	Node *new_node = n->expand_parent_edge(n);
	new_node->add_child(new_leaf);
//	cout << "New Super Tree" << super_tree->str_subtree() << endl;

	super_tree->set_depth(0);
	super_tree->fix_depths();
	super_tree->preorder_number();
	int distance;
	if (RANDOM_TREE)
		distance = 0;
	else {
		if (APPROX) {
			if (UNROOTED)
				distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
			else
				distance = rSPR_total_approx_distance(super_tree, gene_trees);
		}
		else {
			if (UNROOTED)
				distance = rSPR_total_distance_unrooted(super_tree, gene_trees,
						min_distance);
			else
				distance = rSPR_total_distance(super_tree, gene_trees, min_distance);
		}
	}
	if (distance < min_distance) {
		min_distance = distance;
		if (RF_TIES)
			min_tie_distance = rf_total_distance(super_tree, gene_trees);
		*best_sibling = n;
		num_ties = 2;
	}
	else if (distance == min_distance) {
		bool check_tie = true;
		bool tie = false;
		if (RF_TIES) {
			int rf_distance = rf_total_distance(super_tree, gene_trees);
			if (rf_distance < min_tie_distance) {
				min_tie_distance = rf_distance;
				*best_sibling = n;
				num_ties = 2;
				check_tie = false;
			}
			else if (rf_distance == min_tie_distance)
				tie = true;
			else
				check_tie = false;
		}
		if (check_tie) {
			int r = rand();
			if (r < RAND_MAX/num_ties) {
				*best_sibling = n;
			}
			num_ties++;
		}
	}
//	cout << "distance: " << distance;
//	cout << endl;

//	cout << super_tree->str_subtree() << endl;
	new_leaf->cut_parent();
//	cout << super_tree->str_subtree() << endl;
	new_node = new_node->undo_expand_parent_edge();
	delete new_node;

	if (n->parent() != NULL)
		if ((status == 1 && n->parent()->lchild() != n)
				|| (status == 2 && n->parent()->rchild() != n)) {
			Node *rc = n->parent()->lchild();
			rc->cut_parent();
			n->parent()->add_child(rc);
		}
	super_tree->set_depth(0);
	super_tree->fix_depths();
	super_tree->preorder_number();
//	cout << "Reverted: " << super_tree->str_subtree() << endl;

}

vector<Node *> *find_best_siblings(Node *super_tree, vector<Node *> &gene_trees, int label, int num_siblings) {
	Node *new_leaf = new Node(itos(label));
	int min_distance = INT_MAX;
	int min_tie_distance = INT_MAX;
	int num_ties = 0;
	multimap<int, Node*> best_siblings = multimap<int, Node*>();
	find_best_siblings_helper(super_tree, new_leaf, super_tree, gene_trees,
			min_distance, min_tie_distance, num_ties, &best_siblings, num_siblings);
	vector<Node *> *bs = new vector<Node *>();
	int end = best_siblings.size();
	for(int i = 0; i < end; i++) {
		multimap<int,  Node *>::iterator p =
			best_siblings.begin();
//		cout << i << ": " << p->first << endl;
		bs->push_back(p->second);
		best_siblings.erase(p);
	}
	delete new_leaf;
	return bs;
}

void find_best_siblings_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &min_tie_distance,
		int &num_ties, multimap<int, Node*> *best_siblings, int num_siblings) {

/*	vector<Node *> descendants = N->find_descendants();
	descendants.push_back(N);
	int end = descendants.size();
	for(int i = 0; i < end; i++) {
		Node *n = descendants[i];
		*/
	if (n->lchild() != NULL) {
		find_best_siblings_helper(n->lchild(), new_leaf, super_tree,
				gene_trees, min_distance, min_tie_distance, num_ties,
				best_siblings, num_siblings);
	}
	if (n->rchild() != NULL) {
		find_best_siblings_helper(n->rchild(), new_leaf, super_tree,
				gene_trees, min_distance, min_tie_distance, num_ties,
				best_siblings, num_siblings);
	}

	// TODO: modify this to keep a single intermediate node rather
	//than recreating and destroying it
//	cout << "Previous Super Tree" << super_tree->str_subtree() << endl;
	int status = -1;
	if (n->parent() != NULL)
		if (n->parent()->lchild() == n)
			status = 1;
		else
			status = 2;
	Node *new_node = n->expand_parent_edge(n);
	new_node->add_child(new_leaf);
//	cout << "New Super Tree" << super_tree->str_subtree() << endl;

	super_tree->set_depth(0);
	super_tree->fix_depths();
	super_tree->preorder_number();
	int distance;
	distance = rSPR_total_approx_distance(super_tree, gene_trees,
			min_distance);
{
	if (distance < min_distance) {
		best_siblings->insert(make_pair(distance, n));
		if (best_siblings->size() > num_siblings) {
			best_siblings->erase(--(best_siblings->end()));
			min_distance = (--(best_siblings->end()))->first;
		}
		num_ties = 2;
	}
	else if (distance == min_distance) {
		int r = rand();
		if (r < RAND_MAX/num_ties) {
			if (best_siblings->size() > num_siblings) {
				best_siblings->erase(--(best_siblings->end()));
				best_siblings->insert(make_pair(distance, n));
			}
		}
		num_ties++;
	}
}
//	cout << "distance: " << distance;
//	cout << endl;

//	cout << super_tree->str_subtree() << endl;
	new_leaf->cut_parent();
//	cout << super_tree->str_subtree() << endl;
	new_node = new_node->undo_expand_parent_edge();
	delete new_node;

	if (n->parent() != NULL)
		if ((status == 1 && n->parent()->lchild() != n)
				|| (status == 2 && n->parent()->rchild() != n)) {
			Node *rc = n->parent()->lchild();
			rc->cut_parent();
			n->parent()->add_child(rc);
		}
	super_tree->set_depth(0);
	super_tree->fix_depths();
	super_tree->preorder_number();
//	cout << "Reverted: " << super_tree->str_subtree() << endl;
//	}

}


void find_best_spr(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling) {
	int min_distance = INT_MAX;
	int num_ties = 0;
	find_best_spr_helper(super_tree, super_tree, gene_trees,
			best_spr_move, best_sibling, min_distance, num_ties);
}

void find_best_spr_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties) {

	if (n->lchild() != NULL) {
		find_best_spr_helper(n->lchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties);
	}
	if (n->rchild() != NULL) {
		find_best_spr_helper(n->rchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties);
	}

		find_best_spr_helper(n, super_tree, super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties);

}

void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties) {
	// do not consider invalid spr moves to n's subtree
	if (new_sibling == n)
		return;

	// recurse
	if (new_sibling->lchild() != NULL) {
		find_best_spr_helper(n, new_sibling->lchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties);
	}
	if (new_sibling->rchild() != NULL) {
		find_best_spr_helper(n, new_sibling->rchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties);
	}
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "s: " << new_sibling->str_subtree() << endl;

		//cout << "foo1" << endl;
	if (n->parent() != NULL &&
			(new_sibling == n->parent())) {/* ||
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent())) */
//		cout << "rule 1" << endl;
		return;
	}
	if (n->parent() != NULL &&
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent()) {
//		cout << "rule 2" << endl;
		return;
	}
	if (new_sibling == n->get_sibling()) {
//		cout << "rule 3" << endl;
		return;
	}
//		cout << "foo3" << endl;
//	cout << "foo4" << endl;
	if (n != super_tree && new_sibling != n) {
		Node *old_sibling = n->get_sibling();
		//if (new_sibling != old_sibling)


/*
		cout << "SPR Move:" << endl;
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Previous Super Tree: "
		<< super_tree->str_subtree() << endl;
		cout << "Subtree: " << n->str_subtree() << endl;
		cout << "New Sibling: " << new_sibling->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/
	

		int which_sibling = 0;
		Node *undo = n->spr(new_sibling, which_sibling);


		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();

/*
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Proposed Super Tree: " << super_tree->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/


		int distance;
		if (APPROX) {
			if (UNROOTED)
				distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
			else
				distance = rSPR_total_approx_distance(super_tree, gene_trees);
		}
		else {
			if (UNROOTED)
				distance = rSPR_total_distance_unrooted(super_tree, gene_trees,
						min_distance);
			else
				distance = rSPR_total_distance(super_tree, gene_trees, min_distance);
		}
//		cout << "\t" << distance << endl;

		if (distance < min_distance) {
			min_distance = distance;
			best_spr_move = n;
			best_sibling = new_sibling;
			num_ties = 2;
		}
		else if (distance == min_distance) {
			int r = rand();
			if (r < RAND_MAX/num_ties) {
				min_distance = distance;
				best_spr_move = n;
				best_sibling = new_sibling;
			}
			num_ties++;
		}
		// restore the previous tree


		n->spr(undo, which_sibling);

		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
	}

}


void get_support(Node *super_tree, vector<Node *> *gene_trees) {
	get_support(super_tree, super_tree, gene_trees);
}

void get_support(Node *n, Node *super_tree, vector<Node *> *gene_trees) {
	if (!n->is_leaf() ){
		list<Node *>::iterator c;
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			get_support(*c, super_tree, gene_trees);
		}

		if (n != super_tree) {
	
			int count = 0;
			double root_avg_acc = 0;
			Node *temp_tree = new Node(*super_tree);
			temp_tree->reroot(temp_tree->find_by_prenum(n->get_preorder_number()));
			temp_tree->set_depth(0);
			temp_tree->fix_depths();
			temp_tree->preorder_number();
			int end = gene_trees->size();
			#pragma omp parallel for schedule(static) reduction(+: root_avg_acc) reduction(+: count)
			for(int i = 0; i < end; i++) {
				(*gene_trees)[i]->preorder_number();
				double acc = find_best_root_acc(temp_tree, (*gene_trees)[i]);
				if (acc > -1) {
					root_avg_acc += acc;
					count++;
				}
			}
			if (count > 0)
				root_avg_acc /= count;
			root_avg_acc /= 2.0;
			n->set_support(root_avg_acc);
			temp_tree->delete_tree();
		}
	}
}

void get_transfer_support(Node *super_tree, vector<Node *> *gene_trees) {
	vector<Node *> descendants = super_tree->find_descendants();
	descendants.push_back(super_tree);
	int end = descendants.size();
	#pragma omp parallel for
	for(int i = 0; i < end; i++) {
		get_transfer_support(descendants[i], super_tree, gene_trees);
	}
}

void get_transfer_support(Node *n, Node *super_tree, vector<Node *> *gene_trees) {
	if (!n->is_leaf() ){
		Node *temp_tree = new Node(*super_tree);
		Node *temp_n = temp_tree->find_by_prenum(n->get_preorder_number());
		temp_tree->disallow_siblings_subtree();
		temp_n->allow_siblings_subtree();
		int distance = rSPR_total_distance(temp_tree, *gene_trees);
		n->set_support(distance);
		temp_tree->delete_tree();
	}
}

void get_bipartition_support(Node *super_tree, vector<Node *> *gene_trees,
		enum RELAXATION relaxed) {
	vector<Node *> descendants = super_tree->find_descendants();
	descendants.push_back(super_tree);
	for(int i = 0; i < descendants.size(); i++) {
		descendants[i]->set_support(0);
		descendants[i]->set_support_normalization(0);
	}
	int end = gene_trees->size();
	#pragma omp parallel for
	for(int i = 0; i < end; i++) {
		modify_bipartition_support(super_tree, (*gene_trees)[i], relaxed);
	}
}

bool is_taboo(list<Node> taboo_trees, Node *super_tree) {
	list<Node>::iterator i;
	Forest F1 = Forest(super_tree);
	for(i = taboo_trees.begin(); i != taboo_trees.end(); i++) {
		Forest F2 = Forest(&(*i));
		if (rSPR_branch_and_bound(&F1, &F2, 0) >= 0)
			return true;
	}
	return false;
}

/*Joel: Limiting starting point*/
void find_best_distance(Node * n, Node * super_tree, vector<Node *> &gene_trees, vector< pair<Node*, int> > &scores, int &num_zeros, int &best_distance){
	if(n->lchild() != NULL){
		find_best_distance(n->lchild(), super_tree, gene_trees, scores, num_zeros, best_distance);
	}
	if(n->rchild() != NULL){
		find_best_distance(n->rchild(), super_tree, gene_trees, scores, num_zeros, best_distance);
	}
	if(n->parent() != NULL){
		UndoMachine um = UndoMachine();
		vector<Node *> components = vector<Node *>();
		Node *n_parent = n->parent();
		um.add_event(new CutParent(n));
		n->cut_parent();
		ContractEvent(&um, n_parent);
		Node *post_contract = n_parent->contract();
		components.push_back(super_tree);
		components.push_back(n);
		Forest f =Forest(components);
		//f.print_components();
		components.clear();
		//find_best_distance_helper(super_tree, &f, gene_trees, scores);
		int distance;
		if(GREEDY_REFINED){
			distance = rSPR_total_approx_distance(&f, gene_trees);
		}
		else
			distance = rSPR_total_distance(&f,gene_trees);
		if (best_distance - distance == 0)
			num_zeros++;
		f.erase_components();
		um.undo_all();
		scores.push_back(make_pair(n,distance));

	}
		
}

vector<pair<Node *, int> > get_best_scores(vector<pair<Node *, int> > &scores){
	int count = 0;
	std::sort(scores.begin(), scores.end());
	vector< pair<Node*, int> > best_scores = vector< pair<Node*, int> >();
	for(vector<pair<Node *, int> >::const_iterator i = scores.begin(); count != S_NUM && i!=scores.end(); i++){
		best_scores.push_back(make_pair(i->first,i->second));
		count++;
	}
	return best_scores;	
}
bool operator<(const pair<Node *, int> &a, const pair<Node *, int> &b){
	return a.second < b.second;
}
/*end*/

/*Joel: Starting point stats*/
void find_best_spr(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, vector <pair <int, pair<int, int> > > &stats) {
	int min_distance = INT_MAX;
	int num_ties = 0;
	find_best_spr_helper(super_tree, super_tree, gene_trees,
			best_spr_move, best_sibling, min_distance, num_ties, stats);

}

void find_best_spr_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <int, pair<int, int> > > &stats) {
	if (n->lchild() != NULL) {
		find_best_spr_helper(n->lchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, stats);
	}
	if (n->rchild() != NULL) {
		find_best_spr_helper(n->rchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, stats);
	}

		find_best_spr_helper(n, super_tree, super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, stats);

}

void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <int, pair<int, int> > > &stats) {
	// do not consider invalid spr moves to n's subtree
	if (new_sibling == n)
		return;

	// recurse
	if (new_sibling->lchild() != NULL) {
		find_best_spr_helper(n, new_sibling->lchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, stats);
	}
	if (new_sibling->rchild() != NULL) {
		find_best_spr_helper(n, new_sibling->rchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, stats);
	}
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "s: " << new_sibling->str_subtree() << endl;

		//cout << "foo1" << endl;
	if (n->parent() != NULL &&
			(new_sibling == n->parent())) {/* ||
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent())) */
//		cout << "rule 1" << endl;
		return;
	}
	if (n->parent() != NULL &&
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent()) {
//		cout << "rule 2" << endl;
		return;
	}
	if (new_sibling == n->get_sibling()) {
//		cout << "rule 3" << endl;
		return;
	}
//		cout << "foo3" << endl;
//	cout << "foo4" << endl;
	if (n != super_tree && new_sibling != n) {
		Node *old_sibling = n->get_sibling();
		//if (new_sibling != old_sibling)
/*
		cout << "SPR Move:" << endl;
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Previous Super Tree: "
		<< super_tree->str_subtree() << endl;
		cout << "Subtree: " << n->str_subtree() << endl;
		cout << "New Sibling: " << new_sibling->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/

		int which_sibling = 0;
		Node *undo = n->spr(new_sibling, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
/*
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Proposed Super Tree: " << super_tree->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/

		int distance;
		if (APPROX) {
			if (UNROOTED)
				distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
			else
				distance = rSPR_total_approx_distance(super_tree, gene_trees);
		}	
		else {
			if (UNROOTED)
				distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
			else
				distance = rSPR_total_distance(super_tree, gene_trees);
		}
		if(stats.size() == 0)
			stats.push_back(make_pair(0,make_pair(rSPR_total_distance(super_tree, gene_trees),rSPR_total_approx_distance(super_tree, gene_trees))));
		else
			stats.push_back(make_pair(++(stats.back().first),make_pair(rSPR_total_distance(super_tree, gene_trees),rSPR_total_approx_distance(super_tree, gene_trees))));
/*		cout << "\t" << distance << endl;
*/
		if (distance < min_distance) {
			min_distance = distance;
			best_spr_move = n;
			best_sibling = new_sibling;
			num_ties = 2;
		}
		else if (distance == min_distance) {
			int r = rand();
			if (r < RAND_MAX/num_ties) {
				min_distance = distance;
				best_spr_move = n;
				best_sibling = new_sibling;
			}
			num_ties++;
		}
		// restore the previous tree
		n->spr(undo, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
//		cout << "Reverted Super Tree: "
//	<< super_tree->str_subtree() << endl;
	}

}
/*end*/

/*Joel: Greedy Search*/
void find_best_spr(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, vector<pair <pair<Node*,Node*>, int> > &approx_moves) {
	int min_distance = INT_MAX;
	int num_ties = 0;
	find_best_spr_helper(super_tree, super_tree, gene_trees,
			best_spr_move, best_sibling, min_distance, num_ties, approx_moves);

}

void find_best_spr_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <pair<Node*,Node*>, int> > &approx_moves) {
	if (n->lchild() != NULL) {
		find_best_spr_helper(n->lchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, approx_moves);
	}
	if (n->rchild() != NULL) {
		find_best_spr_helper(n->rchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, approx_moves);
	}

//	cout << "n: " << n->str_subtree() << endl;
	vector<Node *> leaves = n->find_leaves();
//	vector<Node *>::iterator leaf;
//	for(leaf = leaves.begin(); leaf != leaves.end(); leaf++) {
//		cout << " " << (*leaf)->get_name() << endl;
//	}
//	cout << endl;
	vector<Node *> current_gene_trees = vector<Node *>();
	int offset = 0;
	for(int i = 0; i < gene_trees.size(); i++) {
		bool include = false;
		for(int j = 0; j < leaves.size(); j++) {
			if (gene_trees[i]->contains_leaf(atoi(leaves[j]->get_name().c_str()))) {
				include = true;
				break;
			}
		}
		//int r = rand();
		if (include) {// && r < RAND_MAX/10) {
				current_gene_trees.push_back(gene_trees[i]);
		}
		else {
			Forest F1 = Forest(super_tree);
			Forest F2 = Forest(gene_trees[i]);
			offset += rSPR_worse_3_approx(&F1, &F2);
		}
	}
//	cout << gene_trees.size() << endl;
//	cout << current_gene_trees.size() << endl;
//	cout << offset << endl;

		find_best_spr_helper(n, super_tree, super_tree,
				current_gene_trees, best_spr_move, best_sibling, min_distance, num_ties, approx_moves, offset);

}

void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <pair<Node*,Node*>, int> > &approx_moves) {
	find_best_spr_helper(n, new_sibling, super_tree, gene_trees,
			best_spr_move, best_sibling, min_distance, num_ties, approx_moves, 0);
}
void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, vector<pair <pair<Node*,Node*>, int> > &approx_moves, int offset) {
	// do not consider invalid spr moves to n's subtree
	if (new_sibling == n)
		return;
	// recurse
	if (new_sibling->lchild() != NULL) {
		find_best_spr_helper(n, new_sibling->lchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, approx_moves, offset);
	}
	if (new_sibling->rchild() != NULL) {
		find_best_spr_helper(n, new_sibling->rchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance, num_ties, approx_moves, offset);
	}
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "s: " << new_sibling->str_subtree() << endl;

		//cout << "foo1" << endl;
	if (n->parent() != NULL &&
			(new_sibling == n->parent())) {/* ||
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent())) */
//		cout << "rule 1" << endl;
		return;
	}
	if (n->parent() != NULL &&
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent()) {
//		cout << "rule 2" << endl;
		return;
	}
	if (new_sibling == n->get_sibling()) {
//		cout << "rule 3" << endl;
		return;
	}
//		cout << "foo3" << endl;
//	cout << "foo4" << endl;
	if (n != super_tree && new_sibling != n) {
		Node *old_sibling = n->get_sibling();
		//if (new_sibling != old_sibling)
/*
		cout << "SPR Move:" << endl;
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Previous Super Tree: "
		<< super_tree->str_subtree() << endl;
		cout << "Subtree: " << n->str_subtree() << endl;
		cout << "New Sibling: " << new_sibling->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/

		int which_sibling = 0;
		Node *undo = n->spr(new_sibling, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
/*
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Proposed Super Tree: " << super_tree->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/

		int distance;
		distance = rSPR_total_approx_distance(super_tree, gene_trees);
		distance += offset;
		approx_moves.push_back(make_pair(make_pair(n,new_sibling), distance));



		/*if (distance < min_distance) {
			min_distance = distance;
			best_spr_move = n;
			best_sibling = new_sibling;
			num_ties = 2;
		}
		else if (distance == min_distance) {
			int r = rand();
			if (r < RAND_MAX/num_ties) {
				min_distance = distance;
				best_spr_move = n;
				best_sibling = new_sibling;
			}
			num_ties++;
		}*/
		// restore the previous tree
		n->spr(undo, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
//		cout << "Reverted Super Tree: "
//	<< super_tree->str_subtree() << endl;
	}

}


bool sort_approx_moves(const pair<pair<Node*, Node*>, int> &a, const pair<pair<Node*, Node*>, int> &b){
	return a.second<b.second;
}
/*end*/

/*Joel: Limit SPR Radius*/
void find_best_spr_r(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, int r) {
	find_best_spr_r(super_tree, gene_trees, best_spr_move, best_sibling,
			r, NULL);
}
void find_best_spr_r(Node *super_tree, vector<Node *> &gene_trees, Node *&best_spr_move, Node *&best_sibling, int r, vector<int> *original_scores) {
	int min_distance = INT_MAX;
	int min_tie_distance = INT_MAX;
	int num_ties = 0;
	find_best_spr_r_helper(super_tree, super_tree, gene_trees,
			best_spr_move, best_sibling, min_distance, min_tie_distance,
			num_ties, r, original_scores);
}

void find_best_spr_r_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, int r) {
	int min_tie_distance = INT_MAX;
	find_best_spr_r_helper(n, super_tree, gene_trees, best_spr_move,
			best_sibling, min_distance, min_tie_distance, num_ties, r, NULL);
}

void find_best_spr_r_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &min_tie_distance,
		int &num_ties, int r, vector<int> *original_scores) {

	if (n->lchild() != NULL) {
		find_best_spr_r_helper(n->lchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance,
				min_tie_distance, num_ties, r, original_scores);
	}
	if (n->rchild() != NULL) {
		find_best_spr_r_helper(n->rchild(), super_tree,
				gene_trees, best_spr_move, best_sibling, min_distance,
				min_tie_distance, num_ties, r, original_scores);
	}
	vector<Node *> current_gene_trees;
	vector<Node *> *gene_trees_p = &gene_trees;
	int offset = 0;
	if (C_SOURCE != 1)
		cout << "\r \r";
	cout << C_SOURCE << "/" << NUM_SOURCE << flush;
	C_SOURCE++;
	
	// TODO: check here if there will be any moves
	// if not, then we do not need to look at the gene trees
	if (original_scores != NULL) {
//		cout << "selecting gene_trees" << endl;
		vector<Node *> leaves = n->find_leaves();
		current_gene_trees = vector<Node *>();
		for(int i = 0; i < gene_trees.size(); i++) {
			bool include = false;
			for(int j = 0; j < leaves.size(); j++) {
				if (gene_trees[i]->contains_leaf(atoi(leaves[j]->get_name().c_str()))) {
					include = true;
					break;
				}
			}
			if (include)
					current_gene_trees.push_back(gene_trees[i]);
			else {
				offset += (*original_scores)[i];
			}
		}
		gene_trees_p = &current_gene_trees;
//		cout << "gene_trees: " << gene_trees.size() << endl;
//		cout << "c_gene_trees: " << current_gene_trees.size() << endl;
	}
// TODO FROM HERE
// TODO: I think this works for the sibling subtree but goes one
//	to far in the grandparent

////	cout << endl;
	if(n->parent() != NULL) {
		if(!R_CONTROL  && R_RAND)
			r = find_r(R_PROB);
		if(n->parent()->lchild() == n)
			find_best_spr_r_helper(n, n->parent(), super_tree, *gene_trees_p,
					best_spr_move, best_sibling, min_distance, min_tie_distance,
					num_ties, r+1, 1, offset);
		else if(n->parent()->rchild() == n)
			find_best_spr_r_helper(n, n->parent(), super_tree, *gene_trees_p,
					best_spr_move, best_sibling, min_distance, min_tie_distance,
					num_ties, r+1, 2, offset);
	}

}
// TODO FROM HERE
void find_best_spr_r_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, int r,
		int origin) {
	int min_tie_distance = INT_MAX;
	find_best_spr_r_helper(n, new_sibling, super_tree, gene_trees,
			best_spr_move, best_sibling, min_distance, min_tie_distance,
			num_ties, r, origin, 0);
}

/*
 origin 1: left child to parent
 origin 2: right child to parent
 origin 3: parent to either child
*/
void find_best_spr_r_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &min_tie_distance,
		int &num_ties, int r, int origin, int offset) {
	// do not consider invalid spr moves to n's subtree
	if (new_sibling == n)
		return;


	// recurse
//		n->numbers_to_labels(&reverse_label_map);
//		new_sibling->numbers_to_labels(&reverse_label_map);
//	cout << n->str_subtree() << "  " << new_sibling->str_subtree() << " origin " << origin << endl;
//		n->labels_to_numbers(&label_map, &reverse_label_map);	
//		new_sibling->labels_to_numbers(&label_map, &reverse_label_map);	
	if(origin == 1 && r > 0 ){
		//cout << "n1: " << n->str_subtree() << "\t\tnew_sibling: " << new_sibling->str_subtree() << endl;
		if(new_sibling->parent() != NULL
				&& (!BIPARTITION_CLUSTER || new_sibling->get_support() < SUPPORT_THRESHOLD
					|| new_sibling->parent()->parent() == NULL)) {//!new_sibling->is_protected())){
			if(new_sibling->parent()->lchild() == new_sibling){
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
						num_ties, --r, 1, offset);	
				//cout << "Origin1, Target: " << new_sibling->parent()->str_subtree() << endl;
			}
			else{
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
						num_ties, --r, 2, offset);	
				//cout << "Origin1, Target: " << new_sibling->parent()->str_subtree() << endl;
			}
		}
		if(new_sibling->rchild() != NULL) {
//				&& (!BIPARTITION_CLUSTER
//						|| new_sibling->rchild()->is_protected())){
			//cout << "Going to Origin 3\n";
			//cout << "n to 3: " << n->str_subtree() << "\tsibling to 3: " << new_sibling->rchild()->str_subtree() << endl;
			find_best_spr_r_helper(n, new_sibling->rchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
					num_ties, --r, 3, offset);
			//cout << "Origin1, Target: " << new_sibling->rchild()->str_subtree() << endl;
		}

	}

	if(origin == 2 && r > 0){
		//cout << "n2: " << n->str_subtree() << "\t\tnew_sibling: " << new_sibling->str_subtree() << endl;
		if(new_sibling->parent() != NULL
				&& (!BIPARTITION_CLUSTER || new_sibling->get_support() < SUPPORT_THRESHOLD
					|| new_sibling->parent()->parent() == NULL)) {//!new_sibling->is_protected())){
			if(new_sibling->parent()->lchild() == new_sibling){
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
						num_ties, --r, 1, offset);
				//cout << "Origin2, Target: " << new_sibling->parent()->str_subtree() << endl;	
			}			
			else{
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
						num_ties, --r, 2, offset);	
				//cout << "Origin2, Target: " << new_sibling->parent()->str_subtree() << endl;	
			}
		}		
			if(new_sibling->lchild() != NULL) {
//				&& (!BIPARTITION_CLUSTER
//						|| new_sibling->lchild()->is_protected())){
			//cout << "Going to Origin 3\n";
			//cout << "n to 3: " << n->str_subtree() << "\tsibling to 3: " << new_sibling->lchild()->str_subtree() << endl;
			//cout << "Origin2, Target: " << new_sibling->lchild()->str_subtree() << endl;
			find_best_spr_r_helper(n, new_sibling->lchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
					num_ties, --r, 3, offset);
		}

	}
		
	
	if(origin == 3 && r > 0){
		//cout << "n3: " << n->str_subtree() << "\t\tnew_sibling: " << new_sibling->str_subtree() << endl;
		if(new_sibling->lchild() != NULL
				&& (!BIPARTITION_CLUSTER
//						|| !new_sibling->is_protected())){
				|| new_sibling->get_support() < SUPPORT_THRESHOLD)) {//!new_sibling->is_protected())){
			find_best_spr_r_helper(n, new_sibling->lchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
					num_ties, --r, 3, offset);
			//cout << "Origin3, Target: " << new_sibling->lchild()->str_subtree() << endl;
		}
		if(new_sibling->rchild() != NULL
				&& (!BIPARTITION_CLUSTER
//						|| !new_sibling->is_protected())){
				|| new_sibling->get_support() < SUPPORT_THRESHOLD)) {//!new_sibling->is_protected())){
			find_best_spr_r_helper(n, new_sibling->rchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, min_tie_distance,
					num_ties, --r, 3, offset);
			//cout << "Origin3, Target: " << new_sibling->rchild()->str_subtree() << endl;
		}
	}
//	cout << n->str_subtree() << "  " << new_sibling->str_subtree() << " origin " << origin << endl;
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "s: " << new_sibling->str_subtree() << endl;

		//cout << "foo1" << endl;
	if (n->parent() != NULL &&
			(new_sibling == n->parent())) {/* ||
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent())) */
//		cout << "rule 1" << endl;
		return;
	}
	if (n->parent() != NULL && new_sibling->parent() != NULL && n->parent()->parent() == new_sibling->parent()) {
//		cout << "rule 2" << endl;
		return;
	}
	
	if (new_sibling == n->get_sibling()) {
//		cout << "rule 3" << endl;
		return;
	}
//		cout << "foo3" << endl;
//		cout << "foo4" << endl;
	if (n != super_tree && new_sibling != n) {
		Node *old_sibling = n->get_sibling();
		//if (new_sibling != old_sibling)

/*		cout << endl;
		cout << "SPR Move:" << endl;
//		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Previous Super Tree: "
		<< super_tree->str_subtree() << endl;
		cout << "Subtree: " << n->str_subtree() << endl;
		cout << "New Sibling: " << new_sibling->str_subtree() << endl;
//		super_tree->labels_to_numbers(&label_map, &reverse_label_map);	
*/

//		cout << endl;
//		cout << "Previous Super Tree: "
//		<< super_tree->str_support_subtree(true) << endl;
		int which_sibling = 0;
		Node *undo = n->spr(new_sibling, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
//		cout << "Proposed Super Tree: "
//		<< super_tree->str_support_subtree(true) << endl;
/*
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Proposed Super Tree: " << super_tree->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/	
		int distance;
		if (APPROX) {
			if (UNROOTED)
				distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
			else
				distance = rSPR_total_approx_distance(super_tree, gene_trees);
		}
		else {
			if (UNROOTED)
				distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
			else
				distance = rSPR_total_distance(super_tree, gene_trees);
		}
//		cout << "\t" << distance << endl;

//		cout << "After SPR Distance Super Tree: "
//		<< super_tree->str_support_subtree(true) << endl;
		distance += offset;
		if (distance < min_distance) {
			if (!TABOO_SEARCH || !is_taboo(taboo_trees, super_tree)) {
				min_distance = distance;
				if (RF_TIES)
					min_tie_distance = rf_total_distance(super_tree, gene_trees);
				best_spr_move = n;
				best_sibling = new_sibling;
				num_ties = 2;
			}
		}
		else if (distance == min_distance) {
			bool check_tie = true;
			bool tie = false;
			if (RF_TIES) {
				int rf_distance = rf_total_distance(super_tree, gene_trees);
				if (rf_distance < min_tie_distance) {
					if (!TABOO_SEARCH || !is_taboo(taboo_trees, super_tree)) {
						min_tie_distance = rf_distance;
						best_spr_move = n;
						best_sibling = new_sibling;
						num_ties = 2;
					}
					check_tie = false;
				}
				else if (rf_distance == min_tie_distance)
					tie = true;
				else
					check_tie = false;
			}
			if (check_tie) {
				int r = rand();
				if (r < RAND_MAX/num_ties) {
					if (!TABOO_SEARCH || !is_taboo(taboo_trees, super_tree)) {
						min_distance = distance;
						best_spr_move = n;
						best_sibling = new_sibling;
					}
				}
				num_ties++;
			}
		}
		// restore the previous tree
		n->spr(undo, which_sibling);

//		cout << "Reverted Super Tree: "
//		<< super_tree->str_support_subtree(true) << endl;

		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
//		cout << "Reverted Super Tree: "
//		<< super_tree->str_subtree() << endl;

	}
}
/*End*/

/*Joel: Refined Greedy Search*/
void find_best_spr_r_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties, int r, int origin, vector<pair <pair<Node*,Node*>, int> > &approx_moves) {
	// do not consider invalid spr moves to n's subtree
	if (new_sibling == n)
		return;


	// recurse
	if(origin == 1 && r > 0 ){
		//cout << "n1: " << n->str_subtree() << "\t\tnew_sibling: " << new_sibling->str_subtree() << endl;
		if(new_sibling->parent() != NULL){
			if(new_sibling->parent()->lchild() == new_sibling){
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 1, approx_moves);	
				//cout << "Origin1, Target: " << new_sibling->parent()->str_subtree() << endl;
			}
			else{
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 2, approx_moves);	
				//cout << "Origin1, Target: " << new_sibling->parent()->str_subtree() << endl;
			}
		}
		if(new_sibling->rchild() != NULL){
			//cout << "Going to Origin 3\n";
			//cout << "n to 3: " << n->str_subtree() << "\tsibling to 3: " << new_sibling->rchild()->str_subtree() << endl;
			find_best_spr_r_helper(n, new_sibling->rchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 3, approx_moves);
			//cout << "Origin1, Target: " << new_sibling->rchild()->str_subtree() << endl;
		}

	}

	if(origin == 2 && r > 0){
		//cout << "n2: " << n->str_subtree() << "\t\tnew_sibling: " << new_sibling->str_subtree() << endl;
		if(new_sibling->parent() != NULL){
			if(new_sibling->parent()->lchild() == new_sibling){
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 1, approx_moves);
				//cout << "Origin2, Target: " << new_sibling->parent()->str_subtree() << endl;	
			}			
			else{
				find_best_spr_r_helper(n, new_sibling->parent(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 2, approx_moves);	
				//cout << "Origin2, Target: " << new_sibling->parent()->str_subtree() << endl;	
			}
		}		
			if(new_sibling->lchild() != NULL){
			//cout << "Going to Origin 3\n";
			//cout << "n to 3: " << n->str_subtree() << "\tsibling to 3: " << new_sibling->lchild()->str_subtree() << endl;
			//cout << "Origin2, Target: " << new_sibling->lchild()->str_subtree() << endl;
			find_best_spr_r_helper(n, new_sibling->lchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 3, approx_moves);
		}

	}
		
	
	if(origin == 3 && r > 0){
		//cout << "n3: " << n->str_subtree() << "\t\tnew_sibling: " << new_sibling->str_subtree() << endl;
		if(new_sibling->lchild() != NULL){
			find_best_spr_r_helper(n, new_sibling->lchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 3, approx_moves);
			//cout << "Origin3, Target: " << new_sibling->lchild()->str_subtree() << endl;
		}
		if(new_sibling->rchild() != NULL){
			find_best_spr_r_helper(n, new_sibling->rchild(), super_tree, gene_trees, best_spr_move, best_sibling, min_distance, num_ties, --r, 3, approx_moves);
			//cout << "Origin3, Target: " << new_sibling->rchild()->str_subtree() << endl;
		}
	}
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "s: " << new_sibling->str_subtree() << endl;

		//cout << "foo1" << endl;
	if (n->parent() != NULL &&
			(new_sibling == n->parent())) {/* ||
			new_sibling->parent() != NULL &&
			n->parent()->parent() == new_sibling->parent())) */
//		cout << "rule 1" << endl;
		return;
	}
	if (n->parent() != NULL && new_sibling->parent() != NULL && n->parent()->parent() == new_sibling->parent()) {
//		cout << "rule 2" << endl;
		return;
	}
	if (new_sibling == n->get_sibling()) {
//		cout << "rule 3" << endl;
		return;
	}
//		cout << "foo3" << endl;
//		cout << "foo4" << endl;
	if (n != super_tree && new_sibling != n) {
		Node *old_sibling = n->get_sibling();

		int which_sibling = 0;
		Node *undo = n->spr(new_sibling, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
/*
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Proposed Super Tree: " << super_tree->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
*/	
		int distance;
		distance = rSPR_total_approx_distance(super_tree, gene_trees);
		approx_moves.push_back(make_pair(make_pair(n,new_sibling), distance));

		// restore the previous tree
		n->spr(undo, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
		super_tree->preorder_number();
//		cout << "Reverted Super Tree: "
//		<< super_tree->str_subtree() << endl;

	}
}
/*end*/

/*Joel: Limit SPR Radius*/
int find_r(){
	if (R_BIAS)
		return find_r(R_PROB);
	else if(R_LIMIT && R_RAND)	{
		int count = 1;
		int r = rand() % 2;
		while(r != 1 ){
			count++;
			r = rand() % 2;
		}
		return count;
	}
	else if(R_LIMIT)
		return R_DISTANCE;
	else 
		return 0;
}

int find_r(double probability){
	int count = 1;
	double r = ((double)rand()/(double)RAND_MAX);
	while (r <= probability){
		count++;
		r = ((double)rand()/(double)RAND_MAX);
	}
	return count;
}
/*end*/

bool supported_spr(Node *source, Node *target) {
	// a supported source can still be moved
	source = source->parent();
	target = target->parent();
	while(source != target) {
		if (source != NULL && (target == NULL ||   source->get_depth() > target->get_depth())) {
			if (source->parent() != NULL && source->parent()->parent() != NULL
					&& source->get_support() >= SUPPORT_THRESHOLD) {
				return false;
			}
			else {
				source = source->parent();
			}
		}
		else {
			if (target->parent() != NULL && target->parent()->parent() != NULL && target->get_support() >= SUPPORT_THRESHOLD) {
				return false;
			}
			else {
				target = target->parent();
			}
		}
	}
	return true;
}

bool pair_comparator (pair<int, int> a, pair<int, int> b) {
	return a.first > b.first;
}
