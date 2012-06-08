/*******************************************************************************
rspr.cpp

Usage: rspr [OPTIONS]
Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of rooted binary trees from STDIN in newick format.
Supports arbitrary labels. See the README for more information.

Copyright 2009-2012 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
May 3, 2012
Version 1.03

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
just that subset of optimizations. See the README for more information.

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
#include <list>
#include "rspr.h"

#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "UndoMachine.h"

using namespace std;


// options to pick default
bool DEFAULT_ALGORITHM=true;
bool DEFAULT_OPTIMIZATIONS=true;


bool FPT = false;
bool QUIET = false;
bool UNROOTED = false;
bool LCA_TEST = false;
bool CLUSTER_TEST = false;
bool TOTAL = false;
bool APPROX = false;
bool LOWER_BOUND = false;
int MULTI_TEST = 0;

string USAGE =
"rspr, version 1.02\n"
"\n"
"usage: rspr [OPTIONS]\n"
"Calculate approximate and exact Subtree Prune and Regraft (rSPR)\n"
"distances and the associated maximum agreement forests (MAFs) between pairs\n"
"of rooted binary trees from STDIN in newick format.\n"
"Supports arbitrary labels. See the README for more information.\n"
"\n"
"Copyright 2009-2011 Chris Whidden\n"
"whidden@cs.dal.ca\n"
"http://kiwi.cs.dal.ca/Software/RSPR\n"
"November 2, 2011\n"
"Version 1.02\n"
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
"            FPT algorithm. This is the default option.\n"
"\n"
"-approx     Calculate just a linear -time 3-approximation of the\n"
"            rSPR distance\n"
"\n"
"\n"
"-cluster_test   Use the cluster reduction to speed up the exact algorithm.\n"
"                This is enabled by default.\n"
"\n"
"-total          Find the total SPR distance from the first input tree to\n"
"                the rest of the list of trees. Uses the other algorithm\n"
"                options as specified (including unrooted options).\n"
"*******************************************************************************\n"
"OPTIMIZATIONS\n"
"*******************************************************************************\n"
"\n"
"These options control the use of optimized branching. All optimizations are\n"
"enabled by default. Specifying any subset of -cob, -cab, and -sc will use\n"
"just that subset of optimizations. See the README for more information.\n"
"\n"
"-allopt     Use -cob -cab -sc. This is the default option\n"
"\n"
"-noopt      Use 3-way branching for all FPT algorithms\n"
"\n"
"-cob        Use \"cut one b\" improved branching\n"
"\n"
"-cab        Use \"cut all b\" improved branching\n"
"\n"
"-sc         Use \"separate components\" improved branching\n"
"\n"
"*******************************************************************************\n"
"UNROOTED COMPARISON OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-unrooted   Compare the first input tree to each other input tree.\n"
"            Output the best found distance and agreement forest\n"
"-unrooted_min_approx\n"
"            Compare the first input tree to each other input tree.\n"
"            Run the exact algorithms on the pair with the\n"
"            minimum approximate rspr distance\n"
"\n"
"*******************************************************************************\n"
"OTHER OPTIONS\n"
"*******************************************************************************\n"
"-cc         Calculate a potentially better approximation with a quadratic\n"
"            time algorithm\n"
"\n"
"-q          Quiet; Do not output the input trees or approximation\n"
"*******************************************************************************\n";

int main(int argc, char *argv[]) {
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
			APPROX_REVERSE_CUT_ONE_B = true;
			APPROX_EDGE_PROTECTION = true;
		}
		else if (strcmp(arg, "-a_cob") == 0) {
			APPROX_CUT_ONE_B = true;
		}
		else if (strcmp(arg, "-a_c2b") == 0) {
			APPROX_CUT_TWO_B = true;
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
		else if (strcmp(arg, "-cut_two_b") == 0 ||
				strcmp(arg, "-c2b") == 0) {
			CUT_TWO_B = true;
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
					SPLIT_APPROX_THRESHOLD = atof(arg2);
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
		else if (strcmp(arg, "--help") == 0) {
			cout << USAGE;
			return 0;
		}
			
	}
	if (DEFAULT_OPTIMIZATIONS) {
		CUT_ALL_B=true;
		CUT_ONE_B = true;
		REVERSE_CUT_ONE_B = true;
		CUT_TWO_B = true;
		CUT_AC_SEPARATE_COMPONENTS = true;
		EDGE_PROTECTION = true;
//		if (ALL_MAFS == false)
//			ABORT_AT_FIRST_SOLUTION = true;
//		PREORDER_SIBLING_PAIRS = true;
		NEAR_PREORDER_SIBLING_PAIRS = true;
		LEAF_REDUCTION = true;
		//LEAF_REDUCTION2 = true;

		APPROX_CUT_ONE_B = true;
		APPROX_CUT_TWO_B = true;
		APPROX_REVERSE_CUT_ONE_B = true;
		APPROX_EDGE_PROTECTION = true;
	}
	PREORDER_SIBLING_PAIRS = true;
	if (DEFAULT_ALGORITHM) {
		BB=true;
		CLUSTER_TEST = true;
		PREFER_RHO = true;
	}

	// Label maps to allow string labels
	map<string, int> label_map= map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();

	// set random seed
	srand(unsigned(time(0)));


	// Normal operation
	if (!UNROOTED && !UNROOTED_MIN_APPROX && !TOTAL) {
		string T1_line = "";
		string T2_line = "";
		while (getline(cin, T1_line) && getline(cin, T2_line)) {
			Node *T1 = build_tree(T1_line);
			Node *T2 = build_tree(T2_line);

			if (MULTI_TEST > 0) {
				vector<Node *> interior = T2->find_interior();
				vector<Node *> remove = random_select(interior, MULTI_TEST);
				vector<Node *>::iterator i;
				for(i = remove.begin(); i != remove.end(); i++) {
					(*i)->contract_node();
				}
			}
			// TODO: should we sync here to prune out additional leaves?
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
				return 0;
			}

			T1->labels_to_numbers(&label_map, &reverse_label_map);
			T2->labels_to_numbers(&label_map, &reverse_label_map);
	
//			ClusterForest F1 = ClusterForest(T1);
//			ClusterForest F2 = ClusterForest(T2);
//			ClusterForest F3 = ClusterForest(T1);
//			ClusterForest F4 = ClusterForest(T2);
			ClusterForest F1 = ClusterForest(T1);
			ClusterForest F2 = ClusterForest(T2);
			Forest F3 = Forest(T1);
			Forest F4 = Forest(T2);

			if (CLUSTER_TEST) {
				int exact_k = rSPR_branch_and_bound_simple_clustering(T1,T2,true, &label_map, &reverse_label_map);
				//int exact_k = rSPR_branch_and_bound_simple_clustering(&F3,&F4,true, &label_map, &reverse_label_map);

				T1->delete_tree();
				T2->delete_tree();
				continue;
			}
	
			// APPROX ALGORITHM
			int approx_spr = rSPR_worse_3_approx(&F1, &F2);
			int min_spr = approx_spr / 3;
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
				int exact_spr = rSPR_branch_and_bound(&F1, &F2);
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
				if (exact_spr == -1)
						if (FPT)
						cout << "exact drSPR=?  " << "k=" << k << " too large"
							<< endl;
						else
						cout << "exact BB drSPR=?  " << "k=" << k << " too large"
							<< endl;
				cout << "\n";
			}
	
			// cleanup
			T1->delete_tree();
			T2->delete_tree();
		}
	}
	// Comparison between a rooted tree and all rootings of an unrooted tree
	else if (!TOTAL && (UNROOTED || UNROOTED_MIN_APPROX)) {
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
		int min_spr = (int)1E9;
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
		if (!getline(cin, line))
			return 0;
		Node *T1 = build_tree(line);
		if (!QUIET) {
			cout << "T1: ";
			T1->print_subtree();
		}
		T1->labels_to_numbers(&label_map, &reverse_label_map);
		while (getline(cin, line)) {
			if (UNROOTED)
				line = root(line);
			Node *T2 = build_tree(line);
			if (!QUIET) {
				cout << "T2: ";
				T2->print_subtree();
			}
			T2->labels_to_numbers(&label_map, &reverse_label_map);
			if (UNROOTED)
				T2->preorder_number();
			trees.push_back(T2);
		}
		cout << endl;

		int distance;
		if (APPROX) {
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
		T1->delete_tree();
		for(vector<Node *>::iterator T2 = trees.begin(); T2 != trees.end(); T2++)
			(*T2)->delete_tree();
	}
	return 0;
}
