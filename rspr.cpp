/*******************************************************************************
rspr.cpp

Usage: rspr [OPTIONS]
Calculate approximate and exact Subtree Prune and Regraft (rSPR)
distances and the associated maximum agreement forests (MAFs) between pairs
of rooted binary trees from STDIN in newick format. By default, computes a
3-approximation of the rSPR distance. Supports arbitrary labels. See the
README for more information.

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

//#define DEBUG 1
#define MAX_SPR 1000

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
#include <deque>
#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"

using namespace std;

int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
		deque<Node *> *sibling_pairs);
int rSPR_3_approx(Forest *T1, Forest *T2);
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
		deque<Node *> *sibling_pairs);
int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync);
int rSPR_worse_3_approx(Forest *T1, Forest *T2);
int rSPR_FPT_hlpr(Forest *T1, Forest *T2, int k, deque<Node *> *sibling_pairs,
		list<Node *> *singletons, bool cut_b_only);
int rSPR_FPT(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound(Forest *T1, Forest *T2);
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k);
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k);
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
		int end_k);
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		deque<Node *> *sibling_pairs, list<Node *> *singletons, bool cut_b_only);
inline void copy_trees(Forest **T1, Forest **T2, deque<Node *> **sibling_pairs,
		Node **T1_a, Node **T1_c, Node **T2_a, Node **T2_c,
		Forest **T1_copy, Forest **T2_copy, deque<Node *> **sibling_pairs_copy,
		Node **T1_a_copy, Node **T1_c_copy, Node **T2_a_copy, Node **T2_c_copy);

// options to pick default
bool DEFAULT_ALGORITHM=true;
bool DEFAULT_OPTIMIZATIONS=true;


bool FPT = false;
bool BB = false;
bool QUIET = false;
bool APPROX_CHECK_COMPONENT = false;
bool UNROOTED = false;
bool UNROOTED_MIN_APPROX = false;
bool CUT_ONE_B = false;
bool CUT_ALL_B = false;
bool CUT_AC_SEPARATE_COMPONENTS = false;
bool CUT_ONE_AB = false;
bool LCA_TEST = false;
bool CLUSTER_TEST = false;
bool CLUSTER_REDUCTION = false;
bool PREFER_RHO = false;

string USAGE =
"rspr, version 1.01\n"
"\n"
"usage: rspr [OPTIONS]\n"
"Calculate approximate and exact Subtree Prune and Regraft (rSPR)\n"
"distances and the associated maximum agreement forests (MAFs) between pairs\n"
"of rooted binary trees from STDIN in newick format. By default, computes a\n"
"3-approximation of the rSPR distance. Supports arbitrary labels. See the\n"
"README for more information.\n"
"\n"
"Copyright 2009-2010 Chris Whidden\n"
"whidden@cs.dal.ca\n"
"http://kiwi.cs.dal.ca/Software/RSPR\n"
"March 22, 2010\n"
"Version 1.01\n"
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
"-cc         Calculate a potentially better approximation with a quadratic time\n"
"            algorithm\n"
"\n"
"-q          Quiet; Do not output the input trees or approximation\n"
"*******************************************************************************\n";

int main(int argc, char *argv[]) {
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
		else if (strcmp(arg, "-cluster") == 0) {
			CLUSTER_REDUCTION = true;
			PREFER_RHO = true;
		}
		else if (strcmp(arg, "-cluster_test") == 0) {
			CLUSTER_TEST = true;
			PREFER_RHO = true;
		}
		else if (strcmp(arg, "-prefer_rho") == 0) {
			PREFER_RHO = true;
		}
		else if (strcmp(arg, "--help") == 0) {
			cout << USAGE;
			return 0;
		}
			
	}
	if (DEFAULT_OPTIMIZATIONS) {
		CUT_ALL_B=true;
		CUT_ONE_B = true;
		CUT_AC_SEPARATE_COMPONENTS = true;
	}
	if (DEFAULT_ALGORITHM) {
		BB=true;
	}

	// Label maps to allow string labels
	map<string, int> label_map= map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();

	// Normal operation
	if (!UNROOTED && !UNROOTED_MIN_APPROX) {
		string T1_line = "";
		string T2_line = "";
		while (getline(cin, T1_line) && getline(cin, T2_line)) {
			Node *T1 = build_tree(T1_line);
			Node *T2 = build_tree(T2_line);
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
				exit(0);
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
				int approx_spr;
				//int approx_spr = rSPR_3_approx(&F3, &F4);
				//	cout << "approx drSPR=" << approx_spr << endl;
				//	cout << "\n";

				sync_twins(&F1, &F2);
				//sync_interior_twins_real(&F1, &F2);
				sync_interior_twins(&F1, &F2);
				list<Node *> *cluster_points = find_cluster_points(&F1);
//				if (!cluster_points->empty()) {
					/*
					cout << "CLUSTER POINTS:" << endl;
					for(list<Node *>::iterator i = cluster_points->begin();
							i != cluster_points->end(); i++) {
						string s1 = (*i)->str_subtree();
						string s2 = (*i)->get_twin()->str_subtree();
						if (s1 != s2)
							cout << s1 << "\t" << s2 << endl;
					}
					cout << endl;
					*/
					// TODO: not quite right because we need to put
					// in the leaf with a ref and not do contractions
					for(list<Node *>::iterator i = cluster_points->begin();
							i != cluster_points->end(); i++) {
						stringstream ss;
						string cluster_name = "X";
						ss << F1.size();
						cluster_name += ss.str();
						int num_labels = label_map.size();
						label_map.insert(make_pair(cluster_name,num_labels));
						reverse_label_map.insert(
								make_pair(num_labels,cluster_name));
						ss.str("");
						ss << num_labels;
						cluster_name = ss.str();

						Node *n = *i;
						Node *n_parent = n->parent();
						Node *twin = n->get_twin();
						Node *twin_parent = twin->parent();

						F1.add_cluster(n,cluster_name);

						F2.add_cluster(twin,cluster_name);

						Node *n_cluster =
								F1.get_cluster_node(F1.num_clusters()-1);
						Node *twin_cluster =
								F2.get_cluster_node(F2.num_clusters()-1);
						n_cluster->set_twin(twin_cluster);
						twin_cluster->set_twin(n_cluster);
					}
					cout << endl << "CLUSTERS" << endl;
					//cout << "F1: " << endl;
					//F1.print_components("\n");
					//cout << endl;
					//cout << "F2: " << endl;
					//F2.print_components("\n");

					// component 0 needs to be done last
					F1.add_component(F1.get_component(0));
					F2.add_component(F2.get_component(0));

					int exact_spr;
					int k;
					int num_clusters = F1.num_components();
					int total_k = 0;
					for(int i = 1; i < num_clusters; i++) {
						Forest f1 = Forest(F1.get_component(i));

						Forest f2 = Forest(F2.get_component(i));
						Forest f1a = Forest(f1);
						Forest f2a = Forest(f2);
						Forest *f1_cluster;
						Forest *f2_cluster;

//						f1.numbers_to_labels(&reverse_label_map);
//						f2.numbers_to_labels(&reverse_label_map);
						cout << "C" << i << "_1: ";
						f1.print_components();
						cout << "C" << i << "_2: ";
						f2.print_components();

						int approx_spr = rSPR_worse_3_approx(&f1a, &f2a);
						if (!(QUIET && (BB || FPT)))
							cout << "cluster approx drSPR=" << approx_spr << endl;

						cout << endl;

						if (BB) {
							int min_spr = approx_spr / 3;
							for(k = min_spr; k <= MAX_SPR; k++) {
								Forest f1t = Forest(f1);
								Forest f2t = Forest(f2);
								f1t.unsync();
								f2t.unsync();
								cout << k << " ";
								cout.flush();
//								cout << "foo1" << endl;
//								f1t.print_components();
//								f1t.print_components();
//								cout << "foo1a" << endl;
								exact_spr = rSPR_branch_and_bound(&f1t, &f2t, k);
//								cout << "foo2" << endl;
								if (exact_spr >= 0) {
//									cout << "foo3" << endl;
									if (i < num_clusters - 1) {
//										cout << "foo3a" << endl;
//										cout << "rho=" << f1t.contains_rho() << endl;
//										cout << i << endl;
//										cout << &f1t << endl;
//										cout << f1t.num_components() << endl;
//										f1t.print_components();
										// TODO: problem things getting
										// deleted?
										F1.join_cluster(i,&f1t);
//										cout << "foo3b" << endl;
//										cout << "rho=" << f2t.contains_rho() << endl;
										//F2.print_components();
//										f2t.print_components();
										F2.join_cluster(i,&f2t);
									}
//									cout << "foo4" << endl;
									//f1.numbers_to_labels(&reverse_label_map);
									//f2.numbers_to_labels(&reverse_label_map);
									cout << endl;
									cout << "F" << i << "_1: ";
									f1t.print_components();
//									cout << "foo5" << endl;
									cout << "F" << i << "_2: ";
									f2t.print_components();
//									cout << "foo6" << endl;
									cout << "cluster exact drSPR=" << exact_spr << endl;
									cout << endl;
									total_k += exact_spr;
									break;
								}
							}
							if (exact_spr == -1)
							cout << "exact drSPR=?  " << "k=" << k << " too large" << endl;
							cout << "\n";
						}
					}
				//}

				if (BB) {
					F1.erase_components(0, num_clusters-1);
					F2.erase_components(0, num_clusters-1);
					cout << "F1: ";
					F1.print_components();
					cout << "F2: ";
					F2.print_components();
					cout << "total exact drSPR=" << total_k << endl;
				}

				delete cluster_points;
				exit(0);
			}
	
			// APPROX ALGORITHM
			int approx_spr = rSPR_worse_3_approx(&F1, &F2);
			int min_spr = approx_spr / 3;
			if (!(QUIET && (BB || FPT))) {
				F1.numbers_to_labels(&reverse_label_map);
				F2.numbers_to_labels(&reverse_label_map);
				cout << "F1: ";
				F1.print_components();
				cout << "F2: ";
				F2.print_components();
				cout << "approx drSPR=" << approx_spr << endl;
				cout << "\n";
			}
	
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
						break;
					}
				}
				if (exact_spr == -1)
					cout << "exact drSPR=?  " << "k=" << k << " too large" << endl;
				cout << "\n";
			}
		
			if (BB) {
				// BRANCH AND BOUND FPT ALGORITHM
				Forest F1 = Forest(F3);
				Forest F2 = Forest(F4);
				exact_spr = rSPR_branch_and_bound(&F1, &F2);
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
					cout << "exact BB drSPR=" << exact_spr << endl;
					break;
				}
				if (exact_spr == -1)
						cout << "exact BB drSPR=?  " << "k=" << k << " too large"
							<< endl;
				cout << "\n";
			}
	
			// cleanup
			delete T1;
			delete T2;
		}
	}
	// Comparison between a rooted tree and all rootings of an unrooted tree
	else if (UNROOTED || UNROOTED_MIN_APPROX) {
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

		// FPT ALGORITHM
		int exact_spr = -1;
		int k = min_spr;
		if (FPT) {
			for(k = min_spr; k <=MAX_SPR;  k++) {
				cout << k << " ";
				cout.flush();
				for (int i = 0; i < trees.size(); i++) {
					Forest F3 = Forest(F1);
					Forest F4 = Forest(trees[i]);
					exact_spr = rSPR_FPT(&F3, &F4, k);
					if (exact_spr >= 0) {
						sync_twins(&F1, &trees[i]);
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
						cout << "exact drSPR=" << exact_spr << endl;
						break;
					}
				}
				if (exact_spr >= 0) {
					break;
				}
			}
			if (exact_spr == -1)
				cout << "exact drSPR=?  " << "k=" << k << " too large" << endl;
			cout << "\n";
		}
	
		if (BB) {
			// BRANCH AND BOUND FPT ALGORITHM
			for(k = min_spr; k <=MAX_SPR;  k++) {
				cout << k << " ";
				cout.flush();
				for (int i = 0; i < trees.size(); i++) {
					Forest F3 = Forest(F1);
					Forest F4 = Forest(trees[i]);
					exact_spr = rSPR_branch_and_bound(&F3, &F4, k);
					if (exact_spr >= 0) {
						sync_twins(&F1, &trees[i]);
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
	return 0;
}


/* rSPR_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_3_approx(Forest *T1, Forest *T2) {
	// find sibling pairs of T1
	// match up nodes of T1 and T2
	sync_twins(T1, T2);
	// find singletons of T2
	deque<Node *> sibling_pairs = T1->find_sibling_pairs();
	list<Node *> singletons = T2->find_singletons();
	return rSPR_3_approx_hlpr(T1, T2, &singletons, &sibling_pairs);
}

// rSPR_3_approx recursive helper function
int rSPR_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
		deque<Node *> *sibling_pairs) {
	int num_cut = 0;
	while(!singletons->empty() || !sibling_pairs->empty()) {
		// Case 1 - Remove singletons
		while(!singletons->empty()) {

			Node *T2_a = singletons->back();
			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();
			// if this is in the first component of T_2 then
			// it is not really a singleton.
			if (T2_a == T2->get_component(0))
				continue;

			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			// cut the edge above T1_a
			T1_a->cut_parent();
			T1->add_component(T1_a);
			//delete(T1_a);

			Node *node = T1_a_parent->contract();
			if (potential_new_sibling_pair && node->is_sibling_pair()){
				sibling_pairs->push_front(node->lchild());
				sibling_pairs->push_front(node->rchild());
			}

		}
		if(!sibling_pairs->empty()) {
			Node *T1_a = sibling_pairs->back();
			sibling_pairs->pop_back();
			Node *T1_c = sibling_pairs->back();
			sibling_pairs->pop_back();
			if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
			}
			Node *T1_ac = T1_a->parent();
			// lookup in T2 and determine the case
			Node *T2_a = T1_a->get_twin();
			Node *T2_c = T1_c->get_twin();

			// Case 2 - Contract identical sibling pair
			if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
				Node *T2_ac = T2_a->parent();
				T1_ac->contract_sibling_pair();
				T2_ac->contract_sibling_pair();
				T1_ac->set_twin(T2_ac);
				T2_ac->set_twin(T1_ac);
				T1->add_deleted_node(T1_a);
				T1->add_deleted_node(T1_c);
				T2->add_deleted_node(T2_a);
				T2->add_deleted_node(T2_c);

				// check if T2_ac is a singleton
				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
					singletons->push_back(T2_ac);
				// check if T1_ac is part of a sibling pair
				if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
					sibling_pairs->push_back(T1_ac->parent()->lchild());
					sibling_pairs->push_back(T1_ac->parent()->rchild());
				}
			}
			// Case 3
			else {
				
				//  ensure T2_a is below T2_c
				if (T2_a->get_depth() < T2_c->get_depth()) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
				}
				else if (T2_a->get_depth() == T2_c->get_depth()) {
					if (T2_a->parent() && T2_c->parent() &&
							(T2_a->parent()->get_depth() <
							T2_c->parent()->get_depth())) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
					}
				}

				// get T2_b
				Node *T2_ab = T2_a->parent();
				Node *T2_b = T2_ab->rchild();
				if (T2_b == T2_a)
					T2_b = T2_ab->lchild();
				// cut T1_a, T1_c, T2_a, T2_b, T2_c

				bool cut_b_only = false;
				if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
					cut_b_only = true;
					sibling_pairs->push_back(T1_a);
					sibling_pairs->push_back(T1_c);
				}

				if (!cut_b_only) {
					T1_a->cut_parent();
					T1_c->cut_parent();
					// contract parents
					Node *node = T1_ac->contract();
					// check for T1_ac sibling pair
					if (node && node->is_sibling_pair()){
						sibling_pairs->push_back(node->lchild());
						sibling_pairs->push_back(node->rchild());
					}
				}

				bool same_component = true;
				if (APPROX_CHECK_COMPONENT)
					same_component = (T2_a->find_root() == T2_c->find_root());

				if (!cut_b_only) {
					T2_a->cut_parent();
					num_cut++;
				}
				bool cut_b = false;
				if (same_component && T2_ab->parent() != NULL) {
					T2_b->cut_parent();
					num_cut++;
					cut_b = true;
				}
				// T2_b will move up after contraction
				else if (T2_ab->parent() == NULL) {
					T2_b = T2_b->parent();
				}
				// check for T2 parents as singletons
				Node *node = T2_ab->contract();
				if (node != NULL && node->is_singleton()
						&& node != T2->get_component(0))
					singletons->push_back(node);

				// if T2_c is gone then its replacement is in singleton list
				// contract might delete old T2_c, see where it is
				bool add_T2_c = true;
				T2_c = T1_c->get_twin();
				// ignore T2_c if it is a singleton
				if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

					Node *T2_c_parent = T2_c->parent();
					T2_c->cut_parent();
					num_cut++;
					node = T2_c_parent->contract();
					if (node != NULL && node->is_singleton()
							&& node != T2->get_component(0))
						singletons->push_back(node);
				}
				else {
					add_T2_c = false;
				}

				
				if (!cut_b_only)
					T1->add_component(T1_a);
				if (!cut_b_only)
					T1->add_component(T1_c);
				// put T2 cut parts into T2
				if (!cut_b_only)
					T2->add_component(T2_a);
				// may have already been added
				if (cut_b)
					T2->add_component(T2_b);
				// problem if c is deleted
				if (add_T2_c)
					T2->add_component(T2_c);

				// may have already been added
				if (T2_b->is_leaf())
					singletons->push_back(T2_b);

			}
		}
	}
		// if the first component of the forests differ then we have to cut p
		if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
			num_cut++;
			T1->add_component(new Node("p"));
			T2->add_component(new Node("p"));
		}
		return num_cut;
}

/* rSPR_worse_3_approx
 * Calculate an approximate maximum agreement forest and SPR distance
 * RETURN At most 3 times the rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_worse_3_approx(Forest *T1, Forest *T2) {
	// match up nodes of T1 and T2
	cout << "foo" << endl;
	sync_twins(T1, T2);
	cout << "foo" << endl;
	// find sibling pairs of T1
	deque<Node *> sibling_pairs = T1->find_sibling_pairs();
	cout << "foo" << endl;
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	cout << "singletons:" << endl;
	for(list<Node *>::iterator i = singletons.begin(); i != singletons.end(); i++) {
		cout << " ";
		(*i)->print_subtree();
	}
	return rSPR_worse_3_approx_hlpr(T1, T2, &singletons, &sibling_pairs);
}

int rSPR_worse_3_approx(Forest *T1, Forest *T2, bool sync) {
	// match up nodes of T1 and T2
	if (sync) {
		sync_twins(T1, T2);
	}
	// find sibling pairs of T1
	deque<Node *> sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();

	return rSPR_worse_3_approx_hlpr(T1, T2, &singletons, &sibling_pairs);
}

// rSPR_worse_3_approx recursive helper function
int rSPR_worse_3_approx_hlpr(Forest *T1, Forest *T2, list<Node *> *singletons,
		deque<Node *> *sibling_pairs) {
	//return rSPR_3_approx_hlpr(T1, T2, singletons, sibling_pairs);
	#ifdef DEBUG
		cout << "rSPR_worse_3_approx_hlpr" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
	#endif
	int num_cut = 0;
	while(!singletons->empty() || !sibling_pairs->empty()) {
		// Case 1 - Remove singletons
		while(!singletons->empty()) {
			#ifdef DEBUG
				cout << "Case 1" << endl;
			#endif

			Node *T2_a = singletons->back();
			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();
			// if this is in the first component of T_2 then
			// it is not really a singleton.
	// TODO: problem when we cluster and have a singleton as the
	//		first comp of T2
	//    NEED TO MODIFY CUTTING?
	// 		HERE AND IN BB?
			if (T2_a == T2->get_component(0))
				continue;

			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			// cut the edge above T1_a
			T1_a->cut_parent();
			T1->add_component(T1_a);
			//delete(T1_a);

			Node *node = T1_a_parent->contract();
			if (potential_new_sibling_pair && node->is_sibling_pair()){
				sibling_pairs->push_front(node->lchild());
				sibling_pairs->push_front(node->rchild());
			}

	#ifdef DEBUG
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
	#endif
		}
		if(!sibling_pairs->empty()) {
			Node *T1_a = sibling_pairs->back();
			sibling_pairs->pop_back();
			Node *T1_c = sibling_pairs->back();
			sibling_pairs->pop_back();
			if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
			}
			Node *T1_ac = T1_a->parent();
			// lookup in T2 and determine the case
			Node *T2_a = T1_a->get_twin();
			Node *T2_c = T1_c->get_twin();

			#ifdef DEBUG
				cout << "Fetching sibling pair" << endl;
				T1_ac->print_subtree();
				cout << "T2_a" << ": ";
				cout << " d=" << T2_a->get_depth() << " ";
				T2_a->print_subtree();
				cout << "T1_c" << ": ";
				T1_c->print_subtree();
				cout << "T2_c" << ": ";
				cout << " d=" << T2_c->get_depth() << " ";
				T2_c->print_subtree();
			#endif

			// Case 2 - Contract identical sibling pair
			if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
				#ifdef DEBUG
					cout << "Case 2" << endl;
					T1->print_components();
					T2->print_components();
				#endif
				Node *T2_ac = T2_a->parent();
				T1_ac->contract_sibling_pair();
				T2_ac->contract_sibling_pair();
				T1_ac->set_twin(T2_ac);
				T2_ac->set_twin(T1_ac);
				T1->add_deleted_node(T1_a);
				T1->add_deleted_node(T1_c);
				T2->add_deleted_node(T2_a);
				T2->add_deleted_node(T2_c);

				// check if T2_ac is a singleton
				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
					singletons->push_back(T2_ac);
				// check if T1_ac is part of a sibling pair
				if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
					sibling_pairs->push_back(T1_ac->parent()->lchild());
					sibling_pairs->push_back(T1_ac->parent()->rchild());
				}
			}
			// Case 3
			else {
				#ifdef DEBUG
					cout << "Case 3" << endl;
				#endif
				
				//  ensure T2_a is below T2_c
				if ((T2_a->get_depth() < T2_c->get_depth()
						&& T2_c->parent() != NULL)
						|| T2_a->parent() == NULL) {
					#ifdef DEBUG
						cout << "swapping" << endl;
					#endif
				
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
				}
				else if (T2_a->get_depth() == T2_c->get_depth()) {
					if (T2_a->parent() && T2_c->parent() &&
							(T2_a->parent()->get_depth() <
							T2_c->parent()->get_depth())) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
					}
				}

				// get T2_b
				Node *T2_ab = T2_a->parent();
				Node *T2_b = T2_ab->rchild();
				#ifdef DEBUG
				cout << "T2_b" << ": ";
				cout.flush();
				T2_b->print_subtree();
			#endif
				if (T2_b == T2_a)
					T2_b = T2_ab->lchild();
				// cut T1_a, T1_c, T2_a, T2_b, T2_c

				bool cut_b_only = false;
				if (T2_a->parent() != NULL && T2_a->parent()->parent() != NULL && T2_a->parent()->parent() == T2_c->parent()) {
					cut_b_only = true;
					sibling_pairs->push_back(T1_a);
					sibling_pairs->push_back(T1_c);
				}

				if (!cut_b_only) {
					T1_a->cut_parent();
					T1_c->cut_parent();
					// contract parents
					Node *node = T1_ac->contract();
					// check for T1_ac sibling pair
					if (node && node->is_sibling_pair()){
						sibling_pairs->push_back(node->lchild());
						sibling_pairs->push_back(node->rchild());
					}
				}

				bool same_component = true;
				if (APPROX_CHECK_COMPONENT)
					same_component = (T2_a->find_root() == T2_c->find_root());

				if (!cut_b_only) {
					T2_a->cut_parent();
				}
				bool cut_b = false;
				if (same_component && T2_ab->parent() != NULL) {
					T2_b->cut_parent();
					cut_b = true;
				}
				// T2_b will move up after contraction
				else if (T2_ab->parent() == NULL) {
					T2_b = T2_b->parent();
				}
				// check for T2 parents as singletons
				Node *node = T2_ab->contract();
				if (node != NULL && node->is_singleton()
						&& node != T2->get_component(0))
					singletons->push_back(node);

				// if T2_c is gone then its replacement is in singleton list
				// contract might delete old T2_c, see where it is
				bool add_T2_c = true;
				T2_c = T1_c->get_twin();
				// ignore T2_c if it is a singleton
				if (T2_c != node && T2_c->parent() != NULL && !cut_b_only) {

					Node *T2_c_parent = T2_c->parent();
					T2_c->cut_parent();
					node = T2_c_parent->contract();
					if (node != NULL && node->is_singleton()
							&& node != T2->get_component(0))
						singletons->push_back(node);
				}
				else {
					add_T2_c = false;
				}

				
				if (!cut_b_only)
					T1->add_component(T1_a);
				if (!cut_b_only)
					T1->add_component(T1_c);
				// put T2 cut parts into T2
				if (!cut_b_only)
					T2->add_component(T2_a);
				// may have already been added
				if (cut_b)
					T2->add_component(T2_b);
				// problem if c is deleted
				if (add_T2_c)
					T2->add_component(T2_c);

				// may have already been added
				if (T2_b->is_leaf())
					singletons->push_back(T2_b);

					num_cut+=3;

			}
		}
	}
		// if the first component of the forests differ then we have to cut p
		if (T1->get_component(0)->get_twin() != T2->get_component(0)) {
			num_cut+=3;
			T1->add_component(new Node("p"));
			T2->add_component(new Node("p"));
		}
		return num_cut;
}

/* rSPR_FPT
 * Calculate a maximum agreement forest and SPR distance
 * RETURN The rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_FPT(Forest *T1, Forest *T2, int k) {
	#ifdef DEBUG
		cout << "rSPR_FPT()" << endl;
	#endif
	// find sibling pairs of T1
	sync_twins(T1, T2);
	deque<Node *> sibling_pairs = T1->find_sibling_pairs();
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	int final_k = 
		rSPR_FPT_hlpr(T1, T2, k, &sibling_pairs, &singletons, false);
	if (final_k >= 0)
		final_k = k - final_k;
	return final_k;
}

// recursive helper method for rSPR_FPT
int rSPR_FPT_hlpr(Forest *T1, Forest *T2, int k, deque<Node *> *sibling_pairs,
		list<Node *> *singletons, bool cut_b_only) {
		#ifdef DEBUG
			cout << "K=" << k << endl;
		cout << "sibling pairs:";
		for (deque<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
			cout << "  ";
			(*i)->print_subtree_hlpr();
		}
		cout << endl;
		#endif
	
	while(true) {//!sibling_pairs->empty() || !singletons->empty()) {
		// Case 1 - Remove singletons
		while(!singletons->empty()) {
			#ifdef DEBUG
				cout << "Case 1" << endl;
			#endif

			Node *T2_a = singletons->back();
			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();

			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			// cut the edge above T1_a
			if (T2_a == T2->get_component(0)) {
				T1->add_rho();
				T2->add_rho();
				k--;
			}
			T1_a->cut_parent();
			T1->add_component(T1_a);

			Node *node = T1_a_parent->contract();
			if (potential_new_sibling_pair && node->is_sibling_pair()){
				sibling_pairs->push_front(node->lchild());
				sibling_pairs->push_front(node->rchild());
			}
		}
		if(!sibling_pairs->empty()) {
			Node *T1_a = sibling_pairs->back();
			sibling_pairs->pop_back();
			Node *T1_c = sibling_pairs->back();
			sibling_pairs->pop_back();
			if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				cut_b_only = false;
				continue;
			}
			Node *T1_ac = T1_a->parent();
			// lookup in T2 and determine the case
			Node *T2_a = T1_a->get_twin();
			Node *T2_c = T1_c->get_twin();

			// Case 2 - Contract identical sibling pair
			if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
				#ifdef DEBUG
					cout << "Case 2" << endl;
					T1_ac->print_subtree();
				#endif
				Node *T2_ac = T2_a->parent();
				T1_ac->contract_sibling_pair();
				T2_ac->contract_sibling_pair();
				T1_ac->set_twin(T2_ac);
				T2_ac->set_twin(T1_ac);
				T1->add_deleted_node(T1_a);
				T1->add_deleted_node(T1_c);
				T2->add_deleted_node(T2_a);
				T2->add_deleted_node(T2_c);

				// check if T2_ac is a singleton
				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
					singletons->push_back(T2_ac);
				// check if T1_ac is part of a sibling pair
				if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
					sibling_pairs->push_back(T1_ac->parent()->lchild());
					sibling_pairs->push_back(T1_ac->parent()->rchild());
				}
				cut_b_only = false;
			}
			/* need to copy trees and lists for branching
			 * use forest copy constructor for T1 and T2 giving T1' and T2'
			 * T1' twins are in T2, and same for T2' and T1.
			 * singleton list will be empty except for maybe above the cut,
			 * so this can be created.
			 * fix one set of twins (T2->T1' or T1->T2' not sure)
			 * exploit chained twin relationship to copy sibling pair list
			 * fix other set of twins
			 * swap T2 and T2' root nodes
			 * now do the cut
			 *
			 * note: don't copy for 3rd cut, is a waste
			 */

			//Note: only do BRANCH 2 if we have
			//just done a branch 2
			// need to clear it when we do a Case 2

			// Case 3
			// note: guaranteed that singleton list is empty
			else {
				if (k <= 0)
					return k-1;
				Forest *best_T1;
				Forest *best_T2;
				int best_answer = -1;
				int answer_a = -1;
				int answer_b = -1;
				int answer_c = -1;
				bool cut_ab_only = false;
				//  ensure T2_a is below T2_c
				if (T2_a->get_depth() < T2_c->get_depth()) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
				}
				else if (T2_a->get_depth() == T2_c->get_depth()) {
					if (T2_a->parent() && T2_c->parent() &&
							(T2_a->parent()->get_depth() <
							T2_c->parent()->get_depth())) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
					}
				}
				// get T2_b
				Node *T2_b = T2_a->parent()->rchild();
				if (T2_b == T2_a)
					T2_b = T2_a->parent()->lchild();

			if (CUT_ONE_B) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL)
					cut_b_only=true;
			}
			if (CUT_ONE_AB) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL)
					cut_ab_only=true;
			}
				#ifdef DEBUG
					cout << "Case 3" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "K=" << k << endl;
					cout << "sibling pairs:";
					for (deque<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif


				// copy elements
				Forest *T1_copy;
				Forest *T2_copy;
				deque<Node *> *sibling_pairs_copy;
				Node *T1_a_copy;
				Node *T1_c_copy;
				Node *T2_a_copy;
				Node *T2_c_copy;

				// make copies for the branching
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);


				// cut T2_a
				Node *T2_ab = T2_a->parent();
				T2_a->cut_parent();
				Node *node = T2_ab->contract();
				if (node != NULL && node->is_singleton()
						&& node != T2->get_component(0))
					singletons->push_back(node);
				singletons->push_back(T2_a);
				T2->add_component(T2_a);
				if (cut_b_only == false)
					answer_a =
						rSPR_FPT_hlpr(T1, T2, k-1, sibling_pairs, singletons, false);
					best_answer = answer_a;
					best_T1 = T1;
					best_T2 = T2;



				//load the copy
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();


				// make copies for the branching
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);


				// get T2_b
				T2_ab = T2_a->parent();
				T2_b = T2_ab->rchild();

				if (T2_b == T2_a)
					T2_b = T2_ab->lchild();

				if (!CUT_AC_SEPARATE_COMPONENTS ||
						T2_a->find_root() == T2_c->find_root()) {
					// cut T2_b
					T2_b->cut_parent();
					node = T2_ab->contract();
					if (node != NULL && node->is_singleton()
							&& node != T2->get_component(0))
						singletons->push_back(node);
					T2->add_component(T2_b);
					if (T2_b->is_leaf())
						singletons->push_back(T2_b);
					sibling_pairs->push_back(T1_a);
					sibling_pairs->push_back(T1_c);
					if (CUT_ALL_B) {
						answer_b =
							rSPR_FPT_hlpr(T1, T2, k-1, sibling_pairs,
							singletons, true);
					}
					else {
						answer_b =
							rSPR_FPT_hlpr(T1, T2, k-1, sibling_pairs,
							singletons, false);
					}
				}
				if (answer_b > best_answer 
						|| (answer_b == best_answer
							&& PREFER_RHO
							&& T2->contains_rho() )) {
					best_answer = answer_b;
					swap(&best_T1, &T1);
					swap(&best_T2, &T2);
				}


				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;

				//load the copy
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();

				if (T2_c->parent() != NULL) {
					Node *T2_c_parent = T2_c->parent();
					T2_c->cut_parent();
					node = T2_c_parent->contract();
					if (node != NULL && node->is_singleton()
							&& node != T2->get_component(0))
						singletons->push_back(node);
					T2->add_component(T2_c);
				}
				else {
					// don't increase k
					k++;
				}
				singletons->push_back(T2_c);
				if (cut_b_only == false && cut_ab_only == false)
					answer_c =
						rSPR_FPT_hlpr(T1, T2, k-1, sibling_pairs, singletons, false);
				if (answer_c > best_answer 
						|| (answer_c == best_answer
							&& PREFER_RHO
							&& T2->contains_rho() )) {
					best_answer = answer_c;
					swap(&best_T1, &T1);
					swap(&best_T2, &T2);
				}
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;

				T1 = best_T1;
				T2 = best_T2;
				return best_answer;
			}
		}
		else {
			break;
		}
	}
	return k;
}

int rSPR_branch_and_bound(Forest *T1, Forest *T2) {
	return rSPR_branch_and_bound_range(T1, T2, MAX_SPR);
}


int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int end_k) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	int approx_spr = rSPR_worse_3_approx(&F1, &F2);
	int min_spr = approx_spr / 3;
	return rSPR_branch_and_bound_range(T1, T2, min_spr, end_k);
}
	
int rSPR_branch_and_bound_range(Forest *T1, Forest *T2, int start_k,
		int end_k) {
	int exact_spr = -1;
	for(int k = start_k; k <= end_k; k++) {
		Forest F1 = Forest(T1);
		F1.print_components();
		Forest F2 = Forest(T2);
		exact_spr = rSPR_branch_and_bound(&F1, &F2, k);
		if (exact_spr >= 0 || k == end_k) {
			F1.swap(T1);
			F2.swap(T2);
			break;
		}
	}
	return exact_spr;
}

/* rSPR_branch_and_bound
 * Calculate a maximum agreement forest and SPR distance
 * Uses a branch and bound optimization to not explore paths
 * guaranteed to be incorrect based on rspr_3_approx
 * RETURN The rSPR distance
 * NOTE: destructive. The computed forests replace T1 and T2.
 */
int rSPR_branch_and_bound(Forest *T1, Forest *T2, int k) {
	// find sibling pairs of T1
	cout << "foo1" << endl;
	sync_twins(T1, T2);
	cout << "foo2" << endl;
	deque<Node *> sibling_pairs = T1->find_sibling_pairs();
	cout << "foo3" << endl;
	// find singletons of T2
	list<Node *> singletons = T2->find_singletons();
	cout << "foo4" << endl;
	int final_k = 
		rSPR_branch_and_bound_hlpr(T1, T2, k, &sibling_pairs, &singletons, false);
		cout << "foo" << endl;
	if (final_k >= 0)
		final_k = k - final_k;
	return final_k;
}

// rSPR_branch_and_bound recursive helper function
int rSPR_branch_and_bound_hlpr(Forest *T1, Forest *T2, int k,
		deque<Node *> *sibling_pairs, list<Node *> *singletons,
		bool cut_b_only) {
//		#ifdef DEBUG
		cout << "rSPR_branch_and_bound_hlpr()" << endl;
		cout << "K=" << k << endl;
		cout << "sibling pairs:";
		for (deque<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
			cout << "  ";
			(*i)->print_subtree_hlpr();
		}
		cout << endl;
//		#endif

		if (CLUSTER_REDUCTION) {
			// clean up singletons
			// TODO: this is duplication
		while(!singletons->empty()) {
			Node *T2_a = singletons->back();
			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();
			// if this is in the first component of T_2 then
			// it is not really a singleton.
			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			if (T2_a == T2->get_component(0)) {
				T1->add_rho();
				T2->add_rho();
				k--;
			}

			// cut the edge above T1_a
			T1_a->cut_parent();
			T1->add_component(T1_a);
			Node *node = T1_a_parent->contract();
			if (potential_new_sibling_pair && node->is_sibling_pair()){
				sibling_pairs->push_front(node->lchild());
				sibling_pairs->push_front(node->rchild());
			}
		}
//			cout << "foo" << endl;
//			cout << "foo2" << endl;
			sync_interior_twins_real(T1, T2);
			list<Node *> *cluster_points = find_cluster_points(T1);

			// TODO: could this be faster by using the approx to allocate
			// a certain amount of the k to different clusters?
			// TODO: write pseudocode for what we need
			// TODO: then implement it
			// NOTE: need to make a list of ClusterInstances and then
			// solve each.
			// TODO: where should we do this? Just before we would
			// normally branch?

			cout << "k=" << k << endl;
			cout << "cp=" << cluster_points->size() << endl;
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
			if (!cluster_points->empty()) {
				cout << "CLUSTERS" << endl;
				for(int j = 0; j < 70; j++) {
					cout << "*";
				}
				cout << endl;
				for(list<Node *>::iterator i = cluster_points->begin();
						i != cluster_points->end(); i++) {
					cout << (*i)->str_subtree() << endl;
					cout << (*i)->get_twin()->str_subtree() << endl;
					for(int j = 0; j < 70; j++) {
						cout << "*";
					}
					cout << endl;
				}
				cout << endl;

				list<ClusterInstance> clusters =
					cluster_reduction(T1, T2, cluster_points);

				// TODO: make it so we don't need this?
				T1->unsync_interior();
				T2->unsync_interior();
				while(!clusters.empty()) {
					ClusterInstance cluster = clusters.front();
					clusters.pop_front();
					cluster.F1->unsync_interior();
					cluster.F2->unsync_interior();
					cout << "CLUSTER_START" << endl;
//					cout << &(*cluster.F1->get_component(0)) << endl;
//					cout << &(*clusters.front().F1->get_component(0)) << endl;
					cout << "\tF1: ";
					cluster.F1->print_components_with_twins();
					cout << "\tF2: ";
					cluster.F2->print_components_with_twins();
					int cluster_spr = -1;
					if (k >= 0)
						cluster_spr = rSPR_branch_and_bound_range(cluster.F1,
								cluster.F2, k);
					cout << "foo" << endl;
					if (cluster_spr >= 0) {
						cout << "cluster k=" << cluster_spr << endl;
						cout << "\tF1: ";
						cluster.F1->print_components();
						cout << "\tF2: ";
						cluster.F2->print_components();
						k -= cluster_spr;
					}
					else {
						k = -1;
					}
					cout << "\tF1: ";
					T1->print_components();
					cout << "\tF2: ";
					T2->print_components();
					if (!cluster.is_original()) {
						cluster.join_cluster(T1, T2);

//						cout << cluster.F1_cluster_node->str_subtree() << endl;
//						cout << cluster.F2_cluster_node->str_subtree() << endl;
//						Node *p = cluster.F1_cluster_node;
//						while (p->parent() != NULL)
//							p = p->parent();
//						cout << &(*p) << endl;
//						cout << p->str_subtree() << endl;


						cout << "\tF1: ";
						T1->print_components();
						cout << "\tF2: ";
						T2->print_components();
					}
					else {
						cout << "original" << endl;
					}
				}
				delete cluster_points;
				cout << "returning k=" << k << endl;
				return k;
			}
			else {
				T1->unsync_interior();
				T2->unsync_interior();
			}
			delete cluster_points;

//			cout << "done" << endl;
		}
	
	while(!singletons->empty() || !sibling_pairs->empty()) {
		// Case 1 - Remove singletons
		while(!singletons->empty()) {
			Node *T2_a = singletons->back();
			#ifdef DEBUG
				cout << "Case 1" << endl;
				cout << "\tT1: ";
				T1->print_components();
				cout << "\tT2: ";
				T2->print_components();
				cout << "a " << T2_a->str() << endl;
			#endif

			singletons->pop_back();
			// find twin in T1
			Node *T1_a = T2_a->get_twin();
			// if this is in the first component of T_2 then
			// it is not really a singleton.
			Node *T1_a_parent = T1_a->parent();
			if (T1_a_parent == NULL)
				continue;
			bool potential_new_sibling_pair = T1_a_parent->is_sibling_pair();
			if (T2_a == T2->get_component(0)) {
				T1->add_rho();
				T2->add_rho();
				k--;
				#ifdef DEBUG
				cout << "adding p element, k=" << k << endl;
				#endif
			}

			// cut the edge above T1_a
			T1_a->cut_parent();
			T1->add_component(T1_a);
			Node *node = T1_a_parent->contract();
			if (potential_new_sibling_pair && node->is_sibling_pair()){
				sibling_pairs->push_front(node->lchild());
				sibling_pairs->push_front(node->rchild());
			}
		}
		if(!sibling_pairs->empty()) {
			Node *T1_a = sibling_pairs->back();
			sibling_pairs->pop_back();
			Node *T1_c = sibling_pairs->back();
			sibling_pairs->pop_back();
			if (T1_a->parent() == NULL || T1_a->parent() != T1_c->parent()) {
				continue;
			}
			Node *T1_ac = T1_a->parent();
			// lookup in T2 and determine the case
			Node *T2_a = T1_a->get_twin();
			Node *T2_c = T1_c->get_twin();

			if (T2_a->parent() != NULL && T2_a->parent() == T2_c->parent()) {
				#ifdef DEBUG
					cout << "Case 2" << endl;
					T1_ac->print_subtree();
				#endif
				Node *T2_ac = T2_a->parent();
				T1_ac->contract_sibling_pair();
				T2_ac->contract_sibling_pair();
				T1_ac->set_twin(T2_ac);
				T2_ac->set_twin(T1_ac);
				T1->add_deleted_node(T1_a);
				T1->add_deleted_node(T1_c);
				T2->add_deleted_node(T2_a);
				T2->add_deleted_node(T2_c);

				// check if T2_ac is a singleton
				if (T2_ac->is_singleton() && !T1_ac->is_singleton() && T2_ac != T2->get_component(0))
					singletons->push_back(T2_ac);
				// check if T1_ac is part of a sibling pair
				if (T1_ac->parent() != NULL && T1_ac->parent()->is_sibling_pair()) {
					sibling_pairs->push_back(T1_ac->parent()->lchild());
					sibling_pairs->push_back(T1_ac->parent()->rchild());
				}
			}
			/* need to copy trees and lists for branching
			 * use forest copy constructor for T1 and T2 giving T1' and T2'
			 * T1' twins are in T2, and same for T2' and T1.
			 * singleton list will be empty except for maybe above the cut,
			 * so this can be created.
			 * fix one set of twins (T2->T1' or T1->T2' not sure)
			 * exploit chained twin relationship to copy sibling pair list
			 * fix other set of twins
			 * swap T2 and T2' root nodes
			 * now do the cut
			 *
			 * note: don't copy for 3rd cut, is a waste
			 */

			// Case 3
			// note: guaranteed that singleton list is empty
			else {
				if (k <= 0) {
					return k-1;
				}
				Forest *best_T1;
				Forest *best_T2;
				int best_k = -1;
				int answer_a = -1;
				int answer_b = -1;
				int answer_c = -1;
				bool cut_ab_only = false;
				//  ensure T2_a is below T2_c
				if ((T2_a->get_depth() < T2_c->get_depth()
						&& T2_c->parent() != NULL)
						|| T2_a->parent() == NULL) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
				}
				else if (T2_a->get_depth() == T2_c->get_depth()) {
					if (T2_a->parent() && T2_c->parent() &&
							(T2_a->parent()->get_depth() <
							T2_c->parent()->get_depth())) {
					swap(&T1_a, &T1_c);
					swap(&T2_a, &T2_c);
					}
				}
				Node *T2_b = T2_a->parent()->rchild();
				if (T2_b == T2_a)
					T2_b = T2_a->parent()->lchild();

			if (CUT_ONE_B) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL)
					cut_b_only=true;
			}
			if (CUT_ONE_AB) {
				if (T2_a->parent()->parent() == T2_c->parent()
					&& T2_c->parent() != NULL)
					cut_ab_only=true;
			}
				#ifdef DEBUG
					cout << "Case 3" << endl;
					cout << "\tT1: ";
					T1->print_components();
					cout << "\tT2: ";
					T2->print_components();
					cout << "K=" << k << endl;
					cout << "sibling pairs:";
					for (deque<Node *>::iterator i = sibling_pairs->begin(); i != sibling_pairs->end(); i++) {
						cout << "  ";
						(*i)->print_subtree_hlpr();
					}
					cout << endl;
					cout << "\tcut_b_only=" << cut_b_only << endl;
					cout << "\tT2_a " << T2_a->str() << " "
						<< T2_a->get_depth() << endl;
					cout << "\tT2_c " << T2_c->str() << " "
						<< T2_c->get_depth() << endl;
					cout << "\tT2_b " << T2_b->str_subtree() << " "
						<< T2_b->get_depth() << endl;
				#endif
			


				// copy elements
				Forest *T1_copy;
				Forest *T2_copy;
				deque<Node *> *sibling_pairs_copy;
				Node *T1_a_copy;
				Node *T1_c_copy;
				Node *T2_a_copy;
				Node *T2_c_copy;
				//list<Node *> *singletons_copy = new list<Node *>();

				// make copies for the approx
				// be careful we do not kill real T1 and T2
				// ie use the copies
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);
				int approx_spr = rSPR_worse_3_approx_hlpr(T1_copy, T2_copy,
						singletons, sibling_pairs_copy);
				delete T1_copy;
				delete T2_copy;
				delete sibling_pairs_copy;
				if (approx_spr  >  3*k){
					return -1;
				}
				

				// make copies for the branching
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);

				// cut T2_a
				Node *T2_ab = T2_a->parent();
				T2_a->cut_parent();
				Node *node = T2_ab->contract();
				if (node != NULL && node->is_singleton() &&
						node != T2->get_component(0))
					singletons->push_back(node);
				singletons->push_back(T2_a);
				T2->add_component(T2_a);
				if (cut_b_only == false)
					answer_a =
						rSPR_branch_and_bound_hlpr(T1, T2, k-1,
								sibling_pairs, singletons, false);
				best_k = answer_a;
				best_T1 = T1;
				best_T2 = T2;



				//load the copy
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();


				// make copies for the branching
				copy_trees(&T1, &T2, &sibling_pairs, &T1_a, &T1_c, &T2_a, &T2_c,
						&T1_copy, &T2_copy, &sibling_pairs_copy,
						&T1_a_copy, &T1_c_copy, &T2_a_copy, &T2_c_copy);

				// get T2_b
				T2_ab = T2_a->parent();
				T2_b = T2_ab->rchild();

				if (T2_b == T2_a)
					T2_b = T2_ab->lchild();

				if (!CUT_AC_SEPARATE_COMPONENTS || T2_a->find_root() == T2_c->find_root()) {
					// cut T2_b
					T2_b->cut_parent();
					node = T2_ab->contract();
					if (node != NULL && node->is_singleton()
							&& node != T2->get_component(0))
						singletons->push_back(node);
					T2->add_component(T2_b);
					if (T2_b->is_leaf())
						singletons->push_back(T2_b);
					sibling_pairs->push_back(T1_a);
					sibling_pairs->push_back(T1_c);
					if (CUT_ALL_B) {
						answer_b =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1,
									sibling_pairs, singletons, true);
					}
					else {
						answer_b =
							rSPR_branch_and_bound_hlpr(T1, T2, k-1,
									sibling_pairs, singletons, false);
					}
				}
				if (answer_b > best_k
						|| (answer_b == best_k
							&& PREFER_RHO
							&& T2->contains_rho() )) {
					best_k = answer_b;
					swap(&best_T1, &T1);
					swap(&best_T2, &T2);
				}
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;

				//load the copy
				T1 = T1_copy;
				T2 = T2_copy;
				T1_a = T1_a_copy;
				T1_c = T1_c_copy;
				T2_a = T2_a_copy;
				T2_c = T2_c_copy;
				sibling_pairs = sibling_pairs_copy;
				singletons = new list<Node *>();

				if (T2_c->parent() != NULL) {
					Node *T2_c_parent = T2_c->parent();
					T2_c->cut_parent();
					node = T2_c_parent->contract();
					if (node != NULL && node->is_singleton()
							&& node != T2->get_component(0))
						singletons->push_back(node);
					T2->add_component(T2_c);
				}
				else {
					// don't increase k
					k++;
				}
				singletons->push_back(T2_c);
				if (cut_b_only == false && cut_ab_only == false)
					answer_c =
						rSPR_branch_and_bound_hlpr(T1, T2, k-1,
								sibling_pairs, singletons, false);
				if (answer_c > best_k
						|| (answer_c == best_k
							&& PREFER_RHO
							&& T2->contains_rho() )) {
					best_k = answer_c;
					swap(&best_T1, &T1);
					swap(&best_T2, &T2);
				}
				delete T1;
				delete T2;
				delete sibling_pairs;
				delete singletons;

				T1 = best_T1;
				T2 = best_T2;

				return best_k;
			}
			cut_b_only = false;
		}
	}

	return k;
}

/* copy_trees
 * copies T1 and T2 into T1_copy and T2_copy
 * to maintain twins among the trees, T2 and T2_copy have
 * to switch components so pointers within the trees have
 * to be updated
 */
inline void copy_trees(Forest **T1, Forest **T2, deque<Node *> **real_sibling_pairs,
		Node **T1_a, Node **T1_c, Node **T2_a, Node **T2_c,
		Forest **T1_copy, Forest **T2_copy, deque<Node *> **real_sibling_pairs_copy,
		Node **T1_a_copy, Node **T1_c_copy, Node **T2_a_copy, Node **T2_c_copy) {
	(*T1_copy) = new Forest(**T1);
	*T2_copy = new Forest(**T2);
	(*T1_copy)->resync();
	// update sibling pair list to point to new trees
	deque<Node *> *sibling_pairs = *real_sibling_pairs;
	*real_sibling_pairs_copy = new deque<Node *>();
	deque<Node *> *sibling_pairs_copy = *real_sibling_pairs_copy;
	for(int i = 0; i < (sibling_pairs)->size(); i++) {
		Node *node = (*sibling_pairs)[i]->get_twin()->get_twin();
		if (node->parent() == NULL) {
			i++;
			continue;
		}
		(sibling_pairs_copy)->push_back(node);
		i++;
		node = (*sibling_pairs)[i]->get_twin()->get_twin();
		if (node->parent() == NULL) {
			(sibling_pairs_copy)->pop_back();
			continue;
		}
		(sibling_pairs_copy)->push_back(node);
	}
	(*T2_copy)->resync();


	// swap T2_copy with T2
	swap(T2, T2_copy);
	*T1_a_copy = (*T2_a)->get_twin();
	*T1_c_copy = (*T2_c)->get_twin();
	*T2_a_copy = *T2_a;
	*T2_c_copy = *T2_c;
	*T2_a = (*T1_a)->get_twin();
	*T2_c = (*T1_c)->get_twin();
}
