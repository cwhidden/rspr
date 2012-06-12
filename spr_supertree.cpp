/*******************************************************************************
spr_supertree.cpp

Copyright 2011 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/SPR_Supertrees
August 22, 2011

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

"*******************************************************************************\n";

*******************************************************************************

*******************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <list>
#include <time.h>
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
bool SIMPLE_UNROOTED = false;
bool APPROX = false;
bool TIMING = false;
int NUM_ITERATIONS = 25;
bool SMALL_TREES = false;
bool CONVERT_LIST = false;
bool VALID_TREES = false;
bool MULTI_TREES = false;
int NUM_LEAVES=-1;
int APPROX_SIBLINGS = 0;

string USAGE =
"spr_supertrees, version 1.00\n"
"\n"
"usage: spr_supertrees [OPTIONS]\n"
"Calculate binary rooted supertrees that minimize the Subtree Prune and\n"
"Regraft (SPR) distance to a set of rooted or unrooted binary trees from\n"
"STDIN in newick format. Supports arbitrary leaf labels. See the\n"
"README for more information.\n"
"\n"
"Copyright 2011 Chris Whidden\n"
"whidden@cs.dal.ca\n"
"http://kiwi.cs.dal.ca/Software/SPR_Supertrees\n"
"October 28, 2011\n"
"Version 1.00\n"
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
"-max k      Calculate the exact rSPR distance if it is k or less and\n"
"            otherwise use the 3-approximation\n"
"\n"
"-i x        Apply x iterations of Global SPR rearrangements\n"
"\n"
"*******************************************************************************\n"
"UNROOTED COMPARISON OPTIONS\n"
"*******************************************************************************\n"
"\n"
"-unrooted   Compare the supertree to each rooting of the input trees.\n"
"            Use the best found distance.\n"
"\n"
"-unrooted_min_approx   Compare the supertree to each rooting of the\n"
"                       input trees. Run the exact algorithm on the\n"
"                       rooting with the minimum approximate rspr distance\n"
"\n"
"*******************************************************************************\n"
"OTHER OPTIONS\n"
"*******************************************************************************\n"
"-cc             Calculate a potentially better approximation with a quadratic time\n"
"                algorithm\n"
"\n"
"-valid_trees    Output the set of trees that appear valid\n"
"\n"
"-multi_trees    Output the set of multifurcating or invalid trees\n"
"\n"
"*******************************************************************************\n";

Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees,
		int label);
Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees,
		vector<Node *> *best_siblings, int label);
void find_best_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &num_ties,
		Node **best_sibling);
vector<Node *> *find_best_siblings(Node *super_tree, vector<Node *> &gene_trees, int label, int num_siblings);
void find_best_siblings_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &num_ties,
		multimap<int, Node*> *best_siblings, int num_siblings);
void test_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &num_ties,
		Node **best_sibling);
void find_best_spr(Node *super_tree, vector<Node *> &gene_trees,
		Node *&best_spr_move, Node *&best_sibling);
void find_best_spr_helper(Node *n, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties);
void find_best_spr_helper(Node *n, Node *new_sibling, Node *super_tree,
		vector<Node *> &gene_trees, Node *&best_spr_move,
		Node *&best_sibling, int &min_distance, int &num_ties);
Node *find_best_root(Node *T1, Node *T2);
double find_best_root_acc(Node *T1, Node *T2);
void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc);
void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties);

	map<string, int> label_map;
	map<int, string> reverse_label_map;

bool is_pow_2(int n) {
	return (n) && !(n & (n - 1));
}

int main(int argc, char *argv[]) {

	// ignore multifurcating trees by default
	IGNORE_MULTI = true;

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
		else if (strcmp(arg, "-simple_unrooted") == 0)
			SIMPLE_UNROOTED = true;
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
				cout << "CLUSTER_MAX_SPR=" << MAX_SPR << endl;
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
/*
		else if (strcmp(arg, "-small_trees") == 0) {
			SMALL_TREES=true;
		}
*/
		else if (strcmp(arg, "-convert_list") == 0) {
			CONVERT_LIST=true;
		}
		else if (strcmp(arg, "-valid_trees") == 0) {
			VALID_TREES=true;
		}
		else if (strcmp(arg, "-multi_trees") == 0) {
			MULTI_TREES=true;
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
		else if (strcmp(arg, "--help") == 0) {
			cout << USAGE;
			return 0;
		}
			
	}
	if (DEFAULT_OPTIMIZATIONS) {
		CUT_ALL_B=true;
		CUT_ONE_B = true;
		REVERSE_CUT_ONE_B = true;
		CUT_AC_SEPARATE_COMPONENTS = true;
		EDGE_PROTECTION = true;
//		if (ALL_MAFS == false)
//			ABORT_AT_FIRST_SOLUTION = true;
//		PREORDER_SIBLING_PAIRS = true;
		NEAR_PREORDER_SIBLING_PAIRS = true;
		LEAF_REDUCTION = true;
//		LEAF_REDUCTION2 = true;

		APPROX_CUT_ONE_B = true;
		APPROX_CUT_TWO_B = true;
		APPROX_REVERSE_CUT_ONE_B = true;
		APPROX_EDGE_PROTECTION = true;
	}
	PREORDER_SIBLING_PAIRS = true;
	if (DEFAULT_ALGORITHM) {
		BB=true;
	}



	// initialize random number generator
	srand((unsigned(time(0))));

	// Label maps to allow string labels
	vector<int> label_counts = vector<int>();
	label_map= map<string, int>();
	reverse_label_map = map<int, string>();

	string T_line = "";
	vector<Node *> gene_trees = vector<Node *>();
	multimap<int, pair<Node*, string> > gene_tree_map
		= multimap<int, pair<Node*, string> >();
	vector<string> gene_tree_names = vector<string>();
	int skipped_multifurcating = 0;
	int skipped_small = 0;
	int skipped_no_bracket = 0;
	int skipped_star = 0;
	while (getline(cin, T_line)) {
		string name = "";
		size_t loc = T_line.find_first_of("(");
		if (loc != string::npos) {
			if (loc != 0) {
				name = T_line.substr(0,loc);
				T_line.erase(0,loc);
			}
			if (UNROOTED || SIMPLE_UNROOTED) {
				T_line = root(T_line);
			}
			Node *T = build_tree(T_line);
			// TODO: check that this works
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

			if (T->is_leaf() && T->str() == "p") {
				if (MULTI_TREES)
					cout << name << T_line << endl;

				skipped_multifurcating++;
				continue;
			}
			int T_size = T->size();
			if ((T_size <= 4)
					|| T_size == 5 && !SMALL_TREES) {
				skipped_small++;
				continue;
			}
			if (!IGNORE_MULTI) {
				int T_depth = T->max_depth();
				if (T_depth <= 2 ||
						((UNROOTED || SIMPLE_UNROOTED) && T_depth <= 3)) {
					skipped_star++;
					continue;
				//cout << T->str_subtree() << endl;
				//cout << T_depth << endl;
				}
			}
			if (UNROOTED || SIMPLE_UNROOTED)
				T->preorder_number();
			//gene_tree_names.push_back(name);
			//gene_trees.push_back(T);
			gene_tree_map.insert(make_pair(T->size(), make_pair(T, name)));
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

	int end = gene_tree_map.size();
	for(int i = 0; i < end; i++) {
		multimap<int,  pair<Node *, string> >::iterator p =
			gene_tree_map.begin();
		gene_trees.push_back(p->second.first);
		gene_tree_names.push_back(p->second.second);
		gene_tree_map.erase(p);
	}

	cout << gene_trees.size() << " gene trees remaining" << endl;

	if (gene_trees.size() <= 0)
		exit(0);

	for(int i = 0; i < gene_trees.size(); i++) {
		if (VALID_TREES) {
			cout << gene_tree_names[i];
			cout << gene_trees[i]->str_subtree() << endl;
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
		labels.insert(make_pair(label_counts[i],i));
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


//	for(auto label = labels.rbegin(); label != labels.rend(); label++) {
//		cout << label->second ;
//		cout << " (" << reverse_label_map.find(label->second)->second << ") ";
//		cout << ": " << label->first << endl;
//	}

	// 4 most common leaves
	multimap<int, int>::reverse_iterator label = labels.rbegin();
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


	// choose the best starting tree
	int best_distance;
	if (APPROX)
		if (UNROOTED)
			best_distance = rSPR_total_approx_distance_unrooted(ST[0], gene_trees); 
		else
			best_distance = rSPR_total_approx_distance(ST[0], gene_trees); 
	else 
		if (UNROOTED)
			best_distance = rSPR_total_distance_unrooted(ST[0], gene_trees); 
		else
			best_distance = rSPR_total_distance(ST[0], gene_trees); 
	int best_tree = 0;
//	cout << best_distance <<  ": ";
//	cout << ST[0]->str_subtree() << endl;
	for(int j = 1; j < ST.size(); j++) {
		//cout << ST[j]->str_subtree() << endl;
		int distance;
		if (APPROX)
			if (UNROOTED)
				distance = rSPR_total_approx_distance_unrooted(ST[j], gene_trees);
			else
				distance = rSPR_total_approx_distance(ST[j], gene_trees);
		else
			if (UNROOTED || SIMPLE_UNROOTED)
				distance = rSPR_total_distance_unrooted(ST[j], gene_trees);
			else
				distance = rSPR_total_distance(ST[j], gene_trees);
//		cout << distance <<  ": ";
//		cout << ST[j]->str_subtree() << endl;
		if (distance < best_distance) {
			best_distance = distance;
			best_tree = j;
		}
	}
	Node *super_tree = ST[best_tree];

	for(int i = 0; i < ST.size(); i++) {
		if (i != best_tree)
			ST[i]->delete_tree();
	}

	cout << endl;
	cout << "Initial Supertree:  " << super_tree->str_subtree() << endl;
	double time;
	double current_time;
	if (TIMING)
		time = clock()/(double)CLOCKS_PER_SEC;

	int x = 0;
	int leaf_num = 5;
	for(; label != labels.rend() &&
			NUM_LEAVES < 0 || leaf_num <= NUM_LEAVES; label++) {
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
		bool APPROX_ROOTING = false;
		if (SIMPLE_UNROOTED) {
			if (is_pow_2(leaf_num-5)) {
				// reroot the supertree based on the balanced accuracy of splits
				vector<Node *> descendants =
					super_tree->find_interior();
					//super_tree->find_descendants();
				Node *best_root = super_tree->lchild();
				double best_root_avg_acc = 0;
				if (APPROX_ROOTING)
					best_root_avg_acc = -INT_MAX;
				int num_ties = 2;
				for(int j = 0; j < descendants.size(); j++) {
					double root_avg_acc = 0;
					int count = 1;
					super_tree->reroot(descendants[j]);
					super_tree->preorder_number();
					if (APPROX_ROOTING) {
						root_avg_acc = -rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
						//root_avg_acc = -rSPR_total_distance_unrooted(super_tree, current_gene_trees);
					}
					else {
						for(int i = 0; i < current_gene_trees.size(); i++) {
							double acc;
							acc = find_best_root_acc(super_tree, current_gene_trees[i]);
		//					cout <<  j << "\t" << i << "\t" << acc << endl;
							if (acc > -1) {
								root_avg_acc += (acc - root_avg_acc) / count;
								count++;
							}
						}
						int lsize = super_tree->lchild()->size_using_prenum();
						int rsize = super_tree->rchild()->size_using_prenum();
						int size = (lsize < rsize) ? lsize : rsize;
						root_avg_acc *= mylog2(size);
					}
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
				super_tree->reroot(best_root);
			}

			// reroot the gene trees based on the balanced accuracy of splits
			super_tree->preorder_number();
			for(int i = 0; i < current_gene_trees.size(); i++) {
				current_gene_trees[i]->preorder_number();
				Node *new_root =
					find_best_root(super_tree, current_gene_trees[i]);
				if (new_root != NULL)
					current_gene_trees[i]->reroot(new_root);
			}
		}
		Node *best_sibling;
		if (APPROX_SIBLINGS > 0) {
			vector<Node *> *best_siblings = find_best_siblings(super_tree,
					gene_trees, label->second, APPROX_SIBLINGS);
			best_sibling = find_best_sibling(super_tree,
					current_gene_trees, best_siblings, label->second);
			delete best_siblings;
		}
		else {
			best_sibling = find_best_sibling(super_tree,
					current_gene_trees, label->second);
		}


		Node *node = best_sibling->expand_parent_edge(best_sibling);

		node->add_child(new Node(itos(label->second)));
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << super_tree->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
	}

	cout << endl;
	super_tree->numbers_to_labels(&reverse_label_map);
	cout << "Initial Supertree: " <<  super_tree->str_subtree() << endl;
	super_tree->labels_to_numbers(&label_map, &reverse_label_map);
	Node *best_supertree = new Node(*super_tree);
	if (APPROX)
		if (UNROOTED)
			best_distance = rSPR_total_approx_distance_unrooted(super_tree, gene_trees);
		else
			best_distance = rSPR_total_approx_distance(super_tree, gene_trees);
	else
		if (UNROOTED)
			best_distance = rSPR_total_distance_unrooted(super_tree, gene_trees);
		else
			best_distance = rSPR_total_distance(super_tree, gene_trees);
	cout << "Total Distance: " << best_distance << endl;
		if (TIMING) {
			current_time = time;
			time = clock()/(double)CLOCKS_PER_SEC;
			current_time = time - current_time;
			cout << "\t" << current_time << "\t" << time << endl;
		}
	//super_tree->numbers_to_labels(&reverse_label_map);

	if (NUM_ITERATIONS < 0)
		NUM_ITERATIONS=labels.size(); 
	for(int i = 0; i < NUM_ITERATIONS; i++) {
		Node *best_subtree_root;
		Node *best_sibling;
		find_best_spr(super_tree, gene_trees, best_subtree_root, best_sibling);
		best_subtree_root->spr(best_sibling);
		super_tree->numbers_to_labels(&reverse_label_map);
		cout << "Current Supertree: " <<  super_tree->str_subtree() << endl;
		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
//		super_tree->labels_to_numbers(&label_map, &reverse_label_map);
		int current_distance;
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
		if (current_distance < best_distance) {
			best_supertree->delete_tree();
			best_supertree = new Node(*super_tree);
			best_distance = current_distance;
		}
		cout << "Total Distance: " << current_distance << endl;
		if (TIMING) {
			current_time = time;
			time = clock()/(double)CLOCKS_PER_SEC;
			current_time = time - current_time;
			cout << "\t" << current_time << "\t" << time << endl;
		}
	}
	super_tree->delete_tree();
	super_tree=best_supertree;
	super_tree->numbers_to_labels(&reverse_label_map);
	cout << "Final Supertree: " <<  super_tree->str_subtree() << endl;
	cout << "Final Distance: " << best_distance << endl;


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
	int num_ties = 0;
	find_best_sibling_helper(super_tree, new_leaf, super_tree, gene_trees,
			min_distance, num_ties, &best_sibling);
	delete new_leaf;
	return best_sibling;
}

void find_best_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &num_ties,
		Node **best_sibling) {

	if (n->lchild() != NULL) {
		find_best_sibling_helper(n->lchild(), new_leaf, super_tree,
				gene_trees, min_distance, num_ties, best_sibling);
	}
	if (n->rchild() != NULL) {
		find_best_sibling_helper(n->rchild(), new_leaf, super_tree,
				gene_trees, min_distance, num_ties, best_sibling);
	}
	test_sibling_helper(n, new_leaf, super_tree,
			gene_trees, min_distance, num_ties, best_sibling);
}

Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees,
		vector<Node *> *best_siblings, int label) {
	Node *best_sibling;
	Node *new_leaf = new Node(itos(label));
	int min_distance = INT_MAX;
	int num_ties = 0;
	for(int i = 0; i < best_siblings->size(); i++) {
		test_sibling_helper((*best_siblings)[i], new_leaf, super_tree, gene_trees,
				min_distance, num_ties, &best_sibling);
	}
	delete new_leaf;
	return best_sibling;
}

void test_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
		vector<Node *> &gene_trees, int &min_distance, int &num_ties,
		Node **best_sibling) {

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
	if (distance < min_distance) {
		min_distance = distance;
		*best_sibling = n;
		num_ties = 2;
	}
	else if (distance == min_distance) {
		int r = rand();
		if (r < RAND_MAX/num_ties) {
			min_distance = distance;
			*best_sibling = n;
		}
		num_ties++;
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
//	cout << "Reverted: " << super_tree->str_subtree() << endl;

}

vector<Node *> *find_best_siblings(Node *super_tree, vector<Node *> &gene_trees, int label, int num_siblings) {
	Node *new_leaf = new Node(itos(label));
	int min_distance = INT_MAX;
	int num_ties = 0;
	multimap<int, Node*> best_siblings = multimap<int, Node*>();
	find_best_siblings_helper(super_tree, new_leaf, super_tree, gene_trees,
			min_distance, num_ties, &best_siblings, num_siblings);
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
		vector<Node *> &gene_trees, int &min_distance, int &num_ties,
		multimap<int, Node*> *best_siblings, int num_siblings) {

	if (n->lchild() != NULL) {
		find_best_siblings_helper(n->lchild(), new_leaf, super_tree,
				gene_trees, min_distance, num_ties, best_siblings, num_siblings);
	}
	if (n->rchild() != NULL) {
		find_best_siblings_helper(n->rchild(), new_leaf, super_tree,
				gene_trees, min_distance, num_ties, best_siblings, num_siblings);
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
	int distance;
	distance = rSPR_total_approx_distance(super_tree, gene_trees,
			min_distance);
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
//	cout << "Reverted: " << super_tree->str_subtree() << endl;

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
	}

}

Node *find_best_root(Node *T1, Node *T2, double *best_root_b_acc) {
	Forest F1 = Forest(T1);
	Forest F2 = Forest(T2);
	Node *t1 = F1.get_component(0);
	Node *t2 = F2.get_component(0);
	//F1->preorder_number();
	int lchild_pre = t1->lchild()->get_preorder_number();
	int rchild_pre = t1->rchild()->get_preorder_number();
	if (!sync_twins(&F1, &F2)) {
		return NULL;
	}
	if (t1->lchild()->get_preorder_number() != lchild_pre ||
			t1->rchild()->get_preorder_number() != rchild_pre)
		return NULL;
	if (F2.get_component(0)->get_children().size() > 2)
		F2.get_component(0)->fixroot();
	// TODO: maybe stop if F2 is too small?
	int pre_separator = t1->lchild()->get_preorder_number();
	int group_1_total;
	int group_2_total;
	if (t1->rchild()->get_preorder_number() > pre_separator) {
		pre_separator = t1->rchild()->get_preorder_number();
		group_1_total = t1->lchild()->find_leaves().size();
		group_2_total = t1->rchild()->find_leaves().size();
	}
	else {
		group_1_total = t1->rchild()->find_leaves().size();
		group_2_total = t1->lchild()->find_leaves().size();
	}
//	cout << "g1_total: " << group_1_total << endl;
//	cout << "g2_total: " << group_2_total << endl;
	Node *best_root = t2->lchild();
	find_best_root_hlpr(t2, pre_separator, group_1_total, group_2_total,
			&best_root, best_root_b_acc);
//	cout << "t1: " << t1->str_subtree() << endl;
//	cout << "t2: " << t2->str_subtree() << endl;
//	cout << "best_root: " << best_root->str_subtree() << endl;
	best_root = T2->find_by_prenum(best_root->get_preorder_number());
//	cout << "T1: " << T1->str_subtree() << endl;
//	cout << "best_root: " << best_root->str_subtree() << endl;
//	T2->reroot(best_root);
//	cout << "T2: " << T2->str_subtree() << endl;
	return best_root;
}
Node *find_best_root(Node *T1, Node *T2) {
	double best_root_b_acc = 0;
	return find_best_root(T1, T2, &best_root_b_acc);
}

double find_best_root_acc(Node *T1, Node *T2) {
	double best_root_b_acc = -1;
	find_best_root(T1, T2, &best_root_b_acc);
	return best_root_b_acc;
}

void find_best_root_hlpr(Node *T2, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
	int num_ties = 2;
	for(c = T2->get_children().begin(); c != T2->get_children().end(); c++) {
		find_best_root_hlpr(*c, pre_separator, group_1_total,
				group_2_total, best_root, best_root_b_acc,
				&group_1_descendants, &group_2_descendants, &num_ties);
	}
}

void find_best_root_hlpr(Node *n, int pre_separator, int group_1_total,
		int group_2_total, Node **best_root, double *best_root_b_acc,
		int *p_group_1_descendants, int *p_group_2_descendants, int *num_ties) {
	list<Node*>::iterator c;
	int group_1_descendants = 0;
	int group_2_descendants = 0;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		find_best_root_hlpr(*c, pre_separator, group_1_total,
				group_2_total, best_root, best_root_b_acc,
				&group_1_descendants, &group_2_descendants, num_ties);
	}
	if (n->is_leaf()) {
		int pre = n->get_twin()->get_preorder_number();
		if (pre < pre_separator)
			group_1_descendants++;
		else
			group_2_descendants++;
	}
	// balanced accuracy
	// don't bother averaging since we only directly compare them
	double tpos = group_1_descendants;
	double fpos = group_2_descendants;
	double fneg = (group_1_total - group_1_descendants);
	double tneg = (group_2_total - group_2_descendants);
	// balanced accuracy for this bipartition
	double b_acc =  tpos / (tpos + fneg)
			+ tneg / (tneg + fpos);
	// balanced accuracy for the opposite bipartition
	double b_acc_opp = fpos / (fpos + tneg)
			+ fneg / (fneg + tpos);
	// use the better bipartition
//	cout << "n: " << n->str_subtree() << endl;
//	cout << "b_acc: " << b_acc << endl;
//	cout << "b_acc_opp: " << b_acc_opp << endl;
	if (b_acc_opp > b_acc)
		b_acc = b_acc_opp;
	if (b_acc > *best_root_b_acc) {
		*best_root = n;
		*best_root_b_acc = b_acc;
		*num_ties = 2;
	}
	else if(b_acc == *best_root_b_acc) {
		int r = rand();
		if (r < RAND_MAX/ *num_ties) {
		*best_root = n;
		*best_root_b_acc = b_acc;
		*num_ties = 2;
		}
	}
	*p_group_1_descendants += group_1_descendants;
	*p_group_2_descendants += group_2_descendants;
}

