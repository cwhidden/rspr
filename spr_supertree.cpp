/*******************************************************************************
spr_supertree.cpp

Copyright 2011 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
August 22, 2011

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

TODO:

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
#include <unordered_map>
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
bool UNROOTED_MIN_APPROX = false;
bool LCA_TEST = false;
bool CLUSTER_TEST = false;
bool APPROX = false;
int NUM_ITERATIONS = -1;

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

Node *find_best_sibling(Node *super_tree, vector<Node *> &gene_trees,
		int label);
void find_best_sibling_helper(Node *n, Node *new_leaf, Node *super_tree,
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
			APPROX=true;
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
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					MAX_CLUSTERS = atoi(arg2);
				cout << "MAX_CLUSTERS=" << MAX_CLUSTERS << endl;
			}
		}
		else if (strcmp(arg, "-cluster_test") == 0) {
			CLUSTER_TEST = true;
			PREFER_RHO = true;
		}
		else if (strcmp(arg, "-prefer_rho") == 0) {
			PREFER_RHO = true;
		}
		else if (strcmp(arg, "-memoize") == 0) {
			MEMOIZE = true;
		}
		else if (strcmp(arg, "-all_mafs") == 0) {
			ALL_MAFS= true;
		}
		else if (strcmp(arg, "-i") == 0) {
			if (max_args > argc) {
				char *arg2 = argv[argc+1];
				if (arg2[0] != '-')
					NUM_ITERATIONS = atoi(arg2);
				cout << "NUM_ITERATIONS=" << NUM_ITERATIONS << endl;
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
		CUT_AC_SEPARATE_COMPONENTS = true;
	}
	if (DEFAULT_ALGORITHM) {
		BB=true;
	}

	// initialize random number generator
	srand((unsigned(time(0))));

	// Label maps to allow string labels
	map<string, int> label_map= map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();
	vector<int> label_counts = vector<int>();

	string T_line = "";
	vector<Node *> gene_trees = vector<Node *>();
	vector<string> gene_tree_names = vector<string>();
	int skipped_multifurcating = 0;
	int skipped_small = 0;
	int skipped_no_bracket = 0;
	while (getline(cin, T_line)) {
		string name = "";
		size_t loc = T_line.find_first_of("(");
		if (loc != string::npos) {
			if (loc != 0) {
			name = T_line.substr(0,loc-1);
			T_line.erase(0,loc-1);
			}
			Node *T = build_tree(T_line);
			if (T->is_leaf() && T->str() == "p") {
				skipped_multifurcating++;
				continue;
			}
			if (T->size() <= 5) {
				cout << "skipped_small" << endl;
				cout << T->str_subtree() << endl;
				skipped_small++;
				continue;
			}
			gene_tree_names.push_back(name);
			gene_trees.push_back(T);
		}
		else
			skipped_no_bracket++;
	}

	cout << "skipped " << skipped_no_bracket << " lines with no opening bracket " << endl;
	cout << "skipped " << skipped_multifurcating << " multifurcating or invalid trees" << endl;
	cout << "skipped " << skipped_small << " trees with less than 4 leaves" << endl;
	cout << gene_trees.size() << " gene trees remaining" << endl;

	for(int i = 0; i < gene_trees.size(); i++) {
//		cout << gene_tree_names[i];
//		cout << gene_trees[i]->str_subtree() << endl;

		gene_trees[i]->labels_to_numbers(&label_map, &reverse_label_map);
		gene_trees[i]->count_numbered_labels(&label_counts);
	}

	// iterate over the taxa by number of occurences
	multimap<int, int> labels = multimap<int, int>();
	for(int i = 0; i < label_counts.size(); i++) {
		labels.insert(make_pair(label_counts[i],i));
	}


//	for(auto label = labels.rbegin(); label != labels.rend(); label++) {
//		cout << label->second ;
//		cout << " (" << reverse_label_map.find(label->second)->second << ") ";
//		cout << ": " << label->first << endl;
//	}

	// 4 most common leaves
	auto label = labels.rbegin();
	vector<int> l = vector<int>(4);
	for(int i = 0; i < 4; i++) {
		l[i] = label->second;
		label++;
		//cout << l[i] << endl;
	}

//	cout << endl;
//	cout << "Starting Trees:" << endl;

	// create all trees on the 4 leaves
	// TODO: really all or all 3 unrooted trees?
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
		best_distance = rSPR_total_approx_distance(ST[0], gene_trees); 
	else 
		best_distance = rSPR_total_distance(ST[0], gene_trees); 
	int best_tree = 0;
//	cout << best_distance <<  ": ";
//	cout << ST[0]->str_subtree() << endl;
	for(int j = 1; j < ST.size(); j++) {
		//cout << ST[j]->str_subtree() << endl;
		int distance;
		if (APPROX)
			distance = rSPR_total_approx_distance(ST[j], gene_trees);
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

	int x = 0;
	for(; label != labels.rend(); label++) {
		cout << "Adding leaf " << label->second << endl;
		Node *best_sibling = find_best_sibling(super_tree, gene_trees, label->second);
		Node *node = best_sibling->expand_parent_edge(best_sibling);

		node->add_child(new Node(itos(label->second)));
		cout << super_tree->str_subtree() << endl;
	}

	cout << endl;
	super_tree->numbers_to_labels(&reverse_label_map);
	cout << "Initial Supertree: " <<  super_tree->str_subtree() << endl;
	super_tree->labels_to_numbers(&label_map, &reverse_label_map);
	Node *best_supertree = new Node(*super_tree);
	if (APPROX)
		best_distance = rSPR_total_approx_distance(super_tree, gene_trees);
	else
		best_distance = rSPR_total_distance(super_tree, gene_trees);
	cout << "Total Distance: " << best_distance << endl;
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
			current_distance = rSPR_total_approx_distance(super_tree, gene_trees);
		else
			current_distance = rSPR_total_distance(super_tree, gene_trees);
		if (current_distance < best_distance) {
			best_supertree->delete_tree();
			best_supertree = new Node(*super_tree);
			best_distance = current_distance;
		}
		cout << "Total Distance: " << current_distance << endl;
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

	// TODO: modify this to keep a single intermediate node rather
	//than recreating and destroying it
	Node *new_node = n->expand_parent_edge(n);
	new_node->add_child(new_leaf);
//	cout << super_tree->str_subtree() << endl;

	super_tree->set_depth(0);
	super_tree->fix_depths();
	int distance;
	if (APPROX) {
		distance = rSPR_total_approx_distance(super_tree, gene_trees);
	}
	else {
		distance = rSPR_total_distance(super_tree, gene_trees);
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
			num_ties = 2;
		}
		else {
			num_ties++;
		}
	}
//	cout << "distance: " << distance;
//	cout << endl;

	new_leaf->cut_parent();
	//cout << super_tree->str_subtree() << endl;
	new_node = new_node->undo_expand_parent_edge();
	delete new_node;
	super_tree->set_depth(0);
	super_tree->fix_depths();
	//cout << super_tree->str_subtree() << endl;

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
		//cout << "foo2" << endl;
	// Note: we might want to allow this one!
	if (new_sibling == n->get_sibling()) {
//		cout << "rule 3" << endl;
		return;
	}
//		cout << "foo3" << endl;
//	cout << "foo4" << endl;
	if (n != super_tree && new_sibling != n) {
		Node *old_sibling = n->get_sibling();
		//if (new_sibling != old_sibling)
//		cout << "SPR Move:" << endl;
//		cout << "Previous Super Tree: "
//		<< super_tree->str_subtree() << endl;
//		cout << "Subtree: " << n->str_subtree() << endl;
//		cout << "New Sibling: " << new_sibling->str_subtree() << endl;

		int which_sibling = 0;
		Node *undo = n->spr(new_sibling, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
//		cout << "Super Tree: " << super_tree->str_subtree() << endl;

		int distance;
		if (APPROX) {
			distance = rSPR_total_approx_distance(super_tree, gene_trees);
		}
		else {
			distance = rSPR_total_distance(super_tree, gene_trees);
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
				num_ties = 2;
			}
			else {
				num_ties++;
			}
		}
		// restore the previous tree
		n->spr(undo, which_sibling);
		super_tree->set_depth(0);
		super_tree->fix_depths();
//		cout << "Reverted Super Tree: "
//	<< super_tree->str_subtree() << endl;
	}

}
