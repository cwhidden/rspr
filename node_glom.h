/*******************************************************************************
node_glom.h

Functions for Node Glom greedy supertre construction

Copyright 2013-2014 Chris Whidden
whidden@cs.dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 3, 2014
Version 1.2.1

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

*******************************************************************************/
#ifndef INCLUDE_NODE_GLOM

#define INCLUDE_NODE_GLOM
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
//#include <list>
//#include <deque>
#include "Node.h"
//#include "LCA.h"
//#include <map>
#include <limits>

using namespace std;

// function prototypes
void count_neighbours(Node *n, SparseCounts<double> *neighbour_counts,
		vector<int> *leaf_counts);
void count_neighbours_hlpr(Node *n, SparseCounts<double> *neighbour_counts,
		vector<int> *leaf_counts);
void count_leaves(Node *n, vector<int> *leaf_counts);
void glom_gene_tree_bottom_up(Node *n, int a, int b);
void glom_gene_tree_top_down(Node *n, Node *glom_root, int a, int b);
vector<vector<int> > find_component_trees(vector<Node *> *gene_trees, vector<Node *> *super_forest, int num_labels);
void append_component_trees(int tree, Node *n, vector<vector<int > > *component_trees);
int max(int a, int b);


// functions

void count_neighbours(Node *n, SparseCounts<double> *neighbour_counts,
		vector<int> *leaf_counts) {

	vector<Node *> leaf_children = vector<Node *>();
	vector<Node *> leaf_grandchildren = vector<Node *>();

	list<Node *>::const_iterator c1;
	for(c1 = n->get_children().begin(); c1 != n->get_children().end(); c1++) {
		if ((*c1)->is_leaf()) {
			leaf_children.push_back(*c1);
		}
		else {
			// recurse on children
			count_neighbours_hlpr(*c1, neighbour_counts, leaf_counts);

			if (n->get_children().size() == 2) {
				// find possible unrooted matches
				list<Node *>::const_iterator c2;
				for(c2 = (*c1)->get_children().begin(); c2 != (*c1)->get_children().end();
						c2++) {
					if ((*c2)->is_leaf()) {
						leaf_grandchildren.push_back(*c2);
					}
				}
			}
		}
	}

	// increment counts for each pair of leaf children
	for(int i = 0; i < leaf_children.size(); i++) {
		for(int j = i+1; j < leaf_children.size(); j++) {
			int num_i = leaf_children[i]->get_name_num();
			int num_j = leaf_children[j]->get_name_num();
			if (leaf_counts == NULL ||
					((*leaf_counts)[num_i] == 1 && (*leaf_counts)[num_j] == 1)) {
				neighbour_counts->increment(num_i, num_j);
				}
//			else {
//				neighbour_counts->increment(num_i, num_j, 1/max((*leaf_counts)[num_i],(*leaf_counts)[num_j]));
//			}
		}
	}

	// increment count for each pair of unrooted leaf children
	for(int i = 0; i < leaf_children.size(); i++) {
		for(int j = 0; j < leaf_grandchildren.size(); j++) {
			int num_i = leaf_children[i]->get_name_num();
			int num_j = leaf_grandchildren[j]->get_name_num();
			if (leaf_counts == NULL ||
					((*leaf_counts)[num_i] == 1 && (*leaf_counts)[num_j] == 1)) {
				neighbour_counts->increment(num_i, num_j);
				}
//			else {
//				neighbour_counts->increment(num_i, num_j, 1/max((*leaf_counts)[num_i],(*leaf_counts)[num_j]));
//			}
		}
	}
}

void count_neighbours_hlpr(Node *n, SparseCounts<double> *neighbour_counts,
		vector<int> *leaf_counts) {

	vector<Node *> leaf_children = vector<Node *>();

	list<Node *>::const_iterator c1;
	for(c1 = n->get_children().begin(); c1 != n->get_children().end(); c1++) {
		if ((*c1)->is_leaf()) {
			leaf_children.push_back(*c1);
		}
		else {
			// recurse on children
			count_neighbours_hlpr(*c1, neighbour_counts, leaf_counts);
		}
	}
	// increment counts for each pair of leaf children
	for(int i = 0; i < leaf_children.size(); i++) {
		for(int j = i+1; j < leaf_children.size(); j++) {
			int num_i = leaf_children[i]->get_name_num();
			int num_j = leaf_children[j]->get_name_num();
			if (leaf_counts == NULL ||
					((*leaf_counts)[num_i] == 1 && (*leaf_counts)[num_j] == 1)) {
				neighbour_counts->increment(num_i, num_j);
				}
//			else {
//				neighbour_counts->increment(num_i, num_j, 1/max((*leaf_counts)[num_i],(*leaf_counts)[num_j]));
//			}
		}
	}
}

void count_leaves(Node *n, vector<int> *leaf_counts) {
	if (n->is_leaf()) {
		(*leaf_counts)[n->get_name_num()]++;
	}
	else {
		list<Node *>::const_iterator c1;
		for(c1 = n->get_children().begin(); c1 != n->get_children().end();
				c1++) {
			count_leaves(*c1, leaf_counts);
		}
	}
}

// join the two chosen components of the super_forest
void glom_super_forest(vector<Node *> *super_forest, int a, int b) {
	Node *new_root = new Node();
	new_root->add_child((*super_forest)[a]);
	new_root->add_child((*super_forest)[b]);
	(*super_forest)[b] = NULL;
	(*super_forest)[a] = new_root;
}

// join the chosen nodes and relabel conflicts with the lower number
// 	move bottom up first and then find unrooted merges top down
void glom_gene_tree(Node *n, int a, int b) {

	glom_gene_tree_bottom_up(n, a, b);

	if (n->get_children().size() == 2) {
		list<Node *>::const_iterator c;
		Node *glom_root = NULL;
		for(c = n->get_children().begin(); c != n->get_children().end();
				c++) {
			if ((*c)->is_leaf()) {
				if (((*c)->get_name_num() == a || (*c)->get_name_num() == b)) {
					glom_root = *c;
					break;
				}
			}
		}

		if (glom_root != NULL) {
			c = n->get_children().begin();
			while(c != n->get_children().end()) {
				Node *child = *c;
				c++;
				if (!child->is_leaf()) {
					glom_gene_tree_top_down(child, glom_root, a, b);
				}
			}
		}
	}
}

void glom_gene_tree_bottom_up(Node *n, int a, int b) {

	list<Node *>::const_iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		if (!(*c)->is_leaf()) {
			glom_gene_tree_bottom_up((*c), a, b);
		}
	}

	vector<Node *> chosen_children = vector<Node *>();

	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		if ((*c)->is_leaf() &&
					((*c)->get_name_num() == a || (*c)->get_name_num() == b)) {
			chosen_children.push_back(*c);
		}
	}
	if (chosen_children.size() > 0) {
		// remove all but 1 of the chosen children
		for(int i = 1; i < chosen_children.size(); i++) {
			chosen_children[i]->cut_parent();
			chosen_children[i]->delete_tree();
		}
		// renumber the last child
		stringstream ss;
		ss << a;
		if (n->get_children().size() == 1) {
			chosen_children[0]->cut_parent();
			chosen_children[0]->delete_tree();
			n->set_name(ss.str());
		}
		else {
			chosen_children[0]->set_name(ss.str());
		}
	}
}

void glom_gene_tree_top_down(Node *n, Node *glom_root, int a, int b) {
	list<Node *>::const_iterator c;

	c = n->get_children().begin();
	bool cut_node = false;
	while(c != n->get_children().end()) {
		Node *child = *c;
		c++;
		if (child->is_leaf() && ((child->get_name_num() == a || child->get_name_num() == b))) {
			child->cut_parent();
			child->delete_tree();
			cut_node = true;
		}
	}
	Node *cont = NULL;
	if (n->get_children().size() == 1 && cut_node) {
//		cont = n->lchild();
		n->contract_node();
		cont = glom_root->get_sibling();
		if (!cont->is_leaf()) {
			glom_gene_tree_top_down(cont, glom_root, a, b);
		}
	}
}

vector<vector<int> > find_component_trees(vector<Node *> *gene_trees, vector<Node *> *super_forest, int num_labels) {
	vector<vector<int> > component_trees = vector<vector<int> >(num_labels);
	for(int i = 0; i < num_labels; i++) {
		
//		if ((*super_forest)[i] != NULL) {
			component_trees[i] = vector<int>();
//		}
	}
	for(int i = 0; i < gene_trees->size(); i++) {
		append_component_trees(i, (*gene_trees)[i],&component_trees);
	}
	return component_trees;
}

void append_component_trees(int tree, Node *n, vector<vector<int > > *component_trees) {
	if (n->is_leaf()) {
		int label = atoi(n->get_name().c_str());
		(*component_trees)[label].push_back(tree);
	}
	else {
		list<Node *>::const_iterator c;
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			append_component_trees(tree, *c, component_trees);
		}
	}
}

int max(int a, int b) {
	if (a > b)
		return a;
	return b;
}

#endif
