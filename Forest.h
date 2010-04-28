/*******************************************************************************
Forest.h

Data structure for a forest of binary trees

Copyright 2009-2010 Chris Whidden
cwhidden@dal.ca
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

*******************************************************************************/
#ifndef INCLUDE_CSTDIO
	#define INCLUDE_CSTDIO
	#include <cstdio>
#endif
#ifndef INCLUDE_CSTDLIB
	#define INCLUDE_CSTDLIB
	#include <cstdlib>
#endif
#ifndef INCLUDE_STRING
	#define INCLUDE_STRING
	#include <string>
#endif
#ifndef INCLUDE_CSTRING
	#define INCLUDE_CSTRING
	#include <cstring>
#endif
#ifndef INCLUDE_IOSTREAM
	#define INCLUDE_IOSTREAM
	#include <iostream>
#endif
#ifndef INCLUDE_SSTREAM
	#define INCLUDE_SSTREAM
	#include <sstream>
#endif
#ifndef INCLUDE_MATH
	#define INCLUDE_MATH
	#include <math.h>
#endif
#ifndef INCLUDE_VECTOR
	#define INCLUDE_VECTOR
	#include <vector>
#endif
#ifndef INCLUDE_LIST
	#define INCLUDE_LIST
	#include <list>
#endif
#ifndef INCLUDE_DEQUE
	#define INCLUDE_DEQUE
	#include <deque>
#endif
#ifndef INCLUDE_NODE
	#define INCLUDE_NODE
	#include "Node.h"
#endif
#ifndef INCLUDE_MAP
	#define INCLUDE_MAP
	#include <map>
#endif
using namespace std;

vector<Node *> find_labels(vector<Node *> components);

class Forest {
	public:
		vector<Node *> components;
		vector<Node *> deleted_nodes;

	public:
	Forest() {
		init(vector<Node *>());
	}
	Forest(vector<Node *> components) {
		init(components);
	}
	Forest(Node *head) {
		components = vector<Node *>();
		components.push_back(new Node(*head));
		deleted_nodes = vector<Node *>();
	}
	Forest(const Forest &f) {
		components = vector<Node *>(f.components.size());
		for(int i = 0; i < f.components.size(); i++) {
			//if (f.components[i] != NULL)
			components[i] = new Node(*f.components[i]);
		}
		deleted_nodes = vector<Node *>();
	}

	void init(vector<Node *> components) {
		this->components = components;
		deleted_nodes = vector<Node *>();
	}
	~Forest() {
		for(int i = 0; i < components.size(); i++) {
			//if (components[i] != NULL)
				components[i]->delete_tree();
		}
		for(int i = 0; i < deleted_nodes.size(); i++) {
			//if (deleted_nodes[i] != NULL)
				delete deleted_nodes[i];
		}
	} 

	// swap the contents of two forests
	void swap(Forest *f) {
		vector<Node *> components_temp = this->components;
		this->components = f->components;
		f->components = components_temp;
		
		/*
		vector<Node *> deleted_nodes_temp = this->deleted_nodes;
		this->deleted_nodes = f->deleted_nodes;
		f->deleted_nodes = deleted_nodes_temp;
		
		*/
	}

	// print the forest
	void print_components() {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else
				root->print_subtree_hlpr();
			cout << " ";
		}
		cout << endl;
	}

	inline void add_component(Node *head) {
		components.push_back(head);
	}
	inline void add_deleted_node(Node *n) {
		deleted_nodes.push_back(n);
	}

	inline Node *get_component(int i) {
		return components[i];
	}
	inline int num_components() {
		return components.size();
	}

	// return a deque of the sibling pairs
	deque<Node *> find_sibling_pairs() {
		deque<Node *> sibling_pairs = deque<Node *>();
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *component = *i;
			vector<Node *> new_sibling_pairs = component->find_sibling_pairs();
			sibling_pairs.insert(sibling_pairs.end(), new_sibling_pairs.begin(), new_sibling_pairs.end());
		}
		return sibling_pairs;
	}

	// return a deque of the singleton leaves
	list<Node *> find_singletons() {
		list<Node *> singletons = list<Node *>();
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *component = *i;
			if (component->is_leaf()) {
				singletons.push_back(component);
			}
		}
		return singletons;
	}

	// make nodes pointed to in the forest point back
	void resync() {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			(*i)->resync();
		}
	}

void labels_to_numbers(map<string, int> *label_map, map<int, string> *reverse_label_map) {
	vector<Node *>::iterator i;
	for(i = components.begin(); i != components.end(); i++) {
		(*i)->labels_to_numbers(label_map, reverse_label_map);
	}
}

void numbers_to_labels(map<int, string> *reverse_label_map) {
	vector<Node *>::iterator i;
	for(i = components.begin(); i != components.end(); i++) {
		(*i)->numbers_to_labels(reverse_label_map);
	}
}

};

// Make the leaves of two forests point to their twin in the other tree
void sync_twins(Forest *T1, Forest *T2) {
	vector<Node *> T1_labels = vector<Node *>();
	vector<Node *> T2_labels = vector<Node *>();
	vector<Node *> T1_components = T1->components;
	vector<Node *> T2_components = T2->components;
	vector<Node *>::iterator i;
	for(i = T1_components.begin(); i != T1_components.end(); i++) {
		Node *component = *i;
		vector<Node *> unsorted_labels = component->find_leaves();
		vector<Node *>::iterator j;
		for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
			Node *leaf = *j;
			string name = leaf->str();
			int number = atoi(name.c_str());
			if (number >= T1_labels.size())
				T1_labels.resize(number+1, NULL);
			T1_labels[number] = leaf;
		}
	}
	for(i = T2_components.begin(); i != T2_components.end(); i++) {
		Node *component = *i;
		vector<Node *> unsorted_labels = component->find_leaves();
		vector<Node *>::iterator j;
		for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
			Node *leaf = *j;
			string name = leaf->str();
			int number = atoi(name.c_str());
			if (number >= T2_labels.size())
				T2_labels.resize(number+1, NULL);
			T2_labels[number] = leaf;
		}
	}
	int size = T1_labels.size();
	if (size > T2_labels.size())
		size = T2_labels.size();
	for(int i = 0; i < size; i++) {
		Node *T1_a = T1_labels[i];
		Node *T2_a = T2_labels[i];
		if (T1_a == NULL && T2_a != NULL) {
			Node *node = T2_a->parent();
			delete T2_a;
			if (node != NULL)
				node->contract();
		}
		else if (T2_a == NULL && T1_a != NULL) {
			Node *node = T1_a->parent();
			delete T1_a;
			if (node != NULL)
				node->contract();
		}
		if (T1_a != NULL && T2_a != NULL) {
			T1_a->set_twin(T2_a);
			T2_a->set_twin(T1_a);
		}
	}
	for(int i = size; i < T1_labels.size(); i++) {
		Node *T1_a = T1_labels[i];
		Node *node = T1_a->parent();
		delete T1_a;
		if (node != NULL)
			node->contract();
	}
	for(int i = size; i < T2_labels.size(); i++) {
		Node *T2_a = T2_labels[i];
		Node *node = T2_a->parent();
		delete T2_a;
		if (node != NULL)
			node->contract();
	}
	return;
}
