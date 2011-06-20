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
#ifndef INCLUDE_FOREST

#define INCLUDE_FOREST
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <list>
#include <deque>
#include "Node.h"
#include "LCA.h"
#include <map>
#include <limits>

using namespace std;

class Forest {
	public:
		vector<Node *> components;
		vector<Node *> deleted_nodes;
		bool rho;

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
		rho = false;
	}
	Forest(const Forest &f) {
		components = vector<Node *>(f.components.size());
		for(int i = 0; i < f.components.size(); i++) {
			//if (f.components[i] != NULL)
			components[i] = new Node(*f.components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f.rho;
	}

	Forest(Forest *f) {
		components = vector<Node *>(f->components.size());
		for(int i = 0; i < f->components.size(); i++) {
			//if (f->components[i] != NULL)
			components[i] = new Node(*f->components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f->rho;
	}

	void init(vector<Node *> components) {
		this->components = vector<Node *>(components);
		deleted_nodes = vector<Node *>();
		rho = false;
	}
	~Forest() {
		for(int i = 0; i < components.size(); i++) {
			//if (components[i] != NULL) {
				components[i]->delete_tree();
				//components[i] = NULL;
			//}
		}
		for(int i = 0; i < deleted_nodes.size(); i++) {
			//if (deleted_nodes[i] != NULL) {
				delete deleted_nodes[i];
				//deleted_nodes[i] = NULL;
			//}
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

		bool rho_temp = this->rho;
		this->rho = f->rho;
		f->rho = rho_temp;
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

	// print the components seperated by s
	void print_components(string s) {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else
				root->print_subtree_hlpr();
			cout << s;
		}
		cout << endl;
	}

	// print the forest showing twins
	void print_components_with_twins() {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else
				root->print_subtree_twin_hlpr();
			cout << " ";
		}
		cout << endl;
	}

	inline void add_component(Node *head) {
		components.push_back(head);
	}
	inline void add_component(int pos, Node *head) {
		components.insert(components.begin() + pos, head);
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
	// clear twin pointers
	void unsync() {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			(*i)->unsync();
		}
	}
	// clear interior twin pointers
	void unsync_interior() {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			(*i)->unsync_interior();
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

int size() {
	return components.size();
}

bool add_rho() {
	if (rho)
		return false;
	Node *T_p = new Node("p");
	add_component(T_p);
	rho = true;
	return true;
}

bool contains_rho() {
	return rho;
}

void erase_components(int start, int end) {
	components.erase(components.begin()+start, components.begin()+end);
}
void erase_components() {
	components.clear();
}
};

// Functions

vector<Node *> find_labels(vector<Node *> components);
void sync_twins(Forest *T1, Forest *T2);
void sync_interior_twins(Forest *T1, Forest *T2);
void sync_interior_twins(Node *n, LCA *twin_LCA);
void sync_interior_twins(Node *n, vector<LCA> *F2_LCAs);
list<Node *> *find_cluster_points(Forest *f);
void find_cluster_points(Node *n, list<Node *> *cluster_points);
void delete_and_merge_LCAs(list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs, list<Node *>:: iterator node1_location,
		list<Node *>:: iterator node2_location);
void delete_and_merge_LCAs(Node *n, list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs);


// return the smallest number in s
int stomini(string s) {
//	cout << "stomini" << endl;
	string number_characters = "+-0123456789";
	int i = 0;
	int min = INT_MAX;
	string current = "";
	for(int i = 0; i < s.size(); i++) {
		if (number_characters.find(s[i]) != string::npos) {
			current += s[i];
		}
		else if (current.size() > 0) {
			int num = atoi(current.c_str());
			if (num < min)
				min = num;
			current = "";
		}
	}
	if (current.size() > 0) {
		int num = atoi(current.c_str());
		if (num < min)
			min = num;
		current = "";
	}
//	cout << "returning " << min << endl;
	return min;
}


// Make the leaves of two forests point to their twin in the other tree
// Note: removes unique leaves
void sync_twins(Forest *T1, Forest *T2) {
	vector<Node *> T1_labels = vector<Node *>();
	vector<Node *> T2_labels = vector<Node *>();
	vector<Node *> T1_components = T1->components;
	vector<Node *> T2_components = T2->components;
	vector<Node *>::iterator i;
	Node *T1_rho = NULL;
	Node *T2_rho = NULL;
	for(i = T1_components.begin(); i != T1_components.end(); i++) {
		Node *component = *i;
		vector<Node *> unsorted_labels = component->find_leaves();
		vector<Node *>::iterator j;
		for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
			Node *leaf = *j;
			cout << "T1: " << leaf->str() << endl;
			if (leaf->str() == "p") {
				T1_rho = leaf;
			}
			else {
				// find smallest number contained in the label
				int number = stomini(leaf->str());
				cout << "\t" << number << endl;
				if (number >= T1_labels.size())
					T1_labels.resize(number+1, NULL);
				T1_labels[number] = leaf;
			}
		}
	}
	for(i = T2_components.begin(); i != T2_components.end(); i++) {
		Node *component = *i;
		vector<Node *> unsorted_labels = component->find_leaves();
		vector<Node *>::iterator j;
		for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
			Node *leaf = *j;
			cout << "T2: " << leaf->str() << endl;
			if (leaf->str() == "p") {
				T2_rho = leaf;
			}
			else {
				// find smallest number contained in the label
				int number = stomini(leaf->str());
				cout << "\t" << number << endl;
				if (number >= T2_labels.size())
					T2_labels.resize(number+1, NULL);
				T2_labels[number] = leaf;
			}
		}
	}
	T1_labels.resize(T1_labels.size()+1);
	T1_labels[T1_labels.size()-1]=T1_rho;
	T2_labels.resize(T2_labels.size()+1);
	T2_labels[T2_labels.size()-1]=T2_rho;

	int size = T1_labels.size();
	if (size > T2_labels.size())
		size = T2_labels.size();
	cout << "Syncing Twins" << endl;
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
			cout << T1_a->str() << endl;
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
//	cout << "Finished Syncing Twins" << endl;
	return;
}

/* make interior nodes point to the lca of their descendants in the other tree
   assumes that sync_twins has already been called
   assumes that component 1 of T1 matches with 1 of T2
      NOTE: this isn't true during the algorithm so this will need to be changed
      if we want to interleave clustering. It should be just component 1 of T1
	  matching multiple components of T2 (The first several components?)
   */
void sync_interior_twins(Forest *T1, Forest *T2) {
	Node  *root1 = T1->get_component(0);
	Node  *root2 = T2->get_component(0);
	LCA T1_LCA = LCA(root1);
	LCA T2_LCA = LCA(root2);
	sync_interior_twins(root1, &T2_LCA);
	sync_interior_twins(root2, &T1_LCA);
}

/* make interior nodes point to the LCA of their descendants in the other
	 forest if there is one unambiguous LCA
	 * This is true for a node n of T1 if all leaves that are a descendant
	 	of T1 either map to a single component of F2 or are from another
	 	component of F2 such that the root of that component maps to
	 	a descendant of n (i.e. a finished component)
   * assumes that sync_twins has already been called
	 */
// TODO: initializing parameters seems to be slow
void sync_interior_twins_real(Forest *T1, Forest *F2) {
	Node  *T1_root = T1->get_component(0);
	LCA T1_LCA = LCA(T1_root);
	int T1_size = T1_root->size_using_prenum();
	// roots of F2
	vector<Node *> F2_roots = vector<Node *>();
	// LCA queries for F2
	vector<LCA> F2_LCAs = vector<LCA>();
	// lists of root nodes that map to a given T1 node
	T1_root->initialize_parameter(ROOT_LCAS, list<Node *>());
	// list of active descendants
	T1_root->initialize_parameter(ACTIVE_DESCENDANTS, list<Node *>());

	// should be fine.
	for(int i = 0; i < F2->num_components(); i++) {
//		cout << "starting i" << endl;
		F2_roots.push_back(F2->get_component(i));
#ifdef DEBUG
		cout << "COMPONENT " << i << ": " << F2_roots[i]->str_subtree() << endl;
		cout << "foo" << endl;
		cout << F2_roots[i]->get_twin() << endl;
		if (F2_roots[i]->get_twin() != NULL) {
		cout << F2_roots[i]->get_twin()->str_subtree() << endl;
		cout << F2_roots[i]->get_twin()->parent() << endl;
		cout << "fooa" << endl;
		if (F2_roots[i]->get_twin()->parent() != NULL) {
		cout << F2_roots[i]->get_twin()->parent()->str_subtree() << endl;
		}
		cout << "foob" << endl;
		}
#endif
		// ignore finished components
//		if (F2_roots[i]->get_twin() != NULL && F2_roots[i]->get_twin()->parent() == NULL) {
//			cout << "fooc" << endl;
		//F2_LCAs.push_back(F2_roots[i]);
//			F2_LCAs.push_back(LCA(F2_roots[i]));
//			cout << "fooc" << endl;
//			continue;
//		}
		cout << "foo" << endl;
		// ignore rho components
		if (F2_roots[i]->str() == "p") {
		//	F2_LCAs.push_back(LCA());
		F2_LCAs.push_back(F2_roots[i]);
		//F2_LCAs.push_back(NULL);//LCA(F2_roots[i]));
			continue;
		}
		cout << "aa" << endl;
		F2_LCAs.push_back(LCA(F2_roots[i]));
		// number the component
		F2_roots[i]->initialize_parameter(COMPONENT_NUMBER,i);
		// list of nodes that get deleted when a component is finished
		F2_roots[i]->initialize_parameter(REMOVABLE_DESCENDANTS,list<list<Node *>::iterator>());
		// sync the component with T1
		if (F2_roots[i]->str() != "p" &&
				!(F2_roots[i]->get_twin() != NULL && F2_roots[i]->get_twin()->parent() == NULL)) {
			sync_interior_twins(F2_roots[i], &T1_LCA);
		}
		// keep reverse pointer for the root's twin
		cout << "a" << endl;
		/*
		cout << F2_roots[i] << endl;
		cout << F2_roots[i]->str_subtree() << endl;
		cout << F2_roots[i]->get_twin() << endl;
		cout << F2_roots[i]->get_twin()->str_subtree() << endl;
		cout << F2_roots[i]->get_twin()->get_parameter_ref(ACTIVE_DESCENDANTS) << endl;
		cout << F2_roots[i]->get_twin()->get_parameter_ref(ROOT_LCAS) << endl;
		cout << boost::any_cast<list<Node *> >(F2_roots[i]->get_twin()->get_parameter_ref(ROOT_LCAS))->size() << endl;
		*/
		if (i > 0)
			boost::any_cast<list<Node *> >(F2_roots[i]->get_twin()->get_parameter_ref(ROOT_LCAS))->push_back(F2_roots[i]);
		else
			boost::any_cast<list<Node *> >(T1_root->get_parameter_ref(ROOT_LCAS))->push_back(F2_roots[i]);
//		cout << "b" << endl;
	}
	cout << "syncing" << endl;
	sync_interior_twins(T1_root, &F2_LCAs); 
}

/* make interior nodes point to the lca of their descendants in the other
 * tree
 * assumes that sync_twins has already been called
 */
void sync_interior_twins(Node *n, LCA *twin_LCA) {
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	if (lc != NULL)
		sync_interior_twins(lc, twin_LCA);
	if (rc != NULL)
		sync_interior_twins(rc, twin_LCA);
	if (lc == NULL && rc != NULL)
		n->set_twin(rc->get_twin());
	else if (lc != NULL && rc == NULL)
		n->set_twin(lc->get_twin());
	else if (lc != NULL && rc != NULL) {
		Node *twin = twin_LCA->get_lca(lc->get_twin(), rc->get_twin());
		n->set_twin(twin);
	}
}

void sync_interior_twins(Node *n, vector<LCA> *F2_LCAs) {
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	list<Node *> *active_descendants = boost::any_cast<list<Node *> >(n->get_parameter_ref(ACTIVE_DESCENDANTS));
	// visit left subtree first
	if (lc != NULL)
		sync_interior_twins(lc, F2_LCAs);
	// visit right subtree first
	if (rc != NULL)
		sync_interior_twins(rc, F2_LCAs);
	cout << "SYNC_INTERIOR_TWINS()" << endl;
	cout << n->str_subtree() << endl;
	if (lc == NULL && rc == NULL) {
		cout << "leaf" << endl;
		active_descendants->push_back(n->get_twin());
		list<Node *>::iterator node_location = active_descendants->end();
		node_location--;
			boost::any_cast<list<list<Node *>::iterator> >(n->get_twin()->get_parameter_ref(REMOVABLE_DESCENDANTS))->push_back(node_location);
	}
	// no lc so propogate up
	if (lc == NULL && rc != NULL) {
		cout << "no lc" << endl;
		n->set_twin(rc->get_twin());
		list<Node *> *rc_active_descendants = boost::any_cast<list<Node *> >(rc->get_parameter_ref(ACTIVE_DESCENDANTS));
		active_descendants->splice(active_descendants->end(),*rc_active_descendants);
	}
	// no rc so propogate up
	if (lc != NULL && rc == NULL) {
		cout << "no rc" << endl;
		n->set_twin(lc->get_twin());
		list<Node *> *lc_active_descendants = boost::any_cast<list<Node *> >(lc->get_parameter_ref(ACTIVE_DESCENDANTS));
		active_descendants->splice(active_descendants->end(),*lc_active_descendants);
	}
	// two children so put their info together
	else if (lc != NULL && rc != NULL) {
		cout << "two children" << endl;
		list<Node *> *lc_active_descendants = boost::any_cast<list<Node *> >(lc->get_parameter_ref(ACTIVE_DESCENDANTS));
		list<Node *> *rc_active_descendants = boost::any_cast<list<Node *> >(rc->get_parameter_ref(ACTIVE_DESCENDANTS));

	cout << "active_descendants lc" << endl;
	for(list<Node *>::iterator i =  lc_active_descendants-> begin(); i != lc_active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
		cout << endl;
	cout << "active_descendants rc" << endl;
	for(list<Node *>::iterator i =  rc_active_descendants-> begin(); i != rc_active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
		cout << endl;

		bool done = false;
		if (lc_active_descendants->empty() || rc_active_descendants->empty())
			done = true;
		cout << active_descendants->size() << endl;
		active_descendants->splice(active_descendants->end(),*lc_active_descendants);
		cout << active_descendants->size() << endl;
		list<Node *>::iterator node1_location = active_descendants->end();
		node1_location--;
		list<Node *>::iterator node2_location = node1_location;
		active_descendants->splice(active_descendants->end(),*rc_active_descendants);
		node2_location++;
		/* check the intersection point to see if we have two
			leaves from the same component
		*/
//		cout << "foo" << endl;
		cout << active_descendants->size() << endl;
		if (!done)
			delete_and_merge_LCAs(active_descendants, F2_LCAs, node1_location,
					node2_location);
		cout << "done first merge" << endl;

		/* check to see if n is twinned by a root of F2
			 if so, then remove each leaf twinned by that component
			 and check each of the new intersection points
		*/
		list<Node *> *root_lcas = boost::any_cast<list<Node *> >(n->get_parameter_ref(ROOT_LCAS));
		while(!root_lcas->empty()) {

			Node *root_lca = root_lcas->front();
			root_lcas->pop_front();
			/* TODO: problem when n is a root
				 We don't care about this but it might mean there is a different
				 problem
				 */
			if (n->parent() != NULL) {
//			cout << "deleting from component " << 
//	boost::any_cast<int>(root_lca->get_parameter(COMPONENT_NUMBER))
//			<< endl;
				delete_and_merge_LCAs(root_lca, active_descendants, F2_LCAs);
			}
		}
//		cout << "done checking component" << endl;

		/* If we have a single element in n's active descendants
			 list then set twin pointers appropriately
		*/

//	cout << "active_descendants done" << endl;
//	for(list<Node *>::iterator i =  active_descendants->begin(); i != active_descendants->end(); i++) {
//		cout << "\t" << (*i)->str_subtree() << endl;
//	}
		if (active_descendants->size() == 1) {
//			cout << "found twin" << endl;
			Node *twin = active_descendants->front();
			n->set_twin(twin);
		}
	}
}

/* merge two nodes from a list into their LCA if they are from
	 the same component
	 */
void delete_and_merge_LCAs(list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs, list<Node *>:: iterator node1_location,
		list<Node *>:: iterator node2_location) {
#ifdef DEBUG
	cout << "DELETE_AND_MERGE FIRST" << endl;
	cout << "active_descendants before" << endl;
	for(list<Node *>::iterator i =  active_descendants-> begin(); i != active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
	cout << endl;
#endif

//	while(active_descendants->size() > 1) {
	Node *n1 = *node1_location;
	Node *n2 = *node2_location;
#ifdef DEBUG
	cout << "n1=" << n1->str_subtree() << endl;
	cout << "n2=" << n2->str_subtree() << endl;
#endif
	int component1 = boost::any_cast<int>(n1->get_parameter(COMPONENT_NUMBER));
	int component2 = boost::any_cast<int>(n2->get_parameter(COMPONENT_NUMBER));
#ifdef DEBUG
	cout << "c1=" << component1 << endl;
	cout << "c2=" << component2 << endl;
#endif
	if (component1 == component2) {
#ifdef DEBUG
		cout << "size=" << (*F2_LCAs).size() << endl;
		for(int i = 0; i < (*F2_LCAs).size(); i++) {
			if ((*F2_LCAs)[i].get_tree() == NULL)
				cout << "\t" << "NULL" << endl;
			else
				cout << "\t" << (*F2_LCAs)[i].get_tree()->str_subtree() << endl;
		}
		cout << (*F2_LCAs)[component1].get_tree()->str_subtree() << endl;
#endif
		Node *lca = (*F2_LCAs)[component1].get_lca(n1,n2);
#ifdef DEBUG
		cout << lca->str_subtree() << endl;
#endif
//		cout << "xa" << endl;
		list<Node *>::iterator lca_location =
			active_descendants->insert(node1_location,lca);
//		cout << "xb" << endl;
		active_descendants->erase(node1_location);
//		cout << "xc" << endl;
		boost::any_cast<list<list<Node *>::iterator> >(n1->get_parameter_ref(REMOVABLE_DESCENDANTS))->clear();
//		cout << "xd" << endl;
		active_descendants->erase(node2_location);
//		cout << "xe" << endl;
		boost::any_cast<list<list<Node *>::iterator> >(n2->get_parameter_ref(REMOVABLE_DESCENDANTS))->clear();
//		cout << "xf" << endl;
		boost::any_cast<list<list<Node *>::iterator> >(lca->get_parameter_ref(REMOVABLE_DESCENDANTS))->push_back(lca_location);
//		cout << "xg" << endl;
	}
//		else {
//			break;
//		}
//	}
#ifdef DEBUG
	cout << "active_descendants after" << endl;
	for(list<Node *>::iterator i =  active_descendants-> begin(); i != active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
#endif
}

/* delete each leaf from the list that is twinned with the component
	 of n. For each such deleted node, merge its predecessor
	 and successor in the list into their LCA if they are from
	 the same component (other than n's component)
	 */
void delete_and_merge_LCAs(Node *n, list<Node *>
		*active_descendants, vector<LCA> *F2_LCAs) {
//	cout << n->str_subtree() << endl;
	int component = boost::any_cast<int>(n->get_parameter(COMPONENT_NUMBER));
	list<list<Node *>::iterator> *removable_descendants	=
			boost::any_cast<list<list<Node *>::iterator> >(n->get_parameter_ref(REMOVABLE_DESCENDANTS));

//	cout << "removable_descendants" << endl;
//	for(list<list<Node *>::iterator>::iterator i =  removable_descendants-> begin(); i != removable_descendants->end(); i++) {
//		cout << "\t" << (**i)->str_subtree() << endl;
//	}
//	cout << "foo" << endl;
	while (!removable_descendants->empty()) {
//		cout << "fooa" << endl;
		list<Node *>::iterator leaf_location = removable_descendants->front();
//		cout << "foob" << endl;
		removable_descendants->pop_front();
//		cout << "fooc" << endl;
//		cout << *leaf_location << endl;
//		cout << (*leaf_location)->str_subtree() << endl;
//		cout << "food" << endl;

		if (active_descendants->front() != *leaf_location
				&& active_descendants->back() != *leaf_location) {
			list<Node *>::iterator node1_location = leaf_location;
//		cout << "fooe" << endl;
			list<Node *>::iterator node2_location = leaf_location;
//		cout << "foof" << endl;
			node1_location--;
//		cout << "foog" << endl;
			node2_location++;
//		cout << "fooh" << endl;
			active_descendants->erase(leaf_location);
//		cout << "fooi" << endl;
			int node1_component = boost::any_cast<int>((*node1_location)->get_parameter(COMPONENT_NUMBER));
//		cout << "fooj" << endl;
			if (component != node1_component)
				delete_and_merge_LCAs(active_descendants, F2_LCAs, node1_location,
						node2_location);
//		cout << "fook" << endl;
		}
		else if (active_descendants->size() > 1){
			active_descendants->erase(leaf_location);
		}

	}
//	cout << "foo end" << endl;
	// TODO: continue to lc and rc?
	if (n->lchild() != NULL)
		delete_and_merge_LCAs(n->lchild(), active_descendants, F2_LCAs);
	if (n->rchild() != NULL)
		delete_and_merge_LCAs(n->rchild(), active_descendants, F2_LCAs);
}

list<Node *> *find_cluster_points(Forest *F) {
	list<Node *> *cluster_points = new list<Node *>();
	find_cluster_points(F->get_component(0), cluster_points);
	//cout << "foo" << endl;
	return cluster_points;
}

// find the cluster points
void find_cluster_points(Node *n, list<Node *> *cluster_points) {
	//cout << "Start: " << n->str_subtree() << endl;
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	if (lc != NULL)
		find_cluster_points(lc, cluster_points);
	if (rc != NULL)
		find_cluster_points(rc, cluster_points);
	/*
	cout << "here" << endl;
	cout << n->get_depth() << endl;
	if (n->get_twin() != NULL)
		cout << n->get_twin()->get_twin()->get_depth() << endl;
	cout << n->parent() << endl;
	if (lc != NULL) {
	cout << "lc= " << lc->str_subtree() << endl;
	cout << lc->get_depth() << endl;
	if (lc->get_twin() != NULL)
	cout << lc->get_twin()->get_twin()->get_depth() << endl;
	}
	if (rc != NULL) {
	cout << "rc= " << rc->str_subtree() << endl;
	cout << rc->get_depth() << endl;
	if (rc->get_twin() != NULL)
	cout << rc->get_twin()->get_twin()->get_depth() << endl;
	}
	*/
	if (n->get_twin() != NULL
			&& n->parent() != NULL
			&& lc != NULL
			&& rc != NULL
			&& n->get_depth() <= n->get_twin()->get_twin()->get_depth()
			&& (lc->get_twin() == NULL
				|| lc->get_depth() > lc->get_twin()->get_twin()->get_depth()
				|| rc->get_twin() == NULL
				|| rc->get_depth() > rc->get_twin()->get_twin()->get_depth())) {
//		cout << "added cluster_point" << endl;
		cluster_points->push_back(n);
	}
	//cout << "End: " << n->str_subtree() << endl;
}

// swap two forests
void swap(Forest **a, Forest **b) {
	(*a)->swap(*b);
}

// expand all contracted nodes
void expand_contracted_nodes(Forest *F) {
	for(int i = 0; i < F->size(); i++) {
		expand_contracted_nodes(F->get_component(i));
	}
}
#endif
