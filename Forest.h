/*******************************************************************************
Forest.h

Data structure for a forest of binary trees

Copyright 2009-2014 Chris Whidden
cwhidden@dal.ca
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
//#include "ClusterInstance.h"

using namespace std;

bool MULTI_CLUSTER = false;

class ClusterInstance;
extern bool LEAF_REDUCTION2;

class Forest {
	public:
		vector<Node *> components;
		vector<Node *> deleted_nodes;
		bool rho;
		Forest *twin;
		ClusterInstance *cluster;

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
		twin = NULL;
		cluster = NULL;
	}
	Forest(const Forest &f) {
		components = vector<Node *>(f.components.size());
		for(int i = 0; i < f.components.size(); i++) {
			//if (f.components[i] != NULL)
			components[i] = new Node(*f.components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f.rho;
		twin = NULL;
		cluster = f.cluster;
		//label_nodes_with_forest();
	}

	Forest(Forest *f) {
		components = vector<Node *>(f->components.size());
		for(int i = 0; i < f->components.size(); i++) {
			//if (f->components[i] != NULL)
			components[i] = new Node(*f->components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f->rho;
		twin = NULL;
		cluster = f->cluster;
		//label_nodes_with_forest();
	}

	Forest(Forest *f, bool b) {
		components = vector<Node *>(f->components.size());
		for(int i = 0; i < f->components.size(); i++) {
			//if (f->components[i] != NULL)
			components[i] = new Node(*f->components[i]);
		}
		deleted_nodes = vector<Node *>();
		rho = f->rho;
		twin = NULL;
		cluster = f->cluster;
		//label_nodes_with_forest();
	}


	void init(vector<Node *> components) {
		this->components = vector<Node *>(components);
		deleted_nodes = vector<Node *>();
		rho = false;
		for(int i = 0; i < components.size(); i++) {
				if (components[i]->str() == "p")
					rho = true;
		}
		twin = NULL;
		cluster = NULL;
	}
	~Forest() {
		for(int i = 0; i < components.size(); i++) {
			//if (components[i] != NULL) {
				components[i]->delete_tree();
				components[i] = NULL;
			//}
		}
		for(int i = 0; i < deleted_nodes.size(); i++) {
			//if (deleted_nodes[i] != NULL) {
				deleted_nodes[i]->delete_tree();
				deleted_nodes[i] = NULL;
			//}
		}
	} 

	// swap the contents of two forests
	void swap(Forest *f) {
		vector<Node *> components_temp = this->components;
		this->components = f->components;
		f->components = components_temp;
		
//		/*
		vector<Node *> deleted_nodes_temp = this->deleted_nodes;
		this->deleted_nodes = f->deleted_nodes;
		f->deleted_nodes = deleted_nodes_temp;
//		*/

		bool rho_temp = this->rho;
		this->rho = f->rho;
		f->rho = rho_temp;

		ClusterInstance *c_temp = this->cluster;
		this->cluster = f->cluster;
		f->cluster = c_temp;
	}

	// print the forest
	void print_components() {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else if (root->is_leaf() && root->str() == "")
				cout << "*";
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

	// print the forest showing edge intervals
	void print_components_with_edge_pre_interval() {
		vector<Node *>::iterator it = components.begin();
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				cout << "!";
			else
				cout << root->str_edge_pre_interval_subtree();
			cout << " ";
		}
		cout << endl;
	}

	// return the string for this forest
	string str() {
		string s = "";
		vector<Node *>::iterator it;
		for(it = components.begin(); it != components.end(); it++) {
			Node *root = *it;
			if (root == NULL)
				s += "!";
			else
				s += root->str_subtree();
			s += " ";
		}
		return s;
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
	void set_component(int i, Node *head) {
		components[i] = head;
	}
	inline int num_components() {
		return components.size();
	}
	void set_twin(Forest *f) {
		twin = f;
	}
	Forest *get_twin() {
		return twin;
	}
	ClusterInstance *get_cluster() {
		return cluster;
	}
	void set_cluster(ClusterInstance *c) {
		cluster = c;
	}
	void set_cluster(ClusterInstance &c) {
		cluster = &c;
	}
	void update_component(Node *old_c, Node *new_c) {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *component = *i;
			if (&(*component) == &(*old_c))
				*i = new_c;
		}
	}

	// return a list of the sibling pairs
	list<Node *> *find_sibling_pairs() {
		list<Node *> *sibling_pairs = new list<Node *>();
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *component = *i;
			component->append_sibling_pairs(sibling_pairs);
		}
		return sibling_pairs;
	}

	// return a deque of the singleton leaves
	list<Node *> find_singletons() {
		list<Node *> singletons = list<Node *>();
		vector<Node *>::iterator i = components.begin();
		// TODO: is this a hack or correct? We don't want the first
		// component to be a singleton because of rho!
		i++;
		for(; i != components.end(); i++) {
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
	// clear twin pointers
	Node *find_by_prenum(int prenum) {
		vector<Node *>::iterator i;
		for(i = components.begin(); i != components.end(); i++) {
			Node *ans = (*i)->find_by_prenum(prenum);
			if (ans != NULL)
				return ans;
		}
		return NULL;
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

void set_rho(bool b) {
	rho = b;
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

// tell the nodes what forest they are in
void label_nodes_with_forest() {
	for(int i = 0; i < size(); i++) {
		get_component(i)->set_forest_rec(this);
	}
}

void move_first_component_to_end() {
	Node *temp = components[0];
	components[0] = components[components.size()-1];
	components[components.size()-1] = temp;
}
void unprotect_edges() {
	for(int i = 0; i < size(); i++) {
		get_component(i)->unprotect_subtree();
	}
}

};

// Functions

vector<Node *> find_labels(vector<Node *> components);
bool sync_twins(Forest *T1, Forest *T2);
void sync_interior_twins(Forest *T1, Forest *T2);
void sync_interior_twins(Node *n, LCA *twin_LCA);
void sync_interior_twins(Node *n, vector<LCA> *F2_LCAs);
list<Node *> *find_cluster_points(Forest *F1, Forest *F2);
void find_cluster_points(Node *n, list<Node *> *cluster_points,
		vector<int> *leaf_counts_F1, vector<int> *leaf_counts_F2);
void delete_and_merge_LCAs(list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs, list<Node *>:: iterator node1_location,
		list<Node *>:: iterator node2_location);
void delete_and_merge_LCAs(Node *n, list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs);



// Make the leaves of two forests point to their twin in the other tree
// Note: removes unique leaves
bool sync_twins(Forest *T1, Forest *T2) {
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
//			cout << "T1: " << leaf->str() << endl;
			if (leaf->str() == "p") {
				T1_rho = leaf;
			}
			else {
				// find smallest number contained in the label
				int number = stomini(leaf->str());
//				cout << "\t" << number << endl;
				if (number < INT_MAX) {
					if (number >= T1_labels.size())
						T1_labels.resize(number+1, 0);
					T1_labels[number] = leaf;
				}
			}
		}
	}
	for(i = T2_components.begin(); i != T2_components.end(); i++) {
		Node *component = *i;
		vector<Node *> unsorted_labels = component->find_leaves();
		vector<Node *>::iterator j;
		for(j = unsorted_labels.begin(); j != unsorted_labels.end(); j++) {
			Node *leaf = *j;
//			cout << "T2: " << leaf->str() << endl;
			if (leaf->str() == "p") {
				T2_rho = leaf;
			}
			else {
				// find smallest number contained in the label
				int number = stomini(leaf->str());
//				cout << "\t" << number << endl;
				if (number < INT_MAX) {
					if (number >= T2_labels.size())
						T2_labels.resize(number+1, 0);
					T2_labels[number] = leaf;
				}
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
//	cout << "Syncing Twins" << endl;
	for(int i = 0; i < size; i++) {
		Node *T1_a = T1_labels[i];
		Node *T2_a = T2_labels[i];
		if (T1_a == NULL && T2_a != NULL) {
			Node *node = T2_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T2_a)
						sibling = node->rchild();
				T2_labels[stomini(sibling->str())] = sibling;
			}
			delete T2_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
		}
		else if (T2_a == NULL && T1_a != NULL) {
			Node *node = T1_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T1_a)
						sibling = node->rchild();
				T1_labels[stomini(sibling->str())] = sibling;
			}
			delete T1_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
			
		}
		if (T1_a != NULL && T2_a != NULL) {
			T1_a->set_twin(T2_a);
			T2_a->set_twin(T1_a);
//			cout << T1_a->str() << endl;
		}
	}
	for(int i = size; i < T1_labels.size(); i++) {
		Node *T1_a = T1_labels[i];
		if (T1_a != NULL) {
			Node *node = T1_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T1_a)
						sibling = node->rchild();
				T1_labels[stomini(sibling->str())] = sibling;
			}
			delete T1_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
			
		}
	}
	for(int i = size; i < T2_labels.size(); i++) {
		Node *T2_a = T2_labels[i];
		if (T2_a != NULL) {
			Node *node = T2_a->parent();
			if (node == NULL)
				return false;
			int numc = node->get_children().size();
			if (node->parent() == NULL && node->lchild()->is_leaf() &&
					(numc == 1 || (numc == 2 && node->rchild()->is_leaf()))) {
				return false;
				Node *sibling = node->lchild();
				if (sibling == T2_a)
						sibling = node->rchild();
				T2_labels[stomini(sibling->str())] = sibling;
			}
			delete T2_a;
			if (node->get_children().size() < 2) {
				if (node->get_children().size() == 1)
					node->lchild()->lost_child();
				node = node->contract(true);
			}
		}
	}
//	if (T1_loss != NULL)
//		*T1_loss = T1->get_component(0)->count_lost_children_subtree();
//		- T1->get_component(0)->num_lost_children();;
//	cout << "T1_lost = " << T1_lost << endl;
//	if (T2_loss != NULL)
//		*T2_loss = T2->get_component(0)->count_lost_children_subtree();
//		- T2->get_component(0)->num_lost_children();;
//	cout << "T2_lost = " << T2_lost << endl;

//	cout << "Finished Syncing Twins" << endl;
	return true;
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
	T1_root->initialize_root_lcas(list<Node *>());
	// list of active descendants
	T1_root->initialize_active_descendants(list<Node *>());

	// should be fine.
	for(int i = 0; i < F2->num_components(); i++) {
//		cout << "starting i" << endl;
		F2_roots.push_back(F2->get_component(i));
#ifdef DEBUG_SYNC
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
//		cout << "foo" << endl;
		// ignore rho components
		if (F2_roots[i]->str() == "p") {
		//	F2_LCAs.push_back(LCA());
		F2_LCAs.push_back(F2_roots[i]);
		//F2_LCAs.push_back(NULL);//LCA(F2_roots[i]));
			continue;
		}
//		cout << "aa" << endl;
		F2_LCAs.push_back(LCA(F2_roots[i]));
		// number the component
		F2_roots[i]->initialize_component_number(i);
		// list of nodes that get deleted when a component is finished
		F2_roots[i]->initialize_removable_descendants(list<list<Node *>::iterator>());
		// sync the component with T1
		if (F2_roots[i]->str() != "p" &&
				!(F2_roots[i]->get_twin() != NULL && F2_roots[i]->get_twin()->parent() == NULL)) {
			sync_interior_twins(F2_roots[i], &T1_LCA);
		}
		// keep reverse pointer for the root's twin
//		cout << "a" << endl;
		/*
		cout << F2_roots[i] << endl;
		cout << F2_roots[i]->str_subtree() << endl;
		cout << F2_roots[i]->get_twin() << endl;
		cout << F2_roots[i]->get_twin()->str_subtree() << endl;
		cout << F2_roots[i]->get_twin()->get_parameter_ref(ACTIVE_DESCENDANTS) << endl;
		cout << F2_roots[i]->get_twin()->get_parameter_ref(ROOT_LCAS) << endl;
		cout << boost::any_cast<list<Node *> >(F2_roots[i]->get_twin()->get_parameter_ref(ROOT_LCAS))->size() << endl;
		*/
		if (i > 0 || T1->contains_rho())
			F2_roots[i]->get_twin()->get_root_lcas()->push_back(F2_roots[i]);
		else
			T1_root->get_root_lcas()->push_back(F2_roots[i]);
//		cout << "b" << endl;
	}
//	cout << "syncing" << endl;
	sync_interior_twins(T1_root, &F2_LCAs); 
}

/* make interior nodes point to the lca of their descendants in the other
 * tree
 * assumes that sync_twins has already been called
 */
void sync_interior_twins(Node *n, LCA *twin_LCA) {
	list<Node *>::iterator c = n->get_children().begin();
	if (c == n->get_children().end())
		return;
	sync_interior_twins(*c, twin_LCA);
	n->set_twin((*c)->get_twin());
	c++;
	while(c != n->get_children().end()) {
		sync_interior_twins(*c, twin_LCA);
		Node *twin = twin_LCA->get_lca(n->get_twin(), (*c)->get_twin());
		n->set_twin(twin);
		c++;
	}
}

void sync_interior_twins(Node *n, vector<LCA> *F2_LCAs) {
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	list<Node *> *active_descendants = n->get_active_descendants();
	// visit children first
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		sync_interior_twins(*c, F2_LCAs);
	}
	#ifdef DEBUG_SYNC
	cout << "SYNC_INTERIOR_TWINS()" << endl;
	cout << n->str_subtree() << endl;
	#endif
	if (n->get_children().size() == 0) {
//		cout << "leaf" << endl;
		active_descendants->push_back(n->get_twin());
		list<Node *>::iterator node_location = active_descendants->end();
		node_location--;
			n->get_twin()->get_removable_descendants()->push_back(node_location);
	}
	// no rc so propogate up
	if (n->get_children().size() == 1) {
		Node *lc = n->get_children().front();
//		cout << "no rc" << endl;
		n->set_twin(lc->get_twin());
		list<Node *> *lc_active_descendants = lc->get_active_descendants();
		active_descendants->splice(active_descendants->end(),*lc_active_descendants);
	}
	// TODO: generalize from here for 2 or more children
	// two children so put their info together
	else if (lc != NULL && rc != NULL) {
//		cout << "two children" << endl;
		list<Node *> *lc_active_descendants = lc->get_active_descendants();
		list<Node *> *rc_active_descendants = rc->get_active_descendants();

/*	#ifdef DEBUG_SYNC
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
	#endif
*/

		vector<list<Node *>::iterator> node_location =
				vector<list<Node *>::iterator>();
		list<Node *>::iterator node1_location;
		int nonempty_active_descendants_count = 0;
		for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
			if (!(*c)->get_active_descendants()->empty()) {
				nonempty_active_descendants_count++;
				if (nonempty_active_descendants_count > 1) {
					node1_location = active_descendants->end();
					node1_location--;
					node_location.push_back(node1_location);
				}
//		cout << active_descendants->size() << endl;
				active_descendants->splice(active_descendants->end(),
						*((*c)->get_active_descendants()));
//		cout << active_descendants->size() << endl;
			}
		}

		/* check the intersection points to see if we have two
			leaves from the same component
		*/
//		cout << "foo" << endl;
		#ifdef DEBUG_SYNC
		cout << active_descendants->size() << endl;
		#endif
		for(int i = 0; i < node_location.size(); i++) {
			list<Node *>::iterator node1_location = node_location[i];
			list<Node *>::iterator node2_location = node1_location;
			node2_location++;
			delete_and_merge_LCAs(active_descendants, F2_LCAs, node1_location,
					node2_location);
		}
		#ifdef DEBUG_SYNC
		cout << "done first merge" << endl;
		#endif

		/* check to see if n is twinned by a root of F2
			 if so, then remove each leaf twinned by that component
			 and check each of the new intersection points
		*/
		list<Node *> *root_lcas = n->get_root_lcas();
		while(!root_lcas->empty()) {

			Node *root_lca = root_lcas->front();
			root_lcas->pop_front();
			/* TODO: problem when n is a root
				 We don't care about this but it might mean there is a different
				 problem
				 */
			if (n->parent() != NULL) {
			#ifdef DEBUG_SYNC
			cout << "deleting from component " << endl;
			#endif
				delete_and_merge_LCAs(root_lca, active_descendants, F2_LCAs);
			}
		}
		#ifdef DEBUG_SYNC
		cout << "done checking component" << endl;
		#endif

		/* If we have a single element in n's active descendants
			 list then set twin pointers appropriately
		*/

	#ifdef DEBUG_SYNC
	cout << "active_descendants done" << endl;
	for(list<Node *>::iterator i =  active_descendants->begin(); i != active_descendants->end(); i++) {
		cout << "\t" << (*i)->str_subtree() << endl;
	}
	#endif
		if (active_descendants->size() == 1) {
			#ifdef DEBUG_SYNC
			cout << "found twin" << endl;
			#endif
			Node *twin = active_descendants->front();
			n->set_twin(twin);
		}
	}
	if (n->parent() == NULL)
		active_descendants->clear();
}

void sync_af_twins(Forest *F1, Forest *F2) {
	F1->unsync();
	F2->unsync();
	sync_twins(F1, F2);
	for(int i = 0; i < F1->num_components(); i++) {
		F1->get_component(i)->sync_af_twins();
	}
}

/* merge two nodes from a list into their LCA if they are from
	 the same component
	 */
void delete_and_merge_LCAs(list<Node *> *active_descendants,
		vector<LCA> *F2_LCAs, list<Node *>:: iterator node1_location,
		list<Node *>:: iterator node2_location) {
#ifdef DEBUG_SYNC
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
#ifdef DEBUG_SYNC
	cout << "n1=" << n1->str_subtree() << endl;
	cout << "n2=" << n2->str_subtree() << endl;
#endif
	int component1 = n1->get_component_number();
	int component2 = n2->get_component_number();
#ifdef DEBUG_SYNC
	cout << "c1=" << component1 << endl;
	cout << "c2=" << component2 << endl;
#endif
	if (component1 == component2) {
#ifdef DEBUG_SYNC
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
#ifdef DEBUG_SYNC
		cout << lca->str_subtree() << endl;
#endif
//		cout << "xa" << endl;
		list<Node *>::iterator lca_location =
			active_descendants->insert(node1_location,lca);
//		cout << "xb" << endl;
//		active_descendants->erase(node1_location);
//		cout << "xc" << endl;

// TODO: could this be faster?
		bool remove = false;
		list<list<Node *>::iterator>::iterator i;
		for(i = n1->get_removable_descendants()->begin(); i != n1->get_removable_descendants()->end(); i++) {
			if (*i == node1_location) {
				//active_descendants->erase(*i);
				remove = true;
				break;
			}
		}
		if (remove) {
			active_descendants->erase(*i);
			n1->get_removable_descendants()->erase(i);
		}
//		n1->get_removable_descendants()->clear();
		// TODO: delete each when clearing?
//		cout << "xd" << endl;
//		active_descendants->erase(node2_location);
//		cout << "xe" << endl;
		remove = false;
		for(i = n2->get_removable_descendants()->begin(); i != n2->get_removable_descendants()->end(); i++) {
			if (*i == node2_location) {
				//active_descendants->erase(*i);
				remove = true;
				break;
			}
		}
		if (remove) {
			active_descendants->erase(*i);
			n2->get_removable_descendants()->erase(i);
		}
//		n2->get_removable_descendants()->clear();
//		cout << "xf" << endl;
		lca->get_removable_descendants()->push_back(lca_location);
//		cout << "xg" << endl;
	}
//		else {
//			break;
//		}
//	}
#ifdef DEBUG_SYNC
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
	int component = n->get_component_number();
	list<list<Node *>::iterator> *removable_descendants	=
			n->get_removable_descendants();

	if (n->lchild() != NULL)
		delete_and_merge_LCAs(n->lchild(), active_descendants, F2_LCAs);
	if (n->rchild() != NULL)
		delete_and_merge_LCAs(n->rchild(), active_descendants, F2_LCAs);
	#ifdef DEBUG_SYNC

	cout << n->str_subtree() << endl;
	cout << "removable_descendants" << endl;
	for(list<list<Node *>::iterator>::iterator i =  removable_descendants-> begin(); i != removable_descendants->end(); i++) {
		cout << "\t" << (**i)->str_subtree() << endl;
	}	
	#endif
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

// TODO: problem here
		if (leaf_location != active_descendants->begin() &&
				leaf_location != active_descendants->end() &&
				leaf_location != -- active_descendants->end()) {
//		if (active_descendants->front() != *leaf_location
//				&& active_descendants->back() != *leaf_location) {
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
			int node1_component = (*node1_location)->get_component_number();
//		cout << "fooj" << endl;
			if (component != node1_component)
				delete_and_merge_LCAs(active_descendants, F2_LCAs, node1_location,
						node2_location);
//		cout << "fook" << endl;
		}
		else {//if (active_descendants->size() > 1){
			active_descendants->erase(leaf_location);
		}

	}
//	cout << "foo end" << endl;
	// TODO: continue to lc and rc?
}

list<Node *> *find_cluster_points(Forest *F1, Forest *F2) {
	list<Node *> *cluster_points = new list<Node *>();
	vector<int> *leaf_counts_F1 = NULL;
	vector<int> *leaf_counts_F2 = NULL;
	if (MULTI_CLUSTER) {
		leaf_counts_F1 = F1->get_component(0)->find_leaf_counts(); 
		leaf_counts_F2 = F2->get_component(0)->find_leaf_counts(); 
	}
	find_cluster_points(F1->get_component(0), cluster_points, leaf_counts_F1,
			leaf_counts_F2);
	if (MULTI_CLUSTER) {
		delete leaf_counts_F1;
//		delete leaf_counts_F2;
	}
	//cout << "foo" << endl;
	return cluster_points;
}

// find the cluster points
void find_cluster_points(Node *n, list<Node *> *cluster_points,
		vector<int> *leaf_counts_F1, vector<int> *leaf_counts_F2) {
//	cout << "Start: " << n->str_subtree() << endl;
	list<Node *>::iterator c;
	for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
		find_cluster_points(*c, cluster_points, leaf_counts_F1,
				leaf_counts_F2);
	}
	/*
	cout << "here" << endl;
	cout << "n= " << n->str_subtree() << endl;
	cout << n->get_depth() << endl;
	if (n->get_twin() != NULL) {
		cout << "n_twin= " << n->get_twin()->str_subtree() << endl;
		cout << "n_twin_twin= " << n->get_twin()->get_twin()->str_subtree() << endl;
		cout << n->get_twin()->get_twin()->get_depth() << endl;
	}
	cout << n->parent() << endl;
	*/
	bool is_cluster = true;
	int num_clustered_children = 0;
	if (n->get_twin() == NULL ||
			n->parent() == NULL ||
			n->get_children().size() < 2 ||
			n->get_depth() > n->get_twin()->get_twin()->get_depth())
		is_cluster = false;
	else {
#ifdef RSPR
		if (!LEAF_REDUCTION2){
#endif
			for(c = n->get_children().begin(); c != n->get_children().end(); c++) {
				if ((*c)->get_twin() != NULL &&
						(*c)->get_depth() <= (*c)->get_twin()->get_twin()->get_depth())
				num_clustered_children++;
			}
			if (num_clustered_children == n->get_children().size())
				is_cluster = false;
#ifdef RSPR
		}
#endif
	}
	if (is_cluster) {
//		cout << "added cluster_point" << endl;
		cluster_points->push_back(n);
	}
	// buggy, needs testing, doesn't seem worth it
	else if (MULTI_CLUSTER && n->get_twin() != NULL && n->parent() != NULL &&
			n->get_children().size() >= 2) {
		// TODO: use find_leaf_counts if this works
		Node *n_twin = n->get_twin();
		int num_leaves = (*leaf_counts_F1)[n->get_preorder_number()];
		vector<Node *> chosen = vector<Node *>();
		int chosen_leaves = 0;
		if (n_twin != NULL && n->get_edge_pre_start() > -1 && n->get_edge_pre_end() > -1 && n_twin->get_children().size() > 2) {
//			cout << "foo" << endl;
//			cout << n->str_subtree() << endl;
//			cout << num_leaves << endl;
//			cout << n->get_edge_pre_start() << endl;
//			cout << n->get_edge_pre_end() << endl;
			for(c = n_twin->get_children().begin(); c != n_twin->get_children().end(); c++) {
//				cout << "\t" << (*c)->str_subtree() << endl;
				int c_num_leaves = (*leaf_counts_F2)[(*c)->get_preorder_number()];
//				cout << "\t" << c_num_leaves << endl;
				int c_twin_pre = (*c)->get_twin()->get_preorder_number();
//				cout << "\t" << c_twin_pre << endl;
				if (c_twin_pre >= n->get_edge_pre_start() &&
						c_twin_pre <= n->get_edge_pre_end()) {
//					cout << "yes" << endl;
					chosen.push_back(*c);
					chosen_leaves += c_num_leaves;
				}
			}
			// PROBLEM: the new node should have its own preorder number
			// and its own size
			if (num_leaves == chosen_leaves) {
				Node *new_child = new Node();
				n_twin->add_child(new_child);
				new_child->set_preorder_number(n_twin->get_preorder_number());
				new_child->set_edge_pre_start(n_twin->get_edge_pre_start());
				new_child->set_edge_pre_end(n_twin->get_edge_pre_end());
				n->set_twin(new_child);
				new_child->set_twin(n);
				cluster_points->push_back(n);
				for(int i = 0; i < chosen.size(); i++) {
					new_child->add_child(chosen[i]);
				}
			}

		}
	}
//	cout << "End: " << n->str_subtree() << endl;
}

// swap two forests
void swap(Forest **a, Forest **b) {
	(*a)->swap(*b);
}

// expand all contracted nodes
void expand_contracted_nodes(Forest *F) {
	for(int i = 0; i < F->num_components(); i++) {
		F->get_component(i)->expand_contracted_nodes();
	}
}

Forest *build_finished_forest(string &name) {
	Forest *new_forest = new Forest();
	string::iterator i = name.begin();
		size_t old_loc = 0;
		size_t loc = 0;
		while ((loc = name.find(" ", old_loc)) != string::npos) {
			//cout << "old_loc=" << old_loc << endl;
			//cout << "loc=" << loc << endl;
			//cout << name.substr(old_loc,loc-old_loc) << endl;
			//new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
			new_forest->add_component(new Node(name.substr(old_loc,loc-old_loc)));
			if (name.substr(old_loc,loc-old_loc) == "p")
				new_forest->set_rho(true);
			//new_forest->print_components();
			old_loc = loc+1;
		}
		new_forest->add_component(new Node(name.substr(old_loc,loc-old_loc)));
		return new_forest;
		//new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
}
Forest *build_forest(string &name) {
	Forest *new_forest = new Forest();
	string::iterator i = name.begin();
		size_t old_loc = 0;
		size_t loc = 0;
		while ((loc = name.find(" ", old_loc)) != string::npos) {
			//cout << "old_loc=" << old_loc << endl;
			//cout << "loc=" << loc << endl;
			//cout << name.substr(old_loc,loc-old_loc) << endl;
			//new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
			new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
			if (name.substr(old_loc,loc-old_loc) == "p")
				new_forest->set_rho(true);
			//new_forest->print_components();
			old_loc = loc+1;
		}
		new_forest->add_component(build_tree(name.substr(old_loc,loc-old_loc)));
		return new_forest;
}
#endif
