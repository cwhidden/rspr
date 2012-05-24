/*******************************************************************************
Node.h

Data structure for a node of a binary tree
Contains methods to recursively work on a node's subtree

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

#ifndef INCLUDE_NODE

#define INCLUDE_NODE

#define COPY_CONTRACTED

#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include "Forest.h"

using namespace std;

// representation of a component with no leaves
#define DEAD_COMPONENT "*"
/*void find_sibling_pairs_hlpr(Node *node, list<Node *> &sibling_pairs);
void find_leaves_hlpr(Node *node, vector<Node *> &leaves);
void str_subtree_hlpr(string *s);
list<Node *> find_sibling_pairs(Node *node);
vector<Node *> find_leaves(Node *node);
*/
int stomini(string s);

class Forest;

class Node {
	private:
	Node *lc;			// left child
	Node *rc;			// right child
	Node *p;			// parent
	Node *twin;			// counterpart in another tree
	string name;		// label
	int depth;			//distance from root
	int pre_num;	// preorder number

	int component_number;
	list <Node *> active_descendants;
	list <Node *> root_lcas;
	list<list <Node *>::iterator> removable_descendants;
	list <Node *>::iterator sibling_pair_loc;
	int sibling_pair_status;
	int num_clustered_children;
	Forest *forest;
	Node *contracted_lc;
	Node *contracted_rc;

	public:
	Node() {
		init(NULL, NULL, NULL, "", 0);
	}
	Node(string n) {
		init(NULL, NULL, NULL, n, 0);
	}
	Node(string n, int d) {
		init(NULL, NULL, NULL, n, d);
	}
	Node(Node *lc, Node *rc, Node *p, string n, int d) {
		init(lc, rc, p, n, d);
	}
	void init(Node *lc, Node *rc, Node *p, string n, int d) {
		this->lc = lc;
		this->rc = rc;
		this->p = p;
		this->name = string(n);
		this->twin = NULL;
		this->depth = d;
		this->pre_num = -1;
		this->component_number = -2;
		this->active_descendants = list <Node *>();
		this->root_lcas = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		this->contracted_lc = NULL;
		this->contracted_rc = NULL;
	}
	// copy constructor
	Node(const Node &n) {
		p = NULL;
		name = n.name;
		twin = n.twin;
//		depth = n.depth;
		depth = 0;
		pre_num = n.pre_num;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		if (n.lc == NULL)
			lc = NULL;
		else
			lc = new Node(*(n.lc), this);
		if (n.rc == NULL)
			rc = NULL;
		else
			rc = new Node(*(n.rc), this);
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else
			contracted_lc = new Node(*(n.contracted_lc), this);
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else
			contracted_rc = new Node(*(n.contracted_rc), this);
#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
#endif
	}
	Node(const Node &n, Node *parent) {
		p = parent;
		name = n.name;
		twin = n.twin;
		if (p != NULL)
			depth = p->depth+1;
		else
			depth = 0;
		pre_num = n.pre_num;
		component_number = n.component_number;
		this->active_descendants = list <Node *>();
		this->removable_descendants = list< list<Node *>::iterator>();
		this->root_lcas = list <Node *>();
		//sibling_pair_loc = n.sibling_pair_loc;
		//sibling_pair_status = n.sibling_pair_status;
		this->sibling_pair_loc = list<Node *>::iterator(); 
		this->sibling_pair_status = 0;
		this->num_clustered_children = 0;
		this->forest = NULL;
		if (n.lc == NULL)
			lc = NULL;
		else
			lc = new Node(*(n.lc), this);
		if (n.rc == NULL)
			rc = NULL;
		else
			rc = new Node(*(n.rc), this);
#ifdef COPY_CONTRACTED
		if (n.contracted_lc == NULL)
			contracted_lc = NULL;
		else
			contracted_lc = new Node(*(n.contracted_lc), this);
		if (n.contracted_rc == NULL)
			contracted_rc = NULL;
		else
			contracted_rc = new Node(*(n.contracted_rc), this);
#else
		this->contracted_lc = n.contracted_lc;
		this->contracted_rc = n.contracted_rc;
#endif
	}
	~Node() {
		if (lc != NULL) {
			lc->p = NULL;
			//lc = NULL;
		}
		if (rc != NULL) {
			rc->p = NULL;
			//rc = NULL;
		}
		if (p != NULL) {
			p->delete_child(this);
			//p = NULL;
		}
		active_descendants.clear();
		root_lcas.clear();
		removable_descendants.clear();
	}
	void fake_delete() {
		if (lc != NULL) {
			lc->p = NULL;
			//lc = NULL;
		}
		if (rc != NULL) {
			rc->p = NULL;
			//rc = NULL;
		}
		if (p != NULL) {
			p->delete_child(this);
			//p = NULL;
		}
		active_descendants.clear();
		root_lcas.clear();
		removable_descendants.clear();
	}
	// delete a subtree
	void delete_tree() {
		Node *lc = this->lc;
		Node *rc = this->rc;
		if (lc != NULL) {
			lc->delete_tree();
		}
		lc = NULL;
		if (rc != NULL) {
			rc->delete_tree();
		}
		rc = NULL;
#ifdef COPY_CONTRACTED
		if (contracted_lc != NULL) {
			contracted_lc->delete_tree();
		}
		contracted_lc = NULL;
		if (contracted_rc != NULL) {
			contracted_rc->delete_tree();
		}
		contracted_rc = NULL;
#endif
		delete this;
	}

	// cut edge between parent and child
	// should really be cut_child, no deleting occurs
	void delete_child(Node *n) {

		if (lc != NULL && lc == n) {
			lc = NULL;
		}
		else if (rc != NULL && rc == n) {
			rc = NULL;
		}
		else {
		}
	}

	// add a child if possible
	int add_child(Node *n) {
		if (lc == NULL) {
			set_lchild(n);
			return 0;
		}
		else if (rc == NULL) {
			set_rchild(n);
			return 0;
		}
		return -1;
	}


	Node *set_lchild(Node *n) {
		this->lc = n;
		if (n != NULL) {
			n->p = this;
			n->depth = depth+1;
		}
		return lc;
	}

	Node *set_lchild_keep_depth(Node *n) {
		this->lc = n;
		if (n != NULL) {
			n->p = this;
		}
		return lc;
	}
	Node *set_rchild(Node *n) {
		this->rc = n;
		if (n != NULL) {
			n->p = this;
			n->depth = depth+1;
		}
		return rc;
	}
	Node *set_rchild_keep_depth(Node *n) {
		this->rc = n;
		if (n != NULL) {
			n->p = this;
		}
		return rc;
	}
	Node *set_parent(Node *n) {
		p = n;
		return p;
	}
	Node *set_twin(Node *n) {
		twin = n;
		return twin;
	}
	Node *set_name(string n) {
		name = string(n);
	}
	int set_depth(int d) {
		depth = d;
		return depth;
	}
	void fix_depths() {
		if (lc != NULL) {
			lc->depth = depth+1;
			lc->fix_depths();
		}
		if (rc != NULL) {
			rc->depth = depth+1;
			rc->fix_depths();
		}
	}
	int set_preorder_number(int p) {
		pre_num = p;
		return pre_num;
	}

	int set_component_number(int c) {
		component_number = c;
	}
	int get_component_number() {
		return component_number;
	}
	void increase_clustered_children() {
		num_clustered_children++;
	}
	void decrease_clustered_children() {
		num_clustered_children--;
	}
	int set_num_clustered_children(int c) {
		num_clustered_children = c;
	}
	int get_num_clustered_children() {
		return num_clustered_children;
	}
	list <Node *> *get_active_descendants() {
		return &active_descendants;
	}
	list <Node *> *get_root_lcas(){
		return &root_lcas;
	}
	int get_sibling_pair_status(){
		return sibling_pair_status;
	}
	int set_sibling_pair_status(int s){
		sibling_pair_status = s;
	}
	void set_forest(Forest *f) {
		forest = f;
	}
	void set_forest_rec(Forest *f) {
		forest = f;
		if (lc != NULL)
			lc->set_forest_rec(f);
		if (rc != NULL)
			rc->set_forest_rec(f);
	}
	Forest *get_forest() {
		return forest;
	}
	list<list <Node *>::iterator> *get_removable_descendants() {
		return &removable_descendants;
	}
	void initialize_component_number(int value) {
		component_number = value;
		if (lc != NULL)
			lc->initialize_component_number(value);
		if (rc != NULL)
			rc->initialize_component_number(value);
	}
	void initialize_active_descendants(list <Node *> value) {
		active_descendants = value;
		if (lc != NULL)
			lc->initialize_active_descendants(value);
		if (rc != NULL)
			rc->initialize_active_descendants(value);
	}
	void initialize_root_lcas(list <Node *> value) {
		root_lcas = value;
		if (lc != NULL)
			lc->initialize_root_lcas(value);
		if (rc != NULL)
			rc->initialize_root_lcas(value);
	}
	void initialize_removable_descendants(list<list <Node *>::iterator> value) {
		removable_descendants = value;
		if (lc != NULL)
			lc->initialize_removable_descendants(value);
		if (rc != NULL)
			rc->initialize_removable_descendants(value);
	}

	/* contract:
	 * if this node has a parent and one child then contract it out
	 * if this node has a parent and no child then contract it out
	 *   parent will have one child so contract it as well.
	 * if this node has no parent and one child then take the
	 *	child's children and contract it out
	 * return the first degree two parent found or NULL if there
	 * was no contraction
	 */
	// TODO: check handling twins for interior nodes
	Node *contract(bool remove) {
		Node *parent = p;
		Node *child;
		Node *ret = NULL;
		// contract out this node and give child to parent
		if (parent != NULL) {
		//cout << p->str_subtree() << endl;
			if (lc && !rc) {
				child = lc;
				if (remove)
					delete this;
				//this->fake_delete();
				parent->delete_child(this);
				parent->add_child(child);
				ret = parent;
			}
			else if (rc && !lc) {
				child = rc;
				if (remove)
					delete this;
				//this->fake_delete();
				parent->delete_child(this);
				parent->add_child(child);
				ret = parent;
			}
			else if (lc == NULL && rc == NULL) {
				parent->delete_child(this);
				if (remove)
					delete this;
				//this->fake_delete();
				ret = parent->contract(remove);
			}
			else
				ret = this;

		}
		// if no parent then take children of single child and remove it
		else {

			// dead component or singleton, will be cleaned up by the forest
			if (lc == NULL && rc == NULL) {
				if (str() == "")
					name = DEAD_COMPONENT;
			}
			else if ((bool)lc xor (bool)rc) {
				child = lc;
				if (rc == NULL) {
					child = lc;
					lc->p = NULL;
					lc = NULL;
				}
				else {
					child = rc;
					rc->p = NULL;
					rc = NULL;
				}
				/* cluster hack - if we delete a cluster node then
				 * we may try to use it later. This only happens once
				 * per cluster so we can spend linear time to update
				 * the forest
				 */
				if (child->num_clustered_children > 0) {
					delete_child(child);
					if (remove)
						delete this;
					//this->fake_delete();
					ret = child;
				}
				else {
					// if child is a leaf then get rid of this so we don't lose refs
					// problem: if the child is not c, then we want to copy
					// otherwise we don't
					// copy other parameters and join the twin
					//to this if the child is a label
					Node *new_lc = child->lchild();
					Node *new_rc = child->rchild();
					if (child->is_leaf()) {
						if (child->get_twin() != NULL) {
							set_twin(child->get_twin());
							child->get_twin()->set_twin(this);
						}
						name = child->str();
					}
					delete_child(child);
					if (remove)
						delete child;
					//child->fake_delete();
					if (new_lc != NULL)
						add_child(new_lc);
					if (new_rc != NULL)
						add_child(new_rc);
					ret = this;
				}
			}
		}

		return ret;
	}

	Node *contract() {
		return contract(false);
	}

	/* contract_sibling_pair:
	 * if this node has two child leaves then contract them out
	 * return true if contracted, otherwise false
	 */
	bool contract_sibling_pair() {
		if (lc != NULL && lc->is_leaf()
				&& rc != NULL && rc->is_leaf()) {
			#ifdef DEBUG
				string new_name = "<" + lc->str() + "," + rc->str() + ">";
			#else
				string new_name = "(" + lc->str() + "," + rc->str() + ")";
			#endif
			set_name(new_name);
			lc->set_parent(NULL);
			rc->set_parent(NULL);
			lc = NULL;
			rc = NULL;
			return true;
		}
		return false;
	}

	bool contract_sibling_pair_undoable() {
		if (lc != NULL && lc->is_leaf()
				&& rc != NULL && rc->is_leaf()) {
			/*
			#ifdef DEBUG
				string new_name = "<" + lc->str() + "," + rc->str() + ">";
			#else
				string new_name = "(" + lc->str() + "," + rc->str() + ")";
			#endif
			set_name(new_name);
			*/
			lc->set_parent(NULL);
			rc->set_parent(NULL);
			contracted_lc = lc;
			contracted_rc = rc;
			lc = NULL;
			rc = NULL;
			return true;
		}
		return false;
	}

	void undo_contract_sibling_pair() {
		lc = contracted_lc;
		contracted_lc = NULL;
		rc = contracted_rc;
		contracted_rc = NULL;
		lc->set_parent(this);
		rc->set_parent(this);
	}


	// cut the edge between this node and its parent
	bool cut_parent() {
		p->delete_child(this);
		p = NULL;
	}

	Node *parent() {
		return p;
	}
	Node *lchild() {
		return lc;
	}
	Node *rchild() {
		return rc;
	}
	Node *get_twin() {
		return twin;
	}
	int get_depth() {
		return depth;
	}
	int get_preorder_number() {
		return pre_num;
	}
	string str() {
		string s = "";
		str_hlpr(&s);
		return s;
	}
	string get_name() {
		return name;
	}

	void str_hlpr(string *s) {
		if (!name.empty())
			*s += name;
		if (contracted_lc != NULL || contracted_rc != NULL) {
			#ifdef DEBUG_CONTRACTED
				*s += "<";
			#else
				*s += "(";
			#endif
			if (contracted_lc != NULL) {
				contracted_lc->str_c_subtree_hlpr(s);
			}
			*s += ",";
			if (contracted_rc != NULL) {
				contracted_rc->str_c_subtree_hlpr(s);
			}
			#ifdef DEBUG_CONTRACTED
				*s += ">";
			#else
				*s += ")";
			#endif
		}
	}

	string str_subtree() {
		string s = "";
		str_subtree_hlpr(&s);
		return s;
	}

	void str_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			*s += "(";
			if (lc != NULL) {
				lc->str_subtree_hlpr(s);
				if (lc->parent() != this)
					cout << "#";
			}
			*s += ",";
			if (rc != NULL) {
				rc->str_subtree_hlpr(s);
				if (rc->parent() != this)
					cout << "#";
			}
			*s += ")";
		}
	}

	void str_c_subtree_hlpr(string *s) {
		str_hlpr(s);
		if (!is_leaf()) {
			#ifdef DEBUG_CONTRACTED
				*s += "<";
			#else
				*s += "(";
			#endif
			if (lc != NULL) {
				lc->str_c_subtree_hlpr(s);
			}
			*s += ",";
			if (rc != NULL) {
				rc->str_c_subtree_hlpr(s);
			}
			#ifdef DEBUG_CONTRACTED
				*s += ">";
			#else
				*s += ")";
			#endif
		}
	}

	string str_subtree_twin() {
		string s = "";
		str_subtree_twin_hlpr(&s);
		return s;
	}

	void str_subtree_twin_hlpr(string *s) {
		*s += name;//str_hlpr(s);
		if (twin != NULL) {
			*s += "{";
				twin->str_subtree_hlpr(s);
			*s += "}";
		}
		if (!is_leaf()) {
			*s += "(";
			if (lc != NULL) {
				lc->str_subtree_twin_hlpr(s);
			}
			*s += ",";
			if (rc != NULL) {
				rc->str_subtree_twin_hlpr(s);
			}
			*s += ")";
		}
	}


	void print() {
		cout << name;
	}
	void print_subtree() {
		cout << str_subtree();
		cout << endl;
	}
	void print_subtree_hlpr() {
		cout << str_subtree();
	}
	void print_subtree_twin_hlpr() {
		cout << str_subtree_twin();
	}

	bool is_leaf() {
		return (lc == NULL && rc == NULL);
	}

	bool is_sibling_pair() {
		return (lc != NULL && lc->is_leaf()
				&& rc != NULL && rc->is_leaf());
	}

	bool is_singleton() {
		return (parent() == NULL && is_leaf());
	}

	void find_sibling_pairs_hlpr(list<Node *> *sibling_pairs) {
		Node *lchild = this->lchild();
		Node *rchild = this->rchild();
		bool lchild_leaf = false;
		bool rchild_leaf = false;
		if (lchild != NULL) {
			if (lchild->is_leaf())
				lchild_leaf = true;
			else
				lchild->find_sibling_pairs_hlpr(sibling_pairs);
		}
		if (rchild != NULL) {
			if (rchild->is_leaf())
				rchild_leaf = true;
			else
				rchild->find_sibling_pairs_hlpr(sibling_pairs);
		}
		if (lchild_leaf && rchild_leaf) {
			sibling_pairs->push_back(lchild);
			sibling_pairs->push_back(rchild);
			//lchild->add_to_sibling_pairs(sibling_pairs, 1);
			//rchild->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	
	// find the sibling pairs in this node's subtree
	void append_sibling_pairs(list<Node *> *sibling_pairs) {
		find_sibling_pairs_hlpr(sibling_pairs);
	}

	// find the sibling pairs in this node's subtree
	list<Node *> find_sibling_pairs() {
		list<Node *> sibling_pairs = list<Node *>();
		find_sibling_pairs_hlpr(&sibling_pairs);
		return sibling_pairs;
	}
	
	void find_leaves_hlpr(vector<Node *> &leaves) {
		Node *lchild = this->lchild();
		Node *rchild = this->rchild();
		bool lchild_leaf = false;
		bool rchild_leaf = false;
		if (lchild != NULL) {
			if (lchild->is_leaf())
				leaves.push_back(lchild);
			else
				lchild->find_leaves_hlpr(leaves);
		}
		if (rchild != NULL) {
			if (rchild->is_leaf())
				leaves.push_back(rchild);
			else
				rchild->find_leaves_hlpr(leaves);
		}
	}
	
	// find the leaves in this node's subtree
	vector<Node *> find_leaves() {
		vector<Node *> leaves = vector<Node *>();
		if (is_leaf())
			leaves.push_back(this);
		else
			find_leaves_hlpr(leaves);
		return leaves;
	}

	// make twins point to this tree in this node's subtree
	void resync() {
		if (lc != NULL)
			lc->resync();
		if (rc != NULL)
			rc->resync();
		if (twin != NULL)
			twin->set_twin(this);
	}

	// remove all twins
	void unsync() {
		if (lc != NULL)
			lc->unsync();
		if (rc != NULL)
			rc->unsync();
		twin = NULL;
	}
	// remove all twins
	void unsync_interior() {
		if (lc != NULL)
			lc->unsync_interior();
		if (rc != NULL)
			rc->unsync_interior();
		if (lc != NULL || rc != NULL)
			twin = NULL;
	}

	// find the root of this node's tree
	Node *find_root() {
		Node *root = this;
		while (root->parent() != NULL)
			root = root->parent();
		return root;
	}

	void labels_to_numbers(map<string, int> *label_map, map<int, string> *reverse_label_map) {
		if (name != "") {
			map<string, int>::iterator i = label_map->find(name);
			stringstream ss;
			if (i != label_map->end()) {
				ss << i->second;
				name = ss.str();
			}
			else {
				int num = label_map->size();
				ss << num;
				label_map->insert(make_pair(name, num));
				reverse_label_map->insert(make_pair(num, name));
				name = ss.str();
			}
		}
		if (lc != NULL)
			lc->labels_to_numbers(label_map, reverse_label_map);
		if (rc != NULL)
			rc->labels_to_numbers(label_map, reverse_label_map);
		if (contracted_lc != NULL)
			contracted_lc->labels_to_numbers(label_map, reverse_label_map);
		if (contracted_rc != NULL)
			contracted_rc->labels_to_numbers(label_map, reverse_label_map);
	}
	
	void numbers_to_labels(map<int, string> *reverse_label_map) {
		if (name != "") {
			string converted_name = "";
			string current_num = "";
			string::iterator i = name.begin();
			size_t old_loc = 0;
			size_t loc = 0;
			while ((loc = name.find_first_of("0123456789", old_loc)) != string::npos) {
				converted_name.append(name.substr(old_loc, loc - old_loc)); 
				old_loc = loc;
				loc = name.find_first_not_of("0123456789", old_loc);
				string label = "";
				if (loc == string::npos)
					loc = name.size();
				label = name.substr(old_loc, loc - old_loc);
				map<int, string>::iterator j = reverse_label_map->find(atoi(label.c_str()));
				if (j != reverse_label_map->end()) {
					stringstream ss;
					ss << j->second;
					label = ss.str();
				}
				converted_name.append(label);
				old_loc = loc;
			}
			converted_name.append(name.substr(old_loc, name.size() - old_loc)); 
			name = converted_name;



		}
		if (lc != NULL)
			lc->numbers_to_labels(reverse_label_map);
		if (rc != NULL)
			rc->numbers_to_labels(reverse_label_map);
		if (contracted_lc != NULL)
			contracted_lc->numbers_to_labels(reverse_label_map);
		if (contracted_rc != NULL)
			contracted_rc->numbers_to_labels(reverse_label_map);
	}

	void count_numbered_labels(vector<int> *label_counts) {
		if (name != "") {
			int label = stomini(name);
			if (label_counts->size() <= label)
				label_counts->resize(label+1,0);
			(*label_counts)[label]++;
		}

		if (lc != NULL)
			lc->count_numbered_labels(label_counts);
		if (rc != NULL)
			rc->count_numbered_labels(label_counts);
		if (contracted_lc != NULL)
			contracted_lc->count_numbered_labels(label_counts);
		if (contracted_rc != NULL)
			contracted_rc->count_numbered_labels(label_counts);
	}

	void preorder_number() {
		preorder_number(0);
	}
	int preorder_number(int next) {
		set_preorder_number(next);
		next++;
		if(lchild() != NULL) {
			next = lchild()->preorder_number(next);
		}
		if(rchild() != NULL) {
			next = rchild()->preorder_number(next);
		}
		return next;
	}

	int size() {
		int s  = 1;
		if(lchild() != NULL)
			s += lchild()->size();
		if(rchild() != NULL)
			s += rchild()->size();
		return s;
	}

	int size_using_prenum() {
		if(rchild() != NULL)
			return rchild()->size_using_prenum();
		else if (lchild() != NULL)
			return lchild()->size_using_prenum();
		else
			return get_preorder_number();
	}

void add_to_front_sibling_pairs(list<Node *> *sibling_pairs, int status) {
	sibling_pairs->push_front(this);
	clear_sibling_pair(sibling_pairs);
	sibling_pair_status = status;
	sibling_pair_loc = sibling_pairs->begin();
}

void add_to_sibling_pairs(list<Node *> *sibling_pairs, int status) {
	sibling_pairs->push_back(this);
	clear_sibling_pair(sibling_pairs);
	sibling_pair_status = status;
	sibling_pair_loc = sibling_pairs->end();
	sibling_pair_loc--;
}

void remove_sibling_pair(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;

		if (sibling_loc != sibling_pairs->end()) {
			Node *old_sibling = *sibling_loc;
			old_sibling->sibling_pair_status = 0;
			sibling_pairs->erase(sibling_loc);
		}
		sibling_pairs->erase(loc);
		sibling_pair_status = 0;
	}
}

void clear_sibling_pair(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;

		if (sibling_loc != sibling_pairs->end()) {
			Node *old_sibling = *sibling_loc;
			old_sibling->sibling_pair_status = 0;
		}
		sibling_pair_status = 0;
	}
}

Node *get_sibling() {
	if (p == NULL)
		return NULL;
	else if (p->lc == this)
		return p->rc;
	else
		return p->lc;
}

Node *get_sibling(list<Node *> *sibling_pairs) {
	if (sibling_pair_status > 0) {
		list<Node *>::iterator loc = sibling_pair_loc;
		list<Node *>::iterator sibling_loc = loc;
		if (sibling_pair_status == 1)
			sibling_loc++;
		else if (sibling_pair_status == 2)
			sibling_loc--;
		return *sibling_loc;
	}
	else
		return NULL;
}

void set_sibling(Node *sibling) {
	if (sibling->sibling_pair_status > 0) {
		sibling_pair_loc = sibling->sibling_pair_loc;
		if (sibling->sibling_pair_status == 1)
			sibling_pair_loc++;
		else if (sibling->sibling_pair_status == 2)
			sibling_pair_loc--;
	}
}

void clear_sibling_pair_status() {
	sibling_pair_status = 0;
}

// fix parents
void fix_parents() {
	if (lc != NULL) {
		if (lc->parent() != this)
			lc->set_parent(this);
		lc->fix_parents();
	}
	if (rc != NULL) {
		if (rc->parent() != this)
			rc->set_parent(this);
		rc->fix_parents();
	}
}


void left_rotate() {
	if (lc->lc != NULL) {
		Node *new_lc = lc->lc;
		Node *new_rc = lc;
		Node *new_rc_lc = lc->rc;
		Node *new_rc_rc = rc;
		lc = new_lc;
		lc->p = this;
		rc = new_rc;
		rc->p = this;
		rc->lc = new_rc_lc;
		rc->lc->p = rc;
		rc->rc = new_rc_rc;
		rc->rc->p = rc;
	}
}

void right_rotate() {
	if (rc->lc != NULL) {
		Node *new_lc = rc->lc;
		Node *new_rc = rc;
		Node *new_rc_lc = rc->rc;
		Node *new_rc_rc = lc;
		lc = new_lc;
		lc->p = this;
		rc = new_rc;
		rc->p = this;
		rc->lc = new_rc_lc;
		rc->lc->p = rc;
		rc->rc = new_rc_rc;
		rc->rc->p = rc;
	}
}

void next_rooting() {
	if (lc->lc != NULL)
		left_rotate();
	else if (rc->lc != NULL)
			right_rotate();
	else
		return;
	if (lc->pre_num < rc->pre_num && (lc->pre_num != 1 || rc->pre_num != 2 || !lc->is_leaf()))
		next_rooting();
}

Node *expand_parent_edge(Node *n) {
	if (p != NULL) {
		Node *old_p = p;
		cut_parent();
		p = new Node();
		p->add_child(n);
		old_p->add_child(p);
		return p;
	}
	else {
		Node *new_child = new Node(name);
		new_child->add_child(lc);
		new_child->add_child(rc);
		new_child->contracted_lc = contracted_lc;
		new_child->contracted_rc = contracted_rc;
		name = "";
		lc = NULL;
		rc = NULL;
		contracted_lc = NULL;
		contracted_rc = NULL;
		add_child(new_child);
		return this;
	}
}

Node *undo_expand_parent_edge() {
	if (p != NULL) {
		Node *old_p = p;
		cut_parent();
		Node *child = lc;
		if (lc == NULL)
			child = rc;
		delete_child(child);
		old_p->add_child(child);
		return this;
	}
	else {
		Node *child = lc;
		if (lc == NULL)
			child = rc;
		name = child->name;
		Node *new_lc = child->lc;
		Node *new_rc = child->rc;
		contracted_lc = child->contracted_lc;
		contracted_rc = child->contracted_rc;
		child->lc = NULL;
		child->rc = NULL;
		child->contracted_lc = NULL;
		child->contracted_rc = NULL;
		child->name = "";
		set_lchild(new_lc);
		set_rchild(new_rc);
		return child;
	}
}

/*  apply an SPR operation to move this subtree to be a
 *	sibling of new_sibling
 *
 *  Note: problems will occur if new_sibling is a descendant of this
 *  Note: does nothing if this is the root
 *
 *	Returns a node that will reverse the spr (old_sibling unless it
 *		was moved to maintain a root, in which case the root is returned,
 *		NULL if no spr occured)
 */
Node *spr(Node *new_sibling, int &which_child) {
	Node *reverse;
	int prev_child_loc = 0;
	if (p == NULL || new_sibling == NULL)
		return NULL;
	Node *old_sibling = get_sibling();
	if (old_sibling == new_sibling)
		return NULL;
	Node *grandparent = p->p;
	if (p->lc == this)
		prev_child_loc = 1;
	else
		prev_child_loc = 2;
	// Prune
	if (grandparent != NULL) {
		p->delete_child(old_sibling);
		grandparent->delete_child(p);
		grandparent->add_child(old_sibling);
		reverse = old_sibling;
	}
	else {
		if (old_sibling->is_leaf())
			return NULL;
		Node *root = p;
		bool leftc = false;
		if (root->lc == this)
			leftc = true;
		root->delete_child(this);
		root->delete_child(old_sibling);
		root->add_child(old_sibling->lc);
		root->add_child(old_sibling->rc);
		old_sibling->delete_child(old_sibling->lc);
		old_sibling->delete_child(old_sibling->rc);
		if (leftc)
			old_sibling->set_lchild(this);
		else
			old_sibling->set_rchild(this);
		reverse = root;
	}


	// Regraft
	if (new_sibling->p != NULL) {
		grandparent = new_sibling->p;
		grandparent->delete_child(new_sibling);
		grandparent->add_child(p);
		p->add_child(new_sibling);
	}
	else {
		Node *root = new_sibling;
		new_sibling = p;
		p->delete_child(this);
		p->add_child(root->lc);
		p->add_child(root->rc);
		root->delete_child(root->lc);
		root->delete_child(root->rc);
		// problem here
		if (which_child == 0)
			which_child = prev_child_loc;
		if (which_child == 1) {
			root->add_child(this);
			root->add_child(new_sibling);
		}
		else {
			root->add_child(new_sibling);
			root->add_child(this);
		}


	}

	which_child = prev_child_loc;
	return reverse;
}

Node *spr(Node *new_sibling) {
	int na = 0;
	spr(new_sibling, na);
}

};

// function prototypes
Node *build_tree(string s);
Node *build_tree(string s, int start_depth);
int build_tree_helper(int start, const string& s, Node *parent,
		bool &valid);
//void preorder_number(Node *node);
//int preorder_number(Node *node, int next);


// build a tree from a newick string
Node *build_tree(string s) {
	return build_tree(s, 0);
}
Node *build_tree(string s, int start_depth) {
	if (s == "")
		return new Node();
	Node *dummy_head = new Node("p", start_depth-1);
	bool valid = true;
	build_tree_helper(0, s, dummy_head, valid);
	Node *head = dummy_head->lchild();
	if (valid) {
		delete dummy_head;
		return head;
	}
	else {
		head->delete_tree();
		return dummy_head;
	}

}

// build_tree recursive helper function
int build_tree_helper(int start, const string& s, Node *parent,
		bool &valid) {
	int loc = s.find_first_of("(,)", start);
	if (loc == string::npos) {
		string name = s.substr(start, s.size() - start);
		Node *node = new Node(name);
		parent->add_child(node);
		loc = s.size()-1;
		return loc;
	}
	while(s[start] == ' ' || s[start] == '\t')
		start++;
	int end = loc;
	while(s[end] == ' ' || s[end] == '\t')
		end--;
	string name = s.substr(start, end - start);
	Node *node = new Node(name);
	parent->add_child(node);
	if (s[loc] == '(') {
			loc = build_tree_helper(loc + 1, s, node, valid);
			loc = build_tree_helper(loc + 1, s, node, valid);
//			int loc_check = s.find_first_of("(,)", loc);
//			if (loc_check != string::npos &&
//					s[loc_check] == ','
			if (s[loc] != ')') {
				valid = false;
					return s.size()-1;
			}
			//loc = build_tree_helper(loc + 1, s, node);
			loc++;
	}
	return loc;
}

// swap two nodes
void swap(Node **a, Node **b) {
	Node *temp = *a;
	*a = *b;
	*b = temp;
}

// expand all contracted nodes of a subtree starting at n
void expand_contracted_nodes(Node *n) {
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	if (lc != NULL)
		expand_contracted_nodes(lc);
	if (rc != NULL)
		expand_contracted_nodes(rc);
	if (n->is_leaf()) {
		Node *subtree = build_tree(n->str(), n->get_depth());
		Node *new_lc = subtree->lchild();
		Node *new_rc = subtree->rchild();
		if (new_lc != NULL) {
			n->set_lchild(new_lc);
			new_lc->set_parent(n);
			n->set_name("");
		}
		if (new_rc != NULL) {
			n->set_rchild(new_rc);
			new_rc->set_parent(n);
			n->set_name("");
		}
		delete subtree;
	}
}


/*
void preorder_number(Node *node) {
	preorder_number(node, 0);
}
int preorder_number(Node *node, int next) {
	node->set_preorder_number(next);
	next++;
	if(node->lchild() != NULL) {
		next = preorder_number(node->lchild(), next);
	}
	if(node->rchild() != NULL) {
		next = preorder_number(node->rchild(), next);
	}
	return next;
}
*/
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

string root(string s) {
	string r = "";
	int i = 0;
	int depth = 0;
	int first_c = -1;
	int second_c = -1;
	int last_bracket = -1;
	for(int i = 0; i < s.size(); i++) {
		if (s[i] == '(')
			depth++;
		else if (s[i] == ')') {
			depth--;
			last_bracket = i;
		}
		else if (depth == 1 && s[i] == ',') {
			if (first_c == -1)
				first_c = i;
			else if  (second_c == -1)
				second_c = i;
		}
	}
	if (second_c == -1 || last_bracket == -1)
		return s;
	else {
		r.append(s.substr(0,first_c+1));
		r.append("(");
		r.append(s.substr(first_c+1,last_bracket-first_c));
		r.append(")");
		r.append(s.substr(last_bracket+1,string::npos));
	}
//	cout << r << endl;
	return r;
}

#endif
