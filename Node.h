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

#ifndef INCLUDE_CSTDIO
	#define INCLUDE_CSTDIO
	#include <cstdio>
#endif
#ifndef INCLUDE_STRING
	#define INCLUDE_STRING
	#include <string>
#endif
#ifndef INCLUDE_IOSTREAM
	#define INCLUDE_IOSTREAM
	#include <iostream>
#endif
#ifndef INCLUDE_SSTREAM
	#define INCLUDE_SSTREAM
	#include <sstream>
#endif
#ifndef INCLUDE_MAP
	#define INCLUDE_MAP
	#include <map>
#endif
using namespace std;

// representation of a component with no leaves
#define DEAD_COMPONENT "*"
/*void find_sibling_pairs_hlpr(Node *node, vector<Node *> &sibling_pairs);
void find_leaves_hlpr(Node *node, vector<Node *> &leaves);
void str_subtree_hlpr(string *s);
vector<Node *> find_sibling_pairs(Node *node);
vector<Node *> find_leaves(Node *node);
*/
class Node {
	private:
	Node *lc;			// left child
	Node *rc;			// right child
	Node *p;			// parent
	Node *twin;			// counterpart in another tree
	string name;		// label
	int depth;			//distance from root
	int pre_num;	// preorder number

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
	}
	// copy constructor
	Node(const Node &n) {
		p = NULL;
		name = n.name;
		twin = n.twin;
		depth = n.depth;
		if (n.lc == NULL)
			lc = NULL;
		else
			lc = new Node(*n.lc, this);
		if (n.rc == NULL)
			rc = NULL;
		else
			rc = new Node(*n.rc, this);
	}
	Node(const Node &n, Node *parent) {
		p = parent;
		name = n.name;
		twin = n.twin;
		depth = n.depth;
		if (n.lc == NULL)
			lc = NULL;
		else
			lc = new Node(*n.lc, this);
		if (n.rc == NULL)
			rc = NULL;
		else
			rc = new Node(*n.rc, this);
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
	}
	// delete a subtree
	void delete_tree() {
		Node *lc = this->lc;
		Node *rc = this->rc;
		delete this;
		if (lc != NULL) {
			lc->delete_tree();
		}
		if (rc != NULL) {
			rc->delete_tree();
		}
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
		n->p = this;
		n->depth = depth+1;
		return lc;
	}
	Node *set_rchild(Node *n) {
		this->rc = n;
		n->p = this;
		n->depth = depth+1;
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
	int set_preorder_number(int p) {
		pre_num = p;
		return pre_num;
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
	Node *contract() {
		Node *parent = p;
		Node *child;
		Node *ret = NULL;
		// contract out this node and give child to parent
		if (p != NULL) {
		//cout << p->str_subtree() << endl;
			if (lc && !rc) {
				child = lc;
				delete this;
				parent->add_child(child);
				ret = parent;
			}
			else if (rc && !lc) {
				child = rc;
				delete this;
				parent->add_child(child);
				ret = parent;
			}
			else if (lc == NULL && rc == NULL) {
				delete this;
				ret = parent->contract();
			}
			//else if (!(lc && rc)) {
			else
				ret = this;
		}
		// if no parent then take children of single child and remove it
		else {
			// dead component or singleton, will be cleaned up by the forest
			if (lc == NULL && rc == NULL) {
				if (name == "")
					name = DEAD_COMPONENT;
				
			}
			//else if (!(lc && rc)) {
			else if ((bool)lc xor (bool)rc) {
				child = lc;
				if (child == NULL)
					child = rc;
		// if child is a leaf then get rid of this so we don't lose refs
				// problem: if the child is not c, then we want to copy
				// otherwise we don't
				// copy other parameters and join the twin
				//to this if the child is a label
				Node *new_lc = child->lchild();
				Node *new_rc = child->rchild();
				if (child->is_leaf()) {
					set_twin(child->get_twin());
					child->get_twin()->set_twin(this);
					name = child->str();
				}
				delete child;
				if (new_lc != NULL)
					add_child(new_lc);
				if (new_rc != NULL)
					add_child(new_rc);
				ret = this;
			}
		}

		return ret;
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
		return name;
	}
	string str_subtree() {
		string s = "";
		str_subtree_hlpr(&s);
		return s;
	}

	void str_subtree_hlpr(string *s) {
		*s += name;
		if (!is_leaf()) {
			*s += "(";
			if (lc != NULL) {
				lc->str_subtree_hlpr(s);
			}
			*s += ",";
			if (rc != NULL) {
				rc->str_subtree_hlpr(s);
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

	void find_sibling_pairs_hlpr(vector<Node *> &sibling_pairs) {
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
			sibling_pairs.push_back(lchild);
			sibling_pairs.push_back(rchild);
		}
	}
	
	// find the sibling pairs in this node's subtree
	vector<Node *> find_sibling_pairs() {
		vector<Node *> sibling_pairs = vector<Node *>();
		find_sibling_pairs_hlpr(sibling_pairs);
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
	}
};

