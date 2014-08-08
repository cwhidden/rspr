/*******************************************************************************
SiblingPair.h

Data structure for a sibling pair of a binary tree

Copyright 2012-2014 Chris Whidden
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

#ifndef INCLUDE_SIBLINGPAIR

#define INCLUDE_SIBLINGPAIR

#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <set>
#include "Node.h"

using namespace std;

class SiblingPair {
	public:
	Node *a;
	Node *c;
	int key;
	int key2;

	SiblingPair() {
		a = NULL;
		c = NULL;
		key = -1;
		key2 = -1;
	}
	SiblingPair(Node *A, Node *C) {
		a = A;
		c = C;
		key = a->get_preorder_number();
		if (c->get_preorder_number() < key)
			key = c->get_preorder_number();
	//	a->get_preorder_number();
	//	if (c->get_preorder_number() < key)
	//		key = c->get_preorder_number();
	}
	bool operator< (const SiblingPair &sp) const {
		return key < sp.key;
	}
};

	// TODO: binary only
	void find_sibling_pairs_set_hlpr(Node *n,
			set<SiblingPair> *sibling_pairs) {
		Node *lchild = n->lchild();
		Node *rchild = n->rchild();
		bool lchild_leaf = false;
		bool rchild_leaf = false;
		if (lchild != NULL) {
			if (lchild->is_leaf())
				lchild_leaf = true;
			else
				find_sibling_pairs_set_hlpr(lchild,sibling_pairs);
		}
		if (rchild != NULL) {
			if (rchild->is_leaf())
				rchild_leaf = true;
			else
				find_sibling_pairs_set_hlpr(rchild,sibling_pairs);
		}
		if (lchild_leaf && rchild_leaf) {
			sibling_pairs->insert(SiblingPair(lchild,rchild));
			//lchild->add_to_sibling_pairs(sibling_pairs, 1);
			//rchild->add_to_sibling_pairs(sibling_pairs, 2);
		}
	}
	
	// find the sibling pairs in this node's subtree
	void append_sibling_pairs_set(Node *n,set<SiblingPair> *sibling_pairs) {
		find_sibling_pairs_set_hlpr(n,sibling_pairs);
	}

	// find the sibling pairs in this node's subtree
	set<SiblingPair> *find_sibling_pairs_set(Node *n) {
		set<SiblingPair> *sibling_pairs = new set<SiblingPair>();
		find_sibling_pairs_set_hlpr(n,sibling_pairs);
		return sibling_pairs;
	}

	// return a set of the sibling pairs
	set<SiblingPair> *find_sibling_pairs_set(Forest *f) {
		set<SiblingPair> *sibling_pairs = new set<SiblingPair>();
		for(int i = 0; i < f->num_components(); i++) {
			Node *component = f->get_component(i);
			append_sibling_pairs_set(component,sibling_pairs);
		}
		return sibling_pairs;
	}
#endif
