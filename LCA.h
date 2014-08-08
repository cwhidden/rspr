/*******************************************************************************
LCA.h

Data structure for LCA computations on a binary tree
Implementation of the RMQ-based methods of Bender and Farach-Colton

Copyright 2010-2014 Chris Whidden
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

#ifndef INCLUDE_LCA

#define INCLUDE_LCA
#include <cstdio>
#include <string>
#include <iostream>
#include "Node.h"
#include <vector>
#include <cmath>
using namespace std;

int mylog2 (int val) {
	if (val <= 0)
		return -1;
	return 31 - __builtin_clz(val);
	cout << "start" << ",";
	cout << val << ",";
	cout <<  31 - __builtin_clz(val) << ",";;
    int ret = -1;
    while (val > 0) {
        val >>= 1;
        ret++;
    }
		cout << ret << endl;
		cout << endl;
    return ret;
}

class LCA {
	private:
	Node *tree;
	vector<int> E;		// preorder numbers of euler tour
	vector<int> L;		// levels of euler tour
	vector<int> H;		// first occurence of a preorder number in E
	vector<int> T;    // real preorder to internal preorder mapping
	vector<Node *> N;	// preorder to node mapping
	vector<vector<int> > RMQ;	// precomputed RMQ values

	public:
	LCA(Node *tree) {
		this->tree = tree;
		T = vector<int>();
		if (tree->get_preorder_number() == -1)
			tree->preorder_number();
		euler_tour(tree, 0);
		precompute_rmq();
	}

	LCA() {
		this->tree = NULL;
	}
	void euler_tour(Node *node, int depth) {
		// First visit
		int preorder_number = N.size();
		int euler_number = E.size();
		N.push_back(node);
		if (T.size() <= node->get_preorder_number())
			T.resize(node->get_preorder_number()+1,-1);
		T[node->get_preorder_number()] = preorder_number;

		//cout << preorder_number << "\t";
		//node->print_subtree();
		H.push_back(euler_number);
		L.push_back(depth);
		E.push_back(preorder_number);

		// TODO: check that this is correct for multifurcating trees
		
		list<Node *>::const_iterator c;
		for(c = node->get_children().begin(); c != node->get_children().end();
				c++) {
			euler_tour(*c, depth+1);
			// Middle/Last visit
			L.push_back(depth);
			E.push_back(preorder_number);
		}

	}

	void precompute_rmq() {
		// E gives queries of length 1
		RMQ.push_back(E);
		for(int length = 1; length < E.size(); length *= 2) {
			vector<int> V = vector<int>();
			for(int start = 0; start < E.size() - length; start++) {
				V.push_back(rmq(start, start + length));
			}
			RMQ.push_back(V);
		}
	}
	//vector<vector<int> > RMQ;	// precomputed RMQ values

	// find the index of the rmq between indices i and j of E
	int rmq(int i, int j) {
		if (i == j)
			return E[i];
		int length = j - i - 1;
		//cout << "  " << i << " " << j << endl;
		//cout << "  " << length << endl;
		//cout << "  " << mylog2(length) << endl;
		length = (mylog2(length));
		//cout << "  " << length << endl;
		//cout << "  " << j - (1 << (length)) << endl;
		int rmq1 = RMQ[length+1][i];
		int rmq2;
		if (length >= 0)
			rmq2 = RMQ[length+1][j - (1 << (length))];
		else
			rmq2 = RMQ[length+1][j];
		if (rmq1 < rmq2)
			return rmq1;
		return rmq2;
	}
	Node *get_lca(Node *a, Node *b) {
		int preorder_a = T[a->get_preorder_number()];
		int preorder_b = T[b->get_preorder_number()];
		int lca_index;
		if (preorder_a <= preorder_b)
			lca_index = rmq(H[preorder_a], H[preorder_b]);
		else
			lca_index = rmq(H[preorder_b], H[preorder_a]);
		return N[lca_index];
	}

	Node *get_tree() {
		return tree;
	}

	/* copy constructor
	LCA(const LCA &n) {
	}
	*/
	void debug() {
		for(int i = 0; i < E.size(); i++) {
			cout << " " << E[i];
		}
		cout << endl;
		cout << endl;
		for(int i = 0; i < L.size(); i++) {
			cout << " " << L[i];
		}
		cout << endl;
		cout << endl;
		for(int i = 0; i < H.size(); i++) {
			cout << " " << H[i];
		}
		cout << endl;
		cout << endl;
		cout << endl;
		for(int i = 0; i < RMQ.size(); i++) {
			for(int j = 0; j < RMQ[i].size(); j++) {
			cout << " " << RMQ[i][j];
			}
			cout << endl;
			cout << endl;
		}
	}
};

#endif
