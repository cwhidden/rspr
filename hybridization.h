/*******************************************************************************
hybridization.h

Data structure for hybridization
Contains methods to detect agreement forest cycles

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

#ifndef INCLUDE_HYBRIDIZATION

#define INCLUDE_HYBRIDIZATION

#include <cstdio>
#include <string>
#include <iostream>
#include <list>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
//#include <boost/tuple/tuple.hpp>
//#include <boost/graph/visitors.hpp>
//#include <boost/graph/graph_utility.hpp>
using namespace std;

/*struct Vertex {
	string name;
	enum COLOUR colour;
	enum COLOUR new_colour;

};
struct Edge {
	int weight;
	enum COLOUR colour;
};
*/
typedef boost::adjacency_list<
	boost::vecS, boost::vecS, boost::directedS,
	boost::no_property, boost::no_property> Graph;

struct cycle_detector : public boost::dfs_visitor<>
  {
    cycle_detector( bool& has_cycle) 
      : _has_cycle(has_cycle) { }

    template <class Edge, class Graph>
    void back_edge(Edge, Graph&) {
      _has_cycle = true;
    }
  protected:
    bool& _has_cycle;
  };

// function prototypes
bool detect_cycle(Node *T1, Node *T2, Forest *AF);
void add_AF_edges(Node *n, Forest *AF, vector<int> *leaf_cnumber,
		vector<int> *node_cnumber, vector<int> *node_lcount,
		vector<list<int> > *node_lists, vector<int> *component_lcount,
		Graph *G);
bool has_cycle(Graph *g);

bool detect_cycle(Node *T1, Node *T2, Forest *AF) {
//	cout << "BEGIN detect_cycle()" << endl;
	// the component number of a leaf/node
	vector<int> leaf_cnumber = vector<int>();
	vector<int> node_cnumber = vector<int>();
	// the number of leaves below a node
	vector<int> node_lcount = vector<int>();
	// the list of descendants of a node that have the same hybrid edge root
	vector<list<int> > node_lists = vector<list<int> >();
	vector<int> component_lcount = vector<int>(AF->size(), 0);
	// number the leaves of T1 with their cnumber and lnumber
	for(int i = 0; i < AF->size(); i++) {
//		cout << "component: " << i << endl;
		vector<Node *> leaves = AF->get_component(i)->find_leaves();
//		cout << leaves.size() << " leaves" << endl;
		vector<Node *>::iterator leaf;
		for(leaf = leaves.begin(); leaf != leaves.end(); leaf++) {
			string name = (*leaf)->str();
//			cout << "\tleaf: " << name << endl;
			int number = atoi(name.c_str());
			if (number >= leaf_cnumber.size())
				leaf_cnumber.resize(number+1);
			leaf_cnumber[number] = i;
			component_lcount[i]++;
		}
		// for each leaf
		// number its component in leaf_cnumber
	}
	// preorder_number the trees
	T1->preorder_number();
	T2->preorder_number();

	// create the graph 
	Graph G = Graph(AF->size());
//	cout << "Adding T1 edges" << endl;
	add_AF_edges(T1, AF, &leaf_cnumber, &node_cnumber, &node_lcount,
		&node_lists, &component_lcount, &G);
//	cout << "Adding T2 edges" << endl;
	add_AF_edges(T2, AF, &leaf_cnumber, &node_cnumber, &node_lcount,
		&node_lists, &component_lcount, &G);
	// check for a cycle

//	cout << "END detect_cycle()" << endl;
	bool found_cycle = has_cycle(&G);
	//write_graphviz(cout, G);
	return found_cycle;
}

void add_AF_edges(Node *n, Forest *AF, vector<int> *leaf_cnumber,
		vector<int> *node_cnumber, vector<int> *node_lcount,
		vector<list<int> > *node_lists, vector<int> *component_lcount,
		Graph *G) {
//	cout << "BEGIN add_AF_edges()" << endl;
//	cout << "1" << endl;
	if (n == NULL)
		return;
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	int n_cnumber, lc_cnumber, rc_cnumber;
	int n_lcount, lc_lcount, rc_lcount;
	list<int> n_list, lc_list, rc_list;
	n_list = list<int>();
	string n_name = n->str();
	int n_name_number = atoi(n_name.c_str());
	int n_number = n->get_preorder_number();
	// find lc's cnumber and lcount
//	cout << "2" << endl;
	if (lc != NULL) {
		// recurse on the lc so it's values will be in node_*
		add_AF_edges(lc, AF, leaf_cnumber, node_cnumber, node_lcount,
			node_lists, component_lcount,  G);
		string lc_name = lc->str();
		int lc_number = lc->get_preorder_number();
		lc_cnumber = (*node_cnumber)[lc_number];
		lc_lcount = (*node_lcount)[lc_number];
		lc_list = (*node_lists)[lc_number];
	}
	else {
		lc_cnumber = -1;
		lc_lcount = 0;
		lc_list = list<int>();
	}
//	cout << "3" << endl;
	// find rc's cnumber and lcount
	if (rc != NULL) {
		// recurse on the rc so it's values will be in node_*
		add_AF_edges(rc, AF, leaf_cnumber, node_cnumber, node_lcount,
			node_lists, component_lcount,  G);
		string rc_name = rc->str();
		int rc_number = rc->get_preorder_number();
		rc_cnumber = (*node_cnumber)[rc_number];
		rc_lcount = (*node_lcount)[rc_number];
		rc_list = (*node_lists)[rc_number];
	}
	else {
		rc_cnumber = -1;
		rc_lcount = 0;
		rc_list = list<int>();
	}
//	cout << "4" << endl;
	// lookup values if this is a leaf
	if (lc == NULL && rc == NULL && leaf_cnumber->size() > n_name_number) {
		n_cnumber = (*leaf_cnumber)[n_name_number];
		n_lcount = 1;
		if ((*component_lcount)[n_cnumber] == 1)
			n_list.push_back(n_number);
	}
	// use other child's values if one is NULL
	else if (lc == NULL) {
		n_cnumber = rc_cnumber;
		n_lcount = rc_lcount;
		n_list.splice(n_list.begin(), rc_list);
	}
	else if (rc == NULL) {
		n_cnumber = lc_cnumber;
		n_lcount = lc_lcount;
		n_list.splice(n_list.begin(), lc_list);
	}
	else {
		// same component
		if (lc_cnumber != -1 && lc_cnumber == rc_cnumber) {
			n_cnumber = lc_cnumber;
			n_lcount = lc_lcount + rc_lcount;
			if (n_lcount == (*component_lcount)[n_cnumber]) {
				n_list.push_back(n_number);
				//cout << "component " << n_cnumber << " " << n_lcount << endl;
			}
//			cout << "\tsame component" << endl;
		}
		// both children are finished components
		else if ((lc_cnumber == -1 ||
				lc_lcount == (*component_lcount)[lc_cnumber]) &&
				(rc_cnumber == -1 ||
				rc_lcount == (*component_lcount)[rc_cnumber])) {
			//cout << "components: " << lc_cnumber << " " << rc_cnumber << endl;
			//cout << "leaves: " << lc_lcount << " " << rc_lcount << endl;
			//cout << n_list.size() << " + " << lc_list.size() << " + " << rc_list.size()
				//<< " = ";
			n_list.splice(n_list.begin(), lc_list);
			n_list.splice(n_list.begin(), rc_list);
			//cout << n_list.size() << "\n";
			n_cnumber = -1;
			n_lcount = 0;
		}
		// otherwise, one of the children must be a finished component
		else if (lc_cnumber == -1 ||
				lc_lcount == (*component_lcount)[lc_cnumber]) {
			// rc's component is the parent of each list component
			list<int>::iterator i;
			//cout << "lc\n";
			for(i = lc_list.begin(); i != lc_list.end(); i++) {
				boost::add_edge(rc_cnumber, (*node_cnumber)[*i], *G);
				//cout << "adding edge (" << rc_cnumber << "," << (*node_cnumber)[*i]
			//		<< ")\n";
			}
			n_cnumber = rc_cnumber;
			n_lcount = rc_lcount;
//			cout << "\trc->lc" << endl;
		}
		else if (rc_cnumber == -1 ||
				rc_lcount == (*component_lcount)[rc_cnumber]) {
			// lc's component is the parent of each list component
			//cout << "rc\n";
			list<int>::iterator i;
			for(i = rc_list.begin(); i != rc_list.end(); i++) {
				boost::add_edge(lc_cnumber, (*node_cnumber)[*i], *G);
				//cout << "adding edge (" << lc_cnumber << "," << (*node_cnumber)[*i]
					//<< ")\n";
			}
			n_cnumber = lc_cnumber;
			n_lcount = lc_lcount;
//			cout << "\tlc->cc" << endl;
		}
		else {
			// error
//			cout << "error" << endl;
		}
	}
//	cout << "5" << endl;
	// add the values to node_cnumber and node_lcount
	if (n_number >= node_cnumber->size())
		node_cnumber->resize(n_number+1);
	(*node_cnumber)[n_number] = n_cnumber;
	if (n_number >= node_lcount->size())
		node_lcount->resize(n_number+1);
	(*node_lcount)[n_number] = n_lcount;
	if (n_number >= node_lists->size())
		node_lists->resize(n_number+1);
	(*node_lists)[n_number] = list<int>(n_list);
	//cout << n_list.size() << endl;
	//cout << (*node_lists)[n_number].size() << endl;

//	cout << "END add_AF_edges()" << endl;
}

bool has_cycle(Graph *g) {
	bool has_cycle = false;
  	cycle_detector vis(has_cycle);
  	boost::depth_first_search(*g, visitor(vis));
	return has_cycle; 	
}

#endif
