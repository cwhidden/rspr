/*******************************************************************************
ClusterForest.h


Copyright 2011 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
March 7, 2010
Version 1.02

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
	#include "Forest.h"
#endif
using namespace std;

class ClusterForest: public Forest {
	public:
		vector<Node *> cluster_nodes;

	public:
	ClusterForest() : Forest(){
		init();
	}
	ClusterForest(vector<Node *> components) : Forest(components) {
		init();
	}
	ClusterForest(Node *head) : Forest(head){
		init();
	}
	ClusterForest(const ClusterForest &f) : Forest(f) {
		cluster_nodes = vector<Node *>(f.cluster_nodes);
	}

	void init() {
		cluster_nodes = vector<Node *>();
		cluster_nodes.push_back(NULL);
	}
	~ClusterForest() {
		/* don't delete cluster nodes because they are, by
		   definition, in the components
		*/
		cluster_nodes.clear();
	} 

	// add a new cluster at n
	void add_cluster(Node *n, string name) {
		Node *n_parent = n->parent();
		n->cut_parent();
		Node *n_child = new Node(name);
		n_parent->add_child(n_child);
		add_component(n);
		cluster_nodes.push_back(n_child);
	}

	// swap the contents of two forests
	void swap(ClusterForest *f) {
		vector<Node *> components_temp = this->components;
		this->components = f->components;
		f->components = components_temp;
		
		/*
		vector<Node *> deleted_nodes_temp = this->deleted_nodes;
		this->deleted_nodes = f->deleted_nodes;
		f->deleted_nodes = deleted_nodes_temp;
		*/

		vector<Node *> cluster_nodes_temp = this->cluster_nodes;
		this->cluster_nodes = f->cluster_nodes;
		f->cluster_nodes = cluster_nodes_temp;
	}

	inline Node *get_cluster_node(int i) {
		return cluster_nodes[i];
	}
	inline int num_clusters() {
		return cluster_nodes.size();
	}
	void join_cluster(int cluster_loc, Forest *solved_cluster) {
		print_components();
		cout << "a" << endl;
		Node *cluster_node = get_cluster_node(cluster_loc);
		cout << "b" << endl;
		Node *cluster_parent = cluster_node->parent();
		cout << "c" << endl;
		cluster_node->cut_parent();
		cout << "d" << endl;
		delete cluster_node;
		delete components[cluster_loc];
		components[cluster_loc] = NULL;
		cout << "e" << endl;
		int start = 0;
		cout << "f" << endl;
		// problem here!!!
		if (solved_cluster->contains_rho()) {
			cluster_parent->contract();
		}
		else {
			cluster_parent->add_child(new Node(*(solved_cluster->get_component(0))));
			start = 1;
		}
		cout << "g" << endl;
		// should we add these to a finished_components or something?
		for(int i = start; i < solved_cluster->num_components(); i++) {
			if (solved_cluster->get_component(i)->str() != "p")
				add_component(new Node(*(solved_cluster->get_component(i))));
		}
		cout << "h" << endl;
		print_components();

	}

};

// Functions

// swap two cluster forests
void swap(ClusterForest **a, ClusterForest **b) {
	(*a)->swap(*b);
}
