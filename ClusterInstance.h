/*******************************************************************************
ClusterInstance.h


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
#ifndef INCLUDE_CLUSTERINSTANCE

#define INCLUDE_CLUSTERINSTANCE
#include "Forest.h"
#include "Node.h"

using namespace std;

class ClusterInstance {
	public:
		Forest *F1;
		Forest *F2;
		Node *F1_cluster_node;
		Node *F2_cluster_node;
		bool F2_has_component_zero;

	public:
	ClusterInstance(){
		init(NULL, NULL, NULL, NULL, false);
	}
	ClusterInstance(Forest *f1, Forest *f2, Node *f1_cluster_node,
			Node *f2_cluster_node,bool f2_has_component_zero) {
		init (f1, f2, f1_cluster_node, f2_cluster_node, f2_has_component_zero);
	}
	ClusterInstance(const ClusterInstance &c) {
		init(c.F1,c.F2,c.F1_cluster_node,c.F2_cluster_node,c.F2_has_component_zero);
	}
	~ClusterInstance() {
	}

	void init(Forest *f1, Forest *f2, Node *f1_cluster_node,
			Node *f2_cluster_node, bool f2_has_component_zero) {
		F1 = f1;
		F2 = f2;
		F1_cluster_node = f1_cluster_node;
		F2_cluster_node = f2_cluster_node;
		F2_has_component_zero = f2_has_component_zero;
	}

	bool is_original() {
		if (F1_cluster_node == NULL && F2_cluster_node == NULL)
			return true;
		else
			return false;
	}

	int join_cluster(Forest *original_F1, Forest *original_F2) {
		int adjustment = 0;
		int rho_1 = -1;
		int rho_2 = -1;
		int start = 0;
		Forest *F1_cluster_node_forest = NULL;
		Forest *F2_cluster_node_forest = NULL;
		#ifdef DEBUG_CLUSTERS
			cout << "join_cluster" << endl;
			cout << "original_F1: ";
			original_F1->print_components();
			cout << "original_F2: ";
			original_F2->print_components();
			cout << "cluster_F1: ";
			F1->print_components();
			cout << "cluster_F2: ";
			F2->print_components();
			if (F1_cluster_node != NULL) {
				cout << "F1_cluster_node: ";
				F1_cluster_node->print_subtree();
				cout << "F1_cluster_node_forest: ";
				F1_cluster_node->get_forest()->print_components();
				F1_cluster_node_forest = F1_cluster_node->get_forest();
			}
			if (F2_cluster_node != NULL) {
				cout << "F2_cluster_node: ";
				F2_cluster_node->print_subtree();
				cout << "F2_cluster_node_forest: ";
				F2_cluster_node->get_forest()->print_components();
				F2_cluster_node_forest = F2_cluster_node->get_forest();
			}
		#endif
		if (F1_cluster_node != NULL) {
			F1_cluster_node->decrease_clustered_children();
			F1_cluster_node_forest = F1_cluster_node->get_forest();
		}
		if (F2_cluster_node != NULL)
			F2_cluster_node->decrease_clustered_children();

		if (F1->contains_rho()) {
			bool contract = true;
			if (F1_cluster_node->is_leaf()) {
				if (F1_cluster_node->get_num_clustered_children() >= 1) {
					contract = false;
				}
				else {
					F1_cluster_node->add_child(F1->get_component(0));
					start = 1;
				}
			}
			else {
				F1_cluster_node->get_forest()->add_component(F1->get_component(0));
				start = 1;
			}
			if (F1_cluster_node->parent() == NULL) {
				if (F1_cluster_node->lchild() != NULL &&
					 F1_cluster_node->lchild()->get_num_clustered_children() > 0) {
					F1_cluster_node_forest->update_component(
							F1_cluster_node, F1_cluster_node->lchild());
				}
				if (F1_cluster_node->rchild() != NULL &&
					 F1_cluster_node->rchild()->get_num_clustered_children() > 0) {
					F1_cluster_node_forest->update_component(
							F1_cluster_node, F1_cluster_node->rchild());
				}
			}
			if (contract)
				F1_cluster_node->contract();
		}
		else {
			F1_cluster_node->add_child(F1->get_component(0));
			if (F1_cluster_node->get_num_clustered_children() <= 0)
				F1_cluster_node->contract();
			start = 1;
		}
		// should we add these to a finished_components or something?
		for(int i = start; i < F1->num_components(); i++) {
			if (F1->get_component(i)->str() != "p")
				original_F1->add_component(F1->get_component(i));
			else {
				rho_1 = i;
			}
		}

		bool skip = false;
		Node *F2_cluster = F1->get_component(0)->get_twin();
		// TODO: this is still not quite right. We should only add the cluster
		//to the front if that is where it came from
//		cout << "F1_rho=" << F1->contains_rho() << endl;
//		cout << "F2_rho=" << F2->contains_rho() << endl;
		if ((F1->contains_rho() || F2->contains_rho()) && !F2_has_component_zero && F2_cluster_node != NULL) {
//			if (F2_has_component_zero) {
//				cout << "PROBLEM!!" << endl;
				//if (!original_F1->contains_rho())
				//	original_F1->add_rho();
				//if (!original_F2->contains_rho())
					//original_F2->add_rho();
//			}
			bool contract = true;
			if (F2_cluster_node->is_leaf()) {
				if (F2_cluster_node->get_num_clustered_children() >= 1) {
					contract = false;
				}
				else {
					//F2_cluster_node->add_child(F2->get_component(0));
					F2_cluster_node->add_child(F2_cluster);
					skip = true;
				}
			}
			if (F2_cluster_node->parent() == NULL) {
				if (F2_cluster_node->lchild() != NULL &&
					 F2_cluster_node->lchild()->get_num_clustered_children() > 0) {
					F2_cluster_node->get_forest()->update_component(
							F2_cluster_node, F2_cluster_node->lchild());
				}
				if (F2_cluster_node->rchild() != NULL &&
					 F2_cluster_node->rchild()->get_num_clustered_children() > 0) {
					F2_cluster_node->get_forest()->update_component(
							F2_cluster_node, F2_cluster_node->rchild());
				}
			}
			if (contract) {
				F2_cluster_node->contract();
				if (F1_cluster_node_forest->get_twin()->get_component(0)->str() == F2_cluster->str() && F1_cluster_node_forest->get_component(0)->str() == F2_cluster->str() && !F1_cluster_node_forest->contains_rho()) {
					F1_cluster_node_forest->add_rho();
					F1_cluster_node_forest->get_twin()->add_rho();
					adjustment++;
				}
			}
		}
		else { // if (F2_cluster != NULL) {
			if (F2_cluster == NULL)
				F2_cluster = F2->get_component(0);
			skip = true;
//			cout << "skip=" << skip << endl;
			// TODO: this is not right yet
			if (F2_has_component_zero) {
//				cout << __LINE__ << endl;
				F1_cluster_node_forest->get_twin()->add_component(0, F2_cluster);
				if (F1->contains_rho() || F2->contains_rho()) {
					if (!F1_cluster_node_forest->contains_rho())
						F1_cluster_node_forest->add_rho();
					if (!F1_cluster_node_forest->get_twin()->contains_rho())
						F1_cluster_node_forest->get_twin()->add_rho();
				}
			}
			else if (F2_cluster_node == NULL) {
				skip = false;
			}
			else {
				// TODO: check this
				F2_cluster_node->add_child(F2->get_component(0));
				F2_cluster = F2->get_component(0);
				if ((F2_cluster_node->lchild() == NULL || F2_cluster_node->rchild() == NULL) && F2_cluster_node->get_num_clustered_children() <= 0) {
					F2_cluster_node->contract();
				}
			}
		}
		// should we add these to a finished_components or something?
		for(int i = 0; i < F2->num_components(); i++) {
			if (F2->get_component(i) == F2_cluster) {
				if (!skip) {
					if (F2->get_component(i)->str() != "p")
						F1_cluster_node_forest->get_twin()->add_component(F2->get_component(i));
					else
						rho_2 = i;
				}
			}
			else {
				if (F2->get_component(i)->str() != "p")
					original_F2->add_component(F2->get_component(i));
				else
					rho_2 = i;
			}
		}

		/*
		cout << "Joined" << endl;
		cout << F1_cluster_node << endl;
		cout << F1->get_component(0) << endl;
		if (!F1->contains_rho() && F1_cluster_node != NULL)
			cout << F1_cluster_node->str_subtree() << endl;
		else {
			cout << F1->get_component(0)->str_subtree() << endl;
		}
		if (!F1->contains_rho() && !F2->contains_rho() && F2_cluster_node != NULL)
			cout << F2_cluster_node->str_subtree() << endl;
		else
			cout << F2_cluster->str_subtree() << endl;
			*/

		if (rho_1 >= 0) {
			F1->get_component(rho_1)->delete_tree();
			F1->set_component(rho_1,NULL);
		}
		if (rho_2 >= 0) {
			F2->get_component(rho_2)->delete_tree();
			F2->set_component(rho_2,NULL);
		}
		F1->erase_components();
		F2->erase_components();

		#ifdef DEBUG_CLUSTERS
			cout << "join_cluster finished" << endl;
			cout << "original_F1: ";
			original_F1->print_components();
			cout << "original_F2: ";
			original_F2->print_components();
			if (F1_cluster_node_forest != NULL) {
				cout << "F1_cluster_node_forest: ";
				F1_cluster_node_forest->print_components();
			}
			if (F2_cluster_node_forest != NULL) {
				cout << "F2_cluster_node_forest: ";
				F2_cluster_node_forest->print_components();
			}
		#endif
		return adjustment;
	}

};

void cluster_reduction_find_components(Node *n,
				vector<bool> *F2_cluster_copy_components,
				vector<bool> *old_F2_keep_components,
				int cluster_component_number) {
	Node *lc = n->lchild();
	Node *rc = n->rchild();
	if (lc != NULL)
		cluster_reduction_find_components(lc, F2_cluster_copy_components,
				old_F2_keep_components, cluster_component_number);
	if (rc != NULL)
		cluster_reduction_find_components(rc, F2_cluster_copy_components,
				old_F2_keep_components, cluster_component_number);
	if (lc == NULL && rc == NULL) {
		int cnumber = n->get_twin()->get_component_number();
		if (cnumber != cluster_component_number) {
			(*F2_cluster_copy_components)[cnumber] = true;
			(*old_F2_keep_components)[cnumber] = false;
		}
	}
}

list<ClusterInstance> cluster_reduction(Forest *old_F1, Forest *old_F2,
		list<Node *> *cluster_points) {
	list<ClusterInstance> clusters = list<ClusterInstance>();
	vector<bool> old_F2_keep_components =
		vector<bool>(old_F2->num_components(), true);
	for(list<Node *>::iterator i = cluster_points->begin();
			i != cluster_points->end(); i++) {
		// Cluster F1
		Node *F1_root_node = *i;
		Node *F1_cluster_node = F1_root_node->parent();
//		Node *p = F1_cluster_node;
//		while (p->parent() != NULL)
//			p = p->parent();
//		cout << "root address=" << &(*p) << endl;

		F1_root_node->cut_parent();
		vector<Node *> cluster_F1_components = vector<Node *>();
		cluster_F1_components.push_back(F1_root_node);
		Forest *F1 = new Forest(cluster_F1_components);

		// Cluster F2
		Node *F2_root_node = F1_root_node->get_twin();
		Node *F2_cluster_node = F2_root_node->parent();
		vector<bool> F2_cluster_copy_components =
			vector<bool>(old_F2->num_components(), false);
		int cnumber = F2_root_node->get_component_number();
		bool F2_has_component_zero = false;
		bool skip_F2_cluster = false;
		if (F2_root_node->parent() != NULL) {
			F2_root_node->cut_parent();
		}
		else {
			if (old_F2_keep_components[cnumber] == false) {
				skip_F2_cluster = true;
			}
			else if (old_F2->get_component(cnumber) == F2_root_node){
				old_F2_keep_components[cnumber] = false;		
				if (cnumber == 0) {
					F2_has_component_zero = true;
				}
			}
			else {
//				if (F2_root_node->get_forest()->get_cluster()->F2_has_component_zero == false) {
					F2_root_node->get_forest()->get_cluster()
							->F2_has_component_zero = true;
//					if (F2_root_node->get_forest()->contains_rho() == false)
//						F2_root_node->get_forest()->add_rho();
					F2_cluster_node = F2_root_node->get_forest()->get_cluster()->F2_cluster_node;
//				}
				skip_F2_cluster = true;
			}
		}

		cluster_reduction_find_components(F1_root_node,
				&F2_cluster_copy_components, &old_F2_keep_components, cnumber);
		vector<Node *> cluster_F2_components = vector<Node *>();
		if (!skip_F2_cluster)
			cluster_F2_components.push_back(F2_root_node);
		for(int i = 0; i < F2_cluster_copy_components.size(); i++) {
			if (F2_cluster_copy_components[i] == true)
				cluster_F2_components.push_back(old_F2->get_component(i));
		}
		Forest *F2 = new Forest(cluster_F2_components);
		clusters.push_back(ClusterInstance(F1, F2, F1_cluster_node,
					F2_cluster_node, F2_has_component_zero));
		F1->set_cluster(clusters.back());
		F2->set_cluster(clusters.back());
		if (F1_cluster_node != NULL)
			F1_cluster_node->increase_clustered_children();
		if (F2_cluster_node != NULL)
			F2_cluster_node->increase_clustered_children();
		F1->label_nodes_with_forest();
		F2->label_nodes_with_forest();
		F1->set_twin(F2);
		F2->set_twin(F1);

		for(int i = 0; i < F2_cluster_copy_components.size(); i++) {
			if (F2_cluster_copy_components[i] == true) {
				cluster_F2_components.push_back(
						(old_F2->get_component(i)));
//				cout << "true " ;
			}
//			else
//				cout << "false " ;
		}
//		cout << endl;
	}
	// remove any clustered components from old_F2
	vector<Node *> old_F2_remaining_components = vector<Node *>();
//	cout << "size=" << old_F2->num_components() << endl;
	for(int i = 0; i < old_F2_keep_components.size(); i++) {
		if (old_F2_keep_components[i] == true) {
			old_F2_remaining_components.push_back(
					(old_F2->get_component(i)));
//			cout << "true " ;
		}
//		else
//			cout << "false " ;
	}
//	cout << endl;
//	cout << endl;
	Forest *replace_old_F2 = new Forest(old_F2_remaining_components);
	replace_old_F2->swap(old_F2);
	replace_old_F2->erase_components();
	delete replace_old_F2;


	clusters.push_back(ClusterInstance(old_F1, old_F2, NULL, NULL, true));
	old_F1->set_cluster(clusters.back());
	old_F2->set_cluster(clusters.back());
	old_F1->label_nodes_with_forest();
	old_F2->label_nodes_with_forest();
	old_F1->set_twin(old_F2);
	old_F2->set_twin(old_F1);

	return clusters;

}


#endif
