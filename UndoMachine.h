/*******************************************************************************
UndoMachine.h

Data structure for recording and undoing tree alterations

Copyright 2012-2014 Chris Whidden
whidden@cs.dal.ca
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
#ifndef DEFINE_UNDOMACHINE

#define DEFINE_UNDOMACHINE

#include "Node.h"
#include "Forest.h"
#include <list>
#include <typeinfo>


class Node;

class Undoable {
	public:
	virtual ~Undoable() {}
	virtual void undo() = 0;
};

class UndoMachine {
	public:
	list<Undoable *> events;
	int size;

	UndoMachine() {
		events = list<Undoable *>();
		size = 0;
	}

	void add_event(Undoable *event) {
		events.push_back(event);
		size++;
	}

	void insert_event(list<Undoable *>::iterator i, Undoable *event) {
		list<Undoable *>::iterator j = i;
		j++;
		events.insert(j, event);
		size++;
	}

	 list<Undoable *>::iterator get_bookmark() {
		 list<Undoable *>::iterator bookmark = events.end();
		 if (size > 0)
		 	bookmark--;
		 else
			 bookmark = events.begin();
		 return bookmark;
	 }

	void undo() {
		if (!events.empty()) {
			Undoable *event = events.back();
#ifdef DEBUG_UNDO
			cout << typeid(*event).name() << endl;
#endif
			events.pop_back();
			event->undo();
			delete event;
			size--;
		}
	}

	void undo(int num) {
		while (num > 0)
			undo();
	}

	void undo_all() {
		undo_to(0);
	}

	void undo_to(int to) {
		while (size > to)
			undo();
	}

	void clear_to(int to) {
		while (size > to) {
			Undoable *event = events.back();
			events.pop_back();
			delete event;
			size--;
		}

	}

	int num_events() {
		return size;
	}
};

void ContractEvent(UndoMachine *um, Node *n);

class AddRho : public Undoable {
	public:
	Forest *F;
	AddRho(Forest *f) {
		F = f;
	}
	void undo() {
		F->set_rho(false);
		F->get_component(F->num_components()-1)->delete_tree();
		F->erase_components(F->num_components()-1,F->num_components());
	}
};

class AddComponent : public Undoable {
	public:
	Forest *F;
	AddComponent(Forest *f) {
		F = f;
	}
	void undo() {
		F->erase_components(F->num_components()-1,F->num_components());
	}
};

class AddComponentToFront : public Undoable {
	public:
	Forest *F;
	AddComponentToFront(Forest *f) {
		F = f;
	}
	void undo() {
		F->erase_components(0,1);
	}
};

// TODO: use new insert_child function with a stored successor sibling
// does the end work? maybe a seperate variable for that?
class CutParent : public Undoable {
	public:
	Node *child;
	Node *parent;
	int branch;
	int depth;
	CutParent(Node *c) {
		child = c;
		parent = c->parent();
		branch = 0;
		depth = c->get_depth();
		if (parent != NULL) {
			if (parent->lchild() == child)
				branch = 1;
			else
				branch = 2;
		}
	}

	void undo() {
		if (branch == 1) {
			if (parent->is_leaf())
				parent->add_child(child);
			else
				parent->insert_child(parent->get_children().front(), child);
		}
		else if (branch == 2)
			parent->add_child(child);
		child->set_depth(depth);
	}
};

class ClearSiblingPair : public Undoable {
	public:
		Node *a;
		Node *c;
		int a_status;
		int c_status;
		ClearSiblingPair(Node *x, Node *y) {
			if (x->get_sibling_pair_status() == 1 ||
					y->get_sibling_pair_status() == 2) {
				a = x;
				c = y;
			}
			else {
				a = y;
				c = x;
			}
			a_status = a->get_sibling_pair_status();
			c_status = c->get_sibling_pair_status();
		}

		void undo() {
			a->set_sibling_pair_status(1);
			c->set_sibling_pair_status(2);
			if (a_status == 0) {
				c->set_sibling(a);
			}
			else if (c_status == 0) {
				a->set_sibling(c);
			}

		}
};

class PopClearedSiblingPair : public Undoable {
	public:
		Node *a;
		Node *c;
		int a_status;
		int c_status;

		list<Node *> *sibling_pairs;
		PopClearedSiblingPair(Node *x, Node *y, list<Node *> *s) {
			sibling_pairs = s;
			a = x;
			c = y;
		}

		void undo() {
			sibling_pairs->push_back(c);
			sibling_pairs->push_back(a);
		}
};

class PopSiblingPair : public Undoable {
	public:
		Node *a;
		Node *c;
		list<Node *> *sibling_pairs;
		PopSiblingPair(Node *x, Node *y, list<Node *> *s) {
			sibling_pairs = s;
				a = x;
				c = y;
		}

		void undo() {
			sibling_pairs->push_back(c);
			sibling_pairs->push_back(a);
		}
};

class ContractSiblingPair : public Undoable {
	public:
		Node *node;
		int c1_depth;
		int c2_depth;
		bool binary_node;
		bool node_protected;
		ContractSiblingPair(Node *n) {
			init(n);
			if (n->is_protected())
				node_protected = true;
			else
				node_protected = false;
		}
		ContractSiblingPair(Node *n, Node *child1, Node *child2,
				UndoMachine *um) {
			if (n->get_children().size() == 2)
				init(n);
			else {
				um->add_event(new CutParent(child1));
				um->add_event(new CutParent(child2));
				binary_node = false;
				node = n;
			}
			if (n->is_protected())
				node_protected = true;
			else
				node_protected = false;
		}

		void init(Node *n) {
			node = n;
			if (n->lchild() != NULL)
				c1_depth = n->lchild()->get_depth();
			else
				c1_depth = -1;
			if (n->rchild() != NULL)
				c2_depth = n->rchild()->get_depth();
			else
				c2_depth = -1;
			binary_node = true;
		}

		void undo() {
			if (binary_node) {
				node->undo_contract_sibling_pair();
				if (c1_depth > -1)
					node->lchild()->set_depth(c1_depth);
				if (c2_depth > -1)
					node->rchild()->set_depth(c2_depth);
			}
			if (node_protected)
				node->protect_edge();
		}
};

class AddToFrontSiblingPairs : public Undoable {
	public:
		list<Node *> *sibling_pairs;
		AddToFrontSiblingPairs(list<Node *> *s) {
			sibling_pairs = s;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				sibling_pairs->pop_front();
				sibling_pairs->pop_front();
			}
		}
};

class AddToSiblingPairs : public Undoable {
	public:
		list<Node *> *sibling_pairs;
		AddToSiblingPairs(list<Node *> *s) {
			sibling_pairs = s;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				sibling_pairs->pop_back();
				sibling_pairs->pop_back();
			}
		}
};

class AddToSetSiblingPairs : public Undoable {
	public:
		set<SiblingPair> *sibling_pairs;
		SiblingPair pair;
		AddToSetSiblingPairs(set<SiblingPair> *sp, SiblingPair p) {
			sibling_pairs = sp;
			//pair = SiblingPair(p->a,p->c);
			pair = p;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				sibling_pairs->erase(pair);
			}
		}
};

class RemoveSetSiblingPairs : public Undoable {
	public:
		set<SiblingPair> *sibling_pairs;
		SiblingPair pair;
		RemoveSetSiblingPairs(set<SiblingPair> *sp, SiblingPair p) {
			sibling_pairs = sp;
			pair = p;
		}
		void undo() {
			sibling_pairs->insert(pair);
		}
};

class AddInSiblingPairs : public Undoable {
	public:
		list<Node *> *sibling_pairs;
		int pos;
		AddInSiblingPairs(list<Node *> *s, int p) {
			sibling_pairs = s;
			pos = p;
		}
		void undo() {
			if (!sibling_pairs->empty()) {
				list<Node *>::iterator c = sibling_pairs->begin();
				for(int i = 0; i <= pos && c != sibling_pairs->end(); i++) {
					c++;
				}
				if (c != sibling_pairs->end()) {
					list<Node *>::iterator rem = c;
					c++;
					sibling_pairs->erase(rem);
					rem = c;
					c++;
					sibling_pairs->erase(rem);
				}
			}
		}
};

class SetTwin : public Undoable {
	public:
		Node *node;
		Node *twin;

		SetTwin(Node *n) {
			node = n;
			twin = n->get_twin();
		}

		void undo() {
			node->set_twin(twin);
		}
};

class ChangeName : public Undoable {
	public:
		Node *node;
		string name;

		ChangeName(Node *n) {
			node = n;
			name = n->get_name();//str();
//			name = n->str();
		}

		void undo() {
			node->set_name(name);
		}
};

class ChangeEdgePreInterval : public Undoable {
	public:
		Node *node;
		int start;
		int end;

		ChangeEdgePreInterval(Node *n) {
			start = n->get_edge_pre_start();
			end = n->get_edge_pre_end();
			node = n;
		}

		void undo() {
			node->set_edge_pre_start(start);
			node->set_edge_pre_end(end);
		}
};

class ChangePreNum : public Undoable {
	public:
		Node *node;
		int prenum;

		ChangePreNum(Node *n) {
			prenum = n->get_preorder_number();
			node = n;
		}

		void undo() {
			node->set_preorder_number(prenum);
		}
};

class ChangeRightChild : public Undoable {
	public:
		Node *node;
		Node *rchild;
		int rchild_depth;

		ChangeRightChild(Node *n) {
			node = n;
			rchild = n->rchild();
			if (rchild != NULL)
				rchild_depth = rchild->get_depth();
		}

		void undo() {
			if (rchild != NULL) {
				//node->add_child_keep_depth(rchild);
				node->add_child(rchild);
				rchild->set_depth(rchild_depth);
			}
			else
				if (node->rchild() != NULL)
					node->rchild()->cut_parent();
		}
};

class ChangeLeftChild : public Undoable {
	public:
		Node *node;
		Node *lchild;
		int lchild_depth;

		ChangeLeftChild(Node *n) {
			node = n;
			lchild = n->lchild();
			if (lchild != NULL)
				lchild_depth = lchild->get_depth();
		}

		void undo() {
			if (lchild != NULL) {
				//node->add_child_keep_depth(lchild);
				node->add_child(lchild);
				lchild->set_depth(lchild_depth);
			}
			else
				if (node->lchild() != NULL)
					node->lchild()->cut_parent();
		}
};

class AddChild : public Undoable {
	public:
		Node *child;
		int depth;

		AddChild(Node *c) {
			child = c;
			if (c != NULL)
				depth = c->get_depth();
		}

		void undo() {
			if (child != NULL) {
				child->cut_parent();
				child->set_depth(depth);
			}
		}
};

class AddContractedLC : public Undoable {
	public:
		Node *node;

		AddContractedLC(Node *n) {
			node = n;
		}

		void undo() {
			node->set_contracted_lc(NULL);
		}
};

class AddContractedRC : public Undoable {
	public:
		Node *node;

		AddContractedRC(Node *n) {
			node = n;
		}

		void undo() {
			node->set_contracted_rc(NULL);
		}
};

class CreateNode : public Undoable {
	public:
		Node *node;

		CreateNode(Node *n) {
			node = n;
		}

		void undo() {
			if (node != NULL)
				delete node;
		}
};

class ProtectEdge : public Undoable {
	public:
		Node *node;

		ProtectEdge(Node *n) {
			node = n;
		}

		void undo() {
			if (node != NULL)
				node->unprotect_edge();
		}
};

class UnprotectEdge : public Undoable {
	public:
		Node *node;

		UnprotectEdge(Node *n) {
			node = n;
		}

		void undo() {
			if (node != NULL)
				node->protect_edge();
		}
};

class ListPushBack : public Undoable {
	public:
		list<Node *> *l;

		ListPushBack(list<Node *> *x) {
			l = x;
		}

		void undo() {
			l->pop_back();
		}
};

class ListPopBack : public Undoable {
	public:
		list<Node *> *l;
		Node *node;

		ListPopBack(list<Node *> *x) {
			l = x;
			if (!l->empty())
				node = l->back();
			else
				node = NULL;
		}

		void undo() {
			if (node != NULL)
				l->push_back(node);
		}
};


void ContractEvent(UndoMachine *um, Node *n, list<Undoable *>::iterator
		bookmark) {
		Node *parent = n->parent();
		Node *child;
		Node *lc = n->lchild();
		Node *rc = n->rchild();
		Node *ret = NULL;
		// contract out this node and give child to parent
		if (parent != NULL) {
			if (lc && !rc) {
				child = lc;
				um->add_event(new ChangeEdgePreInterval(child));
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
				if (n->is_protected() && !child->is_protected())
					um->add_event(new ProtectEdge(child));
			}
			else if (rc && !lc) {
				child = rc;
				um->add_event(new ChangeEdgePreInterval(child));
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
				if (n->is_protected() && !child->is_protected())
					um->add_event(new ProtectEdge(child));
			}
			else if (lc == NULL && rc == NULL) {
				um->insert_event(bookmark, new CutParent(n));
				parent->delete_child(n);
				ContractEvent(um, parent);
				parent->add_child(n);
			}
		}
		// if no parent then take children of single child and remove it
		else {

			// dead component or singleton, will be cleaned up by the forest
			if (n->get_children().empty()) {
				um->add_event(new ChangeName(n));
			}
			else if (n->get_children().size() == 1) {
				child = n->get_children().front();
//				if (rc == NULL) {
//					um->add_event(new ChangeRightChild(n));
//					child = lc;
//				}
//				else {
//					um->add_event(new ChangeLeftChild(n));
//					child = rc;
//				}
				um->add_event(new CutParent(child));
				/* cluster hack - if we delete a cluster node then
				 * we may try to use it later. This only happens once
				 * per cluster so we can spend linear time to update
				 * the forest
				 */
				if (child->get_num_clustered_children() > 0) {
					//um->add_event(new CutParent(n));
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
							um->add_event(new SetTwin(n));
							um->add_event(new SetTwin(child->get_twin()));
						}
						um->add_event(new ChangeName(n));
					}
					um->add_event(new ChangePreNum(n));
					//um->add_event(new CutParent(n));
					list<Node *>::iterator c;
					for(c = child->get_children().begin();
							c != child->get_children().end();
							c++) {
						um->add_event(new CutParent(*c));
					}
					if (child->get_contracted_lc() != NULL)
						um->add_event(new AddContractedLC(n));
					if (child->get_contracted_rc() != NULL)
						um->add_event(new AddContractedRC(n));
				}
			}
		}

	}

void ContractEvent(UndoMachine *um, Node *n) {
	list<Undoable *>::iterator bookmark = um->get_bookmark();
	ContractEvent(um, n, bookmark);
}

#endif
