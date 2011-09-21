#ifndef DEFINE_UNDOMACHINE

#define DEFINE_UNDOMACHINE

#include "Node.h"
#include "Forest.h"

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

	void undo() {
		if (!events.empty()) {
			Undoable *event = events.back();
			//cout << typeid(*event).name() << endl;
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
			else if (parent->rchild() == child)
				branch = 2;
		}
	}

	void undo() {
		if (branch > 0)
		child->set_parent(parent);
		child->set_depth(depth);
		if (branch == 1)
			parent->set_lchild_keep_depth(child);
		else if (branch == 2)
			parent->set_rchild_keep_depth(child);
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
		ContractSiblingPair(Node *n) {
			node = n;
		}

		void undo() {
			node->undo_contract_sibling_pair();
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
			name = n->str();
		}

		void undo() {
			node->set_name(name);
		}
};

class ChangeRightChild : public Undoable {
	public:
		Node *node;
		Node *rchild;

		ChangeRightChild(Node *n) {
			node = n;
			rchild = n->rchild();
		}

		void undo() {
			node->set_rchild_keep_depth(rchild);
		}
};

class ChangeLeftChild : public Undoable {
	public:
		Node *node;
		Node *lchild;

		ChangeLeftChild(Node *n) {
			node = n;
			lchild = n->lchild();
		}

		void undo() {
			node->set_lchild_keep_depth(lchild);
		}
};



void ContractEvent(UndoMachine *um, Node *n) {
		Node *parent = n->parent();
		Node *child;
		Node *lc = n->lchild();
		Node *rc = n->rchild();
		Node *ret = NULL;
		// contract out this node and give child to parent
		if (parent != NULL) {
			if (lc && !rc) {
				child = lc;
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
			}
			else if (rc && !lc) {
				child = rc;
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
			}
			else if (lc == NULL && rc == NULL) {
				um->add_event(new CutParent(n));
				parent->delete_child(n);
				ContractEvent(um, parent);
				parent->add_child(n);
			}
		}
		// if no parent then take children of single child and remove it
		else {

			// dead component or singleton, will be cleaned up by the forest
			if (lc == NULL && rc == NULL) {
				um->add_event(new ChangeName(n));
			}
			else if ((bool)lc xor (bool)rc) {
				child = lc;
				if (rc == NULL) {
					um->add_event(new ChangeRightChild(n));
					child = lc;
				}
				else {
					um->add_event(new ChangeLeftChild(n));
					child = rc;
				}
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
						um->add_event(new SetTwin(n));
						um->add_event(new SetTwin(child->get_twin()));
						um->add_event(new ChangeName(n));
					}
					//um->add_event(new CutParent(n));
					if (new_lc != NULL) {
						um->add_event(new CutParent(new_lc));
					}
					if (new_rc != NULL)
						um->add_event(new CutParent(new_rc));
				}
			}
		}

	}
#endif
