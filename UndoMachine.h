#ifndef DEFINE_UNDOMACHINE

#define DEFINE_UNDOMACHINE

#include "Node.h"
#include "Forest.h"
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

// TODO: this needs to be fixed with the removal of set_parent
// TODO: this needs to be fixed with the removal of set_xchild
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
			else if (parent->rchild() == child)
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
		ContractSiblingPair(Node *n) {
			init(n);
		}
		ContractSiblingPair(Node *n, Node *child1, Node *child2,
				UndoMachine *um) {
			if (n->get_children().size() == 2)
				init(n);
			else {
				um->add_event(new CutParent(child1));
				um->add_event(new CutParent(child2));
				binary_node = false;
			}
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
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
			}
			else if (rc && !lc) {
				child = rc;
				um->add_event(new CutParent(child));
				um->add_event(new CutParent(n));
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

void ContractEvent(UndoMachine *um, Node *n) {
	list<Undoable *>::iterator bookmark = um->get_bookmark();
	ContractEvent(um, n, bookmark);
}

#endif
