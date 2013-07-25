
/*******************************************************************************
sparse_counts.h

Data structure for a sparse matrix of counts

Copyright 2013 Chris Whidden
cwhidden@dal.ca
http://kiwi.cs.dal.ca/Software/RSPR
July 19, 2013
Version 1.2

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
#ifndef INCLUDE_SPARSE_COUNTS

#define INCLUDE_SPARSE_COUNTS
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
//#include <list>
//#include <deque>
//#include "Node.h"
//#include "LCA.h"
#include <map>
#include <limits>

using namespace std;

// TODO: doesn't allow insertion with larger x than initially allocated

void swap_if_larger(int *x, int *y);
	vector<pair<int, int>> find_most_common_pairs_scaled(vector<double> *scale);
int count_intersection(vector<int> *A, vector<int> *B);

template <class T>
class SparseCounts {
	private:
	vector<map<int, T> > sparse_counts; 

	public:

	SparseCounts(int x, int y) {
		sparse_counts = vector<map<int, T> >(x);
		for (int i = 0; i < y; i++) {
			sparse_counts[i] = map<int, T> ();
		}
	}

	void increment(int x, int y) {
		swap_if_larger(&x, &y);
		T val = sparse_get(x, y);
		sparse_set(x, y, val+1);
	}

	T sparse_get(int x, int y) {
		swap_if_larger(&x, &y);
		typename map<int, T>::iterator it;
		it = sparse_counts[x].find(y);
		if (it == sparse_counts[x].end()) {
			return 0;
		}
		else {
			return it->second;
		}
	}

	void sparse_set(int x, int y, int val) {
		swap_if_larger(&x, &y);
		//sparse_counts[x][y] = make_pair<int, T>(y, val);
		sparse_counts[x][y] = val;
	}

	void sparse_print() {
		int end1 = sparse_counts.size();
		for(int i = 0; i < end1; i++) {
			typename map<int, T>::iterator it;
			
			for(it = sparse_counts[i].begin(); it != sparse_counts[i].end();
					it++) {
				cout << i << ", " << it->first << ", " << it->second << endl;
			}
		}
	}

	void sparse_labelled_print(map<int, string> *reverse_label_map) {
		int end1 = sparse_counts.size();
		for(int i = 0; i < end1; i++) {
			typename map<int, T>::iterator it;
			
			for(it = sparse_counts[i].begin(); it != sparse_counts[i].end();
					it++) {
				cout << reverse_label_map->find(i)->second << ", "
						<< reverse_label_map->find(it->first)->second
						<< ", " << it->second << endl;
			}
		}
	}

	void clear() {
		for(int i = 0; i < sparse_counts.size(); i++) {
			sparse_counts[i].clear();
		}
	}

	// return a vector of the most common pairs
	vector<pair<int, int>> find_most_common_pairs() {
		int end1 = sparse_counts.size();
		vector<pair<int, int>> mcp = vector<pair<int, int> >();
		int max_count = 0;
		typename map<int, T>::iterator it;
		for(int i = 0; i < end1; i++) {
			for(it = sparse_counts[i].begin(); it != sparse_counts[i].end();
					it++) {
				if (it->second > max_count) {
					max_count = it->second;
					mcp.clear();
					mcp.push_back(make_pair(i, it->first));
				}
				else if (it->second == max_count) {
					mcp.push_back(make_pair(i, it->first));
				}
			}
		}
		return mcp;
	}

	vector<pair<int, int>> find_most_common_pairs_scaled(vector<double> *scale) {
		int end1 = sparse_counts.size();
		vector<pair<int, int>> mcp = vector<pair<int, int> >();
		double max_count = 0;
		typename map<int, T>::iterator it;
		for(int i = 0; i < end1; i++) {
			for(it = sparse_counts[i].begin(); it != sparse_counts[i].end();
					it++) {
				double score = it->second * it->second;
				score /= (*scale)[it->first];
				score /= (*scale)[i];
				score = sqrt(score);
				if (score > max_count) {
					max_count = score;
					mcp.clear();
					mcp.push_back(make_pair(i, it->first));
				}
				else if (score == max_count) {
					mcp.push_back(make_pair(i, it->first));
				}
			}
		}
		return mcp;
	}

	vector<pair<int, int>> find_most_common_pairs_scaled(vector<double> *scale, vector<vector<int> > *component_trees) {
		int end1 = sparse_counts.size();
		vector<pair<int, int>> mcp = vector<pair<int, int> >();
		double max_count = 0;
		typename map<int, T>::iterator it;
		int end2 = component_trees->size();
		int max_trees = 0;
		for(int i = 0; i < end2; i++) {
			int count = (*component_trees)[i].size();
			if (count > max_trees) {
				max_trees = count;
			}
		}
//		cout << "max: " << max_trees << endl;
		for(int i = 0; i < end1; i++) {
			for(it = sparse_counts[i].begin(); it != sparse_counts[i].end();
					it++) {
				double score = it->second * it->second;
				int common_trees = count_intersection(&(*component_trees)[i], &(*component_trees)[it->first]);
				score /= (double)common_trees;
				score /= (double)common_trees;
				score /= (*scale)[it->first];
				score /= (*scale)[i];
				score = sqrt(score);
//				cout << i << ", " << it->first << ": " << common_trees << endl;
				if (score > max_count) {
					max_count = score;
					mcp.clear();
					mcp.push_back(make_pair(i, it->first));
				}
				else if (score == max_count) {
					mcp.push_back(make_pair(i, it->first));
				}
			}
		}
		return mcp;
	}

};

void swap_if_larger(int *x, int *y) {
	if (*x > *y) {
		int temp = *x;
		*x = *y;
		*y = temp;
	}
}

// count the number of common integers from two sorted vectors
int count_intersection(vector<int> *A, vector<int> *B) {
	vector<int>::iterator a = A->begin();
	vector<int>::iterator b = B->begin();
	int count;
	while(a != A->end() && b != B->end()) {
		if (*a < *b) {
			a++;
		}
		else if (*b > *a) {
			b++;
		}
		else {
			a++;
			b++;
			count++;
		}
	}
	return count;
}

// TODO: functions for union and intersection of sorted int vectors

#endif
