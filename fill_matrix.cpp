// fill an upper right triangular matrix symmetrically
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>

using namespace std;

int main() {

vector<vector<int> > m = vector<vector<int> >();
string line = "";
int size = 0;
while(getline(cin, line)) {
	string token = "";
	m.push_back(vector<int>());
	for(int i = 0; i < line.size(); i++) {
		if (line[i] == ',') {
			int num = -1;
			if (token != "") {
				num = atoi(token.c_str());
			}
			m[size].push_back(num);
			token = "";
		}
		else {
			token.push_back(line[i]);
		}
	}
	int num = -1;
	if (token != "") {
		num = atoi(token.c_str());
	}
	m[size].push_back(num);
	size++;
}


for(int i = 0; i < m.size(); i++) {
	for(int j = i+1; j < m.size(); j++) {
		m[j][i] = m[i][j];
	}
}

for(int i = 0; i < m.size(); i++) {
	cout << m[i][0];
	for(int j = 1; j < m.size(); j++) {
		cout << ",";
		if (m[i][j] >= 0) {
			cout << m[i][j];
		}
	}
	cout << endl;
}


}
