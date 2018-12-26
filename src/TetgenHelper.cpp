#include "TetgenHelper.h"

#include "Node.h"
#include <fstream>
#include <iostream>

using namespace Eigen;
using namespace std;

void TetgenHelper::createNodeFile(vector<shared_ptr<Node>> &nodes, char *filename) {
	ofstream fout(filename);
	//fout << "#node" << endl;
	fout << (int)nodes.size() << "  " << "3" << "  " << "0" << "  " << "0" <<  endl;
	for (int i = 0; i < (int)nodes.size(); ++i) {
		auto node = nodes[i];
		fout << node->i << "  " << node->x.x() << "  " << node->x.y() << "  " << node->x.z() << endl;
	}
	fout.close();
}