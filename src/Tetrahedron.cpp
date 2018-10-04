#include "Tetrahedron.h"
#include "Node.h"

using namespace Eigen;
using namespace std;


Tetrahedron::Tetrahedron()
{

}

Tetrahedron::Tetrahedron(double young, double poisson, double density, const vector<shared_ptr<Node>> &nodes):
m_young(young), m_poisson(poisson), m_density(density), m_nodes(nodes)
{



}

Tetrahedron:: ~Tetrahedron() {

}