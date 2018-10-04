#pragma once
#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Node;

class Tetrahedron
{
public:
	Tetrahedron();
	Tetrahedron(double young, double poisson, double density, const std::vector<std::shared_ptr<Node>> &nodes);

	virtual ~Tetrahedron();

	std::vector<std::shared_ptr<Node>> m_nodes;	// i, j, k, l

private:

	double m_young;
	double m_poisson;
	// Lame coefficients
	double mu;
	double lambda;
	double mass;
	double m_density;
	

	
};