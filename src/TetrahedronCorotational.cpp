#include "TetrahedronCorotational.h"
#include "Node.h"

using namespace std;
using namespace Eigen;

TetrahedronCorotational::TetrahedronCorotational(double young, double poisson, double density, Material material, const vector<shared_ptr<Node>> &nodes) :
	Tetrahedron(young, poisson, density, material, nodes)
{

}

void TetrahedronCorotational::precompute() {
	// Compute volume
	this->W = abs(1.0 / 6.0 * Dm.determinant());
	m_mass = (this->W * this->m_density);

	// Distribute 1/4 mass to each node
	// Initialize the undeformed position vector
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->m += this->m_mass * 0.25;
		this->xx.segment<3>(3 * i) = m_nodes[i]->x0;
	}

	// Compute stiffness matrix
	Matrix6d E = hooke(m_young, m_poisson);
	Matrix3d E1, E2;
	E1 = E.block<3, 3>(0, 0);
	E2 = E.block<3, 3>(3, 3);
	Ke.setZero();

	MatrixXd ys(4, 3); // Four auxiliary vectors y0, y1, y2, y3 
	ys.block<3, 3>(0, 0) = this->Bm;
	ys.row(3) = -ys.row(0) - ys.row(1) - ys.row(2);

	Matrix3d Ni, Nj, NiEn, Si, Sj, SiEs, Kij; // Diagonal Matrix
	Vector3d yi, yj;

	for (int i = 0; i < 4; i++) {

		// Compute vertex i 
		int ii = 3 * i; // Local starting index into the 12x12 Ke matrix
		yi = ys.row(i);
		Ni.setZero();
		for (int idx = 0; idx < 3; idx++) {
			Ni(idx, idx) = yi(idx);
		}

		NiEn = Ni * E1;

		Si << yi(1), 0, yi(2),
			yi(0), yi(2), 0,
			0, yi(1), yi(0);

		SiEs = Si * E2;

		for (int j = i; j < 4; j++) {
			//Compute vertex j
			int jj = 3 * j; // Local starting index into the 12x12 Ke matrix
			yj = ys.row(j);
			Nj.setZero();
			for (int idx = 0; idx < 3; idx++) {
				Nj(idx, idx) = yj(idx);
			}

			Sj << yj(1), 0, yj(2),
				yj(0), yj(2), 0,
				0, yj(1), yj(0);

			Kij = NiEn * Nj + SiEs * (Sj.transpose());

			Ke.block<3, 3>(ii, jj) -= Kij;
			if (i != j) {
				Ke.block<3, 3>(jj, ii) -= Kij.transpose();
			}
		}
	}

	this->Kexx = this->Ke * this->xx;
}

VectorXd TetrahedronCorotational::computeElasticForces(VectorXd f) {
	this->F = computeDeformationGradient();
	this->R = gs3(this->F);

	this->Re.setZero();
	for (int i = 0; i < 4; ++i) {
		this->Re.block<3, 3>(3 * i, 3 * i) = this->R;
	}
	this->RKR.noalias() = this->Re * this->Ke * (this->Re.transpose());

	// Forces
	this->RKxx = this->Re * this->Kexx;

	for (int i = 0; i < 4; ++i) {
		this->x_vec.segment<3>(3 * i) = m_nodes[i]->x;
	}

	this->Kx = this->RKR * this->x_vec;

	for (int i = 0; i < 4; i++) {
		int rowi = m_nodes[i]->idxM;
		f.segment<3>(rowi) = f.segment<3>(rowi) - this->RKxx.segment<3>(3 * i) + this->Kx.segment<3>(3 * i);
	}
	return f;
}

void TetrahedronCorotational::computeForceDifferentials(MatrixXd &K_global) {
	int i = m_nodes[0]->idxM;
	int j = m_nodes[1]->idxM;
	int k = m_nodes[2]->idxM;
	int l = m_nodes[3]->idxM;

	K_global.block<3, 3>(i, i) += this->RKR.block<3, 3>(0, 0);
	K_global.block<3, 3>(i, j) += this->RKR.block<3, 3>(0, 3);
	K_global.block<3, 3>(i, k) += this->RKR.block<3, 3>(0, 6);
	K_global.block<3, 3>(i, l) += this->RKR.block<3, 3>(0, 9);

	K_global.block<3, 3>(j, i) += this->RKR.block<3, 3>(3, 0);
	K_global.block<3, 3>(j, j) += this->RKR.block<3, 3>(3, 3);
	K_global.block<3, 3>(j, k) += this->RKR.block<3, 3>(3, 6);
	K_global.block<3, 3>(j, l) += this->RKR.block<3, 3>(3, 9);

	K_global.block<3, 3>(k, i) += this->RKR.block<3, 3>(6, 0);
	K_global.block<3, 3>(k, j) += this->RKR.block<3, 3>(6, 3);
	K_global.block<3, 3>(k, k) += this->RKR.block<3, 3>(6, 6);
	K_global.block<3, 3>(k, l) += this->RKR.block<3, 3>(6, 9);

	K_global.block<3, 3>(l, i) += this->RKR.block<3, 3>(9, 0);
	K_global.block<3, 3>(l, j) += this->RKR.block<3, 3>(9, 3);
	K_global.block<3, 3>(l, k) += this->RKR.block<3, 3>(9, 6);
	K_global.block<3, 3>(l, l) += this->RKR.block<3, 3>(9, 9);
}

void TetrahedronCorotational::computeForceDifferentialsSparse(vector<T> &K_) {
	int a = m_nodes[0]->idxM;
	int b = m_nodes[1]->idxM;
	int c = m_nodes[2]->idxM;
	int d = m_nodes[3]->idxM;
	for (int idx = 0; idx < 3; idx++) {
		for (int jdx = 0; jdx < 3; jdx++) {
			K_.push_back(T(a + idx, a + jdx, RKR(0 + idx, 0 + jdx)));
			K_.push_back(T(a + idx, b + jdx, RKR(0 + idx, 3 + jdx)));
			K_.push_back(T(a + idx, c + jdx, RKR(0 + idx, 6 + jdx)));
			K_.push_back(T(a + idx, d + jdx, RKR(0 + idx, 9 + jdx)));

			K_.push_back(T(b + idx, a + jdx, RKR(3 + idx, 0 + jdx)));
			K_.push_back(T(b + idx, b + jdx, RKR(3 + idx, 3 + jdx)));
			K_.push_back(T(b + idx, c + jdx, RKR(3 + idx, 6 + jdx)));
			K_.push_back(T(b + idx, d + jdx, RKR(3 + idx, 9 + jdx)));

			K_.push_back(T(c + idx, a + jdx, RKR(6 + idx, 0 + jdx)));
			K_.push_back(T(c + idx, b + jdx, RKR(6 + idx, 3 + jdx)));
			K_.push_back(T(c + idx, c + jdx, RKR(6 + idx, 6 + jdx)));
			K_.push_back(T(c + idx, d + jdx, RKR(6 + idx, 9 + jdx)));

			K_.push_back(T(d + idx, a + jdx, RKR(9 + idx, 0 + jdx)));
			K_.push_back(T(d + idx, b + jdx, RKR(9 + idx, 3 + jdx)));
			K_.push_back(T(d + idx, c + jdx, RKR(9 + idx, 6 + jdx)));
			K_.push_back(T(d + idx, d + jdx, RKR(9 + idx, 9 + jdx)));
		}
	}
}
