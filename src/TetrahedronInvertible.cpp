#include "TetrahedronInvertible.h"
#include "Node.h"

using namespace std;
using namespace Eigen;

#define Fthreshold 0.1

TetrahedronInvertible::TetrahedronInvertible(double young, double poisson, double density, Material material, const vector<shared_ptr<Node>> &nodes):
Tetrahedron(young, poisson, density, material, nodes), m_isInvertible(true)
{

}

void TetrahedronInvertible::precompute() {
	// Compute volume
	this->W = abs(1.0 / 6.0 * Dm.determinant());
	m_mass = (this->W * this->m_density);

	// Distribute 1/4 mass to each node
	// Initialize the undeformed position vector
	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->m += this->m_mass * 0.25;
	}

	m_delta_L = 0.001;
	m_delta_U = 2.0;

}

//Matrix3x4d TetrahedronInvertible::computeAreaWeightedVertexNormals() {
//	Vector3d va, vb, vc, vd;
//	va = m_nodes[0]->x;
//	vb = m_nodes[1]->x;
//	vc = m_nodes[2]->x;
//	vd = m_nodes[3]->x;
//
//	// Computes normals for the four faces: acb, adc, abd, bcd
//	Vector3d acb_normal, adc_normal, abd_normal, bcd_normal;
//	acb_normal = (vc - va).cross(vb - va);
//	adc_normal = (vd - va).cross(vc - va);
//	abd_normal = (vb - va).cross(vd - va);
//	bcd_normal = (vc - vb).cross(vd - vb);
//
//	// if the tet vertices abcd form a positive orientation, no need to correct
//	// if not, flip them
//	double orientation = (vd - va).dot((vb - va).cross(vc - va));
//	if (orientation < 0.0) {
//		acb_normal *= -1.0;
//		adc_normal *= -1.0;
//		abd_normal *= -1.0;
//		bcd_normal *= -1.0;
//	}
//
//	// Computes the area of triangles
//	// area = 0.5 * | u x v |
//	double acb_area, adc_area, abd_area, bcd_area;
//	acb_area = 0.5 * sqrt(acb_normal.dot(acb_normal));
//	adc_area = 0.5 * sqrt(adc_normal.dot(adc_normal));
//	abd_area = 0.5 * sqrt(abd_normal.dot(abd_normal));
//	bcd_area = 0.5 * sqrt(bcd_normal.dot(bcd_normal));
//
//	acb_normal.normalize();
//	adc_normal.normalize();
//	abd_normal.normalize();
//	bcd_normal.normalize();
//
//	//this->Nm.col(0) = -(acb_area * acb_normal + adc_area * adc_normal + abd_area * abd_normal) / 3.0;
//	//this->Nm.col(1) = -(acb_area * acb_normal + abd_area * abd_normal + bcd_area * bcd_normal) / 3.0;
//	//this->Nm.col(2) = -(acb_area * acb_normal + adc_area * adc_normal + bcd_area * bcd_normal) / 3.0;
//	//this->Nm.col(3) = -(adc_area * adc_normal + abd_area * abd_normal + bcd_area * bcd_normal) / 3.0;
//
//	return this->Nm;
//}

bool TetrahedronInvertible::checkNecessityForSVD(double deltaL, double deltaU, Matrix3d F) {
	//Matrix3d C;
	//C.noalias() = F.transpose() * F;
	// C = |c1 c2 c3|
	//     |c2 c4 c5|
	//     |c3 c5 c6|

	// f(lambda) = lambda^3 + a * lambda^2 + b * lambda + c
	// a = -(c1 + c4 + c6)
	// b = (c1*c4 + c4*c6 + c1*c6 - c2^2 - c3^2 - c5^2)
	// c = (c1*c5^2 + c4*c3^2 +c6*c2^2 - c1*c4*c6-2*c2*c3*c5) 
	// 
	// alpha, beta are the solutions of f'(lambda) = 0 (alpha < beta)
	// f'(lambda) = 3*lambda^2 + 2a*lambda + b = 0

	double c1, c2, c3, c4, c5, c6;


	c1 = F(0, 0) * F(0, 0) + F(1, 0) * F(1, 0) + F(2, 0) * F(2, 0);
	c2 = F(0, 0) * F(0, 1) + F(1, 0) * F(1, 1) + F(2, 0) * F(2, 1);
	c3 = F(0, 0) * F(0, 2) + F(1, 0) * F(1, 2) + F(2, 0) * F(2, 2);
	c4 = F(0, 1) * F(0, 1) + F(1, 1) * F(1, 1) + F(2, 1) * F(2, 1);
	c5 = F(0, 1) * F(0, 2) + F(1, 1) * F(1, 2) + F(2, 1) * F(2, 2);
	c6 = F(0, 2) * F(0, 2) + F(1, 2) * F(1, 2) + F(2, 2) * F(2, 2);

	double a, b, c;
	a = -(c1 + c4 + c6);

	double c14, c46, c16, c22, c33, c55;
	c14 = c1 * c4;
	c46 = c4 * c6;
	c16 = c1 * c6;

	c22 = c2 * c2;
	c33 = c3 * c3;
	c55 = c5 * c5;

	b = (c14 + c46 + c16 - c22 - c33 - c55);
	c = c1 * c55 + c4 * c33 + c6 * c22 - c14 * c6 - 2 * c2 * c3 * c5;//13

	double alpha, beta, same;

	same = a * a - 3.0 * b;
	if (same >= -0.00000001) {
		same = abs(same);
		// f'(lambda) = 0 has one or two solutions
		same = sqrt(same);
		alpha = (-a - same) / 3.0;
		beta = (-a + same) / 3.0;
		// If the element has a small deformation, it should satisfy these four conditions:
		// f(deltaL) < 0, 
		// f(deltaU) > 0, 
		// deltaL < alpha, 
		// beta < deltaU

		double f_deltaL = deltaL*deltaL*(deltaL + a) + b * deltaL + c;
		double f_deltaU = deltaU*deltaU*(deltaU + a) + b * deltaU + c;

		if ((f_deltaL < 0) && (f_deltaU > 0) && (deltaL < alpha) && (beta < deltaU)) {
			// This element has a small deformation, no need to do SVD
			m_isSVD = false;
			//cout << "no no " << endl;
		}
		else {
			// large deformation, need to do SVD
			m_isSVD = true;
		}
	}
	else {
		// no solutions
		m_isSVD = true;
		//cout << "Do svd!no sol " << endl;
	}

	return m_isSVD;
}


VectorXd TetrahedronInvertible::computeElasticForces(VectorXd f) {
	bool print = false;

	this->F = computeDeformationGradient();
	this->FTF.noalias() = F.transpose() * F;

	// The deformation gradient is available in this->F
	if (this->F.determinant() <= 0.0) {
		m_isInverted = true;
		m_isSVD = true;
	}
	else {
		// We need an additional step to decide whether the element undergoes large deformation
		m_isInverted = false;
		//checkNecessityForSVD(m_delta_L, m_delta_U, this->F);
		
		// The necessity of SVD is available in m_isSVD
	}
	m_isSVD = true;
	if (m_isSVD) {
		// SVD on the deformation gradient
		int modifiedSVD = 1;

		if (!SVD(this->F, this->U, this->Fhats, this->V, 1e-8, modifiedSVD)) {
			//cout << "error in svd " << endl;
		}
		this->UT = this->U.transpose();
		this->VT = this->V.transpose();
		this->Fhat = Fhats.asDiagonal();

		// SVD result is available in this->U, this->V, Fhat_vec, this->Fhat
		// clamp if below the principal stretch threshold
		clamped = 0;
		for (int i = 0; i < 3; i++)
		{
			if (this->Fhats(i) < Fthreshold)
			{
				this->Fhat(i, i) = Fthreshold;
				clamped |= (1 << i);
			}
		}

		//clamped = 0; // disable clamping

		// Computes the internal forces
		// Computes P first and computes the nodal forces G=PBm in section 4 of [Irving 04]

		// Computes the diagonal P tensor
		this->Phat = computePKStress(this->Fhat, m_mu, m_lambda);

		// P = U * diag(Phat) * V'
		this->P.noalias() = this->U * this->Phat * this->VT;
	}
	else {
		// Don't do SVD
		this->P = computePKStress(this->F, m_mu, m_lambda);
	}

	this->H.noalias() = -W * this->P * this->BmT;

	// Computes the nodal forces by G=PBm=PNm
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		int rowi = m_nodes[i]->idxM;
		f.segment<3>(rowi) += this->H.col(i);
		int row3 = m_nodes[3]->idxM;
		f.segment<3>(row3) -= this->H.col(i);

	}
	return f;
}


void TetrahedronInvertible::computeForceDifferentials(VectorXd dx, VectorXd &df) {

	this->dF = computeDeformationGradientDifferential(dx);
	Matrix3d UTdFV;
	UTdFV.noalias() = this->UT * this->dF * this->V;

	this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);
	this->dP.noalias() = this->U * this->dPhat * this->VT;
	this->dH.noalias() = -W * dP * this->BmT;

	//modify hessian to compute correct values if in the inversion handling regime
	//hessian = clampHessian(hessian, clamped);

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		df.segment<3>(3 * m_nodes[i]->i) += dH.col(i);
		df.segment<3>(3 * m_nodes[3]->i) -= dH.col(i);
	}
}

void TetrahedronInvertible::computeForceDifferentials(VectorXd dx, int row, int col, MatrixXd &K) {

	this->dF = computeDeformationGradientDifferential(dx);
	Matrix3d UTdFV;
	UTdFV.noalias() = this->UT * this->dF * this->V;
	this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);

	this->dP.noalias() = this->U * this->dPhat * this->VT;
	this->dH.noalias() = -W * dP * this->BmT;

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		Vector3d temp = dH.col(i);
		for (int j = 0; j < 3; ++j) {
			K(row + 3 * m_nodes[i]->i + j, col) += temp(j);
			K(row + 3 * m_nodes[3]->i + j, col) += -temp(j);
		}
	}
}

void TetrahedronInvertible::computeForceDifferentials(MatrixXd &K_global) {
	this->K.setZero();
	Matrix3d UTdFV;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			this->dDs.setZero();
			this->dDs(j, i) = 1.0;
			this->dF.noalias() = this->dDs * this->Bm;
			UTdFV.noalias() = this->UT * this->dF * this->V;
			this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);

			this->dP.noalias() = this->U * this->dPhat * this->VT;
			this->dH.noalias() = -W * dP * this->BmT;

			//Matrix3d hessian = this->dP * this->Nm.block<3, 3>(0, 0);
			for (int t = 0; t < 3; t++) {
				Vector3d temp = dH.col(t);
				this->K.block<3, 1>(t * 3, i * 3 + j) = temp;
				this->K.block<3, 1>(9, i * 3 + j) -= temp;
				this->K.block<1, 3>(i * 3 + j, 9) -= temp;
			}
		}
	}

	this->K.block<3, 3>(9, 9) = -this->K.block<3, 3>(0, 9) - this->K.block<3, 3>(3, 9) - this->K.block<3, 3>(6, 9);
	//cout << "Ki" << endl << this->K << endl;
	//VectorXcd ev = this->K.eigenvalues();
	//cout << ev << endl; 

	int i = m_nodes[0]->idxM;
	int j = m_nodes[1]->idxM;
	int k = m_nodes[2]->idxM;
	int l = m_nodes[3]->idxM;

	K_global.block<3, 3>(i, i) += this->K.block<3, 3>(0, 0);
	K_global.block<3, 3>(i, j) += this->K.block<3, 3>(0, 3);
	K_global.block<3, 3>(i, k) += this->K.block<3, 3>(0, 6);
	K_global.block<3, 3>(i, l) += this->K.block<3, 3>(0, 9);

	K_global.block<3, 3>(j, i) += this->K.block<3, 3>(3, 0);
	K_global.block<3, 3>(j, j) += this->K.block<3, 3>(3, 3);
	K_global.block<3, 3>(j, k) += this->K.block<3, 3>(3, 6);
	K_global.block<3, 3>(j, l) += this->K.block<3, 3>(3, 9);

	K_global.block<3, 3>(k, i) += this->K.block<3, 3>(6, 0);
	K_global.block<3, 3>(k, j) += this->K.block<3, 3>(6, 3);
	K_global.block<3, 3>(k, k) += this->K.block<3, 3>(6, 6);
	K_global.block<3, 3>(k, l) += this->K.block<3, 3>(6, 9);

	K_global.block<3, 3>(l, i) += this->K.block<3, 3>(9, 0);
	K_global.block<3, 3>(l, j) += this->K.block<3, 3>(9, 3);
	K_global.block<3, 3>(l, k) += this->K.block<3, 3>(9, 6);
	K_global.block<3, 3>(l, l) += this->K.block<3, 3>(9, 9);
}

void TetrahedronInvertible::computeForceDifferentialsSparse(VectorXd dx, int row, int col, vector<T> &K_) {

	this->dF = computeDeformationGradientDifferential(dx);
	if (m_isSVD) {
		Matrix3d UTdFV = this->UT * this->dF * this->V;
		this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);
		this->dP = this->U * this->dPhat * this->VT;
	}
	else
	{
		this->dP = computePKStressDerivative(this->F, this->dF, m_mu, m_lambda);
	}
	this->dH.noalias() = -W * dP * this->BmT;
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {

		for (int j = 0; j < 3; ++j) {
			K_.push_back(T(row + 3 * m_nodes[i]->i + j, col, dH(j, i)));
			K_.push_back(T(row + 3 * m_nodes[3]->i + j, col, -dH(j, i)));
		}
	}
}

Matrix3d TetrahedronInvertible::clampHessian(Matrix3d &hessian, int clamped) {
	if (clamped & 1) // first lambda was clamped (in inversion handling)
	{
		hessian(0, 0) = 0.0;
		hessian(0, 1) = 0.0;
		hessian(1, 0) = 0.0;
		hessian(2, 0) = 0.0;
		hessian(0, 2) = 0.0;
		//cout << "clamped 1" << endl;
	}

	if (clamped & 2) // second lambda was clamped (in inversion handling)
	{
		hessian(0, 1) = 0.0;
		hessian(1, 0) = 0.0;
		hessian(1, 1) = 0.0;
		hessian(1, 2) = 0.0;
		hessian(2, 1) = 0.0;
		//cout << "clamped 2" << endl;
	}

	if (clamped & 4) // third lambda was clamped (in inversion handling)
	{
		Matrix3d new_hessian = hessian;
		hessian.setZero();
		hessian(1, 1) = new_hessian(1, 1);

		/*hessian(0, 0) = 0.0;
		hessian(0, 1) = 0.0;
		hessian(1, 0) = 0.0;
		hessian(2, 0) = 0.0;
		hessian(0, 2) = 0.0;
		hessian(1, 2) = 0.0;
		hessian(2, 1) = 0.0;
		hessian(2, 2) = 0.0;*/
		//cout << "clamped 3" << endl;
	}
	return hessian;
}