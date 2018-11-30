#include "Tetrahedron.h"
#include "Node.h"
#include <iostream>
#include "svd3.h"
#include <cmath>

#include "MatrixStack.h"
#include "Program.h"

using namespace Eigen;
using namespace std;

#define Fthreshold 0.1
// 0.45 corotated
// 0.65 neo


Tetrahedron::Tetrahedron()
{

}

Tetrahedron::Tetrahedron(double young, double poisson, double density, Material material, const vector<shared_ptr<Node>> &nodes, const Vector108d &dDsdU_vec) :
	m_young(young), m_poisson(poisson), m_density(density), m_material(material), m_nodes(nodes)
{
	m_mu = m_young / (2.0 * (1.0 + m_poisson));
	m_lambda = m_young * m_poisson / ((1.0 + m_poisson) * (1.0 - 2.0 * m_poisson));

	// Compute the inverse of the Dm matrix
	// Dm is a 3x3 matrix where the columns are the edge vectors of a tet in rest configuration

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->Dm.col(i) = m_nodes[i]->x0 - m_nodes[3]->x0;
	}
	this->Bm = this->Dm.inverse();
	
	// Compute volume
	this->W = abs(1.0 / 6.0 * Dm.determinant());
	m_mass = this->W * this->m_density;

	// Distribute 1/4 mass to each node
	for (int i = 0; i < (int)nodes.size(); i++) {
		nodes[i]->m += this->m_mass * 0.25;
	}

	computeAreaWeightedVertexNormals();

	// Compute dFdU = d F / d U
	// constant
	compute_dFdU(dDsdU_vec);

	// set the renumbering indices for conversion from Teran's order to row-major order
	rowMajorMatrixToTeran[0] = 0;
	rowMajorMatrixToTeran[1] = 3;
	rowMajorMatrixToTeran[2] = 5;
	rowMajorMatrixToTeran[3] = 4;
	rowMajorMatrixToTeran[4] = 1;
	rowMajorMatrixToTeran[5] = 7;
	rowMajorMatrixToTeran[6] = 6;
	rowMajorMatrixToTeran[7] = 8;
	rowMajorMatrixToTeran[8] = 2;

	for (int i = 0; i<9; i++)
		teranToRowMajorMatrix[rowMajorMatrixToTeran[i]] = i;


	m_delta_L = 0.25;
	m_delta_U = 5.0;

	m_delta_L = 0.0;
	m_delta_U = 1005.0;
}

void Tetrahedron::compute_dFdU(const Vector108d &dDsdU_vec) {

	for (int index = 0; index<108; index++)
	{
		int n = index % 3;
		int m = (int)(index / 3) % 4;
		int j = (int)(index / 12) % 3;
		int i = (int)(index / 36) % 3;
		double result = 0.0;
	
		for (int k = 0; k < 3; k++)
			result += dDsdU_vec[tensor9x12Index(i, k, m, n)] * this->Bm(k, j);
			//result += dDSdU[tensor9x12Index(i, k, m, n)] * dmInv[k][j];
		this->dFdU[tensor9x12Index(i, j, m, n)] = result;
	}
}

Matrix3d Tetrahedron::computeDeformationGradient() {
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	return this->F;
}

Matrix3d Tetrahedron::computeDeformationGradientDifferential(VectorXd dx) {
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->dDs.col(i) = dx.segment<3>(3 * m_nodes[i]->i) - dx.segment<3>(3 * m_nodes[3]->i);
	}
	this->dF = dDs * Bm;
	return this->dF;
}

Matrix3x4d Tetrahedron::computeAreaWeightedVertexNormals() {
	Vector3d va, vb, vc, vd;
	va = m_nodes[0]->x;
	vb = m_nodes[1]->x;
	vc = m_nodes[2]->x;
	vd = m_nodes[3]->x;

	// Computes normals for the four faces: acb, adc, abd, bcd
	Vector3d acb_normal, adc_normal, abd_normal, bcd_normal;
	acb_normal = (vc - va).cross(vb - va);
	adc_normal = (vd - va).cross(vc - va);
	abd_normal = (vb - va).cross(vd - va);
	bcd_normal = (vc - vb).cross(vd - vb);

	// if the tet vertices abcd form a positive orientation, no need to correct
	// if not, flip them
	double orientation = (vd - va).dot((vb- va).cross(vc - va));
	if (orientation < 0.0) {
		acb_normal *= -1.0;
		adc_normal *= -1.0;
		abd_normal *= -1.0;
		bcd_normal *= -1.0;
	}

	// Computes the area of triangles
	// area = 0.5 * | u x v |
	double acb_area, adc_area, abd_area, bcd_area;
	acb_area = 0.5 * sqrt(acb_normal.dot(acb_normal));
	adc_area = 0.5 * sqrt(adc_normal.dot(adc_normal));
	abd_area = 0.5 * sqrt(abd_normal.dot(abd_normal));
	bcd_area = 0.5 * sqrt(bcd_normal.dot(bcd_normal));

	acb_normal.normalize();
	adc_normal.normalize();
	abd_normal.normalize();
	bcd_normal.normalize();

	this->Nm.col(0) = -(acb_area * acb_normal + adc_area * adc_normal + abd_area * abd_normal) / 3.0;
	this->Nm.col(1) = -(acb_area * acb_normal + abd_area * abd_normal + bcd_area * bcd_normal) / 3.0;
	this->Nm.col(2) = -(acb_area * acb_normal + adc_area * adc_normal + bcd_area * bcd_normal) / 3.0;
	this->Nm.col(3) = -(adc_area * adc_normal + abd_area * abd_normal + bcd_area * bcd_normal) / 3.0;

	return this->Nm;
}

VectorXd Tetrahedron::computeElasticForces(VectorXd f) {
	this->F = computeDeformationGradient();
	// The deformation gradient is available in this->F

	this->P = computePKStress(F, m_mu, m_lambda);
	this->H = -W * P * (Bm.transpose());

	/*if (isInvert && m_isInvertible) {
		this->H = -W * U * P * V.transpose() * (Bm.transpose());
	}*/

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		int rowi = m_nodes[i]->idxM;
		f.segment<3>(rowi) += H.col(i);
		int row3 = m_nodes[3]->idxM;
		f.segment<3>(row3) -= H.col(i);
		//m_nodes[i]->addForce(H.col(i));
		//m_nodes[3]->addForce(-H.col(i));
	}
	return f;
}

bool Tetrahedron::checkNecessityForSVD(double deltaL, double deltaU, Matrix3d F) {
	Matrix3d C = F.transpose() * F;
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

	c1 = F.col(0).transpose() * F.col(0);
	c2 = F.col(0).transpose() * F.col(1);
	c3 = F.col(0).transpose() * F.col(2);
	c4 = F.col(1).transpose() * F.col(1);
	c5 = F.col(1).transpose() * F.col(2);
	c6 = F.col(2).transpose() * F.col(2);

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

		double f_deltaL = pow(deltaL, 3) + a * pow(deltaL, 2) + b * deltaL + c;
		double f_deltaU = pow(deltaU, 3) + a * pow(deltaU, 2) + b * deltaU + c;

		if ((f_deltaL < 0) && (f_deltaU > 0) && (deltaL < alpha) && (beta < deltaU)) {
			// This element has a small deformation, no need to do SVD
			m_isSVD = false;
			//cout << "no no " << endl;
		}
		else {
			// large deformation, need to do SVD
			m_isSVD = true;
			//cout << "Do svd! " << endl;
			//cout << "F" << this->F << endl;
			//cout << "same" << same << endl;
			//cout << "a" << a << endl;
			//cout << "b " << b << endl;
			//cout << "tem" << a * a - 3.0 *b << endl;
			//cout << "fdeltaL: " << f_deltaL << endl;
			//cout << "fdeltaU: " << f_deltaU << endl;
			//cout << "alpha: " << alpha << endl;
			//cout << "beta: " << beta << endl;
		}
	}
	else {
		// no solutions
		m_isSVD = true;
		//cout << "Do svd!no sol " << endl;

	}
	//m_isSVD = true;
	m_isSVD = true;

	return m_isSVD;
}

VectorXd Tetrahedron::computeInvertibleElasticForces(VectorXd f) {
	bool print = false;

	this->F = computeDeformationGradient();
	if (abs(this->F(0, 1) > 3000)) {
		print = true;
	}
	if (print) {
		cout << "F" << this->F << endl;
			cout << "P suppose: " << endl<<computePKStress(F, m_mu, m_lambda);
	}
	
	// The deformation gradient is available in this->F
	if (this->F.determinant() <= 0.0) {
		m_isInverted = true;
		m_isSVD = true;
	}
	else {
		// We need an additional step to decide whether the element undergoes large deformation
		m_isInverted = false;
		checkNecessityForSVD(m_delta_L, m_delta_U, this->F);
		// The necessity of SVD is available in m_isSVD
	}

	if (m_isSVD) {
		// SVD on the deformation gradient
		int modifiedSVD = 1;

		if (!SVD(this->F, this->U, Fhats, this->V, 1e-8, modifiedSVD)) {
			//cout << "error in svd " << endl;
		}
		this->Fhat = Fhats.asDiagonal();
		if (print) {
			cout << "Fhat" << this->Fhat << endl;
		}

		// SVD result is available in this->U, this->V, Fhat_vec, this->Fhat
		// clamp if below the principal stretch threshold
		clamped = 0;
		for (int i = 0; i < 3; i++)
		{
			if (this->Fhat(i, i) < Fthreshold)
			{
				this->Fhat(i, i) = Fthreshold;
				clamped |= (1 << i);
			}
		}
		if (print) {
			cout << "Fhat" << this->Fhat << endl;
		}
		//clamped = 0; // disable clamping

		// Computes the internal forces
		// Computes P first and computes the nodal forces G=PBm in section 4 of [Irving 04]

		// Computes the diagonal P tensor
		this->Phat = computeInvertiblePKStress(this->Fhat, m_mu, m_lambda);

		if (print) {
			cout << "Phat" << this->Phat << endl;
		}
		// P = U * diag(Phat) * V'
		this->P = this->U * this->Phat * this->V.transpose();

		if (print) {
			cout << "P" << this->P << endl;
			Matrix3d ttemp = computePKStress(F, m_mu, m_lambda);
			this->H = -W * ttemp * (Bm.transpose());
			cout << "H" << this->H << endl;
			cout << "PNm" << this->P * this->Nm.block<3, 3>(0, 0) << endl;
		}
		

	}
	else {
		// Don't do SVD
		this->P = computePKStress(this->F, m_mu, m_lambda);
		
	}
	
	// Computes the nodal forces by G=PBm=PNm
	for (int i = 0; i < (int)m_nodes.size()-1; i++) {
		int rowi = m_nodes[i]->idxM;
		f.segment<3>(rowi) += this->P * this->Nm.col(i);
		int row3 = m_nodes[3]->idxM;
		f.segment<3>(row3) -= this->P * this->Nm.col(i);
		//cout << "P" << this->P << endl;
		//cout << "Nm" << this->Nm.col(i) << endl;
		//cout << "f" << endl << this->P * this->Nm.col(i) << endl;

		if (m_isInverted) {
			//cout <<"f"<< endl<< this->P * this->Nm.col(i) << endl;
		}
		 
	}
	return f;
}

void Tetrahedron::computeForceDifferentials(VectorXd dx, VectorXd& df) {
	this->F = computeDeformationGradient();
	this->dF = computeDeformationGradientDifferential(dx);	
	this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);
	//cout << "dP" << endl << this->dP << endl;
	//cout << "df " << endl << this->dP * Nm.block<3, 3>(0, 0) << endl;
	this->dH = -W * dP * (Bm.transpose());
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		df.segment<3>(3 * m_nodes[i]->i) += this->dH.col(i);
		df.segment<3>(3 * m_nodes[3]->i) -= this->dH.col(i);
	}
}

void Tetrahedron::computeForceDifferentials(MatrixXd &K_global) {
	this->F = computeDeformationGradient();
	this->K.setZero();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			this->dDs.setZero();
			this->dDs(j, i) = 1.0;
			this->dF = this->dDs * this->Bm;			
			this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);

			Matrix3d hessian = this->dP * this->Nm.block<3, 3>(0, 0);
			for (int t = 0; t < 3; t++) {
				Vector3d temp = hessian.col(t);
				this->K.block<3, 1>(t * 3, i * 3 + j) = temp;
				this->K.block<3, 1>(9, i * 3 + j) = -temp;
				this->K.block<1, 3>(i * 3 + j, 9) = -temp;
			}
		}
	}

	this->K.block<3, 3>(9, 9) = -this->K.block<3, 3>(0, 9) - this->K.block<3, 3>(3, 9) - this->K.block<3, 3>(6, 9);
	//cout << "Ki" << endl << this->K << endl;
	VectorXcd ev = this->K.eigenvalues();
	cout << ev << endl;

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

void Tetrahedron::computeInvertibleForceDifferentials(VectorXd dx, VectorXd &df) {

	this->dF = computeDeformationGradientDifferential(dx);
	Matrix3d UTdFV = this->U.transpose() * this->dF * this->V;

	this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);	
	//cout << "inv:dPhat" << endl << this->dPhat << endl;
	this->dP = this->U * this->dPhat * this->V.transpose();
	Matrix3d hessian = this->dP * this->Nm.block<3, 3>(0, 0);
	//cout << hessian << endl;
	// modify hessian to compute correct values if in the inversion handling regime
	//hessian = clampHessian(hessian, clamped);

	//cout << "inv:dP" << endl << this->dP << endl;
	//cout << "inv:df " << endl << this->dP * Nm.block<3, 3>(0, 0) << endl;
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		df.segment<3>(3 * m_nodes[i]->i) += hessian.col(i);
		df.segment<3>(3 * m_nodes[3]->i) -= hessian.col(i);
		
		//df.segment<3>(3 * m_nodes[i]->i) += this->dP * this->Nm.col(i);
		//df.segment<3>(3 * m_nodes[3]->i) -= this->dP * this->Nm.col(i);
		if (m_isInverted) {
			//cout << "df:"<< endl<< hessian.col(i) << endl;
		}
	}
}

void Tetrahedron::computeInvertibleForceDifferentials(VectorXd dx, int row, int col, MatrixXd &K) {

	this->dF = computeDeformationGradientDifferential(dx);
	Matrix3d UTdFV = this->U.transpose() * this->dF * this->V;
	this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);
	
	this->dP = this->U * this->dPhat * this->V.transpose();
	Matrix3d hessian = this->dP * this->Nm.block<3, 3>(0, 0);

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		Vector3d temp = hessian.col(i);
		for (int j = 0; j < 3; ++j) {

			K(row + 3 * m_nodes[i]->i + j, col) += temp(j);
			K(row + 3 * m_nodes[3]->i + j, col) += -temp(j);
		}
	}
}

void Tetrahedron::computeInvertibleForceDifferentials(MatrixXd &K_global) {
	this->K.setZero();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			this->dDs.setZero();
			this->dDs(j, i) = 1.0;
			this->dF = this->dDs * this->Bm;
			Matrix3d UTdFV = this->U.transpose() * this->dF * this->V;
			this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);

			this->dP = this->U * this->dPhat * this->V.transpose();
			Matrix3d hessian = this->dP * this->Nm.block<3, 3>(0, 0);
			for (int t = 0; t < 3; t++) {
				Vector3d temp = hessian.col(t);
				this->K.block<3, 1>(t * 3, i * 3 + j) = temp;
				this->K.block<3, 1>(9, i * 3 + j) = -temp;
				this->K.block<1, 3>(i * 3 + j, 9) = -temp;
			}
		}
	}

	this->K.block<3, 3>(9, 9) = -this->K.block<3, 3>(0, 9) - this->K.block<3, 3>(3, 9) - this->K.block<3, 3>(6, 9);
	//cout << "Ki" << endl << this->K << endl;
	VectorXcd ev = this->K.eigenvalues();
	cout << ev << endl; 

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

void Tetrahedron::computeTetK(MatrixXd &K_global) {
	double K[144];
	memset(K, 0, sizeof(double) * 144);
	/*
	dP/dF is a column major matrix, but is stored as a 1D vector

	| dP_11/dF_11  dP_11/dF_12  dP_11/dF_13  dP_11/dF_21 ... dP_11/dF_33 |
	| dP_12/dF_11  dP_12/dF_12  dP_12/dF_13  dP_12/dF_21 ... dP_12/dF_33 |
	|                              ...                                   |
	| dP_33/dF_11  dP_33/dF_12  dP_33/dF_13  dP_33/dF_21 ... dP_33/dF_33 |
	*/
	Vector81d dPdF; //in 9x9 matrix format
	Vector81d dGdF; //in 9x9 matrix format

	compute_dPdF(dPdF, this->clamped);
	//cout << dPdF << endl;

	compute_dGdF(this->Nm.col(0), this->Nm.col(1), this->Nm.col(2), dPdF, dGdF);
	//cout << dGdF << endl;
	for (int row = 0; row<9; row++)
	{
		for (int column = 0; column<12; column++)
		{
			double result = 0;
			for (int inner = 0; inner<9; inner++)
			{
				//dGdF is 9x9, and dFdU is 9x12
				result += dGdF[9 * row + inner] * dFdU[12 * inner + column];
			}
			K[12 * column + row] = result;

		}
	}

	for (int row = 0; row < 12; row++)
	{
		//10th column
		K[12 * row + 9] = -K[12 * row + 0] - K[12 * row + 3] - K[12 * row + 6];
		//11th column
		K[12 * row + 10] = -K[12 * row + 1] - K[12 * row + 4] - K[12 * row + 7];
		//12th column
		K[12 * row + 11] = -K[12 * row + 2] - K[12 * row + 5] - K[12 * row + 8];
	}

	Matrix<double, 12, 12, RowMajor> K_;
	K_ = Map<Matrix<double, 12, 12> > (K);
	
	int ia = m_nodes[0]->idxM;
	int ib = m_nodes[1]->idxM;
	int ic = m_nodes[2]->idxM;
	int id = m_nodes[3]->idxM;

	K_global.block<3, 3>(ia, ia) += K_.block<3, 3>(0, 0);
	K_global.block<3, 3>(ia, ib) += K_.block<3, 3>(0, 3);
	K_global.block<3, 3>(ia, ic) += K_.block<3, 3>(0, 6);
	K_global.block<3, 3>(ia, id) += K_.block<3, 3>(0, 9);
	K_global.block<3, 3>(ib, ia) += K_.block<3, 3>(3, 0);
	K_global.block<3, 3>(ib, ib) += K_.block<3, 3>(3, 3);
	K_global.block<3, 3>(ib, ic) += K_.block<3, 3>(3, 6);
	K_global.block<3, 3>(ib, id) += K_.block<3, 3>(3, 9);
	K_global.block<3, 3>(ic, ia) += K_.block<3, 3>(6, 0);
	K_global.block<3, 3>(ic, ib) += K_.block<3, 3>(6, 3);
	K_global.block<3, 3>(ic, ic) += K_.block<3, 3>(6, 6);
	K_global.block<3, 3>(ic, id) += K_.block<3, 3>(6, 9);
	K_global.block<3, 3>(id, ia) += K_.block<3, 3>(9, 0);
	K_global.block<3, 3>(id, ib) += K_.block<3, 3>(9, 3);
	K_global.block<3, 3>(id, ic) += K_.block<3, 3>(9, 6);
	K_global.block<3, 3>(id, id) += K_.block<3, 3>(9, 9);

}

void Tetrahedron::computeEnergyGradient(Vector3d invariants, Vector3d &gradient) {
	if (m_material == NEO_HOOKEAN) {
		double IIIC = invariants[2];
		gradient[0] = 0.5 * m_mu;
		gradient[1] = 0.0;
		gradient[2] = (-0.5 * m_mu + 0.25 * m_lambda * log(IIIC)) / IIIC;

	}
	else if(m_material == STVK) {
		double IC = invariants[0];
		gradient[0] = 0.25 * m_lambda * (IC - 3.0) - 0.5 * m_mu;
		gradient[1] = 0.25 * m_mu;
		gradient[2] = 0.0;



	}else if(m_material == CO_ROTATED){

	}

}

void Tetrahedron::computeEnergyHessian(Vector3d invariants, Vector6d &hessian) {
	if (m_material == NEO_HOOKEAN) {
		double IIIC = invariants[2];
		// 11
		hessian[0] = 0.0;
		// 12
		hessian[1] = 0.0;
		// 13
		hessian[2] = 0.0;
		// 22
		hessian[3] = 0.0;
		// 23
		hessian[4] = 0.0;
		// 33
		hessian[5] = (0.25 * m_lambda + 0.5 * m_mu - 0.25 *  m_lambda * log(IIIC)) / (IIIC * IIIC);
	}
	else if(m_material == STVK) {
		// 11
		hessian[0] = 0.25 * m_lambda;
		// 12
		hessian[1] = 0.0;
		// 13
		hessian[2] = 0.0;
		// 22
		hessian[3] = 0.0;
		// 23
		hessian[4] = 0.0;
		// 33
		hessian[5] = 0.0;
	}
	else if (m_material == CO_ROTATED) {

	}
}

void Tetrahedron::compute_dPdF(Vector81d &dPdF, int clamped) {
	Vector3d sigma = this->Fhats;
	cout << "sigma" << this->Fhats << endl;

	double sigma1square = sigma[0] * sigma[0];
	double sigma2square = sigma[1] * sigma[1];
	double sigma3square = sigma[2] * sigma[2];

	Vector3d invariants;
	invariants[0] = sigma1square + sigma2square + sigma3square;
	invariants[1] = (sigma1square * sigma1square +
		sigma2square * sigma2square +
		sigma3square * sigma3square);
	invariants[2] = sigma1square * sigma2square * sigma3square;

	//double E[3];
	//E[0] = 0.5 * (Fhats[el][0] * Fhats[el][0] - 1);
	//E[1] = 0.5 * (Fhats[el][1] * Fhats[el][1] - 1);
	//E[2] = 0.5 * (Fhats[el][2] * Fhats[el][2] - 1);
	Vector3d gradient;
	computeEnergyGradient(invariants, gradient);

	/*
	in order (11,12,13,22,23,33)
	| 11 12 13 |   | 0 1 2 |
	| 21 22 23 | = | 1 3 4 |
	| 31 32 33 |   | 2 4 5 |
	*/
	Vector6d hessian;
	computeEnergyHessian(invariants, hessian);

	// modify hessian to compute correct values if in the inversion handling regime
	if (clamped & 1) // first lambda was clamped (in inversion handling)
	{
		hessian[0] = hessian[1] = hessian[2] = 0.0;
	}

	if (clamped & 2) // second lambda was clamped (in inversion handling)
	{
		hessian[1] = hessian[3] = hessian[4] = 0.0;
	}

	if (clamped & 4) // third lambda was clamped (in inversion handling)
	{
		hessian[0] = hessian[1] = hessian[2] = hessian[4] = hessian[5] = 0.0;
	}

	double alpha11 = 2.0 * gradient[0] + 8.0 * sigma1square * gradient[1];
	double alpha22 = 2.0 * gradient[0] + 8.0 * sigma2square * gradient[1];
	double alpha33 = 2.0 * gradient[0] + 8.0 * sigma3square * gradient[1];
	double alpha12 = 2.0 * gradient[0] + 4.0 * (sigma1square + sigma2square) * gradient[1];
	double alpha13 = 2.0 * gradient[0] + 4.0 * (sigma1square + sigma3square) * gradient[1];
	double alpha23 = 2.0 * gradient[0] + 4.0 * (sigma2square + sigma3square) * gradient[1];

	double beta11 = 4.0 * sigma1square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma1square;
	double beta22 = 4.0 * sigma2square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma2square;
	double beta33 = 4.0 * sigma3square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma3square;
	double beta12 = 4.0 * sigma[0] * sigma[1] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[0] * sigma[1]);
	double beta13 = 4.0 * sigma[0] * sigma[2] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[0] * sigma[2]);
	double beta23 = 4.0 * sigma[1] * sigma[2] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[1] * sigma[2]);

	double gamma11 = gammaValue(0, 0, sigma, invariants, gradient, hessian);
	double gamma22 = gammaValue(1, 1, sigma, invariants, gradient, hessian);
	double gamma33 = gammaValue(2, 2, sigma, invariants, gradient, hessian);
	double gamma12 = gammaValue(0, 1, sigma, invariants, gradient, hessian);
	double gamma13 = gammaValue(0, 2, sigma, invariants, gradient, hessian);
	double gamma23 = gammaValue(1, 2, sigma, invariants, gradient, hessian);

	double x1111, x2222, x3333;
	double x2211, x3311, x3322;
	double x2121, x3131, x3232;
	double x2112, x3113, x3223;

	x1111 = alpha11 + beta11 + gamma11;
	x2222 = alpha22 + beta22 + gamma22;
	x3333 = alpha33 + beta33 + gamma33;

	x2211 = gamma12;
	x3311 = gamma13;
	x3322 = gamma23;

	x2121 = alpha12;
	x3131 = alpha13;
	x3232 = alpha23;

	x2112 = beta12;
	x3113 = beta13;
	x3223 = beta23;

	//if (enforceSPD)
	//{
	//	FixPositiveIndefiniteness(x1111, x2211, x3311, x2222, x3322, x3333);
	//	FixPositiveIndefiniteness(x2121, x2112);
	//	FixPositiveIndefiniteness(x3131, x3113);
	//	FixPositiveIndefiniteness(x3232, x3223);
	//}

	//double dPdF_atFhat[81];
	Vector81d dPdF_atFhat;
	dPdF_atFhat.setZero();
	//memset(dPdF_atFhat, 0, sizeof(double) * 81);
	dPdF_atFhat[tensor9x9Index(0, 0, 0, 0)] = x1111;
	dPdF_atFhat[tensor9x9Index(0, 0, 1, 1)] = x2211;
	dPdF_atFhat[tensor9x9Index(0, 0, 2, 2)] = x3311;

	dPdF_atFhat[tensor9x9Index(1, 1, 0, 0)] = x2211;
	dPdF_atFhat[tensor9x9Index(1, 1, 1, 1)] = x2222;
	dPdF_atFhat[tensor9x9Index(1, 1, 2, 2)] = x3322;

	dPdF_atFhat[tensor9x9Index(2, 2, 0, 0)] = x3311;
	dPdF_atFhat[tensor9x9Index(2, 2, 1, 1)] = x3322;
	dPdF_atFhat[tensor9x9Index(2, 2, 2, 2)] = x3333;

	dPdF_atFhat[tensor9x9Index(0, 1, 0, 1)] = x2121;
	dPdF_atFhat[tensor9x9Index(0, 1, 1, 0)] = x2112;

	dPdF_atFhat[tensor9x9Index(1, 0, 0, 1)] = x2112;
	dPdF_atFhat[tensor9x9Index(1, 0, 1, 0)] = x2121;

	dPdF_atFhat[tensor9x9Index(0, 2, 0, 2)] = x3131;
	dPdF_atFhat[tensor9x9Index(0, 2, 2, 0)] = x3113;

	dPdF_atFhat[tensor9x9Index(2, 0, 0, 2)] = x3113;
	dPdF_atFhat[tensor9x9Index(2, 0, 2, 0)] = x3131;

	dPdF_atFhat[tensor9x9Index(1, 2, 1, 2)] = x3232;
	dPdF_atFhat[tensor9x9Index(1, 2, 2, 1)] = x3223;

	dPdF_atFhat[tensor9x9Index(2, 1, 1, 2)] = x3223;
	dPdF_atFhat[tensor9x9Index(2, 1, 2, 1)] = x3232;

	/*
	| P_00 P_01 P_02 |        | F_00 F_01 F_02 |
	if P= | P_10 P_11 P_12 | and F= | F_10 F_11 F_12 |
	| P_20 P_21 P_22 |        | F_20 F_21 F_22 |
	| dP_00/dF_00  dP_00/dF_01 dP_00/dF_02 dP_00/dF_10 ... dP00/dF_22 |
	| dP_01/dF_00  dP_01/dF_01 dP_01/dF_02 dP_01/dF_10 ... dP01/dF_22 |
	| dP_02/dF_00  dP_02/dF_01 dP_02/dF_02 dP_02/dF_10 ... dP02/dF_22 |
	| dP_10/dF_00  dP_10/dF_01 dP_10/dF_02 dP_10/dF_10 ... dP10/dF_22 |
	|                               ...                               |
	| dP_22/dF_00  dP_22/dF_01 dP_22/dF_02 dP_22/dF_10 ... dP22/dF_22 |
	*/

	Matrix3d UT = this->U.transpose(); // trans(*U);
	Matrix3d VT = this->V.transpose(); // trans(*V);


	double eiejVector[9];
	memset(eiejVector, 0, sizeof(double) * 9);

	dPdF.setZero();
	for (int column = 0; column < 9; column++)
	{
		eiejVector[column] = 1.0;
		Matrix<double, 3, 3, RowMajor> ei_ej;
		
		ei_ej = Map<Matrix<double, 3, 3> >(eiejVector);
		//cout << ei_ej << endl;

		Matrix3d ut_eiej_v = UT*ei_ej*(this->V);

		Vector9d ut_eiej_v_TeranVector;	//in Teran order
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[0]] = ut_eiej_v(0, 0);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[1]] = ut_eiej_v(0, 1);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[2]] = ut_eiej_v(0, 2);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[3]] = ut_eiej_v(1, 0);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[4]] = ut_eiej_v(1, 1);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[5]] = ut_eiej_v(1, 2);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[6]] = ut_eiej_v(2, 0);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[7]] = ut_eiej_v(2, 1);
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[8]] = ut_eiej_v(2, 2);

		
		double dPdF_resultVector[9]; // not in Teran order

		for (int innerRow = 0; innerRow < 9; innerRow++)
		{
			double tempResult = 0.0;
			for (int innerColumn = 0; innerColumn < 9; innerColumn++)
			{
				tempResult += dPdF_atFhat[innerRow * 9 + innerColumn] *
					ut_eiej_v_TeranVector[innerColumn];

				//cout << dPdF_atFhat[innerRow * 9 + innerColumn] << endl;
				//cout << ut_eiej_v_TeranVector[innerColumn] << endl;
			}
			dPdF_resultVector[teranToRowMajorMatrix[innerRow]] = tempResult;
			
		}

		Matrix<double, 3, 3, RowMajor> dPdF_resultMatrix;
		dPdF_resultMatrix = Map<Matrix<double, 3, 3> >(dPdF_resultVector);



		Matrix3d u_dpdf_vt = this->U * dPdF_resultMatrix * VT;
		dPdF[column + 0] = u_dpdf_vt(0, 0);
		dPdF[column + 9] = u_dpdf_vt(0, 1);
		dPdF[column + 18] = u_dpdf_vt(0, 2);
		dPdF[column + 27] = u_dpdf_vt(1, 0);
		dPdF[column + 36] = u_dpdf_vt(1, 1);
		dPdF[column + 45] = u_dpdf_vt(1, 2);
		dPdF[column + 54] = u_dpdf_vt(2, 0);
		dPdF[column + 63] = u_dpdf_vt(2, 1);
		dPdF[column + 72] = u_dpdf_vt(2, 2);
		// reset
		eiejVector[column] = 0.0;
	}
}

void Tetrahedron::compute_dGdF(const Vector3d &b0, const Vector3d &b1, const Vector3d &b2, const Vector81d &dPdF, Vector81d &dGdF) {
	//Both G and F are 3x3 matrices, so dGdF has 81 entries
	dGdF.setZero();
	/*
			| ga_x gb_x gc_x |   | 0 1 2 |
	if G =	| ga_y gb_y gc_y | = | 3 4 5 |
			| ga_z gb_z gc_z |   | 6 7 8 |
	where ga, gb, gc are the nodal forces at vertex a,b,c
			| ba_0 bb_0 bc_0 |   | 0 1 2 |
	and B = | ba_1 bb_1 bc_1 | = | 3 4 5 |
			| ba_2 bb_2 bc_2 |   | 6 7 8 |
			| dga_x/dF_00 dga_x/dF_01 dga_x/dF_02 dga_x/dF_10 ... dga_x/dF_22 |
			| dga_y/dF_00 dga_y/dF_01 dga_y/dF_02 dga_y/dF_10 ... dga_y/dF_22 |
			| dga_z/dF_00 dga_z/dF_01 dga_z/dF_02 dga_z/dF_10 ... dga_z/dF_22 |
	dGdF =	| dgb_x/dF_00 dgb_x/dF_01 dgb_x/dF_02 dgb_x/dF_10 ... dgb_x/dF_22 |
			|                                 ...                             |
			| dgc_z/dF_00 dgc_z/dF_01 dgc_z/dF_02 dgc_z/dF_10 ... dgc_z/dF_22 |
	*/
	Matrix3d bVec;
	bVec.col(0) = b0;
	bVec.col(1) = b1;
	bVec.col(2) = b2;
	//dga_x/dF, dga_y/dF, dga_z/dF
	//dgb_x/dF, dgb_y/dF, dgb_z/dF
	//dgc_x/dF, dgc_y/dF, dgc_z/dF
	
	for (int abc = 0; abc<3; abc++)
		for (int i = 0; i<3; i++)
			for (int column = 0; column<9; column++)
				for (int k = 0; k<3; k++)
					dGdF[27 * abc + 9 * i + column] += dPdF[(3 * i + k) * 9 + column] * bVec(k, abc);


}


	/*
	The "i" goes from 0 to 2 inclusively
	The "j" goes from 0 to 2 inclusively
	See [Teran 05].
	*/
double Tetrahedron::gammaValue(int i, int j, Vector3d sigma, Vector3d invariants, Vector3d gradient, Vector6d hessian)
{
	/*
	The hessian is in order (11,12,13,22,23,33)
	| 11 12 13 |   | 0 1 2 |
	| 21 22 23 | = | 1 3 4 |
	| 31 32 33 |   | 2 4 5 |
	*/

	Vector3d tempGammaVec1;
	tempGammaVec1[0] = 2.0 * sigma[i];
	tempGammaVec1[1] = 4.0 * sigma[i] * sigma[i] * sigma[i];
	tempGammaVec1[2] = 2.0 * invariants[2] / sigma[i];

	Vector3d tempGammaVec2;
	tempGammaVec2[0] = 2.0 * sigma[j];
	tempGammaVec2[1] = 4.0 * sigma[j] * sigma[j] * sigma[j];
	tempGammaVec2[2] = 2.0 * invariants[2] / sigma[j];

	Vector3d productResult;
	productResult[0] = (tempGammaVec2[0] * hessian[0] + tempGammaVec2[1] * hessian[1] +
		tempGammaVec2[2] * hessian[2]);
	productResult[1] = (tempGammaVec2[0] * hessian[1] + tempGammaVec2[1] * hessian[3] +
		tempGammaVec2[2] * hessian[4]);
	productResult[2] = (tempGammaVec2[0] * hessian[2] + tempGammaVec2[1] * hessian[4] +
		tempGammaVec2[2] * hessian[5]);

	return (tempGammaVec1[0] * productResult[0] + tempGammaVec1[1] * productResult[1] + tempGammaVec1[2] * productResult[2] + 4.0 * invariants[2] * gradient[2] / (sigma[i] * sigma[j]));
}

int Tetrahedron::tensor9x9Index(int i, int j, int m, int n)
{
	/*
	|  dP_0/dF_0  dP_0/dF_4  dP_0/dF_8  ...  dP_0/dF_5  |
	|  dP_4/dF_0  dP_4/dF_4  dP_4/dF_8  ...  dP_4/dF_5  |
	|                         ...                       |
	|  dP_5/dF_0  dP_5/dF_4  dP_5/dF_8  ...  dP_5/dF_5  |
	*/
	int rowIndex_in9x9Matrix = rowMajorMatrixToTeran[3 * i + j];
	int columnIndex_in9x9Matrix = rowMajorMatrixToTeran[3 * m + n];
	return (9 * rowIndex_in9x9Matrix + columnIndex_in9x9Matrix);
}

void Tetrahedron::computeInvertibleForceDifferentialsSparse(VectorXd dx, int row, int col, vector<T> &K_) {
	
	this->dF = computeDeformationGradientDifferential(dx);
	if (m_isSVD) {
		Matrix3d UTdFV = this->U.transpose() * this->dF * this->V;
		this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);
		this->dP = this->U * this->dPhat * this->V.transpose();

	}
	else
	{
		this->dP = computePKStressDerivative(this->F, this->dF, m_mu, m_lambda);

	}
		
	Matrix3d hessian = this->dP * this->Nm.block<3, 3>(0, 0);
	//

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		Vector3d temp = hessian.col(i);

		for (int j = 0; j < 3; ++j) {
			//cout << j << endl;
			K_.push_back(T(row + 3 * m_nodes[i]->i + j, col, temp(j)));
			//cout << temp(j) << endl;
			K_.push_back(T(row + 3 * m_nodes[3]->i + j, col, -temp(j)));
		}
	}
}

void Tetrahedron::computeForceDifferentialsSparse(VectorXd dx, int row, int col, vector<T> &K_) {
	this->F = computeDeformationGradient();
	this->dF = computeDeformationGradientDifferential(dx);
	this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);
	this->dH = -W * dP * (Bm.transpose());

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		Vector3d temp = this->dH.col(i);
		for (int j = 0; j < 3; ++j) {
			K_.push_back(T(row + 3 * m_nodes[i]->i + j, col, temp(j)));
			K_.push_back(T(row + 3 * m_nodes[3]->i + j, col, -temp(j)));
		}
	}
}

Matrix3d Tetrahedron::computeInvertiblePKStress(Matrix3d F, double mu, double lambda) {
	Vector3d invariants;
	Vector3d lambda1, lambda2;
	lambda1 << F(0, 0), F(1, 1), F(2, 2);

	lambda2 << lambda1(0) * lambda1(0), lambda1(1) * lambda1(1), lambda1(2) * lambda1(2);
	double IC = lambda2(0) + lambda2(1) + lambda2(2);
	double IIC = lambda2(0) * lambda2(0) + lambda2(1) * lambda2(1) + lambda2(2) * lambda2(2);
	double IIIC = lambda2(0) * lambda2(1) * lambda2(2);

	invariants << IC, IIC, IIIC;

	Vector3d dPsidIV;
	dPsidIV << 0.5 * mu, 0.0, (-0.5 * mu + 0.25 * lambda * log(IIIC)) / IIIC;

	// PDiag = [ dI / dlambda ]^T * dPsidI
	Matrix3d matM;
	matM << 2.0 * lambda1(0), 2.0 * lambda1(1), 2.0 * lambda1(2),
		4.0 * lambda1(0) * lambda1(0) * lambda1(0), 4.0 * lambda1(1) * lambda1(1) * lambda1(1), 4.0 * lambda1(2) * lambda1(2) * lambda1(2),
		2.0 * lambda1(0) * lambda2(1) * lambda2(2), 2.0 * lambda1(1) * lambda2(0) * lambda2(2), 2.0 * lambda1(2) * lambda2(0) * lambda2(1);

	Vector3d result;
	result = matM.transpose() * dPsidIV;
	this->Phat = result.asDiagonal();
	return this->Phat;

}

Matrix3d Tetrahedron::computePKStress(Matrix3d F, double mu, double lambda) {

	Matrix3d E = Matrix3d::Zero();
	Matrix3d P = Matrix3d::Zero();
	Matrix3d I = Matrix3d::Identity();

	switch (m_material)
	{
	case LINEAR:
	{
		E = 0.5 * (F + F.transpose()) - I;
		psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		P = 2.0 * mu * E + lambda * E.trace() * I;
		break;
	}

	case NEO_HOOKEAN:
	{
		double I1 = (F.transpose() * F).trace();
		double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		double I3 = (F.transpose() * F).determinant();
		double J = sqrt(I3);
		psi = 1.0 / 2.0 * mu *(I1 - 3.0) - mu * log(J) + 1.0 / 2.0 * lambda * log(J)*log(J);
		P = mu * (F - F.inverse().transpose()) + lambda * log(J)*(F.inverse().transpose());


		break;
	}

	case STVK:
	{
		E = 0.5 * (F.transpose() * F - I);
		psi = mu * E.norm()*E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		P = F * (2.0 * mu * E + lambda * E.trace() * I);
		break;
	}

	case CO_ROTATED:
	{
		// Polar decomposition
		Matrix3d A = F.adjoint() * F;
		SelfAdjointEigenSolver<Matrix3d> es(A);
		Matrix3d S = es.operatorSqrt();
		Matrix3d R = F * S.inverse();

		E = S - I;
		psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		//P = R * (2.0 * mu * E + lambda * E.trace() * I);
		P = 2.0 * mu * (F - R) + lambda * (R.transpose() * F - I).trace() * R;
		break;
	}

	default:
	{
		break;
	}
	}

	return P;
}

Matrix3d Tetrahedron::computePKStressDerivative(Matrix3d F, Matrix3d dF, double mu, double lambda) {
	Matrix3d E = Matrix3d::Zero();
	Matrix3d P = Matrix3d::Zero();
	Matrix3d dE = Matrix3d::Zero();
	Matrix3d dP = Matrix3d::Zero();
	Matrix3d I3 = Matrix3d::Identity();
	switch (m_material) {
	case CO_ROTATED:
	{
		Matrix3d A = F.adjoint() * F;
		SelfAdjointEigenSolver<Matrix3d> es(A);
		Matrix3d S = es.operatorSqrt();
		Matrix3d R = F * S.inverse();
		E = S - I3;
		P = 2.0 * mu *(F - R) + lambda * (R.transpose()*F - I3).trace() * R;
		dP = 2.0 * mu * dF + lambda * (R.transpose()*dF).trace() * R;
		break;
	}

	case STVK:
	{
		E = 1.0 / 2.0 * (F.transpose() * F - I3);
		dE = 1.0 / 2.0 * (dF.transpose() * F + F.transpose() * dF);
		P = F * (2.0 * mu * E + lambda * E.trace() * I3);
		dP = dF * (2.0 * mu * E + lambda * E.trace() * I3) + F * (2.0 * mu * dE + lambda * dE.trace() * I3);
		break;
	}

	case NEO_HOOKEAN:
	{	
		/*double I1 = (F.transpose() * F).trace();
		double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		double I3 = (F.transpose() * F).determinant();
		double J = sqrt(I3);
		psi = 1.0 / 2.0 * mu *(I1 - 3.0) - mu * log(J) + 1.0 / 2.0 * lambda * log(J)*log(J);
		P = mu * (F - F.inverse().transpose()) + lambda * log(J)*(F.inverse().transpose());
*/
		//double I1 = (F.norm()) * (F.norm());
		//double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		MatrixXd FT = F.transpose();
		MatrixXd FIT = F.inverse().transpose();
		//double I3 = (F.transpose() * F).determinant();
		double I3 = (FT * F).determinant();
		double J = sqrt(I3);
		P = mu * (F - FIT) + lambda * log(J) * FIT;
		dP = mu * dF + (mu - lambda * log(J)) * FIT * (dF.transpose()) * FIT + lambda * ((F.inverse() * dF)).trace() * FIT;
		//P = mu * (F - (F.inverse().transpose())) + lambda * log(J) * (F.inverse().transpose());
		//dP = mu * dF + (mu - lambda * log(J)) * (F.inverse().transpose()) * (dF.transpose()) * (F.inverse().transpose()) + lambda * ((F.inverse() * dF)).trace() * (F.inverse().transpose());
		break;
	}
	case LINEAR:
	{
		E = 1.0 / 2.0 * (F + F.transpose()) - I3;
		dE = 1.0 / 2.0 * (dF + dF.transpose());
		P = 2.0 * mu * E + lambda * E.trace() * I3;
		dP = 2.0 * mu * dE + lambda * dE.trace() * I3;
		break;
	}
	default:
		break;
	}

	return dP;
}

void Tetrahedron::diagDeformationGradient(Matrix3d F_) {
	Matrix3f F = F_.cast<float>();
	float a11, a12, a13, a21, a22, a23, a31, a32, a33;

	a11 = F(0, 0); a12 = F(0, 1); a13 = F(0, 2);
	a21 = F(1, 0); a22 = F(1, 1); a23 = F(1, 2);
	a31 = F(2, 0); a32 = F(2, 1); a33 = F(2, 2);

	float u11, u12, u13,
		u21, u22, u23,
		u31, u32, u33;

	float s11, s12, s13,
		s21, s22, s23,
		s31, s32, s33;

	float v11, v12, v13,
		v21, v22, v23,
		v31, v32, v33;

	svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
		u11, u12, u13, u21, u22, u23, u31, u32, u33,
		s11, s12, s13, s21, s22, s23, s31, s32, s33,
		v11, v12, v13, v21, v22, v23, v31, v32, v33);

	Matrix3f U_, V_, S_;
	U_ << u11, u12, u13,
		u21, u22, u23,
		u31, u32, u33;
	V_ << v11, v12, v13, v21, v22, v23, v31, v32, v33;
	S_ << s11, s12, s13, s21, s22, s23, s31, s32, s33;
	this->U = U_.cast<double>();
	this->V = V_.cast<double>();
	this->Fhat = S_.cast<double>();
}

double Tetrahedron::computeEnergy() {
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	this->P = computePKStress(F, m_mu, m_lambda);
	this->m_energy = W * psi;
	return this->m_energy;
}

Matrix3d Tetrahedron::clampHessian(Matrix3d &hessian, int clamped) {
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

void Tetrahedron::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	if (m_isInverted) {
		prog->bind();
		for (int i = 0; i < 4; i++) {
			auto node = m_nodes[i];
			node->draw(MV, prog);
		}
		prog->unbind();
	}

}
