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

Tetrahedron::Tetrahedron(double young, double poisson, double density, Material material, const vector<shared_ptr<Node>> &nodes) :
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
	this->W = (1.0 / 6.0 * Dm.determinant());
	m_mass = (this->W * this->m_density);

	// Distribute 1/4 mass to each node
	for (int i = 0; i < (int)nodes.size(); i++) {
		nodes[i]->m += this->m_mass * 0.25;
	}

	//computeAreaWeightedVertexNormals();

	m_delta_L = 0.001;
	m_delta_U = 2.0;

	//m_delta_L = 0.0;
	//m_delta_U = 1005.0;
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
		this->Phat = computePKStress(this->Fhat, m_mu, m_lambda);

		if (print) {
			cout << "Phat" << this->Phat << endl;
		}
		// P = U * diag(Phat) * V'
		this->P = this->U * this->Phat * this->V.transpose();

		if (print) {
			Matrix3d ttemp = computePKStress(F, m_mu, m_lambda);
			this->H = -W * ttemp * (Bm.transpose());

		}
	}
	else {
		// Don't do SVD
		this->P = computePKStress(this->F, m_mu, m_lambda);
	}
	
	this->H = -W * this->P * (Bm.transpose());

	// Computes the nodal forces by G=PBm=PNm
	for (int i = 0; i < (int)m_nodes.size()-1; i++) {
		int rowi = m_nodes[i]->idxM;
		f.segment<3>(rowi) += this->H.col(i);
		int row3 = m_nodes[3]->idxM;
		f.segment<3>(row3) -= this->H.col(i);

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
			this->dH = -W * dP * (Bm.transpose());
			for (int t = 0; t < 3; t++) {
				Vector3d temp = this->dH.col(t);
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

void Tetrahedron::computeInvertibleForceDifferentials(VectorXd dx, VectorXd &df) {

	this->dF = computeDeformationGradientDifferential(dx);
	Matrix3d UTdFV = this->U.transpose() * this->dF * this->V;

	this->dPhat = computePKStressDerivative(this->Fhat, UTdFV, m_mu, m_lambda);	
	this->dP = this->U * this->dPhat * this->V.transpose();
	this->dH = -W * dP * (Bm.transpose());

	//modify hessian to compute correct values if in the inversion handling regime
	//hessian = clampHessian(hessian, clamped);

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		df.segment<3>(3 * m_nodes[i]->i) += dH.col(i);
		df.segment<3>(3 * m_nodes[3]->i) -= dH.col(i);
		
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
	this->dH = -W * dP * (Bm.transpose());

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		Vector3d temp = dH.col(i);
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
			this->dH = -W * dP * (Bm.transpose());

			//Matrix3d hessian = this->dP * this->Nm.block<3, 3>(0, 0);
			for (int t = 0; t < 3; t++) {
				Vector3d temp = dH.col(t);
				this->K.block<3, 1>(t * 3, i * 3 + j) = temp;
				this->K.block<3, 1>(9, i * 3 + j) -= temp;
				this->K.block<1, 3>(i * 3 + j, 9) -= temp;
			}
		}
	}

	this->K.block<3, 3>(9, 9) = - this->K.block<3, 3>(0, 9) - this->K.block<3, 3>(3, 9) - this->K.block<3, 3>(6, 9);
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
	this->dH = -W * dP * (Bm.transpose());

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		//Vector3d temp = dH.col(i);

		for (int j = 0; j < 3; ++j) {
			K_.push_back(T(row + 3 * m_nodes[i]->i + j, col, dH(j, i)));
			K_.push_back(T(row + 3 * m_nodes[3]->i + j, col, -dH(j, i)));
		}
	}
}

void Tetrahedron::computeForceDifferentialsSparse(VectorXd dx, int row, int col, vector<T> &K_) {
	this->F = computeDeformationGradient();
	this->dF = computeDeformationGradientDifferential(dx);
	this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);	
	this->dH = -W * dP * (Bm.transpose());

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		//Vector3d temp = this->dH.col(i);
		for (int j = 0; j < 3; ++j) {
			K_.push_back(T(row + 3 * m_nodes[i]->i + j, col, dH(j, i)));
			K_.push_back(T(row + 3 * m_nodes[3]->i + j, col, -dH(j, i)));
		}
	}
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
		//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		P = 2.0 * mu * E + lambda * E.trace() * I;
		break;
	}

	case NEO_HOOKEAN:
	{
		//double I1 = (F.transpose() * F).trace();
		//double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		double I3 = (F.transpose() * F).determinant();
		double J = sqrt(I3);
		//psi = 1.0 / 2.0 * mu *(I1 - 3.0) - mu * log(J) + 1.0 / 2.0 * lambda * log(J)*log(J);
		P = mu * (F - F.inverse().transpose()) + lambda * log(J)*(F.inverse().transpose());
		break;
	}

	case STVK:
	{
		E = 0.5 * (F.transpose() * F - I);
		//psi = mu * E.norm()*E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
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

		//E = S - I;
		//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
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
		//E = S - I3;
		//P = 2.0 * mu *(F - R) + lambda * (R.transpose()*F - I3).trace() * R;
		dP = 2.0 * mu * dF + lambda * (R.transpose()*dF).trace() * R;
		break;
	}

	case STVK:
	{
		E = 1.0 / 2.0 * (F.transpose() * F - I3);
		dE = 1.0 / 2.0 * (dF.transpose() * F + F.transpose() * dF);
		//P = F * (2.0 * mu * E + lambda * E.trace() * I3);
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
		//P = mu * (F - FIT) + lambda * log(J) * FIT;
		dP = mu * dF + (mu - lambda * log(J)) * FIT * (dF.transpose()) * FIT + lambda * ((F.inverse() * dF)).trace() * FIT;
		//P = mu * (F - (F.inverse().transpose())) + lambda * log(J) * (F.inverse().transpose());
		//dP = mu * dF + (mu - lambda * log(J)) * (F.inverse().transpose()) * (dF.transpose()) * (F.inverse().transpose()) + lambda * ((F.inverse() * dF)).trace() * (F.inverse().transpose());
		break;
	}
	case LINEAR:
	{
		E = 1.0 / 2.0 * (F + F.transpose()) - I3;
		dE = 1.0 / 2.0 * (dF + dF.transpose());
		//P = 2.0 * mu * E + lambda * E.trace() * I3;
		dP = 2.0 * mu * dE + lambda * dE.trace() * I3;
		break;
	}
	default:
		break;
	}

	return dP;
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
