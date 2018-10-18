#include "Tetrahedron.h"
#include "Node.h"
#include <iostream>
#include "svd3.h"

using namespace Eigen;
using namespace std;


Tetrahedron::Tetrahedron()
{

}

Tetrahedron::Tetrahedron(double young, double poisson, double density, Material material, const vector<shared_ptr<Node>> &nodes):
m_young(young), m_poisson(poisson), m_density(density), m_material(material), m_nodes(nodes)
{
	m_mu = m_young / (2.0 * (1.0 + m_poisson));
	m_lambda = m_young * m_poisson / ((1.0 + m_poisson) * (1.0 - 2.0 * m_poisson));

	assert(m_nodes.size() == 4);
	for (int i = 0; i < m_nodes.size() - 1; i++) {
		this->Dm.col(i) = m_nodes[i]->x0 - m_nodes[3]->x0;
	}
	this->Bm = this->Dm.inverse();
	this->W = 1.0 / 6.0 * Dm.determinant();
	m_mass = this->W * this->m_density;

	// Distribute 1/4 mass to each node
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->m += this->m_mass * 0.25;
	}
}

VectorXd Tetrahedron::computeElasticForces(VectorXd f) {
	isInverted();

	for (int i = 0; i < m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	if (isInvert) {
		this->F = this->Fhat;
	}

	this->P = computePKStress(F, m_mu, m_lambda);

	this->H = -W * P * (Bm.transpose());

	if (isInvert) {
		this->H = -W * U * P * V.transpose() * (Bm.transpose());
	}

	for (int i = 0; i < m_nodes.size() - 1; i++) {
		int rowi = m_nodes[i]->idxM;
		f.segment<3>(rowi) += H.col(i);
		int row3 = m_nodes[3]->idxM;
		f.segment<3>(row3) -= H.col(i);

		//m_nodes[i]->addForce(H.col(i));
		//m_nodes[3]->addForce(-H.col(i));
	}

	return f;
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

bool Tetrahedron::isInverted() {
	for (int i = 0; i < m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	if (this->F.determinant() < -0.000001) { // some threshold todo
		//cout << "tet_" << this->i << " is inverted! " << endl;
		diagDeformationGradient(this->F);
		isInvert = true;
	}
	else {
		isInvert = false;
	}

	return isInvert;
}

void Tetrahedron::diagDeformationGradient(Eigen::Matrix3d F_) {
	Eigen::Matrix3f F = F_.cast<float>();
	float a11, a12, a13, a21, a22, a23, a31, a32, a33;

	a11 = F(0, 0); a12 = F(0, 1); a13 = F(0, 2);
	a21 = F(1, 0); a22 = F(1, 1); a23 = F(1, 2);
	a31 = F(2, 0); a32 = F(2, 1); a33 = F(2, 2);

	float u11, u12, u13,
		u21, u22, u23,
		u31, u32, u33;

	float 	s11, s12, s13,
		s21, s22, s23,
		s31, s32, s33;

	float 	v11, v12, v13,
		v21, v22, v23,
		v31, v32, v33;

	svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
		u11, u12, u13, u21, u22, u23, u31, u32, u33,
		s11, s12, s13, s21, s22, s23, s31, s32, s33,
		v11, v12, v13, v21, v22, v23, v31, v32, v33);

	Eigen::Matrix3f U_, V_, S_;
	U_ << u11, u12, u13,
		u21, u22, u23,
		u31, u32, u33;
	V_ << v11, v12, v13, v21, v22, v23, v31, v32, v33;
	S_ << s11, s12, s13, s21, s22, s23, s31, s32, s33;
	this->U = U_.cast<double>();
	this->V = V_.cast<double>();
	this->Fhat = S_.cast<double>();

}


void Tetrahedron::computeForceDifferentials(VectorXd dx, VectorXd& df) {
	for (int i = 0; i < m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	if (isInverted()) {
		this->F = this->Fhat;
	}


	for (int i = 0; i < m_nodes.size() - 1; i++) {
		this->dDs.col(i) = dx.segment<3>(3 * m_nodes[i]->i) - dx.segment<3>(3 * m_nodes[3]->i);
	}

	this->dF = dDs * Bm;
	this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);
	this->dH = -W * dP * (Bm.transpose());
	for (int i = 0; i < m_nodes.size() - 1; i++) {
		df.segment<3>(3 * m_nodes[i]->i) += this->dH.col(i);
		df.segment<3>(3 * m_nodes[3]->i) -= this->dH.col(i);
	}

	/*MatrixXd dFRow(4, 3);
	for (int i = 0; i < 3; ++i) {
	dFRow.row(i) = Bm.row(i);
	dFRow(3, i) = -Bm(0, i) - Bm(1, i) - Bm(2, i);
	}

	K.setZero();
	for (int row = 0; row < 4; ++row) {
	MatrixXd Kb(12, 3);
	Kb.setZero();
	for (int kk = 0; kk < 3; ++kk) {
	Matrix3d dF;
	dF.setZero();
	dF.row(kk) = dFRow.row(row);
	MatrixXd dP = computePKStressDerivative(F, dF, material, mu, lambda);
	dH = -W * dP * (Bm.transpose());
	for (int ii = 0; ii < 3; ii++) {
	for (int ll = 0; ll < 3; ll++) {
	Kb(ii * 3 + ll, kk) = dH(ll, ii);
	}
	Kb(9 + ii, kk) = -dH(ii, 0) - dH(ii, 1) - dH(ii, 2);
	}
	}
	K.block<12, 3>(0, 3 * row) = Kb;
	}*/
}

double Tetrahedron::computeEnergy() {
	isInverted();

	for (int i = 0; i < m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	if (isInvert) {
		this->F = this->Fhat;
	}

	this->P = computePKStress(F, m_mu, m_lambda);
	this->m_energy = W * psi;
	return this->m_energy;
}

Tetrahedron:: ~Tetrahedron() {

}