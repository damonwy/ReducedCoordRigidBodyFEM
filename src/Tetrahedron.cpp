#include "rmpch.h"
#include "Tetrahedron.h"

#include "Node.h"
#include <cmath>
#include "Body.h"

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
	m_enclosed_color << (float)(rand() % 255) / 255.0f, (float)(rand() % 255) / 255.0f, (float)(rand() % 255) / 255.0f;
	m_mu = m_young / (2.0 * (1.0 + m_poisson));
	m_lambda = m_young * m_poisson / ((1.0 + m_poisson) * (1.0 - 2.0 * m_poisson));
	m_isInverted = false;
	// Compute the inverse of the Dm matrix
	// Dm is a 3x3 matrix where the columns are the edge vectors of a tet in rest configuration

	for (int i = 0; i < 3; i++) {
		this->Dm.col(i) = m_nodes[i]->x0 - m_nodes[3]->x0;
	}
	this->Bm = this->Dm.inverse();
	this->BmT = this->Bm.transpose();

}

void Tetrahedron::precompute() {
	// Compute volume
	this->W = 1.0 / 6.0 * Dm.determinant();
	m_mass = this->W * this->m_density;

	// Distribute 1/4 mass to each node
	// Initialize the undeformed position vector
	for (int i = 0; i < 4; i++) {
		m_nodes[i]->m += this->m_mass * 0.25;
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

void Tetrahedron::computeElasticForces() {
	this->F = computeDeformationGradient();
	// The deformation gradient is available in this->F
	this->FTF.noalias() = F.transpose() * F;

	this->P = computePKStress(F, m_mu, m_lambda);
	this->H = -W * P * this->BmT;

	/*if (isInvert && m_isInvertible) {
		this->H = -W * U * P * V.transpose() * this->BmT);
	}*/
}

void Tetrahedron::computeElasticForces(VectorXd &f) {
	computeElasticForces();
	assembleGlobalForceVector(f);
}

void Tetrahedron::assembleGlobalForceVector(Eigen::VectorXd &f) {
	int rowi;
	auto node3 = m_nodes[3];
	int row3 = node3->idxM;

	for (int i = 0; i < (int)m_nodes.size() - 1; ++i) {
		auto node = m_nodes[i];
		rowi = node->idxM;

		f.segment<3>(rowi) += H.col(i);	
		f.segment<3>(row3) -= H.col(i);		

	}
}

void Tetrahedron::computeForceDifferentials(VectorXd dx, VectorXd& df) {
	this->F = computeDeformationGradient();
	this->dF = computeDeformationGradientDifferential(dx);	
	this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);

	this->dH = -W * dP * this->BmT;
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		df.segment<3>(3 * m_nodes[i]->i) += this->dH.col(i);
		df.segment<3>(3 * m_nodes[3]->i) -= this->dH.col(i);
	}
}

void Tetrahedron::computeForceDifferentials() {
	this->F = computeDeformationGradient();
	this->K.setZero();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			this->dDs.setZero();
			this->dDs(j, i) = 1.0;
			this->dF = this->dDs * this->Bm;
			this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);
			this->dH = -W * dP * this->BmT;
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
}

bool Tetrahedron::checkInverted() {
	if (this->m_isInverted) {		
		//cout << "tet " << i << " is inverted!" << endl;
		return true;
	}
	else {
		return false;
	}
}

void Tetrahedron::assembleGlobalStiffnessMatrixSparse(vector<T> &K_) {
	int a = m_nodes[0]->idxM;
	int b = m_nodes[1]->idxM;
	int c = m_nodes[2]->idxM;
	int d = m_nodes[3]->idxM;

	for (int idx = 0; idx < 3; idx++) {
		for (int jdx = 0; jdx < 3; jdx++) {
			K_.push_back(T(a + idx, a + jdx, this->K(0 + idx, 0 + jdx)));
			K_.push_back(T(a + idx, b + jdx, this->K(0 + idx, 3 + jdx)));
			K_.push_back(T(a + idx, c + jdx, this->K(0 + idx, 6 + jdx)));
			K_.push_back(T(a + idx, d + jdx, this->K(0 + idx, 9 + jdx)));

			K_.push_back(T(b + idx, a + jdx, this->K(3 + idx, 0 + jdx)));
			K_.push_back(T(b + idx, b + jdx, this->K(3 + idx, 3 + jdx)));
			K_.push_back(T(b + idx, c + jdx, this->K(3 + idx, 6 + jdx)));
			K_.push_back(T(b + idx, d + jdx, this->K(3 + idx, 9 + jdx)));

			K_.push_back(T(c + idx, a + jdx, this->K(6 + idx, 0 + jdx)));
			K_.push_back(T(c + idx, b + jdx, this->K(6 + idx, 3 + jdx)));
			K_.push_back(T(c + idx, c + jdx, this->K(6 + idx, 6 + jdx)));
			K_.push_back(T(c + idx, d + jdx, this->K(6 + idx, 9 + jdx)));

			K_.push_back(T(d + idx, a + jdx, this->K(9 + idx, 0 + jdx)));
			K_.push_back(T(d + idx, b + jdx, this->K(9 + idx, 3 + jdx)));
			K_.push_back(T(d + idx, c + jdx, this->K(9 + idx, 6 + jdx)));
			K_.push_back(T(d + idx, d + jdx, this->K(9 + idx, 9 + jdx)));
		}
	}
}

void Tetrahedron::assembleGlobalStiffnessMatrixDense(MatrixXd &K_global) {
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


void Tetrahedron::computeForceDifferentials(MatrixXd &K_global) {
	computeForceDifferentials();
	assembleGlobalStiffnessMatrixDense(K_global);
}


void Tetrahedron::computeForceDifferentialsSparse(VectorXd dx, int row, int col, vector<T> &K_) {
	this->F = computeDeformationGradient();
	this->dF = computeDeformationGradientDifferential(dx);
	this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);	
	this->dH.noalias() = -W * dP * this->BmT;

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		//Vector3d temp = this->dH.col(i);
		for (int j = 0; j < 3; ++j) {
			K_.push_back(T(row + 3 * m_nodes[i]->i + j, col, dH(j, i)));
			K_.push_back(T(row + 3 * m_nodes[3]->i + j, col, -dH(j, i)));
		}
	}
}

Matrix3d Tetrahedron::computePKStress(Matrix3d F, double mu, double lambda) {

	Matrix3d P = Matrix3d::Zero();
	//if (m_isSVD) {
	//	m_material = CO_ROTATED;
	//}
	//else {
	//	//m_material = STVK;
	//}

	switch (m_material)
	{
	case LINEAR:
	{
		Matrix3d I = Matrix3d::Identity();
		Matrix3d E = Matrix3d::Zero();

		E.noalias() = 0.5 * (F + F.transpose()) - I;
		//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		P.noalias() = 2.0 * mu * E + lambda * E.trace() * I;
		break;
	}

	case NEO_HOOKEAN:
	{
		Matrix3d FIT, FT;

		//double I1 = (F.transpose() * F).trace();
		//double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		FT = F.transpose();
		FIT = F.inverse().transpose();

		double I3 = (FT * F).determinant();
		double J = sqrt(I3);
		//psi = 1.0 / 2.0 * mu *(I1 - 3.0) - mu * log(J) + 1.0 / 2.0 * lambda * log(J)*log(J);
		P.noalias() = mu * (F - FIT) + lambda * log(J)*FIT;
		break;
	}

	case STVK:
	{
		Matrix3d I = Matrix3d::Identity();
		Matrix3d E = Matrix3d::Zero();

		//E.noalias() = 0.5 * (F.transpose() * F - I);
		E.noalias() = 0.5 * (this->FTF - I);
		//psi = mu * E.norm()*E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		P.noalias() = F * (2.0 * mu * E + lambda * E.trace() * I);
		break;
	}

	case CO_ROTATED:
	{
		// Polar decomposition
		Matrix3d R = gs3(F);
		Matrix3d I = Matrix3d::Identity();


		//A.noalias() = F.adjoint() * F;
		//SelfAdjointEigenSolver<Matrix3d> es(A);
		//S.noalias() = es.operatorSqrt();
		//R.noalias() = F * S.inverse();

		//E = S - I;
		//psi = mu * E.norm() * E.norm() + 1.0 / 2.0 * lambda * E.trace() * E.trace();
		//P = R * (2.0 * mu * E + lambda * E.trace() * I);
		P.noalias() = 2.0 * mu * (F - R) + lambda * (R.transpose() * F - I).trace() * R;
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
	
	Matrix3d dP = Matrix3d::Zero();

	switch (m_material) {
	case CO_ROTATED:
	{
		Matrix3d R = gs3(F);
		//A.noalias() = F.adjoint() * F;
		//SelfAdjointEigenSolver<Matrix3d> es(A);
		//S = es.operatorSqrt();
		//R.noalias() = F * S.inverse();

		//E = S - I3;
		//P = 2.0 * mu *(F - R) + lambda * (R.transpose()*F - I3).trace() * R;
		dP.noalias() = 2.0 * mu * dF + lambda * (R.transpose() * dF).trace() * R;
		break;
	}

	case STVK:
	{
		Matrix3d E = Matrix3d::Zero();

		Matrix3d dE = Matrix3d::Zero();
		Matrix3d I3 = Matrix3d::Identity();
		E.noalias() = 0.5 * (this->FTF - I3);
		//E.noalias() = 1.0 / 2.0 * (F.transpose() * F - I3);
		dE.noalias() = 0.5 * (dF.transpose() * F + F.transpose() * dF);
		//P = F * (2.0 * mu * E + lambda * E.trace() * I3);
		dP.noalias() = dF * (2.0 * mu * E + lambda * E.trace() * I3) + F * (2.0 * mu * dE + lambda * dE.trace() * I3);
		break;
	}

	case NEO_HOOKEAN:
	{	
		Matrix3d FIT, FT;

		/*double I1 = (F.transpose() * F).trace();
		double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		double I3 = (F.transpose() * F).determinant();
		double J = sqrt(I3);
		psi = 1.0 / 2.0 * mu *(I1 - 3.0) - mu * log(J) + 1.0 / 2.0 * lambda * log(J)*log(J);
		P = mu * (F - F.inverse().transpose()) + lambda * log(J)*(F.inverse().transpose());
*/
		//double I1 = (F.norm()) * (F.norm());
		//double I2 = ((F.transpose() * F) *  (F.transpose() * F)).trace();
		FT = F.transpose();
		FIT = F.inverse().transpose();
		//double I3 = (F.transpose() * F).determinant();
		double I3 = (FT * F).determinant();
		double J = sqrt(I3);
		//P = mu * (F - FIT) + lambda * log(J) * FIT;
		dP.noalias() = mu * dF + (mu - lambda * log(J)) * FIT * (dF.transpose()) * FIT + lambda * ((F.inverse() * dF)).trace() * FIT;
		//P = mu * (F - (F.inverse().transpose())) + lambda * log(J) * (F.inverse().transpose());
		//dP = mu * dF + (mu - lambda * log(J)) * (F.inverse().transpose()) * (dF.transpose()) * (F.inverse().transpose()) + lambda * ((F.inverse() * dF)).trace() * (F.inverse().transpose());
		break;
	}
	case LINEAR:
	{
		Matrix3d E = Matrix3d::Zero();
		Matrix3d dE = Matrix3d::Zero();
		Matrix3d I3 = Matrix3d::Identity();
		E.noalias() = 1.0 / 2.0 * (F + F.transpose()) - I3;
		dE.noalias() = 1.0 / 2.0 * (dF + dF.transpose());
		//P = 2.0 * mu * E + lambda * E.trace() * I3;
		dP.noalias() = 2.0 * mu * dE + lambda * dE.trace() * I3;
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

bool Tetrahedron::checkSameSide(const shared_ptr<Node> &v0, const shared_ptr<Node> &v1, const shared_ptr<Node> &v2, const shared_ptr<Node> &v3, const shared_ptr<Node> &p) {
	Vector3d normal = (v1->x - v0->x).cross(v2->x - v0->x);
	double dotV3 = normal.dot(v3->x - v0->x);
	double dotP = normal.dot(p->x - v0->x);

	if (dotV3 * dotP >= 0) {
		return true;
	}
	else {
		return false;
	}
}

bool Tetrahedron::checkPointInside(const shared_ptr<Node>& p) {
	return checkSameSide(m_nodes[0], m_nodes[1], m_nodes[2], m_nodes[3], p) &&
		checkSameSide(m_nodes[1], m_nodes[2], m_nodes[3], m_nodes[0], p) &&
		checkSameSide(m_nodes[2], m_nodes[3], m_nodes[0], m_nodes[1], p) &&
		checkSameSide(m_nodes[3], m_nodes[0], m_nodes[1], m_nodes[2], p);
}

Vector4d Tetrahedron::computeBarycentricWeightAndSave(const shared_ptr<Node>& p) {
	Vector4d weight;
	Vector3d vap = p->x - m_nodes[0]->x;
	Vector3d vbp = p->x - m_nodes[1]->x;
	Vector3d vcp = p->x - m_nodes[2]->x;
	Vector3d vdp = p->x - m_nodes[3]->x;

	Vector3d vab = m_nodes[1]->x - m_nodes[0]->x;
	Vector3d vac = m_nodes[2]->x - m_nodes[0]->x;
	Vector3d vad = m_nodes[3]->x - m_nodes[0]->x;

	Vector3d vbc = m_nodes[2]->x - m_nodes[1]->x;
	Vector3d vbd = m_nodes[3]->x - m_nodes[1]->x;

	double va = ScalarTripleProduct(vbp, vbd, vbc);
	double vb = ScalarTripleProduct(vap, vac, vad);
	double vc = ScalarTripleProduct(vap, vad, vab);
	double vd = ScalarTripleProduct(vap, vab, vac);

	double v = 1.0 / ScalarTripleProduct(vab, vac, vad);
	weight << va * v, vb * v, vc * v, vd * v;
	m_barycentric_weights.push_back(weight);
	return weight;
}

double Tetrahedron::ScalarTripleProduct(const Vector3d &a, const Vector3d &b, const Vector3d &c) {
	return a.dot(b.cross(c));
}

Vector3d Tetrahedron::computePositionByBarycentricWeight(const Vector4d &weight) {
	Vector3d pos;
	pos.setZero();
	for (int i = 0; i < 4; i++) {
		pos += m_nodes[i]->x * weight(i);
	}

	return pos;
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
	/*if ((int)m_enclosed_points.size() > 0) {
		prog->bind();
		for (int i = 0; i < (int)m_enclosed_points.size(); i++) {
			auto node = m_enclosed_points[i];
			glUniform3fv(prog->getUniform("kd"), 1, m_enclosed_color.data());
			node->draw(MV, prog);
		}
		prog->unbind();
	}*/

}
