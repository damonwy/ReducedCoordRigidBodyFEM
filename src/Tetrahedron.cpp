#include "Tetrahedron.h"
#include "Node.h"
#include <iostream>
#include "svd3.h"


#include "MatrixStack.h"
#include "Program.h"

using namespace Eigen;
using namespace std;

#define Fthreshold 0.7

Tetrahedron::Tetrahedron()
{

}

Tetrahedron::Tetrahedron(double young, double poisson, double density, Material material, const vector<shared_ptr<Node>> &nodes) :
	m_young(young), m_poisson(poisson), m_density(density), m_material(material), m_nodes(nodes)
{
	m_mu = m_young / (2.0 * (1.0 + m_poisson));
	m_lambda = m_young * m_poisson / ((1.0 + m_poisson) * (1.0 - 2.0 * m_poisson));

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->Dm.col(i) = m_nodes[i]->x0 - m_nodes[3]->x0;
	}
	this->Bm = this->Dm.inverse();
	this->W = abs(1.0 / 6.0 * Dm.determinant());
	m_mass = this->W * this->m_density;

	// Distribute 1/4 mass to each node
	for (int i = 0; i < (int)nodes.size(); i++) {
		nodes[i]->m += this->m_mass * 0.25;
	}

	computeAreaWeightedVertexNormals();

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
}



Eigen::Matrix3d Tetrahedron::computeDeformationGradient() {

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	return this->F;
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
	/*if (m_isInvertible) {
		isInverted();
	}
	*/

	this->F = computeDeformationGradient();
	// The deformation gradient is available in this->F


	//if (isInvert && m_isInvertible) {
	//	this->F = this->Fhat; // Use the new F
	//}

	this->P = computePKStress(F, m_mu, m_lambda);
	this->H = -W * P * (Bm.transpose());

	/*if (isInvert && m_isInvertible) {
		this->H = -W * U * P * V.transpose() * (Bm.transpose());
	}*/

	//for (int i = 0; i < 3; i++) {
	//	double force = this->H.col(i).norm();
	//	if (force > 1.0e2) {
	//		this->H *= 1.0e2 / force;
	//	}
	//}
	//

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

VectorXd Tetrahedron::computeInvertibleElasticForces(VectorXd f) {
	bool print = false;

	this->F = computeDeformationGradient();
	if (abs(this->F(0, 1) > 0.3)) {
		print = true;
	}
	if (print) {
		cout << "F" << this->F << endl;
			cout << "P suppose: " << endl<<computePKStress(F, m_mu, m_lambda);
	}
	
	// The deformation gradient is available in this->F
	if (this->F.determinant() < 0.0) {
		isInvert = true;
	}
	else {
		isInvert = false;
	}

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

	}

	Matrix3d ttemp = computePKStress(F, m_mu, m_lambda);
	this->H = -W * ttemp * (Bm.transpose());
	if (print) {
		cout << "H" << this->H << endl;
		cout << "PNm" << this->P * this->Nm.block<3, 3>(0, 0) << endl;
	}
	
	// Computes the nodal forces by G=PBm=PNm

	for (int i = 0; i < (int)m_nodes.size()-1; i++) {
		int rowi = m_nodes[i]->idxM;
		f.segment<3>(rowi) += this->P * this->Nm.col(i);
		int row3 = m_nodes[3]->idxM;
		f.segment<3>(row3) -= this->P * this->Nm.col(i);
	}

	return f;
}

void Tetrahedron::Compute_dGdF(Vec3d * b0, Vec3d * b1, Vec3d * b2,
	double dPdF[81], double dGdF[81])
{
	//Both G and F are 3x3 matrices, so dGdF has 81 entries
	memset(dGdF, 0, sizeof(double) * 81);

	/*
	| ga_x gb_x gc_x |   | 0 1 2 |
	if G = | ga_y gb_y gc_y | = | 3 4 5 |
	| ga_z gb_z gc_z |   | 6 7 8 |
	where ga, gb, gc are the nodal forces at vertex a,b,c
	| ba_0 bb_0 bc_0 |   | 0 1 2 |
	and B = | ba_1 bb_1 bc_1 | = | 3 4 5 |
	| ba_2 bb_2 bc_2 |   | 6 7 8 |
	| dga_x/dF_00 dga_x/dF_01 dga_x/dF_02 dga_x/dF_10 ... dga_x/dF_22 |
	| dga_y/dF_00 dga_y/dF_01 dga_y/dF_02 dga_y/dF_10 ... dga_y/dF_22 |
	| dga_z/dF_00 dga_z/dF_01 dga_z/dF_02 dga_z/dF_10 ... dga_z/dF_22 |
	dGdF = | dgb_x/dF_00 dgb_x/dF_01 dgb_x/dF_02 dgb_x/dF_10 ... dgb_x/dF_22 |
	|                                 ...                             |
	| dgc_z/dF_00 dgc_z/dF_01 dgc_z/dF_02 dgc_z/dF_10 ... dgc_z/dF_22 |
	*/

	Vec3d * bVec[3] = { b0, b1, b2 };

	//dga_x/dF, dga_y/dF, dga_z/dF
	//dgb_x/dF, dgb_y/dF, dgb_z/dF
	//dgc_x/dF, dgc_y/dF, dgc_z/dF
	memset(dGdF, 0, sizeof(double) * 81);
	for (int abc = 0; abc<3; abc++)
		for (int i = 0; i<3; i++)
			for (int column = 0; column<9; column++)
				for (int k = 0; k<3; k++)
					dGdF[27 * abc + 9 * i + column] += dPdF[(3 * i + k) * 9 + column] * (*(bVec[abc]))[k];

	/*
	printf("---- printing dGdF ----\n");
	for (int i=0; i<9; i++)
	{
	for (int j=0; j<9; j++)
	printf("%G ", dGdF[i*9+j]);
	printf("\n");
	}
	*/
}

void Tetrahedron::ComputeTetK(int el, double K[144], int clamped)
{
	/*
	dP/dF is a column major matrix, but is stored as a 1D vector

	| dP_11/dF_11  dP_11/dF_12  dP_11/dF_13  dP_11/dF_21 ... dP_11/dF_33 |
	| dP_12/dF_11  dP_12/dF_12  dP_12/dF_13  dP_12/dF_21 ... dP_12/dF_33 |
	|                              ...                                   |
	| dP_33/dF_11  dP_33/dF_12  dP_33/dF_13  dP_33/dF_21 ... dP_33/dF_33 |
	*/
	double dPdF[81]; //in 9x9 matrix format
	double dGdF[81]; //in 9x9 matrix format

	Compute_dPdF(el, dPdF, clamped);
	Vec3d areaWeightedVertexNormals0(Nm.col(0)(0), Nm.col(0)(1), Nm.col(0)(2));
	Vec3d areaWeightedVertexNormals1(Nm.col(1)(0), Nm.col(1)(1), Nm.col(1)(2));
	Vec3d areaWeightedVertexNormals2(Nm.col(2)(0), Nm.col(2)(1), Nm.col(2)(2));

	Compute_dGdF(&(areaWeightedVertexNormals0), &(areaWeightedVertexNormals1),
		&(areaWeightedVertexNormals2), dPdF, dGdF);
	//dF_dU was already computed by the constructor before calling this function
	double * dFdU; //= &dFdUs[108 * el];

	// K is stored column-major (however, it doesn't matter because K is symmetric)
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

	//The last three columns are combinations of the first nine columns.
	//The reason is that the nodal force of the 4th vertex equals to 
	//the minus of the sum of the 1st, 2nd, and 3rd vertices (see p3 
	//section 4 of [Irving 04]
	for (int row = 0; row < 12; row++)
	{
		//10th column
		K[12 * row + 9] = -K[12 * row + 0] - K[12 * row + 3] - K[12 * row + 6];
		//11th column
		K[12 * row + 10] = -K[12 * row + 1] - K[12 * row + 4] - K[12 * row + 7];
		//12th column
		K[12 * row + 11] = -K[12 * row + 2] - K[12 * row + 5] - K[12 * row + 8];
	}
}

void Compute_dFdU()
{
	
	//double * dFdU = &dFdUs[108 * el];
	//Mat3d & dmInv = dmInverses[el];

	for (int index = 0; index<108; index++)
	{
		int n = index % 3;
		int m = (int)(index / 3) % 4;
		int j = (int)(index / 12) % 3;
		int i = (int)(index / 36) % 3;
		double result = 0.0;
		//for (int k = 0; k<3; k++)
		//result += dDSdU[tensor9x12Index(i, k, m, n)] * dmInv[k][j];
		//dFdU[tensor9x12Index(i, j, m, n)] = result;
	}
}

void ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient, double mu, double lambda) // invariants and gradient are 3-vectors
{
	double IIIC = invariants[2];

	gradient[0] = 0.5 * mu;
	gradient[1] = 0.0;
	gradient[2] = (-0.5 * mu + 0.25 * lambda * log(IIIC)) / IIIC;

}

void ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian, double mu, double lambda) // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
{
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
	hessian[5] = (0.25 * lambda + 0.5 * mu - 0.25 * lambda * log(IIIC)) / (IIIC * IIIC);
}

// gradient of P with respect to F(9x9 matrix, row - major)
// see [Teran 05]
void Tetrahedron::Compute_dPdF(int el, double dPdF[81], int clamped)
{
	double sigma[3] = { Fhats[0], Fhats[1], Fhats[2] };

	double sigma1square = sigma[0] * sigma[0];
	double sigma2square = sigma[1] * sigma[1];
	double sigma3square = sigma[2] * sigma[2];

	double invariants[3];
	invariants[0] = sigma1square + sigma2square + sigma3square;
	invariants[1] = (sigma1square * sigma1square +
		sigma2square * sigma2square +
		sigma3square * sigma3square);
	invariants[2] = sigma1square * sigma2square * sigma3square;

	//double E[3];
	//E[0] = 0.5 * (Fhats[el][0] * Fhats[el][0] - 1);
	//E[1] = 0.5 * (Fhats[el][1] * Fhats[el][1] - 1);
	//E[2] = 0.5 * (Fhats[el][2] * Fhats[el][2] - 1);

	double gradient[3];
	ComputeEnergyGradient(el, invariants, gradient, m_mu, m_lambda);

	/*
	in order (11,12,13,22,23,33)
	| 11 12 13 |   | 0 1 2 |
	| 21 22 23 | = | 1 3 4 |
	| 31 32 33 |   | 2 4 5 |
	*/
	double hessian[6];
	ComputeEnergyHessian(el, invariants, hessian, m_mu, m_lambda);

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

	/*if (enforceSPD)
	{
		FixPositiveIndefiniteness(x1111, x2211, x3311, x2222, x3322, x3333);
		FixPositiveIndefiniteness(x2121, x2112);
		FixPositiveIndefiniteness(x3131, x3113);
		FixPositiveIndefiniteness(x3232, x3223);
	}*/

	double dPdF_atFhat[81];
	memset(dPdF_atFhat, 0, sizeof(double) * 81);
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

	
	//Mat3d UT = trans(Us[el]); // trans(*U);
	//Mat3d VT = trans(Vs[el]); // trans(*V);

							  /*
							  U->print();
							  V->print();
							  UT.print();
							  VT.print();
							  */

	double eiejVector[9];
	memset(eiejVector, 0, sizeof(double) * 9);
	memset(dPdF, 0, sizeof(double) * 81);
	Matrix3d eiejMatrix;
	eiejMatrix.setZero();

	double ut[9];
	memset(ut, 0, sizeof(double) * 9);
	for (int i = 0; i < 9; ++i) {
		ut[i] = U.transpose()(i);
	}

	Mat3d UT(ut);

	double v[9];
	memset(v, 0, sizeof(double) * 9);
	for (int i = 0; i < 9; ++i) {
		v[i] = V(i);
	}

	Mat3d V(v);

	for (int column = 0; column<9; column++)
	{
		eiejVector[column] = 1.0;


		Mat3d ei_ej(eiejVector);
		Mat3d ut_eiej_v = UT * ei_ej * V;
		double ut_eiej_v_TeranVector[9]; //in Teran order
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[0]] = ut_eiej_v[0][0];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[1]] = ut_eiej_v[0][1];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[2]] = ut_eiej_v[0][2];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[3]] = ut_eiej_v[1][0];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[4]] = ut_eiej_v[1][1];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[5]] = ut_eiej_v[1][2];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[6]] = ut_eiej_v[2][0];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[7]] = ut_eiej_v[2][1];
		ut_eiej_v_TeranVector[rowMajorMatrixToTeran[8]] = ut_eiej_v[2][2];
		double dPdF_resultVector[9]; // not in Teran order
		for (int innerRow = 0; innerRow<9; innerRow++)
		{
			double tempResult = 0.0;
			for (int innerColumn = 0; innerColumn<9; innerColumn++)
			{
				tempResult += dPdF_atFhat[innerRow * 9 + innerColumn] *
					ut_eiej_v_TeranVector[innerColumn];
			}
			dPdF_resultVector[teranToRowMajorMatrix[innerRow]] = tempResult;
		}
		Mat3d dPdF_resultMatrix(dPdF_resultVector);
		Mat3d u_dpdf_vt = (UT)*dPdF_resultMatrix*V;
		dPdF[column + 0] = u_dpdf_vt[0][0];
		dPdF[column + 9] = u_dpdf_vt[0][1];
		dPdF[column + 18] = u_dpdf_vt[0][2];
		dPdF[column + 27] = u_dpdf_vt[1][0];
		dPdF[column + 36] = u_dpdf_vt[1][1];
		dPdF[column + 45] = u_dpdf_vt[1][2];
		dPdF[column + 54] = u_dpdf_vt[2][0];
		dPdF[column + 63] = u_dpdf_vt[2][1];
		dPdF[column + 72] = u_dpdf_vt[2][2];
		// reset
		eiejVector[column] = 0.0;
	}

	/*
	printf("---- full dPdF ----\n");
	for (int i=0; i<9; i++)
	{
	for (int j=0; j<9; j++)
	printf("%G ", dPdF[i*9+j]);
	printf(";\n");
	}
	*/
}


double Tetrahedron::gammaValue(int i, int j, double sigma[3], double invariants[3], double gradient[3], double hessian[6])
{
	/*
	The hessian is in order (11,12,13,22,23,33)
	| 11 12 13 |   | 0 1 2 |
	| 21 22 23 | = | 1 3 4 |
	| 31 32 33 |   | 2 4 5 |
	*/

	double tempGammaVec1[3];
	tempGammaVec1[0] = 2.0 * sigma[i];
	tempGammaVec1[1] = 4.0 * sigma[i] * sigma[i] * sigma[i];
	tempGammaVec1[2] = 2.0 * invariants[2] / sigma[i];
	double tempGammaVec2[3];
	tempGammaVec2[0] = 2.0 * sigma[j];
	tempGammaVec2[1] = 4.0 * sigma[j] * sigma[j] * sigma[j];
	tempGammaVec2[2] = 2.0 * invariants[2] / sigma[j];
	double productResult[3];
	productResult[0] = (tempGammaVec2[0] * hessian[0] + tempGammaVec2[1] * hessian[1] +
		tempGammaVec2[2] * hessian[2]);
	productResult[1] = (tempGammaVec2[0] * hessian[1] + tempGammaVec2[1] * hessian[3] +
		tempGammaVec2[2] * hessian[4]);
	productResult[2] = (tempGammaVec2[0] * hessian[2] + tempGammaVec2[1] * hessian[4] +
		tempGammaVec2[2] * hessian[5]);
	return (tempGammaVec1[0] * productResult[0] + tempGammaVec1[1] * productResult[1] +
		tempGammaVec1[2] * productResult[2] + 4.0 * invariants[2] * gradient[2] / (sigma[i] * sigma[j]));
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


void Tetrahedron::computeInvertibleForceDifferentials(VectorXd dx, VectorXd &df) {
	//this->F = computeDeformationGradient();

	/*if (isInvert && m_isInvertible) {
	this->F = this->Fhat;
	}*/

	//// clamp if below the principal stretch threshold
	//int clamped = 0;
	//for (int i = 0; i < 3; i++)
	//{
	//	if (abs(this->F(i, i)) < Fthreshold)
	//	{
	//		//dropBelowThreshold = true;
	//		cout << this->F(i, i) << endl;
	//		if (this->F(i, i) < 0.0) {
	//			this->F(i, i) = -Fthreshold;
	//		}
	//		else {
	//			this->F(i, i) = Fthreshold;
	//		}

	//		clamped |= (1 << i);
	//	}
	//}

	//clamped = 0; // disable clamping

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->dDs.col(i) = dx.segment<3>(3 * m_nodes[i]->i) - dx.segment<3>(3 * m_nodes[3]->i);
	}

	this->dF = dDs * Bm;
	this->dPhat = computePKStressDerivative(this->Fhat, dF, m_mu, m_lambda);
	this->dP = this->U * this->dPhat * this->V.transpose();
	//this->dH = -W * dP * (Bm.transpose());
	//cout << "dP_INVERT" << this->dP << endl;
	//cout << "df" << this->dP * Nm.block<3, 3>(0, 0) << endl;

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		//df.segment<3>(3 * m_nodes[i]->i) += this->dH.col(i);
		//df.segment<3>(3 * m_nodes[3]->i) -= this->dH.col(i);

		//df.segment<3>(3 * m_nodes[i]->i) += this->dP * this->Nm.col(i);
		//df.segment<3>(3 * m_nodes[3]->i) -= this->dP * this->Nm.col(i);
	}

}

void Tetrahedron::computeInvertibleForceDifferentialsSparse(Eigen::VectorXd dx, int row, int col, std::vector<T> &K_) {

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->dDs.col(i) = dx.segment<3>(3 * m_nodes[i]->i) - dx.segment<3>(3 * m_nodes[3]->i);
	}

	this->dF = dDs * Bm;
	this->dPhat = computePKStressDerivative(this->F, dF, m_mu, m_lambda);
	this->dP = this->U * this->dPhat * this->V.transpose();

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		Vector3d temp = this->dP * this->Nm.col(i);

		for (int j = 0; j < 3; ++j) {
			K_.push_back(T(row + 3 * m_nodes[i]->i + j, col, temp(j)));
			K_.push_back(T(row + 3 * m_nodes[3]->i + j, col, -temp(j)));

		}	
	}
}

void Tetrahedron::computeForceDifferentials(VectorXd dx, VectorXd& df) {
	this->F = computeDeformationGradient();

	/*if (isInvert && m_isInvertible) {
	this->F = this->Fhat;
	}*/

	//// clamp if below the principal stretch threshold
	//int clamped = 0;
	//for (int i = 0; i < 3; i++)
	//{
	//	if (abs(this->F(i, i)) < Fthreshold)
	//	{
	//		//dropBelowThreshold = true;
	//		cout << this->F(i, i) << endl;
	//		if (this->F(i, i) < 0.0) {
	//			this->F(i, i) = -Fthreshold;
	//		}
	//		else {
	//			this->F(i, i) = Fthreshold;
	//		}

	//		clamped |= (1 << i);
	//	}
	//}

	//clamped = 0; // disable clamping

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->dDs.col(i) = dx.segment<3>(3 * m_nodes[i]->i) - dx.segment<3>(3 * m_nodes[3]->i);
	}

	this->dF = dDs * Bm;
	this->dP = computePKStressDerivative(F, dF, m_mu, m_lambda);
	Matrix3d temp = computePKStressDerivative(Fhat, dF, m_mu, m_lambda);
	//cout << "dP" << dP << endl;
	//cout << "dP_hat" << temp << endl;
	this->dH = -W * dP * (Bm.transpose());
	//cout << "dH " <<this->dH << endl;
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
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

void Tetrahedron::computeForceDifferentialsSparse(VectorXd dx, int row, int col, std::vector<T> &K_) {
	this->F = computeDeformationGradient();

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->dDs.col(i) = dx.segment<3>(3 * m_nodes[i]->i) - dx.segment<3>(3 * m_nodes[3]->i);
	}

	this->dF = dDs * Bm;
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

Matrix3d Tetrahedron::computeInvertiblePKStressDerivative(Matrix3d F, Matrix3d dF, double mu, double lambda) {


	return F;

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
		// modify hessian to compute correct values if in the inversion handling regime
		if (clamped & 1) // first lambda was clamped (in inversion handling)
		{
			dP(0, 0) = 0.0;
			dP(1, 0) = 0.0;
			dP(0, 1) = 0.0;
			//hessian[0] = hessian[1] = hessian[2] = 0.0;
		}

		if (clamped & 2) // second lambda was clamped (in inversion handling)
		{
			dP(2, 1) = 0.0;
			dP(1, 2) = 0.0;
			dP(1, 0) = 0.0;
			dP(0, 1) = 0.0;
			dP(1, 1) = 0.0;
			//hessian[1] = hessian[3] = hessian[4] = 0.0;
		}

		if (clamped & 4) // third lambda was clamped (in inversion handling)
		{
			dP(0, 0) = 0.0;
			dP(1, 0) = 0.0;
			dP(0, 1) = 0.0;
			dP(2, 1) = 0.0;
			dP(1, 2) = 0.0;
			dP(2, 2) = 0.0;
			//hessian[0] = hessian[1] = hessian[2] = hessian[4] = hessian[5] = 0.0;
		}

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

	int modifiedSVD = 1;
	Vector3d Fhat_vec;
	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	if (this->F.determinant() < 0.0) { // some threshold todo
		SVD(this->F, this->U, Fhat_vec, this->V, 1e-8, modifiedSVD);
		
		// clamp if below the principal stretch threshold
		int clamped = 0;
		for (int i = 0; i < 3; i++)
		{
			if (abs(Fhat_vec(i)) < Fthreshold)
			{
				if (Fhat_vec(i) < 0.0) {
					Fhat_vec(i) = -Fthreshold;
				}
				else {
					Fhat_vec(i) = Fthreshold;
				}
				
				//Fhat_vec(i) = 0.5;
				clamped |= (1 << i);

			}
		}
		
		clamped = 0; // disable clamping

		this->Fhat = Fhat_vec.asDiagonal();
		
		//diagDeformationGradient(this->F);
		isInvert = true;

		m_nodes[0]->r = 0.05;
		m_nodes[1]->r = 0.05;
		m_nodes[2]->r = 0.05;
		m_nodes[3]->r = 0.05;

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




double Tetrahedron::computeEnergy() {
	//isInverted();

	for (int i = 0; i < (int)m_nodes.size() - 1; i++) {
		this->Ds.col(i) = m_nodes[i]->x - m_nodes[3]->x;
	}

	this->F = Ds * Bm;
	/*if (isInvert && m_isInvertible) {
		this->F = this->Fhat;
	}*/

	this->P = computePKStress(F, m_mu, m_lambda);
	this->m_energy = W * psi;
	return this->m_energy;
}


void Tetrahedron::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	if (isInvert) {
		prog->bind();
		for (int i = 0; i < 4; i++) {
			auto node = m_nodes[i];
			node->draw(MV, prog);
		}
		prog->unbind();
	}

}

void Tetrahedron::compute_dPdF() {
	Vector3d sigma, invariants, gradient;
	sigma << Fhat(0, 0), Fhat(1, 1), Fhat(2, 2);

	double sigma1square = sigma[0] * sigma[0];
	double sigma2square = sigma[1] * sigma[1];
	double sigma3square = sigma[2] * sigma[2];

	invariants[0] = sigma1square + sigma2square + sigma3square;
	invariants[1] = (sigma1square * sigma1square +
		sigma2square * sigma2square +
		sigma3square * sigma3square);
	invariants[2] = sigma1square * sigma2square * sigma3square;
	
	gradient << 0.5 * m_mu, 0.0, (-0.5 * m_mu + 0.25 * m_lambda * log(invariants[2])) / invariants[2];

	Vector6d hessian;
	hessian.setZero();
	hessian(5) = (0.25 * m_lambda + 0.5 * m_mu - 0.25 * m_lambda * log(invariants[2])) / (invariants[2] * invariants[2]);

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



}

void Tetrahedron::compute_dGdF() {
	Vector3d b0 = this->Nm.col(0);
	Vector3d b1 = this->Nm.col(1);
	Vector3d b2 = this->Nm.col(2);



}

void Tetrahedron::compute_dFdU() {



}



Tetrahedron:: ~Tetrahedron() {

}