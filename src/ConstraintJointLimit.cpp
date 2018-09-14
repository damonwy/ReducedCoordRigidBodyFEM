#include "ConstraintJointLimit.h"

#include "Joint.h"
#include <json.hpp>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

ConstraintJointLimit::ConstraintJointLimit() {


}

ConstraintJointLimit::ConstraintJointLimit(std::shared_ptr<Joint> joint):
Constraint(0, 0, 0, 1),
m_joint(joint)
{
	m_name =  m_joint->getName() + "-LIMIT";

}

void ConstraintJointLimit::computeJacIneqR_(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr) {
	int row = idxIR;
	int col = m_joint->idxR;
	idxQ = col;
	nQ = m_joint->m_ndof;
	cout << m_joint->m_q(0) << endl;

	if (m_joint->m_q(0) <= m_ql) {
		Cr.block(row, col, nconIR, m_joint->m_ndof).setOnes();
		Cr.block(row, col, nconIR, m_joint->m_ndof) *= -1;
		activeR = true;
	}
	else if (m_joint->m_q(0) >= m_qu) {
		Cr.block(row, col, nconIR, m_joint->m_ndof).setOnes();
		activeR = true;
	}
	else {
		activeR = false;
	}
}