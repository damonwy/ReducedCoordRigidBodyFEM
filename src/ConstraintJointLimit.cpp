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