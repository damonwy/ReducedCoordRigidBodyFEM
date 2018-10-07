#include "World.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Joint.h"
#include "JointRevolute.h"
#include "Body.h"
#include "SoftBody.h"
#include "MatrixStack.h"
#include "Program.h"
#include "SE3.h"
#include "JsonEigen.h"
#include "ConstraintJointLimit.h"
#include "ConstraintNull.h"
#include "ConstraintLoop.h"
#include "ConstraintAttachSpring.h"
#include "ConstraintAttachSoftBody.h"
#include "Spring.h"
#include "SpringSerial.h"
#include "SpringNull.h"
#include "JointNull.h"
#include "SoftBodyNull.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

World::World():
nr(0), nm(0), nem(0), ner(0), ne(0), nim(0), nir(0), m_nbodies(0), m_njoints(0), m_nsprings(0), m_constraints(0), m_countS(0), m_countCM(0),
m_nsoftbodies(0)
{

}

World::World(WorldType type):
m_type(type),
nr(0), nm(0), nem(0), ner(0), ne(0), nim(0), nir(0), m_nbodies(0), m_njoints(0), m_nsprings(0), m_nconstraints(0), m_countS(0), m_countCM(0),
m_nsoftbodies(0)
{
}

World::~World() {
}

void World::load(const std::string &RESOURCE_DIR) {

	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	json js;
	i >> js;
	i.close();

	double density;
	Eigen::Vector3d sides;
	Matrix4d E;
	Vector3d p;

	switch (m_type)
	{
	case SERIAL_CHAIN: 
		{	
			m_h = 1.0e-2;
			density = 1.0;
			m_grav << 0.0, -98, 0.0;
			Eigen::from_json(js["sides"], sides);
			//m_nbodies = 5;
			//m_njoints = 5;
			m_Hexpected = 10000; // todo
			m_tspan << 0.0, 5.0;
			
			// Inits rigid bodies
			for (int i = 0; i < 1; i++) {

				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

				// Inits joints
				if (i == 0) {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				}
				else {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, m_joints[i-1]);
				}
			}
			break;
		}
	case DIFF_REVOLUTE_AXES:
		break;
	case BRANCHING:
		{
			m_h = 1.0e-2;
			m_tspan << 0.0, 50.0;
			density = 1.0;
			m_grav << 0.0, -98, 0.0;
			Eigen::from_json(js["sides"], sides);
			Vector3d sides_0;
			sides_0 << 1.0, 10.0, 1.0;
			Vector3d sides_1;
			sides_1 << 20.0, 1.0, 1.0;

			auto b0 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b1 = addBody(density, sides_1, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box20_1_1.obj");
			auto b2 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b3 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b4 = addBody(density, sides, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
			auto b5 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b6 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b7 = addBody(density, sides, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
			auto b8 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b9 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");

			auto j0 = addJointRevolute(b0, Vector3d::UnitX(), Vector3d(0.0, 15.0, 0.0), Matrix3d::Identity(), 0.0);
			auto j1 = addJointRevolute(b1, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, j0);
			auto j2 = addJointRevolute(b2, Vector3d::UnitX(), Vector3d(-10.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j1);
			auto j3 = addJointRevolute(b3, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j1);
			auto j4 = addJointRevolute(b4, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j2);
			auto j5 = addJointRevolute(b5, Vector3d::UnitX(), Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j4);
			auto j6 = addJointRevolute(b6, Vector3d::UnitY(), Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j4);
			auto j7 = addJointRevolute(b7, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j3);
			auto j8 = addJointRevolute(b8, Vector3d::UnitX(), Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j7);
			auto j9 = addJointRevolute(b9, Vector3d::UnitY(), Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, j7);
		}
		break;
	case SHPERICAL_JOINT:
		break;
	case LOOP:
		{
			m_h = 1.0e-2;
			m_tspan << 0.0, 50.0;
			density = 1.0;
			m_grav << 0.0, -98, 0.0;
			Eigen::from_json(js["sides"], sides);
			Vector3d sides_0;
			sides_0 << 1.0, 10.0, 1.0;
			Vector3d sides_1;
			sides_1 << 20.0, 1.0, 1.0;

			auto b0 = addBody(density, sides_1, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box20_1_1.obj");
			auto b1 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b2 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
			auto b3 = addBody(density, sides_1, Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box20_1_1.obj");
			auto b4 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");

			auto j0 = addJointRevolute(b0, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
			auto j1 = addJointRevolute(b1, Vector3d::UnitZ(), Vector3d(-10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, j0);
			auto j2 = addJointRevolute(b2, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, j0);
			auto j3 = addJointRevolute(b3, Vector3d::UnitZ(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, j1);
			auto j4 = addJointRevolute(b4, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, j3);
			j4->m_qdot(0) = 5.0;

			auto constraint = make_shared<ConstraintLoop>(b2, b3);
			m_constraints.push_back(constraint);
			constraint->setPositions(Vector3d(0.0, -5.0, 0.0), Vector3d(10.0, 0.0, 0.0));
			m_nconstraints++;
	
		}
		break;
	case JOINT_TORQUE:
		break;
	case JOINT_LIMITS:
		{
			m_h = 1.0e-2;
			m_tspan << 0.0, 50.0;
			density = 1.0;
			m_grav << 0.0, -98, 0.0;
			Eigen::from_json(js["sides"], sides);

			for (int i = 0; i < 6; i++) {
				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

				// Inits joints
				if (i == 0) {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				}
				else {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, m_joints[i - 1]);
				}

				// Init constraints
				if (i > 0) {
					addConstraintJointLimit(m_joints[i], -M_PI / 4, M_PI / 4);
				}
			}
		

		}

		break;
	case EQUALITY_CONSTRAINED_ANGLES:
		break;
	case EQUALITY_AND_LOOP:
		break;
	case HYBRID_DYNAMICS:
		break;
	case EXTERNAL_WORLD_FORCE:
		break;
	case JOINT_STIFFNESS:
		break;
	case SPRINGS: 
		{	
			m_h = 1.0e-2;
			m_tspan << 0.0, 50.0;
			density = 1.0;
			m_grav << 0.0, -98, 0.0;
			m_stiffness = 5.0e3;
			Eigen::from_json(js["sides"], sides);
			
			for (int i = 0; i < 2; i++) {
				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

				// Inits joints
				if (i == 0) {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				}
				else {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, m_joints[i - 1]);
				}
			}

			// Init springs
			auto spring0 = addSpringSerial(sides(0)*sides(1)*sides(2)*density, 3, nullptr, Vector3d(10.0 * m_nbodies + 10.0, 10.0, 0.0), m_bodies[m_nbodies - 1], Vector3d(5.0, 0.0, 0.0));
			spring0->setStiffness(m_stiffness);
			auto spring1 = addSpringSerial(sides(0)*sides(1)*sides(2)*density, 2, m_bodies[0], Vector3d(0.0, 0.0, 0.0), m_bodies[m_nbodies - 1], Vector3d(0.0, 0.0, 0.0));
			spring1->setStiffness(m_stiffness);
			for (int i = 0; i < m_springs.size(); i++) {
				m_springs[i]->load(RESOURCE_DIR);
			}
			
		}
		break;
	case SOFT_BODIES:
		{
			m_h = 1.0e-2;
			m_tspan << 0.0, 1.0;
			density = 1.0;
			m_grav << 0.0, -98, 0.0;
			Eigen::from_json(js["sides"], sides);
			double young = 1e3;
			double possion = 0.45;

			for (int i = 0; i < 5; i++) {
				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

				// Inits joints
				if (i == 0) {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				}
				else {
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, m_joints[i - 1]);
				}
			}

			/*auto softbody = make_shared<SoftBody>(0.01 * density, young, possion);
			softbody->load(RESOURCE_DIR, "cylinder");
			softbody->transform(Vector3d(10.0, 0.0, 0.0));*/
			
			//m_softbodies.push_back(softbody);
			//m_nsoftbodies++;

			auto softbody = addSoftBody(0.01 * density, young, possion, RESOURCE_DIR, "cylinder");
			softbody->transform(Vector3d(10.0, 0.0, 0.0));
			/*auto softbody1 = make_shared<SoftBody>(0.01 * density, young, possion);
			softbody1->load(RESOURCE_DIR, "cylinder");
			softbody1->transform(Vector3d(20.0, 0.0, 0.0));

			m_softbodies.push_back(softbody1);
			m_nsoftbodies++;*/

		}
		break;
	default:
		break;
	}
	
}

shared_ptr<SoftBody> World::addSoftBody(double density, double young, double possion, const string &RESOURCE_DIR, string file_name) {
	auto softbody = make_shared<SoftBody>(density, young, possion);
	softbody->load(RESOURCE_DIR, file_name);
	//softbody->transform(Vector3d(10.0, 0.0, 0.0));
	m_softbodies.push_back(softbody);
	m_nsoftbodies++;
	return softbody;
}


shared_ptr<Body> World::addBody(double density, Vector3d sides, Vector3d p, Matrix3d R, const string &RESOURCE_DIR, string file_name) {
	auto body = make_shared<Body>(density, sides);
	Matrix4d E = SE3::RpToE(R, p);
	body->setTransform(E);
	body->load(RESOURCE_DIR, file_name);
	m_bodies.push_back(body);
	m_nbodies++;
	return body;
}

shared_ptr<JointRevolute> World::addJointRevolute(shared_ptr<Body> body, Vector3d axis, Vector3d p, Matrix3d R, double q, shared_ptr<Joint> parent) {
	auto joint = make_shared<JointRevolute>(body, axis, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->m_q(0) = q;
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<ConstraintJointLimit> World::addConstraintJointLimit(shared_ptr<Joint> joint, double ql, double qu) {
	auto constraint = make_shared<ConstraintJointLimit>(joint);
	m_constraints.push_back(constraint);
	constraint->setLimits(ql, qu);
	m_nconstraints++;
	return constraint;
}

shared_ptr<SpringSerial> World::addSpringSerial(double mass, int n_points, shared_ptr<Body> body0, Vector3d r0, shared_ptr<Body> body1, Vector3d r1) {

	auto spring = make_shared<SpringSerial>(n_points, m_countS, m_countCM);
	m_springs.push_back(spring);
	spring->setStiffness(m_stiffness);
	spring->setMass(mass); 
	spring->setAttachments(body0, r0, body1, r1);
	m_nsprings++;
	return spring;

}


shared_ptr<ConstraintNull> World::addConstraintNull() {

	auto constraint = make_shared<ConstraintNull>();
	m_nconstraints++;
	m_constraints.push_back(constraint);
	return constraint;

}

shared_ptr<JointNull> World::addJointNull() {
	auto joint = make_shared<JointNull>();
	m_njoints++;
	m_joints.push_back(joint);
	return joint;
}

shared_ptr<SpringNull> World::addSpringNull() {

	auto spring = make_shared<SpringNull>();
	m_nsprings++;
	m_springs.push_back(spring);
	return spring;
}

void World::init() {
	for (int i = 0; i < m_nbodies; i++) {
		
		m_bodies[i]->init(nm);
		if (i < m_nbodies-1) {
			m_bodies[i]->next = m_bodies[i + 1];
		}
	}
	
	nm = 0;
	/*for (int i = 0; i < m_njoints; i++) {
		m_joints[i]->init(nm, nr);
	}*/

	//joint ordering
	// todo

	for (int i = m_njoints - 1; i > -1; i--) {
		m_joints[i]->init(nm, nr);
	}

	for (int i = 0; i < m_njoints; i++) {
		if (i < m_njoints - 1) {
			m_joints[i]->next = m_joints[i + 1];
		}
		if (i > 0) {
			m_joints[i]->prev = m_joints[i - 1];
		}
	}

	if (m_njoints == 0) {
		addJointNull();
	}

	m_joints[0]->update();
	
	for (int i = 0; i < m_nsprings; i++) {
		m_springs[i]->countDofs(nm, nr);
		
		m_springs[i]->init();
		// Create attachment constraints
		auto constraint = make_shared<ConstraintAttachSpring>(m_springs[i]);
		m_constraints.push_back(constraint);
		m_nconstraints++;
		if (i < m_nsprings - 1) {
			m_springs[i]->next = m_springs[i + 1];
		}
	}
	
	if (m_nsprings == 0) {
		addSpringNull();
	}

	if (m_type == SOFT_BODIES) {
		m_softbodies[0]->setAttachments(0, m_bodies[0]);
		m_softbodies[0]->setAttachments(3, m_bodies[0]);
		m_softbodies[0]->setAttachments(6, m_bodies[0]);
		m_softbodies[0]->setAttachments(9, m_bodies[0]);
		m_softbodies[0]->setAttachments(12, m_bodies[0]);
		m_softbodies[0]->setAttachments(19, m_bodies[0]);
		m_softbodies[0]->setAttachments(25, m_bodies[0]);
		//m_softbodies[0]->setAttachments(30, m_bodies[0]);

		m_softbodies[0]->setAttachments(60, m_bodies[1]);
		m_softbodies[0]->setAttachments(63, m_bodies[1]);
		m_softbodies[0]->setAttachments(67, m_bodies[1]);
		m_softbodies[0]->setAttachments(69, m_bodies[1]);
		m_softbodies[0]->setAttachments(72, m_bodies[1]);

	

		/*m_softbodies[1]->setAttachments(0, m_bodies[1]);
		m_softbodies[1]->setAttachments(3, m_bodies[1]);
		m_softbodies[1]->setAttachments(6, m_bodies[1]);
		m_softbodies[1]->setAttachments(9, m_bodies[1]);
		m_softbodies[1]->setAttachments(12, m_bodies[1]);
		m_softbodies[1]->setAttachments(19, m_bodies[1]);

		m_softbodies[1]->setAttachments(60, m_bodies[2]);
		m_softbodies[1]->setAttachments(63, m_bodies[2]);
		m_softbodies[1]->setAttachments(67, m_bodies[2]);
		m_softbodies[1]->setAttachments(69, m_bodies[2]);
		m_softbodies[1]->setAttachments(72, m_bodies[2]);*/

		//m_softbodies[1]->setAttachments(0, m_bodies[1]);
	//	m_softbodies[1]->setAttachments(3, m_bodies[1]);
		//m_softbodies[1]->setAttachments(6, m_bodies[2]);

	}

	for (int i = 0; i < m_nsoftbodies; i++) {
		m_softbodies[i]->countDofs(nm, nr);
		m_softbodies[i]->init();
		// Create attachment constraints
		auto constraint = make_shared<ConstraintAttachSoftBody>(m_softbodies[i]);
		m_constraints.push_back(constraint);
		m_nconstraints++;


		if (i < m_nsoftbodies - 1) {
			m_softbodies[i]->next = m_softbodies[i + 1];
		}
	}

	if (m_nsoftbodies == 0) {
		//todo
		auto softbody = make_shared<SoftBodyNull>();
		m_softbodies.push_back(softbody);
		m_nsoftbodies++;
	}

	// init constraints

	for (int i = 0; i < m_nconstraints; i++) {
		m_constraints[i]->countDofs(nem, ner, nim, nir);
		if (i < m_nconstraints - 1) {
			m_constraints[i]->next = m_constraints[i + 1];
		}
	}
	
	if (m_nconstraints == 0) {
		addConstraintNull();
	}
}

void World::update() {
	for (int i = 0; i < m_nbodies; i++) {
		m_bodies[i]->update();
	}
}

int World::getNsteps() {
	// Computes the number of results
	int nsteps = (m_tspan(1) - m_tspan(0)) / m_h;
	return nsteps;
}

shared_ptr<Body> World::getBody(int uid) {
	MapBodyUID::const_iterator it = m_bodyUID.find(uid);
	return (it == m_bodyUID.end() ? NULL : it->second);
}

shared_ptr<Body> World::getBody(const string &name) {
	MapBodyName::const_iterator it = m_bodyName.find(name);
	return (it == m_bodyName.end() ? NULL : it->second);
}

shared_ptr<Joint> World::getJoint(int uid) {
	MapJointUID::const_iterator it = m_jointUID.find(uid);
	return (it == m_jointUID.end() ? NULL : it->second);
}

shared_ptr<Joint> World::getJoint(const string &name) {
	MapJointName::const_iterator it = m_jointName.find(name);
	return (it == m_jointName.end() ? NULL : it->second);
}

void World::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, const shared_ptr<Program> progSoft, shared_ptr<MatrixStack> P) {
	// Draw rigid bodies
	for (int i = 0; i < m_nbodies; i++) {
		m_bodies[i]->draw(MV, prog, P);
	}

	// Draw joints
	for (int i = 0; i < m_njoints; i++) {
		m_joints[i]->draw(MV, progSimple, P);
	}

	// Draw constraints
	for (int i = 0; i < 1; i++) {

	}

	// Draw springs
	for (int i = 0; i < m_nsprings; i++) {
		m_springs[i]->draw(MV, prog, progSimple, P);
	}

	// Draw soft bodies
	for (int i = 0; i < m_nsoftbodies; i++) {
		m_softbodies[i]->draw(MV, prog, progSimple, P);
	}
}