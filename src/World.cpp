#include "World.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Joint.h"
#include "JointNull.h"
#include "JointFixed.h"
#include "JointRevolute.h"
#include "JointRevoluteHyperReduced.h"
#include "JointSphericalExp.h"
#include "JointUniversal.h"
#include "JointUniversalYZ.h"
#include "JointUniversalXY.h"

#include "JointFree.h"
#include "JointTranslational.h"
#include "JointSplineCurve.h"
#include "JointSplineSurface.h"

#include "Node.h"
#include "Body.h"
#include "BodyCuboid.h"
#include "SoftBodyNull.h"
#include "SoftBodyInvertibleFEM.h"
#include "SoftBodyCorotationalLinear.h"
#include "MeshEmbedding.h"
#include "MeshEmbeddingNull.h"

#include "SoftBody.h"
#include "FaceTriangle.h"

#include "MatrixStack.h"
#include "Program.h"
#include "SE3.h"
#include "JsonEigen.h"

#include "ConstraintJointLimit.h"
#include "ConstraintNull.h"
#include "ConstraintLoop.h"
#include "ConstraintAttachSpring.h"
#include "ConstraintAttachSoftBody.h"
#include "ConstraintPrescBody.h"
#include "ConstraintPrescJoint.h"

#include "Deformable.h"
#include "DeformableSpring.h"
#include "DeformableNull.h"

#include "Spring.h"
#include "SpringNull.h"
#include "SpringDamper.h"

#include "Comp.h"
#include "CompNull.h"
#include "CompSphere.h"
#include "CompCylinder.h"
#include "CompDoubleCylinder.h"

#include "WrapNull.h"
#include "WrapSphere.h"
#include "WrapCylinder.h"
#include "WrapDoubleCylinder.h"
#include "Vector.h"

#include "TetgenHelper.h"
#include "Line.h"
#include "Surface.h"

#include <json\writer.h>
#include <json\json.h>
#include <json\value.h>

#define D90 M_PI/2.0
#define D180 M_PI
#define D30 M_PI/6.0
#define D45 M_PI/4.0
#define D15 M_PI/12.0
#define D10 M_PI/18.0
#define D60 M_PI/3.0
#define D360 2 * M_PI

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

World::World() :
	nr(0), nm(0), nR(0), nem(0), ner(0), ne(0), nim(0), nir(0), m_nbodies(0), m_njoints(0), m_ndeformables(0), m_constraints(0), m_countS(0), m_countCM(0),
	m_nsoftbodies(0), m_ncomps(0), m_nwraps(0), m_nsprings(0), m_nmeshembeddings(0)
{
	m_energy.K = 0.0;
	m_energy.V = 0.0;

	export_part = 0;
}

World::World(WorldType type) :
	m_type(type),
	nr(0), nm(0), nR(0),nem(0), ner(0), ne(0), nim(0), nir(0), m_nbodies(0), m_njoints(0), m_ndeformables(0), m_nconstraints(0), m_countS(0), m_countCM(0),
	m_nsoftbodies(0), m_ncomps(0), m_nwraps(0), m_nsprings(0), m_nmeshembeddings(0)
{
	m_energy.K = 0.0;
	m_energy.V = 0.0;
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
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 3; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");


			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		//auto con0 = make_shared<ConstraintPrescJoint>(m_joints[0], REDMAX_EULER);
		//m_nconstraints++;
		//m_constraints.push_back(con0);
		break;
	}
	case DIFF_REVOLUTE_AXES:
		break;
	case BRANCHING:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
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

		auto j0 = addJointRevolute(b0, Vector3d::UnitX(), Vector3d(0.0, 15.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
		auto j1 = addJointRevolute(b1, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j0);
		auto j2 = addJointRevolute(b2, Vector3d::UnitX(), Vector3d(-10.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j1);
		auto j3 = addJointRevolute(b3, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j1);
		auto j4 = addJointRevolute(b4, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j2);
		auto j5 = addJointRevolute(b5, Vector3d::UnitX(), Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j4);
		auto j6 = addJointRevolute(b6, Vector3d::UnitY(), Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j4);
		auto j7 = addJointRevolute(b7, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j3);
		auto j8 = addJointRevolute(b8, Vector3d::UnitX(), Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j7);
		auto j9 = addJointRevolute(b9, Vector3d::UnitY(), Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j7);
	}
	break;
	case SHPERICAL_JOINT:
		break;
	case LOOP:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
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

		auto j0 = addJointRevolute(b0, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
		auto j1 = addJointRevolute(b1, Vector3d::UnitZ(), Vector3d(-10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j0);
		auto j2 = addJointRevolute(b2, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j0);
		auto j3 = addJointRevolute(b3, Vector3d::UnitZ(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j1);
		auto j4 = addJointRevolute(b4, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j3);
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
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 10; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
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
	{	
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, 0.0, 0.0;
		Eigen::from_json(js["sides"], sides);
		
		m_stiffness = 1.0e4;
		m_damping = 1.0e3;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 3; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
			m_joints[i]->setStiffness(m_stiffness);
			m_joints[i]->setDamping(m_damping);

		}

		m_joints[0]->m_qdot(0) = 1.0;

		break;

	}

		break;
	case SPRINGS:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		m_stiffness = 5.0e3;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		// Init springs
		auto deformable0 = addDeformableSpring(sides.x()*sides.y()*sides.z() * density, 3, nullptr, Vector3d(10.0 * m_nbodies + 10.0, 10.0, 0.0), m_bodies[m_nbodies - 1], Vector3d(5.0, 0.0, 0.0));
		deformable0->setStiffness(m_stiffness);
		auto deformable1 = addDeformableSpring(sides.x()*sides.y()*sides.z() * density, 2, m_bodies[0], Vector3d(0.0, 0.0, 0.0), m_bodies[m_nbodies - 1], Vector3d(0.0, 0.0, 0.0));
		deformable1->setStiffness(m_stiffness);
		for (int i = 0; i < (int)m_deformables.size(); i++) {
			m_deformables[i]->load(RESOURCE_DIR);
		}

	}
	break;
	case SOFT_BODIES:
	{
		m_h = 1.0e-4;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e2;
		double possion = 0.4;

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);

				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				
			}
			if (i > 0) {
			//addConstraintJointLimit(m_joints[i], -M_PI / 4, M_PI / 4);
			}
		}

		m_joints[0]->m_qdot(0) = 10.0;
		m_joints[1]->m_qdot(0) = -40.0;
		// Init constraints
		
		auto softbody = addSoftBody( 0.001 * density, young, possion, NEO_HOOKEAN, RESOURCE_DIR, "pqz", "box10_1_1");
		softbody->transform(Vector3d(10.0, 0.0, 0.0));
		softbody->setColor(Vector3f(255.0, 204.0, 153.0) / 255.0);

		// auto softbody1 = addSoftBody(0.01 * density, young, possion, RESOURCE_DIR, "cylinder");
		// softbody1->transform(Vector3d(20.0, 0.0, 0.0));

	}
	break;

	case COMPONENT:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);

				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				joint->m_qdot(0) = -5.0;
			}
		}
		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(3.0, 1.0, 0.0));

		auto compSphere = addCompSphere(1.0, m_bodies[0], E, RESOURCE_DIR);
		auto compCylinder = addCompCylinder(1.0, m_bodies[1], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0), RESOURCE_DIR, "obstacle.obj");
		auto compDoubleCylinder = addCompDoubleCylinder(
			0.5, m_bodies[0], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0),
			0.5, m_bodies[2], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0),
			RESOURCE_DIR, "obstacle.obj", "obstacle.obj");
	}
	break;

	case WRAP_SPHERE:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				joint->m_qdot(0) = -5.0;
			}
		}

		auto compSphere = addCompSphere(1.0, m_bodies[1], Matrix4d::Identity(), RESOURCE_DIR);
		auto wrapSphere = addWrapSphere(m_bodies[0], Vector3d(1.0, 0.0, 0.0), m_bodies[2], Vector3d(1.0, 0.0, 0.0), compSphere, 20, RESOURCE_DIR);

	}
	break;

	case WRAP_CYLINDER:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				//joint->m_qdot(0) = -5.0;
			}
		}

		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(0.0, 1.0, 0.0));
		auto compCylinder = addCompCylinder(0.5, m_bodies[1], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0), RESOURCE_DIR, "obstacle.obj");
		auto wrapCylinder = addWrapCylinder(m_bodies[0], Vector3d(1.0, 0.0, 0.0), m_bodies[2], Vector3d(1.0, 0.0, 0.0), compCylinder, 20, RESOURCE_DIR);

	}
	break;

	case WRAP_DOUBLECYLINDER:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				joint->m_qdot(0) = -5.0;
			}
		}

		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(3.0, 1.0, 0.0));
		auto compDoubleCylinder = addCompDoubleCylinder(
			0.5, m_bodies[0], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0),
			0.5, m_bodies[2], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0), 
			RESOURCE_DIR, "obstacle.obj", "obstacle.obj");

		auto wrapDoubleCylinder = addWrapDoubleCylinder(
			m_bodies[0], Vector3d(-5.0, 0.5, 0.0), 
			m_bodies[2], Vector3d(5.0, 0.5, 0.0), 
			Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 0.0), 
			Vector3d(0.0, 0.0, -1.0), Vector3d(0.0, 0.0, 1.0), 
			compDoubleCylinder, 20, RESOURCE_DIR);

	}
	break;

	case SPLINE_CURVE_JOINT:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		auto body0 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_y_9.obj");
		auto joint0 = addJointRevolute(body0, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
		auto body1 = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
		auto joint1 = make_shared<JointSplineCurve>(body1, joint0);
		m_joints.push_back(joint1);
		m_njoints++;
		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(0.0, -10.0, 0.0));
		joint1->setJointTransform(E);
		joint1->load(RESOURCE_DIR, "joint_spline_curve2.obj");
		Matrix4d cf0 = SE3::RpToE(SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), M_PI), Vector3d(-10.0, 0.0, 0.0));
		Matrix4d cf1 = SE3::RpToE(SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), M_PI / 2.0), Vector3d(0.0, 2.0, 0.0));
		Matrix4d cf2 = SE3::RpToE(SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), 0.0), Vector3d(10.0, 0.0, 0.0));
		Matrix4d cf3 = SE3::RpToE(SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), -M_PI / 2.0), Vector3d(0.0, -2.0, 0.0));

		joint1->addControlFrame(cf0);
		joint1->addControlFrame(cf1);
		joint1->addControlFrame(cf2);
		joint1->addControlFrame(cf3);

		auto body2 = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
		auto joint2 = addJointRevolute(body2, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint1);
		E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(10.0, 0.0, 0.0));

		joint2->setJointTransform(E);
		joint1->m_q(0) = - M_PI/4.0;
		joint2->m_q(0) = 15.0 * M_PI / 16.0;

	}
	break;
	case SPLINE_SURFACE_JOINT:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		auto body0 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto joint0 = addJointRevolute(body0, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
		auto body1 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");

		auto joint1 = make_shared<JointSplineSurface>(body1, joint0);
		m_joints.push_back(joint1);
		m_njoints++;
		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(0.0, -14.0, 0.0));
		joint1->setJointTransform(E);

		double t0 = 15.0;
		double r0 = M_PI * 0.25;
		double x, y, z, a, b, c, s1, s2;
		for (int i = 0; i < 4; ++i) {
			s1 = i / 3.0;
			x = (1 - s1)*(-t0) + s1 * t0;
			a = (1 - s1) * (-r0) + s1 * r0;
			for (int j = 0; j < 4; ++j) {
				s2 = j / 3.0;
				y = (1 - s2) * (-t0) + s2 * t0;
				z = 0.05 * (x * x + y * y);
				b = (1 - s1) * (-r0) + s1 * r0;
				c = 0;
				Vector6d ctf;
				ctf << x, z, y, a, c, b;
				joint1->addControlFrame(i,j,ctf);
			}
		}

		auto body2 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto joint2 = addJointRevolute(body2, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint1);
		E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(0.0, -10.0, 0.0));
		joint2->setJointTransform(E);
		joint0->m_q(0) = M_PI / 8.0;
		joint1->m_q(0) = 0.5;
		joint1->m_q(1) = 0.5;
		joint2->m_q(0) = M_PI / 4.0;

	}
	break;
	case SOFT_BODIES_CUBE_INVERTIBLE:
	{
		m_h = 1.0e-3;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e3;
		double possion = 0.40;

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);

				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		m_joints[0]->m_qdot(0) = 10.0;
		m_joints[1]->m_qdot(0) = -40.0;

		auto softbody = addSoftBody(0.001 * density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqz", "box10_1_1");
		softbody->transform(Vector3d(10.0, 0.0, 0.0));
		softbody->setColor(Vector3f(255.0, 204.0, 153.0) / 255.0);

	}
	break;
	case SOFT_BODIES_CYLINDER_INVERTIBLE:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e3;
		double possion = 0.35;

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);

			}
		}

		//m_joints[0]->m_qdot(0) = 10.0;
		//m_joints[1]->m_qdot(0) = -40.0;

		auto softbody = addSoftBodyInvertibleFEM(0.001 * density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqz", "muscle_cyc_cyc"); //
		softbody->transform(Vector3d(10.0, 0.0, 0.0));
		softbody->setColor(Vector3f(255.0, 204.0, 153.0) / 255.0);

	}
	break;
	case SOFT_BODIES_CYLINDER_COROTATIONAL_LINEAR:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e3;
		double possion = 0.35;

		for (int i = 0; i < 2; i++) {
			auto body = addBody( density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);

				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);

			}
		}

		//m_joints[0]->m_qdot(0) = 10.0;
		//m_joints[1]->m_qdot(0) = -40.0;

		auto softbody = addSoftBodyCorotationalLinearFEM(0.001 *density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqz", "muscle_cyc_cyc"); //
		softbody->transform(Vector3d(10.0, 0.0, 0.0));
		softbody->setColor(Vector3f(255.0, 204.0, 153.0) / 255.0);

	}
	break;

	case SOFT_BODIES_CUBE_COROTATIONAL_LINEAR:
	{
		m_h = 1.0e-3;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e3;
		double possion = 0.40;

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);

				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		m_joints[0]->m_qdot(0) = 10.0;
		m_joints[1]->m_qdot(0) = -40.0;

		auto softbody = addSoftBodyCorotationalLinearFEM(0.001 * density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqz", "box10_1_1");
		softbody->transform(Vector3d(10.0, 0.0, 0.0));
		softbody->setColor(Vector3f(255.0, 204.0, 153.0) / 255.0);

	}
	break;
	case SPRING_DAMPER:
	{
		m_h = 1.0e-1;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e2;
		double possion = 0.40;

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			// Inits joints
			if (i == 0) {
				addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		m_joints[1]->m_q(0) = - M_PI / 2.0;
		m_joints[2]->m_q(0) = M_PI / 2.0;

		auto spring = make_shared<SpringDamper>(m_bodies[0], Vector3d(-2.0, -0.5, 0.0), m_bodies[1], Vector3d(2.0, -0.5, 0.0));
		m_springs.push_back(spring);
		m_nsprings++;
		spring->setStiffness(1.0e5);
		spring->setDamping(1.0e3);
		spring->load(RESOURCE_DIR);

	}
	break;

	case MESH_EMBEDDING:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e3;
		double possion = 0.35;
		m_stiffness = 1.0e4;
		m_damping = 1.0e3;

		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(10.0, 0.0, 0.0));

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			if (i == 0) {
				addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);

				//addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);

			}

			// Init constraints
			if (i > 0) {
				//addConstraintJointLimit(m_joints[i], -M_PI / 2.0, M_PI / 2.0);
			}

			//m_joints[i]->setStiffness(m_stiffness);
			m_joints[i]->setDamping(m_damping);

		}
		//m_joints[0]->m_qdot(0) = 10.0;
		m_joints[1]->m_qdot(0) = -10.0;

		auto mesh_embedding = addMeshEmbedding(0.001 *density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqz", "pqz", "embedded", "muscle_cyc_cyc", SOFT_INVERTIBLE); //
		mesh_embedding->transformDenseMesh(E);
		mesh_embedding->transformCoarseMesh(E);
		mesh_embedding->precomputeWeights();

	}
	break;
	case HUMAN_BODY:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e3;
		double possion = 0.35;
		m_stiffness = 1.0e4;
		m_damping = 1.0e3;
		isleftleg = false;
		isrightleg = false;

		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(10.0, 0.0, 0.0));

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}

			// Init constraints
			if (i > 0) {
				//addConstraintJointLimit(m_joints[i], -M_PI / 4.0, M_PI / 4.0);
			}

			m_joints[i]->setStiffness(m_stiffness);
			m_joints[i]->setDamping(m_damping);
		}

		auto torso = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "torso.obj"); // body2
		torso->setColor(Vector3f(143.0f / 255.0f, 188.0f / 255.0f, 143.0f/255.0f));
		addJointFixed(torso, Vector3d(0, 5.0, 0.0), Matrix3d::Identity(), 0.0);
		auto head = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "head.obj");// body3
		head->setColor(Vector3f(245.0f / 255.0f, 222.0f / 255.0f, 179.0f / 255.0f));
		addJointFixed(head, Vector3d(-7.5, 12.0, 0.0), Matrix3d::Identity(), 0.0);

		auto invisible_fixed_left_arm = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");// body4
		invisible_fixed_left_arm->toggleDrawing(false);
		
		auto joint_left_arm_fixed = addJointFixed(invisible_fixed_left_arm, Vector3d(-15, 0.0, 0.0), Matrix3d::Identity(), 0.0);
		
		auto left_arm_0 = addBody(density, sides, Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");// body5
		auto joint_left_arm_0 = addJointRevolute(left_arm_0, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint_left_arm_fixed);
		auto left_arm_1 = addBody(density, sides, Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");// body6
		auto joint_left_arm_1 = addJointRevolute(left_arm_1, Vector3d::UnitZ(), Vector3d(-10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint_left_arm_0);
		joint_left_arm_0->setDamping(m_damping);
		joint_left_arm_1->setDamping(m_damping);

		//addConstraintJointLimit(joint_left_arm_0, -M_PI / 2.0, M_PI / 2.0);
		//addConstraintJointLimit(joint_left_arm_1, -M_PI / 2.0, M_PI / 2.0);


		Vector3f muscle_color = Vector3f(255.0f, 228.0f, 196.0f);
		muscle_color /= 255.0f;	

		if (isleftleg) {
			auto fix_left_leg = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_y_9.obj");
			fix_left_leg->toggleDrawing(false);

			auto joint_left_leg_fixed = addJointFixed(fix_left_leg, Vector3d(-12.0, -13.0, 0.0), Matrix3d::Identity(), 0.0);
			auto left_leg_0 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_y_9.obj");
			auto joint_left_leg0 = addJointRevolute(left_leg_0, Vector3d::UnitX(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint_left_leg_fixed);
			auto left_leg1 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_y_9.obj");
			auto joint_left_leg1 = addJointRevolute(left_leg1, Vector3d::UnitX(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint_left_leg0);


			joint_left_leg1->m_qdot(0) = 20.0;
			addConstraintJointLimit(joint_left_leg0, -M_PI / 4.0, M_PI / 4.0);
			addConstraintJointLimit(joint_left_leg1, -M_PI / 4.0, M_PI / 4.0);
			auto left_leg_mesh = addMeshEmbedding(0.001 *density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqzi", "pqz", "left_leg_coarse", "left_leg_dense", SOFT_INVERTIBLE); //
			left_leg_mesh->precomputeWeights();
			left_leg_mesh->getDenseMesh()->setColor(muscle_color);
		}

		if (isrightleg) {
			auto fix_right_leg = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_y_9.obj");
			fix_right_leg->toggleDrawing(false);
			auto joint_right_leg_fixed = addJointFixed(fix_right_leg, Vector3d(-3.0, -13.0, 0.0), Matrix3d::Identity(), 0.0);
			auto right_leg_0 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_y_9.obj");
			auto joint_right_leg_0 = addJointRevolute(right_leg_0, Vector3d::UnitX(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint_right_leg_fixed);
			auto right_leg_1 = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_y_9.obj");
			auto joint_right_leg_1 = addJointRevolute(right_leg_1, Vector3d::UnitX(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, joint_right_leg_0);

			joint_right_leg_1->m_qdot(0) = -20.0;

			addConstraintJointLimit(joint_right_leg_0, -M_PI / 4.0, M_PI / 4.0);
			addConstraintJointLimit(joint_right_leg_1, -M_PI / 4.0, M_PI / 4.0);
			Matrix4d E_l_r = Matrix4d::Identity();
			E_l_r.block<3, 1>(0, 3) = Vector3d(9.0, 0.0, 0.0);
			auto right_leg_mesh = addMeshEmbedding(0.001 *density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqzi", "pqz", "left_leg_coarse", "left_leg_dense", SOFT_INVERTIBLE); //
			right_leg_mesh->transformCoarseMesh(E_l_r);
			right_leg_mesh->transformDenseMesh(E_l_r);
			right_leg_mesh->precomputeWeights();
			right_leg_mesh->getDenseMesh()->setColor(muscle_color);
		}

		//m_joints[0]->m_qdot(0) = 10.0;
		m_joints[1]->m_qdot(0) = 20*1.2;
		joint_left_arm_1->m_qdot(0) = 20.0;

		//m_joints[1]->m_qdot(0) = 5;

		//joint_left_arm_1->m_qdot(0) = -5;
		//joint_left_leg1->m_qdot(0) = 5.0;
		//joint_right_leg_1->m_qdot(0) = -5.0;

		auto right_arm_mesh = addMeshEmbedding(0.001 *density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqzi", "pqz", "right_arm_coarse_1", "right_arm_dense_1", SOFT_INVERTIBLE); //
		right_arm_mesh->precomputeWeights();
		
		right_arm_mesh->getDenseMesh()->setColor(muscle_color);

		//auto torso = addMeshEmbedding(0.001 *density, young, possion, CO_ROTATED, RESOURCE_DIR, "torso_embedded_w_hole", "torso"); //
		//torso->precomputeWeights();

		auto left_arm_mesh = addMeshEmbedding(0.001 *density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqzi", "pqz", "left_arm_coarse_1", "left_arm_dense_1", SOFT_INVERTIBLE); //
		left_arm_mesh->precomputeWeights();
		left_arm_mesh->getDenseMesh()->setColor(muscle_color);
	
		// Init springs
		auto deformable0 = addDeformableSpring(sides(0)*sides(1)*sides(2)*density, 2, nullptr, Vector3d(25, 0.0, 0.0), m_bodies[1], Vector3d(5.0, 0.0, 0.0));
		deformable0->setStiffness(m_stiffness);
		auto deformable1 = addDeformableSpring(sides(0)*sides(1)*sides(2)*density, 2, nullptr, Vector3d(-45.0, 0.0, 0.0), m_bodies[6], Vector3d(-5.0, 0.0, 0.0));
		deformable1->setStiffness(m_stiffness);

		for (int i = 0; i < (int)m_deformables.size(); i++) {
			m_deformables[i]->load(RESOURCE_DIR);
		}

	}
	break;

	case WORM:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e4;
		double possion = 0.4;
		m_stiffness = 1.0e4;
		m_damping = 1.0e2;

		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(10.0, 0.0, 0.0));

		for (int i = 0; i < 8; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
			//m_joints[i]->setStiffness(m_stiffness);
			m_joints[i]->setDamping(m_damping);
		}

		Vector3d start_pt = Vector3d::Zero();
		Vector3d end_pt = Vector3d(80.0, 0.0, 0.0);
		int nsegments = 8;
		int nsamples = 5;

		vector<std::shared_ptr<Node>> additional_nodes;
		
		for (int i = 0; i < nsegments; ++i) {
			Vector3d n0 = Vector3d(i * 10.0, 0.0, 0.0);
			Vector3d n1 = n0;
			n1.x() += 10.0;
			auto line = make_shared<Line>(n0, n1);
			line->setBody(m_bodies[i]);
			m_lines.push_back(line);

			for (int s = 0; s < nsamples; s++) {
				double f = double(s) / double(nsamples);
				auto node = make_shared<Node>();
				node->i = nsamples * i + s;
				node->x = f * n1 + (1.0 - f) * n0;
				additional_nodes.push_back(node);
			}
		}
	
		TetgenHelper::createNodeFile(additional_nodes, (char *)(RESOURCE_DIR + "worm8_coarse.a.node").c_str());

		Vector3f worm_color = Vector3f(102.0f, 204.0f, 0.0f);
		worm_color /= 255.0f;

		m_joints[0]->m_qdot[0] = -11;

		//m_joints[1]->m_qdot[0] = 5;

		//m_joints[2]->m_qdot[0] = 2;
		//m_joints[3]->m_qdot[0] = 2;

		m_joints[1]->m_qdot[0] = 10;
		m_joints[3]->m_qdot[0] = 2;
		m_joints[4]->m_qdot[0] = 2;
		m_joints[5]->m_qdot[0] = 2;

		m_joints[6]->m_qdot[0] = 3;
		m_joints[7]->m_qdot[0] = 3;
		
		Floor f0(-31.0f, Vector2f(-80.0f, 80.0f), Vector2f(-80.0f, 80.0f));
		m_floors.push_back(f0);

		auto worm = addMeshEmbedding(density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqzi", "pqz", "worm8_coarse", "worm8_dense", SOFT_INVERTIBLE);
		//worm->transformCoarseMesh(SE3::RpToE(Matrix3d::Identity(), Vector3d(40.0, 0.0, 0.0)));
		worm->transformDenseMesh(SE3::RpToE(Matrix3d::Identity(), Vector3d(40.0, 0.0, 0.0)));
		worm->precomputeWeights();
		worm->setDamping(10.0);
		worm->getDenseMesh()->setColor(worm_color);
		worm->getCoarseMesh()->setFloor(-31.0);
		worm->getCoarseMesh()->setYthreshold(0.0);
	}
	break;

	case CROSS: {
		m_h = 1.0e-1;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98.0*5.0, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e3;
		double possion = 0.4;
		m_stiffness = 1.0e4;
		m_damping = 1.0e1;

		nlegs = 4;
		nsegments = 4;
		double len_segment = 10.0;
		double rotation = 360.0 / nlegs;

		for (int k = 0; k < nlegs; ++k) {
			for (int i = 0; i < nsegments; i++) {
				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
				if (i == 0) {
					//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
					addJointRevolute(body, Vector3d::UnitZ(), Vector3d::Zero(), SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), (0.0 + rotation * k) / 180.0 * M_PI), 0.0, RESOURCE_DIR);
				}
				else {
					auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(len_segment, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[nsegments * k + i - 1]);
				}
				//m_joints[i]->setStiffness(m_stiffness);
				m_joints[i]->setDamping(m_damping);
			}
		}
		//auto con0 = make_shared<ConstraintPrescBody>(m_bodies[3], Vector3d(1, 3, 5), REDMAX_EULER);
		//auto con0 = make_shared<ConstraintPrescJoint>(m_joints[3], REDMAX_EULER);
		//m_nconstraints++;
		//m_constraints.push_back(con0);
		//scene.constraints{ end + 1 } = reduced.ConstraintPrescBody(scene.bodies{ end }, [2 4 6], vel);
		//scene.constraints{ end + 1 } = reduced.ConstraintPrescJoint(scene.joints{ 3 }, vel);

		int nsamples = 4;		
		vector<std::shared_ptr<Node>> additional_nodes;
		int idx_body = 0;
		addSkeleton(Vector3d::Zero(), Vector3d(40.0, 0.0, 0.0), nsegments, nsamples, idx_body, additional_nodes);
		addSkeleton(Vector3d::Zero(), Vector3d(0.0, 0.0, -40.0), nsegments, nsamples, idx_body, additional_nodes);
		addSkeleton(Vector3d::Zero(), Vector3d(-40.0, 0.0, 0.0), nsegments, nsamples, idx_body, additional_nodes);
		addSkeleton(Vector3d::Zero(), Vector3d(0.0, 0.0, 40.0), nsegments, nsamples, idx_body, additional_nodes);

		TetgenHelper::createNodeFile(additional_nodes, (char *)(RESOURCE_DIR + "cross_coarse.a.node").c_str());

		Vector3f cross_color = Vector3f(102.0f, 204.0f, 0.0f);
		cross_color /= 255.0f;

		m_joints[0]->m_qdot[0] = 5.0;
		m_joints[1]->m_qdot[0] = 5;
		m_joints[2]->m_qdot[0] = 5;
		m_joints[3]->m_qdot[0] = 5;

		m_joints[4]->m_qdot[0] = 5.0;
		m_joints[5]->m_qdot[0] = 6.0;
		m_joints[6]->m_qdot[0] = 3.0;
		m_joints[7]->m_qdot[0] = 3.0;
		m_joints[8]->m_qdot[0] = 2.0;
		//m_joints[12]->m_qdot[0] = 5;
		double y_floor = -21.0;

		Floor f0(float(y_floor), Vector2f(-80.0f, 80.0f), Vector2f(-80.0f, 80.0f));
		m_floors.push_back(f0);

		auto cross = addMeshEmbedding(density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqzi", "pqz", "cross_coarse", "cross_dense", SOFT_INVERTIBLE);
		cross->precomputeWeights();
		cross->setDamping(15.0);
		cross->getDenseMesh()->setColor(cross_color);
		cross->getCoarseMesh()->setFloor(y_floor);

	}
		break;

	case STARFISH:
	{
		m_h = 1.0e-1;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -980.0, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e4;
		double possion = 0.35;
		m_stiffness = 1.0e4;
		m_damping = 1.0e3;
		double mesh_damping = 210.0;
		double y_floor = -40.0;
		nlegs = 5;
		nsegments = 4;
		double rotation = 360.0 / nlegs;
		Vector3f starfish_color = Vector3f(255.0f, 99.0f, 71.0f);
		starfish_color /= 255.0f;
		std::string coarse_file_name = "thin_coarse";//starfishco  starfish_coarse3
		std::string dense_file_name = "starfish2";
		//std::string coarse_file_name = "starfish_coarse_new0";//starfishcostarfish_extrude_0
		//std::string dense_file_name = "SFdense";

		double len_segment = 20.0;
		Vector3d root_sides;
		Eigen::from_json(js["cube_sides"], root_sides);
		auto root_body = addBody(density, root_sides, Vector3d::Zero(), Matrix3d::Identity(), RESOURCE_DIR, "box_1_1_1.obj");
		/*	VectorXi dof(1);
		dof << 4;
		addConstraintPrescBody(root_body, dof);*/
		//body->setDrawingOption(false);
		Vector6d qdot0;
		qdot0 << 0.0, 0.0, 0.0, 0.0, 6.0, 0.0;
		auto root_joint = addJointFree(root_body, Vector3d::Zero(), Matrix3d::Identity(), Vector6d::Zero(), qdot0, RESOURCE_DIR);

		for (int k = 0; k < nlegs; ++k) {
			for (int i = 0; i < nsegments; i++) {
				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
				body->setDrawingOption(false);

				shared_ptr<Joint> joint;
				if (i == 0) {
					//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
					Vector6d q0;
					q0 << 0.0, 0.0, 0.0, 0.0, 5.0, 0.0;
					//joint = addJointFree(body, Vector3d::Zero(), SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), (18.0 + rotation * k) / 180.0 * M_PI),  Vector6d::Zero(), q0, RESOURCE_DIR);
					joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d::Zero(), SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), (18.0 + rotation * k)/180.0 * M_PI), 0.0, RESOURCE_DIR, root_joint);
					//addConstraintPrescJoint(joint);
				
				}
				else {
					joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(len_segment, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[nsegments * k + i - 1 + 1]);
					//addConstraintPrescJoint(joint);

				}
				addConstraintPrescJoint(joint);
				//m_joints[i]->setStiffness(m_stiffness);
				m_joints[i]->setDamping(2.0);
			}
		}
		VectorXi dof(6);
		dof << 0, 1, 2, 3, 4, 5;
		VectorXi dof3(3);
		dof3 << 3, 4, 5;

		VectorXi dof_(1);
		dof_ << 4;

		VectorXi dofx(1);
		dofx << 3;

		VectorXi dofxy(2);
		dofxy << 3, 4;
		addConstraintPrescBody(m_bodies[8], dof3);
		addConstraintPrescBody(m_bodies[16], dof_);
		addConstraintPrescBody(m_bodies[0], dof);
		addConstraintPrescBody(m_bodies[20], dof3);
		addConstraintPrescBody(m_bodies[4], dof3);
		addConstraintPrescBody(m_bodies[12], dofxy);
		addConstraintPrescBody(m_bodies[1], dofx);


		double len_skeleton = nsegments * len_segment;
		int nsamples = 4;

		vector<std::shared_ptr<Node>> additional_nodes;
		int idx_body = 1;
		Vector3d end_pt;
		for (int i = 0; i < nlegs; ++i) {
			double theta = (18.0 + rotation * i) / 180.0 * M_PI;
			end_pt = len_skeleton * Vector3d(cos(theta), 0.0, -sin(theta));
			addSkeleton(Vector3d::Zero(), end_pt, nsegments, nsamples, idx_body, additional_nodes);
		}

		TetgenHelper::createNodeFile(additional_nodes, (char *)(RESOURCE_DIR + coarse_file_name + ".a.node").c_str());

		//m_joints[0]->m_qdot[0] = -5.0;
		//m_joints[8]->m_qdot[0] = -5.0;

		//m_joints[16]->m_qdot[0] = -5.0;
		////m_joints[16]->m_qdot[0] = -0.7;
		//m_joints[4]->m_qdot[0] = -7.0;
			
		Floor f0(float(y_floor), Vector2f(-150.0f, 150.0f), Vector2f(-150.0f, 150.0f));
		m_floors.push_back(f0);

		auto starfish = addMeshEmbedding(density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqziV", "pq1.4a60.0z", coarse_file_name, dense_file_name, SOFT_INVERTIBLE);//pqz
		starfish->precomputeWeights();
		starfish->setDamping(mesh_damping);
		//starfish->toggleDrawingDenseMesh(false);
		starfish->getDenseMesh()->setColor(starfish_color);
		//starfish->getDenseMesh()->setFloor(y_floor);
		starfish->getCoarseMesh()->setFloor(y_floor );

	}
	break;
	case FREEJOINT: 
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -1.0, 0.0;
		Eigen::from_json(js["cube_sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 1; i++) {

			auto body = addBody(density, sides, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box_1_1_1.obj");
			// Inits joints
			if (i == 0) {
				Vector6d qdot0;
				qdot0 << 0.2, 0.4, 0.6, 0.0, 0.0, 3.0;
				addJointFree(body, Vector3d::Zero(), Matrix3d::Identity(), Vector6d::Zero(),qdot0, RESOURCE_DIR);
				
			}
			else {
				
			}
		}
		break;
	}
	case STARFISH_2:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, 0.0, 0.0;
		Eigen::from_json(js["sides"], sides);
		Vector3d root_sides;
		Eigen::from_json(js["cube_sides"], root_sides);
		double young = 1e4;
		double possion = 0.35;
		m_stiffness = 1.0e1;
		m_damping = 1.0e3;
		double mesh_damping = 1.0;
		double y_floor = -40.0;
		nlegs = 5;
		nsegments = 8;
		double rotation = 360.0 / nlegs;
		Vector3f starfish_color = Vector3f(255.0f, 99.0f, 71.0f);
		starfish_color /= 255.0f;
		std::string coarse_file_name = "starfish_coarse3";//starfishco
		std::string dense_file_name = "starfish2";

		double len_segment = 10.0;
		auto root_body = addBody(density, root_sides, Vector3d::Zero(), Matrix3d::Identity(), RESOURCE_DIR, "box_1_1_1.obj");
	/*	VectorXi dof(1);
		dof << 4;
		addConstraintPrescBody(root_body, dof);*/
		//body->setDrawingOption(false);
		Vector6d qdot0;
		qdot0 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
		auto root_joint = addJointFree(root_body, Vector3d::Zero(), Matrix3d::Identity(), Vector6d::Zero(), qdot0, RESOURCE_DIR);
		//auto root_joint = addJointFixed(root_body, Vector3d::Zero(), Matrix3d::Identity(), 0.0);
		for (int k = 0; k < nlegs; ++k) {
			for (int i = 0; i < nsegments; i++) {
				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
				//body->setDrawingOption(false);

				shared_ptr<Joint> joint;
				if (i == 0) {
					//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
					joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d::Zero(), SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), (18.0 + rotation * k)/180.0 * M_PI), 0.0, RESOURCE_DIR, root_joint);
				}
				else {
					joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(len_segment, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[nsegments * k + i - 1 + 1]);
					//addConstraintPrescJoint(joint);

				}
				
				if (k * nsegments + i == 7) {
					VectorXi dof(1);
					dof << 4;
					addConstraintPrescBody(body, dof);
					body->setDrawingOption(false);
				}
				//addConstraintPrescJoint(joint);
				//m_joints[i]->setStiffness(m_stiffness);
				//m_joints[i]->setDamping(m_damping);
			}
		}

		double len_skeleton = nsegments * len_segment;
		int nsamples = 4;

		vector<std::shared_ptr<Node>> additional_nodes;
		int idx_body = 1;
		Vector3d end_pt;
		for (int i = 0; i < nlegs; ++i) {
			double theta = (18.0 + rotation * i) / 180.0 * M_PI;
			end_pt = len_skeleton * Vector3d(cos(theta), 0.0, -sin(theta));
			addSkeleton(Vector3d::Zero(), end_pt, nsegments, nsamples, idx_body, additional_nodes);
		}

		TetgenHelper::createNodeFile(additional_nodes, (char *)(RESOURCE_DIR + coarse_file_name + ".a.node").c_str());

		//m_joints[0]->m_qdot[0] = -5.0;
		//m_joints[8]->m_qdot[0] = -5.0;

		//m_joints[16]->m_qdot[0] = -5.0;
		////m_joints[16]->m_qdot[0] = -0.7;
		//m_joints[4]->m_qdot[0] = -7.0;

		Floor f0(float(y_floor), Vector2f(-80.0f, 80.0f), Vector2f(-80.0f, 80.0f));
		m_floors.push_back(f0);

		auto starfish = addMeshEmbedding(density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqziYV", "pqz", dense_file_name, dense_file_name, SOFT_INVERTIBLE);
		starfish->precomputeWeights();
		starfish->setDamping(mesh_damping);
		starfish->getDenseMesh()->setColor(starfish_color);
		//starfish->getCoarseMesh()->setFloor(y_floor);

	}
	break;
	case TEST_REDUCED_HYBRID_DYNAMICS:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		auto con0 = make_shared<ConstraintPrescJoint>(m_joints[0], REDMAX_EULER);
		m_nconstraints++;
		m_constraints.push_back(con0);


		break;
	}
	case TEST_MAXIMAL_HYBRID_DYNAMICS:
	{
		m_h = 5.0e-2;
		density = 1.0;
		m_grav << 0.0, -980, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 4; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), -M_PI / 2.0), 0.0, RESOURCE_DIR);
			}
			else if (i == 1) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), M_PI / 2.0), 0.0, RESOURCE_DIR, m_joints[i - 1]);

			}
			else if (i == 2) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), M_PI / 2.0), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), -M_PI / 2.0), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}
		Vector3i dof;
		dof << 2, 3, 4;
		auto con0 = make_shared<ConstraintPrescBody>(m_bodies[3], dof, REDMAX_EULER);
		m_nconstraints++;
		m_constraints.push_back(con0);
		break;
	}
	case FINGERS:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -0.0, 0.0;
		Eigen::from_json(js["sides"], sides);

		double y_floor = 0.0;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		
		double scale = 4.0;
		double len_thumb0 = 11.5 * scale;
		double len_thumb1 = 6.0 * scale;
		double len_thumb2 = 4.0 * scale;

		double len_index_0 = 12.0 * scale;
		double len_index_1 = 7.0 * scale;
		double len_index_2 = 4.0 * scale;
		double len_index_3 = 3.5 *scale;

		double len_middle_0 = 11.0 * scale;
		double len_middle_1 = 7.0 * scale;
		double len_middle_2 = 4.5 * scale;
		double len_middle_3 = 3.5 * scale;

		double len_ring_0 = 10.0 * scale;
		double len_ring_1 = 6.5 * scale;
		double len_ring_2 = 4.2 * scale;
		double len_ring_3 = 3.3 * scale;

		double len_pinky_0 = 9.0 * scale;
		double len_pinky_1 =  4.0* scale;
		double len_pinky_2 = 3.0 * scale;
		double len_pinky_3 =  2.8* scale;

		double len_elbow = 36.0* scale;

		double short_side = 4.0;

		Floor f0(float(y_floor), Vector2f(-50.0f, 500.0f), Vector2f(-150.0f, 500.0f));
		m_floors.push_back(f0);
		Matrix3d changeXY = SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), D90) * SE3::aaToMat(Vector3d(1, 0, 0), D90);
		Matrix3d postXY = SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), -D90) * SE3::aaToMat(Vector3d(1, 0, 0), -D90);
		Matrix3d preXY = SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), D90);

		auto elbow = addBody(density, Vector3d(len_elbow, short_side, short_side), Vector3d(len_elbow / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "36.obj");
		auto j_elbow = addJointRevolute(elbow, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(),0.0), 0.0, RESOURCE_DIR);

		auto wrist = addBody(density, Vector3d(short_side, short_side, short_side), Vector3d(0.0, 0.0, 0.0), preXY, RESOURCE_DIR, "wrist.obj");
		auto j_wrist = addJointUniversalXY(wrist, Vector3d(len_elbow, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), D30) * changeXY, RESOURCE_DIR, j_elbow);
//
		auto thumb_0 = addBody(density, Vector3d(len_thumb0, short_side, short_side), Vector3d(0.0, 0.0, len_thumb0 / 2.0), preXY, RESOURCE_DIR, "11_5.obj");
		auto thumb_1 = addBody(density, Vector3d(len_thumb1, short_side, short_side), Vector3d(len_thumb1 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "6.obj");
		auto thumb_2 = addBody(density, Vector3d(len_thumb2, short_side, short_side), Vector3d(len_thumb2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4.obj");

		auto index_finger_0 = addBody(density, Vector3d(len_index_0, short_side, short_side), Vector3d(len_index_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "12.obj");
		auto index_finger_1 = addBody(density, Vector3d(len_index_1, short_side, short_side), Vector3d(0.0, 0.0, len_index_1 / 2.0), preXY, RESOURCE_DIR, "7.obj");
		auto index_finger_2 = addBody(density, Vector3d(len_index_2, short_side, short_side), Vector3d(len_index_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4.obj");
		auto index_finger_3 = addBody(density, Vector3d(len_index_3, short_side, short_side), Vector3d(len_index_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3_5.obj");
		
		auto middle_finger_0 = addBody(density, Vector3d(len_middle_0, short_side, short_side), Vector3d(len_middle_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "11.obj");
		auto middle_finger_1 = addBody(density, Vector3d(len_middle_1, short_side, short_side), Vector3d(0.0, 0.0, len_middle_1 / 2.0), preXY, RESOURCE_DIR, "7.obj");
		auto middle_finger_2 = addBody(density, Vector3d(len_middle_2, short_side, short_side), Vector3d(len_middle_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4_5.obj");
		auto middle_finger_3 = addBody(density, Vector3d(len_middle_3, short_side, short_side), Vector3d(len_middle_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3_5.obj");
		
		auto ring_finger_0 = addBody(density, Vector3d(len_ring_0, short_side, short_side), Vector3d(len_ring_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "10.obj");
		auto ring_finger_1 = addBody(density, Vector3d(len_ring_1, short_side, short_side), Vector3d(0.0, 0.0, len_ring_1 / 2.0), preXY, RESOURCE_DIR, "6_5.obj");
		auto ring_finger_2 = addBody(density, Vector3d(len_ring_2, short_side, short_side), Vector3d(len_ring_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4_2.obj");
		auto ring_finger_3 = addBody(density, Vector3d(len_ring_3, short_side, short_side), Vector3d(len_ring_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3_3.obj");
		
		auto pinky_finger_0 = addBody(density, Vector3d(len_pinky_0, short_side, short_side), Vector3d(len_pinky_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "9.obj");
		auto pinky_finger_1 = addBody(density, Vector3d(len_pinky_1, short_side, short_side), Vector3d(0.0, 0.0, len_pinky_1 / 2.0), preXY, RESOURCE_DIR, "4.obj");
		auto pinky_finger_2 = addBody(density, Vector3d(len_pinky_2, short_side, short_side), Vector3d(len_pinky_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3.obj");
		auto pinky_finger_3 = addBody(density, Vector3d(len_pinky_3, short_side, short_side), Vector3d(len_pinky_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "2_8.obj");

		auto j_thumb_0 = addJointUniversalXY(thumb_0, Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), 0.0) * SE3::aaToMat(Vector3d::UnitX(), -M_PI/10.0)  * SE3::aaToMat(Vector3d::UnitY(), -D15), RESOURCE_DIR, j_wrist);
		auto j_thumb_1 = addJointRevolute(thumb_1, Vector3d::UnitZ(), Vector3d(0.0, 0.0, len_thumb0), postXY*SE3::aaToMat(Vector3d::UnitZ(), -D30/1.5), 0.0, RESOURCE_DIR, j_thumb_0);
		auto j_thumb_2 = addJointRevolute(thumb_2, Vector3d::UnitZ(), Vector3d(len_thumb1, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D30), 0.0, RESOURCE_DIR, j_thumb_1);

		auto j_index_finger_0 = addJointFixed(index_finger_0, Vector3d(0.0, 0.0, 0.0), postXY, 0.0, j_wrist);
		auto j_index_finger_1 = addJointUniversalXY(index_finger_1, Vector3d(len_index_0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D30)*changeXY, RESOURCE_DIR, j_index_finger_0);
		auto j_index_finger_2 = addJointRevolute(index_finger_2, Vector3d::UnitZ(), Vector3d(0.0, 0.0, len_index_1), postXY*SE3::aaToMat(Vector3d::UnitZ(), -D45), 0.0, RESOURCE_DIR, j_index_finger_1);
		auto j_index_finger_3 = addJointRevoluteHyperReduced(index_finger_3, Vector3d::UnitZ(), j_index_finger_2, 0.5, Vector3d(len_index_2, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D45), 0.0, RESOURCE_DIR, j_index_finger_2);

		auto j_middle_finger_0 = addJointFixed(middle_finger_0, Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitX(), M_PI / 10.0)* postXY, 0.0, j_wrist);
		auto j_middle_finger_1 = addJointUniversalXY(middle_finger_1, Vector3d(len_middle_0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D30 / 1.3)*changeXY, RESOURCE_DIR, j_middle_finger_0);
		auto j_middle_finger_2 = addJointRevolute(middle_finger_2, Vector3d::UnitZ(), Vector3d(0.0, 0.0, len_middle_1), postXY*SE3::aaToMat(Vector3d::UnitZ(), -D45), 0.0, RESOURCE_DIR, j_middle_finger_1);
		auto j_middle_finger_3 = addJointRevoluteHyperReduced(middle_finger_3, Vector3d::UnitZ(), j_middle_finger_2, 0.5, Vector3d(len_middle_2, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D45), 0.0, RESOURCE_DIR, j_middle_finger_2);
		
		auto j_ring_finger_0 = addJointFixed(ring_finger_0, Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), 0.0) * SE3::aaToMat(Vector3d::UnitX(), 2 * M_PI / 10.0) * postXY, 0.0,j_wrist);
		auto j_ring_finger_1 = addJointUniversalXY(ring_finger_1, Vector3d(len_ring_0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D30/1.5)*changeXY, RESOURCE_DIR, j_ring_finger_0);
		auto j_ring_finger_2 = addJointRevolute(ring_finger_2, Vector3d::UnitZ(), Vector3d(0.0, 0.0, len_ring_1), postXY*SE3::aaToMat(Vector3d::UnitZ(), -D30), 0.0 , RESOURCE_DIR, j_ring_finger_1);
		auto j_ring_finger_3 = addJointRevoluteHyperReduced(ring_finger_3, Vector3d::UnitZ(), j_ring_finger_2, 0.5, Vector3d(len_ring_2, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D45), 0.0, RESOURCE_DIR, j_ring_finger_2);
		
		auto j_pinky_finger_0 = addJointFixed(pinky_finger_0, Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), 0.0) * SE3::aaToMat(Vector3d::UnitX(), 3 * M_PI / 10.0)* postXY, 0.0, j_wrist);
		auto j_pinky_finger_1 = addJointUniversalXY(pinky_finger_1, Vector3d(len_pinky_0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D10)*changeXY, RESOURCE_DIR, j_pinky_finger_0);
		auto j_pinky_finger_2 = addJointRevolute(pinky_finger_2, Vector3d::UnitZ(), Vector3d(0.0, 0.0, len_pinky_1), postXY*SE3::aaToMat(Vector3d::UnitZ(), -D30), 0.0, RESOURCE_DIR, j_pinky_finger_1);
		auto j_pinky_finger_3 = addJointRevoluteHyperReduced(pinky_finger_3, Vector3d::UnitZ(), j_pinky_finger_2, 0.5, Vector3d(len_pinky_2, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), -D15), 0.0, RESOURCE_DIR, j_pinky_finger_2);


		for (int i = 0; i < (int)m_joints.size(); ++i) {
			m_joints[i]->setDrawRadius(2.0);
			m_joints[i]->setStiffness(7e8);
			m_joints[i]->setDamping(1e1);
		}
		m_joints[eBone_Thumb0]->setStiffness(1e8);
		m_joints[eBone_Thumb1]->setStiffness(1e8);
		m_joints[eBone_Thumb2]->setStiffness(1e8);


		VectorXi dof(3);
		dof << 3, 4, 5; // y direction 
		// stay on table
		addConstraintPrescBody(thumb_2, dof);
		addConstraintPrescBody(index_finger_3, dof);
		addConstraintPrescBody(middle_finger_3, dof);
		addConstraintPrescBody(ring_finger_3, dof);
		addConstraintPrescBody(pinky_finger_3, dof);

		addConstraintPrescJoint(j_elbow); // fixed
		addConstraintPrescJoint(j_wrist); // control z angle

		break;
	}
	case STARFISH3:
	{
		m_h = 1.0e-1;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -980.0, 0.0;
		Eigen::from_json(js["sides"], sides);
		double young = 1e4;
		double possion = 0.35;
		m_stiffness = 1.0e4;
		m_damping = 1.0e3;
		double mesh_damping = 10.0;
		double y_floor = -40.0;
		nlegs = 5;
		nsegments = 4;
		double rotation = 360.0 / nlegs;
		Vector3f starfish_color = Vector3f(255.0f, 99.0f, 71.0f);
		starfish_color /= 255.0f;

		std::string dense_file_name = "SFdense";

		double len_segment = 20.0;
		Vector3d root_sides;
		Eigen::from_json(js["cube_sides"], root_sides);
		auto root_body = addBody(density, root_sides, Vector3d::Zero(), Matrix3d::Identity(), RESOURCE_DIR, "box_1_1_1.obj");
		/*	VectorXi dof(1);
		dof << 4;
		addConstraintPrescBody(root_body, dof);*/
		//body->setDrawingOption(false);
		Vector6d qdot0;
		qdot0 << 0.0, 0.0, 0.0, 0.0, 6.0, 0.0;
		auto root_joint = addJointFree(root_body, Vector3d::Zero(), Matrix3d::Identity(), Vector6d::Zero(), qdot0, RESOURCE_DIR);

		for (int k = 0; k < nlegs; ++k) {
			for (int i = 0; i < nsegments; i++) {
				auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");
				body->setDrawingOption(false);

				shared_ptr<Joint> joint;
				if (i == 0) {
					//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
					Vector6d q0;
					q0 << 0.0, 0.0, 0.0, 0.0, 5.0, 0.0;
					//joint = addJointFree(body, Vector3d::Zero(), SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), (18.0 + rotation * k) / 180.0 * M_PI),  Vector6d::Zero(), q0, RESOURCE_DIR);
					joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d::Zero(), SE3::aaToMat(Vector3d(0.0, 1.0, 0.0), (15.0 + rotation * k) / 180.0 * M_PI), 0.0, RESOURCE_DIR, root_joint);
					//addConstraintPrescJoint(joint);

				}
				else {
					joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(len_segment, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[nsegments * k + i - 1 + 1]);
					//addConstraintPrescJoint(joint);

				}
				addConstraintPrescJoint(joint);
				m_joints[i]->setStiffness(m_stiffness);
				m_joints[i]->setDamping(4.0);
			}
		}
		VectorXi dof(6);
		dof << 0, 1, 2, 3, 4, 5;
		VectorXi dof_(1);
		dof_ << 4;
		addConstraintPrescBody(m_bodies[8], dof_);
		addConstraintPrescBody(m_bodies[16], dof_);
		addConstraintPrescBody(m_bodies[0], dof);
		addConstraintPrescBody(m_bodies[20], dof_);
		addConstraintPrescBody(m_bodies[4], dof_);
		addConstraintPrescBody(m_bodies[12], dof_);

		double len_skeleton = nsegments * len_segment;
		int nsamples = 4;

		vector<std::shared_ptr<Node>> additional_nodes;
		int idx_body = 1;
		Vector3d end_pt;
		for (int i = 0; i < nlegs; ++i) {
			double theta = (15.0 + rotation * i) / 180.0 * M_PI;
			end_pt = len_skeleton * Vector3d(cos(theta), 0.0, -sin(theta));
			addSkeleton(Vector3d::Zero(), end_pt, nsegments, nsamples, idx_body, additional_nodes);
		}

		TetgenHelper::createNodeFile(additional_nodes, (char *)(RESOURCE_DIR + dense_file_name + ".a.node").c_str());

		Floor f0(float(y_floor), Vector2f(-80.0f, 80.0f), Vector2f(-80.0f, 80.0f));
		m_floors.push_back(f0);
		auto starfish = addSoftBodyInvertibleFEM(density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqziV", dense_file_name);
		//auto starfish = addMeshEmbedding(density, young, possion, CO_ROTATED, RESOURCE_DIR, "pqziVY", "pqz", dense_file_name, dense_file_name, SOFT_INVERTIBLE);
		//starfish->precomputeWeights();
		starfish->setDamping(mesh_damping);
		//starfish->toggleDrawingDenseMesh(false);
		starfish->setColor(starfish_color);
		//starfish->getDenseMesh()->setFloor(y_floor);
		starfish->setFloor(y_floor);

	}

	break;
	case TEST_HYPER_REDUCED_COORDS:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -9.8, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 2; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			//// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				//addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				addJointRevoluteHyperReduced(body, Vector3d::UnitZ(), m_joints[i-1], 0.8, Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}
		VectorXi dof(1);
		dof << 2;
		addConstraintPrescJoint(m_joints[0]);		
		addConstraintPrescBody(m_bodies[1], dof);
	}
	break;
	case TEST_JOINT_UNIVERSAL: {

		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -9.8, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		for (int i = 0; i < 3; i++) {

			auto body = addBody(density, sides, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");

			//// Inits joints
			if (i == 0) {
				addJointUniversalYZ(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR);
			}
			else {
				addJointUniversalYZ(body, Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, m_joints[i-1]);

			}
		}
		m_joints[0]->m_q(0) = M_PI / 8.0;
		m_joints[1]->m_q(1) = M_PI / 8.0;
		m_joints[2]->m_q(0) = M_PI / 8.0;

	}
		break;
default:
		break;
	}
}

shared_ptr<ConstraintPrescJoint> World::addConstraintPrescJoint(shared_ptr<Joint> j) {
	auto con = make_shared<ConstraintPrescJoint>(j, REDMAX_EULER);
	m_nconstraints++;
	m_constraints.push_back(con);
	return con;
}

shared_ptr<ConstraintPrescBody> World::addConstraintPrescBody(shared_ptr<Body> b, Eigen::VectorXi dof) {
	auto con = make_shared<ConstraintPrescBody>(b, dof, REDMAX_EULER);
	m_nconstraints++;
	m_constraints.push_back(con);
	return con;
}

shared_ptr<SoftBody> World::addSoftBody(double density, double young, double possion, Material material, const string &RESOURCE_DIR, const string &TETGEN_FLAGS, string file_name) {
	auto softbody = make_shared<SoftBody>(density, young, possion, material);
	softbody->load(RESOURCE_DIR, file_name, TETGEN_FLAGS);
	m_softbodies.push_back(softbody);
	m_nsoftbodies++;
	return softbody;
}

shared_ptr<SoftBodyInvertibleFEM> World::addSoftBodyInvertibleFEM(double density, double young, double possion, Material material, const string &RESOURCE_DIR, const string &TETGEN_FLAGS, string file_name) {
	auto softbody = make_shared<SoftBodyInvertibleFEM>(density, young, possion, material);
	softbody->load(RESOURCE_DIR, file_name, TETGEN_FLAGS);
	m_softbodies.push_back(softbody);
	m_nsoftbodies++;
	return softbody;
}

shared_ptr<SoftBodyCorotationalLinear> World::addSoftBodyCorotationalLinearFEM(double density, double young, double possion, Material material, const string &RESOURCE_DIR, const string &TETGEN_FLAGS, string file_name) {
	auto softbody = make_shared<SoftBodyCorotationalLinear>(density, young, possion, material);
	softbody->load(RESOURCE_DIR, file_name, TETGEN_FLAGS);
	m_softbodies.push_back(softbody);
	m_nsoftbodies++;
	return softbody;
}

shared_ptr<Body> World::addBody(double density, Vector3d sides, Vector3d p, Matrix3d R, const string &RESOURCE_DIR, string file_name) {
	auto body = make_shared<BodyCuboid>(density, sides);
	Matrix4d E = SE3::RpToE(R, p);
	body->setTransform(E);
	body->load(RESOURCE_DIR, file_name);
	m_bodies.push_back(body);
	m_nbodies++;
	return body;
}

shared_ptr<JointRevolute> World::addJointRevolute(shared_ptr<Body> body, 
	Vector3d axis, 
	Vector3d p, 
	Matrix3d R, 
	double q, 
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) {
	auto joint = make_shared<JointRevolute>(body, axis, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->m_q(0) = q;
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointUniversalYZ> World::addJointUniversalYZ(shared_ptr<Body> body,
	Vector3d p,
	Matrix3d R,
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) {
	auto joint = make_shared<JointUniversalYZ>(body, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointUniversalXY> World::addJointUniversalXY(shared_ptr<Body> body,
	Vector3d p,
	Matrix3d R,
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) {
	auto joint = make_shared<JointUniversalXY>(body, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointRevoluteHyperReduced> World::addJointRevoluteHyperReduced(shared_ptr<Body> body,
	Vector3d axis,
	shared_ptr<Joint> friend_joint,
	double scalar,
	Vector3d p,
	Matrix3d R,
	double q,
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) {
	auto joint = make_shared<JointRevoluteHyperReduced>(body, axis, friend_joint, scalar, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->m_q(0) = q;
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointFree> World::addJointFree(
	shared_ptr<Body> body,
	Vector3d p,
	Matrix3d R,
	Vector6d q0,
	Vector6d qdot0,
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) 
{
	auto joint = make_shared<JointFree>(body, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->m_q = q0;
	joint->m_qdot = qdot0;
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointFixed> World::addJointFixed(shared_ptr<Body> body, Vector3d p, Matrix3d R, double q, shared_ptr<Joint> parent) {
	auto joint = make_shared<JointFixed>(body, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	
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

shared_ptr<DeformableSpring> World::addDeformableSpring(double mass, int n_points, shared_ptr<Body> body0, Vector3d r0, shared_ptr<Body> body1, Vector3d r1) {

	auto deformable = make_shared<DeformableSpring>(n_points, m_countS, m_countCM);
	m_deformables.push_back(deformable);
	deformable->setStiffness(m_stiffness);
	deformable->setMass(mass);
	deformable->setAttachments(body0, r0, body1, r1);
	m_ndeformables++;
	return deformable;
}

shared_ptr<CompSphere> World::addCompSphere(double r, shared_ptr<Body> parent, Matrix4d E, const string &RESOURCE_DIR) {
	auto comp = make_shared<CompSphere>(parent, r);
	m_comps.push_back(comp);
	comp->setTransform(E);
	comp->load(RESOURCE_DIR);
	m_ncomps++;
	return comp;
}

shared_ptr<CompCylinder> World::addCompCylinder(double r, shared_ptr<Body> parent, Matrix4d E, Vector3d z, Vector3d o, const string &RESOURCE_DIR, string shape) {
	auto comp = make_shared<CompCylinder>(parent, r);
	m_comps.push_back(comp);
	comp->setTransform(E);
	auto z_axis = make_shared<Vector>();
	z_axis->dir0 = z;
	comp->setZAxis(z_axis);
	auto origin = make_shared<Node>();
	origin->x0 = o;
	comp->setOrigin(origin);
	comp->load(RESOURCE_DIR, shape);

	m_ncomps++;
	return comp;
}

shared_ptr<CompDoubleCylinder> World::addCompDoubleCylinder(
	double rA, shared_ptr<Body> parentA, Matrix4d EA, Vector3d z_a, Vector3d o_a,
	double rB, shared_ptr<Body> parentB, Matrix4d EB, Vector3d z_b, Vector3d o_b,
	const string &RESOURCE_DIR, string shapeA, string shapeB) {
	auto comp = make_shared<CompDoubleCylinder>(parentA, rA, parentB, rB);
	m_comps.push_back(comp);
	comp->setTransformA(EA);
	comp->setTransformB(EB);
	auto  za = make_shared<Vector>();
	za->dir0 = z_a;
	comp->setZAxisA(za);
	auto zb = make_shared<Vector>();
	zb->dir0 = z_b;
	comp->setZAxisB(zb);
	auto originA = make_shared<Node>();
	originA->x0 = o_a;
	auto originB = make_shared<Node>();
	originB->x0 = o_b;
	comp->setOriginA(originA);
	comp->setOriginB(originB);
	comp->load(RESOURCE_DIR, shapeA, shapeB);
	m_ncomps++;
	return comp;
}

shared_ptr<MeshEmbedding> World::addMeshEmbedding(
	double density,
	double young,
	double possion,
	Material material,
	const string &RESOURCE_DIR,
	const string &TETGEN_FLAGS_0,
	const string &TETGEN_FLAGS_1,
	string coarse_mesh,
	string dense_mesh,
	SoftBodyType soft_body_type)
{
	
	auto mesh_embedding = make_shared<MeshEmbedding>(density, young, possion, material, soft_body_type);
	mesh_embedding->load(RESOURCE_DIR, coarse_mesh, TETGEN_FLAGS_0, dense_mesh, TETGEN_FLAGS_1);

	m_nmeshembeddings++;
	m_meshembeddings.push_back(mesh_embedding);
	return mesh_embedding;
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

shared_ptr<DeformableNull> World::addDeformableNull() {

	auto deformable = make_shared<DeformableNull>();
	m_ndeformables++;
	m_deformables.push_back(deformable);
	return deformable;
}

shared_ptr<CompNull> World::addCompNull() {
	auto comp = make_shared<CompNull>();
	m_ncomps++;
	m_comps.push_back(comp);
	return comp;
}

shared_ptr<WrapObst> World::addWrapNull() {
	auto wrap = make_shared<WrapNull>();
	m_nwraps++;
	m_wraps.push_back(wrap);
	return wrap;

}

shared_ptr<SpringNull> World::addSpringNull() {
	auto spring = make_shared<SpringNull>();
	m_nsprings++;
	m_springs.push_back(spring);
	return spring;
}

shared_ptr<MeshEmbeddingNull> World::addMeshEmbeddingNull() {
	auto mesh_null = make_shared<MeshEmbeddingNull>();
	m_nmeshembeddings++;
	m_meshembeddings.push_back(mesh_null);
	return mesh_null;
}

shared_ptr<WrapSphere> World::addWrapSphere(shared_ptr<Body> b0, Vector3d r0, shared_ptr<Body> b1, Vector3d r1, shared_ptr<CompSphere> compSphere, int num_points, const string &RESOURCE_DIR) {
	auto P = make_shared<Node>();
	P->x0 = r0;
	P->setParent(b0);
	auto S = make_shared<Node>();
	S->x0 = r1;
	S->setParent(b1);

	auto wrapSphere = make_shared<WrapSphere>(P, S, compSphere, num_points);
	m_nwraps++;
	wrapSphere->load(RESOURCE_DIR);
	m_wraps.push_back(wrapSphere);
	return wrapSphere;
}

shared_ptr<WrapCylinder> World::addWrapCylinder(shared_ptr<Body> b0, Vector3d r0, shared_ptr<Body> b1, Vector3d r1, shared_ptr<CompCylinder> compCylinder, int num_points, const string &RESOURCE_DIR) {
	auto P = make_shared<Node>();
	P->x0 = r0;
	P->setParent(b0);

	auto S = make_shared<Node>();
	S->x0 = r1;
	S->setParent(b1);

	auto wrapCylinder = make_shared<WrapCylinder>(P, S, compCylinder, num_points);
	m_nwraps++;
	wrapCylinder->load(RESOURCE_DIR);
	m_wraps.push_back(wrapCylinder);

	return wrapCylinder;
}

shared_ptr<WrapDoubleCylinder> World::addWrapDoubleCylinder(shared_ptr<Body> b0, Vector3d r0, shared_ptr<Body> b1, Vector3d r1, Vector3d u, Vector3d v, Vector3d z_u, Vector3d z_v, shared_ptr<CompDoubleCylinder> compDoubleCylinder, int num_points, const string &RESOURCE_DIR) {
	auto z_axis_U = make_shared<Vector>();
	z_axis_U->dir0 = z_u;
	auto z_axis_V = make_shared<Vector>();
	z_axis_V->dir0 = z_v;
	compDoubleCylinder->setZAxisA(z_axis_U);
	compDoubleCylinder->setZAxisB(z_axis_V);

	auto origin_U = make_shared<Node>();
	origin_U->x0 = u;
	compDoubleCylinder->setOriginA(origin_U);
	auto origin_V = make_shared<Node>();
	origin_V->x0 = v;
	compDoubleCylinder->setOriginB(origin_V);

	auto P = make_shared<Node>();
	P->x0 = r0;
	P->setParent(b0);
	auto S = make_shared<Node>();
	S->x0 = r1;
	S->setParent(b1);
	auto wrapDoubleCylinder = make_shared<WrapDoubleCylinder>(P, S, compDoubleCylinder, num_points);
	wrapDoubleCylinder->load(RESOURCE_DIR);
	m_wraps.push_back(wrapDoubleCylinder);
	m_nwraps++;
	return wrapDoubleCylinder;
}

void World::addSkeleton(Vector3d x0, Vector3d x1, int nlines, int nsamples, int &idx_start_body, vector<shared_ptr<Node>> &nodes) {
	double len_skeleton = (x1 - x0).norm();
	double len_line = len_skeleton / nlines;
	Vector3d dir = (x1 - x0) / len_skeleton;
	for (int i = 0; i < nlines; ++i) {
		Vector3d n0 = x0 + i * len_line * dir;
		Vector3d n1 = n0 + len_line * dir;
		
		auto line = make_shared<Line>(n0, n1);
		line->setBody(m_bodies[i + idx_start_body]);
		m_lines.push_back(line);
		line->addSampleNodes(nsamples, nodes);
	}

	idx_start_body = idx_start_body + nlines;
}

void World::init() {
	for (int i = 0; i < m_nbodies; i++) {
		m_bodies[i]->init(nm);
		if (i < m_nbodies - 1) {
			m_bodies[i]->next = m_bodies[i + 1];
		}
	}
	//cout << m_bodies[m_nbodies - 1]->idxM << endl;

	nm = 0;
	/*for (int i = 0; i < m_njoints; i++) {
	m_joints[i]->init(nm, nr);
	}*/

	//joint ordering
	// todo
	for (int i = m_njoints - 1; i > -1; i--) {
		m_joints[i]->countHRDofs(nR);
	}

	for (int i = m_njoints - 1; i > -1; i--) {		
		m_joints[i]->init(nm, nr);
	}

	m_dense_nm = nm;
	m_dense_nr = nr;
	// Until now, save nm, nr for later Sparse Jacobian Computation, 
	// only dense_nm, dense_nr 


	for (int i = 0; i < m_njoints; i++) {
		if (i < m_njoints - 1) {
			m_joints[i]->next = m_joints[i + 1];
		}
		if (i > 0) {
			m_joints[i]->prev = m_joints[i - 1];
		}
	}

	for (int i = 0; i < m_ncomps; ++i) {
		m_comps[i]->init();
		if (i < m_ncomps - 1) {
			m_comps[i]->next = m_comps[i + 1];
		}
	}

	for (int i = 0; i < m_nsprings; ++i) {
		m_springs[i]->init();
		if (i < m_nsprings - 1) {
			m_springs[i]->next = m_springs[i + 1];
		}
	}

	for (int i = 0; i < m_nwraps; ++i) {
		m_wraps[i]->init();
		if (i < m_nwraps - 1) {
			m_wraps[i]->next = m_wraps[i + 1];
		}
	}

	if (m_njoints == 0) {
		addJointNull();
	}

	if (m_ncomps == 0) {
		addCompNull();
	}

	if (m_nwraps == 0) {
		addWrapNull();
	}

	if (m_nsprings == 0) {
		addSpringNull();
	}

	m_joints[0]->update();
	m_comps[0]->update();
	m_wraps[0]->update();
	m_springs[0]->update();
	
	for (int i = 0; i < m_ndeformables; i++) {
		m_deformables[i]->countDofs(nm, nr);

		m_deformables[i]->init();
		// Create attachment constraints
		auto constraint = make_shared<ConstraintAttachSpring>(m_deformables[i]);
		m_constraints.push_back(constraint);
		m_nconstraints++;
		if (i < m_ndeformables - 1) {
			m_deformables[i]->next = m_deformables[i + 1];
		}
	}

	
	if (m_ndeformables == 0) {
		addDeformableNull();
	}

	if (m_type == SOFT_BODIES) {

		//m_softbodies[0]->setAttachmentsByYZCircle(5.0, 0.0001, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		//m_softbodies[0]->setSlidingNodesByYZCircle(7.5, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		m_softbodies[0]->setSlidingNodesByYZCircle(7.5, 0.7, Vector2d(0.0, 0.0), 0.705, m_bodies[0]);

		m_softbodies[0]->setSlidingNodesByYZCircle(10, 0.7, Vector2d(0.0, 0.0), 0.705, m_bodies[0]);
		m_softbodies[0]->setSlidingNodesByYZCircle(15.0, 0.7, Vector2d(0.0, 0.0), 0.705, m_bodies[1]);

	}

	if (m_type == SOFT_BODIES_CUBE_INVERTIBLE || m_type == SOFT_BODIES_CUBE_COROTATIONAL_LINEAR) {
		m_softbodies[0]->setAttachmentsByYZCircle(5.0, 0.0001, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);

		m_softbodies[0]->setAttachmentsByYZCircle(7.5, 1.4, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		m_softbodies[0]->setAttachmentsByYZCircle(13.5, 2.5, Vector2d(0.0, 0.0), 0.5, m_bodies[1]);

		/*m_softbodies[0]->setSlidingNodesByYZCircle(7.5, 0.7,Vector2d(0.0, 0.0), 0.705, m_bodies[0]);

		m_softbodies[0]->setSlidingNodesByYZCircle(10, 0.7, Vector2d(0.0, 0.0), 0.705, m_bodies[0]);
		m_softbodies[0]->setSlidingNodesByYZCircle(15.0, 0.7, Vector2d(0.0, 0.0), 0.705, m_bodies[1]);
		m_softbodies[0]->setSlidingNodesByYZCircle(12.5, 0.7, Vector2d(0.0, 0.0), 0.705, m_bodies[1]);*/
	}

	if (m_type == SOFT_BODIES_CYLINDER_INVERTIBLE || m_type == SOFT_BODIES_CYLINDER_COROTATIONAL_LINEAR) {
		//m_softbodies[0]->setAttachmentsByYZCircle(5.0, 0.0001, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		m_softbodies[0]->setAttachmentsByYZCircle(7.5, 1.4, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		m_softbodies[0]->setAttachmentsByYZCircle(13.5, 2.5, Vector2d(0.0, 0.0), 0.5, m_bodies[1]);
		//m_softbodies[0]->setSlidingNodesByYZCircle(7.5, 1.5, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		//m_softbodies[0]->setSlidingNodesByYZCircle(12.5, 1.5, Vector2d(0.0, 0.0), 0.5, m_bodies[1]);
	}

	if (m_type == MESH_EMBEDDING) {
		m_meshembeddings[0]->setAttachmentsByYZCircle(7.5, 1.4, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		m_meshembeddings[0]->setAttachmentsByYZCircle(7.5, 1.4, Vector2d(0.0, 0.0), 0.8, m_bodies[0]);

		m_meshembeddings[0]->setAttachmentsByYZCircle(13.5, 2.5, Vector2d(0.0, 0.0), 0.5, m_bodies[1]);
		m_meshembeddings[0]->setAttachmentsByYZCircle(13.5, 2.5, Vector2d(0.0, 0.0), 0.8, m_bodies[1]);

	}

	if (m_type == HUMAN_BODY) {
		m_meshembeddings[0]->setAttachmentsByYZCircle(5.0, 4.5, Vector2d(0.0, 0.0), 0.5, m_bodies[0]);
		//m_meshembeddings[0]->setAttachmentsByYZCircle(5.0, 4.5, Vector2d(0.0, 0.0), 0.8, m_bodies[0]);
		m_meshembeddings[0]->setAttachmentsByYZCircle(0.0, 0.1, Vector2d(0.0, 0.0), 1.414*2, m_bodies[0]);
		m_meshembeddings[0]->setAttachmentsByYZCircle(20.2, 1.0, Vector2d(0.0, 0.0), 2.05, m_bodies[1]);

		m_meshembeddings[0]->setAttachmentsByYZCircle(15.0, 4.5, Vector2d(0.0, 0.0), 0.5, m_bodies[1]);
		//m_meshembeddings[0]->setAttachmentsByYZCircle(15.0, 4.5, Vector2d(0.0, 0.0), 0.8, m_bodies[1]);
		//m_meshembeddings[1]->setAttachmentsByXZCircle(0.0, 4.5, Vector2d(-7.5, 0.0), 0.5, m_bodies[2]);
	//	m_meshembeddings[1]->setAttachmentsByXZCircle(-10.0, 4.5, Vector2d(-7.5, 0.0), 0.5, m_bodies[3]);

		m_meshembeddings[1]->setAttachmentsByYZCircle(-20.0, 4.5, Vector2d(0.0, 0.0), 0.5, m_bodies[5]);
		m_meshembeddings[1]->setAttachmentsByYZCircle(-30.0, 4.5, Vector2d(0.0, 0.0), 0.5, m_bodies[6]);
		m_meshembeddings[1]->setAttachmentsByYZCircle(-15.0, 1.0, Vector2d(0.0, 0.0), 2.05, m_bodies[5]);
		m_meshembeddings[1]->setAttachmentsByYZCircle(-35.0, 1.0, Vector2d(0.0, 0.0), 2.05, m_bodies[6]);

		if (isrightleg) {
			m_meshembeddings[2]->setAttachmentsByXZCircle(-18.0, 3.0, Vector2d(-12.0, 0.0), 0.5, m_bodies[8]);
			m_meshembeddings[2]->setAttachmentsByXZCircle(-28.0, 3.0, Vector2d(-12.0, 0.0), 0.5, m_bodies[9]);
		}
		
		if(isrightleg){
			m_meshembeddings[3]->setAttachmentsByXZCircle(-18.0, 3.0, Vector2d(-3.0, 0.0), 0.5, m_bodies[11]);
			m_meshembeddings[3]->setAttachmentsByXZCircle(-28.0, 3.0, Vector2d(-3.0, 0.0), 0.5, m_bodies[12]);
		}
		
	}


	if (m_type == WORM) {

		for (int i = 0; i < (int)m_lines.size(); ++i) {
			m_meshembeddings[0]->setAttachmentsByLine(m_lines[i]);
		}

	}

	if (m_type == CROSS) {	
		for (int i = 0; i < (int)m_lines.size(); ++i) {
			m_meshembeddings[0]->setAttachmentsByLine(m_lines[i]);
		}
	}

	if (m_type == STARFISH) {
		for (int i = 0; i < (int)m_lines.size(); ++i) {
			m_meshembeddings[0]->setAttachmentsByLine(m_lines[i]);
		}
		//for (int i = 0; i < nlegs; i++) {
		//	double r_ =5.0;
		//	for (int j = 0; j < nsegments; j++) {
		//		double theta = (18.0 + i * 72.0) / 180.0 * M_PI;
		//		double x_, y_, z_;
		//		x_ = (j * 20.0 + 10.0) * cos(theta);
		//		y_ = - 0.5 * j;
		//		z_ = - (j * 20.0 + 10.0) * sin(theta);
		//		
		//		m_meshembeddings[0]->setAttachmentsByYZCircle(x_, 9.5, Vector2d(y_, z_), r_, m_bodies[nsegments * i + j]);
		//		//r_ *= 0.8;

		//	}
		//}
	}

	if (m_type == STARFISH_2) {
		for (int i = 0; i < (int)m_lines.size(); ++i) {
			m_meshembeddings[0]->setAttachmentsByLine(m_lines[i]);
		}
	}
	if (m_type == STARFISH3) {
		for (int i = 0; i < (int)m_lines.size(); ++i) {
			m_softbodies[0]->setAttachmentsByLine(m_lines[i]);
		}
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

	int tet = nm;
	for (int i = 0; i < m_nmeshembeddings; i++) {
		m_meshembeddings[i]->countDofs(nm, nr);
		m_meshembeddings[i]->init();
		// Create attachment constraints
		auto constraint = make_shared<ConstraintAttachSoftBody>(m_meshembeddings[i]->getCoarseMesh());
		m_constraints.push_back(constraint);
		m_nconstraints++;

		if (i < m_nmeshembeddings - 1) {
			m_meshembeddings[i]->next = m_meshembeddings[i + 1];
		}
	}

	tet = nm - tet;
	tet = tet / 3;
	m_ntets = tet;

	if (m_nmeshembeddings == 0) {
		auto mesh_null = make_shared<MeshEmbeddingNull>();
		m_meshembeddings.push_back(mesh_null);
		m_nmeshembeddings++;
	}

	// init constraints

	for (int i = 0; i < m_nconstraints; i++) {
		m_constraints[i]->countDofs(nem, ner, nim, nir);
		m_constraints[i]->init();
		if (i < m_nconstraints - 1) {
			m_constraints[i]->next = m_constraints[i + 1];
		}
	}

	if (m_nconstraints == 0) {
		addConstraintNull();
	}
}

void World::update() {

	m_comps[0]->update();
	m_wraps[0]->update();
	m_springs[0]->update();
}

int World::getNsteps() {
	// Computes the number of results
	int nsteps = int((m_tspan(1) - m_tspan(0)) / m_h);
	return nsteps;
}

void World::drawFloor(Floor f, shared_ptr<MatrixStack> MV, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) {
	// Draw grid
	progSimple->bind();
	glUniformMatrix4fv(progSimple->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glLineWidth(2.0f);
	float x0 = f.xrange(0);
	float x1 = f.xrange(1);
	float z0 = f.zrange(0);
	float z1 = f.zrange(1);
	int gridSize = 10;
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	for (int i = 1; i < gridSize; ++i) {
		if (i == gridSize / 2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		}
		else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float x = x0 + i / (float)gridSize * (x1 - x0);
		glVertex3f(x, f.y, z0);
		glVertex3f(x, f.y, z1);
	}
	for (int i = 1; i < gridSize; ++i) {
		if (i == gridSize / 2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		}
		else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float z = z0 + i / (float)gridSize * (z1 - z0);
		glVertex3f(x0, f.y, z);
		glVertex3f(x1, f.y, z);
	}
	glEnd();
	glColor3f(0.4f, 0.4f, 0.4f);
	glBegin(GL_LINE_LOOP);
	glVertex3f(x0, f.y, z0);
	glVertex3f(x1, f.y, z0);
	glVertex3f(x1, f.y, z1);
	glVertex3f(x0, f.y, z1);
	glEnd();
	progSimple->unbind();
}

void World::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, const shared_ptr<Program> progSoft, shared_ptr<MatrixStack> P) {
	m_bodies[0]->draw(MV, prog, P);
	m_joints[0]->draw(MV, prog, progSimple, P);
	m_deformables[0]->draw(MV, prog, progSimple, P);
	m_comps[0]->draw(MV, prog, P);
	m_wraps[0]->draw(MV, prog, progSimple, P);
	m_springs[0]->draw(MV, prog, progSimple, P);
	m_softbodies[0]->draw(MV, prog, progSimple, P);
	m_meshembeddings[0]->draw(MV, prog, progSimple, P);
	for (int i = 0; i < (int)m_floors.size(); ++i) {
		drawFloor(m_floors[i], MV, progSimple, P);
	}
}

Energy World::computeEnergy() {
	m_energy.K = 0.0;
	m_energy.V = 0.0;

	m_joints[0]->computeEnergies(m_grav, m_energy);
	m_deformables[0]->computeEnergies(m_grav, m_energy);
	m_energy = m_softbodies[0]->computeEnergies(m_grav, m_energy);

	if (m_t == 0.0) {
		m_energy0 = m_energy;
	}

	m_energy.V -= m_energy0.V;

	return m_energy;
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

void World::sceneCross(double t) {
	//m_joints[0]->presc->setActive();
	//m_joints[0]->presc->m_q[0] = 0.0;
	//m_joints[0]->presc->m_qdot[0] = 0.0;
	//m_joints[0]->presc->m_qddot[0] = 0.0;
}

void World::sceneStarFish(double t) {

	double sinTheta =  M_PI * sin( t);
	double cosTheta = M_PI * cos( t);
	double d30 = 1.0 / 6.0;
	double d45 = 1.0 / 4.0;
	double d15 = 1.0 / 12.0;
	double d60 = 1.0 / 3.0;
	double d90 = 1.0 / 2.0;
	double d10 = 1.0 / 18.0;
	double d12 = 1.0 / 15.0;
	double d18 = 1.0 / 10.0;
	double d20 = 1.0 / 9.0;
	double d22 = 1.0 / 8.0;
	//
	//for (int i = 0; i < nlegs; i++) {
	//	if (i > -1) {
	//		//m_joints[0 + nsegments * i]->presc->m_q[0] = d60 * sinTheta;
	//		//m_joints[0 + nsegments * i]->presc->m_qdot[0] = d60 * cosTheta;
	//		//m_joints[0 + nsegments * i]->presc->m_qddot[0] = -d60 * sinTheta;


	//		m_joints[1 + nsegments * i]->presc->m_q[0] = -d30 * sinTheta;
	//		m_joints[1 + nsegments * i]->presc->m_qdot[0] = -d30 * cosTheta;
	//		m_joints[1 + nsegments * i]->presc->m_qddot[0] = d30 * sinTheta;

	//		m_joints[2 + nsegments * i]->presc->m_q[0] = -d30 * sinTheta;
	//		m_joints[2 + nsegments * i]->presc->m_qdot[0] = -d30 * cosTheta;
	//		m_joints[2 + nsegments * i]->presc->m_qddot[0] = d30 * sinTheta;

	//		m_joints[3 + nsegments * i]->presc->m_q[0] = -d22 * sinTheta;
	//		m_joints[3 + nsegments * i]->presc->m_qdot[0] = -d22 * cosTheta;
	//		m_joints[3 + nsegments * i]->presc->m_qddot[0] = d22 * sinTheta;

	//		m_joints[4 + nsegments * i]->presc->m_q[0] = -d30 * sinTheta;
	//		m_joints[4 + nsegments * i]->presc->m_qdot[0] = -d30 * cosTheta;
	//		m_joints[4 + nsegments * i]->presc->m_qddot[0] = d30 * sinTheta;


	//		m_joints[5 + nsegments * i]->presc->m_q[0] = -d30 * sinTheta;
	//		m_joints[5 + nsegments * i]->presc->m_qdot[0] = -d30 * cosTheta;
	//		m_joints[5 + nsegments * i]->presc->m_qddot[0] = d30 * sinTheta;

	//		m_joints[6 + nsegments * i]->presc->m_q[0] = -d30 * sinTheta;
	//		m_joints[6 + nsegments * i]->presc->m_qdot[0] = -d30 * cosTheta;
	//		m_joints[6 + nsegments * i]->presc->m_qddot[0] = d30 * sinTheta;

	//		m_joints[7 + nsegments * i]->presc->m_q[0] = -d30 * sinTheta;
	//		m_joints[7 + nsegments * i]->presc->m_qdot[0] = -d30 * cosTheta;
	//		m_joints[7 + nsegments * i]->presc->m_qddot[0] = d30 * sinTheta;

	//	}
	//	else {
	//		m_joints[0 + 8 * i]->presc->m_q[0] = -1.0 / 6.0 * sinTheta;
	//		m_joints[0 + 8 * i]->presc->m_qdot[0] = -1.0 / 6.0 * cosTheta;
	//		m_joints[0 + 8 * i]->presc->m_qddot[0] = 1.0 / 6.0 * sinTheta;


	//		m_joints[1 + 8 * i]->presc->m_q[0] = 1.0 / 18.0 * sinTheta;
	//		m_joints[1 + 8 * i]->presc->m_qdot[0] = 1.0 / 18.0 * cosTheta;
	//		m_joints[1 + 8 * i]->presc->m_qddot[0] = -1.0 / 18.0 * sinTheta;

	//		m_joints[2 + 8 * i]->presc->m_q[0] = 1.0 / 18.0 * sinTheta;
	//		m_joints[2 + 8 * i]->presc->m_qdot[0] = 1.0 / 18.0 * cosTheta;
	//		m_joints[2 + 8 * i]->presc->m_qddot[0] = -1.0 /18.0 * sinTheta;

	//		m_joints[3 + 8 * i]->presc->m_q[0] = 1.0 / 15.0 * sinTheta;
	//		m_joints[3 + 8 * i]->presc->m_qdot[0] = 1.0 / 15.0 * cosTheta;
	//		m_joints[3 + 8 * i]->presc->m_qddot[0] = -1.0 / 15.0 * sinTheta;

	//		m_joints[4 + 8 * i]->presc->m_q[0] = 1.0 / 12.0 * sinTheta;
	//		m_joints[4 + 8 * i]->presc->m_qdot[0] = 1.0 / 12.0 * cosTheta;
	//		m_joints[4 + 8 * i]->presc->m_qddot[0] = -1.0 / 12.0 * sinTheta;


	//		m_joints[5 + 8 * i]->presc->m_q[0] =1.0 / 15.0 * sinTheta;
	//		m_joints[5 + 8 * i]->presc->m_qdot[0] = 1.0 / 15.0 * cosTheta;
	//		m_joints[5 + 8 * i]->presc->m_qddot[0] =- 1.0 / 15.0 * sinTheta;

	//		m_joints[6 + 8 * i]->presc->m_q[0] = 1.0 / 18.0 * sinTheta;
	//		m_joints[6 + 8 * i]->presc->m_qdot[0] = 1.0 /18.0 * cosTheta;
	//		m_joints[6 + 8 * i]->presc->m_qddot[0] =- 1.0 / 18.0 * sinTheta;

	//		m_joints[7 + 8 * i]->presc->m_q[0] = 1.0 / 18.0 * sinTheta;
	//		m_joints[7 + 8 * i]->presc->m_qdot[0] = 1.0 / 18.0 * cosTheta;
	//		m_joints[7 + 8 * i]->presc->m_qddot[0] =-1.0 / 18.0 * sinTheta;

	//	}	
	//}
}

void World::sceneStarFish2(double t) {
	cout << t << endl;
	double d30 = -1.0 / 6.0;
	double d45 = -1.0 / 4.0;
	double d15 = -1.0 / 12.0;
	double d60 = -1.0 / 3.0;
	double d90 = -1.0 / 2.0;
	double d10 = -1.0 / 18.0;
	double d12 = -1.0 / 15.0;
	double d18 = -1.0 / 10.0;
	double d20 = -1.0 / 9.0;
	double d22 = -1.0 / 8.0;
	Vector3d vt_w, wt_i, vtdot_w, wtdot_i;
	
	vt_w.setZero();
	vtdot_w.setZero();
	wt_i.setZero();
	wtdot_i.setZero();
	Vector3d Zero = Vector3d::Zero();

	double alpha = 1.0 / 7.9617 / 2.0;

	double beta = 1.4;
	double q, dq;

	if (t < 6.0) {
		// Hold one arm still and let the starfish falls down to floor and collides with floor
		setMaximalPrescStates(m_bodies[8], Zero, Zero, Zero, Zero);

		// Active:
		// mcon: 8
		// rcon: null
	}
	else if (t < 16.0) {
		// Lift an arm of the starfish to make it off the floor, velocity increases

		double t_ = t - 6.0;
		vt_w << 0.0, beta * t_, 0.0;
		vtdot_w << 0.0, beta, 0.0;

		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, Zero, Zero);

		// Active:
		// mcon: 8
		// rcon: null

	}
	else if (t < 26.0) {
		// Still lifting, velocity decreases to zero

		double t_ = t - 26.0; // t_ = [-10, 0]
		vt_w << 0.0, -beta * t_, 0.0; 
		vtdot_w << 0.0, -beta, 0.0;
		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, Zero, Zero);

		// Active:
		// mcon: 8
		// rcon: null
	}
	else if (t < 100.0) {
		// Starfish starts to wiggle and then stop wiggling

		// Hold one point
		setMaximalPrescStates(m_bodies[8], Zero, Zero, Zero, Zero);
		
		double t0 = 26.0;
		double t_ = t - t0; // t_ = [0, 20]
		double t1 = 100.0;

		//double alpha = 6 * M_PI / (t1 - t0);
		double alpha = 1.0;
		double sinTheta = M_PI * sin(alpha * t_);
		double cosTheta = M_PI * cos(alpha * t_);
		
		computeTargetQ(t0, t1, t, 0, 0, q, dq);
		setReducedPrescStates(m_joints[13], q, dq);

		//setReducedPrescStates(m_joints[5], -d30 * sinTheta, -d30 * alpha * cosTheta);
		//setReducedPrescStates(m_joints[9], d20 * sinTheta, d20 * alpha * cosTheta);
		//setReducedPrescStates(m_joints[17], d30 * sinTheta, d30 * alpha * cosTheta);
		setReducedPrescStates(m_joints[5], -d30 * sinTheta, -d30 * alpha * cosTheta);
		setReducedPrescStates(m_joints[9], d30 * sinTheta, d30 * alpha * cosTheta);
		setReducedPrescStates(m_joints[17], d30 * sinTheta, d30 * alpha * cosTheta);
		// Active:
		// mcon: 8
		// ncon: 5 9 13 17 
	}
	else if (t < 100.0) {
		// Settle down on the floor again
		VectorXi mcon(1);	
		Vector4i ncon;
		mcon << 8;
		ncon << 5, 9, 13, 17;
		deactivateListOfPrescConstraints(mcon, ncon);

		// Active:
		// mcon: null
		// ncon: null
	}else if (t < 110) {

		// Try to flip around 

		double t0 = 100.0;
		double t_ = t - t0;	// t_ = [0, 10]
		double t1 = 110.0;

		computeTargetQ(t0, t1, t, d45 * M_PI, 0, q, dq);
		setReducedPrescStates(m_joints[1], q, dq);
		computeTargetQ(t0, t1, t, d45 * M_PI, 0, q, dq);
		setReducedPrescStates(m_joints[2], q, dq);
		computeTargetQ(t0, t1, t, d45 * M_PI, 0, q, dq);
		setReducedPrescStates(m_joints[3], q, dq);
		computeTargetQ(t0, t1, t, d45 * M_PI, 0, q, dq);
		setReducedPrescStates(m_joints[4], q, dq);

		computeTargetQ(t0, t1, t, d45 * M_PI, 0, q, dq);
		setReducedPrescStates(m_joints[6], q, dq);
		computeTargetQ(t0, t1, t, d45 * M_PI, 0, q, dq);
		setReducedPrescStates(m_joints[7], q, dq);

		double p = 0.1;
		vt_w << 0.0, p * t_, 0.0;
		vtdot_w << 0.0, p, 0.0;
		setMaximalPrescStates(m_bodies[0], vt_w, vtdot_w, Zero, Zero);

		// Active: 
		// mcon: 0
		// ncon: 1 2 3 4 6 7
	}
	else if (t < 120.0) {
		// Flip around

		VectorXi ncon(1);
		VectorXi mcon(1);
		ncon << -2;
		mcon << 0;
		deactivateListOfPrescConstraints(mcon, ncon);

		double t0 = 110.0;
		double t_ = t - t0;	// t_ = [0, 10]
		double t1 = 120.0;


		// Drag the highest body in x +
		double p = 0.2;
		vt_w << p * t_, p * t_, 0.0;
		vtdot_w << p, p, 0.0;

		computeTargetQ(t0, t1, t, -d45 * M_PI, d45 * M_PI, q, dq);
		setReducedPrescStates(m_joints[1], q, dq);
		computeTargetQ(t0, t1, t, -d45 * M_PI, d45 * M_PI, q, dq);
		setReducedPrescStates(m_joints[2], q, dq);
		computeTargetQ(t0, t1, t, -d45 * M_PI, d45 * M_PI, q, dq);
		setReducedPrescStates(m_joints[3], q, dq);
		computeTargetQ(t0, t1, t, -d45 * M_PI, d45 * M_PI, q, dq);
		setReducedPrescStates(m_joints[6], q, dq);
		//computeTargetQ(t0, t1, t, -d45 * M_PI, d45 * M_PI, q, dq);
		setReducedPrescStates(m_joints[7], q, dq);


		setMaximalPrescStates(m_bodies[4], vt_w, vtdot_w, Zero, Zero);

		// Fix one end so it won't be dragged away

		p = 0.3;
		vt_w << -p * t_, 0.0, 0.0;
		vtdot_w << -p, 0.0, 0.0;
		setMaximalPrescStates(m_bodies[12], Zero, Zero, Zero, Zero);
		//setMaximalPrescStates(m_bodies[1], vt_w, vtdot_w, Zero, Zero);


		// Active: 
		// mcon: 4 12
		// ncon: 1 2 3 4 6 7 
	}
	else if(t < 130.0)
	{
		double t0 = 120.0;
		double t_ = t - t0;	// t_ = [0, 10]
		double t1 = 130.0;

		VectorXi ncon(5);
		VectorXi mcon(2);
		ncon << 1, 2, 3, 6, 7;
		mcon << 4, 12;
		deactivateListOfPrescConstraints(mcon, ncon);

		double p = 0.02;
		wt_i << 0.0, 0.0, -p * t_;
		wtdot_i << 0.0, 0.0, -p;

		setMaximalPrescStates(m_bodies[0], Zero, Zero, wt_i, wtdot_i);


		computeTargetQ(t0, t1, t, -d45 * M_PI, d45 * M_PI, q, dq);
		setReducedPrescStates(m_joints[4], q, dq);
		// Active:
		// mcon: 0 
		// ncon: 4

	}
	else if (t < 135.0) {
		VectorXi ncon(1);
		VectorXi mcon(1);
		ncon << 4;
		mcon << 0;
		deactivateListOfPrescConstraints(mcon, ncon);

		double t0 = 130.0;
		double t_ = t - t0;	// t_ = [0, 10]
		double t1 = 135.0;



/*
		setReducedPrescStates(m_joints[1], q, dq);
		setReducedPrescStates(m_joints[9], q, dq);*/

		//computeTargetQ(t0, t1, t, -d20 * M_PI, 0, q, dq);

		//setReducedPrescStates(m_joints[5], q, dq);
		//setReducedPrescStates(m_joints[13], q, dq);

		setMaximalPrescStates(m_bodies[12], Zero, Zero, Zero, Zero);
		setMaximalPrescStates(m_bodies[20], Zero, Zero, Zero, Zero);

		// Active:
		// mcon: 12 20
		// ncon: null

	}
	else if (t < 150.0) {
		VectorXi ncon(1);
		VectorXi mcon(2);
		ncon << -2;
		mcon << 12, 20;
		deactivateListOfPrescConstraints(mcon, ncon);

		// Active:
		// mcon: null
		// ncon: null
	}
	else {

	}
}

void World::sceneStarFishJump(double t) {
	cout << t << endl;
	double d30 = -1.0 / 6.0;
	double d45 = -1.0 / 4.0;
	double d15 = -1.0 / 12.0;
	double d60 = -1.0 / 3.0;
	double d50 = -5.0 / 18.0;
	double d90 = -1.0 / 2.0;
	double d10 = -1.0 / 18.0;
	double d12 = -1.0 / 15.0;
	double d18 = -1.0 / 10.0;
	double d20 = -1.0 / 9.0;
	double d22 = -1.0 / 8.0;
	Vector3d vt_w, wt_i, vtdot_w, wtdot_i;
	double q, dq;

	vt_w.setZero();
	vtdot_w.setZero();
	wt_i.setZero();
	wtdot_i.setZero();
	double beta = 5.0;
	double gamma = 4.0;
	double sigma = 0.25;
	if (t < 5.0) {
		double t0 = 0.0;
		double t_ = t - t0;
		double t1 = 5.0;

		double alpha = M_PI / (t1 - t0);
		double sinTheta = M_PI * sin(alpha * t_);
		double cosTheta = M_PI * cos(alpha * t_);

		m_bodies[0]->presc->setActive();

		vt_w << -gamma * sinTheta, beta * sinTheta, 0.0;
		vtdot_w << - gamma * alpha * cosTheta, beta * alpha * cosTheta, 0.0;
		wt_i << 0.0, sigma * sinTheta, 0.0;
		wtdot_i << 0.0, sigma * alpha *cosTheta, 0.0;

		//wt_i.setZero();
		//wtdot_i.setZero();
		setMaximalPrescStates(m_bodies[0], vt_w, vtdot_w, wt_i, wtdot_i);

		m_bodies[4]->presc->setActive();
		m_bodies[8]->presc->setActive();
		m_bodies[16]->presc->setActive();
		m_bodies[20]->presc->setActive();
		m_bodies[12]->presc->setActive();

		vt_w.setZero();
		vtdot_w.setZero();
		wt_i << 0.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 0.0;

		//wt_i.setZero();
		//wtdot_i.setZero();

		setMaximalPrescStates(m_bodies[4], vt_w, vtdot_w, wt_i, wtdot_i);
		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);
		setMaximalPrescStates(m_bodies[16], vt_w, vtdot_w, wt_i, wtdot_i);
		setMaximalPrescStates(m_bodies[20], vt_w, vtdot_w, wt_i, wtdot_i);
		setMaximalPrescStates(m_bodies[12], vt_w, vtdot_w, wt_i, wtdot_i);




		//m_joints[1]->presc->setActive();
		//m_joints[5]->presc->setActive();
		//m_joints[9]->presc->setActive();
		//m_joints[13]->presc->setActive();
		//m_joints[17]->presc->setActive();

		//q = d30 * sinTheta;
		//dq = d30 * alpha * cosTheta;

		////computeTargetQ(t0, t1, t, d45 * M_PI, 0.0, q, dq);

		//m_joints[1]->presc->m_q[0] = q;
		//m_joints[1]->presc->m_qdot[0] = dq;

		//m_joints[5]->presc->m_q[0] = q;
		//m_joints[5]->presc->m_qdot[0] = dq;

		//m_joints[9]->presc->m_q[0] = q;
		//m_joints[9]->presc->m_qdot[0] = dq;

		//m_joints[13]->presc->m_q[0] = q;
		//m_joints[13]->presc->m_qdot[0] = dq;

		//m_joints[17]->presc->m_q[0] = q;
		//m_joints[17]->presc->m_qdot[0] = dq;

		//m_joints[2]->presc->setActive();
		//m_joints[6]->presc->setActive();
		//m_joints[10]->presc->setActive();
		//m_joints[14]->presc->setActive();
		//m_joints[18]->presc->setActive();

		////computeTargetQ(t0, t1, t, d30 * M_PI, 0.0, q, dq);



		//m_joints[2]->presc->m_q[0] = q;
		//m_joints[2]->presc->m_qdot[0] = dq;
		//m_joints[6]->presc->m_q[0] = q;
		//m_joints[6]->presc->m_qdot[0] = dq;		

		//m_joints[10]->presc->m_q[0] = d20 * sinTheta;
		//m_joints[10]->presc->m_qdot[0] = d20 * alpha * cosTheta;


		//m_joints[14]->presc->m_q[0] = q;
		//m_joints[14]->presc->m_qdot[0] = dq;


		//m_joints[18]->presc->m_q[0] = q;
		//m_joints[18]->presc->m_qdot[0] = dq;

		//m_joints[3]->presc->setActive();
		//m_joints[7]->presc->setActive();
		//m_joints[11]->presc->setActive();
		//m_joints[15]->presc->setActive();
		//m_joints[19]->presc->setActive();

		////computeTargetQ(t0, t1, t, d10 * M_PI, 0.0, q, dq);

		//m_joints[3]->presc->m_q[0] = q;
		//m_joints[3]->presc->m_qdot[0] = dq;
		//m_joints[7]->presc->m_q[0] = q;
		//m_joints[7]->presc->m_qdot[0] = dq;
		//m_joints[11]->presc->m_q[0] = q;
		//m_joints[11]->presc->m_qdot[0] = dq;
		//m_joints[15]->presc->m_q[0] = q;
		//m_joints[15]->presc->m_qdot[0] = dq;
		//m_joints[19]->presc->m_q[0] = q;
		//m_joints[19]->presc->m_qdot[0] = dq;

		//m_joints[4]->presc->setActive();
		//m_joints[8]->presc->setActive();
		//m_joints[12]->presc->setActive();
		//m_joints[16]->presc->setActive();
		//m_joints[20]->presc->setActive();
		////computeTargetQ(t0, t1, t, d10 * M_PI, 0.0, q, dq);

		//m_joints[4]->presc->m_q[0] = q;
		//m_joints[4]->presc->m_qdot[0] = dq;

		//m_joints[8]->presc->m_q[0] = q;
		//m_joints[8]->presc->m_qdot[0] = dq;
		//m_joints[12]->presc->m_q[0] = q;
		//m_joints[12]->presc->m_qdot[0] = dq;
		//m_joints[16]->presc->m_q[0] = q;
		//m_joints[16]->presc->m_qdot[0] = dq;
		//m_joints[20]->presc->m_q[0] = q;
		//m_joints[20]->presc->m_qdot[0] = dq;



	}
	else if (t < 10.0) {

		m_bodies[4]->presc->setInactive();
		m_bodies[8]->presc->setInactive();
		m_bodies[16]->presc->setInactive();
		m_bodies[20]->presc->setInactive();
		m_bodies[12]->presc->setInactive();


		double t0 = 5.0;
		double t1 = 10.0;
		double t_ = t - t0; 
		
		double alpha = M_PI / (t1 - t0);
		double sinTheta = M_PI * sin(alpha * t_);
		double cosTheta = M_PI * cos(alpha * t_);

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				m_joints[1 + i*4 + j]->presc->m_q[0] = 0;
				m_joints[1 + i * 4 + j]->presc->m_qdot[0] = 0;
				//m_joints[5 + i]->presc->m_qdot[0] = 0;
				//m_joints[9 + i]->presc->setInactive();
				//m_joints[13 + i]->presc->setInactive();
				//m_joints[17 + i]->presc->setInactive();
			}

		}

		vt_w << -gamma * sinTheta, - 1.3* beta * sinTheta, 0.0;
		vtdot_w << -gamma * alpha * cosTheta, - 1.3 *beta * alpha * cosTheta, 0.0;
		wt_i << 0.0, sigma * sinTheta, 0.0;
		wtdot_i << 0.0, sigma * alpha *cosTheta, 0.0;
		//m_bodies[0]->presc->setInactive();
		setMaximalPrescStates(m_bodies[0], vt_w, vtdot_w, wt_i, wtdot_i);


		/*vt_w << 0.0, beta * sinTheta, 0.0;
		vtdot_w << 0.0, beta * cosTheta, 0.0;
		wt_i.setZero();
		wtdot_i.setZero();
		setMaximalPrescStates(m_bodies[0], vt_w, vtdot_w, wt_i, wtdot_i);

		computeTargetQ(t0, t1, t, -d45 * M_PI, d45 * M_PI, q, dq);
		m_joints[1]->presc->m_q[0] = q;
		m_joints[1]->presc->m_qdot[0] = dq;
		m_joints[5]->presc->m_q[0] = q;
		m_joints[5]->presc->m_qdot[0] = dq;
		m_joints[9]->presc->m_q[0] = q;
		m_joints[9]->presc->m_qdot[0] = dq;
		m_joints[13]->presc->m_q[0] = q;
		m_joints[13]->presc->m_qdot[0] = dq;
		m_joints[17]->presc->m_q[0] = q;
		m_joints[17]->presc->m_qdot[0] = dq;

		computeTargetQ(t0, t1, t, -d30 * M_PI, d30 * M_PI, q, dq);
		m_joints[2]->presc->m_q[0] = q;
		m_joints[2]->presc->m_qdot[0] = dq;
		m_joints[6]->presc->m_q[0] = q;
		m_joints[6]->presc->m_qdot[0] = dq;
		m_joints[10]->presc->m_q[0] = q;
		m_joints[10]->presc->m_qdot[0] = dq;
		m_joints[14]->presc->m_q[0] = q;
		m_joints[14]->presc->m_qdot[0] = dq;
		m_joints[18]->presc->m_q[0] = q;
		m_joints[18]->presc->m_qdot[0] = dq;

		computeTargetQ(t0, t1, t, -d10 * M_PI, d10 * M_PI, q, dq);
		m_joints[3]->presc->m_q[0] = q;
		m_joints[3]->presc->m_qdot[0] = dq;
		m_joints[7]->presc->m_q[0] = q;
		m_joints[7]->presc->m_qdot[0] = dq;
		m_joints[11]->presc->m_q[0] = q;
		m_joints[11]->presc->m_qdot[0] = dq;
		m_joints[15]->presc->m_q[0] = q;
		m_joints[15]->presc->m_qdot[0] = dq;
		m_joints[19]->presc->m_q[0] = q;
		m_joints[19]->presc->m_qdot[0] = dq;*/

	}
	else if (t < 15.0) {
		double t0 = 10.0;
		double t1 = 15.0;
		double t_ = t - t0;
		/*vt_w << 0.0, 0.0, 0.0;
		vtdot_w << 0.0, 0.0, 0.0;
		wt_i << 0.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 0.0;
		setMaximalPrescStates(m_bodies[0], vt_w, vtdot_w, wt_i, wtdot_i);*/
		m_bodies[0]->presc->setInactive();


		//m_bodies[4]->presc->setActive();
		//m_bodies[8]->presc->setActive();
		//m_bodies[16]->presc->setActive();
		//m_bodies[20]->presc->setActive();
		//m_bodies[12]->presc->setActive();

		vt_w << 0.0, 10 * beta * t_, 0.0;
		vtdot_w << 0.0, 10 * beta, 0.0;
		wt_i << 0.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 0.0;

		//wt_i.setZero();
		//wtdot_i.setZero();

		//setMaximalPrescStates(m_bodies[4], vt_w, vtdot_w, wt_i, wtdot_i);
		//setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);
		//setMaximalPrescStates(m_bodies[16], vt_w, vtdot_w, wt_i, wtdot_i);
		//setMaximalPrescStates(m_bodies[20], vt_w, vtdot_w, wt_i, wtdot_i);
		//setMaximalPrescStates(m_bodies[12], vt_w, vtdot_w, wt_i, wtdot_i);


		//m_bodies[0]->presc->setInactive();
		//setMaximalPrescStates(m_bodies[0], vt_w, vtdot_w, wt_i, wtdot_i);
		vt_w << 0.0, -beta * t_, 0.0;
		vtdot_w << 0.0, -beta, 0.0;
		wt_i.setZero();
		wtdot_i.setZero();
		//setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);

		// vt_w = 0 0 0
	}
	else if (t < 20.0) {
		// Starfish starts to wiggle and then stop wiggling
		double t0 = 15.0;
		double t_ = t - t0; // t_ = [0, 20]

		double t1 = 20.0;

		double alpha = 6 * M_PI / (t1 - t0);
		double sinTheta = M_PI * sin(alpha * t_);
		double cosTheta = M_PI * cos(alpha * t_);
		vt_w << 0.0, beta * t_, 0.0;
		vtdot_w << 0.0, beta, 0.0;
		wt_i << 0.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 0.0;

		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);

		/*m_joints[1]->presc->setActive();
		m_joints[5]->presc->setActive();
		m_joints[9]->presc->setActive();
		m_joints[13]->presc->setActive();
		m_joints[17]->presc->setActive();

		m_joints[5]->presc->m_q[0] = -d20 * sinTheta;
		m_joints[5]->presc->m_qdot[0] = -d20 * alpha * cosTheta;
		m_joints[5]->presc->m_qddot[0] = d20 * alpha * alpha * sinTheta;

		m_joints[9]->presc->m_q[0] = d20 * sinTheta;
		m_joints[9]->presc->m_qdot[0] = d20 * alpha * cosTheta;
		m_joints[9]->presc->m_qddot[0] = -d20 * alpha * alpha * sinTheta;

		m_joints[13]->presc->m_q[0] = 0.0;
		m_joints[13]->presc->m_qdot[0] = 0.0;
		m_joints[13]->presc->m_qddot[0] = 0.0;

		m_joints[17]->presc->m_q[0] = d20 * sinTheta;
		m_joints[17]->presc->m_qdot[0] = d20 * alpha * cosTheta;
		m_joints[17]->presc->m_qddot[0] = -d20 * alpha * alpha * sinTheta;*/

	}
	else if (t < 25.0) {
		// Use two handles to put it down on the floor again

		double t_ = t - 20.0;	// t_ = [0, 5]


		vt_w << 0.0, beta * t_, 0.0;
		vtdot_w << 0.0, beta, 0.0;
		wt_i << 0.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 0.0;

		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);
		/*vt_w << 0.0, -beta * t_, 0.0;
		vtdot_w << 0.0, -beta, 0.0;
		wt_i.setZero();
		wtdot_i.setZero();
		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);*/

		/*m_bodies[16]->presc->setActive();
		vt_w << -beta/2.0 * t_, beta/2.0 * t_, 0.0;
		vtdot_w << -beta/2.0, beta/2.0, 0.0;
		wt_i.setZero();
		wtdot_i.setZero();
		setMaximalPrescStates(m_bodies[16], vt_w, vtdot_w, wt_i, wtdot_i);*/

		m_joints[1]->presc->setInactive();
		m_joints[5]->presc->setInactive();
		m_joints[9]->presc->setInactive();
		m_joints[13]->presc->setInactive();
		m_joints[17]->presc->setInactive();

		//m_bodies[8]->presc->setInactive();

	}
	else if (t < 30.0) {

		double t_ = t - 25.0;	// t_ = [0, 5]


		vt_w << 0.0, beta * t_, 0.0;
		vtdot_w << 0.0, beta, 0.0;
		wt_i << 0.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 0.0;

		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);


		//double t_ = t - 46.0;	// t_ = [0, 10]
		//vt_w << 0.0, -beta * t_, 0.0;
		//vtdot_w << 0.0, -beta, 0.0;
		//wt_i.setZero();
		//wtdot_i.setZero();
		//setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, wt_i, wtdot_i);

		//m_bodies[16]->presc->setActive();
		//vt_w << -beta / 2.0 * t_, beta / 2.0 * t_, 0.0;
		//vtdot_w << -beta / 2.0, beta / 2.0, 0.0;
		//wt_i.setZero();
		//wtdot_i.setZero();
		//setMaximalPrescStates(m_bodies[16], vt_w, vtdot_w, wt_i, wtdot_i);

		//setMaximalPrescStates(m_bodies[16], vt_w, vtdot_w, wt_i, wtdot_i);
	}
	else {
		m_bodies[8]->presc->setInactive();
		vt_w.setZero();
		wt_i.setZero();
		vtdot_w.setZero();
		wtdot_i.setZero();

	}
}


void World::sceneTestReducedHD(double t) {
	m_joints[0]->presc->m_q[0] = 0.0;
	m_joints[0]->presc->m_qdot[0] = 0.0;
	m_joints[0]->presc->m_qddot[0] = 0.0;

}

void World::sceneTestMaximalHD(double t) {
	Matrix4d E = m_bodies[3]->E_wi;
	Matrix3d R = E.topLeftCorner(3, 3);
	Vector6d phi = m_bodies[3]->phi;
	Vector3d vt_w, wt_i, vtdot_w, wtdot_i;

	if (t < 2.0) {
		vt_w.setZero();
		wt_i << 0.0, 0.0, t;
		vtdot_w.setZero();
		wtdot_i << 0.0, 0.0, 1.0;

	}
	else if (t < 4.0) {
		double t_ = t - 4.0;
		vt_w.setZero();
		wt_i << 0.0, 0.0, -t_;
		vtdot_w.setZero();
		wtdot_i << 0.0, 0.0, -1.0;

	}
	else if (t < 6.0) {
		double t_ = t - 4.0;
		vt_w << -2 * t_, 0.0, 0.0;
		wt_i << 0.0, 0.0, -t_;
		vtdot_w << -2.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, -1.0;
	}
	else if (t < 8.0) {
		double t_ = t - 8.0;
		vt_w << 2 * t_, 0.0, 0.0;
		wt_i << 0.0, 0.0, t_;
		vtdot_w << 2.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 1.0;
	}
	else {
		vt_w.setZero();
		wt_i.setZero();
		vtdot_w.setZero();
		wtdot_i.setZero();
	}

	Vector3d vt_i = R.transpose() * vt_w;
	m_bodies[3]->presc->m_qdot.segment<3>(0) = wt_i;
	m_bodies[3]->presc->m_qdot.segment<3>(3) = vt_i;
	m_bodies[3]->presc->m_qddot.segment<3>(0) = wtdot_i;
	m_bodies[3]->presc->m_qddot.segment<3>(3) = R.transpose() * vtdot_w - phi.segment<3>(0).cross(vt_i);
}

void World::sceneFingers(double t) {
	Vector3d vt_w, wt_i, vtdot_w, wtdot_i;
	Vector3d vt_wb, wt_ib;
	Vector3d zero = Vector3d::Zero();
	vt_w.setZero();
	vtdot_w.setZero();
	cout << t << endl;
	wt_i.setZero();
	wtdot_i.setZero();

	double beta = 0.3;
	
	if (t < 5.0) {
		vt_w << -3 * beta * t, 0.0, 0.0;
		vt_wb = m_bodies[eBone_IndexFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_IndexFinger3, vt_wb, wt_i);

		vt_w << -3 * beta * t, 0.0, 3 * beta * t;		
		vt_wb = m_bodies[eBone_MiddleFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_MiddleFinger3, vt_wb, zero);

		vt_w << -3 * beta * t, 0.0, 4 * beta * t;
		vt_wb = m_bodies[eBone_RingFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_RingFinger3, vt_wb, zero);

		vt_w << -3 *0.5 * beta * t, 0.0, 4 * beta * t;
		vt_wb = m_bodies[eBone_PinkyFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_PinkyFinger3, vt_wb, zero);

		vt_w << -3 * beta * t, 0.0, - 5 * beta * t;
		vt_wb = m_bodies[eBone_Thumb2]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_Thumb2, 0.5 * vt_wb, zero);
	}
	else if (t < 10.0) {
		double t_ = t - 5.0;
		vt_w << 3 * beta * t_, 0.0, 0.0;
		vt_wb = m_bodies[eBone_IndexFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_IndexFinger3, vt_wb, zero);

		vt_w << 1.2 * beta * t, 0.0, -1.4 * beta * t;
		vt_wb = m_bodies[eBone_MiddleFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_MiddleFinger3, vt_wb, zero);

		vt_w << 1 * beta * t, 0.0, -1* beta * t;
		vt_wb = m_bodies[eBone_RingFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_RingFinger3, vt_wb, zero);

		vt_w << 1 * 0.5 * beta * t, 0.0, - 1.5 * beta * t;
		vt_wb = m_bodies[eBone_PinkyFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_PinkyFinger3, vt_wb, zero);

		vt_w << 1 * beta * t, 0.0, 5 * beta * t;
		vt_wb = m_bodies[eBone_Thumb2]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_Thumb2, 0.2 * vt_wb, zero);

	}
	else if (t < 15.0) {
		double t_ = t - 10.0;
		double alpha = 1.7;

		vt_w << - 0.3 *alpha * beta * t, 0.0, 0.0;
		vt_wb = m_bodies[eBone_IndexFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_IndexFinger3, vt_wb, zero);

		vt_w << -0.3* alpha *beta * t, 0.0, 1 * beta * t;
		vt_wb = m_bodies[eBone_MiddleFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_MiddleFinger3, vt_wb, zero);

		vt_w << -0.2 * alpha  * beta * t, 0.0, 1 * beta * t;
		vt_wb = m_bodies[eBone_RingFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_RingFinger3, vt_wb, zero);

		vt_w << - 0.6 * alpha * beta * t, 0.0, 1.9 * beta * t;
		vt_wb = m_bodies[eBone_PinkyFinger3]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_PinkyFinger3, 0.5 * vt_wb , zero);

		vt_w << - 0.3 * alpha * beta * t, 0.0, - 0.6 * beta * t;
		vt_wb = m_bodies[eBone_Thumb2]->getBodyVelocityByEndPointVelocity(vt_w);
		setMaximalPrescStates(eBone_Thumb2, vt_wb, zero);
	}
	else if (t < 20.0) {
		// Scene 1:
		//double t_ = t - 15.0;

		//vt_w << 0.5 * beta * t, 0.0, 0.0;
		//setMaximalPrescStates(eBone_IndexFinger3, vt_w, zero);

		//vt_w << 0.3 * beta * t, 0.0, -0.7 * beta * t;
		//setMaximalPrescStates(eBone_MiddleFinger3, vt_w, zero);

		//vt_w << 0.3 *  beta * t, 0.0, -1 * beta * t;
		//setMaximalPrescStates(eBone_RingFinger3, vt_w, zero);

		//vt_w << 0.1 * 0.1 * beta * t, 0.0, -1 * beta * t;

		//setMaximalPrescStates(eBone_PinkyFinger3, 0.5 * vt_w, zero);

		//vt_w << 0.05 * beta * t, 0.0, 1 * beta * t;
		//wt_i << 0.0, 0.0, -1.0;
		//setMaximalPrescStates(eBone_Thumb2, 0.8 * vt_w, wt_i);

		// Scene 2:

		vt_wb = m_bodies[eBone_IndexFinger3]->getBodyVelocityByEndPointVelocity(zero);
		//cout << vt_wb << endl;
		wtdot_i << 0.0, 0.0, 1.0;
		setMaximalPrescStates(eBone_IndexFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_MiddleFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_MiddleFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_RingFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_RingFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_PinkyFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_PinkyFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_Thumb2]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_Thumb2, vt_wb, zero);

	}
	else if (t < 25.0) {
		vt_wb = m_bodies[eBone_IndexFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_IndexFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_MiddleFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_MiddleFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_RingFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_RingFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_PinkyFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_PinkyFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_Thumb2]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_Thumb2, vt_wb, zero);


	}
	else if (t < 30.0) {
		vt_wb = m_bodies[eBone_IndexFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_IndexFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_MiddleFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_MiddleFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_RingFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_RingFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_PinkyFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_PinkyFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_Thumb2]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_Thumb2, vt_wb, zero);

	}
	else if (t < 35.0) {
		vt_wb = m_bodies[eBone_IndexFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_IndexFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_MiddleFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_MiddleFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_RingFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_RingFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_PinkyFinger3]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_PinkyFinger3, vt_wb, zero);

		vt_wb = m_bodies[eBone_Thumb2]->getBodyVelocityByEndPointVelocity(zero);
		setMaximalPrescStates(eBone_Thumb2, vt_wb, zero);



	}
	else if (t < 40.0) {
		double t_ = t - 35.0;

		vt_w << 0.0, 3 * beta * t_, 0.0;

		//setMaximalPrescStates(eBone_IndexFinger3, vt_w, zero);



	}
	else if (t < 45.0) {



	}
	else {
		vt_w.setZero();
		wt_i.setZero();

		setMaximalPrescStates(eBone_IndexFinger3, zero, zero);
		vt_w << 1.5 * beta * t, 0.0, -2 * beta * t;

		setMaximalPrescStates(eBone_MiddleFinger3, zero, zero);

		setMaximalPrescStates(eBone_RingFinger3, zero, zero);

		setMaximalPrescStates(eBone_PinkyFinger3, zero, zero);

		setMaximalPrescStates(eBone_Thumb2, zero, zero);
		setReducedPrescStates(m_joints[eBone_Elbow], 0.0, 0.0);

	}

	double q, dq;

	// Reduced Hybrid Dynamics
	if (t < 5.0) {
		computeTargetQ(0.0, 5.0, t, D10, 0.0, q, dq);
		Vector2d qvec = Vector2d::Zero();
		Vector2d dqvec = Vector2d::Zero();
		qvec.y() = q;
		dqvec.y() = dq;

		//setReducedPrescStates(m_joints[eBone_Wrist], qvec, dqvec);
		setReducedPrescStates(m_joints[eBone_Elbow], 0.0, 0.0);
	}
	else if (t < 10.0) {
		computeTargetQ(5.0, 10.0, t, -D10, D10, q, dq);
		Vector2d qvec = Vector2d::Zero();
		Vector2d dqvec = Vector2d::Zero();
		qvec.y() = q;
		dqvec.y() = dq;
		//setReducedPrescStates(m_joints[eBone_Wrist], qvec, dqvec);
		setReducedPrescStates(m_joints[eBone_Elbow], 0.0, 0.0);

	}
	else if (t < 15.0) {
		computeTargetQ(10.0, 15.0, t, D10, 0, q, dq);
		Vector2d qvec = Vector2d::Zero();
		Vector2d dqvec = Vector2d::Zero();
		qvec.y() = q;
		dqvec.y() = dq;
		//setReducedPrescStates(m_joints[eBone_Wrist], qvec, dqvec);
		setReducedPrescStates(m_joints[eBone_Elbow], 0.0, 0.0);

	}
	else if (t < 20.0) {
		computeTargetQ(15.0, 20.0, t, -D90, D10, q, dq);
		Vector2d qvec = Vector2d::Zero();
		Vector2d dqvec = Vector2d::Zero();
		qvec.y() = q;
		dqvec.y() = dq;
		//setReducedPrescStates(m_joints[eBone_Wrist], qvec, dqvec);


		//computeTargetQ(15.0, 20.0, t, D10, 0, q, dq);
		//setReducedPrescStates(m_joints[eBone_Elbow], 0.0, 0.0);
		computeTargetQ(15.0, 20.0, t, D10*1.5, 0, q, dq);
		setReducedPrescStates(m_joints[eBone_Elbow], q, dq);

		//m_joints[eBone_Elbow]->presc->setInactive();

	}
	else if (t < 25.0) {
		computeTargetQ(20.0, 25.0, t, D90, D10 - D90, q, dq);
		Vector2d qvec = Vector2d::Zero();
		Vector2d dqvec = Vector2d::Zero();
		qvec.y() = q;
		dqvec.y() = dq;
		//setReducedPrescStates(m_joints[eBone_Wrist], qvec, dqvec);

		computeTargetQ(20.0, 25.0, t, -D10*1.5, D10*1.5, q, dq);
		setReducedPrescStates(m_joints[eBone_Elbow], q, dq);

	}
	else if (t < 30.0) {
		computeTargetQ(25.0, 30.0, t, D10*1.5, 0, q, dq);
		setReducedPrescStates(m_joints[eBone_Elbow], q, dq);

	}
	else if (t < 35.0) {		
		computeTargetQ(30.0, 35.0, t, -D10*1.5, D10*1.5, q, dq);
		setReducedPrescStates(m_joints[eBone_Elbow], q, dq);

	}
	else if (t < 40.0) {

		setReducedPrescStates(m_joints[eBone_Elbow], 0, 0);
	//	setReducedPrescStates(m_joints[eBone_Wrist], 0, 0);


	}
	else {
		//setReducedPrescStates(m_joints[eBone_Thumb1], 0.0, 0.0);
		//setReducedPrescStates(m_joints[eBone_Thumb2], 0.0, 0.0);
		//setReducedPrescStates(m_joints[eBone_Thumb3], 0.0, 0.0);

		//setReducedPrescStates(m_joints[eBone_IndexFinger3], 0.0, 0.0);
		//setReducedPrescStates(m_joints[eBone_MiddleFinger3], 0.0, 0.0);
		//setReducedPrescStates(m_joints[eBone_RingFinger3], 0.0, 0.0);
		//setReducedPrescStates(m_joints[eBone_PinkyFinger3], 0.0, 0.0);
	}
}

void World::setMaximalPrescStates(shared_ptr<Body> b, Vector3d vt_w, Vector3d vtdot_w, Vector3d wt_i, Vector3d wtdot_i) {
	b->presc->setActive();
	Matrix4d E = b->E_wi;
	Matrix3d R = E.topLeftCorner(3, 3);
	Vector6d phi = b->phi;

	Vector3d vt_i = R.transpose() * vt_w;
	b->presc->m_qdot.segment<3>(0) = wt_i;
	b->presc->m_qdot.segment<3>(3) = vt_i;
	b->presc->m_qddot.segment<3>(0) = wtdot_i;
	b->presc->m_qddot.segment<3>(3) = R.transpose() * vtdot_w - phi.segment<3>(0).cross(vt_i);
}

void World::setMaximalPrescStates(int index_body, Vector3d vt_w, Vector3d wt_i) {
	auto b = m_bodies[index_body];
	b->presc->setActive();
	Matrix4d E = b->E_wi;
	Matrix3d R = E.topLeftCorner(3, 3);
	Vector6d phi = b->phi;

	Vector3d vt_i = R.transpose() * vt_w;
	b->presc->m_qdot.segment<3>(0) = wt_i;
	b->presc->m_qdot.segment<3>(3) = vt_i;
}

void World::setReducedPrescStates(shared_ptr<Joint> j, VectorXd q, VectorXd dq) {
	j->presc->setActive();
	j->presc->m_q = q;
	j->presc->m_qdot = dq;

}

void World::setReducedPrescStates(shared_ptr<Joint> j, double q, double dq) {
	j->presc->setActive();
	j->presc->m_q[0] = q;
	j->presc->m_qdot[0] = dq;
}

void World::computeTargetQ(double t0, double t1, double t, double angle, double q0, double &q, double &dq) {

	double a = 7.0;
	double b = angle;
	double s = 2 * ((t - t0) / (t1 - t0) - 0.5);
	q = q0 + b / (1 + exp(-a * s));
	double T = t - t0;
	double TT = t0 - t1;
	double w = 2 * T;
	double P = w / TT + 1;
	double e = exp(a * P);
	double Q = T / TT + 1;
	double f = exp(a * Q);
	dq = -(2 * a*b*e) / (TT * (f + 1) *(f + 1));
}

void World::deactivateListOfPrescConstraints(Eigen::VectorXi mcon, Eigen::VectorXi rcon) {

	for (int i = 0; i < (int)mcon.size(); ++i) {
		if (mcon(i) > -1) {
			m_bodies[mcon(i)]->presc->setInactive();
		}	
	}

	for (int i = 0; i < (int)rcon.size(); ++i) {
		if (rcon(i) > -1) {
			m_joints[rcon(i)]->presc->setInactive();
		}	
	}
}

void World::setListOfReducedPrescStates(Eigen::VectorXi rcon, VectorXd q, VectorXd dq) {
	for (int i = 0; i < (int)rcon.size(); ++i) {
		setReducedPrescStates(m_joints[rcon(i)], q, dq);
	}
}


void World::setListOfMaximalPrescStates(Eigen::VectorXi mcon, Vector3d vt_w, Vector3d vtdot_w, Vector3d wt_i, Vector3d wtdot_i) {
	for (int i = 0; i < (int)mcon.size(); ++i) {
		setMaximalPrescStates(m_bodies[mcon(i)], vt_w, vtdot_w, wt_i, wtdot_i);
	}
}


// Final Version of Starfish
void World::sceneStarFish3(double t) {
	//cout << t << endl;
	double d30 = -1.0 / 6.0;
	double d45 = -1.0 / 4.0;
	double d15 = -1.0 / 12.0;
	double d60 = -1.0 / 3.0;
	double d90 = -1.0 / 2.0;
	double d10 = -1.0 / 18.0;
	double d12 = -1.0 / 15.0;
	double d18 = -1.0 / 10.0;
	double d20 = -1.0 / 9.0;
	double d22 = -1.0 / 8.0;
	Vector3d vt_w, wt_i, vtdot_w, wtdot_i;

	vt_w.setZero();
	vtdot_w.setZero();
	wt_i.setZero();
	wtdot_i.setZero();
	Vector3d Zero = Vector3d::Zero();

	double alpha = 1.0 / 7.9617 / 2.0;

	double beta = 1.4;
	double q, dq;

	if (t < 6.0) {
		// Hold one arm still and let the starfish falls down to floor and collides with floor
		setMaximalPrescStates(m_bodies[8], Zero, Zero, Zero, Zero);

		// Active:
		// mcon: 8
		// rcon: null
	}
	else if (t < 16.0) {
		// Lift an arm of the starfish to make it off the floor, velocity increases

		double t_ = t - 6.0;
		vt_w << 0.0, beta * t_, 0.0;
		vtdot_w << 0.0, beta, 0.0;

		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, Zero, Zero);

		// Active:
		// mcon: 8
		// rcon: null

	}
	else if (t < 26.0) {
		// Still lifting, velocity decreases to zero

		double t_ = t - 26.0; // t_ = [-10, 0]
		vt_w << 0.0, -beta * t_, 0.0;
		vtdot_w << 0.0, -beta, 0.0;
		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, Zero, Zero);

		if (t > 20.0 && t < 25.0) {
			double t0 = 20.0;
			double t_ = t - t0; // t_ = [0, 5]
			double t1 = 25.0;

			double alpha = 2 * M_PI / (t1 - t0);
			double sinTheta = M_PI * sin(alpha * t_);
			double cosTheta = M_PI * cos(alpha * t_);

			setReducedPrescStates(m_joints[1], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[2], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[3], d45 * sinTheta, d45 * alpha * cosTheta);
			setReducedPrescStates(m_joints[4], d45 * sinTheta, d45 * alpha * cosTheta);

		}
		else if (t > 25.0) {
			m_joints[1]->presc->setInactive();
			m_joints[2]->presc->setInactive();
			m_joints[3]->presc->setInactive();
			m_joints[4]->presc->setInactive();

		}

		// Active:
		// mcon: 8
		// rcon: null
	}
	else if (t < 56.0) {
		// Starfish starts to wiggle and then stop wiggling

		// Hold one point
		setMaximalPrescStates(m_bodies[8], Zero, Zero, Zero, Zero);

		if (t > 26.0 && t < 36.0) {
			double t0 = 26.0;
			double t_ = t - t0; // t_ = [0, 20]
			double t1 = 36.0;

			double alpha = 1.0;
			double sinTheta = M_PI * sin(alpha * t_);
			double cosTheta = M_PI * cos(alpha * t_);

			setReducedPrescStates(m_joints[1], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[2], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[3], -d30 * sinTheta, -d30 * alpha * cosTheta);
			setReducedPrescStates(m_joints[4], -d30 * sinTheta, -d30 * alpha * cosTheta);

		}
		else if (t > 36.0) {
			m_joints[1]->presc->setInactive();
			m_joints[2]->presc->setInactive();
			m_joints[3]->presc->setInactive();
			m_joints[4]->presc->setInactive();

		}

		if (t > 38.0 && t < 45.0) {
			double t0 = 38.0;
			double t_ = t - t0; // t_ = [0, 20]

			double alpha = 1.0;
			double sinTheta = M_PI * sin(alpha * t_);
			double cosTheta = M_PI * cos(alpha * t_);

			setReducedPrescStates(m_joints[9], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[10], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[11], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[12], d20 * sinTheta, d20 * alpha * cosTheta);

		}
		else if (t > 45.0) {
			m_joints[9]->presc->setInactive();
			m_joints[10]->presc->setInactive();
			m_joints[11]->presc->setInactive();
			m_joints[12]->presc->setInactive();


		}


		if (t > 50.0 && t < 55.0) {
			double t0 = 50.0;
			double t_ = t - t0; // t_ = [0, 20]
								//double sinTheta = M_PI * sin(alpha * t_);
								//double cosTheta = M_PI * cos(alpha * t_);


			double t1 = 55.0;

			double alpha = M_PI / (t1 - t0);
			double sinTheta = M_PI * sin(alpha * t_);
			double cosTheta = M_PI * cos(alpha * t_);

			setReducedPrescStates(m_joints[9], d45 * sinTheta, d45 * alpha * cosTheta);
			setReducedPrescStates(m_joints[10], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[11], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[12], -d20 * sinTheta, -d20 * alpha * cosTheta);

			setReducedPrescStates(m_joints[17], -d10 * sinTheta, -d10 * alpha * cosTheta);
			setReducedPrescStates(m_joints[18], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[19], -d10 * sinTheta, -d10 * alpha * cosTheta);
			setReducedPrescStates(m_joints[20], -d10 * sinTheta, -d10 * alpha * cosTheta);

		}
		else if (t > 55.0) {
			m_joints[9]->presc->setInactive();
			m_joints[10]->presc->setInactive();
			m_joints[11]->presc->setInactive();
			m_joints[12]->presc->setInactive();

			m_joints[17]->presc->setInactive();
			m_joints[18]->presc->setInactive();
			m_joints[19]->presc->setInactive();
			m_joints[20]->presc->setInactive();

		}


		//if (t > 50.0 && t < 55.0) {
		//	double t0 = 50.0;
		//	double t_ = t - t0; // t_ = [0, 20]
		//	//double sinTheta = M_PI * sin(alpha * t_);
		//	//double cosTheta = M_PI * cos(alpha * t_);

		//
		//	double t1 = 55.0;

		//	double alpha = 2.0;
		//	double sinTheta = M_PI * sin(alpha * t_);
		//	double cosTheta = M_PI * cos(alpha * t_);

		//	computeTargetQ(t0, t1, t, 0, 0, q, dq);
		//	setReducedPrescStates(m_joints[13], q, dq);
		//	setReducedPrescStates(m_joints[5], -d10 * sinTheta, -d10 * alpha * cosTheta);
		//	setReducedPrescStates(m_joints[9], d10 * sinTheta, d10 * alpha * cosTheta);
		//	setReducedPrescStates(m_joints[17], d20 * sinTheta, d20 * alpha * cosTheta);


		//}
		//else if (t > 55.0) {
		//	m_joints[5]->presc->setInactive();
		//	m_joints[13]->presc->setInactive();
		//	m_joints[9]->presc->setInactive();
		//	m_joints[17]->presc->setInactive();

		//}
		//setReducedPrescStates(m_joints[9], d20 * sinTheta, d20 * alpha * cosTheta);
		//setReducedPrescStates(m_joints[17], d20 * sinTheta, d20 * alpha * cosTheta);
		// Active:
		// mcon: 8
		// ncon: 5 9 13 17 
	}
	else if (t < 64.0) {
		// Settle down on the floor again

		double t_ = t - 56.0;
		vt_w << 0.0, -0.5 *beta * t_, 0.0;
		vtdot_w << 0.0, -0.5 *beta, 0.0;
		setMaximalPrescStates(m_bodies[8], vt_w, vtdot_w, Zero, Zero);

		if (t > 58.0 && t < 62.0) {
			double t0 = 58.0;
			double t_ = t - t0; // t_ = [0, 20]
			double t1 = 62.0;

			double alpha = 1.0;
			double sinTheta = M_PI * sin(alpha * t_);
			double cosTheta = M_PI * cos(alpha * t_);

			setReducedPrescStates(m_joints[13], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[14], d20 * sinTheta, d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[15], d45 * sinTheta, d45 * alpha * cosTheta);
			setReducedPrescStates(m_joints[16], d30 * sinTheta, d30 * alpha * cosTheta);

		}
		else if (t > 62.0) {
			m_joints[13]->presc->setInactive();
			m_joints[14]->presc->setInactive();
			m_joints[15]->presc->setInactive();
			m_joints[16]->presc->setInactive();
		}

		if (t > 62.0 && t < 64.0) {
			double t0 = 62.0;
			double t_ = t - t0; // t_ = [0, 20]
			double t1 = 64.0;

			double alpha = 1.0;
			double sinTheta = M_PI * sin(alpha * t_);
			double cosTheta = M_PI * cos(alpha * t_);

			setReducedPrescStates(m_joints[17], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[18], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[19], -d45 * sinTheta, -d45 * alpha * cosTheta);
			setReducedPrescStates(m_joints[20], -d30 * sinTheta, -d30 * alpha * cosTheta);

			setReducedPrescStates(m_joints[13], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[14], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[15], -d45 * sinTheta, -d45 * alpha * cosTheta);
			setReducedPrescStates(m_joints[16], -d30 * sinTheta, -d30 * alpha * cosTheta);
		}

		//VectorXi mcon(1);
		//Vector4i ncon;
		//mcon << -2;
		//ncon << 5, 9, 13, 17;
		//deactivateListOfPrescConstraints(mcon, ncon);

		//double t0 = 100.0;
		//double t_ = t - t0; // t_ = [0, 20]
		//double t1 = 120.0;
		//double p = 0.01;
	
		//wt_i << 0.0, 0.0, p * t_;
		//wtdot_i << 0.0, 0.0, p;
		//setMaximalPrescStates(m_bodies[0], Zero, Zero, wt_i, wtdot_i);
		//// wt_i  0 0 0.1
		//// Active:
		//// mcon: null
		//// ncon: null
	}
	else if (t < 67) {
		// Settle down on the floor again
		m_bodies[8]->presc->setInactive();
		m_joints[17]->presc->setInactive();
		m_joints[18]->presc->setInactive();
		m_joints[19]->presc->setInactive();
		m_joints[20]->presc->setInactive();

		m_joints[13]->presc->setInactive();
		m_joints[14]->presc->setInactive();
		m_joints[15]->presc->setInactive();
		m_joints[16]->presc->setInactive();
		// Active: 
		// mcon: 0
		// ncon: 1 2 3 4 6 7
	}
	else if (t < 77.0) {
		// Flip around

		if (t > 70.0 && t < 75.0) {
			double t0 = 70.0;
			double t_ = t - t0; // t_ = [0, 20]
			double t1 = 75.0;

			double alpha = 1.0;
			double sinTheta = M_PI * sin(alpha * t_);
			double cosTheta = M_PI * cos(alpha * t_);

			setReducedPrescStates(m_joints[17], -d10 * sinTheta, -d10 * alpha * cosTheta);
			setReducedPrescStates(m_joints[18], -d20 * sinTheta, -d20 * alpha * cosTheta);
			setReducedPrescStates(m_joints[19], -d30 * sinTheta, -d30 * alpha * cosTheta);
			setReducedPrescStates(m_joints[20], -d30 * sinTheta, -d30 * alpha * cosTheta);

		}
		else if (t > 75.0) {
			m_joints[17]->presc->setInactive();
			m_joints[18]->presc->setInactive();
			m_joints[19]->presc->setInactive();
			m_joints[20]->presc->setInactive();
		}
		//double t0 = 130.0;
		//double t_ = t - t0; // t_ = [0, 20]
		//double t1 = 190;

		//double alpha = 1.0;
		//double sinTheta = M_PI * sin(alpha * t_);
		//double cosTheta = M_PI * cos(alpha * t_);

		//computeTargetQ(t0, t1, t, 0, 0, q, dq);
		//setReducedPrescStates(m_joints[13], q, dq);
		//setReducedPrescStates(m_joints[5], -d20 * sinTheta, -d20 * alpha * cosTheta);
		//setReducedPrescStates(m_joints[9], d20 * sinTheta, d20 * alpha * cosTheta);
		//setReducedPrescStates(m_joints[17], d20 * sinTheta, d20 * alpha * cosTheta);
	
	}
	else if (t < 90.0)
	{
		
		// Active:
		// mcon: 0 
		// ncon: 4

	}
	else if (t < 95.0) {

		double t0 = 90.0;
		double t_ = t - t0;	// t_ = [0, 10]
		double t1 = 95.0;

		double alpha = 2 * M_PI / (t1 - t0);
		double sinTheta = M_PI * sin(alpha * t_);
		double cosTheta = M_PI * cos(alpha * t_);

		setReducedPrescStates(m_joints[1], d20 * sinTheta, d20 * alpha * cosTheta);
		setReducedPrescStates(m_joints[2], d20 * sinTheta, d20 * alpha * cosTheta);
		setReducedPrescStates(m_joints[3], d45 * sinTheta, d45 * alpha * cosTheta);
		setReducedPrescStates(m_joints[4], d45 * sinTheta, d45 * alpha * cosTheta);
		// Active:
		// mcon: 12 20
		// ncon: null

	}
	else if (t < 117.0) {
		m_joints[1]->presc->setInactive();
		m_joints[2]->presc->setInactive();
		m_joints[3]->presc->setInactive();
		m_joints[4]->presc->setInactive();
		// Active:
		// mcon: null
		// ncon: null
	}
	else if (t < 122.0) {
		double t0 = 117.0;
		double t_ = t - t0; // t_ = [0, 20]
		double t1 = 122.0;

		double alpha = 1.0;
		double sinTheta = M_PI * sin(alpha * t_);
		double cosTheta = M_PI * cos(alpha * t_);

		setReducedPrescStates(m_joints[5], d20 * sinTheta, d20 * alpha * cosTheta);
		setReducedPrescStates(m_joints[6], d30 * sinTheta, d30 * alpha * cosTheta);
		setReducedPrescStates(m_joints[7], d30 * sinTheta,d30 * alpha * cosTheta);
		setReducedPrescStates(m_joints[8], d30 * sinTheta,d30 * alpha * cosTheta);


	}
	else if (t < 130.0) {
		m_joints[5]->presc->setInactive();
		m_joints[6]->presc->setInactive();
		m_joints[7]->presc->setInactive();
		m_joints[8]->presc->setInactive();

		double t0 = 122.0;
		double t_ = t - t0; // t_ = [0, 20]
		double t1 = 133.0;

		double alpha = 1.0;
		double sinTheta = M_PI * sin(alpha * t_);
		double cosTheta = M_PI * cos(alpha * t_);

		setReducedPrescStates(m_joints[13], d20 * sinTheta, d20 * alpha * cosTheta);
		setReducedPrescStates(m_joints[14], d30 * sinTheta, d30 * alpha * cosTheta);
		setReducedPrescStates(m_joints[15], d30 * sinTheta, d30 * alpha * cosTheta);
		setReducedPrescStates(m_joints[16], d30 * sinTheta, d30 * alpha * cosTheta);

		setReducedPrescStates(m_joints[1], d20 * sinTheta, d20 * alpha * cosTheta);
		setReducedPrescStates(m_joints[2], d20 * sinTheta, d20 * alpha * cosTheta);
		setReducedPrescStates(m_joints[3], d20 * sinTheta, d20 * alpha * cosTheta);
		setReducedPrescStates(m_joints[4], d30 * sinTheta, d30 * alpha * cosTheta);
	}
	else {
		m_joints[13]->presc->setInactive();
		m_joints[14]->presc->setInactive();
		m_joints[15]->presc->setInactive();
		m_joints[16]->presc->setInactive();

		m_joints[1]->presc->setInactive();
		m_joints[2]->presc->setInactive();
		m_joints[3]->presc->setInactive();
		m_joints[4]->presc->setInactive();
	}
}

void World::sceneTestHyperReduced(double t) {
	double q, dq;
	if (t < 2.0) {
		
		//computeTargetQ(0.0, 2.0, t, M_PI / 8.0, 0.0, q, dq);
		//setReducedPrescStates(m_joints[0], q, dq);
		Vector3d vt_w, wt_i, vtdot_w, wtdot_i;
		vt_w.setZero();
		setMaximalPrescStates(m_bodies[1], vt_w, vt_w, vt_w, vt_w);
	}
	if (t < 4.0) {

		//computeTargetQ(2.0, 4.0, t, -M_PI / 8.0, M_PI / 8.0, q, dq);
		//setReducedPrescStates(m_joints[0], q, dq);

	}

}

// Export
Json::Value jsonworld;
Json::Value frames(Json::arrayValue);

int World::getBrenderCount() const
{
	return 1;
}

vector<string> World::getBrenderNames() const
{
	vector<string> names;
	names.push_back("World");

	return names;
}

vector<string> World::getBrenderExtensions() const
{
	vector<string> extensions;
	string json = "json";
	extensions.push_back(json);

	return extensions;
}

vector<int> World::getBrenderTypes() const
{
	vector<int> types;
	int mytype;
	mytype = Brenderable::ResetAppend;
	types.push_back(mytype);

	//types.push_back(Brenderable::Truncate);

	return types;
}

void World::exportBrender(vector< shared_ptr< ofstream > > outfiles, int frame, double time) const
{
	ofstream &outfile = *outfiles[0];

	Json::Value states(Json::arrayValue);
	vector<string> mybodyname_vec;
	string mybody = "elbow";
	mybodyname_vec.push_back(mybody);
	mybody = "wrist";
	mybodyname_vec.push_back(mybody);
	mybody = "thumb_0";
	mybodyname_vec.push_back(mybody);
	mybody = "thumb_1";
	mybodyname_vec.push_back(mybody);
	mybody = "thumb_2";
	mybodyname_vec.push_back(mybody);
	mybody = "index_finger_0";
	mybodyname_vec.push_back(mybody);
	mybody = "index_finger_1";
	mybodyname_vec.push_back(mybody);
	mybody = "index_finger_2";
	mybodyname_vec.push_back(mybody);
	mybody = "index_finger_3";
	mybodyname_vec.push_back(mybody);
	mybody = "middle_finger_0";
	mybodyname_vec.push_back(mybody);
	mybody = "middle_finger_1";
	mybodyname_vec.push_back(mybody);
	mybody = "middle_finger_2";
	mybodyname_vec.push_back(mybody);
	mybody = "middle_finger_3";
	mybodyname_vec.push_back(mybody);
	mybody = "ring_finger_0";
	mybodyname_vec.push_back(mybody);
	mybody = "ring_finger_1";
	mybodyname_vec.push_back(mybody);
	mybody = "ring_finger_2";
	mybodyname_vec.push_back(mybody);
	mybody = "ring_finger_3";
	mybodyname_vec.push_back(mybody);
	mybody = "pinky_finger_0";
	mybodyname_vec.push_back(mybody);
	mybody = "pinky_finger_1";
	mybodyname_vec.push_back(mybody);
	mybody = "pinky_finger_2";
	mybodyname_vec.push_back(mybody);
	mybody = "pinky_finger_3";
	mybodyname_vec.push_back(mybody);
	std::string resource_dir = "D:/Research/Muscles/Projects/ReducedCoordRigidBodyFEM/resources/";
	
	if (export_part == 0)
	{

		//auto elbow = addBody(density, Vector3d(len_elbow, short_side, short_side), Vector3d(len_elbow / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "36.obj");
		//auto j_elbow = addJointRevolute(elbow, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), 0.0), 0.0, RESOURCE_DIR);

		//auto wrist = addBody(density, Vector3d(short_side, short_side, short_side), Vector3d(0.0, 0.0, 0.0), preXY, RESOURCE_DIR, "wrist.obj");
		//auto j_wrist = addJointUniversalXY(wrist, Vector3d(len_elbow, 0.0, 0.0), SE3::aaToMat(Vector3d::UnitZ(), D30) * changeXY, RESOURCE_DIR, j_elbow);
		////
		//auto thumb_0 = addBody(density, Vector3d(len_thumb0, short_side, short_side), Vector3d(0.0, 0.0, len_thumb0 / 2.0), preXY, RESOURCE_DIR, "11_5.obj");
		//auto thumb_1 = addBody(density, Vector3d(len_thumb1, short_side, short_side), Vector3d(len_thumb1 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "6.obj");
		//auto thumb_2 = addBody(density, Vector3d(len_thumb2, short_side, short_side), Vector3d(len_thumb2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4.obj");

		//auto index_finger_0 = addBody(density, Vector3d(len_index_0, short_side, short_side), Vector3d(len_index_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "12.obj");
		//auto index_finger_1 = addBody(density, Vector3d(len_index_1, short_side, short_side), Vector3d(0.0, 0.0, len_index_1 / 2.0), preXY, RESOURCE_DIR, "7.obj");
		//auto index_finger_2 = addBody(density, Vector3d(len_index_2, short_side, short_side), Vector3d(len_index_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4.obj");
		//auto index_finger_3 = addBody(density, Vector3d(len_index_3, short_side, short_side), Vector3d(len_index_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3_5.obj");

		//auto middle_finger_0 = addBody(density, Vector3d(len_middle_0, short_side, short_side), Vector3d(len_middle_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "11.obj");
		//auto middle_finger_1 = addBody(density, Vector3d(len_middle_1, short_side, short_side), Vector3d(0.0, 0.0, len_middle_1 / 2.0), preXY, RESOURCE_DIR, "7.obj");
		//auto middle_finger_2 = addBody(density, Vector3d(len_middle_2, short_side, short_side), Vector3d(len_middle_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4_5.obj");
		//auto middle_finger_3 = addBody(density, Vector3d(len_middle_3, short_side, short_side), Vector3d(len_middle_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3_5.obj");

		//auto ring_finger_0 = addBody(density, Vector3d(len_ring_0, short_side, short_side), Vector3d(len_ring_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "10.obj");
		//auto ring_finger_1 = addBody(density, Vector3d(len_ring_1, short_side, short_side), Vector3d(0.0, 0.0, len_ring_1 / 2.0), preXY, RESOURCE_DIR, "6_5.obj");
		//auto ring_finger_2 = addBody(density, Vector3d(len_ring_2, short_side, short_side), Vector3d(len_ring_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "4_2.obj");
		//auto ring_finger_3 = addBody(density, Vector3d(len_ring_3, short_side, short_side), Vector3d(len_ring_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3_3.obj");

		//auto pinky_finger_0 = addBody(density, Vector3d(len_pinky_0, short_side, short_side), Vector3d(len_pinky_0 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "9.obj");
		//auto pinky_finger_1 = addBody(density, Vector3d(len_pinky_1, short_side, short_side), Vector3d(0.0, 0.0, len_pinky_1 / 2.0), preXY, RESOURCE_DIR, "4.obj");
		//auto pinky_finger_2 = addBody(density, Vector3d(len_pinky_2, short_side, short_side), Vector3d(len_pinky_2 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "3.obj");
		//auto pinky_finger_3 = addBody(density, Vector3d(len_pinky_3, short_side, short_side), Vector3d(len_pinky_3 / 2.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "2_8.obj");
		//
		Json::Value objs(Json::arrayValue);
		objs[0] = resource_dir + "36.obj";
		objs[1] = resource_dir + "wrist.obj";
		objs[2] = resource_dir + "11_5.obj";
		objs[3] = resource_dir + "6.obj";
		objs[4] = resource_dir + "4.obj";
		objs[5] = resource_dir + "12.obj";
		objs[6] = resource_dir + "7.obj";
		objs[7] = resource_dir + "4.obj";

		objs[8] = resource_dir + "3_5.obj";
		objs[9] = resource_dir + "11.obj";
		objs[10] = resource_dir + "7.obj";
		objs[11] = resource_dir + "4_5.obj";
		objs[12] = resource_dir + "3_5.obj";

		objs[13] = resource_dir + "10.obj";
		objs[14] = resource_dir + "6_5.obj";
		objs[15] = resource_dir + "4_2.obj";
		objs[16] = resource_dir + "3_3.obj";
		objs[17] = resource_dir + "9.obj";
		objs[18] = resource_dir + "4.obj";
		objs[19] = resource_dir + "3.obj";
		objs[20] = resource_dir + "2_8.obj";

		//v0["obj"] = 0;
		//v0["name"] = "elbow";
		//states.append(v0);
		//v0["obj"] = 1;
		//v0["name"] = "wrist";
		//states.append(v0);
		//v0["obj"] = 2;
		//v0["name"] = "thumb_0";
		//states.append(v0);
		//v0["obj"] = 3;
		//v0["name"] = "thumb_1";
		//states.append(v0);
		//v0["obj"] = 4;
		//v0["name"] = "thumb_2";
		//states.append(v0);
		//v0["obj"] = 5;
		//v0["name"] = "index_finger_0";
		//states.append(v0);
		//v0["obj"] = 6;
		//v0["name"] = "index_finger_1";
		//states.append(v0);
		//v0["obj"] = 7;
		//v0["name"] = "index_finger_2";
		//states.append(v0);
		//v0["obj"] = 8;
		//v0["name"] = "index_finger_3";
		//states.append(v0);
		//v0["obj"] = 9;
		//v0["name"] = "middle_finger_0";
		//states.append(v0);
		//v0["obj"] = 10;
		//v0["name"] = "middle_finger_1";
		//states.append(v0);
		//v0["obj"] = 11;
		//v0["name"] = "middle_finger_2";
		//states.append(v0);
		//v0["obj"] = 12;
		//v0["name"] = "middle_finger_3";
		//states.append(v0);
		//v0["obj"] = 13;
		//v0["name"] = "ring_finger_0";
		//states.append(v0);
		//v0["obj"] = 14;
		//v0["name"] = "ring_finger_1";
		//states.append(v0);
		//v0["obj"] = 15;
		//v0["name"] = "ring_finger_2";
		//states.append(v0);
		//v0["obj"] = 16;
		//v0["name"] = "ring_finger_3";
		//states.append(v0);
		//v0["obj"] = 17;
		//v0["name"] = "pinky_finger_0";
		//states.append(v0);
		//v0["obj"] = 18;
		//v0["name"] = "pinky_finger_1";
		//states.append(v0);
		//v0["obj"] = 19;
		//v0["name"] = "pinky_finger_2";
		//states.append(v0);
		//v0["obj"] = 20;
		//v0["name"] = "pinky_finger_3";
		//states.append(v0);

		jsonworld["header"]["objs"] = objs;
		
		for (int i = 0; i < m_bodies.size(); ++i)
		{
			Json::Value vi(Json::objectValue);
			vi["obj"] = i;
			vi["name"] = mybodyname_vec[i];
			states.append(vi);
		}
		jsonworld["header"]["states"] = states;

		Json::Value v(Json::objectValue);
		v["frame"] = frame;

		for (int i = 0; i < m_bodies.size(); ++i)
		{
			v[mybodyname_vec[i]] = m_bodies[i]->exportJson();
		}

		frames.append(v);

	}
	// print frame
	else if (export_part == 1)
	{
		Json::Value v(Json::objectValue);
		v["frame"] = frame;

		for (int i = 0; i < m_bodies.size(); ++i)
		{
			v[mybodyname_vec[i]] = m_bodies[i]->exportJson();
		}

		frames.append(v);
	}
	else {
		jsonworld["body"] = frames;
		*(outfiles[0]) << jsonworld;
	}
}

