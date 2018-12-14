#pragma once

#ifndef MUSCLEMASS_SRC_WORLD_H_
#define MUSCLEMASS_SRC_WORLD_H_
#define EIGEN_USE_MKL_ALL

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>
#include "MLCommon.h"

class Joint;
class JointRevolute;
class Body;
class SoftBody;
class SoftBodyInvertibleFEM;
class SoftBodyCorotationalLinear;
class MatrixStack;
class Program;
class Constraint;
class ConstraintJointLimit;
class ConstraintNull;
class Spring;
class SpringNull;
class SpringDamper;
class Deformable;
class DeformableSpring;
class DeformableNull;
class JointNull;
class JointFixed;
class Comp;
class CompNull;
class CompSphere;
class CompCylinder;
class CompDoubleCylinder;
class WrapObst;
class WrapSphere;
class WrapCylinder;
class WrapDoubleCylinder;
class MeshEmbedding;
class MeshEmbeddingNull;

enum WorldType {
	SERIAL_CHAIN,
	DIFF_REVOLUTE_AXES,
	BRANCHING,
	SHPERICAL_JOINT,
	LOOP,
	JOINT_TORQUE,
	JOINT_LIMITS,
	EQUALITY_CONSTRAINED_ANGLES,
	EQUALITY_AND_LOOP, HYBRID_DYNAMICS,
	EXTERNAL_WORLD_FORCE,
	JOINT_STIFFNESS,
	SPRINGS,
	SOFT_BODIES,
	COMPONENT,
	WRAP_SPHERE,
	WRAP_CYLINDER,
	WRAP_DOUBLECYLINDER,
	SPLINE_CURVE_JOINT,
	SPLINE_SURFACE_JOINT,
	SOFT_BODIES_CUBE_INVERTIBLE,
	SOFT_BODIES_CYLINDER_INVERTIBLE,
	SOFT_BODIES_CUBE_COROTATIONAL_LINEAR,
	SOFT_BODIES_CYLINDER_COROTATIONAL_LINEAR,
	SPRING_DAMPER,
	MESH_EMBEDDING,
	HUMAN_BODY
};


class World
{
public:
	World();
	World(WorldType type);
	virtual ~World() {}

	std::shared_ptr<Body> addBody(
		double density, 
		Vector3d sides, 
		Vector3d p, 
		Matrix3d R, 
		const std::string &RESOURCE_DIR, 
		std::string file_name);

	std::shared_ptr<JointRevolute> addJointRevolute(
		std::shared_ptr<Body> body, 
		Vector3d axis, 
		Vector3d p, 
		Matrix3d R, 
		double q, 
		const std::string &RESOURCE_DIR,
		std::shared_ptr<Joint> parent=nullptr);
	
	std::shared_ptr<JointFixed> addJointFixed(
		std::shared_ptr<Body> body, 
		Vector3d p, 
		Matrix3d R, 
		double q, 
		std::shared_ptr<Joint> parent = nullptr);
	
	std::shared_ptr<ConstraintJointLimit> addConstraintJointLimit(
		std::shared_ptr<Joint> joint, 
		double ql, 
		double qu);
	
	std::shared_ptr<DeformableSpring> addDeformableSpring(
		double mass, 
		int n_points, 
		std::shared_ptr<Body> body0, 
		Vector3d r0, 
		std::shared_ptr<Body> body1, 
		Vector3d r1);
	
	std::shared_ptr<SoftBody> addSoftBody(
		double density, 
		double young, 
		double possion, 
		Material material, 
		const std::string &RESOURCE_DIR, 
		std::string file_name);

	std::shared_ptr<SoftBodyInvertibleFEM> addSoftBodyInvertibleFEM(
		double density,
		double young,
		double possion,
		Material material,
		const std::string &RESOURCE_DIR,
		std::string file_name);

	std::shared_ptr<SoftBodyCorotationalLinear> addSoftBodyCorotationalLinearFEM(
		double density,
		double young,
		double possion,
		Material material,
		const std::string &RESOURCE_DIR,
		std::string file_name);

	std::shared_ptr<MeshEmbedding> addMeshEmbedding(
		double density,
		double young,
		double possion,
		Material material, 
		const std::string &RESOURCE_DIR,
		std::string dense_mesh,
		std::string coarse_mesh
	);

	std::shared_ptr<ConstraintNull> addConstraintNull();
	std::shared_ptr<DeformableNull> addDeformableNull();
	std::shared_ptr<JointNull> addJointNull();
	std::shared_ptr<CompNull> addCompNull();
	std::shared_ptr<WrapObst> addWrapNull();
	std::shared_ptr<SpringNull> addSpringNull();
	std::shared_ptr<MeshEmbeddingNull> addMeshEmbeddingNull();

	std::shared_ptr<CompSphere> addCompSphere(
		double r, 
		std::shared_ptr<Body> parent, 
		Matrix4d E, 
		const std::string &RESOURCE_DIR);
	
	std::shared_ptr<CompCylinder> addCompCylinder(
		double r, 
		std::shared_ptr<Body> parent, 
		Matrix4d E, 
		Vector3d z, 
		Vector3d o, 
		const std::string &RESOURCE_DIR, 
		std::string shape);
	
	std::shared_ptr<CompDoubleCylinder> addCompDoubleCylinder(
		double rA, 
		std::shared_ptr<Body> parentA, 
		Matrix4d EA, 
		Vector3d z_a, 
		Vector3d o_a,
		double rB, 
		std::shared_ptr<Body> parentB, 
		Matrix4d EB,
		Vector3d z_b, 
		Vector3d o_b,
		const std::string &RESOURCE_DIR, 
		std::string shapeA, 
		std::string shapeB);
	
	std::shared_ptr<WrapSphere> addWrapSphere(
		std::shared_ptr<Body> b0, 
		Vector3d r0, 
		std::shared_ptr<Body> b1, 
		Vector3d r1, 
		std::shared_ptr<CompSphere> compSphere, 
		int num_points, 
		const std::string &RESOURCE_DIR);

	std::shared_ptr<WrapCylinder> addWrapCylinder(
		std::shared_ptr<Body> b0, 
		Vector3d r0, 
		std::shared_ptr<Body> b1, 
		Vector3d r1, 
		std::shared_ptr<CompCylinder> compCylinder, 
		int num_points, 
		const std::string &RESOURCE_DIR);

	std::shared_ptr<WrapDoubleCylinder> addWrapDoubleCylinder(
		std::shared_ptr<Body> b0,
		Vector3d r0,
		std::shared_ptr<Body> b1,
		Vector3d r1,
		Vector3d u,
		Vector3d v,
		Vector3d z_u,
		Vector3d z_v,
		std::shared_ptr<CompDoubleCylinder> compDoubleCylinder,
		int num_points,
		const std::string &RESOURCE_DIR);


	Energy computeEnergy();

	void load(const std::string &RESOURCE_DIR);
	void init();
	void update();
	
	void draw(
		std::shared_ptr<MatrixStack> MV, 
		const std::shared_ptr<Program> prog, 
		const std::shared_ptr<Program> progSimple, 
		const std::shared_ptr<Program> progSoft, 
		std::shared_ptr<MatrixStack> P);

	void setTime(double t) { m_t = t; }
	double getTime() const { return m_t; }
	double getH() const { return m_h; }
	void incrementTime() { m_t += m_h; }

	void setGrav(Eigen::Vector3d grav) { m_grav = grav; }
	Vector3d getGrav() const { return m_grav; }

	std::shared_ptr<Body> getBody(int uid);
	std::shared_ptr<Body> getBody(const std::string &name);
	std::shared_ptr<Joint> getJoint(int uid);
	std::shared_ptr<Joint> getJoint(const std::string &name);


	std::shared_ptr<Body> getBody0() const { return m_bodies[0]; }
	std::shared_ptr<Joint> getJoint0() const { return m_joints[0]; }
	std::shared_ptr<Deformable> getDeformable0() const { return m_deformables[0]; }
	std::shared_ptr<SoftBody> getSoftBody0() const { return m_softbodies[0]; }
	std::shared_ptr<Constraint> getConstraint0() const { return m_constraints[0]; }
	std::shared_ptr<Spring> getSpring0() const { return m_springs[0]; }
	std::shared_ptr<MeshEmbedding> getMeshEmbedding0() const { return m_meshembeddings[0]; }

	Vector2d getTspan() const { return m_tspan; }
	int getNsteps();

	int nm;
	int nr;
	int nem;
	int ner;
	int ne;
	int nim;
	int nir;
	int m_countS;
	int m_countCM;

	int m_nbodies;
	int m_nsoftbodies;
	int m_njoints;
	int m_ndeformables;
	int m_nsprings;
	int m_nconstraints;
	int m_ncomps;
	int m_nwraps;
	int m_nmeshembeddings;
	
	int m_dense_nm;
	int m_dense_nr;
	int m_ntets;
private:
	Energy m_energy;		// the energy in current state
	Energy m_energy0;		// the energy in initial state

	WorldType m_type;
	Eigen::Vector3d m_grav;
	double m_stiffness;
	double m_damping;
	double m_t;
	double m_h;
	Eigen::Vector2d m_tspan;	
	bool isleftleg;
	bool isrightleg;

	double m_Hexpected;		// used to check correctness

	std::vector<std::shared_ptr<Body>> m_bodies;
	std::vector<std::shared_ptr<Comp>> m_comps;
	std::vector<std::shared_ptr<WrapObst>> m_wraps;
	std::vector<std::shared_ptr<SoftBody>> m_softbodies;
	std::vector<std::shared_ptr<MeshEmbedding>> m_meshembeddings;
	std::vector<std::shared_ptr<Joint>> m_joints;
	std::vector<std::shared_ptr<Deformable>> m_deformables;
	std::vector<std::shared_ptr<Spring>> m_springs;
	std::vector<std::shared_ptr<Constraint>> m_constraints;

	typedef std::map<std::string, std::shared_ptr<Body>> MapBodyName;
	typedef std::map<int, std::shared_ptr<Body>> MapBodyUID;

	typedef std::map<std::string, std::shared_ptr<Joint>> MapJointName;
	typedef std::map<int, std::shared_ptr<Joint>> MapJointUID;

	MapBodyName m_bodyName;
	MapBodyUID m_bodyUID;
	MapJointName m_jointName;
	MapJointUID m_jointUID;
};

#endif // MUSCLEMASS_SRC_WORLD_H_