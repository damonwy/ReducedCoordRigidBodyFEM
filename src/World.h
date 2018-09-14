#pragma once

#ifndef MUSCLEMASS_SRC_WORLD_H_
#define MUSCLEMASS_SRC_WORLD_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>
#include "MLCommon.h"

class Joint;
class Body;
class MatrixStack;
class Program;
class Constraint;

class World
{
public:
	World();
	World(WorldType type);
	virtual ~World();

	std::shared_ptr<Body> addBody(double density, Eigen::Vector3d sides, Eigen::Vector3d p, Eigen::Matrix3d R, const std::string &RESOURCE_DIR, std::string file_name);
	std::shared_ptr<Joint> addJointRevolute(std::shared_ptr<Body> body, Eigen::Vector3d axis, Eigen::Vector3d p, Eigen::Matrix3d R, double q, std::shared_ptr<Joint> parent=nullptr);
	std::shared_ptr<Constraint> addConstraint();

	void load(const std::string &RESOURCE_DIR);
	void init();
	void update();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P);

	void setTime(double t) { m_t = t; }
	double getTime() const { return m_t; }
	double getH() const { return m_h; }
	void incrementTime(double dt) { m_t += dt; }

	void setGrav(Eigen::Vector3d grav) { m_grav = grav; }
	Eigen::Vector3d getGrav() const { return m_grav; }

	std::shared_ptr<Body> getBody(int uid);
	std::shared_ptr<Body> getBody(const std::string &name);
	std::shared_ptr<Joint> getJoint(int uid);
	std::shared_ptr<Joint> getJoint(const std::string &name);


	std::shared_ptr<Body> getBody0() const { return m_bodies[0]; }
	std::shared_ptr<Joint> getJoint0() const { return m_joints[0]; }
	std::shared_ptr<Constraint> getConstraint0() const { return m_constraints[0]; }

	Eigen::Vector2d getTspan() const { return m_tspan; }
	int getNsteps();

	int nm;
	int nr;
	int nem;
	int ner;
	int ne;
	int nim;
	int nir;


private:
	WorldType m_type;
	Eigen::Vector3d m_grav;
	double m_t;
	double m_h;
	Eigen::Vector2d m_tspan;

	int m_nbodies;
	int m_njoints;
	int m_nconstraints;

	double m_Hexpected;

	std::vector<std::shared_ptr<Body>> m_bodies;
	std::vector<std::shared_ptr<Joint>> m_joints;
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