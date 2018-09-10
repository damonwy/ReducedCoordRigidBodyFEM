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

class World
{
public:
	World();
	World(WorldType type);
	virtual ~World();


	void addBody(std::shared_ptr<Body> body);
	void addJoint(std::shared_ptr<Joint> joint);

	void load(const std::string &RESOURCE_DIR);
	void init();
	void update();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P);

	void setTime(double t) { m_t = t; }
	double getTime() const { return m_t; }
	void incrementTime(double dt) { m_t += dt; }

	void setGrav(Eigen::Vector3d grav) { m_grav = grav; }
	Eigen::Vector3d getGrav() const { return m_grav; }

	std::shared_ptr<Body> getBody(int uid);
	std::shared_ptr<Body> getBody(const std::string &name);
	std::shared_ptr<Joint> getJoint(int uid);
	std::shared_ptr<Joint> getJoint(const std::string &name);

	void updateQ();
	void updateQDot();

	int nm;
	int nr;
private:
	WorldType m_type;
	Eigen::Vector3d m_grav;
	double m_t;
	int m_nbodies;
	int m_njoints;
	double m_Hexpected;



	// 
	std::vector<std::shared_ptr<Body>> m_bodies;
	std::vector<std::shared_ptr<Joint>> m_joints;

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