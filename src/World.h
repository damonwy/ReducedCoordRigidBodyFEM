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

class Rigid;
class Joint;
class Body;

class World
{
public:
	World();
	virtual ~World();

	void load(const std::string &RESOURCE_DIR);

	void addBody(Rigid *body);
	void addJoint(Joint *joint);


	void init();
	void update();
	void draw();

	void setTime(double t) { m_t = t; }
	double getTime() const { return m_t; }
	void incrementTime(double dt) { m_t += dt; }

	void setGrav(Eigen::Vector3d grav) { m_grav = grav; }
	Eigen::Vector3d getGrav() const { return m_grav; }

	Rigid * getBody(int uid);
	Joint * getJoint(int uid);



private:
	Eigen::Vector3d m_grav;
	double m_t;
	int m_nlinks;

	// 
	std::vector<Rigid *> m_bodies;
	std::vector<Joint *> m_joints;



};

#endif // MUSCLEMASS_SRC_WORLD_H_