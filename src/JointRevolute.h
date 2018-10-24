#pragma once

#ifndef MUSCLEMASS_SRC_JOINTREVOLUTE_H_
#define MUSCLEMASS_SRC_JOINTREVOLUTE_H_

#include "Joint.h"

class SE3;
class Body;

class JointRevolute : public Joint {

public:
	JointRevolute();
	JointRevolute(std::shared_ptr<Body> body, Eigen::Vector3d axis, std::shared_ptr<Joint> parent = nullptr);
	virtual ~JointRevolute();

	virtual void updateSelf();

	//Eigen::Vector3d m_axis;

};



#endif MUSCLEMASS_SRC_JOINTREVOLUTE_H_