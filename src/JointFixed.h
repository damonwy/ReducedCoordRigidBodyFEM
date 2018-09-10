#pragma once

#ifndef MUSCLEMASS_SRC_JOINTFIXED_H_
#define MUSCLEMASS_SRC_JOINTFIXED_H_
#include <vector>
#include <memory>
#define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Dense>
#include <iostream>
#include "Joint.h"
#include "MLCommon.h"

class SE3;
class Body;


class JointFixed : public Joint {

public:
	JointFixed();
	JointFixed(std::shared_ptr<Body> body, std::shared_ptr<Joint> parent = nullptr);
	virtual ~JointFixed();

	virtual void updateSelf();
	virtual void draw();

	Eigen::Vector3d m_axis;

};



#endif MUSCLEMASS_SRC_JOINTFIXED_H_