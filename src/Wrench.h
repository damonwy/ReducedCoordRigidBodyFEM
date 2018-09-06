#pragma once
#ifndef MUSCLEMASS_SRC_WRENCH_H_
#define MUSCLEMASS_SRC_WRENCH_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class Body;

class Wrench
{
public:
	Wrench(Body *body);
	virtual ~Wrench();

	void setWrenchLocal(Vector6d wrenchLocal);
	void setWrenchWorld(Vector6d wrenchWorld);

private:
	Body *m_body;
	Vector6d m_wrenchLocal;

};



#endif // MUSCLEMASS_SRC_WRENCH_H_