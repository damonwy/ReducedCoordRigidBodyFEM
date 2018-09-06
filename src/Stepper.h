#pragma once
#ifndef MUSCLEMASS_SRC_STEPPER_H_
#define MUSCLEMASS_SRC_STEPPER_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>
#include "MLCommon.h"

class World;

class Stepper
{
public:

	Stepper();
	virtual ~Stepper();

	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step(World *world);

protected:
	double m_dt;


};



#endif // MUSCLEMASS_SRC_STEPPER_H_

