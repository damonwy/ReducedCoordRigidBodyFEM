#pragma once

#ifndef MUSCLEMASS_SRC_WRAPNULL_H_
#define MUSCLEMASS_SRC_WRAPNULL_H_
#define EIGEN_USE_MKL_ALL

#include "WrapObst.h"


class WrapNull : public WrapObst {

public:
	WrapNull() : WrapObst() {

	}
	virtual ~WrapNull() {}
};

#endif // MUSCLEMASS_SRC_WRAPNULL_H_
