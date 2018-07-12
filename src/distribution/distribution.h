#pragma once

#include "distribution/domain.h"
#include <stdexcept>
#include <string>

class Distribution {
public:
	virtual int createPoints(float* positions, Domain* dom) = 0;
};