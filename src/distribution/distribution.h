#pragma once

#include "distribution/domain.h"
#include <stdexcept>
#include <string>

class Distribution {
public:
	virtual int createPoints(float* positions, Domain* dom) {
		throw std::runtime_error(std::string("Calling createPoints on Distribution with abstract type."));
	}
};