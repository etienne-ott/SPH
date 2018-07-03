#pragma once

#include "distribution/distribution.h"

class SpherePacking : public Distribution {
private:
	int N;
	float r;

public:
	SpherePacking(int N, float r);

	int createPoints(float* positions, Domain* dom) override;
};