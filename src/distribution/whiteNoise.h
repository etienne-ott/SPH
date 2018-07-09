#pragma once

#include "distribution/distribution.h"

class WhiteNoise : public Distribution {
private:
	long seed;
	int N;

public:
	WhiteNoise(int N, long seed);

	int createPoints(float* positions, Domain* dom) override;
};