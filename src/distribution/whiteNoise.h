#pragma once

#include "distribution/distribution.h"

class WhiteNoise : public Distribution {
private:
	long seed;
	int N;

public:
	WhiteNoise(long seed, int N);

	int createPoints(float* positions, Domain* dom) override;
};