#pragma once

#include "distribution/distribution.h"

class FastPoissonDisk : public Distribution {
private:
	long seed;
	int N;
    float disk_radius;
    int disk_tries;

public:
	FastPoissonDisk(int N, long seed, float disk_radius, int disk_tries);

	int createPoints(float* positions, Domain* dom) override;
};