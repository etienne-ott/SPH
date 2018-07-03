#pragma once

#include "distribution/distribution.h"

class Halton : public Distribution {
private:
	int N;
    int pb_1, pb_2, pb_3;

    float phi(int base, int arg);

public:
	Halton(int N, int pb_1, int pb_2, int pb_3);

	int createPoints(float* positions, Domain* dom) override;
};