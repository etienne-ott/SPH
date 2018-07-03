#pragma once

#include "distribution/distribution.h"

class Hammersley : public Distribution {
private:
	int N;
    int pb_1, pb_2;

    float phi(int base, int arg);

public:
	Hammersley(int N, int pb_1, int pb_2);

	int createPoints(float* positions, Domain* dom) override;
};