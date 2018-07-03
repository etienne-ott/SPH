#pragma once

#include "distribution/distribution.h"

class VolumeGrid : public Distribution {
private:
	int N;

public:
	VolumeGrid(int N);

	int createPoints(float* positions, Domain* dom) override;
};