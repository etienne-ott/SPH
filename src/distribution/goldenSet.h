#pragma once

#include "distribution/distribution.h"
#include <vector>

class GoldenSet : public Distribution {
private:
	long seed;
	int N;

    void calculateMinFibFracs(int N, int* result);

    int goldenSequence(int M, std::vector<float>* samples, float a, float b, float s);

    std::vector<float>* sortByFibFrac(std::vector<float>* pos, int minIndex, int M, float scale);

public:
	GoldenSet(int N, long seed);

	int createPoints(float* positions, Domain* dom) override;
};