#include "distribution/halton.h"
#include "distribution/domain.h"
#include <cmath>

Halton::Halton(int N, int pb_1, int pb_2, int pb_3) {
	this->N = N;
    this->pb_1 = pb_1;
    this->pb_2 = pb_2;
	this->pb_3 = pb_3;
}

float Halton::phi(int base, int arg) {
	int runBase = base;
	int runArg = arg;
	float phi = 0.f;

	while (runArg > 0) {
		int a = runArg % base;
		phi += (float)a / runBase;
		runArg = std::floor(runArg / base);
		runBase *= base;
	}

	return phi;
}

int Halton::createPoints(float* positions, Domain* dom) {
	int counter = 0;

	float* box = dom->getBoundingBox();

	// we need N particles, but can't determine how many of our particles
	// are within the domain until we sampled them all and for sampling
	// we need to know how many particles there are... therefore we
	// calculate as many as N multiplied by the ratio of
	// bounding box volume over domain volume. Note that due to the nature
	// of randomness we might be missing a few points towards the domain
	// boundary, but the relative number of these goes down as N increases
	float boxVol = fabs(box[3] - box[0]) * fabs(box[4] - box[1]) * fabs(box[5] - box[2]);
	int M = std::ceil(this->N * boxVol / dom->getVolume());

	float* pos = new float[3 * M];

	// Create M points within the bounding box
	for (int k = 0; k < M; k++) {
		pos[k * 3] = box[0] + (box[3] - box[0]) * phi(this->pb_1, k);
		pos[k * 3 + 1] = box[1] + (box[4] - box[1]) * phi(this->pb_2, k);
		pos[k * 3 + 2] = box[2] + (box[5] - box[2]) * phi(this->pb_3, k);
	}

	// Select the first N points that are within the domain
	float* x = new float[3];
	for (int k = 0; k < M && counter < this->N; k++) {
		x[0] = pos[k * 3];
		x[1] = pos[k * 3 + 1];
		x[2] = pos[k * 3 + 2];
		if (dom->pointIsInDomain(x)) {
			positions[counter * 3] = x[0];
			positions[counter * 3 + 1] = x[1];
			positions[counter * 3 + 2] = x[2];
			counter++;
		}
	}

	delete[] pos;
	delete[] x;
	return counter;
}