#include "distribution/whiteNoise.h"
#include "distribution/domain.h"
#include "util/random_pool.h"

WhiteNoise::WhiteNoise(int N, long seed) {
	this->seed = seed;
	this->N = N;
}

int WhiteNoise::createPoints(float* positions, Domain* dom) {
	RandomPool pool = RandomPool(this->seed);
	int counter = 0;

	float* box = dom->getBoundingBox();
	float* x = new float[3];

	while (counter < this->N) {
		x[0] = box[0] + (box[3] - box[0]) * pool.nextFloat(0.5f, 1.0f);
		x[1] = box[1] + (box[4] - box[1]) * pool.nextFloat(0.5f, 1.0f);
		x[2] = box[2] + (box[5] - box[2]) * pool.nextFloat(0.5f, 1.0f);

		if (dom->pointIsInDomain(x)) {
			positions[counter * 3] = x[0];
			positions[counter * 3 + 1] = x[1];
			positions[counter * 3 + 2] = x[2];
			counter++;
		}
	}

	delete[] x;

	return counter;
}