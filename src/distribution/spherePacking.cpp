#include "distribution/spherePacking.h"
#include "distribution/domain.h"
#include <cmath>

SpherePacking::SpherePacking(int N, float r) {
	this->N = N;
	this->r = r;
}

int SpherePacking::createPoints(float* positions, Domain* dom) {
	int counter = 0;

    float* x = new float[3];
	float* box = dom->getBoundingBox();

    float density = M_PI / sqrt(2.f) / 3.f,
        sq3 = sqrt(3.f),
        third = 1.f / 3.f,
        twoSq6Third = 2.f * sqrt(6.f) * third,
        r = std::pow(0.75f * dom->getVolume() * density / this->N / M_PI, third);

    x[0] = box[0];
	for (int i = 0; x[0] < box[3]; i++) {

        x[1] = box[1];
		for (int j = 0; x[1] < box[4]; j++) {

            x[2] = box[2];
            for (int k = 0; x[2] < box[5]; k++) {
                x[0] = box[0] + r * (2 * i + ((j + k) % 2));
                x[1] = box[1] + r * (sq3 * (j + (k % 2) * third));
                x[2] = box[2] + r * (k * twoSq6Third);

                if (dom->pointIsInDomain(x)) {
                    positions[counter * 3] = x[0];
                    positions[counter * 3 + 1] = x[1];
                    positions[counter * 3 + 2] = x[2];
                    counter++;

                    if (counter >= this->N) {
                        // We've reached the number of desired points before the
                        // end of the nested loops, so we need to exit
                        x[0] = box[3];
                        x[1] = box[4];
                        x[2] = box[5];
                    }
                }
			}
		}
	}

	delete[] x;

	return counter;
}