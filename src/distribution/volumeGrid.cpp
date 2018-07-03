#include "distribution/volumeGrid.h"
#include "distribution/domain.h"
#include <cmath>

VolumeGrid::VolumeGrid(int N) {
	this->N = N;
}

int VolumeGrid::createPoints(float* positions, Domain* dom) {
	int counter = 0;

	float deltaLength = std::pow(dom->getVolume() / this->N, 1.0f / 3.0f);

    float* x = new float[3];
	float* box = dom->getBoundingBox();

	for (int i = 0; i < 3; i++) {
		x[i] = box[i];
	}

	for (int i = 0; x[0] < box[3]; i++) {
		x[1] = box[1];

		for (int j = 0; x[1] < box[4]; j++) {
            x[2] = box[2];

            for (int k = 0; x[2] < box[5]; k++) {

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

                x[2] += deltaLength;
			}

			x[1] += deltaLength;
		}

		x[0] += deltaLength;
	}

	delete[] x;

	return counter;
}