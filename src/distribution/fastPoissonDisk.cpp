#include "distribution/fastPoissonDisk.h"
#include "distribution/domain.h"
#include "util/random_pool.h"
#include <vector>
#include <cmath>
#include "data/grid.h"

FastPoissonDisk::FastPoissonDisk(int N, long seed, float disk_radius, int disk_tries) {
	this->seed = seed;
	this->N = N;
    this->disk_radius = disk_radius;
    this->disk_tries = disk_tries;
}

int FastPoissonDisk::createPoints(float* positions, Domain* dom) {
	RandomPool pool = RandomPool(this->seed);
	int counter = 0;
	float* box = dom->getBoundingBox();

	// init grid
	float gridCellLength = this->disk_radius / sqrt(3.0f);
	int gridSize[3];
	gridSize[0] = ceil(fabs(box[3] - box[0]) / gridCellLength);
	gridSize[1] = ceil(fabs(box[4] - box[1]) / gridCellLength);
	gridSize[2] = ceil(fabs(box[5] - box[2]) / gridCellLength);

	Grid3D<int>* grid = new Grid3D<int>(gridSize[0], gridSize[1], gridSize[2], -1);

	std::vector<int>* activeList = new std::vector<int>();

	float sample[3];
	int gridCoords[3];
	int activeID = counter, listID;

	// choose starting point randomly from entire domain
	do {
		sample[0] = pool.nextFloat(0.5f * (box[3] + box[0]), abs(box[3] - box[0]));
		sample[1] = pool.nextFloat(0.5f * (box[4] + box[1]), abs(box[4] - box[1]));
		sample[2] = pool.nextFloat(0.5f * (box[5] + box[2]), abs(box[5] - box[2]));
	} while (!dom->pointIsInDomain(sample));

	positions[0] = sample[0];
	positions[1] = sample[1];
	positions[2] = sample[2];
	activeList->push_back(counter);
	gridCoords[0] = floor((sample[0] - box[0]) / gridCellLength);
	gridCoords[1] = floor((sample[1] - box[1]) / gridCellLength);
	gridCoords[2] = floor((sample[2] - box[2]) / gridCellLength);
	grid->set(gridCoords[0], gridCoords[1], gridCoords[2], counter);
	counter++;

	while (activeList->size() > 0 && counter < this->N) {
		// Select random new sample from active list
		listID = pool.nextInt(0, activeList->size() - 1);
		activeID = activeList->at(listID);
		activeList->erase(activeList->begin() + listID);

		for (int k = 0; k < this->disk_tries && counter < this->N; k++) {
			// choose random point in poisson shell around active sample
			// while making sure we're within the domain
			do {
				float phi = pool.nextFloat(M_PI, 2 * M_PI);
				float theta = pool.nextFloat(0.5f * M_PI, M_PI);
				float rd = pool.nextFloat(1.5f * this->disk_radius, this->disk_radius);
				sample[0] = positions[activeID * 3] + rd * sin(theta) * cos(phi);
				sample[1] = positions[activeID * 3 + 1] + rd * sin(theta) * sin(phi);
				sample[2] = positions[activeID * 3 + 2] + rd * cos(theta);
			} while (!dom->pointIsInDomain(sample));

			// check if other samples sufficiently distant
			gridCoords[0] = floor((sample[0] - box[0]) / gridCellLength);
			gridCoords[1] = floor((sample[1] - box[1]) / gridCellLength);
			gridCoords[2] = floor((sample[2] - box[2]) / gridCellLength);
			bool reject = false;

			for (int di = -1; di < 2 && !reject; di++) {
				for (int dj = -1; dj < 2 && !reject; dj++) {
					for (int dk = -1; dk < 2 && !reject; dk++) {
						int otherId = grid->at(
							std::min(gridSize[0] - 1, std::max(0, gridCoords[0] + di)),
							std::min(gridSize[1] - 1, std::max(0, gridCoords[1] + dj)),
							std::min(gridSize[2] - 1, std::max(0, gridCoords[2] + dk))
						);
						if (otherId > -1) {
							float dist = sqrt(
								(positions[otherId * 3] - sample[0]) * (positions[otherId * 3] - sample[0])
								+ (positions[otherId * 3 + 1] - sample[1]) * (positions[otherId * 3 + 1] - sample[1])
								+ (positions[otherId * 3 + 2] - sample[2]) * (positions[otherId * 3 + 2] - sample[2])
							);
							if (dist < this->disk_radius) {
								reject = true;
							}
						}
					}
				}
			}

			// found suitable candidate, add it to the active list, the samples
			// list and the background grid
			if (!reject) {
				activeList->push_back(counter);
				grid->set(gridCoords[0], gridCoords[1], gridCoords[2], counter);
				positions[counter * 3] = sample[0];
				positions[counter * 3 + 1] = sample[1];
				positions[counter * 3 + 2] = sample[2];
				counter++;
			}
		}
	}

	delete activeList;
	delete grid;

	return counter;
}