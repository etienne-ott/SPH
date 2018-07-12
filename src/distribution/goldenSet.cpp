#include "distribution/goldenSet.h"
#include "distribution/domain.h"
#include "util/random_pool.h"
#include <cmath>

GoldenSet::GoldenSet(int N, long seed) {
	this->seed = seed;
	this->N = N;
}

void GoldenSet::calculateMinFibFracs(int N, int* result) {
	result[0] = 0;
	result[1] = 1;
	while (result[0] + result[1] < N) {
		result[1] = result[1] + result[0];
		result[0] = result[1] - result[0];
	}
}

int GoldenSet::goldenSequence(int M, std::vector<float>* samples, float a, float b, float s) {
	int minIndex = 0;
	float minVal = b + 1;
	int seq_counter = 0;
	float x;
	float gold_rat_conj = 0.5f * (1 + std::sqrt(5)) - 1;

	while (seq_counter < M) {
		x = s + seq_counter * gold_rat_conj * fabs(b - a);
		if (x > b) {
			int multiple = std::floor((x - a) / fabs(b - a));
			x -= multiple * fabs(b - a);
		}

		if (x < minVal) {
			minVal = x;
			minIndex = seq_counter;
		}

		samples->at(seq_counter) = x;
		seq_counter++;
	}

	return minIndex;
}

std::vector<float>* GoldenSet::sortByFibFrac(std::vector<float>* pos, int minIndex, int M, float scale) {
	int* fibFrac = new int[2];
	this->calculateMinFibFracs(M, fibFrac);

	std::vector<float>* newPos = new std::vector<float>(M);

	int sigma = minIndex + 1;
	bool once;
	newPos->at(0) = pos->at(minIndex);

	for (int i = 1; i < M; i++) {
		once = false;
		while (sigma > M || !once) {
			if (sigma <= fibFrac[0]) {
				sigma += fibFrac[1];
			} else {
				sigma -= fibFrac[0];
			}
			once = true;
		}
		newPos->at(i) = pos->at(sigma - 1) * scale;
	}

	delete[] fibFrac;
	return newPos;
}

int GoldenSet::createPoints(float* positions, Domain* dom) {
	int counter = 0;
	int minIndex;
	float* box = dom->getBoundingBox();
	float* s = new float[3];
	float* x = new float[3];
	RandomPool pool = RandomPool(this->seed);

	// init s as random point within in the domain
	do {
		s[0] = box[0] + (box[3] - box[0]) * pool.nextFloat(0.5, 1.0);
		s[1] = box[1] + (box[4] - box[1]) * pool.nextFloat(0.5, 1.0);
		s[2] = box[2] + (box[5] - box[2]) * pool.nextFloat(0.5, 1.0);
	} while(!dom->pointIsInDomain(s));

	// we need N particles, but can't determine how many of our particles
	// are within the domain until we sampled them all and for sampling
	// we need to know how many particles there are... therefore we
	// calculate as many as N multiplied by the ratio of
	// bounding box volume over domain volume. Note that due to the nature
	// of randomness we might be missing a few points towards the domain
	// boundary, but the relative number of these goes down as N increases
	float boxVol = fabs(box[3] - box[0]) * fabs(box[4] - box[1]) * fabs(box[5] - box[2]);
	int M = std::ceil(this->N * boxVol / dom->getVolume());

	std::vector<float>* samplesX = new std::vector<float>(M);
	std::vector<float>* samplesY;
	std::vector<float>* samplesZ = new std::vector<float>(M);

	// first calculate x coords as golden sequence with the correct scale
	// but no offset
	minIndex = this->goldenSequence(M, samplesX, 0.0, box[3] - box[0], s[0]);

	// The y/z coords will be a permutation of the coords scaled by the ratio
	// of bounding box side lengths. The permutation is the one that sorts
	// the previously calculated x values. There is a clever algorithm that can
	// be used in this specific case that has a runtime complexity of O(N)
	samplesY = this->sortByFibFrac(
		samplesX,
		minIndex,
		M,
		abs(box[4] - box[1]) / fabs(box[3] - box[0])
	);

	// the z coordinates will be simply another golden sequences using
	// the third seed value
    minIndex = this->goldenSequence(M, samplesZ, 0.0, box[5] - box[2], s[2]);

    // for some reason, the z values need to be shuffled
    // @TODO: Find out why
    for (int i = 0; i < M; i++) {
        int otherId = pool.nextInt(0, M - 1);
        float tmp = samplesZ->at(i);
        samplesZ->at(i) = samplesZ->at(otherId);
        samplesZ->at(otherId) = tmp;
    }

	// apply offset
	for (int i = 0; i < M; i++) {
		samplesX->at(i) += box[0];
		samplesY->at(i) += box[1];
		samplesZ->at(i) += box[2];
	}

	// now simply select the first N elements that are within the domain
	for (int i = 0; i < M && counter < this->N; i++) {
		x[0] = samplesX->at(i);
		x[1] = samplesY->at(i);
		x[2] = samplesZ->at(i);

		if (dom->pointIsInDomain(x)) {
			positions[counter * 3] = x[0];
			positions[counter * 3 + 1] = x[1];
			positions[counter * 3 + 2] = x[2];
			counter++;
		}
	}

    delete samplesX;
	delete samplesY;
	delete samplesZ;
	delete[] x;
	delete[] s;

	return counter;
}