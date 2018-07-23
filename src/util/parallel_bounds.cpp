#include "util/parallel_bounds.h"
#include <cmath>

ParallelBounds::ParallelBounds(int nrOfThreads, int N) {
    this->nrOfThreads = nrOfThreads;
    this->N = N;
}

int ParallelBounds::lower(int threadNum) {
    int partial = std::floor(this->N / this->nrOfThreads);
    return threadNum * partial;
}

int ParallelBounds::upper(int threadNum) {
    int partial = std::floor(this->N / this->nrOfThreads);
    if (threadNum == this->nrOfThreads - 1) {
        return this->N;
    } else {
        return (threadNum + 1) * partial;
    }
}