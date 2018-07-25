#include "util/parallel_bounds.h"
#include <cmath>

ParallelBounds::ParallelBounds(int nrOfThreads, int N) {
    this->nrOfThreads = nrOfThreads;
    this->N = N;
}

void ParallelBounds::setNrOfThreads(int newVal) {
    this->nrOfThreads = newVal;
}

int ParallelBounds::getNrOfThreads() {
    return nrOfThreads;
}

void ParallelBounds::setN(int newVal) {
    this->N = newVal;
}

int ParallelBounds::getN() {
    return N;
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