#pragma once

class ParallelBounds {
private:
    /// @var nrOfThreads int The total number of threads
    int nrOfThreads;

    /// @var N int THe total number of particles
    int N;

public:
    /// Constructor.
    ///
    /// @param nrOfThreads int The number of threads
    /// @param N int The number of particles
    ParallelBounds(int nrOfThreads, int N);

    /// Returns the lower bound (inclusive) of the particle indices
    /// the thread with the given thread number should cover. This method
    /// is supposed to be used in for loops.
    ///
    /// @param threadNum int The thread number of the calling thread
    /// @return int The inclusive lower bound of the particle indices to cover
    int lower(int threadNum);

    /// Returns the upper bound (exclusive) of the particle indices
    /// the thread with the given thread number should cover. This method
    /// is supposed to be used in for loops.
    ///
    /// @param threadNum int The thread number of the calling thread
    /// @return int The exclusive upper bound of the particle indices to cover
    int upper(int threadNum);
};
