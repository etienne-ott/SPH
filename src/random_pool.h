#pragma once

#include <random>

class RandomPool {
public:
    /// Constructor. Initializes the random pool with the current unix
    /// time as seed.
    RandomPool();

    /// Constructor. Initializes the random pool with the given seed.
    RandomPool(long seed);

    /// Returns a randomly selected float value between 0.0 and 1.0
    ///
    /// @return float A random float between 0.0 and 1.0
    float NextFloat();

    /// Returns a randomly selected float value between
    /// (mean - 0.5 * scale) and (mean + 0.5 * scale).
    ///
    /// @param mean float The mean value of the distribution
    /// @param scale float The spread of the distribution
    /// @return float A random float in a range determined by the mean and
    ///     scale parameters.
    float NextFloat(float mean, float scale);

private:
    /// @var _seed long The seed used by the random number generator.
    long _seed;

    /// @var _generator std::default_random_engine The random engine used.
    std::default_random_engine _generator;
};
