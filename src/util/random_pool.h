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
    float nextFloat();

    /// Returns a randomly selected float value between
    /// (mean - 0.5 * scale) and (mean + 0.5 * scale).
    ///
    /// @param mean float The mean value of the distribution
    /// @param scale float The spread of the distribution
    /// @return float A random float in a range determined by the mean and
    ///     scale parameters.
    float nextFloat(float mean, float scale);

    /// Returns a randomly selected int value between and including the
    /// min and max values.
    ///
    /// @param min int The minimum value
    /// @param max int The maximum value
    /// @return int A random int value between and including the given limits
    int nextInt(int min, int max);

    /// Returns a randomly selected double value between 0.0 and 1.0
    ///
    /// @return double A random double between 0.0 and 1.0
    double nextDouble();

    /// Returns a randomly selected double value between
    /// (mean - 0.5 * scale) and (mean + 0.5 * scale).
    ///
    /// @param mean double The mean value of the distribution
    /// @param scale double The spread of the distribution
    /// @return double A random double in a range determined by the mean and
    ///     scale parameters.
    double nextDouble(double mean, double scale);

private:
    /// @var _seed long The seed used by the random number generator.
    long _seed;

    /// @var _generator std::default_random_engine The random engine used.
    std::default_random_engine _generator;
};
