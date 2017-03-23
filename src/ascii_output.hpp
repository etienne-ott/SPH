#include "parameter.hpp"
#include <string>

#ifndef __ASCII_OUTPUT
#define __ASCII_OUTPUT

using namespace std;

class ASCIIOutput {
public:

    /// Constructor. Creates a new instance working in the given path.
    ///
    /// @param path string The path where the output files will be stored.
    ASCIIOutput(string path);

    /// Writes the given particle status out to a CSV file, that is ASCII
    /// encoded. Each call since instantiation will increment the file
    /// names by one, so the files will be called field_1.dat, field_2.dat etc.
    ///
    /// @param density double* The particle densities
    /// @param position double* The particle positions
    /// @param param Parameter* The parameter object holding the simulation
    ///   parameters.
    void WriteParticleStatus(double* density, double* position, Parameter* param);

private:
    /// @var _path string The path where the data files will be stored.
    string _path;

    /// @var _count int A counter for the number of files since instantiation.
    int _count;
};
#endif
