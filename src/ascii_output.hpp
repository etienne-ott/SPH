#include <string>
#include <yaml-cpp/yaml.h>

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
    /// @param density float* The particle densities
    /// @param position float* The particle positions
    /// @param pressure float* The particle pressures
    /// @param param YAML::Node& The parameter object holding the simulation
    ///   parameters.
    void WriteParticleStatus(float* density, float* position, float* pressure, YAML::Node& param);

private:
    /// @var _path string The path where the data files will be stored.
    string _path;

    /// @var _count int A counter for the number of files since instantiation.
    int _count;
};
#endif
