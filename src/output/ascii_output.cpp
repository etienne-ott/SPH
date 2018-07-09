#include "output/ascii_output.h"
#include <string>

using namespace std;

ASCIIOutput::ASCIIOutput(string path) {
    _path = path;
    _count = 1;
}

void ASCIIOutput::WriteParticleStatus(float* density, float* position, float* pressure, YAML::Node& param) {
    char* filename = new char[255];
    sprintf(filename, "%sfield_%i.dat", _path.c_str(), _count);
    FILE* handle = fopen(filename, "w");
    delete filename;

    fprintf(handle, "x\ty\tz\tdensity\tpressure\n");

    for (int i = 0; i < param["N"].as<int>(); i++) {
        fprintf(
            handle,
            "%f\t%f\t%f\t%f\t%f\n",
            position[i * 3],
            position[i * 3 + 1],
            position[i * 3 + 2],
            density[i],
            pressure[i]
        );
    }
    fclose(handle);

    _count++;
}