#include "ascii_output.hpp"
#include "parameter.hpp"
#include <string>

using namespace std;

ASCIIOutput::ASCIIOutput(string path) {
    _path = path;
    _count = 1;
}

void ASCIIOutput::WriteParticleStatus(double* density, double* position, Parameter* param) {
    char* filename = new char[255];
    sprintf(filename, "%sfield_%i.dat", _path.c_str(), _count);
    FILE* handle = fopen(filename, "w");
    delete filename;

    fprintf(handle, "x\ty\tz\tdensity\tsmoothing length\tmass\n");

    for (int i = 0; i < param->N; i++) {
        fprintf(
            handle,
            "%f\t%f\t%f\t%f\t%f\t%f\n",
            position[i * 3],
            position[i * 3 + 1],
            position[i * 3 + 2],
            density[i],
            param->h,
            param->mass
        );
    }
    fclose(handle);

    _count++;
}