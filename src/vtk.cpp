#include "vtk.hpp"
#include <iostream>
#include <cmath>

using namespace std;

VTK::VTK(string path, const Kernel* kernel, int size) {
    _path = path;
    _kernel = kernel;
    _size = size;
    _count = 1;
}

void VTK::WriteDensity(double* density, double* position) {
    char* filename = new char[255];
    sprintf(filename, "%sfield_%i.vts", _path.c_str(), _count);

    FILE* handle = fopen(filename, "w");
    delete filename;

    double ds = 1.0 / _size;

    fprintf(handle, "<?xml version=\"1.0\"?>\n");
    fprintf(handle, "<VTKFile type=\"StructuredGrid\">\n");
    fprintf(handle, "<StructuredGrid WholeExtent=\"0 %i 0 %i 0 %i \">\n", _size, _size, _size);
    fprintf(handle, "<Piece Extent=\"0 %i 0 %i 0 %i \">\n", _size, _size, _size);
    fprintf(handle, "<Points>\n");
    fprintf(handle, "<DataArray type=\"Float64\" format=\"ascii\" NumberOfComponents=\"3\">\n");

    for (int z = 0; z <= _size; ++z) {
        for (int y = 0; y <= _size; ++y) {
            for (int x = 0; x <= _size; ++x) {
                fprintf(handle, "%le %le %le\n",
                    (double)x * ds,
                    (double)y * ds,
                    (double)z * ds
                );
            }
        }
    }

    fprintf(handle, "</DataArray>\n");
    fprintf(handle, "</Points>\n");
    fprintf(handle, "<PointData>\n");

    fprintf(handle,
    "<DataArray Name=\"%s\" type=\"Float64\" format=\"ascii\">\n", "density");

    for (int z = 0; z <= _size; ++z) {
        for (int y = 0; y <= _size; ++y) {
            for (int x = 0; x <= _size; ++x) {
                double d = _kernel->InterpolateDensity(x*ds, y*ds, z*ds, density, position);
                fprintf(handle, "%le ", (d > 0.5 ? 1.0 : 0.0));
            }
            fprintf(handle, "\n");
        }
        fprintf(handle, "\n");
    }

    fprintf(handle, "</DataArray>\n");
    fprintf(handle, "</PointData>\n");
    fprintf(handle, "</Piece>\n");
    fprintf(handle, "</StructuredGrid>\n");
    fprintf(handle, "</VTKFile>\n");

    fclose(handle);

    _count++;
}