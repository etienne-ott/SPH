#include "kernel.hpp"
#include <string>

#ifndef __VTK_HPP
#define __VTK_HPP

using namespace std;

class VTK {
public:
    /// Constructor. Creates a new VTK instance working with the given
    /// parameters and saving the output in the given path.
    ///
    /// @param path string Where the VTK output will be saved. The path should
    ///   end with a slash if a directory is targeted, as no directory detection
    ///   takes place.
    /// @param kernel Kernel* The kernel to be used to interpolate the
    ///   density between particles
    /// @param mass double The mass of a single particle
    /// @param size int The size of the regular grid, onto which the density
    ///   values are interpolated
    VTK(string path, const Kernel* kernel, double mass, int size);

    /// Writes the density out as a VTK file, interpolating the values on a
    /// regular grid.
    ///
    /// @param density double* The density values of the particles
    /// @param position double* The position values of the particles
    /// @param N int The number of particles
    void WriteDensity(double* density, double* position, int N);

private:
    /// @var _path string The path where the files are saved
    string _path;

    /// @var _count int Counts up since construction and increments on each
    ///   call of WriteDensity, so the values of each call are saved in files
    ///   of incrementing names.
    int _count;

    /// @var _size int The size of the regular grid, onto which the density
    ///   values are interpolated.
    int _size;

    /// @var _mass double The particle mass.
    double _mass;

    /// @var _kernel Kernel The kernel used for interpolation.
    const Kernel* _kernel;
};
#endif //__VTK__HPP
