# Requirements
* Compiler for C++11 or greater, must be openmp 2.0+ compatible
* git
* SDL2 library installed on OS level

The following libraries will be installed and compiled together with SPH. You will need repositories to the source code of these libraries. The following sources can be used, but be aware that this might change in the future without this document being changed. In that case you will need to find the source in other ways.

* yaml-cpp: https://github.com/jbeder/yaml-cpp
* bitmap: https://github.com/ArashPartow/bitmap


# Installation
The following instructions will install and build a working copy of SPH. Please note that you should never perform instructions you found on the internet without knowing what they do. If you are unsure, please read up on the commands used before proceeding.

1. ```git clone path/to/repo```
2. ```cd SPH```
3. ```mkdir include```
3. ```mkdir build```
4. ```cd include```
5. ```git clone path/to/yaml-cpp-repo```
5. ```git clone path/to/bitmap-repo```
6. ```cd ../build```
7. ```mkdir output```
7. ```mkdir output/ascii```
7. ```mkdir output/vtk```
7. ```mkdir output/bmp```
7. ```cmake -DCMAKE_BUILD_TYPE=Release -DPARALLEL_BUILD=True ..```
8. ```make -j```
9. ```cp ../default_parameter.yaml default_parameter.yaml```

# Using SPH
After installation you can execute the SPH binary in the build subdirectory of the project directory. Note that filepaths used in the program are written relative to the project directory and files might be saved accordingly.

## Build types
There two build types supported: "Release" and "Debug". They differ in the compiler flags used, which also effects the speed of execution or the presence of debugging structures in the compiled source. By default a debug build is made, but this behavious can be changed by using the cmake flag ```-DCMAKE_BUILD_TYPE=Type```, where ```Type``` can be replaced with the desired build type.

## Parallel Build
During calling cmake for a build we can specify if the build is a serial build or a parallel build by using the cmake flag ```-DPARALLEL_BUILD=Flag```, where ```Flag``` can be replaced with the strings ```True``` or ```False```. Please note that in either case due to technical limitations the compiler used must be openmp compatible. However the serial build can still be used in an architecture where multithreading is not possible.

In the parameter file we can specify the number of threads to be used with the parameter ```nr_of_threads```. This need not be the number of physical processors, but usually is chosen as such.

If you switch the cmake flag ```-DPARALLEL_BUILD=Flag``` you must afterwards run ```make clean``` and recompile the project with ```make -j``` as make will not recognize the flag as a reason to recompile.

## Parameters
The parameter file "default_parameter.yaml" contains all parameters that are intended to be changed without recompiling the project. You can find short descriptions within the file and more detailed ones in this document.

### Details on the parameters
@TODO Coming soon(ish)

## Creating a video of the simulation
This software does not provide a good renderer to create nice image from the simulation data. Nevertheless you can quickly make video of a very simple rendering for debug and testing purposes.

First you need to run a simulation and enable the BMP output option. During the simulation the output of the debug renderering will be written as BMP files in the output folder.

Once the simulation is done we can now create a video from the individual frames. We will need the ffmpeg software and libraries installed to do this as the SPH program does not provide any video functionality.

A sample command to make the video would be:

```ffmpeg -framerate 24 -i %d.bmp output.mp4```