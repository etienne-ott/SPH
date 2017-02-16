# add possibility to add a debug visu
vars = Variables('custom.py')

env = Environment(variables=vars)

# do debug build?
debug = ARGUMENTS.get('debug', 0)

# set the compiler.
# For using clang in parallel you have to set all flags by hand or define a
# macro similar to mpic++

# parallel
# env.Replace(CXX='mpicxx.mpich')

# serial
env.Replace(CXX='g++')

# define some general compiler flags
env.Append(
    CXXFLAGS=[
        "-Wall",
        "-Wextra",
        "-pedantic",
        "-std=c++11",
    ],
    LIBS=[
        "SDL2",
    ]
)

# add flags for debug and release build
if debug == 0:
    env['CXXFLAGS'] += ["-O3"]
    # Enable flto for fast stuff, but slightly imprecise calculation
    # env['CXXFLAGS'] += ["-flto"]
else:
    env['CXXFLAGS'] += [
        "-g3",
        "-O0",
    ]

# call SConscript to actually build the project after setting up the environment
env.SConscript("./SConscript", exports='env', variant_dir='./build', duplicate=0)
