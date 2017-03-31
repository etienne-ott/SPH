import fnmatch
import os

# import exported environment
Import('env')

# collect sources
srcsPaths = []
srcs = []
for root, dirnames, filenames in os.walk('../src'):
    for filename in fnmatch.filter(filenames, '*.cpp'):
        srcsPaths.append(os.path.join(root, filename))

# strip the ../ from the output of the walk
for idx in srcsPaths:
    srcs.append(idx[3:])

# give the program a name
name = 'SPH'

# build it
env.Program(name, srcs)
