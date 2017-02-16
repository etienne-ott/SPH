import glob

# import exported environment
Import('env')

# collect sources
srcs = glob.glob('../src/*.cpp')

# give the program a name
name = 'SPH'

# build it
env.Program(name, srcs)
