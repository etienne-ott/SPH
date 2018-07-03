cmake_minimum_required(VERSION 2.8.9)
project(SPH)

include_directories(include/yaml-cpp/include src)

# SPH project
set(SOURCES
    "src/output/ascii_output.cpp"
    "src/simulation/compute.cpp"
    "src/simulation/initialization.cpp"
    "src/main.cpp"
    "src/data/mesh.cpp"
    "src/data/neighbors.cpp"
    "src/distribution/domain.cpp"
    "src/distribution/fastPoissonDisk.cpp"
    "src/distribution/goldenSet.cpp"
    "src/distribution/halton.cpp"
    "src/distribution/hammersley.cpp"
    "src/distribution/spherePacking.cpp"
    "src/distribution/volumeGrid.cpp"
    "src/distribution/whiteNoise.cpp"
    "src/util/random_pool.cpp"
    "src/output/renderer.cpp"
    "src/output/vtk.cpp"
    "src/kernel/kernel.cpp"
    "src/kernel/gaussian.cpp"
    "src/kernel/bidomain.cpp"
    "src/kernel/cubic_spline.cpp"
)

# yaml-cpp library
#file(GLOB YAML_SOURCES_LIST "${PROJECT_SOURCE_DIR}/include/yaml-cpp/src/*.cpp")
#set(YAML_SOURCES "${YAML_SOURCES} ${YAML_SOURCES_LIST}")
set(YAML_SOURCES
    "include/yaml-cpp/src/binary.cpp"
    "include/yaml-cpp/src/collectionstack.h"
    "include/yaml-cpp/src/contrib"
    "include/yaml-cpp/src/convert.cpp"
    "include/yaml-cpp/src/directives.cpp"
    "include/yaml-cpp/src/directives.h"
    "include/yaml-cpp/src/emit.cpp"
    "include/yaml-cpp/src/emitfromevents.cpp"
    "include/yaml-cpp/src/emitter.cpp"
    "include/yaml-cpp/src/emitterstate.cpp"
    "include/yaml-cpp/src/emitterstate.h"
    "include/yaml-cpp/src/emitterutils.cpp"
    "include/yaml-cpp/src/emitterutils.h"
    "include/yaml-cpp/src/exceptions.cpp"
    "include/yaml-cpp/src/exp.cpp"
    "include/yaml-cpp/src/exp.h"
    "include/yaml-cpp/src/indentation.h"
    "include/yaml-cpp/src/memory.cpp"
    "include/yaml-cpp/src/nodebuilder.cpp"
    "include/yaml-cpp/src/nodebuilder.h"
    "include/yaml-cpp/src/node.cpp"
    "include/yaml-cpp/src/node_data.cpp"
    "include/yaml-cpp/src/nodeevents.cpp"
    "include/yaml-cpp/src/nodeevents.h"
    "include/yaml-cpp/src/null.cpp"
    "include/yaml-cpp/src/ostream_wrapper.cpp"
    "include/yaml-cpp/src/parse.cpp"
    "include/yaml-cpp/src/parser.cpp"
    "include/yaml-cpp/src/ptr_vector.h"
    "include/yaml-cpp/src/regeximpl.h"
    "include/yaml-cpp/src/regex_yaml.cpp"
    "include/yaml-cpp/src/regex_yaml.h"
    "include/yaml-cpp/src/scanner.cpp"
    "include/yaml-cpp/src/scanner.h"
    "include/yaml-cpp/src/scanscalar.cpp"
    "include/yaml-cpp/src/scanscalar.h"
    "include/yaml-cpp/src/scantag.cpp"
    "include/yaml-cpp/src/scantag.h"
    "include/yaml-cpp/src/scantoken.cpp"
    "include/yaml-cpp/src/setting.h"
    "include/yaml-cpp/src/simplekey.cpp"
    "include/yaml-cpp/src/singledocparser.cpp"
    "include/yaml-cpp/src/singledocparser.h"
    "include/yaml-cpp/src/streamcharsource.h"
    "include/yaml-cpp/src/stream.cpp"
    "include/yaml-cpp/src/stream.h"
    "include/yaml-cpp/src/stringsource.h"
    "include/yaml-cpp/src/tag.cpp"
    "include/yaml-cpp/src/tag.h"
    "include/yaml-cpp/src/token.h"
    "include/yaml-cpp/src/contrib/graphbuilder.cpp"
    "include/yaml-cpp/src/contrib/graphbuilderadapter.cpp"
    "include/yaml-cpp/src/contrib/graphbuilderadapter.h"
)

add_library(yaml ${YAML_SOURCES})

# build main project
add_executable(SPH ${SOURCES})
target_link_libraries(SPH yaml SDL2)