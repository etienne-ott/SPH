cmake_minimum_required(VERSION 2.8.9)
project(SPH)

include_directories(include/yaml-cpp/include src)

# SPH project
set(SOURCES
    "src/ascii_output.cpp"
    "src/compute.cpp"
    "src/initialization.cpp"
    "src/main.cpp"
    "src/random_pool.cpp"
    "src/renderer.cpp"
    "src/vtk.cpp"
    "src/kernel/kernel.cpp"
    "src/kernel/gaussian.cpp"
    "src/kernel/bidomain.cpp"
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