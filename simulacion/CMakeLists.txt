cmake_minimum_required(VERSION 3.16)
project(ComptonSpectrum)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(compton_spectrum compton_spectrum.cc ${sources})
target_link_libraries(compton_spectrum ${Geant4_LIBRARIES})

set(SCRIPTS run.mac vis.mac co60.mac cs137.mac na22.mac run_vis.sh)
foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY
    )
endforeach()
