cmake_minimum_required(VERSION 3.20)
project(scint_sim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# g4data — exports atomic and photon-attenuation data from Geant4's own NIST/EPICS
# databases, so that the Python screening layer and the transport simulation share
# a single source of truth without vendoring licence-restricted data files.
add_executable(g4data tools/g4data.cc)
target_link_libraries(g4data PRIVATE ${Geant4_LIBRARIES})
