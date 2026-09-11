cmake_minimum_required(VERSION 3.20)
project(scint_sim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Shared simulation components: geometry, surfaces and material handling used by
# both the optical application and the data-export tool.
add_library(scintsim STATIC
    src/SurfaceCatalogue.cc
    src/MaterialLibrary.cc
    src/DetectorConstruction.cc
    src/PhysicsList.cc
    src/Actions.cc
)
target_include_directories(scintsim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(scintsim PUBLIC ${Geant4_LIBRARIES})

# g4data — exports atomic data, photon attenuation and the surface catalogue from
# Geant4's own databases, so the Python layer and the simulation share one source
# of truth without vendoring licence-restricted data files.
add_executable(g4data tools/g4data.cc)
target_link_libraries(g4data PRIVATE scintsim)

# scint_optical -- the light-collection application driven entirely by a macro.
add_executable(scint_optical scint_optical.cc)
target_link_libraries(scint_optical PRIVATE scintsim)
