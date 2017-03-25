#---------------------------------------------------------------------------
# Setup the Project
#

# Requires cmake of Version 2.6, otherwise won't work
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PHASE2SIM)

#include_directories(/usr/include/GL)
#add_library(libGL STATIC IMPORTED)
#set_property(TARGET libGL PROPERTY IMPORTED_LOCATION /usr/lib/i386-linx-gnu/mesa/libGL.so)

# Find Geant4 package
#
SET(CMAKE_MODULE_PATH 
	#/home/vitya/Soft/G4/geant410-install/lib64/Geant4-10.0.1/Modules)
	#/share/data/users/bbilki/geant410/lib64/Geant4-10.0.0/Modules)
	/Users/vkhristenko/G4/geant410-install/lib64/Geant4-10.0.1/Modules)
find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIR})
find_package(Geant4 REQUIRED ui_all vis_all)

# Setup pGeant4 include directories and compiler definitions
#
include (${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Set Pythia8 Include and Lib Directories
#set(PYTHIA8 /share/data/users/vkhristenko/pythia/pythia8)
#set(PYTHIA8_LIBRARIES ${PYTHIA8}/lib/libpythia8.so)
#set(PYTHIA8_INCLUDE ${PYTHIA8}/include)
#include_directories(${PYTHIA8_INCLUDE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#include_directories(${PYTHIA8})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(PHASE2SIM Phase2Sim.cc ${sources} ${headers})
#target_link_libraries(HF_Simulation ${ROOT_LIBRARIES})
target_link_libraries(PHASE2SIM ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${PYTHIA8_LIBRARIES})

#set(wls_scripts
#  init.in
#  Sr90.mac
#  vis.mac
#  wls.in
#  wls.out
#)

#foreach(_script ${wls_scripts})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_script}
#    ${PROJECT_BINARY_DIR}/${_script}
#    COPYONLY
#)
#endforeach()


#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(wls DEPENDS wls)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS PHASE2SIM  DESTINATION bin)



