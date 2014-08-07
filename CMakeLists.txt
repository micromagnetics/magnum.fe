##############################################################
###    General                                             ###
##############################################################

# general settings
cmake_minimum_required (VERSION 2.8)
project (MagnumFE)

# The version number.
set(MAGNUMFE_VERSION_MAJOR 0)
set(MAGNUMFE_VERSION_MINOR 1)

# set cmake module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Add configuration
#configure_file(
#  "${PROJECT_BINARY_DIR}/config.h"
#)

##############################################################
###    Macros                                              ###
##############################################################

# Macro for adding src files from within subdirectories
macro(append_magnumfe_src DIR FILES)
  set(LIST)
  foreach(FILE ${FILES})
    list(APPEND LIST ${DIR}/${FILE})
  endforeach(FILE)
  set(MAGNUMFE_SRC ${MAGNUMFE_SRC};${LIST} PARENT_SCOPE)
endmacro(append_magnumfe_src)

##############################################################
###    Libraries   TODO: check versions                    ###
##############################################################

# Dolfin
find_package(DOLFIN)

# GMSH
find_package(Gmsh)
if (NOT GMSH_FOUND)
  message (FATAL_ERROR "Cannot find libgmsh!")
endif()

# SWIG + Python
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
find_package(PythonLibs)
find_package(Python 2 REQUIRED)

# Numpy
find_package(NumPy)

##############################################################
###    Subdirectories                                      ###
##############################################################

add_subdirectory(magnumfe)
add_subdirectory(magnumfe/expressions)
add_subdirectory(magnumfe/swig)

##############################################################
###    Set Compiler Flags                                  ###
##############################################################

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Compiler flags
set(CMAKE_CXX_FLAGS "-O3 -Wall ${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

##############################################################
###    Compile Library                                     ###
##############################################################

include_directories(
  ${PROJECT_SOURCE_DIR}/magnumfe
  ${PROJECT_BINARY_DIR}
  ${PYTHON_INCLUDE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${DOLFIN_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIRS}
)
include_directories(SYSTEM
  ${DOLFIN_3RD_PARTY_INCLUDE_DIRS}
)

link_directories(${DOLFIN_LIBRARY_DIRS})

add_library(magnumfe SHARED ${MAGNUMFE_SRC})

target_link_libraries(magnumfe
  ${DOLFIN_LIBRARIES}
  ${DOLFIN_3RD_PARTY_LIBRARIES}
  ${GMSH_LIBRARY}
)

##############################################################
###    Swig Stuff                                          ###
##############################################################

set(SWIG_MODULE_NAME cpp)

# Disable warnings for SWIG wrapper code
string(REGEX REPLACE "-Wall" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Set NumPy definitions
set(MAGNUMFE_PYTHON_DEFINITIONS -DNUMPY_VERSION_MAJOR=${NUMPY_VERSION_MAJOR} -DNUMPY_VERSION_MINOR=${NUMPY_VERSION_MINOR} -DNUMPY_VERSION_MICRO=${NUMPY_VERSION_MICRO})
if(NUMPY_VERSION VERSION_GREATER 1.6.2)
  set(MAGNUMFE_PYTHON_DEFINITIONS ${MAGNUMFE_PYTHON_DEFINITIONS} -DNPY_NO_DEPRECATED_API=NPY_${NUMPY_VERSION_MAJOR}_${NUMPY_VERSION_MINOR}_API_VERSION)
endif()

add_definitions(${MAGNUMFE_PYTHON_DEFINITIONS})
set(CMAKE_SWIG_FLAGS ${MAGNUMFE_PYTHON_DEFINITIONS})

set(INTERFACE_FILES ${PROJECT_SOURCE_DIR}/magnumfe/swig/magnumfe.i)

set_source_files_properties (${INTERFACE_FILES} PROPERTIES CPLUSPLUS ON )
set_source_files_properties (${INTERFACE_FILES} PROPERTIES SWIG_FLAGS -O)

swig_add_module(cpp python ${INTERFACE_FILES})
swig_link_libraries(cpp magnumfe)


##############################################################
###    Install Rule                                        ###
##############################################################

message ("-- Found Python ${PYTHON_VERSION}, will install at ${PYTHON_SITE_MODULES}")

# install directory magnum, exclude the symlinks
install(DIRECTORY site-packages/magnumfe DESTINATION "${PYTHON_SITE_MODULES}" PATTERN "site-packages/magnumfe/cpp.py site-packages/magnumfe/_cpp.so" EXCLUDE)
install(FILES ${CMAKE_BINARY_DIR}/cpp.py ${CMAKE_BINARY_DIR}/_cpp.so DESTINATION "${PYTHON_SITE_MODULES}/magnumfe")
