##############################################################
###    General                                             ###
##############################################################

# general settings
cmake_minimum_required (VERSION 2.8)
project (MagnumFE)

# The version number.
set(MAGNUMFE_VERSION_MAJOR 0)
set(MAGNUMFE_VERSION_MINOR 1)

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
include( ${CMAKE_SOURCE_DIR}/cmake/FindGmsh.cmake )
if (NOT GMSH_FOUND)
  message (FATAL_ERROR "Cannot find libgmsh!")
endif()

# SWIG + Python
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
find_package(PythonLibs)

# Numpy
#include( ${CMAKE_SOURCE_DIR}/cmake/FindNumPy.cmake )

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
  #  ${NUMPY_INCLUDE_DIRS}
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

#string(REGEX REPLACE "-Wall" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(SWIG_MODULE_NAME cpp)

set(CMAKE_SWIG_FLAGS "")
set(INTERFACE_FILES ${PROJECT_SOURCE_DIR}/magnumfe/swig/magnumfe.i)

set_source_files_properties (${INTERFACE_FILES} PROPERTIES CPLUSPLUS ON )
set_source_files_properties (${INTERFACE_FILES} PROPERTIES SWIG_FLAGS -O)

swig_add_module(cpp python ${INTERFACE_FILES})
swig_link_libraries(cpp magnumfe)
