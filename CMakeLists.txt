##############################################################
###    General                                             ###
##############################################################

# general settings
cmake_minimum_required (VERSION 2.6)
project (MagnumFE)

# The version number.
set(MAGNUMFE_VERSION_MAJOR 0)
set(MAGNUMFE_VERSION_MINOR 1)

# Add configuration
#configure_file(
#  "${PROJECT_SOURCE_DIR}/magneto/config.h.in"
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

# PkgConfig
include( ${CMAKE_SOURCE_DIR}/cmake/FindPkgConfig.cmake )
if (NOT PKG_CONFIG_FOUND)
  message (FATAL_ERROR "Cannot find PkgConfig!")
endif()

# Dolfin via PkgCOnfig
pkg_search_module(DOLFIN REQUIRED dolfin)
if (NOT DOLFIN_FOUND)
  message (FATAL_ERROR "Cannot find libdolfin!")
endif()

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

if(CMAKE_COMPILER_IS_GNUCXX)
  STRING(REPLACE ";" " " DOLFIN_CXXFLAGS "${DOLFIN_CFLAGS_OTHER}")
  SET(CMAKE_CXX_FLAGS "-O3 -Wall ${DOLFIN_CXXFLAGS}")

  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -march=native -msse3 -fomit-frame-pointer -ffast-math -mfused-madd -std=gnu++0x -pthread")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -march=native -msse3 -fomit-frame-pointer -ffast-math -mfused-madd")
endif(CMAKE_COMPILER_IS_GNUCXX)

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

#message(${DOLFIN_CFLAGS_OTHER})

add_library(magnumfe SHARED ${MAGNUMFE_SRC})
target_link_libraries(magnumfe ${DOLFIN_LIBRARIES})
target_link_libraries(magnumfe ${GMSH_LIBRARY})

##############################################################
###    Swig Stuff                                          ###
##############################################################

#string(REGEX REPLACE "-Wall" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(SWIG_MODULE_NAME cpp)

set(CMAKE_SWIG_FLAGS "")
set(INTERFACE_FILES ${PROJECT_SOURCE_DIR}/magnumfe/swig/magnumfe.i)

set_source_files_properties (${INTERFACE_FILES} PROPERTIES CPLUSPLUS ON )
set_source_files_properties (${INTERFACE_FILES} PROPERTIES SWIG_FLAGS 
  -O
  )

#swig_add_module(cpp python ${INTERFACE_FILES} ${MAGNUMFE_SRC})
#swig_link_libraries(cpp ${DOLFIN_LIBRARIES})
swig_add_module(cpp python ${INTERFACE_FILES})
swig_link_libraries(cpp magnumfe)
