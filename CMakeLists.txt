PROJECT(splotch)
cmake_minimum_required(VERSION 2.8)
  
#-----------------------------------------------------------------------------
# Setup output Directories
#-----------------------------------------------------------------------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Executables."
)
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Libraries"
)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all static libraries."
)

#-----------------------------------------------------------------------------
# OpenMP
#-----------------------------------------------------------------------------
OPTION(SPLOTCH_USE_OPENMP "Compile with OpenMP support" ON)
IF (SPLOTCH_USE_OPENMP)
  ADD_DEFINITIONS(-DHAVE_OPENMP)
  SET(OpenMP_CXX_FLAGS "")
  IF (WIN32 AND MSVC) 
    SET(OpenMP_CXX_FLAGS "/openmp")
  ENDIF (WIN32 AND MSVC) 
ENDIF (SPLOTCH_USE_OPENMP)

#-----------------------------------------------------------------------------
# MPI
#-----------------------------------------------------------------------------
OPTION(SPLOTCH_USE_MPI "Compile with MPI support" ON)

SET(SPLOTCH_MPI_LIBS "")
IF (SPLOTCH_USE_MPI)
  FIND_PACKAGE(MPI)
  ADD_DEFINITIONS(-DUSE_MPI)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH}) 
  SET(SPLOTCH_MPI_LIBS "${MPI_LIBRARY}")
  IF (MPI_EXTRA_LIBRARY)
    SET(SPLOTCH_MPI_LIBS "${MPI_LIBRARY};${MPI_EXTRA_LIBRARY}")
  ENDIF (MPI_EXTRA_LIBRARY)
  ADD_DEFINITIONS(-DSUPPRESS_MPI_INIT)
ENDIF (SPLOTCH_USE_MPI)

#-----------------------------------------------------------------------------
# Include Directories
#-----------------------------------------------------------------------------
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  kernel
  cxxsupport
  reader
  writer
  splotch
)

#-----------------------------------------------------------------------------
# Source files for splotch core library
#-----------------------------------------------------------------------------
SET(SPLOTCH_CORE_LIB_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/kernel/transform.cxx
  ${CMAKE_CURRENT_SOURCE_DIR}/splotch/splotchutils.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/splotch/scenemaker.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/splotch/splotch_host.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cxxsupport/error_handling.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cxxsupport/walltimer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cxxsupport/mpi_support.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cxxsupport/cxxutils.cc
)

ADD_LIBRARY(splotchcore STATIC ${SPLOTCH_CORE_LIB_SRCS})
SET_TARGET_PROPERTIES(splotchcore PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
  
#-----------------------------------------------------------------------------
# Source files for splotch IO library
#-----------------------------------------------------------------------------
SET(SPLOTCH_IO_LIB_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/reader/mesh_reader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/reader/gadget_reader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/reader/millenium_reader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/reader/bin_reader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/reader/bin_reader_mpi.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/writer/write_tga.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/writer/write_bmp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/writer/write_ppm.cc
)

ADD_LIBRARY(splotchio STATIC ${SPLOTCH_IO_LIB_SRCS})
SET_TARGET_PROPERTIES(splotchio PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
TARGET_LINK_LIBRARIES(splotchio splotchcore ${SPLOTCH_MPI_LIBS})
  
#-----------------------------------------------------------------------------
# ParaView plugin
#-----------------------------------------------------------------------------
OPTION(SPLOTCH_BUILD_PARAVIEW_PLUGIN "Build a plugin for rendering with splotch inside paraview" OFF)

IF (SPLOTCH_BUILD_PARAVIEW_PLUGIN)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
  INCLUDE(${PARAVIEW_VTK_DIR}/VTKConfig.cmake)

  IF (ParaView_FOUND)
    ADD_SUBDIRECTORY(paraview)
  ENDIF (ParaView_FOUND)
ENDIF (SPLOTCH_BUILD_PARAVIEW_PLUGIN)

#-----------------------------------------------------------------------------
# Define the Test Executable and set properties
#-----------------------------------------------------------------------------
SET(SPLOTCH_SRCS
  splotch/splotch.cc
)

ADD_EXECUTABLE(splotchdemo ${SPLOTCH_SRCS})
SET_TARGET_PROPERTIES(splotchdemo PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
TARGET_LINK_LIBRARIES(splotchdemo splotchcore splotchio)
