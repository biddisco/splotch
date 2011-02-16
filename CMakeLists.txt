PROJECT(splotch)
cmake_minimum_required(VERSION 2.8)
  
#-----------------------------------------------------------------------------
# OpenMP
#-----------------------------------------------------------------------------
OPTION(SPLOTCH_USE_OPENMP "Compile with OpenMP support" ON)
IF (SPLOTCH_USE_OPENMP)
  ADD_DEFINITIONS(-DHAVE_OPENMP)
  SET(OpenMP_CXX_FLAGS "")
  IF (WIN32 AND MSVC) 
    SET(OpenMP_CXX_FLAGS "/openmp")
    MESSAGE("I just set flags here")
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
  SET(SPLOTCH_MPI_LIBS "${MPI_LIBRARY};${MPI_EXTRA_LIBRARY}")
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
# Source files
#-----------------------------------------------------------------------------
SET(SPLOTCH_SRCS
	kernel/transform.cxx
	cxxsupport/error_handling.cc
	cxxsupport/walltimer.cc
	reader/mesh_reader.cc
	cxxsupport/mpi_support.cc
	cxxsupport/cxxutils.cc
	reader/gadget_reader.cc
	reader/millenium_reader.cc
	reader/bin_reader.cc
	reader/bin_reader_mpi.cc
	writer/write_tga.cc
	writer/write_bmp.cc
	writer/write_ppm.cc
	splotch/splotchutils.cc
	splotch/splotch.cc
	splotch/scenemaker.cc
	splotch/splotch_host.cc
)

#-----------------------------------------------------------------------------
# Define the Executable and any properties
#-----------------------------------------------------------------------------
ADD_EXECUTABLE(splotch ${SPLOTCH_SRCS})
TARGET_LINK_LIBRARIES(splotch ${SPLOTCH_MPI_LIBS})
SET_TARGET_PROPERTIES(splotch PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
MESSAGE("The flags are ${OpenMP_CXX_FLAGS}")
