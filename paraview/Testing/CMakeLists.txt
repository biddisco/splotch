#--------------------------------------------------
# Test splotch raytracer
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../pv-meshless
)

# --------------------------------------------------
# Optional Trilinos support
# This needs to be fixed so that pv-meshless pulls in 
# the libs and other plugins automatically get linked correctly
# --------------------------------------------------
OPTION(USE_TRILINOS "Use Trilinos library" OFF)
IF(USE_TRILINOS)
  FIND_PACKAGE(Trilinos)
  IF(${Trilinos_FOUND})
    INCLUDE_DIRECTORIES("${Trilinos_INCLUDE_DIRS}")
    LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS})
    SET(TRILINOS_CXX ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.cxx)
    SET(TRILINOS_XML ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.xml)
  ENDIF(${Trilinos_FOUND})
ENDIF(USE_TRILINOS)

# --------------------------------------------------
# Optional Trilinos support 
# --------------------------------------------------
IF (USE_HDF5)
  ADD_EXECUTABLE(TestSplotchMapper TestSplotchMapper.cxx)
  TARGET_LINK_LIBRARIES(TestSplotchMapper
    splotchio
    pv_splotch
    pv_meshless
    )
  ADD_TEST(TestSplotchMapper.cxx ${PROJECT_BINARY_DIR}/../bin/TestSplotchMapper.cxx)
ENDIF (USE_HDF5)
