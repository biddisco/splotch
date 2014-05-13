#--------------------------------------------------
# Tests here are added using the 
# ADD_TEST(NAME ... COMMAND ...)
# form and will only run when ctest is run using
# ctest -C Release / ctest -C Debug / etc
#--------------------------------------------------

#--------------------------------------------------
# Test splotch raytracer
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../pv-meshless
)


# --------------------------------------------------
# Optional Splotch+HDF5 test
# --------------------------------------------------
IF (USE_HDF5)
  ADD_EXECUTABLE(TestSplotchMapper TestSplotchMapper.cxx)
  TARGET_LINK_LIBRARIES(TestSplotchMapper
    pv_meshless
    pv_splotch
    splotchio
    )

  ADD_TEST(
    NAME TestSplotchMapper 
    COMMAND 
#      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
      $<TARGET_FILE:TestSplotchMapper> 
      -D ${PROJECT_SOURCE_DIR}/Testing/data/ 
      -F splotch.h5part 
      -T ${PLUGIN_TEST_DIR}
      -R
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/TestSplotchMapper.png
  )  
ENDIF (USE_HDF5)
