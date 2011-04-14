
#--------------------------------------------------
# Trilinos
#--------------------------------------------------
FIND_PACKAGE(Trilinos)
MESSAGE("Trilinos_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS}")
INCLUDE_DIRECTORIES("${Trilinos_INCLUDE_DIRS}")

LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS}/${CMAKE_CFG_INTDIR})
MESSAGE("${Trilinos_LIBRARY_DIRS}/${CMAKE_CFG_INTDIR}")

#--------------------------------------------------
# Test partitioner/raytracer
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../../../pv-meshless
)

ADD_EXECUTABLE(TestParticlePartition TestParticlePartition.cxx)
TARGET_LINK_LIBRARIES(TestParticlePartition
  splotchio
  pv_splotch
  pv_meshless
  ${Trilinos_LIBRARIES}
  )
ADD_TEST(TestParticlePartition.cxx ${PROJECT_BINARY_DIR}/../bin/TestParticlePartition.cxx)
