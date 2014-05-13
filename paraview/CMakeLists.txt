PROJECT("SPLOTCH_plugin")

#--------------------------------------------------
# cmake old/new compatibility
#--------------------------------------------------
cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0017 NEW)
endif(COMMAND cmake_policy)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

#
# For the Painter we need access to some IceT viewport info
# It's awful, but there's no other way to get the image reduction factor
#
INCLUDE_DIRECTORIES(${ParaView_BINARY_DIR}/Utilities/IceT/src/include)

#------------------------------------------------------
# MPI support if ParaView compiled with MPI
#------------------------------------------------------
IF (PARAVIEW_USE_MPI)
  FIND_PACKAGE(MPI)   
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  ADD_DEFINITIONS(-DUSE_MPI)
ENDIF (PARAVIEW_USE_MPI)

IF(PARAVIEW_USE_ICE_T)
  ADD_DEFINITIONS(-DUSE_ICET)
ENDIF(PARAVIEW_USE_ICE_T)
#--------------------------------------------------
# If not called from the splotch source directory
# then build splotch
#--------------------------------------------------
IF (NOT TARGET splotchcore)
  #-----------------------------------------------------------------------------
  # OpenMP
  #-----------------------------------------------------------------------------
  OPTION(SPLOTCH_USE_OPENMP "Compile with OpenMP support" ON)
  IF (SPLOTCH_USE_OPENMP)
    ADD_DEFINITIONS(-DHAVE_OPENMP)
    if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
      SET(OpenMP_CXX_FLAGS "-fopenmp" )
    endif (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
    IF (WIN32 AND MSVC) 
      SET(OpenMP_CXX_FLAGS "/openmp")
    ENDIF (WIN32 AND MSVC) 
    if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
      SET(OpenMP_CXX_FLAGS "-fopenmp" )
    endif (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
  ENDIF (SPLOTCH_USE_OPENMP)

  #-----------------------------------------------------------------------------
  # Include Directories
  #-----------------------------------------------------------------------------
  INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/../
    ${CMAKE_CURRENT_SOURCE_DIR}/../kernel
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader
    ${CMAKE_CURRENT_SOURCE_DIR}/../writer
    ${CMAKE_CURRENT_SOURCE_DIR}/../splotch
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster          
  )

  #-----------------------------------------------------------------------------
  # Source files for splotch IO library
  #-----------------------------------------------------------------------------
  SET(SPLOTCH_IO_LIB_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader/mesh_reader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader/gadget_reader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader/millenium_reader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader/bin_reader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader/bin_reader_mpi.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader/tipsy_reader.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../reader/ramses_reader.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/../writer/write_tga.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/../writer/write_bmp.cc
#    ${CMAKE_CURRENT_SOURCE_DIR}/../writer/write_ppm.cc
  )
  add_definitions(-DPARAVIEW_SPLOTCH)

  ADD_LIBRARY(splotchio STATIC ${SPLOTCH_IO_LIB_SRCS})
  SET_TARGET_PROPERTIES(splotchio PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
  TARGET_LINK_LIBRARIES(splotchio ${SPLOTCH_MPI_LIBS})
  
  #-----------------------------------------------------------------------------
  # Source files for splotch core library
  #-----------------------------------------------------------------------------
  SET(SPLOTCH_CORE_LIB_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/../kernel/transform.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/../splotch/splotchutils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../splotch/scenemaker.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../splotch/splotch_host.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport/error_handling.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport/walltimer.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport/mpi_support.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport/paramfile.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport/string_utils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/m_rotation.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/mesh_creator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/mesh_creator_vs1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/p_selector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/randomizer.cc 
  )

#  MESSAGE("Splotch Flags ${OpenMP_CXX_FLAGS}")
  ADD_LIBRARY(splotchcore STATIC ${SPLOTCH_CORE_LIB_SRCS})
  TARGET_LINK_LIBRARIES(splotchcore splotchio ${SPLOTCH_MPI_LIBS})
  SET_TARGET_PROPERTIES(splotchcore PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

ENDIF (NOT TARGET splotchcore)

#--------------------------------------------------
# Set project include directories 
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

#--------------------------------------------------
# Source files, that are wrapped by Client/Server
#--------------------------------------------------
SET( splotch_plugin_WRAPPED_SRCS
  vtkSplotchRepresentation.cxx  
  vtkSplotchDefaultPainter.cxx
  vtkSplotchPainter.cxx
)

IF(PARAVIEW_BUILD_QT_GUI)
  #--------------------------------------------------
  # Qt .ui files that are wrapped and added to plugin
  #--------------------------------------------------
  set(UI_FORMS
    pqSplotchDisplayPanelDecorator.ui
  )
      
  qt4_wrap_ui(UI_SRCS 
    ${UI_FORMS}
  )

  qt4_wrap_cpp(MOC_SRCS 
    pqSplotchDisplayPanelDecorator.h
    pqSplotchColorScaleEditor.h
    )
     
  set(Plugin_SRCS
    pqSplotchDisplayPanelDecorator.cxx
    pqSplotchColorScaleEditor.cxx
    )

  add_paraview_display_panel_decorator(
    IFACES
    IFACE_SRCS
    CLASS_NAME pqSplotchDisplayPanelDecorator
    PANEL_TYPES pqDisplayProxyEditor)

ENDIF(PARAVIEW_BUILD_QT_GUI)

#--------------------------------------------------
# Splotch uses MIP for LOD
#--------------------------------------------------
IF (TARGET "pv_astro" AND EXISTS "${PROJECT_SOURCE_DIR}/../../pv-astro/")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../pv-astro)
  SET(DEP_LIBS pv_astro)
  ADD_DEFINITIONS(-DPV_SPLOTCH_WITH_MIP)
ENDIF (TARGET "pv_astro" AND EXISTS "${PROJECT_SOURCE_DIR}/../../pv-astro/")

#--------------------------------------------------
# Define Plugin
#--------------------------------------------------
SET(LIBRARY_TARGET_NAME pv_splotch)

ADD_PARAVIEW_PLUGIN(
  ${LIBRARY_TARGET_NAME}
  "1.0" 

  SERVER_MANAGER_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSplotchRepresentation.xml
  SERVER_MANAGER_SOURCES
    ${splotch_plugin_WRAPPED_SRCS}
  SERVER_SOURCES
    ${splotch_plugin_SRCS} 
  GUI_RESOURCE_FILES
    ${splotch_plugin_CLIENT_XML}
  GUI_INTERFACES 
    ${IFACES} 
  GUI_SOURCES 
    ${Plugin_SRCS} 
    ${MOC_SRCS} 
    ${UI_SRCS} 
    ${IFACE_SRCS}
)

TARGET_LINK_LIBRARIES(pv_splotch splotchcore splotchio ${DEP_LIBS})
SET_TARGET_PROPERTIES(pv_splotch PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

#--------------------------------------------------------
# Create the UsePackage configuration for other projects
# that will be using this library
#--------------------------------------------------------
SET(splotch_plugin_INCLUDE_DIRS "${PROJECT_SOURCE_DIR};${PROJECT_BINARY_DIR}")

#--------------------------------------------------
# Install
#--------------------------------------------------
SET(INSTALL_PATH 
  "${CMAKE_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}"
)

INSTALL(
  TARGETS 
    pv_splotch
  DESTINATION ${INSTALL_PATH}
)

#--------------------------------------------------
# Testing
#--------------------------------------------------
IF (BUILD_TESTING)
  SET(PLUGIN_TEST_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)
  MAKE_DIRECTORY(${PLUGIN_TEST_DIR})
  ADD_SUBDIRECTORY(Testing)
ENDIF (BUILD_TESTING)  

#--------------------------------------------------
# Tools : Splotch format converter
#--------------------------------------------------
IF(EXISTS "${PROJECT_SOURCE_DIR}/../../pv-meshless/")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../pv-meshless)
  ADD_EXECUTABLE(splotch2H5Part splotch2H5Part.cpp)
  TARGET_LINK_LIBRARIES(splotch2H5Part splotchcore splotchio pv_meshless)
ENDIF(EXISTS "${PROJECT_SOURCE_DIR}/../../pv-meshless/")


