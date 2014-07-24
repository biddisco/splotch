PROJECT("pv_splotch")

#--------------------------------------------------
# cmake old/new compatibility
#--------------------------------------------------
cmake_minimum_required(VERSION 2.8)

#--------------------------------------------------
# cmake policies
# see http://www.cmake.org/cmake/help/v2.8.12/cmake.html#section_Policies
#--------------------------------------------------
foreach(p
    CMP0003 
    CMP0017
    CMP0022 
    CMP0023 
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF (ParaView_SOURCE_DIR)
  IF(PARAVIEW_BUILD_QT_GUI)
    INCLUDE(${QT_USE_FILE})
  ENDIF(PARAVIEW_BUILD_QT_GUI)
  INCLUDE_DIRECTORIES(
    ${HDF5_INCLUDE_DIR}
    ${PARAVIEW_INCLUDE_DIRS}
    ${PARAVIEW_GUI_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

#
# For the Painter we need access to some IceT viewport info
# It's awful, but there's no other way to get the image reduction factor
#
INCLUDE_DIRECTORIES(${ParaView_BINARY_DIR}/Utilities/IceT/src/include)

#------------------------------------------------------
# MPI support if ParaView compiled with MPI
#------------------------------------------------------
if(PARAVIEW_USE_MPI)
  FIND_PACKAGE(MPI)   
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  ADD_DEFINITIONS(-DUSE_MPI)
  set(SPLOTCH_MPI_LIBS ${MPI_LIBRARY} ${MPI_C_LIBRARIES})
endif(PARAVIEW_USE_MPI)

if(PARAVIEW_USE_ICE_T)
  ADD_DEFINITIONS(-DUSE_ICET)
endif(PARAVIEW_USE_ICE_T)

ADD_DEFINITIONS(-DPARAVIEW_SPLOTCH)

# --------------------------------------------------
# We must build sub libraries with -fPIC
# --------------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#--------------------------------------------------
# if not called from the splotch source directory
# then build splotch
#--------------------------------------------------
if(NOT TARGET splotchcore)
  #-----------------------------------------------------------------------------
  # OpenMP
  #-----------------------------------------------------------------------------
  option(PV_SPLOTCH_USE_OPENMP "Compile with OpenMP support" ON)
  if(PV_SPLOTCH_USE_OPENMP)
    ADD_DEFINITIONS(-DHAVE_OPENMP)
    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
      SET(OpenMP_CXX_FLAGS "-fopenmp" )
    endif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
    if(WIN32 AND MSVC) 
      SET(OpenMP_CXX_FLAGS "/openmp")
    endif(WIN32 AND MSVC) 
    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
      SET(OpenMP_CXX_FLAGS "-fopenmp" )
    endif(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
  endif(PV_SPLOTCH_USE_OPENMP)

  #-----------------------------------------------------------------------------
  # CUDA
  #-----------------------------------------------------------------------------
  #--------------------------------------------------
  # pv-cuda-piston required for cuda support
  #--------------------------------------------------
  find_package(pv_cuda_piston NO_MODULE)
  if(NOT pv_cuda_piston_FOUND)
    set(CUDA_DEFAULT OFF)
  else()
    set(CUDA_DEFAULT ON)
    message("pv cuda piston was found")
  endif()

  option(PV_SPLOTCH_USE_CUDA "Compile with Cuda support" ${CUDA_DEFAULT})

  if(PV_SPLOTCH_USE_CUDA)
    # these should have already been found by pv-cuda-piston
    find_package(CUDA REQUIRED)

    #
    # use /MDd flag in CXX flags to solve some link issues on windows
    #
    if (NOT TARGET pv_cuda_piston)
      message(FATAL_ERROR "Cannot enable Cuda without cuda-piston plugin")
    endif()

    set(PV_SPLOTCH_DEP_LIBS pv_cuda_piston)

    add_definitions(-DPV_SPLOTCH_USE_PISTON)

    INCLUDE_DIRECTORIES(
      ${CMAKE_CURRENT_SOURCE_DIR}/piston
      ${pv_cuda_piston_INCLUDE_DIR}
    )

    #--------------------------------------------------
    # #defines for splotch
    #--------------------------------------------------
    add_definitions(-DCUDA -DHYPERQ -DSPLOTCH_PARAVIEW)

    #--------------------------------------------------
    # include path for Thrust
    #--------------------------------------------------
    include_directories(${CUDATHRUST_INCLUDE_DIRS})

    INCLUDE_DIRECTORIES(
      ${CMAKE_CURRENT_SOURCE_DIR}/../
      ${CMAKE_CURRENT_SOURCE_DIR}/../cuda
    )
      
    set(PV_SPLOTCH_CUDA_CSRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/cuda_splotch.cc
    )
    set(PV_SPLOTCH_CUDA_SRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/cuda_utils.cu
      ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/cuda_device_query.cu
      ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/cuda_policy.cu
#      ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/cuda_kernel.cu
      ${CMAKE_CURRENT_SOURCE_DIR}/../cuda/cuda_render.cu
    )
    set(PV_SPLOTCH_CUDA_PLUGIN_SRCS
      ${CMAKE_CURRENT_SOURCE_DIR}/piston/vtkCUDASplotchPainter.cxx
      ${CMAKE_CURRENT_SOURCE_DIR}/piston/vtkDataSetToSplotch.cxx
    )

    #--------------------------------------------------
    # Define CUDA compiled source files
    #--------------------------------------------------
    if(BUILD_SHARED_LIBS)
      set(GPGPU_BUILD_TYPE SHARED)
    endif()

    # set this, otherwise C_SOURCES is empty in Visual Studio  
    set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF) 

    # call the cuda wrapper
    cuda_compile(PV_SPLOTCH_CUDA_C_SRCS ${PV_SPLOTCH_CUDA_SRCS} ${GPGPU_BUILD_TYPE})
    message("Cuda srcs generated from ${PV_SPLOTCH_CUDA_SRCS} are ${PV_SPLOTCH_CUDA_C_SRCS}|")
  endif()

  #-----------------------------------------------------------------------------
  # Include Directories
  #-----------------------------------------------------------------------------
  INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/../
    ${CMAKE_CURRENT_SOURCE_DIR}/../kernel
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport
    ${CMAKE_CURRENT_SOURCE_DIR}/../splotch
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster
    ${CMAKE_CURRENT_SOURCE_DIR}/../c_utils
  )

  option(PV_SPLOTCH_DISABLE_READERS "Turn off compilation of readers ifparaview ones are used instead" ON)
  if(PV_SPLOTCH_DISABLE_READERS)
    message("adding def PARAVIEW_SPLOTCH_DISABLE_READERS")
    add_definitions(-DPARAVIEW_SPLOTCH_DISABLE_READERS)
  else (PV_SPLOTCH_DISABLE_READERS)
    INCLUDE_DIRECTORIES(
      ${CMAKE_CURRENT_SOURCE_DIR}/../reader
      ${CMAKE_CURRENT_SOURCE_DIR}/../writer
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
    )
    set(SPLOTCH_IO_LIBRARY splotchio)
    ADD_LIBRARY(${SPLOTCH_IO_LIBRARY} STATIC ${SPLOTCH_IO_LIB_SRCS})
    SET_TARGET_PROPERTIES(${SPLOTCH_IO_LIBRARY} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
    TARGET_LINK_LIBRARIES(${SPLOTCH_IO_LIBRARY} ${SPLOTCH_MPI_LIBS})
  endif()
  
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
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport/ls_image.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../cxxsupport/announce.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/m_rotation.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/mesh_creator.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/mesh_creator_vs1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/p_selector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/../booster/randomizer.cc
    ${PV_SPLOTCH_CUDA_CSRCS}
    ${PV_SPLOTCH_CUDA_C_SRCS}
  )

#  MESSAGE("Splotch Flags ${OpenMP_CXX_FLAGS}")
  ADD_LIBRARY(splotchcore STATIC ${SPLOTCH_CORE_LIB_SRCS})

  TARGET_LINK_LIBRARIES(splotchcore 
    ${SPLOTCH_IO_LIBRARY} 
    ${SPLOTCH_MPI_LIBS}     
    ${CUDA_LIBRARIES}
  )
  SET_TARGET_PROPERTIES(splotchcore PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

  add_executable(splotch 
    ${CMAKE_CURRENT_SOURCE_DIR}/../splotch/splotch.cc
  )

  target_link_libraries(splotch 
    splotchcore 
    ${SPLOTCH_IO_LIBRARY} 
    ${SPLOTCH_MPI_LIBS}     
    ${CUDA_LIBRARIES}
  )

endif(NOT TARGET splotchcore)

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
  ${PV_SPLOTCH_CUDA_PLUGIN_SRCS}
)

if(PARAVIEW_BUILD_QT_GUI)
  #--------------------------------------------------
  # Qt .ui files that are wrapped and added to plugin
  #--------------------------------------------------
#  set(UI_FORMS
#    pqSplotchDisplayPanelDecorator.ui
#  )
      
 # qt4_wrap_ui(UI_SRCS 
#    ${UI_FORMS}
#  )

#  qt4_wrap_cpp(MOC_SRCS 
#    pqSplotchDisplayPanelDecorator.h
#    )
     
#  set(Plugin_SRCS
#    pqSplotchDisplayPanelDecorator.cxx
#    )

#  add_paraview_display_panel_decorator(
#    IFACES
 #   IFACE_SRCS
 #   CLASS_NAME pqSplotchDisplayPanelDecorator
 #   PANEL_TYPES pqDisplayProxyEditor)

   qt4_wrap_cpp(MOC_SRCS
    pqSplotchPropertyWidgetDecorator.h)

  set(Plugin_SRCS
    pqSplotchPropertyWidgetDecorator.cxx
    )

  add_paraview_property_widget_decorator(IFACES IFACE_SRCS
    TYPE "splotch_decorator"
    CLASS_NAME pqSplotchPropertyWidgetDecorator)

endif(PARAVIEW_BUILD_QT_GUI)

#--------------------------------------------------
# Splotch uses MIP for LOD
#--------------------------------------------------
if(TARGET "pv_astro" AND EXISTS "${PROJECT_SOURCE_DIR}/../../pv-astro/")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../pv-astro)
  SET(PV_SPLOTCH_DEP_LIBS ${PV_SPLOTCH_DEP_LIBS} pv_astro)
  ADD_DEFINITIONS(-DPV_SPLOTCH_WITH_MIP)
endif(TARGET "pv_astro" AND EXISTS "${PROJECT_SOURCE_DIR}/../../pv-astro/")

#--------------------------------------------------
# Define Plugin
#--------------------------------------------------
SET(PLUGIN_NAME pv_splotch)

ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0" 

  SERVER_MANAGER_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSplotchRepresentation.xml
  SERVER_MANAGER_SOURCES
    ${splotch_plugin_WRAPPED_SRCS}
  SERVER_SOURCES
    ${splotch_plugin_SRCS} 
#    ${SPLOTCH_CORE_LIB_SRCS}
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

TARGET_LINK_LIBRARIES(${PLUGIN_NAME} 
  PUBLIC
    splotchcore 
    ${SPLOTCH_IO_LIBRARY} 
    ${PV_SPLOTCH_DEP_LIBS}
    ${MPI_LIBRARY}
)

SET_TARGET_PROPERTIES(${PLUGIN_NAME} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

#--------------------------------------------------------
# Create the UsePackage configuration for other projects
# that will be using this library
#--------------------------------------------------------
SET(splotch_plugin_INCLUDE_DIRS "${PROJECT_SOURCE_DIR};${PROJECT_BINARY_DIR}")

#--------------------------------------------------------
# Create the module import/export configuration header.
#--------------------------------------------------------
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/module_configure.h.in
               ${PROJECT_BINARY_DIR}/pv_splotch_configure.h
               @ONLY IMMEDIATE)

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
if(BUILD_TESTING)
  SET(PLUGIN_TEST_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)
  MAKE_DIRECTORY(${PLUGIN_TEST_DIR})
  ADD_SUBDIRECTORY(Testing)
endif(BUILD_TESTING)  

#--------------------------------------------------
# Tools : Splotch format converter
#--------------------------------------------------
if(NOT PV_SPLOTCH_DISABLE_READERS)
  if(EXISTS "${PROJECT_SOURCE_DIR}/../../pv-meshless/")
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../pv-meshless)
    ADD_EXECUTABLE(splotch2H5Part splotch2H5Part.cpp)
    TARGET_LINK_LIBRARIES(splotch2H5Part splotchcore splotchio pv_meshless)
  endif(EXISTS "${PROJECT_SOURCE_DIR}/../../pv-meshless/")
endif()

