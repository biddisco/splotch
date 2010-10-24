PROJECT("SPLOTCH_plugin")

CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF(NOT PARAVIEW_LOADED)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
  SET(PARAVIEW_LOADED 1)
ELSE(NOT PARAVIEW_LOADED)
  IF (PARAVIEW_BUILD_QT_GUI)
    INCLUDE(${QT_USE_FILE})
  ENDIF (PARAVIEW_BUILD_QT_GUI)
  INCLUDE_DIRECTORIES(
    ${PARAVIEW_INCLUDE_DIRS}
    ${PARAVIEW_GUI_INCLUDE_DIRS}
    ${PARAVIEW_KWSYS_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
  )
ENDIF(NOT PARAVIEW_LOADED)

#--------------------------------------------------
# Set project include directories 
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

#--------------------------------------------------
# Source files, not wrapped by Client/Server
#--------------------------------------------------
SET(splotch_plugin_SRCS
#  vtkSplotchRaytracePainter.cxx
#  vtkSplotchRaytracePainter.h
#  vtkSplotchRaytracePreprocessor.cxx
#  vtkSplotchRaytracePreprocessor.h
#  vtkSplotchRaytraceDefaultPainter.cxx
#  vtkSplotchRaytraceDefaultPainter.h
#  vtkSplotchRaytraceProperty.cxx
#  vtkSplotchRaytraceProperty.h
#  vtkSMSplotchRepresentation.cxx

  framebufferObject.cpp          
  glErrorUtil.cpp                
  renderbuffer.cpp               
  vtkPointSpriteMapper.cxx            
  vtkPointSpriteDepthSort.cxx  
)

#--------------------------------------------------
# Source files, that are wrapped by Client/Server
#--------------------------------------------------
SET( splotch_plugin_WRAPPED_SRCS
  vtkSplotchRaytraceMapper.cxx
  vtkSplotchRaytraceRepresentation.cxx
)
#--------------------------------------------------
# Define Plugin
#--------------------------------------------------
SET(LIBRARY_TARGET_NAME splotch_plugin)

ADD_PARAVIEW_PLUGIN(
  ${LIBRARY_TARGET_NAME}
  "1.0" 

  SERVER_MANAGER_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/splotch_plugin_server.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSplotchRaytraceRepresentation.xml
  SERVER_MANAGER_SOURCES
    ${splotch_plugin_WRAPPED_SRCS}
  SERVER_SOURCES
    ${splotch_plugin_SRCS}  
  GUI_RESOURCE_FILES
    ${splotch_plugin_CLIENT_XML}
)

TARGET_LINK_LIBRARIES(splotch_plugin splotchlib)

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
    splotch_plugin
  DESTINATION ${INSTALL_PATH}
)

