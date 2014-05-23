/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCUDASplotchPainter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCUDASplotchPainter.h"

#include "vtkgl.h"
#include "vtkOpenGLExtensionManager.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkTransform.h"
#include "vtkPainterDeviceAdapter.h"
//
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkScalarsToColorsPainter.h"
//
#include "vtkObjectFactory.h"
#include "vtkTimerLog.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositeDataIterator.h"
//
#include <thrust/version.h>
#include <thrust/copy.h>
#include <stdexcept>
//
#include "cuda/splotch_cuda2.h"
//
#include <piston/choose_container.h>
//
#include "vtkPistonDataObject.h"
#include "vtkPistonDataWrangling.h"
#include "vtkPistonReference.h"

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCUDASplotchPainter);
//-----------------------------------------------------------------------------
#define NUM_INTEROP_BUFFERS 4
bool vtkCUDASplotchPainter::CudaGLInitted = false;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class vtkCUDASplotchPainter::InternalInfo {
public:
  InternalInfo() {
    this->BufferSize = 0;
    this->CellCount = 0;
    this->DataObjectMTimeCache = 0;
    this->clearBuffers();
  }
  ~InternalInfo() {
    this->clearBuffers();
  }
  void clearBuffers() {
    if (this->BufferSize != 0 && this->vboBuffers[0] != -1) {
      vtkgl::DeleteBuffers(NUM_INTEROP_BUFFERS, this->vboBuffers);
    }
    for (int i=0; i<NUM_INTEROP_BUFFERS; i++) {
      vboBuffers[i] = -1;
      vboResources[i] = NULL;
    }
  }

  int     BufferSize;
  int     CellCount;
  GLuint  vboBuffers[NUM_INTEROP_BUFFERS];
  struct  cudaGraphicsResource* vboResources[NUM_INTEROP_BUFFERS];
  //
  unsigned long DataObjectMTimeCache;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vtkCUDASplotchPainter::vtkCUDASplotchPainter()
{
  this->DataSetToPiston  = vtkSmartPointer<vtkDataSetToSplotch>::New();
  //
  this->Internal = new vtkCUDASplotchPainter::InternalInfo();
}
//-----------------------------------------------------------------------------
vtkCUDASplotchPainter::~vtkCUDASplotchPainter()
{
  this->PrepareDirectRenderBuffers(0, 0);
  delete this->Internal;
}
//-----------------------------------------------------------------------------
void vtkCUDASplotchPainter::PrepareForRendering(vtkRenderer* renderer, vtkActor* actor)
{
  // call superclass function
  this->vtkSplotchPainter::PrepareForRendering(renderer, actor);
  //
  // if we are not using CUDA exit immediately
  if (!this->EnableCUDA) {
    return;
  }
  // if the input dataset is invalid exit immediately
  vtkDataObject* input = this->GetInput();
  if (!input)
  {
    vtkErrorMacro("No input present.");
    return;
  }

  // Now if we have composite data, we need to MapScalars for all leaves.
  if (input->IsA("vtkCompositeDataSet"))
  {
    throw std::runtime_error("Not supported");
  }
  else
  {
    this->DataSetToPiston->SetInputData(input);
    this->DataSetToPiston->SetRadiusArrayName(this->GetRadiusScalars(0));
    this->DataSetToPiston->SetIntensityArrayName(this->GetIntensityScalars(0));
    this->DataSetToPiston->SetScalarArrayName(this->ArrayName);
    this->DataSetToPiston->Update();
  }
  //
  this->Camera = renderer->GetActiveCamera();
  this->Actor  = actor;
}
//-----------------------------------------------------------------------------
void vtkCUDASplotchPainter::PrepareDirectRenderBuffers(int nPoints, int nCells)
{
  if (nPoints==this->Internal->BufferSize && nCells==this->Internal->CellCount) {
    return;
  }
  if (this->Internal->BufferSize != 0) {
    this->Internal->clearBuffers();
  }

  this->Internal->BufferSize = nPoints;
  this->Internal->CellCount  = nCells;
  if (this->Internal->BufferSize == 0) {
    return;
  }

  // Prep shared mem buffer between gl and cuda
  vtkgl::GenBuffers(NUM_INTEROP_BUFFERS, this->Internal->vboBuffers);

  // points 3*n float {x,y,z}
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[0]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*3*sizeof(float), 0,
    vtkgl::DYNAMIC_DRAW);

  // normals 3*n float {n0,n1,n2}
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[1]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*3*sizeof(float), 0,
    vtkgl::DYNAMIC_DRAW);

  // colors 4*n uchar {R,G,B,A} 
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[2]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*4*sizeof(unsigned char), 0,
    vtkgl::DYNAMIC_DRAW);

  // indexes 3*nCells int : triangles assumed {a,b,c} 
  vtkgl::BindBuffer(vtkgl::ELEMENT_ARRAY_BUFFER,
    this->Internal->vboBuffers[3]);
  vtkgl::BufferData(vtkgl::ELEMENT_ARRAY_BUFFER,
    this->Internal->CellCount*3*sizeof(int), 0,
    vtkgl::DYNAMIC_DRAW);

  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[0],
    this->Internal->vboBuffers[0]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[1],
    this->Internal->vboBuffers[1]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[2],
    this->Internal->vboBuffers[2]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[3],
    this->Internal->vboBuffers[3]);
}
// ---------------------------------------------------------------------------
void vtkCUDASplotchPainter::RenderInternal(vtkRenderer* ren, vtkActor* actor, 
  unsigned long typeflags, bool forceCompileOnly)
{
  if (!this->EnableCUDA) {
    this->vtkSplotchPainter::RenderInternal(ren, actor, typeflags, forceCompileOnly);
    return;
  }

  //
  int mydevID = 0;
  int nTasksDev = 1;

  paramfile params;
  params.find("ptypes", this->NumberOfParticleTypes);
  params.find("xres", X);
  params.find("yres", Y);
  for (int i=0; i<this->NumberOfParticleTypes; i++) {
    std::string name;
    name = "intensity_log" + NumToStrSPM<int>(i);
    params.find(name, (this->LogIntensity[i]!=0));
    name = "brightness" + NumToStrSPM<int>(i);
    params.find(name, this->Brightness[i]);
  }
  params.find("gray_absorption", this->GrayAbsorption);
  params.find("zmin", zmin);
  params.find("zmax", zmax);
  params.find("fov",  splotchFOV);
  params.find("projection", true);
  params.find("minrad_pix", 1);
  params.find("a_eq_e", true);
  params.find("colorbar", false);
  params.find("quality_factor", 0.001);
  params.find("boost", false);

//  params.find("intensity_min0", -11.8784);
//  params.find("intensity_max0",  -1.44456);

  params.find("color_min0",  0.152815);
  params.find("color_max0",  6.29244);

  COLOURMAP c1;
  c1.addVal(0.0, COLOUR(0.0, 0.0, 1.0));
  c1.addVal(0.5, COLOUR(0.5, 1.0, 0.5));
  c1.addVal(1.0, COLOUR(1.0, 0.0, 0.0));

  std::vector<COLOURMAP> amap;
  amap.push_back(c1);
  //  amap.push_back(c1);

  vtkPistonDataObject *id = this->DataSetToPiston->GetPistonDataObjectOutput(0);
  if (!id) {
    return;
  }
  vtkPistonReference *tr = id->GetReference();
  
  vtkpiston::vtk_polydata *pd = (vtkpiston::vtk_polydata*)(tr->data);
  if (pd->userPointer) {
    // raw data passed into splotch renderer internals
    cuda_paraview_rendering(mydevID, nTasksDev, pic, particle_data, campos, lookat, sky, amap, 1.0, params);
  }

  //
  this->PostRenderCompositing(ren, actor);
}
//-----------------------------------------------------------------------------
void vtkCUDASplotchPainter::RenderOnGPU(vtkCamera *cam, vtkActor *act)
{
  vtkPistonDataObject *id = this->DataSetToPiston->GetPistonDataObjectOutput(0);
  if (!id) {
    return;
  }

  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  int nPoints = vtkpiston::QueryNumVerts(id);
  int nCells = vtkpiston::QueryNumCells(id);
  this->PrepareDirectRenderBuffers(nPoints, nCells);

  // Transfer what is in tdo to buffer and render it directly on card
  bool hasNormals = false;
  bool hasColors = false;
  bool useindexbuffers = false;
/*
  double cameravec[3], origin[3];
  this->ComputeProjectionVector(cam, act, cameravec, origin);
  if (this->Direction>=0) {
//    vtkpiston::DepthSortPolygons(id, cameravec, this->Direction);
  }
*/

/*
  double scalarrange[2];
  this->ScalarsToColors->GetScalarRange(scalarrange);
  vtkpiston::CudaTransferToGL(
    id, this->Internal->DataObjectMTimeCache,
    this->Internal->vboResources, 
    this->ScalarsToColors->GetRGBAPointer(),
    scalarrange, 
    act->GetProperty()->GetOpacity(),
    hasNormals, hasColors, useindexbuffers);
*/
  glPushAttrib(GL_ENABLE_BIT);

  // Draw the result
  glEnableClientState(GL_VERTEX_ARRAY);
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[0]);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  if (hasNormals)
  {
    glEnableClientState(GL_NORMAL_ARRAY);
    vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[1]);
    glNormalPointer(GL_FLOAT, 0, 0);
  }

  if (hasColors)
  {
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
//    glEnable(GL_BLEND);
//    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
//    glEnable(GL_COLOR_MATERIAL);
    /*
    // because we have used premultiple_with_alpha
    // otherwise we'd use glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    glDepthFunc( GL_LEQUAL );
    glEnable( GL_DEPTH_TEST );
*/

    glEnable(GL_DEPTH_TEST);
    glEnableClientState(GL_COLOR_ARRAY);
    vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[2]);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);

  }
  else {
    //    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    //    glEnable(GL_COLOR_MATERIAL);
    //    glDepthFunc( GL_LEQUAL );
    //    glEnable( GL_DEPTH_TEST );
    //    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //    glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glEnable(GL_BLEND);
    //glEnableClientState(GL_COLOR_ARRAY);
    //vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[2]);
    //glColorPointer(4, GL_FLOAT, 0, 0);
  }

  if (useindexbuffers) {
    //
    int vertsPer = vtkpiston::QueryVertsPer(id);
    vtkgl::BindBuffer(vtkgl::ELEMENT_ARRAY_BUFFER, this->Internal->vboBuffers[3]);
    switch (vertsPer) {
    case 4:
      glDrawElements(GL_QUADS, nCells*4, GL_UNSIGNED_INT, (GLvoid*)0);
      break;
    case 3:
      glDrawElements(GL_TRIANGLES, nCells*3, GL_UNSIGNED_INT, (GLvoid*)0);
      break;
    default:
      glDrawElements(GL_POINTS, nCells*1, GL_UNSIGNED_INT, (GLvoid*)0);
    }
  }
  else {
    int vertsPer = vtkpiston::QueryVertsPer(id);
    switch (vertsPer) {
    case 4:
      glDrawArrays(GL_QUADS, 0, nPoints);
      break;
    case 3:
      glDrawArrays(GL_TRIANGLES, 0, nPoints);
      break;
    default:
      glDrawArrays(GL_POINTS, 0, nPoints);
    }
  }

  glDisableClientState(GL_VERTEX_ARRAY);
  if (hasNormals) glDisableClientState(GL_NORMAL_ARRAY);
  if (hasColors) glDisableClientState(GL_COLOR_ARRAY);

  glPopAttrib();
  // Update object modified time
  this->Internal->DataObjectMTimeCache = id->GetMTime();

  timer->StopTimer();
  double rendertime = timer->GetElapsedTime();
  //  std::cout << setprecision(6) << "RenderTime : << " <<  rendertime << std::endl;
}
//-----------------------------------------------------------------------------

