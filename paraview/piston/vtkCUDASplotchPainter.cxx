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

  params.find("intensity_min0", -11.8784);
  params.find("intensity_max0",  -1.44456);

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
    cuda_paraview_rendering(mydevID, nTasksDev, pic, particle_data, campos, lookat, sky, amap, 1.0, params, pd->userPointer);
  }

  //
  this->PostRenderCompositing(ren, actor);
}
//-----------------------------------------------------------------------------

