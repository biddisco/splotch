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

  this->RenderSplotchParams(ren, actor);

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

  int N = particle_data.size();
  
  vtkpiston::vtk_polydata *pd = (vtkpiston::vtk_polydata*)(tr->data);
  if (pd->userPointer) {
    particle_data.clear();
    particle_data.resize(N);
    for (int i=0; i<N; i++) {
      particle_data[i].e = COLOUR(0,0,0);
      particle_data[i].x = 0;
      particle_data[i].y = 0;
      particle_data[i].z = 0;
      particle_data[i].r = 0;
      particle_data[i].I = 1;
      particle_data[i].type = 0;
    }
    // raw data passed into splotch renderer internals
    cuda_paraview_rendering(mydevID, nTasksDev, pic, particle_data, campos, lookat, sky, amap, 1.0, params, pd->userPointer);
  }

  //
  this->PostRenderCompositing(ren, actor);

   {
    thrust::device_ptr<cu_particle_sim> devPtr((cu_particle_sim*)(pd->userPointer));
    thrust::device_vector<cu_particle_sim> dt_a(devPtr, devPtr + N);
    thrust::copy(dt_a.begin(), dt_a.end(), (cu_particle_sim*)(particle_data.data()));
   }

}
//-----------------------------------------------------------------------------

