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
    if (this->particle_compute) {
      this->DataSetToPiston->Modified();
      // before writing particles to the GPU, re-initialize splotch GPU vars
      cuda_paraview_init(pic, particle_data, campos, lookat, sky, 1.0, params);
    }
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


//  COLOURMAP c1;
//  c1.addVal(0.0, COLOUR(0.0, 0.0, 1.0));
//  c1.addVal(0.5, COLOUR(0.5, 1.0, 0.5));
//  c1.addVal(1.0, COLOUR(1.0, 0.0, 0.0));

//  std::vector<COLOURMAP> amap;
//  amap.push_back(c1);

  vtkPistonDataObject *id = this->DataSetToPiston->GetPistonDataObjectOutput(0);
  if (!id) {
    return;
  }
  vtkPistonReference *tr = id->GetReference();

  vtkpiston::vtk_polydata *pd = (vtkpiston::vtk_polydata*)(tr->data);
  if (pd->userPointer) {
    // raw data passed into splotch renderer internals
    cuda_paraview_rendering(mydevID, nTasksDev, pic, particle_data, campos, lookat, sky, 1.0, params, pd->userPointer);
  }

  //
  this->PostRenderCompositing(ren, actor);
}
//-----------------------------------------------------------------------------

