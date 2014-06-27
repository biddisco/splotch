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

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vtkCUDASplotchPainter::vtkCUDASplotchPainter()
{
  this->DataSetToPiston  = vtkSmartPointer<vtkDataSetToSplotch>::New();
  //
 }
//-----------------------------------------------------------------------------
vtkCUDASplotchPainter::~vtkCUDASplotchPainter()
{
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
//    this->DataSetToPiston->SetTypeArrayName(this->TypeScalars);
    this->DataSetToPiston->SetBrightness(&(this->Brightness[0]));
    this->DataSetToPiston->SetRadiusMultiplier(this->RadiusMultiplier);
    this->DataSetToPiston->SetParticleData(this->particle_data);
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

