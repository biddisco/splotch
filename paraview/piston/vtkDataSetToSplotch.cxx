/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDataSetToSplotch.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkDataSetToSplotch.h"

//
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/version.h>
#include <vector>
//
#include "vtkPistonDataObject.h"
//
#include "vtkPistonDataWrangling.h"
#include "vtkPistonReference.h"
#include "vtkCUDAPiston.h"
//
//#include "cuda/splotch_cuda.h"
struct cu_color
  {
  float r,g,b;
  };

struct cu_particle_sim
  {
    cu_color e;
    float x,y,z,r,I;
    unsigned short type;
    bool active;
  };
//
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkDataSetToSplotch);

//----------------------------------------------------------------------------
vtkDataSetToSplotch::vtkDataSetToSplotch()
{
  this->ScalarArrayName    = NULL;
  this->RadiusArrayName    = NULL;
  this->IntensityArrayName = NULL;
}

//----------------------------------------------------------------------------
vtkDataSetToSplotch::~vtkDataSetToSplotch()
{
  delete []this->ScalarArrayName;
  delete []this->RadiusArrayName;
  delete []this->IntensityArrayName;
}

//----------------------------------------------------------------------------
void vtkDataSetToSplotch::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkDataSetToSplotch::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//-----------------------------------------------------------------------------
bool CheckDirty(vtkDataSet *ds, vtkPistonReference *tr)
{
  unsigned long int dstime = ds->GetMTime();
  if (dstime != tr->mtime)
  {
    tr->mtime = dstime;
    return true;
  }
  return false;
}

//------------------------------------------------------------------------------
int vtkDataSetToSplotch::RequestData(vtkInformation *request,
                                     vtkInformationVector** inputVector,
                                     vtkInformationVector* outputVector)
{
  vtkPistonDataObject *od = vtkPistonDataObject::GetData(outputVector);

  vtkDataObject *ido = this->GetInputDataObject(0,0);
  vtkDataSet *ds = vtkDataSet::SafeDownCast(ido);
  //
  if (ido->GetDataObjectType()!=VTK_POLY_DATA) {
    vtkErrorMacro("Only polydata supported for splotch rendering");
    return 0;
  }
  //
  od->SetBounds(ds->GetBounds());
  vtkPolyData *id = vtkPolyData::GetData(inputVector[0]);
  int nPoints = id->GetNumberOfPoints();
  //
  //
  //
  vtkpiston::AllocGPU(id,od);
  vtkPistonReference *tr = od->GetReference();
  vtkpiston::vtk_polydata *newD = (vtkpiston::vtk_polydata*)(tr->data);
  //
  // Scalars
  //
  vtkDataArray        *inscalars = id->GetPointData()->GetArray(this->ScalarArrayName);
  vtkDataArray         *inradius = id->GetPointData()->GetArray(this->RadiusArrayName);
  vtkDataArray      *inintensity = id->GetPointData()->GetArray(this->IntensityArrayName);
  vtkUnsignedCharArray *incolors = vtkUnsignedCharArray::SafeDownCast(
    id->GetPointData()->GetArray("Color"));
  //
  thrust::host_vector<cu_particle_sim> cuda_particles(nPoints);
  for (int i=0; i<nPoints; i++) {
    double *nextP = id->GetPoint(i);
    cuda_particles[i].x = (float)nextP[0];
    cuda_particles[i].y = (float)nextP[1];
    cuda_particles[i].z = (float)nextP[2];
    //
    if (incolors) {
      double *nextC = incolors->GetTuple4(i);
      cuda_particles[i].e.r = (float)nextC[0];
      cuda_particles[i].e.g = (float)nextC[1];
      cuda_particles[i].e.b = (float)nextC[2];
    }
    else {
      cuda_particles[i].e.r = 0.5;
      cuda_particles[i].e.g = 0.5;
      cuda_particles[i].e.b = 0.5;
    }
    //
    if (inradius) {
      cuda_particles[i].r = (float)inradius->GetTuple1(i);
    }
    else {
      cuda_particles[i].r = 0.01;
    }
    //
    if (inintensity) {
      cuda_particles[i].I = (float)inintensity->GetTuple1(i);
    }
    else {
      cuda_particles[i].I = 1;
    }
    cuda_particles[i].type   = 0;
    cuda_particles[i].active = 1;
  }
  // allocate enough space for an array of cu_particle_sim
  cudaMalloc((void **) &newD->userPointer, nPoints*sizeof(cu_particle_sim));
  // copy from host vector to device, first wrap raw pointer with a device_ptr
  thrust::device_ptr<cu_particle_sim> dev_ptr = thrust::device_pointer_cast<cu_particle_sim>((cu_particle_sim*)newD->userPointer);
  cudaMemcpy(newD->userPointer, &cuda_particles[0], nPoints*sizeof(cu_particle_sim), cudaMemcpyHostToDevice);
  //thrust::copy(cuda_particles.begin(), cuda_particles.end(), dev_ptr);
  //
  return 1;
}

