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

#include <thrust/copy.h>
#include <piston/choose_container.h>
#include <piston/choose_container.h>
#include <piston/image3d.h>
#include <piston/vtk_image3d.h>
#include <vector>
//
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPistonDataObject.h"
#include "vtkPistonDataWrangling.h"
#include "vtkPistonReference.h"
#include "vtkCUDAPiston.h"
//
#include "cuda/splotch_cuda.h"
//
/*
  //
  // the structures we need to pass into cuda for splotch use
  //
  struct cu_color {
    float r,g,b;
    };

  struct cu_particle_sim {
      cu_color e;
      float x,y,z,r,I;
      unsigned short type;
      bool active;
    };

typedef struct
{
  //GPU side representation of a vtkPolyData
  //this is the sibling of vtk_image3D in piston
  int nPoints;
  int vertsPer;
  int nCells;
  thrust::device_vector<float3>        *points;
  thrust::device_vector<uint3>         *cells;
  thrust::device_vector<uint3>         *originalcells;
  thrust::device_vector<float>          distances;
  thrust::device_vector<float>         *scalars;
  thrust::device_vector<float>         *opacities;
  thrust::device_vector<float>         *normals;
  thrust::device_vector<uchar4>        *colors;
  thrust::device_vector<uchar4>         lutRGBA;
} vtk_polydata;

*/
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
  //
  //
  //
  // -------------------
  // Following code taken from vtkPistonConverters and tweaked for splotch compatibility
  // -------------------
  vtkPistonReference *tr = od->GetReference();
  tr->type = VTK_POLY_DATA;
  if (!CheckDirty(id, tr)) {
    return 1;
  }
  //
  // clean previous state
  //
  vtkpiston::DeleteData(tr);
  //
  // allocate a new polydata device object
  //
  vtkpiston::vtk_polydata *newD = new vtkpiston::vtk_polydata;
  tr->data = (void*)newD;
  newD->points        = NULL;
  newD->cells         = NULL;
  newD->originalcells = NULL;
  newD->scalars       = NULL;
  newD->opacities     = NULL;
  newD->normals       = NULL;
  newD->colors        = NULL;
  newD->userPointer   = NULL;
  //
  //
  //
  int nPoints = id->GetNumberOfPoints();
  newD->nPoints = nPoints;


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
  thrust::copy(cuda_particles.begin(), cuda_particles.end(), dev_ptr); 
  //
  return 1;
}

