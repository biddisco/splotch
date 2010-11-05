/*=========================================================================

  Program:   ParaView
  Module:    $RCSfile$

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSplotchRaytraceRepresentation.h"

#include "vtkCompositePolyDataMapper2.h"
#include "vtkDataObject.h"
#include "vtkSplotchRaytraceMapper.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMPIMoveData.h"
#include "vtkObjectFactory.h"
#include "vtkPVArrowSource.h"
#include "vtkPVLODActor.h"
#include "vtkPVRenderView.h"
#include "vtkQuadricClustering.h"
#include "vtkRenderer.h"
#include "vtkUnstructuredDataDeliveryFilter.h"
#include "vtkOrderedCompositeDistributor.h"
#include "vtkPVCacheKeeper.h"

vtkStandardNewMacro(vtkSplotchRaytraceRepresentation);
//----------------------------------------------------------------------------
vtkSplotchRaytraceRepresentation::vtkSplotchRaytraceRepresentation()
{
  this->SplotchMapper = vtkSplotchRaytraceMapper::New();
  this->LODSplotchMapper = vtkSplotchRaytraceMapper::New();

  this->SplotchMapper->SetInputConnection(0,this->Mapper->GetInputConnection(0, 0));
  this->LODSplotchMapper->SetInputConnection(0,this->LODMapper->GetInputConnection(0, 0));
  // does the same as above, just checking
  this->Mapper->SetInputConnection(this->Distributor->GetOutputPort());
  this->LODMapper->SetInputConnection(this->LODDeliveryFilter->GetOutputPort());

  this->Actor->SetMapper(this->SplotchMapper);
  this->Actor->SetLODMapper(this->LODSplotchMapper);
  this->Actor->SetProperty(this->Property);

  // override some settings made in GeometryRepresentation
  this->DeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
  this->LODDeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
  this->Decimator->SetCopyCellData(0);
  //  we don't want the MultiBlockMaker used
  this->CacheKeeper->SetInputConnection(this->GeometryFilter->GetOutputPort());

  this->ColorArrayName = 0;
  this->ColorAttributeType = POINT_DATA;
  this->Representation = SURFACE;

  // Not insanely thrilled about this API on vtkProp about properties, but oh
  // well. We have to live with it.
  vtkInformation* keys = vtkInformation::New();
  this->Actor->SetPropertyKeys(keys);
  keys->Delete();

}
//----------------------------------------------------------------------------
vtkSplotchRaytraceRepresentation::~vtkSplotchRaytraceRepresentation()
{
  this->SplotchMapper->Delete();
  this->LODSplotchMapper->Delete();
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetVisibility(bool val)
{
  this->Superclass::SetVisibility(val);
  this->Actor->SetVisibility(val);
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetIntensityScalars(const char *s)
{
  if (this->SplotchMapper) this->SplotchMapper->SetIntensityScalars(s);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetIntensityScalars(s);
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetRadiusScalars(const char *s)
{
  if (this->SplotchMapper) this->SplotchMapper->SetRadiusScalars(s);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetRadiusScalars(s);
}
//----------------------------------------------------------------------------
int vtkSplotchRaytraceRepresentation::FillInputPortInformation(int port,
  vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkSplotchRaytraceRepresentation::RequestData(vtkInformation* request,
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  return this->Superclass::RequestData(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//**************************************************************************
// Forwarded to vtkPointSpriteMapper
//----------------------------------------------------------------------------
