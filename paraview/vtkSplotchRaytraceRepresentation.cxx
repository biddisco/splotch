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
  this->SplotchMapper    = vtkSplotchRaytraceMapper::New();
  this->Mapper           = this->SplotchMapper;
  this->LODSplotchMapper = NULL;
  this->LODMapper        = vtkPolyDataMapper::New();

  this->Mapper->SetInputConnection(this->Distributor->GetOutputPort());
  this->LODMapper->SetInputConnection(this->LODDeliveryFilter->GetOutputPort());

  this->Actor->SetMapper(this->Mapper);
  this->Actor->SetLODMapper(this->LODMapper);
  this->Actor->SetProperty(this->Property);

  // override some settings made in GeometryRepresentation
  this->DeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
  this->LODDeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
  this->Decimator->SetCopyCellData(0);
  //  we don't want the MultiBlockMaker used
  this->CacheKeeper->SetInputConnection(this->GeometryFilter->GetOutputPort());

  this->ColorArrayName = 0;
  this->ColorAttributeType = POINT_DATA;
  this->Representation = POINTS;

  // Not insanely thrilled about this API on vtkProp about properties, but oh
  // well. We have to live with it.
  vtkInformation* keys = vtkInformation::New();
  this->Actor->SetPropertyKeys(keys);
  keys->Delete();

}
//----------------------------------------------------------------------------
vtkSplotchRaytraceRepresentation::~vtkSplotchRaytraceRepresentation()
{
  // Geometry Representation base class will delete the Mapper and LODMapper which point to our classes
  this->SplotchMapper = NULL;
  this->LODSplotchMapper = NULL;
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
void vtkSplotchRaytraceRepresentation::SetTypeScalars(const char *s)
{
  if (this->SplotchMapper) this->SplotchMapper->SetTypeScalars(s);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetTypeScalars(s);
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetActiveScalars(const char *s)
{
  if (this->SplotchMapper) this->SplotchMapper->SetActiveScalars(s);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetActiveScalars(s);
}
//----------------------------------------------------------------------------



//----------------------------------------------------------------------------
