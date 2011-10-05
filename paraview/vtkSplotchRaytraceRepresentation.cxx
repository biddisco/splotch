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

#include "vtksys/ios/sstream"

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

//  this->GrayAbsorption = 0.0001;
//  this->Brightness = 10.5;
//  this->LogIntensity = 1;
  this->ActiveParticleType = 0;

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

  Settings = vtkSmartPointer<vtkStringArray>::New();

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
void vtkSplotchRaytraceRepresentation::SetActiveParticleType(int p)
{
  // this only allocates space in the mapper, it does not actually set the max
  if (this->SplotchMapper) this->SplotchMapper->SetNumberOfParticleTypes(p+1);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetNumberOfParticleTypes(p+1);
  // this is the active one
  this->ActiveParticleType = p;
}
//----------------------------------------------------------------------------
template <typename T>
std::string NumToStr(T data) {
  vtksys_ios::ostringstream oss;
  oss.setf(0,ios::floatfield);
  oss.precision(5);  
  oss << data;
  return oss.str();
}
//----------------------------------------------------------------------------
vtkStringArray *vtkSplotchRaytraceRepresentation::GetActiveParticleSettings()
{
  this->Settings->Initialize();
  this->Settings->SetNumberOfComponents(1);
  this->Settings->SetNumberOfTuples(5);

  this->Settings->SetValue(0, NumToStr<int>(this->ActiveParticleType).c_str());
  this->Settings->SetValue(1, NumToStr<double>(this->GetBrightness()).c_str());
  this->Settings->SetValue(2, NumToStr<int>(this->GetLogIntensity()).c_str());
  this->Settings->SetValue(3, this->GetIntensityScalars());
  this->Settings->SetValue(4, this->GetRadiusScalars());
  //

  return this->Settings;
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetBrightness(double b)
{
  double value = pow(10,(b/100.0));
  if (this->SplotchMapper) this->SplotchMapper->SetBrightness(this->ActiveParticleType, value);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetBrightness(this->ActiveParticleType, value);
}
//----------------------------------------------------------------------------
double vtkSplotchRaytraceRepresentation::GetBrightness()
{
  return this->SplotchMapper->GetBrightness(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetLogIntensity(int l)
{
  if (this->SplotchMapper) this->SplotchMapper->SetLogIntensity(this->ActiveParticleType, l);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetLogIntensity(this->ActiveParticleType, l);
}
//----------------------------------------------------------------------------
int vtkSplotchRaytraceRepresentation::GetLogIntensity()
{
  return this->SplotchMapper->GetLogIntensity(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetGrayAbsorption(double g)
{
  double value = pow(10,(g/100.0));
  if (this->SplotchMapper) this->SplotchMapper->SetGrayAbsorption(value);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetGrayAbsorption(value);
}
//----------------------------------------------------------------------------
double vtkSplotchRaytraceRepresentation::GetGrayAbsorption()
{
  return this->SplotchMapper->GetGrayAbsorption();
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetInputArrayToProcess(
  int idx, int port, int connection, int fieldAssociation, const char *name)
{
  switch (idx) {
    case 0: this->SetIntensityScalars(name); break;
    case 1: this->SetRadiusScalars(name); break;
    case 2: this->SetTypeScalars(name); break;
    case 3: this->SetActiveScalars(name); break;
  }
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetIntensityScalars(const char *s)
{
  if (this->SplotchMapper) this->SplotchMapper->SetIntensityScalars(this->ActiveParticleType, s);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetIntensityScalars(this->ActiveParticleType, s);
}
//----------------------------------------------------------------------------
void vtkSplotchRaytraceRepresentation::SetRadiusScalars(const char *s)
{
  if (this->SplotchMapper) this->SplotchMapper->SetRadiusScalars(this->ActiveParticleType, s);
  if (this->LODSplotchMapper) this->LODSplotchMapper->SetRadiusScalars(this->ActiveParticleType, s);
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
const char *vtkSplotchRaytraceRepresentation::GetIntensityScalars()
{
  if (this->SplotchMapper) return this->SplotchMapper->GetIntensityScalars(this->ActiveParticleType);
  return NULL;
}
//----------------------------------------------------------------------------
const char *vtkSplotchRaytraceRepresentation::GetRadiusScalars()
{
  if (this->SplotchMapper) return this->SplotchMapper->GetRadiusScalars(this->ActiveParticleType);
  return NULL;
}
//----------------------------------------------------------------------------
const char *vtkSplotchRaytraceRepresentation::GetTypeScalars()
{
  if (this->SplotchMapper) return this->SplotchMapper->GetTypeScalars();
  return NULL;
}
//----------------------------------------------------------------------------
const char *vtkSplotchRaytraceRepresentation::GetActiveScalars()
{
  if (this->SplotchMapper) return this->SplotchMapper->GetActiveScalars();
  return NULL;
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
