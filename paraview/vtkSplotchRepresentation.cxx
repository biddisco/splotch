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
#include "vtkSplotchRepresentation.h"
#include "vtkSplotchDefaultPainter.h"
//
#include "vtksys/ios/sstream"
//
#include "vtkDataObject.h"
#include "vtkDefaultPainter.h"
#include "vtkSplotchPainter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
// we inherit changes to these filters from GeometryRepresentation
#include "vtkPainterPolyDataMapper.h"
#include "vtkPVCacheKeeper.h"
#include "vtkPVUpdateSuppressor.h"
#include "vtkPVLODActor.h"
#include "vtkQuadricClustering.h"

vtkStandardNewMacro(vtkSplotchRepresentation);
//----------------------------------------------------------------------------
vtkSplotchRepresentation::vtkSplotchRepresentation()
{
  this->SplotchDefaultPainter    = vtkSplotchDefaultPainter::New();
  this->LODSplotchDefaultPainter = vtkSplotchDefaultPainter::New();
  this->SplotchPainter           = this->SplotchDefaultPainter->GetSplotchPainter();
  this->LODSplotchPainter        = this->LODSplotchDefaultPainter->GetSplotchPainter();
  this->SplotchPainter->Register(this);
  this->LODSplotchPainter->Register(this);
  this->ActiveParticleType   = 0;
  this->ColorArrayName       = 0;
  this->ColorAttributeType   = POINT_DATA;
  this->Representation       = POINTS;
  this->Settings             = vtkSmartPointer<vtkStringArray>::New();
  //
  // The default Painter based Mapper : vtkCompositePolyDataMapper2 does not
  // pass the ComputeBounds through to the individual painters, so our screenspace
  // compositing from IceT is not handled well. 
  // Since we can't handle multiblock data anyway, use a PolyDataPainter mapper
  //
  this->Mapper->Delete();
  this->LODMapper->Delete();
  this->Mapper = vtkPainterPolyDataMapper::New();
  this->LODMapper = vtkPainterPolyDataMapper::New();
  //
  this->SetupDefaults();
}
//----------------------------------------------------------------------------
vtkSplotchRepresentation::~vtkSplotchRepresentation()
{
  this->SplotchDefaultPainter->Delete();
  this->LODSplotchDefaultPainter->Delete();
  this->SplotchPainter->Delete();
  this->LODSplotchPainter->Delete();
}

//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetupDefaults()
{
  // we changed the default Mapper so we must modify the connections affected
//  this->Mapper->SetInputConnection(this->UpdateSuppressor->GetOutputPort());
//  this->LODMapper->SetInputConnection(this->LODUpdateSuppressor->GetOutputPort());
  // Actors
  this->Actor->SetMapper(this->Mapper);
  this->Actor->SetLODMapper(this->LODMapper);

  // override some settings made in GeometryRepresentation to ensure we get points
  // as output and don't bother copying stuff we don't need.
//  this->DeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
//  this->LODDeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
  this->Decimator->SetCopyCellData(0);
  // We don't want the MultiBlockMaker as we don't support multiblock
  // connect the GeometryFilter to the CacheKeeper and bypass multiblockmaker.
  // The SplotchDefaultPainter removes the composite painter from the painter chain

  this->CacheKeeper->SetInputConnection(this->GeometryFilter->GetOutputPort());

  // Setup painters
  vtkPainterPolyDataMapper* painterMapper = vtkPainterPolyDataMapper::SafeDownCast(this->Mapper);
  this->SplotchDefaultPainter->SetDelegatePainter(painterMapper->GetPainter()->GetDelegatePainter());
  painterMapper->SetPainter(this->SplotchDefaultPainter);
  painterMapper->SetInterpolateScalarsBeforeMapping(0);

  // Setup LOD painters
  painterMapper = vtkPainterPolyDataMapper::SafeDownCast(this->LODMapper);
  this->LODSplotchDefaultPainter->SetDelegatePainter(painterMapper->GetPainter()->GetDelegatePainter());
  painterMapper->SetPainter(this->LODSplotchDefaultPainter);
}

//----------------------------------------------------------------------------
int vtkSplotchRepresentation::FillInputPortInformation(int port,
  vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkSplotchRepresentation::RequestData(vtkInformation* request,
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  return this->Superclass::RequestData(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
void vtkSplotchRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetActiveParticleType(int p)
{
  // this only allocates space in the mapper, it does not actually set the max
  if (this->SplotchPainter) this->SplotchPainter->SetNumberOfParticleTypes(p+1);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetNumberOfParticleTypes(p+1);
  // this is the active one
  this->ActiveParticleType = p;
}
//----------------------------------------------------------------------------
template <typename T>
std::string NumToStr(T data) {
  vtksys_ios::ostringstream oss;
//  oss.setf(0,ios::floatfield);
  oss.precision(5);  
  oss << data;
  return oss.str();
}
//----------------------------------------------------------------------------
vtkStringArray *vtkSplotchRepresentation::GetActiveParticleSettings()
{
  this->Settings->Initialize();
  this->Settings->SetNumberOfComponents(1);
  this->Settings->SetNumberOfTuples(6);

  this->Settings->SetValue(0, NumToStr<int>(this->ActiveParticleType).c_str());
  this->Settings->SetValue(1, NumToStr<double>(this->GetBrightness()).c_str());
  this->Settings->SetValue(2, NumToStr<int>(this->GetLogIntensity()).c_str());
  this->Settings->SetValue(3, this->GetIntensityScalars());
  this->Settings->SetValue(4, this->GetRadiusScalars());
  this->Settings->SetValue(5, NumToStr<int>(this->GetTypeActive()).c_str());
  //
  return this->Settings;
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetBrightness(double b)
{
  double value = pow(10,(b/100.0));
  if (this->SplotchPainter) this->SplotchPainter->SetBrightness(this->ActiveParticleType, value);
//  if (this->LODSplotchPainter) this->LODSplotchPainter->SetBrightness(this->ActiveParticleType, value);
}
//----------------------------------------------------------------------------
double vtkSplotchRepresentation::GetBrightness()
{
  return this->SplotchPainter->GetBrightness(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetBrightnessLOD(double b)
{
  double value = pow(10,(b/100.0));
//  if (this->SplotchPainter) this->SplotchPainter->SetBrightnessLOD(this->ActiveParticleType, value);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetBrightness(this->ActiveParticleType, value);
}
//----------------------------------------------------------------------------
double vtkSplotchRepresentation::GetBrightnessLOD()
{
  return this->LODSplotchPainter->GetBrightness(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetRadiusMultiplier(double r)
{
  double value = pow(10,(r/50.0));
  if (this->SplotchPainter) this->SplotchPainter->SetRadiusMultiplier(value);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetRadiusMultiplier(value);
}
//----------------------------------------------------------------------------
double vtkSplotchRepresentation::GetRadiusMultiplier()
{
  return this->SplotchPainter->GetRadiusMultiplier();
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetLODMIP(int l)
{
  if (this->SplotchPainter) this->SplotchPainter->SetLogIntensity(this->ActiveParticleType, l);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetLogIntensity(this->ActiveParticleType, l);
}
//----------------------------------------------------------------------------
int vtkSplotchRepresentation::GetLODMIP()
{
  return this->SplotchPainter->GetLogIntensity(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetLogIntensity(int l)
{
  if (this->SplotchPainter) this->SplotchPainter->SetLogIntensity(this->ActiveParticleType, l);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetLogIntensity(this->ActiveParticleType, l);
}
//----------------------------------------------------------------------------
int vtkSplotchRepresentation::GetLogIntensity()
{
  return this->SplotchPainter->GetLogIntensity(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetTypeActive(int l)
{
  if (this->SplotchPainter) this->SplotchPainter->SetTypeActive(this->ActiveParticleType, l);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetTypeActive(this->ActiveParticleType, l);
}
//----------------------------------------------------------------------------
int vtkSplotchRepresentation::GetTypeActive()
{
  return this->SplotchPainter->GetTypeActive(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetGrayAbsorption(double g)
{
  double value = pow(10,(g/100.0));
  if (this->SplotchPainter) this->SplotchPainter->SetGrayAbsorption(value);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetGrayAbsorption(value);
}
//----------------------------------------------------------------------------
double vtkSplotchRepresentation::GetGrayAbsorption()
{
  return this->SplotchPainter->GetGrayAbsorption();
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetInputArrayToProcess(
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
void vtkSplotchRepresentation::SetIntensityScalars(const char *s)
{
  if (this->SplotchPainter) this->SplotchPainter->SetIntensityScalars(this->ActiveParticleType, s);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetIntensityScalars(this->ActiveParticleType, s);
  }
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetRadiusScalars(const char *s)
{
  if (this->SplotchPainter) this->SplotchPainter->SetRadiusScalars(this->ActiveParticleType, s);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetRadiusScalars(this->ActiveParticleType, s);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetTypeScalars(const char *s)
{
  if (this->SplotchPainter) this->SplotchPainter->SetTypeScalars(s);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetTypeScalars(s);
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetActiveScalars(const char *s)
{
  if (this->SplotchPainter) this->SplotchPainter->SetActiveScalars(s);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetActiveScalars(s);
}
//----------------------------------------------------------------------------
const char *vtkSplotchRepresentation::GetIntensityScalars()
{
  if (this->SplotchPainter) return this->SplotchPainter->GetIntensityScalars(this->ActiveParticleType);
  return NULL;
}
//----------------------------------------------------------------------------
const char *vtkSplotchRepresentation::GetRadiusScalars()
{
  if (this->SplotchPainter) return this->SplotchPainter->GetRadiusScalars(this->ActiveParticleType);
  return NULL;
}
//----------------------------------------------------------------------------
const char *vtkSplotchRepresentation::GetTypeScalars()
{
  if (this->SplotchPainter) return this->SplotchPainter->GetTypeScalars();
  return NULL;
}
//----------------------------------------------------------------------------
const char *vtkSplotchRepresentation::GetActiveScalars()
{
  if (this->SplotchPainter) return this->SplotchPainter->GetActiveScalars();
  return NULL;
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
