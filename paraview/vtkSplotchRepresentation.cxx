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
#include "vtkUnstructuredDataDeliveryFilter.h"
#include "vtkQuadricClustering.h"

vtkStandardNewMacro(vtkSplotchRepresentation);
//----------------------------------------------------------------------------
vtkSplotchRepresentation::vtkSplotchRepresentation()
{
  this->SplotchDefaultPainter    = vtkSplotchDefaultPainter::New();
  this->LODSplotchDefaultPainter = vtkSplotchDefaultPainter::New();
  this->SplotchPainter           = vtkSplotchPainter::New();
  this->LODSplotchPainter        = vtkSplotchPainter::New();
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
  this->Mapper->SetInputConnection(this->UpdateSuppressor->GetOutputPort());
  this->LODMapper->SetInputConnection(this->LODUpdateSuppressor->GetOutputPort());
  // Actors
  this->Actor->SetMapper(this->Mapper);
  this->Actor->SetLODMapper(this->LODMapper);

  // override some settings made in GeometryRepresentation to ensure we get points
  // as output and don't bother copying stuff we don't need.
  this->DeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
  this->LODDeliveryFilter->SetOutputDataType(VTK_POLY_DATA);
  this->Decimator->SetCopyCellData(0);
  // We don't want the MultiBlockMaker as we don't support multiblock
  // connect the GeometryFilter to the CacheKeeper and bypass multiblockmaker.
  // The SplotchDefaultPainter removes the composite painter from the painter chain

  this->CacheKeeper->SetInputConnection(this->GeometryFilter->GetOutputPort());

  // Setup painters
  vtkPainterPolyDataMapper* painterMapper = vtkPainterPolyDataMapper::SafeDownCast(this->Mapper);
  this->SplotchDefaultPainter->SetDelegatePainter(painterMapper->GetPainter()->GetDelegatePainter());
  painterMapper->SetPainter(this->SplotchDefaultPainter);
  this->SplotchDefaultPainter->SetSplotchPainter(this->SplotchPainter);
  // Setup LOD painters
  painterMapper = vtkPainterPolyDataMapper::SafeDownCast(this->LODMapper);
  this->LODSplotchDefaultPainter->SetDelegatePainter(painterMapper->GetPainter()->GetDelegatePainter());
  painterMapper->SetPainter(this->LODSplotchDefaultPainter);
  this->LODSplotchDefaultPainter->SetSplotchPainter(this->LODSplotchPainter);
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
  this->Settings->SetNumberOfTuples(2);
  //
  this->Settings->SetValue(0, NumToStr<int>(this->ActiveParticleType).c_str());
  this->Settings->SetValue(1, NumToStr<int>(this->GetTypeActive()).c_str());
  //
  return this->Settings;
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
void vtkSplotchRepresentation::SetInputArrayToProcess(
  int idx, int port, int connection, int fieldAssociation, const char *name)
{
  switch (idx) {
    case 0: this->SetTypeScalars(name); break;
    case 1: this->SetActiveScalars(name); break;
  }
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetTypeScalars(const char *s)
{
  if (this->SplotchPainter) this->SplotchPainter->SetTypeScalars(s);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetTypeScalars(s);
}
//----------------------------------------------------------------------------
const char *vtkSplotchRepresentation::GetTypeScalars()
{
  if (this->SplotchPainter) return this->SplotchPainter->GetTypeScalars();
  return NULL;
}
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetActiveScalars(const char *s)
{
  if (this->SplotchPainter) this->SplotchPainter->SetActiveScalars(s);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetActiveScalars(s);
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
