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
#include <vtksys/SystemTools.hxx>
#include "vtkSplotchRepresentation.h"
#include "vtkSplotchDefaultPainter.h"
//
#include "vtksys/ios/sstream"
//
#include "vtkMPI.h"
#include "vtkMPIController.h"
#include "vtkMPICommunicator.h"
//
#include "vtkDataObject.h"
#include "vtkDefaultPainter.h"
#include "vtkSplotchPainter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
// we inherit changes to these filters from GeometryRepresentation
#include "vtkPainterPolyDataMapper.h"
#include "vtkQuadricClustering.h"
#include "vtkPVCacheKeeper.h"
#include "vtkPVUpdateSuppressor.h"
#include "vtkPVLODActor.h"
#include "vtkMultiProcessController.h"
#include "vtkPVRenderView.h"
#include "vtkPVCacheKeeper.h"
#include "vtkPVGeometryFilter.h"

#include <vtksys/SystemInformation.hxx>
#include <vtksys/RegularExpression.hxx>

//#ifdef PARAVIEW_USE_MPI
//#endif

#ifdef PV_SPLOTCH_USE_PISTON
#include "vtkCUDASplotchPainter.h"
#endif

// Otherwise
#include "vtkMultiProcessController.h"
#include "vtkMIPDefaultPainter.h"
#include "vtkMIPPainter.h"

vtkStandardNewMacro(vtkSplotchRepresentation);
//----------------------------------------------------------------------------
vtkSplotchRepresentation::vtkSplotchRepresentation()
{
  this->CheckMPIController();
  //
  this->SplotchDefaultPainter    = vtkSplotchDefaultPainter::New();
  this->LODSplotchDefaultPainter = vtkSplotchDefaultPainter::New();
  this->SplotchPainter           = this->SplotchDefaultPainter->GetSplotchPainter();
  this->LODSplotchPainter        = this->LODSplotchDefaultPainter->GetSplotchPainter();

  this->LODMIPDefaultPainter = vtkMIPDefaultPainter::New();
  this->LODMIPPainter        = this->LODMIPDefaultPainter->GetMIPPainter();

  this->SplotchPainter->Register(this);
  this->LODSplotchPainter->Register(this);
  this->LODMIPPainter->Register(this);
  this->ActiveParticleType   = 0;
//  this->ColorArrayName       = 0;
//  this->ColorAttributeType   = POINT_DATA;
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
void vtkSplotchRepresentation::CheckMPIController()
{
  if (vtkMPIController::GetGlobalController()->IsA("vtkDummyController"))
  {
    vtkDebugMacro("Running vtkDummyController : replacing it");
    int flag = 0;
    MPI_Initialized(&flag);
    if (flag == 0)
    {
      vtkDebugMacro("Running without MPI, attempting to initialize ");
      //int argc = 1;
      //const char *argv = "D:\\cmakebuild\\pv-shared\\bin\\RelWithDebInfo\\paraview.exe";
      //char **_argv = (char**) &argv;
      int provided, rank, size;
      MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      //
      if (rank == 0) {
        if (provided != MPI_THREAD_MULTIPLE) {
          std::cout << "MPI_THREAD_MULTIPLE not set, you may need to recompile your "
          << "MPI distribution with threads enabled" << std::endl;
        }
        else {
          std::cout << "MPI_THREAD_MULTIPLE is OK (splotch(rep) override)" << std::endl;
        }
      }
    }
    //
    vtkDebugMacro("Setting Global MPI controller");

    MPI_Comm *ocomm = new MPI_Comm(MPI_COMM_WORLD);
    vtkMPICommunicatorOpaqueComm comm(ocomm);
    // create a vtkCommunicator and initialize it with the external communicator
    vtkMPICommunicator *communicator = vtkMPICommunicator::New();
    communicator->InitializeExternal(&comm);
    vtkMPICommunicator::SetWorldCommunicator(communicator);

    vtkMPIController *controller = vtkMPIController::New();
    controller->Register(NULL); // to stop it going out of context
//    vtkProcessModule::GlobalController = controller;
    controller->SetCommunicator(communicator);

    if (flag == 0) controller->Initialize(NULL, NULL, 1);
    vtkMPIController::SetGlobalController(controller);
  }
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

  this->MultiBlockMaker->SetInputConnection(this->GeometryFilter->GetOutputPort());
  this->CacheKeeper->SetInputConnection(this->MultiBlockMaker->GetOutputPort());
  this->Decimator->SetInputConnection(this->CacheKeeper->GetOutputPort());

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
//  this->LODSplotchDefaultPainter->SetDelegatePainter(painterMapper->GetPainter()->GetDelegatePainter());
//  painterMapper->SetPainter(this->LODSplotchDefaultPainter);
  this->LODMIPDefaultPainter->SetDelegatePainter(painterMapper->GetPainter()->GetDelegatePainter());
  painterMapper->SetPainter(this->LODMIPDefaultPainter);
  this->Actor->SetEnableLOD(0);
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
void vtkSplotchRepresentation::SetEnableCUDA(int mode)
{
#ifdef PV_SPLOTCH_USE_PISTON
  if (!vtkpiston::IsEnabledCudaGL()) {
    mode = 0;
  }
  this->SplotchDefaultPainter->SetEnableCUDA(mode);
  this->MarkModified();
#endif
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void vtkSplotchRepresentation::SetActiveParticleType(int p)
{
  // this only allocates space in the mapper, it does not actually set the max
  // if (this->SplotchPainter) this->SplotchPainter->SetNumberOfParticleTypes(p+1);
  // if (this->LODSplotchPainter) this->LODSplotchPainter->SetNumberOfParticleTypes(p+1);
  // this is the active one
  this->ActiveParticleType = p;
  // update gui here
  
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
void   vtkSplotchRepresentation::SetMaxRadius(double r)
{
  if (this->SplotchPainter) this->SplotchPainter->SetMaxRadius(this->ActiveParticleType, r);
  if (this->LODSplotchPainter) this->LODSplotchPainter->SetMaxRadius(this->ActiveParticleType, r);
}
//----------------------------------------------------------------------------
double vtkSplotchRepresentation::GetMaxRadius()
{
  return this->SplotchPainter->GetMaxRadius(this->ActiveParticleType);
}
//----------------------------------------------------------------------------
// void vtkSplotchRepresentation::SetLODMIP(int l)
// {
//   if (this->SplotchPainter) this->SplotchPainter->SetLogIntensity(this->ActiveParticleType, l);
//   if (this->LODSplotchPainter) this->LODSplotchPainter->SetLogIntensity(this->ActiveParticleType, l);
// }
//----------------------------------------------------------------------------
// int vtkSplotchRepresentation::GetLODMIP()
// {
//   return this->SplotchPainter->GetLogIntensity(this->ActiveParticleType);
// }
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
/*
void vtkSplotchRepresentation::SetInputArrayToProcess(
  int idx, int port, int connection, int fieldAssociation, const char *name)
{
  std::cout << "SetInputArrayToProcess " << idx << " " << name << std::endl;
  switch (idx) {
    case 0: this->SetIntensityScalars(name); break;
    case 1: this->SetRadiusScalars(name); break;
    case 2: this->SetTypeScalars(name); break;
    case 3: this->SetActiveScalars(name); break;
  }
  vtkGeometryRepresentation::SetInputArrayToProcess(
      idx, port, connection, fieldAssociation, name
  );
}
*/
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
bool vtkSplotchRepresentation::AddToView(vtkView* view)
{
  vtkPVRenderView* rview = vtkPVRenderView::SafeDownCast(view);
#ifdef PV_SPLOTCH_USE_PISTON
  if (rview && !vtkpiston::IsEnabledCudaGL()) {
    
    // if we have no hints as to the display VAR, we'll use the rank (modulus GPU count)
    vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
    int rank = controller->GetLocalProcessId();

    // see if a DISPLAY env var was set for us to get the GPU from
    // match DISPLAY vars of the form :0.0, :0.1, :0.2, :0.N - and get N
    std::string display;
    vtksys::SystemTools::GetEnv("DISPLAY", display);
    vtksys::RegularExpression regex(".*:.*\\.([0-9]+)");
    int displaynum = -1;
    if (regex.find(display.c_str())) {
      try {
        if (regex.match(1).size()>0) {
          displaynum = atoi(regex.match(1).c_str());
        }
      }
      catch (...)
      {
      }
    }
    else {
#ifndef _WIN32
      vtkWarningMacro("DISPLAY environment variable should conform to \":0.0\" format");
#endif
    }
    //
    bool ok = vtkpiston::InitCudaGL(rview->GetRenderWindow(), rank, displaynum);
    // 
    vtksys::SystemInformation sysInfo;
    std::string hostname = sysInfo.GetHostname();
    if (ok) {
      std::cout <<"Rank " << std::setw(3) <<  rank << " Hostname " << hostname.c_str() << " Initialized CudaGL with GPU " << displaynum << std::endl;
    }
    else {
      std::cout <<"Rank " << std::setw(3) <<  rank << " Hostname " << hostname.c_str() << " CudaGL Failed " << std::endl;
    }
  }
#endif
  return this->Superclass::AddToView(view);
}

//----------------------------------------------------------------------------
int vtkSplotchRepresentation::ProcessViewRequest(
  vtkInformationRequestKey* request_type,
  vtkInformation* inInfo, vtkInformation* outInfo)
{
  if (this->GetVisibility() == false)
    {
    return 0;
    }

  if (request_type == vtkPVView::REQUEST_UPDATE_LOD())
    {
    // Called to generate and provide the LOD data to the view.
    // If SuppressLOD is true, we tell the view we have no LOD data to provide,
    // otherwise we provide the decimated data.
    if (!this->SuppressLOD)
      {
      if (inInfo->Has(vtkPVRenderView::USE_OUTLINE_FOR_LOD()))
        {
        // HACK to ensure that when Decimator is next employed, it delivers a
        // new geometry.
        this->Decimator->Modified();

        this->LODOutlineFilter->Update();
        // Pass along the LOD geometry to the view so that it can deliver it to
        // the rendering node as and when needed.
        vtkPVRenderView::SetPieceLOD(inInfo, this,
          this->LODOutlineFilter->GetOutputDataObject(0));
        }
      else
        {
        // HACK to ensure that when Decimator is next employed, it delivers a
        // new geometry.
/*
        this->LODOutlineFilter->Modified();

        if (inInfo->Has(vtkPVRenderView::LOD_RESOLUTION()))
          {
          int division = static_cast<int>(150 *
            inInfo->Get(vtkPVRenderView::LOD_RESOLUTION())) + 10;
          division = 10; 
          this->Decimator->SetNumberOfDivisions(division, division, division);
          }

        this->Decimator->Update();
*/
        // Pass along the LOD geometry to the view so that it can deliver it to
        // the rendering node as and when needed.
        vtkPVRenderView::SetPieceLOD(inInfo, this,
          this->CacheKeeper->GetOutputDataObject(0));
        }
      }
    }
  else if (!this->Superclass::ProcessViewRequest(request_type, inInfo, outInfo))
    {
    // i.e. this->GetVisibility() == false, hence nothing to do.
    return 0;
    }

  return 1;
}


//----------------------------------------------------------------------------
