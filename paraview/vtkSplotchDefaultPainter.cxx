/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSplotchDefaultPainter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSplotchDefaultPainter.h"

#include "vtkGarbageCollector.h"
#include "vtkSplotchPainter.h"
#include "vtkObjectFactory.h"
#include "vtkScalarsToColorsPainter.h"
#include "vtkInformation.h"
#include "vtkDataSet.h"

// If the MIP renderer is available
#ifdef PV_SPLOTCH_WITH_MIP
 #include "vtkMIPPainter.h"
#endif
#ifdef PV_SPLOTCH_USE_PISTON
  #include "vtkCUDASplotchPainter.h"
#endif

vtkStandardNewMacro(vtkSplotchDefaultPainter);
vtkCxxSetObjectMacro(vtkSplotchDefaultPainter, SplotchPainter, vtkSplotchPainter);
//----------------------------------------------------------------------------
vtkSplotchDefaultPainter::vtkSplotchDefaultPainter()
{
#ifndef PV_SPLOTCH_USE_PISTON
  this->SplotchPainter = vtkSplotchPainter::New();
#else
  this->SplotchPainter = vtkCUDASplotchPainter::New();
#endif
  this->MIPPainter = NULL;
  this->UseMIP = 0;
#ifdef PV_SPLOTCH_WITH_MIP
  this->MIPPainter = vtkMIPPainter::New();
#endif
}

//----------------------------------------------------------------------------
vtkSplotchDefaultPainter::~vtkSplotchDefaultPainter()
{
  this->SetSplotchPainter(0);
  if (this->MIPPainter) this->MIPPainter->Delete();
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::SetUseMIP(int m)
{
#ifdef PV_SPLOTCH_WITH_MIP
  if (!this->MIPPainter) {
    this->MIPPainter = vtkMIPPainter::New();
  }
#endif
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::SetEnableCUDA(int m)
{
  this->EnableCUDA = m;
  this->SplotchPainter->SetEnableCUDA(m);
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::UpdateBounds(double bounds[6])
{
  this->SplotchPainter->UpdateBounds(bounds, vtkDataSet::SafeDownCast(this->GetInput()));
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::BuildPainterChain()
{
  // Override painters we don't want.
  this->SetClipPlanesPainter(NULL);
  this->SetDisplayListPainter(NULL);
  this->SetCompositePainter(NULL);
  // Lighting painter aborts render if no input, which locks up our collectives
  this->SetLightingPainter(NULL);
  this->SetRepresentationPainter(NULL);
  this->SetCoincidentTopologyResolutionPainter(NULL);
  //
  // allow superclass to set its internals (nothing left except scalars to colours)
  this->Superclass::BuildPainterChain();
  //
  vtkScalarsToColorsPainter *scalarsToColors = this->GetScalarsToColorsPainter();
  // 
  scalarsToColors->SetDelegatePainter(this->SplotchPainter);


//  vtkPainter *next = scalarsToColors->GetDelegatePainter();
//  scalarsToColors->SetDelegatePainter(NULL);

  // We need the ScalarsToColors Painter as the MIP handles scalar mapping specially
//  this->SplotchPainter->SetDelegatePainter(scalarsToColors);
  this->SplotchPainter->SetScalarsToColorsPainter(scalarsToColors); 
}

//-----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::ProcessInformation(vtkInformation* info)
{
  info->Set(vtkScalarsToColorsPainter::INTERPOLATE_SCALARS_BEFORE_MAPPING(), 0);
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::ReportReferences(vtkGarbageCollector *collector)
{
  this->Superclass::ReportReferences(collector);
  vtkGarbageCollectorReport(collector, this->SplotchPainter, "SplotchPainter");
}

//----------------------------------------------------------------------------
