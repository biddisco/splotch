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

vtkStandardNewMacro(vtkSplotchDefaultPainter);
vtkCxxSetObjectMacro(vtkSplotchDefaultPainter, SplotchPainter, vtkSplotchPainter);
//----------------------------------------------------------------------------
vtkSplotchDefaultPainter::vtkSplotchDefaultPainter()
{
  this->SplotchPainter = vtkSplotchPainter::New();
}

//----------------------------------------------------------------------------
vtkSplotchDefaultPainter::~vtkSplotchDefaultPainter()
{
  this->SetSplotchPainter(0);
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::UpdateBounds(double bounds[6])
{
  if (this->SplotchPainter->GetInput()!=this->GetInput()) {
    this->SplotchPainter->SetInput(this->GetInput());
  }
  this->SplotchPainter->UpdateBounds(bounds);
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::BuildPainterChain()
{
  // Override painters we don't want.
  this->SetDisplayListPainter(NULL);
  this->SetCompositePainter(NULL);
  this->SetCoincidentTopologyResolutionPainter(NULL);
  this->SetRepresentationPainter(NULL);
  // and set ours at the end of the chain
  this->SetDefaultPainterDelegate(this->SplotchPainter);
  // allow superclass to pieces everything together
  this->Superclass::BuildPainterChain();
  // We need the ScalarsToColors Painter as the MIP handles scalar mapping specially
  this->SplotchPainter->SetScalarsToColorsPainter(this->GetScalarsToColorsPainter());
}

//----------------------------------------------------------------------------
void vtkSplotchDefaultPainter::ReportReferences(vtkGarbageCollector *collector)
{
  this->Superclass::ReportReferences(collector);
  vtkGarbageCollectorReport(collector, this->SplotchPainter, "SplotchPainter");
}

//----------------------------------------------------------------------------
