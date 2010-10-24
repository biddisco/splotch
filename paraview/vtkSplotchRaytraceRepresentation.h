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
// .NAME vtkSplotchRaytraceRepresentation
// .SECTION Description
// vtkSplotchRaytraceRepresentation is a representation that uses the vtkSplotchRaytraceMapper
// for rendering glyphs.

#ifndef __vtkSplotchRaytraceRepresentation_h
#define __vtkSplotchRaytraceRepresentation_h

#include "vtkGeometryRepresentation.h"

class vtkSplotchRaytraceMapper;

class VTK_EXPORT vtkSplotchRaytraceRepresentation : public vtkGeometryRepresentation
{
public:
  static vtkSplotchRaytraceRepresentation* New();
  vtkTypeMacro(vtkSplotchRaytraceRepresentation, vtkGeometryRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get/Set the visibility for this representation. When the visibility of
  // representation of false, all view passes are ignored.
  virtual void SetVisibility(bool);

  //**************************************************************************
  // Forwarded to vtkSplotchRaytraceMapper

//BTX
protected:
  vtkSplotchRaytraceRepresentation();
  ~vtkSplotchRaytraceRepresentation();

  // Description:
  // Fill input port information.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  virtual int RequestData(vtkInformation*,
    vtkInformationVector**, vtkInformationVector*);

  vtkSplotchRaytraceMapper              *SplotchMapper;
  vtkSplotchRaytraceMapper              *LODSplotchMapper;

private:
  vtkSplotchRaytraceRepresentation(const vtkSplotchRaytraceRepresentation&); // Not implemented
  void operator=(const vtkSplotchRaytraceRepresentation&); // Not implemented
//ETX
};

#endif
