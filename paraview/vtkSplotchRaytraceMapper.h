/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSplotchRaytraceMapper.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSplotchRaytraceMapper - vtkSplotchRaytrace.
// .SECTION Description
// .SECTION Implementation
//
// .SECTION See Also

#ifndef __vtkSplotchRaytraceMapper_h
#define __vtkSplotchRaytraceMapper_h

#include "vtkPainterPolyDataMapper.h"

class VTK_EXPORT vtkSplotchRaytraceMapper : public vtkPainterPolyDataMapper
{
public:
  static vtkSplotchRaytraceMapper* New();
  vtkTypeMacro(vtkSplotchRaytraceMapper, vtkMapper);

  void Render(vtkRenderer *, vtkActor *);

  vtkSetStringMacro(IntensityScalars);
  vtkGetStringMacro(IntensityScalars);

  vtkSetStringMacro(RadiusScalars);
  vtkGetStringMacro(RadiusScalars);

protected:
   vtkSplotchRaytraceMapper();
  ~vtkSplotchRaytraceMapper();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  char *IntensityScalars;
  char *RadiusScalars;

private:
  vtkSplotchRaytraceMapper(const vtkSplotchRaytraceMapper&); // Not implemented.
  void operator=(const vtkSplotchRaytraceMapper&); // Not implemented.
};

#endif
