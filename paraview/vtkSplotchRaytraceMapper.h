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

#include "vtkToolkits.h" // For VTK_USE_MPI
#include "vtkPainterPolyDataMapper.h"

class VTK_EXPORT vtkSplotchRaytraceMapper : public vtkPainterPolyDataMapper
{
public:
  static vtkSplotchRaytraceMapper* New();
  vtkTypeMacro(vtkSplotchRaytraceMapper, vtkMapper);

  void Render(vtkRenderer *, vtkActor *);

  vtkSetMacro(Brightness, double);
  vtkGetMacro(Brightness, double);

  vtkBooleanMacro(LogIntensity, int);
  vtkGetMacro(LogIntensity, int);
  vtkSetMacro(LogIntensity, int);

  vtkBooleanMacro(LogColour, int);
  vtkGetMacro(LogColour, int);
  vtkSetMacro(LogColour, int);

  vtkSetMacro(GrayAbsorption, double);
  vtkGetMacro(GrayAbsorption, double);

  vtkSetStringMacro(IntensityScalars);
  vtkGetStringMacro(IntensityScalars);

  vtkSetStringMacro(RadiusScalars);
  vtkGetStringMacro(RadiusScalars);

  vtkSetStringMacro(TypeScalars);
  vtkGetStringMacro(TypeScalars);

  vtkSetStringMacro(ActiveScalars);
  vtkGetStringMacro(ActiveScalars);

  // we need to override the bouds so that IceT composites the whole image 
  // and not only the piece bounds
  void GetBounds(double *bounds);
  double *GetBounds();

protected:
   vtkSplotchRaytraceMapper();
  ~vtkSplotchRaytraceMapper();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  double  Brightness;
  double  GrayAbsorption;
  int     LogIntensity;
  int     LogColour;
  char   *IntensityScalars;
  char   *RadiusScalars;
  char   *TypeScalars;
  char   *ActiveScalars;

private:
  vtkSplotchRaytraceMapper(const vtkSplotchRaytraceMapper&); // Not implemented.
  void operator=(const vtkSplotchRaytraceMapper&); // Not implemented.
};

#endif
