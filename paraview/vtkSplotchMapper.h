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

#include <vtkstd/vector> // needed for our arrays
#include <vtkstd/string> // needed for our arrays

class VTK_EXPORT vtkSplotchRaytraceMapper : public vtkPainterPolyDataMapper
{
public:
  static vtkSplotchRaytraceMapper* New();
  vtkTypeMacro(vtkSplotchRaytraceMapper, vtkMapper);

  void Render(vtkRenderer *, vtkActor *);

  // Description:
  // Each particle may have a type, this must be an integer
  // usuall, dark matter, gas, star, etc are denoted by their type Id
  // one array maps all particles
  vtkSetStringMacro(TypeScalars);
  vtkGetStringMacro(TypeScalars);

  // Description:
  // Each particle may be visible or invisible, the active array is used to pass
  // a yes/no value. usually an int array (or char) is used.
  // one array maps all particles
  vtkSetStringMacro(ActiveScalars);
  vtkGetStringMacro(ActiveScalars);

  // Description:
  // There may be N particle types. Each type has its own colour table,
  // Intensity array, brigthness etc. The remaining values are defined per
  // particle type.
  void SetNumberOfParticleTypes(int N);
  vtkGetMacro(NumberOfParticleTypes, int);

  void SetIntensityScalars(int ptype, const char *s);
  const char *GetIntensityScalars(int ptype);

  void SetRadiusScalars(int ptype, const char *s);
  const char *GetRadiusScalars(int ptype);

  void SetBrightness(int ptype, double);
  double GetBrightness(int ptype);

  void SetLogIntensity(int ptype, int);
  int GetLogIntensity(int ptype);

  // don't need this?
  void SetLogColour(int ptype, int);
  int GetLogColour(int ptype);

  void SetTypeActive(int ptype, int);
  int GetTypeActive(int ptype);

  vtkGetMacro(GrayAbsorption,double);
  vtkSetMacro(GrayAbsorption,double);

  // we need to override the bounds so that IceT composites the whole image 
  // and not only the projected piece bounds from each process
  void GetBounds(double *bounds);
  double *GetBounds();

protected:
   vtkSplotchRaytraceMapper();
  ~vtkSplotchRaytraceMapper();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  char   *TypeScalars;
  char   *ActiveScalars;
  double  GrayAbsorption;
  int     NumberOfParticleTypes;
  std::vector<std::string> IntensityScalars;
  std::vector<std::string> RadiusScalars;
  std::vector<double>      Brightness;
  std::vector<int>         LogIntensity;
  std::vector<int>         LogColour;
  std::vector<int>         TypeActive;

private:
  vtkSplotchRaytraceMapper(const vtkSplotchRaytraceMapper&); // Not implemented.
  void operator=(const vtkSplotchRaytraceMapper&); // Not implemented.
};

#endif
