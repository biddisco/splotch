/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSplotchPainter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSplotchPainter - vtkSplotchRaytrace.
// .SECTION Description
// .SECTION Implementation
//
// .SECTION See Also

#ifndef __vtkSplotchPainter_h
#define __vtkSplotchPainter_h

#include "vtkPolyDataPainter.h"

#include <vtkstd/vector> // needed for our arrays
#include <vtkstd/string> // needed for our arrays

class vtkMultiProcessController;
class vtkScalarsToColorsPainter;

class VTK_EXPORT vtkSplotchPainter : public vtkPolyDataPainter
{
public:
  static vtkSplotchPainter* New();
  vtkTypeMacro(vtkSplotchPainter, vtkPolyDataPainter);

  // Description:
  void Render(vtkRenderer* renderer, vtkActor* actor, 
    unsigned long typeflags, bool forceCompileOnly);

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

  // Description:
  // The MIP Renderer needs to manually convert scalars to colours
  // so we must have a copy of the painter used by the rest of the rendering pipeline
  virtual void SetScalarsToColorsPainter(vtkScalarsToColorsPainter* ScalarsToColorsPainter);
  vtkGetObjectMacro(ScalarsToColorsPainter, vtkScalarsToColorsPainter);

  vtkSetMacro(ArrayAccessMode, int);
  vtkSetMacro(ArrayId, int);
  vtkSetStringMacro(ArrayName);
  vtkSetMacro(ArrayComponent, int);
  vtkSetMacro(ScalarMode, int);

  // Description:
  // The MIP painter must return the complete bounds of the whole dataset
  // not just the local 'piece', otherwise the compositing blanks out parts it thinks
  // are not covered by any geometry.
  void UpdateBounds(double bounds[6]);

  // Description:
  // Set/Get the controller used for coordinating parallel writing
  // (set to the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);

protected:
   vtkSplotchPainter();
  ~vtkSplotchPainter();

  // Description:
  // Called before RenderInternal() if the Information has been changed
  // since the last time this method was called.
  virtual void ProcessInformation(vtkInformation*);

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

  vtkMultiProcessController *Controller;
  vtkScalarsToColorsPainter *ScalarsToColorsPainter;

  int   ArrayAccessMode;
  int   ArrayComponent;
  int   ArrayId;
  char *ArrayName;
  int   ScalarMode;

private:
  vtkSplotchPainter(const vtkSplotchPainter&); // Not implemented.
  void operator=(const vtkSplotchPainter&); // Not implemented.
};

#endif
