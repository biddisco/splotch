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
#include "vtkStringArray.h"
#include "vtkSmartPointer.h"

class vtkSplotchRaytraceMapper;

class VTK_EXPORT vtkSplotchRaytraceRepresentation : public vtkGeometryRepresentation
{
public:
  static vtkSplotchRaytraceRepresentation* New();
  vtkTypeMacro(vtkSplotchRaytraceRepresentation, vtkGeometryRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // To simplify the GUI interaction, we make one particle type active
  // so that SetBrightness, SetActiveScalars etc act on that type.
  void SetActiveParticleType(int p);
  vtkGetMacro(ActiveParticleType, int);

  //**************************************************************************
  // Forwarded to vtkSplotchRaytraceMapper
  //**************************************************************************
  virtual void SetInputArrayToProcess(int idx, int port, int connection,
                              int fieldAssociation,
                              const char *name);

  void SetIntensityScalars(const char *);
  void SetRadiusScalars(const char *);
  void SetTypeScalars(const char *);
  void SetActiveScalars(const char *);
  const char *GetIntensityScalars();
  const char *GetRadiusScalars();
  const char *GetTypeScalars();
  const char *GetActiveScalars();

  void   SetBrightness(double b);
  double GetBrightness();

  void   SetLogIntensity(int l);
  int    GetLogIntensity();

  void   SetGrayAbsorption(double b);
  double GetGrayAbsorption();

  // Gather all the settings in one call for feeding back to the gui display
  vtkStringArray *GetActiveParticleSettings();

//BTX
protected:
  vtkSplotchRaytraceRepresentation();
  ~vtkSplotchRaytraceRepresentation();

  // Description:
  // Fill input port information.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Description:
  // Execute the pipeline, not used, but can instantiate filters for extra processing
  // in here.
  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  //
  vtkSplotchRaytraceMapper       *SplotchMapper;
  vtkSplotchRaytraceMapper       *LODSplotchMapper;
  int                             ActiveParticleType;
  vtkSmartPointer<vtkStringArray> Settings;

private:
  vtkSplotchRaytraceRepresentation(const vtkSplotchRaytraceRepresentation&); // Not implemented
  void operator=(const vtkSplotchRaytraceRepresentation&); // Not implemented
//ETX
};

#endif
