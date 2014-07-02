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
// .NAME vtkSplotchRepresentation
// .SECTION Description
// vtkSplotchRepresentation is a representation that uses the splotch
// for rendering particles.

#ifndef __vtkSplotchRepresentation_h
#define __vtkSplotchRepresentation_h

#include "pv_splotch_configure.h"
#include "vtkGeometryRepresentation.h"
#include "vtkStringArray.h"
#include "vtkSmartPointer.h"

class vtkSplotchPainter;
class vtkSplotchDefaultPainter;

class pv_splotch_EXPORT vtkSplotchRepresentation : public vtkGeometryRepresentation
{
public:
  static vtkSplotchRepresentation* New();
  vtkTypeMacro(vtkSplotchRepresentation, vtkGeometryRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // To simplify the GUI interaction, we make one particle type active
  // so that SetBrightness, SetActiveScalars etc act on that type.
  void SetActiveParticleType(int p);
  vtkGetMacro(ActiveParticleType, int);

  //**************************************************************************
  // Forwarded to vtkSplotchMapper
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

  void   SetBrightnessLOD(double b);
  double GetBrightnessLOD();

  void   SetRadiusMultiplier(double r);
  double GetRadiusMultiplier();

  void   SetMaxRadius(double r);
  double GetMaxRadius();

  // void   SetLODMIP(int l);
  // int    GetLODMIP();

  void   SetLogIntensity(int l);
  int    GetLogIntensity();

  void   SetTypeActive(int l);
  int    GetTypeActive();

  void   SetGrayAbsorption(double b);
  double GetGrayAbsorption();

  void   SetEnableCUDA(int mode);

  // Gather all the settings in one call for feeding back to the gui display
  vtkStringArray *GetActiveParticleSettings();

  bool AddToView(vtkView* view);

//BTX
protected:
  vtkSplotchRepresentation();
  ~vtkSplotchRepresentation();

  void CheckMPIController();
  
  // Description:
  // This method is called in the constructor. If the subclasses override any of
  // the iVar vtkObject's of this class e.g. the Mappers, GeometryFilter etc.,
  // they should call this method again in their constructor. It must be totally
  // safe to call this method repeatedly.
  virtual void SetupDefaults();

  // Description:
  // Fill input port information.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Description:
  // Execute the pipeline, not used, but can instantiate filters for extra processing
  // in here.
  virtual int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  //
  vtkSplotchPainter         *SplotchPainter;
  vtkSplotchPainter         *LODSplotchPainter;
  vtkSplotchDefaultPainter  *SplotchDefaultPainter;
  vtkSplotchDefaultPainter  *LODSplotchDefaultPainter;
  //
  int                    ActiveParticleType;
  vtkSmartPointer<vtkStringArray> Settings;

private:
  vtkSplotchRepresentation(const vtkSplotchRepresentation&); // Not implemented
  void operator=(const vtkSplotchRepresentation&); // Not implemented
//ETX
};

#endif
