/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSplotchDefaultPainter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSplotchDefaultPainter - vtkDefaultPainter replacement that
//  inserts the vtkSplotchPainter at the correct position in the painter
//  chain.
//
// .SECTION Description
//  vtkSplotchDefaultPainter is a vtkDefaultPainter replacement
//  that inserts the vtkSplotchPainter at the correct position in the painter
//  chain. also It removes a few other painters that interfere with operation.
//
// .SECTION See Also
//  vtkDefaultPainter vtkSplotchPainter

#ifndef __vtkSplotchDefaultPainter_h
#define __vtkSplotchDefaultPainter_h

#include "pv_splotch_configure.h"
#include "vtkDefaultPainter.h"

class vtkSplotchPainter;

// If the MIP renderer is available
#ifdef PV_SPLOTCH_WITH_MIP
 class vtkMIPPainter;
#else 
 typedef vtkObject vtkMIPPainter;
#endif

// if CUDA rendering is enabled
#ifdef PV_SPLOTCH_USE_PISTON
 class vtkCUDASplotchPainter;
#else 
 typedef vtkObject vtkCUDASplotchPainter;
#endif

class pv_splotch_EXPORT vtkSplotchDefaultPainter : public vtkDefaultPainter
{
public:
  static vtkSplotchDefaultPainter* New();
  vtkTypeMacro(vtkSplotchDefaultPainter, vtkDefaultPainter);

  // Description:
  // Get/Set the Splotch painter.
  void SetSplotchPainter(vtkSplotchPainter*);
  vtkGetObjectMacro(SplotchPainter, vtkSplotchPainter);

  // Description:
  // The MIP painter must return the complete bounds of the whole dataset
  // not just the local 'piece', otherwise the compositing blanks out parts it thinks
  // are not covered by any geometry.
  void UpdateBounds(double bounds[6]);

  // Description:
  // Enable the use of MIP (if this is a LOD renderer)
  void SetUseMIP(int m);
  vtkGetMacro(UseMIP, int);

  // Description:
  // Enable/disble the use of GPU rendering with Piston
  void SetEnableCUDA(int m);
  vtkGetMacro(EnableCUDA, int);

//BTX
protected:
   vtkSplotchDefaultPainter();
  ~vtkSplotchDefaultPainter();

  // Description:
  // Setups the the painter chain.
  virtual void BuildPainterChain();

  // Description:
  // Take part in garbage collection.
  virtual void ReportReferences(vtkGarbageCollector *collector);
  virtual void ProcessInformation(vtkInformation* info);

  // Generic CPU based splotch painter
  vtkSplotchPainter     *SplotchPainter;
  // Special painter using MIP
  vtkMIPPainter         *MIPPainter;
  int                    UseMIP; 
  // CUDA version of splotch painter
  int                    EnableCUDA;
  vtkCUDASplotchPainter *CUDAPainter;

private:
  vtkSplotchDefaultPainter(const vtkSplotchDefaultPainter&); // Not implemented.
  void operator=(const vtkSplotchDefaultPainter&); // Not implemented.
//ETX
};

#endif
