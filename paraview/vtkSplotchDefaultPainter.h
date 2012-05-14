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

#include "vtkDefaultPainter.h"

class vtkSplotchPainter;

class VTK_EXPORT vtkSplotchDefaultPainter : public vtkDefaultPainter
{
public:
  static vtkSplotchDefaultPainter* New();
  vtkTypeMacro(vtkSplotchDefaultPainter, vtkDefaultPainter);

  // Description:
  // Get/Set the Surface LIC painter.
  void SetSplotchPainter(vtkSplotchPainter*);
  vtkGetObjectMacro(SplotchPainter, vtkSplotchPainter);

  // Description:
  // The MIP painter must return the complete bounds of the whole dataset
  // not just the local 'piece', otherwise the compositing blanks out parts it thinks
  // are not covered by any geometry.
  void UpdateBounds(double bounds[6]);

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

  vtkSplotchPainter *SplotchPainter;
private:
  vtkSplotchDefaultPainter(const vtkSplotchDefaultPainter&); // Not implemented.
  void operator=(const vtkSplotchDefaultPainter&); // Not implemented.
//ETX
};

#endif
