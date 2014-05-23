/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDataSetToSplotch.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDataSetToSplotch - converts a DataSet to a PistonDataObject
// .SECTION Description
// Converts vtkDataSets that reside on the CPU into piston data that
// resides on the GPU. Afterward vtkPistonAlgorithms will processed
// it there.
//
// .SECTION See Also
// vtkPistonToDataSet

#ifndef __vtkDataSetToSplotch_h
#define __vtkDataSetToSplotch_h

#include "pv_splotch_configure.h"
#include "vtkDataSetToPiston.h"

class vtkDataSet;

class pv_splotch_EXPORT vtkDataSetToSplotch : public vtkDataSetToPiston
{
public:
  static vtkDataSetToSplotch *New();
  vtkTypeMacro(vtkDataSetToSplotch,vtkDataSetToPiston);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the name of the color array 
  vtkSetStringMacro(ScalarArrayName);
  vtkGetStringMacro(ScalarArrayName);

  // Description:
  // Set/Get the name of the radius array 
  vtkSetStringMacro(RadiusArrayName);
  vtkGetStringMacro(RadiusArrayName);

  // Description:
  // Set/Get the name of the intensity array 
  vtkSetStringMacro(IntensityArrayName);
  vtkGetStringMacro(IntensityArrayName);

protected:
  vtkDataSetToSplotch();
  ~vtkDataSetToSplotch();

  // Description:
  // Method that does the actual calculation. Funnels down to ExecuteData.
  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

  // Description:
  // Overridden to say that we require vtkDataSet inputs
  virtual int FillInputPortInformation(int, vtkInformation*);

  char *ScalarArrayName;
  char *RadiusArrayName;
  char *IntensityArrayName;

private:
  vtkDataSetToSplotch(const vtkDataSetToSplotch&);  // Not implemented.
  void operator=(const vtkDataSetToSplotch&);  // Not implemented.
};

#endif
