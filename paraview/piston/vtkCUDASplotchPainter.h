/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCUDASplotchPainter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCUDASplotchPainter - this painter paints polygons.
// .SECTION Description
// This painter renders Polys in vtkPolyData. It can render the polys
// in any representation (VTK_POINTS, VTK_WIREFRAME, VTK_SURFACE).

#ifndef __vtkCUDASplotchPainter_h
#define __vtkCUDASplotchPainter_h

#include "pv_splotch_configure.h"
#include "vtkCUDAPiston.h"
#include "vtkSplotchPainter.h"
#include "vtkSmartPointer.h"
#include "vtkDataSetToPiston.h"

class vtkRenderWindow;
class vtkCamera;
class vtkActor;

#define VTK_DIRECTION_NO_SORT      -1
#define VTK_DIRECTION_BACK_TO_FRONT 0
#define VTK_DIRECTION_FRONT_TO_BACK 1

class pv_splotch_EXPORT vtkCUDASplotchPainter : public vtkSplotchPainter
{
public:
  static vtkCUDASplotchPainter* New();
  vtkTypeMacro(vtkCUDASplotchPainter, vtkSplotchPainter);

  // Description:
  // Manually call this before any cuda filters are created
  // to use direct GPU rendering.
  static int InitCudaGL(vtkRenderWindow *rw, int rank, int &displayId);

  // Description:
  // Return true if using cuda interop feature otherwise false.
  inline static bool IsEnabledCudaGL()
    {
    return CudaGLInitted;
    }

  // Description:
  // Release any graphics resources that are being consumed by this mapper.
  // The parameter window could be used to determine which graphic
  // resources to release.
  virtual void ReleaseGraphicsResources(vtkWindow *) {};

  void RenderOnGPU(vtkCamera *cam, vtkActor *act);

protected:
  vtkCUDASplotchPainter();
  ~vtkCUDASplotchPainter();

  void PrepareForRendering(vtkRenderer* renderer, vtkActor* actor);

  vtkSmartPointer<vtkDataSetToPiston>            DataSetToPiston;
  vtkCamera                                     *Camera;
  vtkActor                                      *Actor;

protected:
  // Description:
  // Allocates buffers that are shared between CUDA and GL
  void PrepareDirectRenderBuffers(int nPoints, int nCells);

  static bool CudaGLInitted;

  class InternalInfo;
  InternalInfo *Internal;

private:
  vtkCUDASplotchPainter(const vtkCUDASplotchPainter&); // Not implemented.
  void operator=(const vtkCUDASplotchPainter&); // Not implemented.

};


#endif

