/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSplotchRaytraceMapper.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSplotchRaytraceMapper.h"

#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCamera.h"

#include "vtkBitArray.h"
#include "vtkBoundingBox.h"
#include "vtkCompositeDataIterator.h"
#include "vtkCompositeDataSet.h"
#include "vtkDataArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkGraphicsFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTimerLog.h"
#include "vtkTransform.h"

#include <assert.h>
#include <vtkstd/vector>

//
#include "splotch/scenemaker.h"
#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"

vtkInstantiatorNewMacro(vtkSplotchRaytraceMapper);

template <class T>
static T vtkClamp(T val, T min, T max)
{
  val = val < min? min : val;
  val = val > max? max : val;
  return val;
}
//----------------------------------------------------------------------------
// return the correct type of vtkSplotchRaytraceMapper 
vtkSplotchRaytraceMapper *vtkSplotchRaytraceMapper::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = new vtkSplotchRaytraceMapper();
  return static_cast<vtkSplotchRaytraceMapper *>(ret);
}

// ---------------------------------------------------------------------------
// Construct object with scaling on, scaling mode is by scalar value,
// scale factor = 1.0, the range is (0,1), orient geometry is on, and
// orientation is by vector. Clamping and indexing are turned off. No
// initial sources are defined.
vtkSplotchRaytraceMapper::vtkSplotchRaytraceMapper()
{
}

// ---------------------------------------------------------------------------
vtkSplotchRaytraceMapper::~vtkSplotchRaytraceMapper()
{
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
int vtkSplotchRaytraceMapper::FillInputPortInformation(int port,
  vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::Render(vtkRenderer *ren, vtkActor *)
{
  int X = ren->GetSize()[0];
  int Y = ren->GetSize()[1];
  vtkPointSet *input = this->GetInput();
  vtkPoints *pts = input->GetPoints();
  //
  std::vector<particle_sim> particle_data; // raw data 
  vec3 campos, lookat, sky(0.0, 0.0, 0.0);
  ren->GetActiveCamera()->GetPosition(&campos.x);
  ren->GetActiveCamera()->GetFocalPoint(&lookat.x);
  std::vector<COLOURMAP> amap;
  //
  int N = pts->GetNumberOfPoints();
  particle_data.assign(N, particle_sim());
  for (int i=0; i<N; i++) {
    double *p = pts->GetPoint(i);
    particle_data[i].x = p[0];
    particle_data[i].y = p[1];
    particle_data[i].z = p[2];
  }
  arr2<COLOUR> pic(X,Y);
  paramfile params;

  int nColours = 3;
  double step = 1.0/(nColours-1.0);
  vec3 col[3] = { vec3(0, 0, 255)/255.0, vec3(128, 255, 128)/255.0, vec3(255, 0, 0)/255.0 };
  amap[0].addVal(0*step,COLOUR(col[0].x,col[0].y,col[0].z));
  amap[0].addVal(1*step,COLOUR(col[1].x,col[1].y,col[1].z));
  amap[0].addVal(2*step,COLOUR(col[2].x,col[2].y,col[2].z));

  if(particle_data.size()>0) {
    host_rendering(params, particle_data, pic, campos, lookat, sky, amap);
  }


  cerr << "Calling wrong render method!!\n";
}
