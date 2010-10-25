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
#include "vtkgl.h"

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
  vec3 campos, lookat, sky;
  ren->GetActiveCamera()->GetPosition(&campos.x);
  ren->GetActiveCamera()->GetFocalPoint(&lookat.x);
  ren->GetActiveCamera()->GetViewUp(&sky.x);
  double FOV = ren->GetActiveCamera()->GetViewAngle();
  std::vector<COLOURMAP> amap;
  amap.resize(1);
  //
  double N = pts->GetNumberOfPoints();
  particle_data.assign(N, particle_sim());
  double imax = 0;
  for (double i=0; i<N; i++) {
    double *p = pts->GetPoint(i);
    particle_data[i].x = p[0];
    particle_data[i].y = p[1];
    particle_data[i].z = p[2];
    particle_data[i].type = 0;
    particle_data[i].e.r = 0.25 + 0.5*(i/N);
    particle_data[i].e.g = 0.25 + 0.75*(i/N);
    particle_data[i].e.b = 0.25 + 0.25*(i/N);
    particle_data[i].I   = 0.25 + 0.5*(i/N);
    particle_data[i].r   = 0.05;

    if (particle_data[i].I>imax) imax = particle_data[i].I;
  }
  arr2<COLOUR> pic(X,Y);
  paramfile params;

  int nColours = 3;
  double step = 1.0/(nColours-1.0);
  vec3 col[3] = { vec3(0, 0, 255)/255.0, vec3(128, 255, 128)/255.0, vec3(255, 0, 0)/255.0 };
  amap[0].addVal(0*step,COLOUR(col[0].x,col[0].y,col[0].z));
  amap[0].addVal(1*step,COLOUR(col[1].x,col[1].y,col[1].z));
  amap[0].addVal(2*step,COLOUR(col[2].x,col[2].y,col[2].z));
  amap[0].sortMap();

/*
Parser: intensity_log0 = T <default>
Parser: color_log0 = T <default>
Parser: color_asinh0 = F <default>
Parser: color_is_vector0 = F <default>
*/

  params.find("intensity_log0", false);
  params.find("color_log0", false);
  params.find("color_asinh0", false);
  params.find("color_is_vector0", false);
  params.find("xres", X);
  params.find("yres", Y);

  params.find("fov", FOV);
  params.find("projection", true);
  params.find("minrad_pix", 1);
  params.find("a_eq_e", true);
  params.find("zmin", 0.1);
  params.find("zmax", 50);
  params.find("brightness0", 1);

  params.find("intensity_max0", imax);
  params.find("intensity_min0", 0.0);

  if(particle_data.size()>0) {
    host_rendering(params, particle_data, pic, campos, lookat, sky, amap);
  }

  bool master = true;
  bool a_eq_e = true;
//    mpiMgr.allreduceRaw
//      (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres,MPI_Manager::Sum);

    exptable<float32> xexp(-20.0);
    if (master && a_eq_e)
      for (int ix=0;ix<X;ix++)
        for (int iy=0;iy<Y;iy++)
          {
          pic[ix][iy].r=-xexp.expm1(pic[ix][iy].r);
          pic[ix][iy].g=-xexp.expm1(pic[ix][iy].g);
          pic[ix][iy].b=-xexp.expm1(pic[ix][iy].b);
          }
/*
  for (int i=0; i<X; i++) {
    for (int j=0; j<Y; j++) {
      if (i==j) pic(i,j) = COLOUR(1,1,1);
    }
  }
*/

    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(viewport[0], viewport[2], viewport[1], viewport[3], -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  for (int i=0; i<X; i++) {
    glRasterPos2i(X-1-i, 0);
    COLOUR *ptr = &pic[i][0];
    float *x0 = &ptr->r;
    glDrawPixels(1, Y, (GLenum)(GL_RGB), (GLenum)(GL_FLOAT), (GLvoid*)(x0));
  }

    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );   
    glPopMatrix();



 	

//  cerr << "Calling wrong render method!!\n";
}
