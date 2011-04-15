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
vtkSplotchRaytraceMapper::vtkSplotchRaytraceMapper()
{
  this->ValueScalars     = NULL;
  this->IntensityScalars = NULL;
  this->RadiusScalars    = NULL;
  this->TypeScalars      = NULL;
  this->ActiveScalars    = NULL;
}

// ---------------------------------------------------------------------------
vtkSplotchRaytraceMapper::~vtkSplotchRaytraceMapper()
{
  delete []this->ValueScalars;
  delete []this->IntensityScalars;
  delete []this->RadiusScalars;
  delete []this->TypeScalars;
  delete []this->ActiveScalars;
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
void vtkSplotchRaytraceMapper::Render(vtkRenderer *ren, vtkActor *act)
{
  int X = ren->GetSize()[0];
  int Y = ren->GetSize()[1];
  vtkPointSet *input = this->GetInput();
  vtkPoints *pts = input->GetPoints();
  //
  vtkDataArray *RadiusArray = this->RadiusScalars ? 
    input->GetPointData()->GetArray(this->RadiusScalars) : NULL;
  //
  vtkDataArray *IntensityArray = this->IntensityScalars ? 
    input->GetPointData()->GetArray(this->IntensityScalars) : NULL;  
  //
  vtkDataArray *TypeArray = this->TypeScalars ? 
    input->GetPointData()->GetArray(this->TypeScalars) : NULL;  
  //
  vtkDataArray *ActiveArray = this->ActiveScalars ? 
    input->GetPointData()->GetArray(this->ActiveScalars) : NULL;  

  // For vertex coloring, this sets this->Colors as side effect.
  // For texture map coloring, this sets ColorCoordinates
  // and ColorTextureMap as a side effect.
  this->MapScalars( act->GetProperty()->GetOpacity() );


  //
  double N = pts->GetNumberOfPoints();
  double bounds[6];
  input->GetBounds(bounds);
  double length = input->GetLength();
  double radius = length/N;
  //
  std::vector<particle_sim> particle_data; // raw data 
  vec3 campos, lookat, sky;
  double zmin,zmax;
  ren->GetActiveCamera()->GetPosition(&campos.x);
  ren->GetActiveCamera()->GetFocalPoint(&lookat.x);
  ren->GetActiveCamera()->GetViewUp(&sky.x);
  ren->GetActiveCamera()->GetClippingRange(zmin, zmax);
  double FOV = ren->GetActiveCamera()->GetViewAngle();
  std::vector<COLOURMAP> amap;
//  amap.resize(1);
  //
  unsigned char *cdata = this->Colors ? this->Colors->GetPointer(0) : NULL;
  particle_data.assign(N, particle_sim());
  double imax = VTK_DOUBLE_MIN;
  double imin = VTK_DOUBLE_MAX;

  double brightness[2] = {10.5, 1.5};

  for (int i=0; i<N; i++) {
    double *p = pts->GetPoint(i);
    particle_data[i].x      = p[0];
    particle_data[i].y      = p[1];
    particle_data[i].z      = p[2];
    particle_data[i].type   = TypeArray ? TypeArray->GetTuple1(i) : 0;
    particle_data[i].r      = RadiusArray ? RadiusArray->GetTuple1(i) : radius;
    particle_data[i].active = ActiveArray ? (ActiveArray->GetTuple1(i)!=0) : 0;
    particle_data[i].I      = IntensityArray ? IntensityArray->GetTuple1(i) : 1.0;
    if (cdata) {      
      particle_data[i].e.r = (cdata[i*4+0]/255.0);
      particle_data[i].e.g = (cdata[i*4+1]/255.0);
      particle_data[i].e.b = (cdata[i*4+2]/255.0);
    }
    else { // we don't support any other mode
      particle_data[i].e.r = 0.1;
      particle_data[i].e.g = 0.1;
      particle_data[i].e.b = 0.1;
      particle_data[i].I   = 1.0;
    }
    if (particle_data[i].I>imax) imax = particle_data[i].I;
    if (particle_data[i].I<imin) imin = particle_data[i].I;
  }
  arr2<COLOUR> pic(X,Y);
  paramfile params;

  params.find("ptypes", 2);
  params.find("intensity_log0", true);
  params.find("color_log0", true);
  params.find("color_asinh0", false);
  params.find("color_is_vector0", false);
  params.find("xres", X);
  params.find("yres", Y);

  params.find("fov", FOV);
  params.find("projection", true);
  params.find("minrad_pix", 1);
  params.find("a_eq_e", true);
  params.find("zmin", zmin);
  params.find("zmax", zmax);
  params.find("brightness0", 10.5);

  params.find("intensity_max0", imax);
  params.find("intensity_min0", 0.0);
  params.find("gray_absorption", 0.0001);

  if(particle_data.size()>0) {
    particle_project(params, particle_data, campos, lookat, sky);
  }
  for (int i=0; i<N; i++) {
    if (particle_data[i].active) {
      if (cdata) {      
        double b = brightness[particle_data[i].type];
        particle_data[i].e.r *= particle_data[i].I*b;
        particle_data[i].e.g *= particle_data[i].I*b;
        particle_data[i].e.b *= particle_data[i].I*b;
      }
      else { // we don't support any other mode
        particle_data[i].e.r = 0.1;
        particle_data[i].e.g = 0.1;
        particle_data[i].e.b = 0.1;
        particle_data[i].I   = 1.0;
      }
    }
  }

  // ------------------------------------
  // -- Eliminating inactive particles --
  // ------------------------------------
  tsize npart_all;
  particle_eliminate(params, particle_data, npart_all);

  // --------------------------------
  // ----------- Sorting ------------
  // --------------------------------
  bool a_eq_e = params.find<bool>("a_eq_e",true);
  if (!a_eq_e) {
    int sort_type = params.find<int>("sort_type",1);
    particle_sort(particle_data,sort_type,true);
    }

  // ------------------------------------
  // ----------- Rendering ---------------
  // ------------------------------------
  float32 grayabsorb = params.find<float32>("gray_absorption",0.2);
  render_new (particle_data,pic,a_eq_e,grayabsorb);

  mpiMgr->allreduceRaw
    (reinterpret_cast<float *>(&pic[0][0]),3*X*Y,MPI_Manager::Sum);

  exptable<float32> xexp(-20.0);
  if (mpiMgr->master() && a_eq_e) {
    for (int ix=0;ix<X;ix++) {
      for (int iy=0;iy<Y;iy++) {
        pic[ix][iy].r=-xexp.expm1(pic[ix][iy].r);
        pic[ix][iy].g=-xexp.expm1(pic[ix][iy].g);
        pic[ix][iy].b=-xexp.expm1(pic[ix][iy].b);
      }
    }
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
