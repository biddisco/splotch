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

#include "vtksys/ios/sstream"

#include "vtkgl.h"
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

#ifdef VTK_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif

#include "splotch/scenemaker.h"
#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"
#include "cxxsupport/string_utils.h"

vtkInstantiatorNewMacro(vtkSplotchRaytraceMapper);

//----------------------------------------------------------------------------
vtkSplotchRaytraceMapper *vtkSplotchRaytraceMapper::New()
{
  vtkObject* ret = new vtkSplotchRaytraceMapper();
  return static_cast<vtkSplotchRaytraceMapper *>(ret);
}
// ---------------------------------------------------------------------------
vtkSplotchRaytraceMapper::vtkSplotchRaytraceMapper()
{
  this->TypeScalars      = NULL;
  this->ActiveScalars    = NULL;
  this->NumberOfParticleTypes = 0;
  this->SetNumberOfParticleTypes(1); 
  this->GrayAbsorption = 0.001;
  MPI_Manager::GetInstance();
}
// ---------------------------------------------------------------------------
vtkSplotchRaytraceMapper::~vtkSplotchRaytraceMapper()
{
  delete []this->TypeScalars;
  delete []this->ActiveScalars;
}
// ---------------------------------------------------------------------------
int vtkSplotchRaytraceMapper::FillInputPortInformation(int port,
  vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
// ---------------------------------------------------------------------------
double *vtkSplotchRaytraceMapper::GetBounds()
{
  this->GetBounds(this->Bounds);
  return this->Bounds;
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::GetBounds(double *bounds)
{
  //
  // Define box...
  //
  this->GetInput()->GetBounds(bounds);
  //
#ifdef VTK_USE_MPI
  vtkMPICommunicator *communicator = vtkMPICommunicator::SafeDownCast(
    vtkMultiProcessController::GetGlobalController()->GetCommunicator());
  if (communicator)
  {
    double mins[3] = {bounds[0], bounds[2], bounds[4]};
    double maxes[3] = {bounds[1], bounds[3], bounds[5]};
    double globalMins[3], globalMaxes[3];
    communicator->AllReduce(mins, globalMins, 3, vtkCommunicator::MIN_OP);
    communicator->AllReduce(maxes, globalMaxes, 3, vtkCommunicator::MAX_OP);
    bounds[0] = globalMins[0];  bounds[1] = globalMaxes[0];
    bounds[2] = globalMins[1];  bounds[3] = globalMaxes[1];
    bounds[4] = globalMins[2];  bounds[5] = globalMaxes[2];
  }
#endif
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::SetNumberOfParticleTypes(int N)
{
  this->NumberOfParticleTypes = std::max(N,this->NumberOfParticleTypes);
  this->IntensityScalars.resize(this->NumberOfParticleTypes,"");
  this->RadiusScalars.resize(this->NumberOfParticleTypes,"");
  this->Brightness.resize(this->NumberOfParticleTypes,10.5);
  this->LogIntensity.resize(this->NumberOfParticleTypes,0);
  this->TypeActive.resize(this->NumberOfParticleTypes,0);
  this->LogColour.resize(this->NumberOfParticleTypes,0);
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::SetIntensityScalars(int ptype, const char *s)
{
  if (std::string(s)!=this->IntensityScalars[ptype]) {
    this->IntensityScalars[ptype] = s;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
const char *vtkSplotchRaytraceMapper::GetIntensityScalars(int ptype)
{
  return this->IntensityScalars[ptype].c_str();
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::SetRadiusScalars(int ptype, const char *s)
{
  if (std::string(s)!=this->RadiusScalars[ptype]) {
    this->RadiusScalars[ptype] = s;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
const char *vtkSplotchRaytraceMapper::GetRadiusScalars(int ptype)
{
  return this->RadiusScalars[ptype].c_str();
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::SetBrightness(int ptype, double b)
{
  if (b!=this->Brightness[ptype]) {
    this->Brightness[ptype] = b;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
double vtkSplotchRaytraceMapper::GetBrightness(int ptype)
{
  return this->Brightness[ptype];
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::SetLogIntensity(int ptype, int l)
{
  if (l!=this->LogIntensity[ptype]) {
    this->LogIntensity[ptype] = l;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
int vtkSplotchRaytraceMapper::GetLogIntensity(int ptype)
{
  return this->LogIntensity[ptype];
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::SetTypeActive(int ptype, int a)
{
  if (a!=this->TypeActive[ptype]) {
    this->TypeActive[ptype] = a;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
int vtkSplotchRaytraceMapper::GetTypeActive(int ptype)
{
  return this->TypeActive[ptype];
}
// ---------------------------------------------------------------------------
// don't need this?
void vtkSplotchRaytraceMapper::SetLogColour(int ptype, int l)
{
  if (l!=this->LogColour[ptype]) {
    this->LogColour[ptype] = l;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
int vtkSplotchRaytraceMapper::GetLogColour(int ptype)
{
  return this->LogColour[ptype];
}
// ---------------------------------------------------------------------------
template <typename T>
std::string NumToStrSPM(T data) {
  vtksys_ios::ostringstream oss;
  oss.setf(0,ios::floatfield);
  oss.precision(5);  
  oss << data;
  return oss.str();
}
// ---------------------------------------------------------------------------
void vtkSplotchRaytraceMapper::Render(vtkRenderer *ren, vtkActor *act)
{
  int X = ren->GetSize()[0];
  int Y = ren->GetSize()[1];
  vtkPointSet *input = this->GetInput();
  vtkPoints *pts = input->GetPoints();
  //
  std::vector<vtkDataArray *> radiusarrays(this->NumberOfParticleTypes);
  std::vector<vtkDataArray *> intensityarrays(this->NumberOfParticleTypes);
  for (int i=0; i<this->NumberOfParticleTypes; i++) {
    radiusarrays[i] = (this->RadiusScalars[i].size()>0) ? 
      input->GetPointData()->GetArray(this->RadiusScalars[i].c_str()) : NULL;
    intensityarrays[i] = (this->IntensityScalars[i].size()>0) ? 
      input->GetPointData()->GetArray(this->IntensityScalars[i].c_str()) : NULL;  
  }
  //
  vtkDataArray *TypeArray = this->TypeScalars ? 
    input->GetPointData()->GetArray(this->TypeScalars) : NULL;  
  //
  vtkDataArray *ActiveArray = this->ActiveScalars ? 
    input->GetPointData()->GetArray(this->ActiveScalars) : NULL;  

  // For vertex coloring, this sets this->Colors.
  this->MapScalars( act->GetProperty()->GetOpacity() );

  // if one process has no points, pts will be NULL
  vtkIdType N = pts ? pts->GetNumberOfPoints() : 0;
  double bounds[6];
  input->GetBounds(bounds);
  double length = input->GetLength();
  double radius = N>0 ? length/N : length/1000.0;
  //
  std::vector<particle_sim> particle_data; // raw data 
  vec3 campos, lookat, sky;
  double zmin,zmax;
  ren->GetActiveCamera()->GetPosition(&campos.x);
  ren->GetActiveCamera()->GetFocalPoint(&lookat.x);
  ren->GetActiveCamera()->GetViewUp(&sky.x);
  ren->GetActiveCamera()->GetClippingRange(zmin, zmax);
  double FOV = ren->GetActiveCamera()->GetViewAngle();
  double newFOV = tan(vtkMath::RadiansFromDegrees(FOV/2.0))*X/Y;
  double splotchFOV = vtkMath::DegreesFromRadians(2.0*atan(newFOV));
  //
  unsigned char *cdata = this->Colors ? this->Colors->GetPointer(0) : NULL;
  particle_data.assign(N, particle_sim());

  double *brightness = &this->Brightness[0];

  vtkIdType activeParticles = 0;
  for (vtkIdType i=0; i<N; i++) {
    int ptype               = TypeArray ? TypeArray->GetTuple1(i) : 0;
    ptype = ptype<this->NumberOfParticleTypes ? ptype : this->NumberOfParticleTypes-1;
    particle_data[activeParticles].type   = ptype;
    bool active = this->TypeActive[ptype] && (ActiveArray ? (ActiveArray->GetTuple1(i)!=0) : 0);
    if (!active) continue;
    particle_data[activeParticles].active = active;
    //
    double *p = pts->GetPoint(i);
    particle_data[activeParticles].x      = p[0];
    particle_data[activeParticles].y      = p[1];
    particle_data[activeParticles].z      = p[2];
    particle_data[activeParticles].r      = radiusarrays[ptype] ? radiusarrays[ptype]->GetTuple1(i) : radius;
    particle_data[activeParticles].I      = intensityarrays[ptype] ? intensityarrays[ptype]->GetTuple1(i) : 1.0;
    if (cdata) {      
      particle_data[activeParticles].e.r = (cdata[i*4+0]/255.0);
      particle_data[activeParticles].e.g = (cdata[i*4+1]/255.0);
      particle_data[activeParticles].e.b = (cdata[i*4+2]/255.0);
    }
    else { // we don't support any other mode
      particle_data[activeParticles].e.r = 0.1;
      particle_data[activeParticles].e.g = 0.1;
      particle_data[activeParticles].e.b = 0.1;
      particle_data[activeParticles].I   = 1.0;
    }
    activeParticles++;
  }
  particle_data.resize(activeParticles);

  paramfile params;
  params.find("ptypes", this->NumberOfParticleTypes);
  params.find("xres", X);
  params.find("yres", Y);
  for (int i=0; i<this->NumberOfParticleTypes; i++) {
    std::string name;
    name = "intensity_log" + NumToStrSPM<int>(i);
    params.find(name, (this->LogIntensity[i]!=0));
    name = "brightness" + NumToStrSPM<int>(i);
    params.find(name, this->Brightness[i]);
  }
  params.find("gray_absorption", this->GrayAbsorption);

  params.find("zmin", zmin);
  params.find("zmax", zmax);
//  params.find("color_log0", true);
//  params.find("color_asinh0", false);
//  params.find("color_is_vector0", false);

  params.find("fov", splotchFOV);
  params.find("projection", true);
  params.find("minrad_pix", 1);
  params.find("a_eq_e", true);
  params.find("colorbar", false);

  params.find("quality_factor", 0.001);
  params.find("boost", true);

  particle_normalize(params, particle_data, true);

  arr2<COLOUR> pic(X,Y);

  if(particle_data.size()>0) {
    particle_project(params, particle_data, campos, lookat, sky);
  }
  for (int i=0; i<activeParticles; i++) {
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

  float32 grayabsorb;
  for (int i=0; i<this->NumberOfParticleTypes; i++) {
    grayabsorb = params.find<float32>("gray_absorption",this->GrayAbsorption);
  }
  render_new (particle_data, pic, a_eq_e, this->GrayAbsorption);

  MPI_Manager::GetInstance()->allreduceRaw
    (reinterpret_cast<float *>(&pic[0][0]),3*X*Y,MPI_Manager::Sum);

  exptable<float32> xexp(-20.0);
  if (MPI_Manager::GetInstance()->master() && a_eq_e) {
    std::cout << "Ïmage dimensions are " << X << "," << Y << std::endl;
    for (int ix=0;ix<X;ix++) {
      for (int iy=0;iy<Y;iy++) {
        pic[ix][iy].r = -xexp.expm1(pic[ix][iy].r);
        pic[ix][iy].g = -xexp.expm1(pic[ix][iy].g);
        pic[ix][iy].b = -xexp.expm1(pic[ix][iy].b);
      }
    }
  }

  if (!MPI_Manager::GetInstance()->master()) {
    for (int ix=0;ix<X;ix++) {
      for (int iy=0;iy<Y;iy++) {
        pic[ix][iy].r = 0.0;
        pic[ix][iy].g = 0.0;
        pic[ix][iy].b = 0.0;
      }
    }
  }



  if (MPI_Manager::GetInstance()->master()) {
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
  }
}
