/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSplotchPainter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkSplotchPainter.h"

#include "vtksys/ios/sstream"

#include "vtkgl.h"
#include "vtkMapper.h"
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
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
#include "vtkTransform.h"
#include "vtkScalarsToColorsPainter.h"
#include "vtkCellArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
//
#ifdef VTK_USE_MPI
#include "vtkMPICommunicator.h"
#endif
#include "vtkMultiProcessController.h"

#include "splotch/scenemaker.h"
#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"
#include "cxxsupport/string_utils.h"

#undef min
#undef max
#include <algorithm>
#include <limits>

#include "vtkOpenGL.h"
#include "vtkgl.h"
#include "IceTConfig.h"

//----------------------------------------------------------------------------
vtkInstantiatorNewMacro(vtkSplotchPainter);
vtkCxxSetObjectMacro(vtkSplotchPainter, Controller, vtkMultiProcessController);
vtkCxxSetObjectMacro(vtkSplotchPainter, ScalarsToColorsPainter, vtkScalarsToColorsPainter);
//----------------------------------------------------------------------------
vtkSplotchPainter *vtkSplotchPainter::New()
{
  vtkObject* ret = new vtkSplotchPainter();
  return static_cast<vtkSplotchPainter *>(ret);
}
// ---------------------------------------------------------------------------
vtkSplotchPainter::vtkSplotchPainter()
{
  this->TypeScalars            = NULL;
  this->ActiveScalars          = NULL;
  // 
  this->NumberOfParticleTypes  = 1;
  this->SetNumberOfParticleTypes(0); 
  this->GrayAbsorption         = 0.001;
  this->RadiusMultiplier       = 1.0;
  this->ScalarsToColorsPainter = NULL;
  this->Controller             = NULL;
  this->EnableCUDA             = 0;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  //
  this->ArrayName = NULL;
  this->ArrayId = -1;
  this->ArrayComponent = 0;
  this->ArrayAccessMode = VTK_GET_ARRAY_BY_ID;
  //
  this->lastX = -1;
  this->lastY = -1;
  //
  MPI_Manager::GetInstance();
}
// ---------------------------------------------------------------------------
vtkSplotchPainter::~vtkSplotchPainter()
{
  delete []this->ArrayName;
  delete []this->TypeScalars;
  delete []this->ActiveScalars;
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::UpdateBounds(double bounds[6], vtkDataSet *input)
{
  if (!input) return;
  input->GetBounds(bounds);
  //
  if (this->Controller) {
    double mins[3]  = {bounds[0], bounds[2], bounds[4]};
    double maxes[3] = {bounds[1], bounds[3], bounds[5]};
    double globalMins[3], globalMaxes[3];
    this->Controller->AllReduce(mins, globalMins, 3, vtkCommunicator::MIN_OP);
    this->Controller->AllReduce(maxes, globalMaxes, 3, vtkCommunicator::MAX_OP);
    bounds[0] = globalMins[0];  bounds[1] = globalMaxes[0];
    bounds[2] = globalMins[1];  bounds[3] = globalMaxes[1];
    bounds[4] = globalMins[2];  bounds[5] = globalMaxes[2];
  }
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetNumberOfParticleTypes(int N)
{
  this->NumberOfParticleTypes = std::max(N,this->NumberOfParticleTypes);
  this->IntensityScalars.resize(this->NumberOfParticleTypes,"");
  this->RadiusScalars.resize(this->NumberOfParticleTypes,"");
  this->Brightness.resize(this->NumberOfParticleTypes,10.5);
  this->LogIntensity.resize(this->NumberOfParticleTypes,0);
  this->TypeActive.resize(this->NumberOfParticleTypes,1);
  this->LogColour.resize(this->NumberOfParticleTypes,0);
  this->MaxRadius.resize(this->NumberOfParticleTypes,0.0);
  intnorm_min.resize(this->NumberOfParticleTypes,0.0);
  intnorm_max.resize(this->NumberOfParticleTypes,0.0);
  colnorm_min.resize(this->NumberOfParticleTypes,0.0);
  colnorm_max.resize(this->NumberOfParticleTypes,0.0);  

  // Reinitialise array names 
  vtkDataObject *indo = this->GetInput();
  vtkPointSet *input = vtkPointSet::SafeDownCast(indo);

  radiusarrays.resize(this->NumberOfParticleTypes,NULL);
  intensityarrays.resize(this->NumberOfParticleTypes,NULL);
  // TODO: int i = previousNumberOfParticleTypes
  for (int i=0; i<this->NumberOfParticleTypes; i++) {
    radiusarrays[i] = (this->RadiusScalars[i].size()>0) ? 
      input->GetPointData()->GetArray(this->RadiusScalars[i].c_str()) : NULL;
    intensityarrays[i] = (this->IntensityScalars[i].size()>0) ? 
      input->GetPointData()->GetArray(this->IntensityScalars[i].c_str()) : NULL;  
  }
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetTypeActive(int ptype, int a)
{
  if (a!=this->TypeActive[ptype]) {
    this->TypeActive[ptype] = a;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
int vtkSplotchPainter::GetTypeActive(int ptype)
{
  return this->TypeActive[ptype];
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetIntensityScalars(int ptype, const char *s)
{
  if (std::string(s)!=this->IntensityScalars[ptype]) {
    this->IntensityScalars[ptype] = s;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
const char *vtkSplotchPainter::GetIntensityScalars(int ptype)
{
  return this->IntensityScalars[ptype].c_str();
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetRadiusScalars(int ptype, const char *s)
{
  if (std::string(s)!=this->RadiusScalars[ptype]) {
    this->RadiusScalars[ptype] = s;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
const char *vtkSplotchPainter::GetRadiusScalars(int ptype)
{
  return this->RadiusScalars[ptype].c_str();
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetRGBScalars(int ptype, const char *s)
{
  if (std::string(s)!=this->RGBScalars[ptype]) {
    this->RGBScalars[ptype] = s;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
const char *vtkSplotchPainter::GetRGBScalars(int ptype)
{
  return this->RGBScalars[ptype].c_str();
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetBrightness(int ptype, double b)
{
  if (b!=this->Brightness[ptype]) {
    this->Brightness[ptype] = b;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
double vtkSplotchPainter::GetBrightness(int ptype)
{
  return this->Brightness[ptype];
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetLogIntensity(int ptype, int l)
{
  if (l!=this->LogIntensity[ptype]) {
    this->LogIntensity[ptype] = l;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
int vtkSplotchPainter::GetLogIntensity(int ptype)
{
  return this->LogIntensity[ptype];
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// don't need this?
void vtkSplotchPainter::SetLogColour(int ptype, int l)
{
  if (l!=this->LogColour[ptype]) {
    this->LogColour[ptype] = l;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
int vtkSplotchPainter::GetLogColour(int ptype)
{
  return this->LogColour[ptype];
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::SetMaxRadius(int ptype, double r)
{
  if (r!=this->MaxRadius[ptype]) {
    this->MaxRadius[ptype] = r;
    this->Modified();
  }
}
// ---------------------------------------------------------------------------
double vtkSplotchPainter::GetMaxRadius(int ptype)
{
  return this->MaxRadius[ptype];
}

//-----------------------------------------------------------------------------
void vtkSplotchPainter::ProcessInformation(vtkInformation* info)
{
  info->Set(vtkScalarsToColorsPainter::INTERPOLATE_SCALARS_BEFORE_MAPPING(), 0);

  if (info->Has(vtkScalarsToColorsPainter::SCALAR_MODE()))
    {
    this->SetScalarMode(info->Get(vtkScalarsToColorsPainter::SCALAR_MODE()));
    }

  if (info->Has(vtkScalarsToColorsPainter::ARRAY_ACCESS_MODE()))
    {
    this->SetArrayAccessMode(info->Get(vtkScalarsToColorsPainter::ARRAY_ACCESS_MODE()));
    }

  if (info->Has(vtkScalarsToColorsPainter::ARRAY_ID()))
    {
    this->SetArrayId(info->Get(vtkScalarsToColorsPainter::ARRAY_ID()));
    }

  if (info->Has(vtkScalarsToColorsPainter::ARRAY_NAME()))
    {
    this->SetArrayName(info->Get(vtkScalarsToColorsPainter::ARRAY_NAME()));
    }

  if (info->Has(vtkScalarsToColorsPainter::ARRAY_COMPONENT()))
    {
    this->SetArrayComponent(info->Get(vtkScalarsToColorsPainter::ARRAY_COMPONENT()));
    }
  }
// ---------------------------------------------------------------------------
template <typename T>
std::string vtkSplotchPainter::NumToStrSPM(T data) {
  vtksys_ios::ostringstream oss;
  oss.precision(5);  
  oss << data;
  return oss.str();
}
//----------------------------------------------------------------------------
void FloatOrDoubleArrayPointer(vtkDataArray *dataarray, float *&F, double *&D) {
  if (dataarray && vtkFloatArray::SafeDownCast(dataarray)) {
    F = vtkFloatArray::SafeDownCast(dataarray)->GetPointer(0);
    D = NULL;
  }
  if (dataarray && vtkDoubleArray::SafeDownCast(dataarray)) {
    D = vtkDoubleArray::SafeDownCast(dataarray)->GetPointer(0);
    F = NULL;
  }
  //
  if (dataarray && !F && !D) {
    vtkGenericWarningMacro(<< dataarray->GetName() << "must be float or double");
  }
}
//----------------------------------------------------------------------------
#define FloatOrDouble(F, D, index) F ? F[index] : D[index]
#define FloatOrDoubleSet(F, D) ((F!=NULL) || (D!=NULL))
//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// IceT is not exported by paraview, so rather than force lots of include dirs
// and libs, just manually set some defs which will keep the compiler happy
//-----------------------------------------------------------------------------
typedef IceTUnsignedInt32       IceTEnum;
typedef IceTInt32               IceTInt;
typedef void *                  IceTContext;
//#ifdef USE_ICET
 extern "C" ICET_EXPORT void icetGetIntegerv(IceTEnum pname, IceTInt *params);
 extern "C" ICET_EXPORT IceTContext icetGetContext(void);
 extern "C" ICET_EXPORT ICET_EXPORT void icetCompositeMode(IceTEnum mode);
 extern "C" ICET_EXPORT void icetSetColorFormat(IceTEnum color_format);
 extern "C" ICET_EXPORT void icetSetDepthFormat(IceTEnum depth_format);
//#endif

#define ICET_STATE_ENGINE_START (IceTEnum)0x00000000
#define ICET_NUM_TILES          (ICET_STATE_ENGINE_START | (IceTEnum)0x0010)
#define ICET_TILE_VIEWPORTS     (ICET_STATE_ENGINE_START | (IceTEnum)0x0011)
#define ICET_COMPOSITE_MODE_Z_BUFFER    (IceTEnum)0x0301
#define ICET_COMPOSITE_MODE_BLEND       (IceTEnum)0x0302
#define ICET_IMAGE_COLOR_RGBA_UBYTE     (IceTEnum)0xC001
#define ICET_IMAGE_COLOR_RGBA_FLOAT     (IceTEnum)0xC002
#define ICET_IMAGE_COLOR_NONE           (IceTEnum)0xC000

#define ICET_IMAGE_DEPTH_FLOAT          (IceTEnum)0xD001
#define ICET_IMAGE_DEPTH_NONE           (IceTEnum)0xD000

//-----------------------------------------------------------------------------
void vtkSplotchPainter::PrepareForRendering(vtkRenderer* ren, vtkActor* actor)
{
//       icetSetColorFormat(ICET_IMAGE_COLOR_RGBA_UBYTE);
//      icetSetDepthFormat(ICET_IMAGE_DEPTH_NONE);
//     icetCompositeMode(ICET_COMPOSITE_MODE_BLEND);
  //
  // Unneeded because we use iceT window sizes below
  //X = ren->GetSize()[0];
  //Y = ren->GetSize()[1];
  // We arent using cuda
  params.setParam("cuda_paraview_splotch", false);
  // Get input dataset
  vtkDataObject *indo = this->GetInput();
  vtkPointSet *input = vtkPointSet::SafeDownCast(indo);
  vtkPoints *pts = input->GetPoints();
  
  // Init arrays (can use different radius array per ptype etc)
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
   // Active array should not be changeable, as this is used internally by splotch
  // for clipping.(and now also for type filtering in splotch-paraview)
  // vtkDataArray *ActiveArray = this->ActiveScalars ? 
  //   input->GetPointData()->GetArray(this->ActiveScalars) : NULL;  

  //
  // Get the LUT and scalar array
  //
  colourspresent = false;
  int cellFlag=0;
  vtkDataSet* ds = static_cast<vtkDataSet*>(input);
  vtkDataArray* scalars = vtkAbstractMapper::GetScalars(ds,
    this->ScalarMode, this->ArrayAccessMode, this->ArrayId,
    this->ArrayName, cellFlag);

  // if scalars are a RGB colour table, then we don't need to map them.
  cdata3 = NULL;
  cdata4 = NULL;
  if (scalars && scalars->GetDataType()==VTK_UNSIGNED_CHAR && scalars->GetNumberOfComponents()==3) {
    cdata3 = static_cast<unsigned char*>(scalars->GetVoidPointer(0));
    colourspresent = true;
  }
  else {
    vtkScalarsToColors *lut = this->ScalarsToColorsPainter->GetLookupTable();
    //
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkUnsignedCharArray::SafeDownCast(input->GetPointData()->GetScalars());
    cdata4 = colors ? colors->GetPointer(0) : NULL;
    colourspresent = (cdata4!=NULL);
  }

  // We need the viewport/viewsize scaled by the Image Reduction Factor when downsampling
  // with client server. This is a nasty hack because we can't access this information
  // directly.
  // This is the reported size of final image, (which may be wrong)
  int viewsize[2], vieworigin[2];
  ren->GetTiledSizeAndOrigin( &viewsize[0],   &viewsize[1], 
                              &vieworigin[0], &vieworigin[1] );
  // Query IceT for the actual size
  IceTInt ids, vp[32*4] = {0, 0, viewsize[0], viewsize[1],};
#ifdef USE_ICET
  if (icetGetContext()!=NULL) {
    icetGetIntegerv(ICET_NUM_TILES,&ids);
    // when running on a single core, this returns nonsense
    if (ids>0 && ids<32) {
      icetGetIntegerv(ICET_TILE_VIEWPORTS,vp);
    }
  }
#endif 
  // Here we compute the actual viewport scaling factor with the correct adjusted sizes.
  double viewPortRatio[2];
  double *viewPort = ren->GetViewport();
  viewPortRatio[0] = (vp[2]*(viewPort[2]-viewPort[0])) / 2.0 + viewsize[0]*viewPort[0];
  viewPortRatio[1] = (vp[3]*(viewPort[3]-viewPort[1])) / 2.0 + viewsize[1]*viewPort[1];
  // Oops, we must use the IceT sizes not the renderwindow sizes.
  X = vp[2];
  Y = vp[3];

  //
  // We need the transform that reflects the transform point coordinates according to actor's transformation matrix
  //
  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->DeepCopy(
    ren->GetActiveCamera()->GetCompositeProjectionTransformMatrix(ren->GetTiledAspectRatio(),
    0,1));

  //
  // watch out, if one process has no points, pts array will be NULL
  //
  N = pts ? pts->GetNumberOfPoints() : 0;
  float *pointsF = NULL;
  double *pointsD = NULL;
  if (N>0) {
    FloatOrDoubleArrayPointer(pts->GetData(), pointsF, pointsD);
  }

  ren->GetActiveCamera()->GetPosition(&campos.x);
  ren->GetActiveCamera()->GetFocalPoint(&lookat.x);
  ren->GetActiveCamera()->GetViewUp(&sky.x);
  ren->GetActiveCamera()->GetClippingRange(zmin, zmax);
  FOV = ren->GetActiveCamera()->GetViewAngle();
  newFOV = tan(vtkMath::RadiansFromDegrees(FOV/2.0))*X/Y;
  splotchFOV = vtkMath::DegreesFromRadians(2.0*atan(newFOV));

  double bounds[6];
  input->GetBounds(bounds);
  double length = input->GetLength();
  double radius = N>0 ? length/N : length/1000.0;

  this->particle_compute = false;
  if (X!=lastX || Y!=lastY) {
    pic.alloc(X,Y);
    this->particle_compute = true;
    lastX = X;
    lastY = Y;
  }

  if (!this->EnableCUDA || this->GetMTime()>ParticleDataComputeTime.GetMTime()) {
    std::cout << "Modified - need to recompute particle data for " << N << " pts" << std::endl;
    ParticleDataComputeTime.Modified();
    this->particle_compute = true;
  }

//  brightness = &this->Brightness[0];

  if (this->particle_compute) {
    particle_data.resize(N, particle_sim());
    // for openmp, disable activeparticles
    //  vtkIdType activeParticles = 0;
  #define activeParticles i
  #pragma omp parallel for
    for (vtkIdType i=0; i<N; i++) {
      // Check particle type
      int ptype =  TypeArray ? TypeArray->GetTuple1(i) : 0;

      // If this is a new type we must update to handle this
      if(ptype >= this->NumberOfParticleTypes)
      {
        
        #pragma omp critical (update_type)
        {
          std::cout << "Adding new type" << std::endl;
          // Double check in case another thread modified it before we entered critical section
          if(ptype >= this->NumberOfParticleTypes)
          {
            this->SetNumberOfParticleTypes(ptype+1);
          }
        }
      }

      // Set type and active status (if filtered by type)
      particle_data[i].type   = ptype;
      particle_data[i].active = this->TypeActive[ptype];

      // Skip this particle if inactive, otherwise continue to setup parameters
      if (!particle_data[i].active) continue;
      
      particle_data[activeParticles].active = true; 

      if (pointsF) {
        particle_data[activeParticles].x = pointsF[i*3+0];
        particle_data[activeParticles].y = pointsF[i*3+1];
        particle_data[activeParticles].z = pointsF[i*3+2];
      }
      else {
        particle_data[activeParticles].x = pointsD[i*3+0];
        particle_data[activeParticles].y = pointsD[i*3+1];
        particle_data[activeParticles].z = pointsD[i*3+2];
      }

      double radiusdata[1];
      // Limit to max radius (mr) for this type if mr > 0
      double r;
      double mr = this->MaxRadius[particle_data[activeParticles].type];
      if (radiusarrays[ptype]) {
        radiusarrays[ptype]->GetTuple(i,radiusdata);
        r = radiusdata[0]*this->RadiusMultiplier;
      }
      else {
        r = radius*this->RadiusMultiplier;
      }
      particle_data[activeParticles].r = ( r > mr && mr > 0) ? mr : r;

      double intensitydata[1];
      if (intensityarrays[ptype]) {
        intensityarrays[ptype]->GetTuple(i,intensitydata);
        // if (this->LogIntensity[ptype]) {
        //   particle_data[activeParticles].I = pow(10,intensitydata[0]);
        // }
        // else {
          particle_data[activeParticles].I = intensitydata[0];
        // }
      }
      else {
        particle_data[activeParticles].I = 1.0;
      }
      if (cdata3) { // RGB values have been generated by scalars to colours
        particle_data[activeParticles].e.r = (cdata3[i*3+0]/255.0);
        particle_data[activeParticles].e.g = (cdata3[i*3+1]/255.0);
        particle_data[activeParticles].e.b = (cdata3[i*3+2]/255.0);
      }
      else if (cdata4) { // RGBA values have been generated by scalars to colours
        particle_data[activeParticles].e.r = (cdata4[i*4+0]/255.0);
        particle_data[activeParticles].e.g = (cdata4[i*4+1]/255.0);
        particle_data[activeParticles].e.b = (cdata4[i*4+2]/255.0);
      }
      else { // we don't support any other mode for now
        particle_data[activeParticles].e.r = 0.1;
        particle_data[activeParticles].e.g = 0.1;
        particle_data[activeParticles].e.b = 0.1;
      }
    }

    int N = this->NumberOfParticleTypes;
    MPI_Manager::GetInstance()->allreduceRaw<int>(&N, 1, MPI_Manager::Max);
    this->SetNumberOfParticleTypes(N);

    std::cout << "Particle data recomputed successfully" << std::endl; 
  }

  this->RenderSplotchParams(ren, actor);
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::RenderSplotchParams(vtkRenderer* ren, vtkActor* actor)
{
  params.setParam("ptypes", this->NumberOfParticleTypes);
  params.setParam("xres", X);
  params.setParam("yres", Y);
  for (int i=0; i<this->NumberOfParticleTypes; i++) {
    std::string name;
    name = "intensity_log" + NumToStrSPM<int>(i);
    params.setParam(name, (this->LogIntensity[i]!=0));
    name = "brightness" + NumToStrSPM<int>(i);
    params.setParam(name, this->Brightness[i]);
  }
  params.setParam("gray_absorption", this->GrayAbsorption);
  params.setParam("zmin", 0.0); // zmin - (zmax-zmin)/1.0);
  params.setParam("zmax", 1.e23); //zmax + (zmax-zmin)/1.0);
  params.setParam("fov",  splotchFOV);
  params.setParam("projection", true);
  params.setParam("minrad_pix", 1);
  params.setParam("a_eq_e", true);
  params.setParam("colorbar", false);
  params.setParam("quality_factor", 0.001);
  params.setParam("boost", false);

//  params.setParam("color_min0",  0.0);
//  params.setParam("color_max0",  1.0);

  // If particles are recomputed, normalization values are in param struct
  // Otherwise we write the ones we already have to the param struct for use later
  if (this->particle_compute) {
    host_funct::particle_normalize2(params, particle_data, true);
    for(int t = 0; t < this->NumberOfParticleTypes; t++)
    {
      if(params.param_present("intensity_min"+dataToString(t)))
        intnorm_min[t] = params.find<float>("intensity_min"+dataToString(t));
      if(params.param_present("intensity_max"+dataToString(t)))
        intnorm_max[t] = params.find<float>("intensity_max"+dataToString(t));
      if(params.param_present("color_min"+dataToString(t)))
        colnorm_min[t] = params.find<float>("color_min"+dataToString(t));
      if(params.param_present("color_max"+dataToString(t)))
        colnorm_max[t] = params.find<float>("color_max"+dataToString(t));
    }
  }
  else {
    for(int t = 0; t < this->NumberOfParticleTypes; t++)
    {
      params.setParam("intensity_min"+dataToString(t), intnorm_min[t]);
      params.setParam("intensity_max"+dataToString(t), intnorm_max[t]);
      params.setParam("color_min"+dataToString(t), colnorm_min[t]);
      params.setParam("color_max"+dataToString(t), colnorm_max[t]);
    }
  }
}

// ---------------------------------------------------------------------------
void vtkSplotchPainter::RenderInternal(vtkRenderer* ren, vtkActor* actor, 
  unsigned long typeflags, bool forceCompileOnly)
{
  if (this->particle_compute) {
    if(particle_data.size()>0) {
      host_funct::particle_project(params, particle_data, campos, lookat, sky);
    }

  //#pragma omp parallel for
    for (int i=0; i<N /*activeParticles*/; i++) {
      if (colourspresent) {
        double b = this->Brightness[particle_data[i].type];
        particle_data[i].e.r *= particle_data[i].I*b;
        particle_data[i].e.g *= particle_data[i].I*b;
        particle_data[i].e.b *= particle_data[i].I*b;
      }
      else { // we don't support any other mode
        particle_data[i].e.r = 0.9;
        particle_data[i].e.g = 0.9;
        particle_data[i].e.b = 0.9;
        particle_data[i].I   = 1.0;
      }
    }

    // ------------------------------------
    // -- Eliminating inactive particles --
    // ------------------------------------
//    tsize npart_all;
//    host_funct::particle_eliminate(params, particle_data, npart_all);

    // --------------------------------
    // ----------- Sorting ------------
    // --------------------------------
    a_eq_e = params.find<bool>("a_eq_e",true);
    if (!a_eq_e) {
      int sort_type = params.find<int>("sort_type",1);
      host_funct::particle_sort(particle_data,sort_type,true);
    }
  }

  // ------------------------------------
  // ----------- Rendering ---------------
  // ------------------------------------

  float32 grayabsorb;
  for (int i=0; i<this->NumberOfParticleTypes; i++) {
    grayabsorb = params.find<float32>("gray_absorption",this->GrayAbsorption);
  }

  if (particle_data.size()>0) {
    host_funct::render_new(&(particle_data[0]), particle_data.size(), pic, a_eq_e, this->GrayAbsorption);
  }

  this->PostRenderCompositing(ren, actor);
}
// ---------------------------------------------------------------------------
void vtkSplotchPainter::PostRenderCompositing(vtkRenderer* ren, vtkActor* actor)
{
  // 
  // MPI_Manager::GetInstance()->allreduceRaw
  //  (reinterpret_cast<float *>(&pic[0][0]),3*X*Y,MPI_Manager::Sum);

  //if (MPI_Manager::GetInstance()->master() && a_eq_e) {
    if(a_eq_e) {
    // std::cout << "Image dimensions are " << X << "," << Y << std::endl;
    //
    float global_min=std::numeric_limits<double>::max();
    float global_max=std::numeric_limits<double>::min();
#pragma omp parallel 
    {
      float local_min=std::numeric_limits<double>::max();
      float local_max=std::numeric_limits<double>::min();
#pragma omp for nowait
      for (int ix=0;ix<X;ix++) {
        for (int iy=0;iy<Y;iy++) {
          local_min = std::min(pic[ix][iy].r, local_min);
          local_min = std::min(pic[ix][iy].g, local_min);
          local_min = std::min(pic[ix][iy].b, local_min);
          local_max = std::max(pic[ix][iy].r, local_max);
          local_max = std::max(pic[ix][iy].g, local_max);
          local_max = std::max(pic[ix][iy].b, local_max);
        }
      }
#pragma omp critical 
      {
        global_min = std::min(global_min, local_min);
        global_max = std::max(global_max, local_max);
      }
    }

    // std::cout << "global_min, global_max are {" << global_min << "," << global_max << "}" << std::endl;
    exptable<float32> xexp(global_min);
    //
#pragma omp parallel for
    for (int ix=0;ix<X;ix++) {
      for (int iy=0;iy<Y;iy++) {
        pic[ix][iy].r = -xexp.expm1(pic[ix][iy].r);
        pic[ix][iy].g = -xexp.expm1(pic[ix][iy].g);
        pic[ix][iy].b = -xexp.expm1(pic[ix][iy].b);
      }
    }
  }

//   if (!MPI_Manager::GetInstance()->master()) {
// #pragma omp parallel for
//     for (int ix=0;ix<X;ix++) {
//       for (int iy=0;iy<Y;iy++) {
//         pic[ix][iy].r = 0.0;
//         pic[ix][iy].g = 0.0;
//         pic[ix][iy].b = 0.0;
//       }
//     }
//   }

 // if (MPI_Manager::GetInstance()->master()) {
    //
    // copy to OpenGL image buffer
    //
    int viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(viewport[0], viewport[2], viewport[1], viewport[3], -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // we draw our image just in front of the back clipping plane, 
    // so all other geometry will appear in front of it. (z = -0.99)
    // glColor4f(0.0, 0.0, 0.0, 0.5);
    //glPixelTransferf( GL_RED_SCALE, 1.0);
    //glPixelTransferf( GL_RED_BIAS,  0.0);
    //glPixelTransferf( GL_GREEN_SCALE, 1.0);
    //glPixelTransferf( GL_GREEN_BIAS,  0.0);
    //glPixelTransferf( GL_BLUE_SCALE, 1.0);
    //glPixelTransferf( GL_BLUE_BIAS,  0.0);

    // option 1
    // us glPixelTransfer to multiply pixels as we copy them and thus set alpha
    glPixelTransferf( GL_ALPHA_SCALE, 0.0025);
    glPixelTransferf( GL_ALPHA_BIAS,  0.0);

    // option 2
    // create e temporary copy of the image and set alpha by hand
/*
    struct RGBA_TUPLE {
      float r,g,b,a;
    };

    arr2<RGBA_TUPLE> rgba_pic;
    rgba_pic.alloc(X,Y);

#pragma omp parallel for
    for (int ix=0;ix<X;ix++) {
      for (int iy=0;iy<Y;iy++) {
        rgba_pic[ix][iy].r = pic[ix][iy].r;
        rgba_pic[ix][iy].g = pic[ix][iy].g;
        rgba_pic[ix][iy].b = pic[ix][iy].b;
        // lower than this and everythingbecomes invisible.
        // @todo, check if iceT is skipping the composite when alpha is low
        rgba_pic[ix][iy].a = 0.0025;
      }
    }
*/
    GLboolean on = glIsEnabled(GL_BLEND);
    glDisable(GL_BLEND);

    for (int i=0; i<X; i++) {
      glRasterPos3f(X-1-i, 0, -0.99);

///* original RGB
      COLOUR *ptr = &pic[i][0];
      float *x0 = &ptr->r;
      glDrawPixels(1, Y, (GLenum)(GL_RGB), (GLenum)(GL_FLOAT), (GLvoid*)(x0));
//*/
/*
      RGBA_TUPLE *ptr = &rgba_pic[i][0];
      float *x0 = &ptr->r;
      glDrawPixels(1, Y, (GLenum)(GL_RGBA), (GLenum)(GL_FLOAT), (GLvoid*)(x0));
*/
    }
    glPixelTransferf( GL_ALPHA_SCALE, 1.0);
    glPixelTransferf( GL_ALPHA_BIAS,  0.0);

    if (on) {
      glEnable(GL_BLEND);
    }

    glMatrixMode( GL_MODELVIEW );   
    glPopMatrix();
    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
 // }
}

