/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCUDASplotchPainter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCUDASplotchPainter.h"

#include "vtkgl.h"
#include "vtkOpenGLExtensionManager.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkTransform.h"
#include "vtkPainterDeviceAdapter.h"
//
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
//
#include "vtkObjectFactory.h"
#include "vtkTimerLog.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositeDataIterator.h"
//
#include <thrust/version.h>
//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkCUDASplotchPainter);
//-----------------------------------------------------------------------------
#define NUM_INTEROP_BUFFERS 4
bool vtkCUDASplotchPainter::CudaGLInitted = false;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
class vtkCUDASplotchPainter::InternalInfo {
public:
  InternalInfo() {
    this->BufferSize = 0;
    this->CellCount = 0;
    this->DataObjectMTimeCache = 0;
    this->clearBuffers();
  }
  ~InternalInfo() {
    this->clearBuffers();
  }
  void clearBuffers() {
    if (this->BufferSize != 0 && this->vboBuffers[0] != -1) {
      vtkgl::DeleteBuffers(NUM_INTEROP_BUFFERS, this->vboBuffers);
    }
    for (int i=0; i<NUM_INTEROP_BUFFERS; i++) {
      vboBuffers[i] = -1;
      vboResources[i] = NULL;
    }
  }

  int     BufferSize;
  int     CellCount;
  GLuint  vboBuffers[NUM_INTEROP_BUFFERS];
  struct  cudaGraphicsResource* vboResources[NUM_INTEROP_BUFFERS];
  //
  unsigned long DataObjectMTimeCache;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vtkCUDASplotchPainter::vtkCUDASplotchPainter()
{
  this->DataSetToPiston  = vtkSmartPointer<vtkDataSetToPiston>::New();
  //
  this->Internal = new vtkCUDASplotchPainter::InternalInfo();
}
//-----------------------------------------------------------------------------
vtkCUDASplotchPainter::~vtkCUDASplotchPainter()
{
  this->PrepareDirectRenderBuffers(0, 0);
  delete this->Internal;
}
//-----------------------------------------------------------------------------
int device_binding(int mpi_rank)
{
  int local_rank = mpi_rank;
  int dev_count, use_dev_count, my_dev_id;
  char *str;

  if ((str = getenv ("MV2_COMM_WORLD_LOCAL_RANK")) != NULL)
  {
    local_rank = atoi (str);
    printf ("MV2_COMM_WORLD_LOCAL_RANK %s\n", str);
  }

  if ((str = getenv ("MPISPAWN_LOCAL_NPROCS")) != NULL)
  {
    //num_local_procs = atoi (str);
    printf ("MPISPAWN_LOCAL_NPROCS %s\n", str);
  }

  dev_count = vtkpiston::GetCudaDeviceCount();
  if ((str = getenv ("NUM_GPU_DEVICES")) != NULL)
  {
    use_dev_count = atoi (str);
    printf ("NUM_GPU_DEVICES %s\n", str);
  }
  else
  {
    use_dev_count = dev_count;
  }

  my_dev_id = (use_dev_count>0) ? (local_rank % use_dev_count) : 0;
  printf ("local rank = %d dev id = %d\n", local_rank, my_dev_id);
  return my_dev_id;
}
//-----------------------------------------------------------------------------
int vtkCUDASplotchPainter::InitCudaGL(vtkRenderWindow *rw, int rank, int &displayId)
{
  if (!vtkCUDASplotchPainter::CudaGLInitted)
  {
    int major = THRUST_MAJOR_VERSION;
    int minor = THRUST_MINOR_VERSION;
    std::cout << "Thrust v" << major << "." << minor << std::endl;
    //
    vtkOpenGLExtensionManager *em = vtkOpenGLExtensionManager::New();
    em->SetRenderWindow(rw);
    em->Update();
    if (!em->LoadSupportedExtension("GL_VERSION_1_5"))
    {
      std::cout << "WARNING: GL_VERSION_1_5 unsupported Can not use direct piston rendering" << endl;
      std::cout << em->GetExtensionsString() << std::endl;
      em->FastDelete();
      return 0;
    }
    em->FastDelete();
    if (displayId<0 || displayId>=vtkpiston::GetCudaDeviceCount()) {
      // try another method to get the device ID
      displayId = device_binding(rank);
    }
    vtkCUDASplotchPainter::CudaGLInitted = true;
    vtkpiston::CudaGLInit(displayId);
  }
  return 1;
}
//-----------------------------------------------------------------------------
void vtkCUDASplotchPainter::PrepareDirectRenderBuffers(int nPoints, int nCells)
{
  if (nPoints==this->Internal->BufferSize && nCells==this->Internal->CellCount) {
    return;
  }
  if (this->Internal->BufferSize != 0) {
    this->Internal->clearBuffers();
  }

  this->Internal->BufferSize = nPoints;
  this->Internal->CellCount  = nCells;
  if (this->Internal->BufferSize == 0) {
    return;
  }

  // Prep shared mem buffer between gl and cuda
  vtkgl::GenBuffers(NUM_INTEROP_BUFFERS, this->Internal->vboBuffers);

  // points 3*n float {x,y,z}
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[0]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*3*sizeof(float), 0,
    vtkgl::DYNAMIC_DRAW);

  // normals 3*n float {n0,n1,n2}
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[1]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*3*sizeof(float), 0,
    vtkgl::DYNAMIC_DRAW);

  // colors 4*n uchar {R,G,B,A} 
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[2]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*4*sizeof(unsigned char), 0,
    vtkgl::DYNAMIC_DRAW);

  // indexes 3*nCells int : triangles assumed {a,b,c} 
  vtkgl::BindBuffer(vtkgl::ELEMENT_ARRAY_BUFFER,
    this->Internal->vboBuffers[3]);
  vtkgl::BufferData(vtkgl::ELEMENT_ARRAY_BUFFER,
    this->Internal->CellCount*3*sizeof(int), 0,
    vtkgl::DYNAMIC_DRAW);

  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[0],
    this->Internal->vboBuffers[0]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[1],
    this->Internal->vboBuffers[1]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[2],
    this->Internal->vboBuffers[2]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[3],
    this->Internal->vboBuffers[3]);
}
//-----------------------------------------------------------------------------
void vtkCUDASplotchPainter::RenderOnGPU(vtkCamera *cam, vtkActor *act)
{
  vtkPistonDataObject *id = vtkPistonDataObject::SafeDownCast(this->DataSetToPiston->GetOutputDataObject(0));
  if (!id) {
    return;
  }

  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  int nPoints = vtkpiston::QueryNumVerts(id);
  int nCells = vtkpiston::QueryNumCells(id);
  this->PrepareDirectRenderBuffers(nPoints, nCells);

  // Transfer what is in tdo to buffer and render it directly on card
  bool hasNormals = false;
  bool hasColors = false;
  bool useindexbuffers = false;
/*
  double cameravec[3], origin[3];
  this->ComputeProjectionVector(cam, act, cameravec, origin);
  if (this->Direction>=0) {
//    vtkpiston::DepthSortPolygons(id, cameravec, this->Direction);
  }
*/

/*
  double scalarrange[2];
  this->ScalarsToColors->GetScalarRange(scalarrange);
  vtkpiston::CudaTransferToGL(
    id, this->Internal->DataObjectMTimeCache,
    this->Internal->vboResources, 
    this->ScalarsToColors->GetRGBAPointer(),
    scalarrange, 
    act->GetProperty()->GetOpacity(),
    hasNormals, hasColors, useindexbuffers);
*/
  glPushAttrib(GL_ENABLE_BIT);

  // Draw the result
  glEnableClientState(GL_VERTEX_ARRAY);
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[0]);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  if (hasNormals)
  {
    glEnableClientState(GL_NORMAL_ARRAY);
    vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[1]);
    glNormalPointer(GL_FLOAT, 0, 0);
  }

  if (hasColors)
  {
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
//    glEnable(GL_BLEND);
//    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
//    glEnable(GL_COLOR_MATERIAL);
    /*
    // because we have used premultiple_with_alpha
    // otherwise we'd use glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    glDepthFunc( GL_LEQUAL );
    glEnable( GL_DEPTH_TEST );
*/

    glEnable(GL_DEPTH_TEST);
    glEnableClientState(GL_COLOR_ARRAY);
    vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[2]);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, 0);

  }
  else {
    //    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    //    glEnable(GL_COLOR_MATERIAL);
    //    glDepthFunc( GL_LEQUAL );
    //    glEnable( GL_DEPTH_TEST );
    //    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //    glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glEnable(GL_BLEND);
    //glEnableClientState(GL_COLOR_ARRAY);
    //vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[2]);
    //glColorPointer(4, GL_FLOAT, 0, 0);
  }

  if (useindexbuffers) {
    //
    int vertsPer = vtkpiston::QueryVertsPer(id);
    vtkgl::BindBuffer(vtkgl::ELEMENT_ARRAY_BUFFER, this->Internal->vboBuffers[3]);
    switch (vertsPer) {
    case 4:
      glDrawElements(GL_QUADS, nCells*4, GL_UNSIGNED_INT, (GLvoid*)0);
      break;
    case 3:
      glDrawElements(GL_TRIANGLES, nCells*3, GL_UNSIGNED_INT, (GLvoid*)0);
      break;
    default:
      glDrawElements(GL_POINTS, nCells*1, GL_UNSIGNED_INT, (GLvoid*)0);
    }
  }
  else {
    int vertsPer = vtkpiston::QueryVertsPer(id);
    switch (vertsPer) {
    case 4:
      glDrawArrays(GL_QUADS, 0, nPoints);
      break;
    case 3:
      glDrawArrays(GL_TRIANGLES, 0, nPoints);
      break;
    default:
      glDrawArrays(GL_POINTS, 0, nPoints);
    }
  }

  glDisableClientState(GL_VERTEX_ARRAY);
  if (hasNormals) glDisableClientState(GL_NORMAL_ARRAY);
  if (hasColors) glDisableClientState(GL_COLOR_ARRAY);

  glPopAttrib();
  // Update object modified time
  this->Internal->DataObjectMTimeCache = id->GetMTime();

  timer->StopTimer();
  double rendertime = timer->GetElapsedTime();
  //  std::cout << setprecision(6) << "RenderTime : << " <<  rendertime << std::endl;
}
//-----------------------------------------------------------------------------
void vtkCUDASplotchPainter::PrepareForRendering(vtkRenderer* renderer, vtkActor* actor)
{
  vtkDataObject* input = this->GetInput();
  if (!input)
  {
    vtkErrorMacro("No input present.");
    return;
  }

  //
  // does a shallow copy of input to output
  //
  this->Superclass::PrepareForRendering(renderer, actor);

  // Now if we have composite data, we need to MapScalars for all leaves.
  if (input->IsA("vtkCompositeDataSet"))
  {
    throw std::runtime_error("Not supported");
  }
  else
  {
    this->DataSetToPiston->SetInputData(input);
    this->DataSetToPiston->SetOpacityArrayName(NULL);
//    this->DataSetToPiston->SetScalarArrayName(this->ScalarsToColors->GetArrayName());
    this->DataSetToPiston->Update();
  }
  //
  this->Camera = renderer->GetActiveCamera();
  this->Actor  = actor;
}
//-----------------------------------------------------------------------------

