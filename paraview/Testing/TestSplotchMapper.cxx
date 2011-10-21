// TestHDF5PartWriter -I -D D:/ -F test.h5

// windows mpi command line
// cd D:\cmakebuild\csviz\bin\relwithdebinfo
// mpiexec --localonly -n 4 TestH5PartParallelWriter -I -D D:\ -F sph-test.h5 -R -X -N 1000000

// horus slurm command line for auto allocation of nodes
// mpirun -prot -srun -N 6 -n 6 TestH5PartParallelWriter -R -D . -N 1000000

//
#include <vector>
#include <algorithm>
#include "splotch/scenemaker.h"
#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"

#ifdef _WIN32
  #include <windows.h>
#else 
  #include <sys/time.h>
#endif

#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkCamera.h"
#include "vtkPointSource.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkMPIController.h"
#include "vtkParallelFactory.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "Testing/Cxx/vtkTestUtilities.h"
#include "Testing/Cxx/vtkRegressionTestImage.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkDebugLeaks.h"
#include "vtkElevationFilter.h"
#include "vtkH5PartWriter.h"
#include "vtkH5PartReader.h"
#include "vtkMaskPoints.h"
#include "vtkProperty.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkTimerLog.h"
#include "vtkColorTransferFunction.h"
//
#include <vtksys/SystemTools.hxx>
#include <sstream>
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#define _USE_MATH_DEFINES
#include <math.h>
//
#include "vtkSplotchRaytraceMapper.h"
//----------------------------------------------------------------------------
std::string usage = "\n"\
"\t-D path to use for h5part file \n" \
"\t-F name to use for h5part file \n" \
"\t-R Render. Displays points using vtk renderwindow \n" \
"\t-I Interactive (waits until user closes render window if -R selected) \n" \
"\t\n";

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Just pick a tag which is available
static const int RMI_TAG=300; 
//----------------------------------------------------------------------------
struct ParallelArgs_tmp
{
  int* retVal;
  int    argc;
  char** argv;
};
//----------------------------------------------------------------------------
struct ParallelRMIArgs_tmp
{
  vtkMultiProcessController* Controller;
};
//----------------------------------------------------------------------------
// call back to set the iso surface value.
void SetStuffRMI(void *localArg, void* vtkNotUsed(remoteArg), 
                    int vtkNotUsed(remoteArgLen), int vtkNotUsed(id))
{ 
  ParallelRMIArgs_tmp* args = (ParallelRMIArgs_tmp*)localArg;
  vtkMultiProcessController* contrl = args->Controller;
}
//----------------------------------------------------------------------------
// This will be called by all processes
void MyMain( vtkMultiProcessController *controller, void *arg )
{
  // Obtain the id of the running process and the total
  // number of processes
  vtkTypeInt64 myId = controller->GetLocalProcessId();
  vtkTypeInt64 numProcs = controller->GetNumberOfProcesses();

  if (myId==0) {
    std::cout << usage.c_str() << std::endl;
  }
  controller->Barrier();

  //--------------------------------------------------------------
  // command line params : Setup testing utilities/args etc
  //--------------------------------------------------------------
  ParallelArgs_tmp* args = reinterpret_cast<ParallelArgs_tmp*>(arg);
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<args->argc; c++ ) {
    test->AddArgument(args->argv[c]);
  }
  // Get test filename etc
  char *filename = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-F", args->argc, args->argv, "DUMMY_ENV_VAR", "temp.h5");
  char* fullname = vtkTestUtilities::ExpandDataFileName(args->argc, args->argv, filename);
  if (myId==0) {
    std::cout << "Process Id : " << myId << " FileName : " << fullname << std::endl;
  }

  //--------------------------------------------------------------
  // processes 1-N doing nothing for now
  //--------------------------------------------------------------
  if (myId != 0)
    {
    // If I am not the root process
    ParallelRMIArgs_tmp args2;
    args2.Controller = controller;

    // We are not using any RMI's yet
    controller->AddRMI(SetStuffRMI, (void *)&args2, RMI_TAG);
    controller->ProcessRMIs();
    
    }
  //--------------------------------------------------------------
  // Read back all particles on process zero
  //--------------------------------------------------------------
  else {
    //
    // Read particles on N processes
    //
    vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
    reader->SetController(NULL);
    reader->SetFileName(fullname);
    reader->SetGenerateVertexCells(1);
    reader->Update();
    vtkTypeInt64 ReadPoints = reader->GetOutput()->GetNumberOfPoints();
    std::cout << "Process Id : " << myId << " Read : " << ReadPoints << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;

    bool doRender = false;
    if (test->IsFlagSpecified("-R")) {
     std::cout << "Process Id : " << myId << " Rendering" << std::endl;
     doRender = true;
    }

    if (doRender) {
      vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
      timer->StartTimer();
      //
      vtkSmartPointer<vtkSplotchRaytraceMapper> mapper = vtkSmartPointer<vtkSplotchRaytraceMapper>::New();
      mapper->SetNumberOfParticleTypes(2);
      mapper->SetTypeScalars("Type");
      mapper->SetActiveScalars("Active");      
      // type 0
      mapper->SetRadiusScalars(0,"Radius");
      mapper->SetIntensityScalars(0,"Intensity");
      mapper->SetBrightness(0,10.5);
      mapper->SetLogIntensity(0,1);
      mapper->SetLogColour(0,0);
      mapper->SetTypeActive(0,1);
      // type 1
      mapper->SetRadiusScalars(1,"Radius");
      mapper->SetIntensityScalars(1,"Intensity");
      mapper->SetBrightness(1,10.5);
      mapper->SetLogIntensity(1,1);
      mapper->SetLogColour(1,0);
      mapper->SetTypeActive(1,1);
      //
      mapper->SetColorModeToMapScalars();
      mapper->SetScalarModeToUsePointFieldData();
      mapper->SelectColorArray("Scalars");
      //
      vtkSmartPointer<vtkColorTransferFunction>    lut = vtkSmartPointer<vtkColorTransferFunction>::New();
      vtkSmartPointer<vtkActor>                  actor = vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkRenderer>                 ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>       renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor>  iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize( 400, 400);
      mapper->SetInputConnection(reader->GetOutputPort());

      lut->AddRGBPoint(0.0, 0.0, 0.0, 1.0);
      lut->AddRGBPoint(0.5, 0.5, 1.0, 0.5);
      lut->AddRGBPoint(1.0, 1.0, 0.0, 0.0);

      lut->Build();
      mapper->SetLookupTable(lut);
      mapper->SetScalarRange(0,1);

      actor->SetMapper(mapper);
      ren->AddActor(actor);
      renWindow->AddRenderer(ren);

      vtkCamera *cam = ren->GetActiveCamera();
      cam->SetPosition(3244.4, 25289.3, 4764.7);
      cam->SetFocalPoint(-2000, 5289.3, 4764.7);
      cam->SetViewUp(0,0,-1.0);
      ren->ResetCameraClippingRange();

      std::cout << "Process Id : " << myId << " About to Render" << std::endl;
      renWindow->Render();

      // 
      timer->StopTimer();
      controller->Barrier();

      *(args->retVal) = 
        vtkRegressionTester::Test(args->argc, args->argv, renWindow, 10);

      if ( *(args->retVal) == vtkRegressionTester::DO_INTERACTOR)
        {
        iren->Start();
        }
      std::cout << "Process Id : " << myId << " Rendered" << std::endl;
    }

    // Tell the other processors to stop processing RMIs.
    for (int i = 1; i < numProcs; ++i)
    {
      controller->TriggerRMI(i, vtkMultiProcessController::BREAK_RMI_TAG); 
    }
  }

  if (myId==0 && test->IsFlagSpecified("-X")) {
   std::cout << "Process Id : " << myId << " About to Delete file" << std::endl;
//   vtksys::SystemTools::RemoveFile(fullname);
  }
  delete []fullname;
  delete []filename;
}
//----------------------------------------------------------------------------
int main (int argc, char* argv[])
{
  // Check if MPI is already initialized because MPI_Manager in sploth calls it
  int initialized = false;
  if (MPI_Initialized(&initialized)==MPI_SUCCESS && !initialized) {
    MPI_Init(&argc, &argv);
  }

  // Note that this will create a vtkMPIController if MPI
  // is configured, vtkThreadedController otherwise.
  vtkMPIController* controller = vtkMPIController::New();

  controller->Initialize(&argc, &argv, 1);

  vtkParallelFactory* pf = vtkParallelFactory::New();
  vtkObjectFactory::RegisterFactory(pf);
  pf->Delete();
 
  // Added for regression test.
  // ----------------------------------------------
  int retVal = 1;
  ParallelArgs_tmp args;
  args.retVal = &retVal;
  args.argc = argc;
  args.argv = argv;
  // ----------------------------------------------

  controller->SetSingleMethod(MyMain, &args);
  controller->SingleMethodExecute();

  controller->Barrier();
  controller->Finalize();
  controller->Delete();

  return !retVal;
}
//----------------------------------------------------------------------------
