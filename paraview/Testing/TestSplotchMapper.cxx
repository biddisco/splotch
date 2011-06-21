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

#define PART_COUNT 1024
#define NUM_PARTITIONS 32
#define NUM_PARTITIONSS "32"
//----------------------------------------------------------------------------
#define ROWS 50
#define POINTS_PER_PROCESS ROWS*ROWS*ROWS
//----------------------------------------------------------------------------
std::string usage = "\n"\
"\t-D path to use for temp h5part file \n" \
"\t-F name to use for temp h5part file \n" \
"\t-C collective IO for MPI/HDF write (default independent IO) \n" \
"\t-N Generate exactly N points per process (no cube structure) \n" \
"\t-R Render. Displays points using vtk renderwindow \n" \
"\t-I Interactive (waits until user closes render window if -R selected) \n" \
"\t-X delete h5part file on completion \n" \
"\t\n" \
"\tWIN32 mpi example   : mpiexec -n 4 -localonly TestH5PartParallelWriter.exe -D D:\\ -F sph-test.h5 -X\n"\
"\tLinux mpich example : /project/csvis/biddisco/build/mpich2-1.0.7/bin/mpiexec -n 2 ./TestH5PartParallelWriter -D /project/csvis/hdf-scratch -C -N 1000000 -X\n" \
"\tHorus Slurm example : mpirun -prot -srun mpirun -e DISPLAY=horus11:0 -e MPI_IB_CARD_ORDER=0:0 -e MPI_IB_MTU=1024 -e MPI_IC_ORDER=ibv:vapi:udapl:itapi:TCP -srun -N 10 -n 20 bin/TestH5PartParallelWriter -D /project/csvis/hdf-scratch -C -N 100000 -R -I -X\n";

//----------------------------------------------------------------------------
#ifdef _WIN32
unsigned long int random_seed()
{
  LARGE_INTEGER lpPerformanceCount;
  QueryPerformanceCounter(&lpPerformanceCount);
  long int seed = lpPerformanceCount.LowPart + lpPerformanceCount.HighPart;
  srand(seed);
  return seed;
}
#else
unsigned long int random_seed()
{
  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;
  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
  } 
  else {
    if (fread(&seed,sizeof(seed),1,devrandom) == 1) {
      fclose(devrandom);
    } 
    else {
      gettimeofday(&tv,0);
      seed = tv.tv_sec + tv.tv_usec;
    }
  }
  srandom(seed);
  return seed;
}
#endif


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
void SpherePoints(int n, float radius, float X[]) {
  double x, y, z, w, t;
  for(int i=0; i< n; i++ ) {
    #ifdef WIN32
     double r1 = 0.5 + 0.5*double(rand())/RAND_MAX;
     double r2 = 0.5 + 0.5*double(rand())/RAND_MAX;
    #else
     double r1 = drand48();
     double r2 = drand48();
    #endif
    z = 2.0 * r1 - 1.0;
    t = 2.0 * M_PI * r2;
    w = radius * sqrt( 1 - z*z );
    x = w * cos( t );
    y = w * sin( t );
    X[3*i+0] = x;
    X[3*i+1] = y;
    X[3*i+2] = z*3.0*radius;
  }
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

  vtkTypeInt64 numPoints = POINTS_PER_PROCESS;
  double       rows      = ROWS;

  char *number = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-N", args->argc, args->argv, "DUMMY_ENV_VAR", "");
  if (std::string(number)!=std::string("")) {
    vtkstd::stringstream temp;
    temp << number;
    temp >> numPoints;
    rows = floor(pow(numPoints,1.0/3.0)+0.5);
    numPoints = static_cast<vtkTypeInt64>(pow(rows,3));
    if (myId==0) {
      std::cout << "Process Id : " << myId << " Requested Particles : " << numPoints << std::endl;
    }
  }

  //--------------------------------------------------------------
  // splotch particle read
  //--------------------------------------------------------------
  char *splotchparams = vtkTestUtilities::GetArgOrEnvOrDefault(
    "-s", args->argc, args->argv, "DUMMY_ENV_VAR", "");
  MPI_Manager::Instance = new MPI_Manager(false);
  bool master = MPI_Manager::GetInstance()->master();
  paramfile params (splotchparams,false);
  //
  std::vector<particle_sim> particle_data;
  vec3 campos, lookat, sky;
  sceneMaker sMaker(params);
  std::string outfile;
  sMaker.getNextScene (particle_data, campos, lookat, sky, outfile);

  numPoints = particle_data.size();

  //--------------------------------------------------------------
  // allocate scalar arrays
  //--------------------------------------------------------------
  vtkSmartPointer<vtkPolyData>   Sprites = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>      points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>    verts = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkDoubleArray> Values = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray>  Sizes = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> Bright = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray>       Ids = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray>     Ranks = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray>     Parts = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray>     Types = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray>    Active = vtkSmartPointer<vtkIntArray>::New();
  //
  points->SetNumberOfPoints(numPoints);
  //
  verts->Allocate(numPoints,numPoints);
  Sprites->SetPoints(points);
  Sprites->SetVerts(verts);
  //
  Sizes->SetNumberOfTuples(numPoints);
  Sizes->SetNumberOfComponents(1);
  Sizes->SetName("Radius");
  Sprites->GetPointData()->AddArray(Sizes);
  //
  Ids->SetNumberOfTuples(numPoints);
  Ids->SetNumberOfComponents(1);
  Ids->SetName("PointIds");
  Sprites->GetPointData()->AddArray(Ids);  
  //
  Ranks->SetNumberOfTuples(numPoints);
  Ranks->SetNumberOfComponents(1);
  Ranks->SetName("Rank");
  Sprites->GetPointData()->AddArray(Ranks);  
  //
  Parts->SetNumberOfTuples(numPoints);
  Parts->SetNumberOfComponents(1);
  Parts->SetName("Partition");
  Sprites->GetPointData()->AddArray(Parts);  
  //
  Types->SetNumberOfTuples(numPoints);
  Types->SetNumberOfComponents(1);
  Types->SetName("Type");
  Sprites->GetPointData()->AddArray(Types);  
  //
  Bright->SetNumberOfTuples(numPoints);
  Bright->SetNumberOfComponents(1);
  Bright->SetName("Intensity");
  Sprites->GetPointData()->AddArray(Bright);  
  //
  Values->SetNumberOfTuples(numPoints);
  Values->SetNumberOfComponents(1);
  Values->SetName("Scalars");
  Sprites->GetPointData()->AddArray(Values);  
  //
  Active->SetNumberOfTuples(numPoints);
  Active->SetNumberOfComponents(1);
  Active->SetName("Active");
  Sprites->GetPointData()->AddArray(Active);  
/*
  //--------------------------------------------------------------
  // Create default scalar arrays
  //--------------------------------------------------------------
  double radius  = 0.0001;
  radius  = 0.08;
  double spacing = radius*2.0;
  double offset  = myId*spacing*rows;
  const double a = 0.9;
  double index = 0;
  for (vtkTypeInt64 z=0; z<rows; z++) {
    for (vtkTypeInt64 y=0; y<rows; y++) {
      for (vtkTypeInt64 x=0; x<rows; x++) {
        vtkIdType Id = static_cast<vtkIdType>(z*rows*rows + y*rows + x);
        points->SetPoint(Id, x*spacing, y*spacing, z*spacing + offset);
        Sizes->SetValue(Id, radius);
        Ids->SetTuple1(Id, Id + myId*numPoints);
        Ranks->SetTuple1(Id, myId);
        Parts->SetTuple1(Id, myId);
        Bright->SetTuple1(Id, 2.0*sin(2.0*M_PI*index/(numPoints/5.0)));
        verts->InsertNextCell(1,&Id);
        index++;
      }
    }
  }
  SpherePoints(numPoints, (251.0+myId)*0.5/numProcs, vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0));
*/
  for (int i=0; i<numPoints; i++) {
    vtkIdType Id = static_cast<vtkIdType>(i);
    points->SetPoint(Id, particle_data[i].x, particle_data[i].y, particle_data[i].z);
    Sizes->SetValue(Id, particle_data[i].r);
    Ids->SetTuple1(Id, Id);
    Ranks->SetTuple1(Id, myId);
    Parts->SetTuple1(Id, myId);
    Bright->SetTuple1(Id, particle_data[i].I);
    Types->SetTuple1(Id, particle_data[i].type);
    Values->SetTuple1(Id, particle_data[i].e.r);
    Active->SetTuple1(Id, particle_data[i].active);
    verts->InsertNextCell(1,&Id);
  }

  //--------------------------------------------------------------
  // 
  //--------------------------------------------------------------

  bool collective = false;
  if (test->IsFlagSpecified("-C")) {
    if (myId==0) {
     std::cout << "Process Id : " << myId << " Collective IO requested" << std::endl;
    }
   collective = true;
  }


  //--------------------------------------------------------------
  // Create writer on all processes
  //--------------------------------------------------------------
  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetFileModeToWrite();
  writer->SetFileName(fullname);
  writer->SetInput(Sprites);
  writer->SetCollectiveIO(collective);
  writer->SetDisableInformationGather(1);
  writer->SetVectorsWithStridedWrite(0);

/*
  // Randomly give some processes zero points to improve test coverage
  random_seed();
  if (numProcs>1 && rand()%2==3) {
    numPoints = 0;
    Sprites = vtkSmartPointer<vtkPolyData>::New();
    writer->SetInput(Sprites);
  }
  else {
    writer->SetInputConnection(elev->GetOutputPort());
  }
*/

  controller->Barrier();
  if (myId==0) {
    std::cout << "Process Id : " << myId << " Generated N Points : " << numPoints << std::endl;
  }
  //
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  //if (numProcs>1 && myId==0) {
  //  char ch;  
  //  std::cin >> ch;
  //}

  //--------------------------------------------------------------
  // Write in parallel
  //--------------------------------------------------------------
  writer->SetTimeStep(0);
  writer->SetTimeValue(0.5);
  writer->Write();

  // 
  // make sure they have all finished writing before going on to the read part
  //
  timer->StopTimer();
  writer->CloseFile();
  controller->Barrier();

  // memory usage - Ids(int) Size(double) Elevation(float) Verts(double*3)
  double bytes = numPoints*(sizeof(int) + sizeof(double) + sizeof(float) + 3*sizeof(float));
  double MBytes = bytes/(1024*1024);
  double elapsed = timer->GetElapsedTime();
  std::cout << "Process Id : " << myId << " File Written in " << elapsed << " seconds" << std::endl;
  std::cout << "Process Id : " << myId << " IO-Speed " << MBytes/timer->GetElapsedTime() << " MB/s" << std::endl;
  //
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
  else
    {
    std::cout << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
    std::cout << "Process Id : " << myId << " Expected " << static_cast<vtkTypeInt64>(numPoints*numProcs) << std::endl;

    // Read the file we just wrote on N processes
    vtkSmartPointer<vtkH5PartReader> reader = vtkSmartPointer<vtkH5PartReader>::New();
    // we want to read all the particles on this node, so don't use MPI/Parallel
    reader->SetController(NULL);
    reader->SetFileName(fullname);
    reader->Update();
    vtkTypeInt64 ReadPoints = reader->GetOutput()->GetNumberOfPoints();
    std::cout << "Process Id : " << myId << " Read : " << ReadPoints << std::endl;
    std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * " << std::endl;
    //
    // Validate the point Ids to make sure nothing went wrong in the writing
    //
    vtkIntArray *pointIds = vtkIntArray::SafeDownCast(reader->GetOutput()->GetPointData()->GetArray("PointIds"));
    bool IdsGood = true;
    vtkTypeInt64 valid = 0;
    for (vtkTypeInt64 i=0; IdsGood && pointIds && i<ReadPoints; i++) {
      if (pointIds->GetValue(i)==i) {
        valid++;
      }
      else {
        IdsGood = false;
      }
    }
    if (IdsGood && ReadPoints==numPoints*numProcs) {
      *(args->retVal) = 0;
      std::cout << " " << std::endl;
      std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * "   << std::endl;
      std::cout << " All Points read back and Validated OK "               << std::endl;
      std::cout << " * * * * * * * * * * * * * * * * * * * * * * * * * \n" << std::endl;
      //
      if (myId==0) {
        unsigned long size = vtksys::SystemTools::FileLength(fullname);
        double filesize = size/(1024.0*1024.0);
        std::cout << "Process Id : " << myId << " Total IO/Disk-Speed " << filesize/elapsed << " MB/s" << std::endl;
      }
    }
    else {
      *(args->retVal) = 1;
      std::cout << " " << std::endl;
      std::cout << " # # # # # # # # # # # # # # # # # # # # # # # # # " << std::endl;
      std::cout << " FAIL "                                              << std::endl;
      std::cout << " Valid Ids "                                << valid << std::endl;
      std::cout << " # # # # # # # # # # # # # # # # # # # # # # # # # " << std::endl;
      std::cout << " " << std::endl;
    }

    bool doRender = false;
    if (test->IsFlagSpecified("-R")) {
     std::cout << "Process Id : " << myId << " Rendering" << std::endl;
     doRender = true;
    }

    if (doRender) {
      // Generate vertices from the points
      vtkSmartPointer<vtkMaskPoints> verts = vtkSmartPointer<vtkMaskPoints>::New();
      verts->SetGenerateVertices(1);
      verts->SetOnRatio(1);
      verts->SetMaximumNumberOfPoints(numPoints*numProcs);
      verts->SetInputConnection(reader->GetOutputPort());
      verts->Update();

      // Display something to see if we got a good output
      vtkSmartPointer<vtkPolyData> polys = vtkSmartPointer<vtkPolyData>::New();
      polys->ShallowCopy(verts->GetOutput());
      vtkSmartPointer<vtkDataArray> scalars = polys->GetPointData()->GetScalars();
      polys->GetPointData()->SetScalars(polys->GetPointData()->GetArray("Scalars"));
      polys->GetPointData()->AddArray(scalars);
      std::cout << "Process Id : " << myId << " Created vertices : " << polys->GetNumberOfPoints() << std::endl;
      //
#if 0
      vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#else
      vtkSmartPointer<vtkSplotchRaytraceMapper> mapper = vtkSmartPointer<vtkSplotchRaytraceMapper>::New();
//      mapper->SetRadiusScalars("Radius");
      mapper->SetRadiusScalars("Radius");
      mapper->SetIntensityScalars("Intensity");
      mapper->SetTypeScalars("Type");
      mapper->SetActiveScalars("Active");      
#endif



      vtkSmartPointer<vtkColorTransferFunction>    lut = vtkSmartPointer<vtkColorTransferFunction>::New();
      vtkSmartPointer<vtkActor>                  actor = vtkSmartPointer<vtkActor>::New();
      vtkSmartPointer<vtkRenderer>                 ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>       renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor>  iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      iren->SetRenderWindow(renWindow);
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize( 400, 400);
      mapper->SetInput(polys);

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

/*
      vtkSmartPointer<vtkPolyData> polys2 = vtkSmartPointer<vtkPolyData>::New();
      polys2->ShallowCopy(verts->GetOutput());
      polys2->GetPointData()->SetScalars(polys->GetPointData()->GetArray("Rank"));
      vtkSmartPointer<vtkPolyDataMapper>       mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
      vtkSmartPointer<vtkActor>                 actor2 = vtkSmartPointer<vtkActor>::New();
      mapper2->SetInput(polys2);
      mapper2->SetScalarRange(0,numProcs-1);
      actor2->SetMapper(mapper2);
      actor2->GetProperty()->SetPointSize(2);
      actor2->SetPosition(1.0, 0.0, 0.0);
      ren->AddActor(actor2);
*/
      std::cout << "Process Id : " << myId << " About to Render" << std::endl;
      renWindow->Render();

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

//  writer = NULL;

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

  return EXIT_SUCCESS; // !retVal;
}
//----------------------------------------------------------------------------
