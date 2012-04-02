// DGraham Data
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.plymouth.cfg -f H:/ParticleData/DGraham/DAT_t_6_2.ascii 
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.plymouth.cfg -f H:/ParticleData/DGraham/Graham_wave_data_New.ascii 
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.plymouth.cfg -f H:/ParticleData/DGraham/DATA_2.ascii 
//
// EDF Data
// Ascii2H5Part -c C:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.edf.cfg -f C:/data/ParticleData/EDF/test_case_output95.dat
// Ascii2H5Part -c D:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.edf.cfg -f H:\ParticleData\EDF\test_case\test_case_output1.dat
//
// UMan Data
// Ascii2H5Part -c C:/cmakebuild/plugins/bin/ASCII2H5Part.manchester.cfg  -f C:/Data/ParticleData/UMan/Case6/PART_0001 -a C:/Data/ParticleData/UMan/Case6/PART_0000_point_attribs.txt   
//
// ECN Data
// Ascii2H5Part -c C:/code/CSCS/pv-meshless/CSCS/vtkH5Part/Tools/Ascii2H5Part.ecn.cfg -f C:/Data/ParticleData/ECN/BGL_SPhere0POSX.dat
// Ascii2H5Part -c D:/Code/CSCS/csviz/Shared/vtkCSCS/vtkH5Part/Tools/ASCII2H5Part.ecn.cfg -f D:/data/ecn/MillenISOP0POSX.dat 
// 
// Ascii2H5Part -c C:/Code/CSCS/cscs_plugins/vtkCSCSMeshless/h5part/tools/ASCII2H5Part.ens-lyon.cfg -d C:/data/AsciiData-ens-lyon -f dust_sph_030.out
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib> 
#include <vtksys/SystemTools.hxx>
#include <vtksys/Glob.hxx>
#include <vtksys/RegularExpression.hxx>
#include <vtksys/Process.h>
#include <sstream>

#include "vtkSmartPointer.h"
#include "vtkTesting.h"
#include "Testing/Cxx/vtkTestUtilities.h"
#include "vtkAppendPolyData.h"
#include "vtkPointData.h"
#include "vtkH5PartWriter.h"
#include "vtkFloatArray.h"
//
#include "scenemaker.h"
#include "reader.h"
#include "string_utils.h"
//
typedef vtkstd::vector<vtkstd::string> stringlist;
typedef vtkstd::pair< vtkstd::string, vtkstd::string > maptype;
//
vtkstd::map<vtkstd::string, vtkstd::string> RegExMap;
stringlist patterns, regexmatches, Tstrings, Bstrings, Vstrings;
vtkstd::string regex;
//
int Tindex = -1, Tnum = 0;
int Bindex = -1, Bnum = 0;
int Vindex = -1, Vnum = 0;
vtkstd::string timeform, blockform, varform, filepattern;

//----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  vtkstd::cout << "Usage : ASCII2H5Part "  
    << "-c full/path/to/.cfg "
    << "-a full/path/to/IdFile for optional point Id flagging "
    << "-t full/path/to/timefile.txt " 
    << "-f full/path/to/input/file " << vtkstd::endl;
  vtkSmartPointer<vtkTesting> test = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    test->AddArgument(argv[c]);
  }

  std::string asciipath, asciiname, errormsg, asciifile, hdf5file, ascii2h5config, ascii2h5IdFile;

  // input file
  char *filename = vtkTestUtilities::GetArgOrEnvOrDefault("-f", argc, argv, "DUMMY_ENV_VAR", "");
  asciipath = vtksys::SystemTools::GetFilenamePath(filename);
  asciiname = vtksys::SystemTools::GetFilenameName(filename);
  asciifile = asciipath + "/" + asciiname;
  delete []filename;
  if (!vtksys::SystemTools::FileExists(asciifile.c_str()))
  {
    std::cout << "Can't find input file " << asciifile.c_str() << "\n";
    return 0;
  }
  vtkstd::cout << "Input file found     : "  << asciifile.c_str() << "\n";;

  // Point Id file
  char *idf_name = vtkTestUtilities::GetArgOrEnvOrDefault("-a", argc, argv, "DUMMY_ENV_VAR", "");
  vtkstd::map<long int, double> idMap;
  if (vtksys::SystemTools::FileExists(idf_name)) {
    ascii2h5IdFile = idf_name;
    vtkstd::cout << "Using Id file        : "  << ascii2h5IdFile.c_str() << "\n";
    vtkstd::ifstream idfile(idf_name);
    long int low, high;
    double flag;
    while (idfile.good()) {
      idfile >> low >> high >> flag;
      idMap.insert( vtkstd::pair<long int,double>(high, flag));
    }
  }
  delete []idf_name;

  // time override 
  bool   OverrideTime = false;
  double OverrideTimeStep = 0.0;
  char *time_steps = vtkTestUtilities::GetArgOrEnvOrDefault("-t", argc, argv, "DUMMY_ENV_VAR", "");
  std::string timeoverride = time_steps;
  if (timeoverride.size()>0) {
    OverrideTime = true;
    OverrideTimeStep = atof(timeoverride.c_str());
    vtkstd::cout << "Time override set to : " << OverrideTimeStep << " per step \n";;
  }
  delete []time_steps;

  //
  // generate new h5part file name
  //
  hdf5file  = vtksys::SystemTools::GetFilenameWithoutExtension(asciifile);
  hdf5file  = asciipath + "/" + hdf5file + ".h5part";
  vtkstd::cout << "Output HDF5 filename : "  << hdf5file.c_str() << "\n";;


  std::vector<particle_sim> particle_data;
  double boxsize;
  int interpol_mode = 0;

  paramfile params(asciifile,true);
  params.find<bool>("AnalyzeSimulationOnly", false);
/*
  params.find<int>("ptypes", 2);
  params.find<int>("ptype0", 0);
  params.find<int>("ptype1", 4);
*/
  //
  MPI_Manager singleton(true);
  if (MPI_Manager::GetInstance()->master()) {}
  double dummy;
  
  int fidx = 0;
  int simtype = params.find<int>("simtype");
  int spacing = params.find<double>("snapshot_spacing",1);
  int snr1 = int(fidx/spacing)*spacing, snr2=snr1+spacing;
  double frac=(fidx-snr1)/spacing;

// only used if interpol_mode>0
    std::vector<particle_sim> p1,p2;
    std::vector<MyIDType> id1,id2;
    std::vector<uint32> idx1,idx2;
    int snr1_now = -1;
    int snr2_now = -1;
    double time1,time2;
// only used if interpol_mode>1
    std::vector<vec3f> vel1,vel2;

  switch (simtype)
    {
    case 0:
      bin_reader_tab(params,particle_data);
      break;
    case 1:
      bin_reader_block(params,particle_data);
      break;
    case 2:
      if (interpol_mode>0) // Here only the two data sets are prepared, interpolation will be done later
        {
        cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fraction: " << frac << endl;
        cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl;
        if (snr1==snr2_now)
          {
          cout << " old2 = new1!" << endl;
          p1.swap(p2);
          id1.swap(id2);
          idx1.swap(idx2);
          vel1.swap(vel2);

          snr1_now = snr1;
          time1 = time2;
          }
        if (snr1_now!=snr1)
          {
          cout << " reading new1 " << snr1 << endl;
          gadget_reader(params,interpol_mode,p1,id1,vel1,snr1,time1,boxsize);
	  wallTimers.stop("read");
	  wallTimers.start("buildindex");
          buildIndex(id1.begin(),id1.end(),idx1);
	  wallTimers.stop("buildindex");
	  wallTimers.start("read");
          snr1_now = snr1;
          }
        if (snr2_now!=snr2)
          {
          cout << " reading new2 " << snr2 << endl;
          gadget_reader(params,interpol_mode,p2,id2,vel2,snr2,time2,boxsize);
	  wallTimers.stop("read");
	  wallTimers.start("buildindex");
          buildIndex(id2.begin(),id2.end(),idx2);
	  wallTimers.stop("buildindex");
	  wallTimers.start("read");
          snr2_now = snr2;
          }
        }
      else
        {
        double dummy;
        gadget_reader(params,interpol_mode,particle_data,id1,vel1,0,dummy,boxsize);
        }
      break;
    case 3:
#if 0
      enzo_reader(params,particle_data);
#else
      planck_fail("Enzo reader not available in this version!");
#endif
      break;
    case 4:
      {
      double dummy;
      gadget_millenium_reader(params,particle_data,0,&dummy);
      break;
      }
    case 5:
#if defined(USE_MPIIO)
      {
      float maxr, minr;
      bin_reader_block_mpi(params,particle_data, &maxr, &minr, MPI_Manager::GetInstance()->rank(), MPI_Manager::GetInstance()->num_ranks());
      }
#else
      planck_fail("mpi reader not available in non MPI compiled version!");
#endif
      break;
    case 6:
      mesh_reader(params,particle_data);
      break;
#ifdef HDF5
    case 7:
      hdf5_reader(params,particle_data);
      break;
    case 8:
      // GADGET HDF5 READER (this block was initially copied from case 2)
      //
      if (interpol_mode>0) // Here only the two data sets are prepared, interpolation will be done later
        {
        cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fraction: " << frac << endl;
        cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl;
        if (snr1==snr2_now)
          {
          cout << " old2 = new1!" << endl;
          p1.swap(p2);
          id1.swap(id2);
          idx1.swap(idx2);
          vel1.swap(vel2);

          snr1_now = snr1;
          time1 = time2;
          }
        if (snr1_now!=snr1)
          {
          cout << " reading new1 " << snr1 << endl;
          gadget_hdf5_reader(params,interpol_mode,p1,id1,vel1,snr1,time1,boxsize);
          tstack_replace("Input","Particle index generation");
          buildIndex(id1.begin(),id1.end(),idx1);
          tstack_replace("Particle index generation","Input");
          snr1_now = snr1;
          }
        if (snr2_now!=snr2)
          {
          cout << " reading new2 " << snr2 << endl;
          gadget_hdf5_reader(params,interpol_mode,p2,id2,vel2,snr2,time2,boxsize);
          tstack_replace("Input","Particle index generation");
          buildIndex(id2.begin(),id2.end(),idx2);
          tstack_replace("Particle index generation","Input");
          snr2_now = snr2;
          }
        }
      else
        {
        double dummy;
        gadget_hdf5_reader(params,interpol_mode,particle_data,id1,vel1,0,dummy,boxsize);
        }
      break;
#endif
#ifdef SPLVISIVO
    case 10:
      if(!visivo_reader(params,particle_data,opt))
	planck_fail("Invalid read data ...");
	break;
#endif
    default:
      planck_fail("No valid file type given ...");
      break;
    }

  wallTimers.stop("read");

  //
  //
  //
  vtkIdType N = particle_data.size();
  //
  // Create points
  //
  vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(N);
  poly->SetPoints(points);
  // 
  vtkSmartPointer<vtkFloatArray> intensity = vtkSmartPointer<vtkFloatArray>::New();
  intensity->SetNumberOfComponents(1);
  intensity->SetNumberOfTuples(N);
  intensity->SetName("Intensity");
  poly->GetPointData()->AddArray(intensity);
  //
  vtkSmartPointer<vtkFloatArray> radius = vtkSmartPointer<vtkFloatArray>::New();
  radius->SetNumberOfComponents(1);
  radius->SetNumberOfTuples(N);
  radius->SetName("Radius");
  poly->GetPointData()->AddArray(radius);  
  //
  vtkSmartPointer<vtkUnsignedCharArray> color = vtkSmartPointer<vtkUnsignedCharArray>::New();
  color->SetNumberOfComponents(3);
  color->SetNumberOfTuples(N);
  color->SetName("Color");
  poly->GetPointData()->AddArray(color);  
  //
  vtkSmartPointer<vtkUnsignedCharArray> ptype = vtkSmartPointer<vtkUnsignedCharArray>::New();
  ptype->SetNumberOfComponents(1);
  ptype->SetNumberOfTuples(N);
  ptype->SetName("Type");
  poly->GetPointData()->AddArray(ptype);  
  //
  vtkSmartPointer<vtkUnsignedCharArray> active = vtkSmartPointer<vtkUnsignedCharArray>::New();
  active->SetNumberOfComponents(1);
  active->SetNumberOfTuples(N);
  active->SetName("Active");
  poly->GetPointData()->AddArray(active);  
  //
  for (vtkIdType i=0; i<N; i++) {
    particle_sim &pref = particle_data[i];
    points->SetPoint(i, &pref.x);
    intensity->SetValue(i,pref.I);
    radius->SetValue(i,pref.r);
    ptype->SetValue(i,pref.type);
    unsigned char RGB[3] = {
      static_cast<unsigned char>(0.5+(pref.e.r*255.0)), 
      static_cast<unsigned char>(0.5+(pref.e.g*255.0)), 
      static_cast<unsigned char>(0.5+(pref.e.b*255.0))};
    color->SetTupleValue(i,RGB);
    active->SetValue(i,pref.active);
  }

  vtkSmartPointer<vtkH5PartWriter> writer = vtkSmartPointer<vtkH5PartWriter>::New();
  writer->SetInput(poly);
  writer->SetFileName(hdf5file.c_str());
  writer->Write();

  std::cout << "Closing hdf file " << std::endl;
  std::cout << "Written (hopefully) : " << hdf5file.c_str() << std::endl;
  writer->CloseFile();
  std::cout << "Done" << std::endl;
}

