
#include <iostream>
#include <cmath>
#include <algorithm>

#ifdef SPLVISIVO
#include "optionssetter.h"
#endif

#include "splotch/scenemaker.h"
#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"
#include "cxxsupport/walltimer.h"
#include "cxxsupport/ls_image.h"
#include "cxxsupport/announce.h"

#ifdef OPENCL
#include "opencl/splotch_cuda2.h"
#endif

using namespace std;
#ifdef SPLVISIVO
int splotchMain (VisIVOServerOptions opt)
#else
int main (int argc, const char **argv)
#endif
  {
  tstack_push("Splotch total time");
  tstack_push("Setup");
  bool master = mpiMgr.master();
#ifdef SPLVISIVO
//  module_startup ("splotch",argc,argv,2,"<parameter file>",master);
  if(opt.splotchpar.empty())
  {	
	std::cout<<"usage: --splotch <parameter file>"<<std::endl;
	return -1;
   }
  paramfile params (opt.splotchpar.c_str(),false);
#else  
  module_startup ("splotch",argc,argv,2,"<parameter file>",master);
  paramfile params (argv[1],false);
#endif
  

#ifndef OPENCL
  vector<particle_sim> particle_data; //raw data from file
  vec3 campos, lookat, sky;
  vector<COLOURMAP> amap;
#else //ifdef CUDA they will be global vars
  ptypes = params.find<int>("ptypes",1);
  g_params =&params;

  int myID = mpiMgr.rank();
  int nDevNode = check_device(myID);     // number of GPUs available per node
  if (nDevNode < 1)   mpiMgr.abort();

  int nDevProc = params.find<int>("gpu_number",1);  // number of GPU required per process
  int mydevID = 0;
  // We assume a geometry where
  // a) either each process uses only one gpu
  if (nDevProc == 1)
    {
    mydevID = myID;
    if (mydevID >= nDevNode) mydevID = myID%nDevNode;
    if ( mydevID >= nDevNode)
      {
      cout << "There isn't a gpu available for process = " << myID << endl;
      cout << "Configuration supported is 1 gpu for each mpi process" <<endl;
      mpiMgr.abort();
      }
    }
  // b) or processes run on different nodes and use a number of GPUs >= 1 and <= nDevNode
  else if (nDevNode < nDevProc)
    {
    cout << "Number of GPUs available = " << nDevNode << " is lower than the number of GPUs required = " << nDevProc << endl;
    mpiMgr.abort();
    }

  bool gpu_info = params.find<bool>("gpu_info",false);
  if (gpu_info) print_device_info(myID, mydevID);
#endif // CUDA

#ifdef SPLVISIVO
  get_colourmaps(params,amap,opt); 
#else
  get_colourmaps(params,amap); 
#endif
  tstack_pop("Setup");
  string outfile;

#ifdef SPLVISIVO
  sceneMaker sMaker(params,opt);  
  while (sMaker.getNextScene (particle_data, campos, lookat, sky, outfile,opt))
#else
  sceneMaker sMaker(params);  
  while (sMaker.getNextScene (particle_data, campos, lookat, sky, outfile))
#endif
    {
    bool a_eq_e = params.find<bool>("a_eq_e",true);
    int xres = params.find<int>("xres",800),
        yres = params.find<int>("yres",xres);
    arr2<COLOUR> pic(xres,yres);

#ifndef OPENCL
    if(particle_data.size()>0)
#ifdef SPLVISIVO
    host_rendering(params, particle_data, pic, campos, lookat, sky, amap,opt);
#else
      host_rendering(params, particle_data, pic, campos, lookat, sky, amap); 
#endif
#else
    cuda_rendering(mydevID, nDevProc, pic);
#endif

    tstack_push("Post-processing");
    mpiMgr.allreduceRaw
      (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres,MPI_Manager::Sum);

    exptable<float32> xexp(-20.0);
    if (mpiMgr.master() && a_eq_e)
      for (int ix=0;ix<xres;ix++)
        for (int iy=0;iy<yres;iy++)
          {
          pic[ix][iy].r=-xexp.expm1(pic[ix][iy].r);
          pic[ix][iy].g=-xexp.expm1(pic[ix][iy].g);
          pic[ix][iy].b=-xexp.expm1(pic[ix][iy].b);
          }

    tstack_replace("Post-processing","Output");

    if (master && params.find<bool>("colorbar",false))
      {
      cout << endl << "creating color bar ..." << endl;
      add_colorbar(params,pic,amap);
      }

    if(!params.find<bool>("AnalyzeSimulationOnly"))
      {
      if (master)
        {
        cout << endl << "saving file ..." << endl;

        LS_Image img(pic.size1(),pic.size2());

        for (tsize i=0; i<pic.size1(); ++i)
          for (tsize j=0; j<pic.size2(); ++j)
            img.put_pixel(i,j,Colour(pic[i][j].r,pic[i][j].g,pic[i][j].b));
        int pictype = params.find<int>("pictype",0);
        switch(pictype)
          {
          case 0:
            img.write_TGA(outfile);
            break;
          case 1:
            planck_fail("ASCII PPM no longer supported");
            break;
          case 2:
            img.write_PPM(outfile);
            break;
          case 3:
            img.write_TGA_rle(outfile);
            break;
          default:
            planck_fail("No valid image file type given ...");
            break;
          }
        }
      }

    tstack_pop("Output");

#ifdef OPENCL
    cuda_timeReport(params);
#else
    timeReport();
#endif
    }

#ifdef VS
  //Just to hold the screen to read the messages when debugging
  cout << endl << "Press any key to end..." ;
  getchar();
#endif
  }