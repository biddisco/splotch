/*
 * Copyright (c) 2004-2008
 *              Martin Reinecke (1), Klaus Dolag (1)
 *               (1) Max-Planck-Institute for Astrophysics
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include <iostream>
#include <cmath>
#include <algorithm>

#include "splotch/scenemaker.h"
#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"
#include "writer/writer.h"
#include "cxxsupport/walltimer.h"

#ifdef CUDA
#include "cuda/splotch_cuda2.h"
#endif

using namespace std;

int main (int argc, const char **argv)
  {
  wallTimers.start("full");
  wallTimers.start("setup");
  bool master = mpiMgr.master();
  module_startup ("splotch",argc,argv,2,"<parameter file>",master);
  paramfile params (argv[1],false);

#ifndef CUDA
  vector<particle_sim> particle_data; //raw data from file
  vec3 campos, lookat, sky;
  vector<COLOURMAP> amap;
#else //ifdef CUDA they will be global vars
  ptypes = params.find<int>("ptypes",1);
  g_params =&params;
  int nDevProc = 1; // number of gpus per process
#ifdef USE_MPI
  // We assume a geometry where each mpi-process uses one gpu
  int myID = mpiMgr.rank();
  int nDevNode = check_device(myID);     // number of gpus per node
  int mydevID = myID;
  if (myID >= nDevNode) mydevID = myID%nDevNode;
  if (nDevNode == 0 || mydevID >= nDevNode)
  {
      cout << "There isn't a gpu available for process = " << myID << endl;
      cout << "Configuration supported is 1 gpu for each mpi process" <<endl;
      mpiMgr.abort();
  }
  else printf("Rank %d: my device %d\n",myID, mydevID);
#else
  if (nDevNode == 0) exit(EXIT_FAILURE);
  int mydevID = 0;
#ifndef NO_WIN_THREAD
  nDevProc = g_params->find<int>("gpu_number",1);  // number of GPU per process
  if (nDevNode < nDevProc )
  {
      cout << "Number of GPUs available = " << nDevNode << " is lower than the number of GPUs required = " << nDevProc << endl;
      exit(EXIT_FAILURE);
  }
#endif // NO_WIN_THREAD
#endif // USE_MPI
  bool gpu_info = params.find<bool>("gpu_info",false);
  if (gpu_info) device_info(myID, mydevID);
#endif // CUDA 

  get_colourmaps(params,amap);

  wallTimers.stop("setup");

  sceneMaker sMaker(params);
  string outfile;
  while (sMaker.getNextScene (particle_data, campos, lookat, sky, outfile))
    {
    long npart_all = particle_data.size();
    mpiMgr.allreduce (npart_all,MPI_Manager::Sum);

    bool a_eq_e = params.find<bool>("a_eq_e",true);
    int res = params.find<int>("resolution",200);
    arr2<COLOUR> pic(res,res);

#ifndef CUDA
    host_processing(master, params, npart_all, particle_data,
                   campos, lookat, sky, amap);
 
// ----------- Rendering ---------------
    long nsplotch = particle_data.size();
    long nsplotch_all = nsplotch;
    mpiMgr.allreduce (nsplotch_all,MPI_Manager::Sum);
    if (master)
      cout << endl << "host: rendering (" << nsplotch_all << "/" << npart_all << ")..." << endl;

    float64 grayabsorb = params.find<float>("gray_absorption",0.2);

    wallTimers.start("render");
    render_new (particle_data,pic,a_eq_e,grayabsorb,false);
#else
    if (mydevID < nDevNode) cuda_rendering(mydevID, nDevProc, res, pic, npart_all);

  int xres = pic.size1(), yres=pic.size2();
  mpiMgr.allreduceRaw
    (reinterpret_cast<float *>(&pic[0][0]),3*xres*yres,MPI_Manager::Sum);

  exptable xexp(MAX_EXP);
//if (!nopostproc)
  if (mpiMgr.master() && a_eq_e)
      for (int ix=0;ix<xres;ix++)
        for (int iy=0;iy<yres;iy++)
          {
          pic[ix][iy].r=-xexp.expm1(pic[ix][iy].r);
          pic[ix][iy].g=-xexp.expm1(pic[ix][iy].g);
          pic[ix][iy].b=-xexp.expm1(pic[ix][iy].b);
          }
#endif

    wallTimers.stop("render");
 
// ------------ Writing -------------
    wallTimers.start("write");

    if (master && params.find<bool>("colorbar",false))
      {
      cout << endl << "creating color bar ..." << endl;
      add_colorbar(params,pic,amap);
      }

    if (master)
      cout << endl << "saving file ..." << endl;

    int pictype = params.find<int>("pictype",0);

    switch(pictype)
      {
      case 0:
        if (master) write_tga(params,pic,res,outfile);
        break;
      default:
        planck_fail("No valid image file type given ...");
        break;
      }

    wallTimers.stop("write");

    timeReport();
    }

#ifdef VS
  //Just to hold the screen to read the messages when debugging
  cout << endl << "Press any key to end..." ;
  getchar();
#endif
  }
