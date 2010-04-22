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
#else
  ptypes = params.find<int>("ptypes",1);
  g_params =&params;
#endif  //ifdef CUDA they will be global vars

  get_colourmaps(params,amap);

  wallTimers.stop("setup");

  sceneMaker sMaker(params);
  string outfile;
  while (sMaker.getNextScene (particle_data, campos, lookat, sky, outfile))
    {
    long npart=particle_data.size();
    long npart_all=npart;
    mpiMgr.allreduce (npart_all,MPI_Manager::Sum);

    int res = params.find<int>("resolution",200);
    arr2<COLOUR> pic(res,res);

#ifndef CUDA
    host_rendering(master, params, npart_all, pic, particle_data,
                   campos, lookat, sky, amap);
#else
    cuda_rendering(res, pic, npart_all);
#endif

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
