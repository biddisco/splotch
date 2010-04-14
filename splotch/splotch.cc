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
  vector<COLOURMAP> amap,emap;
#else
  ptypes = params.find<int>("ptypes",1);
  g_params =&params;
#endif  //ifdef CUDA they will be global vars

  get_colourmaps(params,amap,emap);

  wallTimers.stop("setup");

  sceneMaker sMaker(params);
  string outfile;
  while (sMaker.getNextScene (particle_data, campos, lookat, sky, outfile))
    {
    long npart=particle_data.size();
    long npart_all=npart;
    mpiMgr.allreduce (npart_all,MPI_Manager::Sum);

#ifdef CUDA
    int res;
    arr2<COLOUR> pic;
    render_cuda(params,res,pic);
#else
    wallTimers.start("range");
    if (master)
      cout << endl << "ranging values (" << npart_all << ") ..." << endl;
    particle_normalize(params,particle_data,true); ///does log calculations and clamps data
    wallTimers.stop("range");

    wallTimers.start("transform");
    if (master)
      cout << endl << "applying geometry (" << npart_all << ") ..." << endl;
    particle_project(params, particle_data, campos, lookat, sky);
    wallTimers.stop("transform");

    wallTimers.start("sort");
    if (master)
      (mpiMgr.num_ranks()>1) ?
        cout << endl << "applying local sort ..." << endl :
        cout << endl << "applying sort (" << npart << ") ..." << endl;
    int sort_type = params.find<int>("sort_type",1);
    particle_sort(particle_data,sort_type,true);
    wallTimers.stop("sort");

    wallTimers.start("coloring");
    if (master)
      cout << endl << "calculating colors (" << npart_all << ") ..." << endl;
    particle_colorize(params, particle_data, amap, emap);
    wallTimers.stop("coloring");

    int res = params.find<int>("resolution",200);
    long nsplotch=particle_data.size();
    long nsplotch_all=nsplotch;
    mpiMgr.allreduce (nsplotch_all,MPI_Manager::Sum);
    if (master)
      cout << endl << "rendering (" << nsplotch_all << "/" << npart_all << ")..." << endl;
    arr2<COLOUR> pic(res,res);
    float64 grayabsorb = params.find<float>("gray_absorption",0.2);
    bool a_eq_e = params.find<bool>("a_eq_e",true);
    wallTimers.start("render");
    bool new_renderer = params.find<bool>("new_renderer",false);
    new_renderer ? render_new    (particle_data,pic,a_eq_e,grayabsorb,false)
                 : render_classic(particle_data,pic,a_eq_e,grayabsorb,false);
    wallTimers.stop("render");
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
