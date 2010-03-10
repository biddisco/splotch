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
#include <fstream>
#include <algorithm>

#include "splotch/splotchutils.h"
#include "writer/writer.h"
#include "reader/reader.h"
#include "cxxsupport/walltimer.h"

#ifdef CUDA
#include "cuda/splotch_cuda2.h"
#endif

using namespace std;

int main (int argc, const char **argv)
  {
  wallTimer.start("full");
  wallTimer.start("setup");
  bool master = mpiMgr.master();
  module_startup ("splotch",argc,argv,2,"<parameter file>",master);

#ifdef INTERPOLATE
#ifndef GEOMETRY_FILE
#error Splotch: interpolation without geometry file makes no sense!
#endif
#endif

// -----------------------------------
// ----------- Needed Data -----------
// -----------------------------------

  //paramfile params (argv[1],master);
  paramfile params (argv[1],false);
#ifndef CUDA_THREADS
  vector<particle_sim> particle_data; ///row data from file
  vec3 campos, lookat, sky; ///A 3D vector class, designed for high efficiency.
  vector<COLOURMAP> amap,emap;
  int ptypes = params.find<int>("ptypes",1); ///each particle type has a color map
#else //if CUDA_THREADS defined
  ptypes = params.find<int>("ptypes",1); ///each particle type has a color map
  g_params =&params;
#endif  //if def CUDA_THREADS they will be a global vars

#ifdef INTERPOLATE
  vector<particle_sim> particle_data1,particle_data2;
  int snr_start = params.find<int>("snap_start",10);
  int snr1=snr_start,snr2=snr_start+1,snr1_now=-1,snr2_now=-1;
  double time1,time2;
#endif
  double time;

// ----------------------------------------------
// ----------- Loading Color Maps ---------------
// ----------------------------------------------

  amap.resize(ptypes);

  if (master)
    cout << "building color maps (" << ptypes << ")..." << endl;
  for (int itype=0;itype<ptypes;itype++)
    {
    if (params.find<bool>("color_is_vector"+dataToString(itype),false))
      {
      if (master)
        cout << " color of ptype " << itype << " is vector, so no colormap to load ..." << endl;
      }
    else
      {
      ifstream infile (params.find<string>("palette"+dataToString(itype)).c_str());
      planck_assert (infile,"could not open palette file  <" +
        params.find<string>("palette"+dataToString(itype)) + ">");
      string dummy;
      int nColours;
      infile >> dummy >> dummy >> nColours;
      if (master)
        cout << " loading " << nColours << " entries of color table of ptype " << itype << endl;
      double step = 1./(nColours-1);
      for (int i=0; i<nColours; i++)
        {
        float rrr,ggg,bbb;
        infile >> rrr >> ggg >> bbb;
        amap[itype].addVal(i*step,COLOUR(rrr/255,ggg/255,bbb/255));
        }
      }
    }
  emap=amap;

  wallTimer.stop("setup");
  wallTimer.start("read");

#ifdef GEOMETRY_FILE

// -------------------------------------
// -- Looping over a flight path -------
// -------------------------------------
  vector<particle_sim> p_orig;

  ifstream inp(params.find<string>("geometry_file").c_str());
  int linecount=0,ninterpol=0,nextfile=0;
  int geometry_skip = params.find<int>("geometry_start",0);
  int geometry_incr = params.find<int>("geometry_incr",1);

  string line;
  for(int i=0; i<geometry_skip; i++, linecount++)
    {
    getline(inp, line);
#ifdef INTERPOLATE
    sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i",
           &campos.x,&campos.y,&campos.z,
           &lookat.x,&lookat.y,&lookat.z,
           &sky.x,&sky.y,&sky.z,&ninterpol);
    if (linecount==nextfile)
      {
      nextfile=linecount+ninterpol;
      snr1=snr2;
      snr2++;
      }
#endif
    }

  while (getline(inp, line))
    {
    sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i",
           &campos.x,&campos.y,&campos.z,
           &lookat.x,&lookat.y,&lookat.z,
           &sky.x,&sky.y,&sky.z,&ninterpol);
    if (master)
      {
      cout << endl << "Next entry <" << linecount << "> in geometry file ..." << endl;
      cout << " Camera:    " << campos << endl;
      cout << " Lookat:    " << lookat << endl;
      cout << " Sky:       " << sky << endl;
#ifdef INTERPOLATE
      cout << " ninterpol: " << ninterpol << endl;
#endif
      }
#ifdef INTERPOLATE
    if(linecount == 0 && nextfile == 0)
      nextfile=linecount+ninterpol;
#endif
#else
    campos.x=params.find<double>("camera_x");
    campos.y=params.find<double>("camera_y");
    campos.z=params.find<double>("camera_z");
    lookat.x=params.find<double>("lookat_x");
    lookat.y=params.find<double>("lookat_y");
    lookat.z=params.find<double>("lookat_z");
    sky.x=params.find<double>("sky_x",0);
    sky.y=params.find<double>("sky_y",0);
    sky.z=params.find<double>("sky_z",0);
#endif


// -----------------------------------
// ----------- Reading ---------------
// -----------------------------------

#if defined(GEOMETRY_FILE) && !defined(INTERPOLATE)
    if (linecount==geometry_skip) // read only once if no interpolation is chosen
      {
#endif
      if (master)
        cout << endl << "reading data ..." << endl;
      int simtype = params.find<int>("simtype"); // 2:Gadget2
      float maxr, minr;
#ifdef INTERPOLATE
      double frac=(linecount-(nextfile-ninterpol))/double(ninterpol);
#endif
      switch (simtype)
        {
        case 0:
          bin_reader_tab(params,particle_data, maxr, minr);
          break;
        case 1:
          bin_reader_block(params,particle_data, maxr, minr);
          break;
        case 2:
#ifdef INTERPOLATE // Here only the two datasets are prepared, interpolation will be done later
          cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fac: " << frac << endl;
          cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl;
          cout << " (pos: " << linecount << " , " << nextfile << " , " << ninterpol << ")" << endl;
          if (snr1==snr2_now)
            {
            cout << " old2 = new1!" << endl;
            particle_data1=particle_data2;
            snr1_now = snr1;
            time1 = time2;
            }
          if (snr1_now!=snr1)
            {
            cout << " reading new1 " << snr1 << endl;
            gadget_reader(params,particle_data1,snr1,&time1);
            snr1_now = snr1;
            }
          if (snr2_now!=snr2)
            {
            cout << " reading new2 " << snr2 << endl;
            gadget_reader(params,particle_data2,snr2,&time2);
            snr2_now = snr2;
            }
#else
          gadget_reader(params,particle_data,0,&time);
#ifdef GEOMETRY_FILE
          p_orig = particle_data;
#endif
#endif
          break;
#if 0
        case 3:
          enzo_reader(params,particle_data);
          break;
#endif
        case 4:
          gadget_millenium_reader(params,particle_data,0,&time);
          break;
        case 5:
#if defined(USE_MPIIO)
          bin_reader_block_mpi(params,particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
#else
          planck_fail("mpi reader not available in non MPI compiled version !");
#endif
          break;
        case 6:
          mesh_reader(params,particle_data, maxr, minr);
          break;
#ifdef HDF5
        case 7:
          hdf5_reader(params,particle_data, maxr, minr);
          break;
#endif
        default:
          planck_fail("No valid file type given ...");
          break;
        }
#if defined(GEOMETRY_FILE) && !defined(INTERPOLATE)
      }
    else
      {
      particle_data = p_orig;
      }
#endif

#ifdef INTERPOLATE
    if (master)
      cout << "Interpolating between " << particle_data1.size() << " and " <<
        particle_data2.size() << " particles ..." << endl;
    particle_interpolate(params,particle_data,particle_data1,particle_data2,frac,time1,time2);
#endif

    long npart=particle_data.size();
    long npart_all=npart;
    mpiMgr.allreduce (npart_all,MPI_Manager::Sum);
    wallTimer.stop("read");

#ifndef CUDA
#ifndef NO_HOST_RANGING
// -----------------------------------
// ----------- Ranging ---------------
// -----------------------------------
    wallTimer.start("range");
    if (master)
      cout << endl << "ranging values (" << npart_all << ") ..." << endl;
    particle_normalize(params,particle_data,true); ///does log calculations and clamps data
    wallTimer.stop("range");
#endif

#ifndef NO_HOST_TRANSFORM
// -------------------------------------
// ----------- Transforming ------------
// -------------------------------------
    wallTimer.start("transform");
    if (master)
      cout << endl << "applying geometry (" << npart_all << ") ..." << endl;
    particle_project(params, particle_data, campos, lookat, sky);
    wallTimer.stop("transform");
#endif

// --------------------------------
// ----------- Sorting ------------
// --------------------------------
    wallTimer.start("sort");
    if (master)
      (mpiMgr.num_ranks()>1) ?
        cout << endl << "applying local sort ..." << endl :
        cout << endl << "applying sort (" << npart << ") ..." << endl;
    int sort_type = params.find<int>("sort_type",1);
    particle_sort(particle_data,sort_type,true);
    wallTimer.stop("sort");

#ifndef NO_HOST_COLORING
// ------------------------------------
// ----------- Coloring ---------------
// ------------------------------------
    wallTimer.start("coloring");
    if (master)
      cout << endl << "calculating colors (" << npart_all << ") ..." << endl;
    particle_colorize(params, particle_data, amap, emap);
    wallTimer.stop("coloring");
#endif

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------
    int res = params.find<int>("resolution",200);
    long nsplotch=particle_data.size();
    long nsplotch_all=nsplotch;
    mpiMgr.allreduce (nsplotch_all,MPI_Manager::Sum);
    if (master)
      cout << endl << "rendering (" << nsplotch_all << "/" << npart_all << ")..." << endl;
    arr2<COLOUR> pic(res,res);
    float64 grayabsorb = params.find<float>("gray_absorption",0.2);
    bool a_eq_e = params.find<bool>("a_eq_e",true);
#ifndef NO_HOST_RENDER
    wallTimer.start("render");
    render(particle_data,pic,a_eq_e,grayabsorb);
    wallTimer.stop("render");
#endif//NO_HOST_RENDER
#endif //if not def CUDA

#ifdef CUDA
    int res;
    arr2<COLOUR> pic;
    render_cuda(params,res,pic);
#endif

    wallTimer.start("write");

// ---------------------------------
// ----------- Colorbar ------------
// ---------------------------------
    if (master && params.find<bool>("colorbar",false))
      {
      cout << endl << "creating color bar ..." << endl;
      add_colorbar(params,pic,amap);
      }

// -------------------------------
// ----------- Saving ------------
// -------------------------------
    if (master)
      cout << endl << "saving file ..." << endl;

    int pictype = params.find<int>("pictype",0);
    string outfile = params.find<string>("outfile");
#ifdef GEOMETRY_FILE
    outfile = outfile + intToString(linecount,4) + ".tga";
    linecount++;
#ifdef INTERPOLATE
    if (linecount==nextfile)
      {
      nextfile=linecount+ninterpol;
      snr1=snr2;
      snr2++;
      }
#endif
#endif

    switch(pictype)
      {
      case 0:
        if (master) write_tga(params,pic,res,outfile);
        break;
      default:
        planck_fail("No valid image file type given ...");
        break;
      }

    wallTimer.stop("write");

#ifdef GEOMETRY_FILE
    for (int i=1; i<geometry_incr; i++)
      {
      getline(inp, line);
      linecount++;
#ifdef INTERPOLATE
      sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i",
             &campos.x,&campos.y,&campos.z,
             &lookat.x,&lookat.y,&lookat.z,
             &sky.x,&sky.y,&sky.z,&ninterpol);
      if (linecount==nextfile)
        {
        nextfile=linecount+ninterpol;
        snr1=snr2;
        snr2++;
        }
#endif
      }
    }
#endif



// -------------------------------
// ----------- Timings -----------
// -------------------------------
  wallTimer.stop("full");
  if (master)
    {
    cout << endl << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "Setup Data (secs)          : " << wallTimer.acc("setup") << endl;
    cout << "Read Data (secs)           : " << wallTimer.acc("read") << endl;
    cout << "Ranging Data (secs)        : " << wallTimer.acc("range") << endl;
    cout << "Transforming Data (secs)   : " << wallTimer.acc("transform") << endl;
    cout << "Sorting Data (secs)        : " << wallTimer.acc("sort") << endl;
    cout << "Coloring Sub-Data (secs)   : " << wallTimer.acc("coloring") << endl;
    cout << "Rendering Sub-Data (secs)  : " << wallTimer.acc("render") << endl;
    cout << "Write Data (secs)          : " << wallTimer.acc("write") << endl;
    cout << "Total (secs)               : " << wallTimer.acc("full") << endl;
    }

#ifdef VS
  //Just to hold the screen to read the messages when debugging
  cout << endl << "Press any key to end..." ;
  getchar();
#endif
  }
