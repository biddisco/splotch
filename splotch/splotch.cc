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
#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>

#ifdef USE_MPI
#include "mpi.h"
#else
#include <sys/time.h>
#endif

#include "mpi_support.h"
#include "arr.h"
#include "cxxutils.h"
#include "paramfile.h"
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"
#include "reader/bin_reader.h"

using namespace std;
using namespace RAYPP;

double times[100];

double myTime()
  {
#ifdef USE_MPI
  return MPI_Wtime();
#else
  using namespace std;
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1e-6*t.tv_usec;
//  return time(0);
#endif
  }


int main (int argc, char **argv)
{
  mpiMgr.startup (argc, argv);
  times[0] = myTime();
  bool master = mpiMgr.master();

  module_startup ("splotch",argc,2,"splotch <parameter file>",master);
  if (mpiMgr.num_ranks()>1)
    {
    cout << "Application was compiled with MPI support," << endl;
    cout << "running with " << mpiMgr.num_ranks() << " MPI tasks." << endl;
    }


// -----------------------------------
// ----------- Needed Data -----------
// -----------------------------------

  //  paramfile params (argv[1],master);
  paramfile params (argv[1],false);
  vector<particle_sim> particle_data;
  vector<particle_splotch> particle_col;
  VECTOR campos, lookat, sky;
  COLOURMAP amap,emap;


// ----------------------------------------------
// ----------- Loading Color Maps ---------------
// ----------------------------------------------

  if (master)
    cout << "building color maps ..." << endl;
  float rrr,ggg,bbb,rrr_old,ggg_old,bbb_old;
  ifstream infile (params.find<string>("palette0").c_str());
  string dummy;
  int nColours;
  infile >> dummy >> dummy >> nColours;
  cout << "Loading " << nColours << " entries of color table"  << endl;
  infile >> rrr_old >> ggg_old >> bbb_old;
  double step = 1./(nColours-1);
  for (int i=1; i<nColours; i++)
    {
    infile >> rrr >> ggg >> bbb;
    amap.Add_Entry(new LINEAR_CMAP_ENTRY((i-1)*step,i*step,
                                         COLOUR(rrr_old/255,ggg_old/255,bbb_old/255),
                                         COLOUR(rrr/255,ggg/255,bbb/255)));
    rrr_old=rrr; ggg_old=ggg; bbb_old=bbb;
    }
  emap=amap;


#ifdef NEVER_COMPILE
// ----------------------------------------------
// ------- How to build Parameter structre ------
// ------- and Color Maps without files ---------
// ----------------------------------------------

  map<string,string> par;
  par["infile"]="snap_92";
  par["simtype"]="1";
// and so on ...
  paramfile params (par);

  COLOUR c1,c2,c3
  c1=COLOUR(1,0,0);           // red
  c2=COLOUR(0.66,0.66,0.66);  // light gray
  c1=COLOUR(0,0,1);           // blue
  amap.Add_Entry(new LINEAR_CMAP_ENTRY( 0,.5,c1,c2));
  amap.Add_Entry(new LINEAR_CMAP_ENTRY(.5,1.,c2,c1));
  emap=amap;

#endif


  times[1] = myTime();

// -----------------------------------
// ----------- Reading ---------------
// -----------------------------------
  if (master)
    cout << "reading data ..." << endl;
  int simtype = params.find<int>("simtype");
  float maxr, minr;
  switch(simtype)
    {
    case 0:
      bin_reader_tab(particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
      break;
    case 1: 
      bin_reader_block(particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
      break;
    case 2: gadget_reader(params,particle_data);
      break;
    case 3: //enzo_reader(params,particle_data);
      break;
    default:
      planck_fail("No valid file type given ...");
      break;
    }
  long npart=particle_data.size();
  long npart_all=npart;
  mpiMgr.allreduce_sum (npart_all);
  times[2] = myTime();


// -----------------------------------
// ----------- Ranging ---------------
// -----------------------------------
  if (master)
    cout << "ranging values (" << npart_all << ") ..." << endl;
  particle_normalize(params,particle_data,true);
  times[3] = myTime();


#ifdef GEOMETRY_FILE

// -------------------------------------
// -- Looping over a flight path -------
// -------------------------------------
  vector<particle_sim> p_orig=particle_data;

  double cam_x,cam_y,cam_z,lat_x,lat_y,lat_z,sky_x,sky_y,sky_z;
  string line;
  ifstream inp(params.find<string>("geometry_file").c_str());
  int linecount=10000;
  int geometry_skip = params.find<int>("geometry_start",0);
  int geometry_incr = params.find<int>("geometry_incr",1);

  for(int i=0; i<geometry_skip; i++, linecount++)
    getline(inp, line);

  while (getline(inp, line))
    {
    particle_data = p_orig;
    sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &campos.x,&campos.y,&campos.z,
           &lookat.x,&lookat.y,&lookat.z,
           &sky.x,&sky.y,&sky.z);
      cout << "Camera: " << campos << endl;
      cout << "Lookat: " << lookat << endl;
      cout << "Sky:    " << sky << endl;
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


// -------------------------------------
// ----------- Transforming ------------
// -------------------------------------
      if (master)
	cout << "applying geometry (" << npart_all << ") ..." << endl;
      paticle_project(params, particle_data, campos, lookat, sky);
      times[4] = myTime();


// --------------------------------
// ----------- Sorting ------------
// --------------------------------
#ifdef USE_MPI
      if (master)
        cout << "applying local sort ..." << endl;
#else
      cout << "applying sort (" << npart << ") ..." << endl;
#endif
      int sort_type = params.find<int>("sort_type",1);
      particle_sort(particle_data,sort_type,true);
      times[5] = myTime();


// ------------------------------------
// ----------- Coloring ---------------
// ------------------------------------
      if (master)                        
        cout << "calculating colors (" << npart_all << ") ..." << endl;
      particle_colorize(params, particle_data, particle_col, amap, emap);
      times[6] = myTime();


// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------
      int res = params.find<int>("resolution",200);
      long nsplotch=particle_col.size();
      long nsplotch_all=nsplotch;
      mpiMgr.allreduce_sum (nsplotch_all);
      if (master)
        cout << "rendering (" << nsplotch_all << "/" << npart << ")..." << endl;
      arr2<COLOUR> pic(res,res);
      float64 grayabsorb = params.find<float>("gray_absorption",0.2);
      bool a_eq_e = params.find<bool>("a_eq_e",true);
      render(particle_col,pic,a_eq_e,grayabsorb);
      times[7] = myTime();


// ---------------------------------
// ----------- Colorbar ------------
// ---------------------------------
      if (master)
	{
	  bool colbar = params.find<bool>("colorbar",false);
	  if (colbar)
	    {
	      cout << "adding colour bar ..." << endl;
	      add_colorbar(pic,amap);
	    }
	}


// -------------------------------
// ----------- Saving ------------
// -------------------------------
      if (master)                             
        cout << "saving file ..." << endl;

      int pictype = params.find<int>("pictype",0);
      string outfile = params.find<string>("outfile");
#ifdef GEOMETRY_FILE
      outfile = outfile + "." + dataToString(linecount);
      linecount++;
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

#ifdef GEOMETRY_FILE
      for(int i=1; i<geometry_incr; i++, linecount++)
        getline(inp, line);
    }
#endif

// -------------------------------
// ----------- Timings -----------
// -------------------------------
  times[8] = myTime();
  if (master)
    {
    cout << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "Setup Data (secs)          : " << times[1]-times[0] << endl;
    cout << "Read Data (secs)           : " << times[2]-times[1] << endl;
    cout << "Ranging Data (secs)        : " << times[3]-times[2] << endl;
    cout << "Transforming Data (secs)   : " << times[4]-times[3] << endl;
    cout << "Sorting Data (secs)        : " << times[5]-times[4] << endl;
    cout << "Coloring Sub-Data (secs)   : " << times[6]-times[5] << endl;
    cout << "Rendering Sub-Data (secs)  : " << times[7]-times[6] << endl;
    cout << "Write Data (secs)          : " << times[8]-times[7] << endl;
    cout << "Total (secs)               : " << times[8]-times[0] << endl;
    }

  mpiMgr.shutdown();
}
