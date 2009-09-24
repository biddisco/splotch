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

template<typename T> void get_minmax (T &minv, T &maxv, T val)
  {
  minv=min(minv,val);
  maxv=max(maxv,val);
  }

double my_asinh (double val)
  {
  return log(val+sqrt(1.+val*val));
  }

template<typename T> void my_normalize (T minv, T maxv, T &val)
  {
  if (minv!=maxv) val =  (val-minv)/(maxv-minv);
  }

template<typename T> void clamp (T minv, T maxv, T &val)
  {
  val = min(maxv, max(minv, val));
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
  paramfile params (argv[1],master);

  if (master)                             // ----------- Reading ---------------
    cout << "reading data ..." << endl;

  int simtype = params.find<int>("simtype");

  float maxr, minr;
  vector<particle_sim> p;
  switch(simtype)
    {
    case 0: 
      bin_reader_tab(p, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
      break;
    case 1: 
      bin_reader_block(p, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
      break;
    case 2: 
      gadget_reader(params,p);
      break;
    case 3: //enzo_reader(params,p);
      break;
    default:
      planck_fail("No valid file type given ...");
      break;
    }

  long npart=p.size();
  long npart_all=npart;
  mpiMgr.allreduce_sum (npart_all);

  if (master)                             // ----------- Ranging ---------------
    cout << "ranging values (" << npart_all << ") ..." << endl;
  times[1] = myTime();

  bool log_int = params.find<bool>("log_intensity0",true);
  bool log_col = params.find<bool>("log_color0",true);
  bool asinh_col = params.find<bool>("asinh_color0",false);
  bool col_vector = params.find<bool>("color_is_vector0",false);
  float32 mincol=1e30, maxcol=-1e30,minint=1e30, maxint=-1e30;

  for (int m=0; m<npart; ++m)
    {
    if (log_int)
      p[m].I = log(p[m].I);
    get_minmax(minint, maxint, p[m].I);
    if (log_col)
      p[m].C1 = log(p[m].C1);
    if(asinh_col)
      p[m].C1 = my_asinh(p[m].C1);
    get_minmax(mincol, maxcol, p[m].C1);
    if (col_vector)
      {
      if (log_col)
        {
        p[m].C2 = log(p[m].C2);
        p[m].C3 = log(p[m].C3);
        }
      if (asinh_col)
        {
        p[m].C2 = my_asinh(p[m].C2);
        p[m].C3 = my_asinh(p[m].C3);
        }
      get_minmax(mincol, maxcol, p[m].C2);
      get_minmax(mincol, maxcol, p[m].C3);
      }
    }

  mpiMgr.allreduce_min(minint);
  mpiMgr.allreduce_min(mincol);
  mpiMgr.allreduce_max(maxint);
  mpiMgr.allreduce_max(maxcol);

  float minval_int = params.find<float>("min_int0",minint);
  float maxval_int = params.find<float>("max_int0",maxint);
  float minval_col = params.find<float>("min_col0",mincol);
  float maxval_col = params.find<float>("max_col0",maxcol);

  cout << "From data: " << endl;
  cout << "Color Range:     " << mincol << " (min) , " <<
                                 maxcol << " (max) " << endl;
  cout << "Intensity Range: " << minint << " (min) , " <<
                                 maxint << " (max) " << endl;
  cout << "Restricted to: " << endl;
  cout << "Color Range:     " << minval_col << " (min) , " <<
                                 maxval_col << " (max) " << endl;
  cout << "Intensity Range: " << minval_int << " (min) , " <<
                                 maxval_int << " (max) " << endl;

  for (int m=0; m<npart; ++m)
    {
    my_normalize(minval_int,maxval_int,p[m].I);
    my_normalize(minval_col,maxval_col,p[m].C1);
    if (col_vector)
      {
      my_normalize(minval_col,maxval_col,p[m].C2);
      my_normalize(minval_col,maxval_col,p[m].C3);
      }
    }


  if (master)                  // ----------- Loading Color Maps ---------------
    cout << "building color maps ..." << endl;

  times[2] = myTime();

  COLOURMAP amap;

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

  COLOURMAP emap=amap;

  if (master)                           // ----------- Transforming ------------
    cout << "applying geometry (" << npart_all << ") ..." << endl;
  times[3] = myTime();

  int res = params.find<int>("resolution",200);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);

#ifdef GEOMETRY_FILE
  vector<particle_sim> p_orig=p;

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
    p = p_orig;

    VECTOR campos, lookat, sky;
    sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
           &campos.x,&campos.y,&campos.z,
           &lookat.x,&lookat.y,&lookat.z,
           &sky.x,&sky.y,&sky.z);
      cout << "Camera: " << campos << endl;
      cout << "Lookat: " << lookat << endl;
      cout << "Sky:    " << sky << endl;
#else
      VECTOR campos(params.find<double>("camera_x"),
                    params.find<double>("camera_y"),
                    params.find<double>("camera_z"));
      VECTOR lookat(params.find<double>("lookat_x"),
                    params.find<double>("lookat_y"),
                    params.find<double>("lookat_z"));
      VECTOR sky(params.find<double>("sky_x",0),
                 params.find<double>("sky_y",1),
                 params.find<double>("sky_z",0));
#endif
      sky.Normalize();
      VECTOR zaxis = (lookat-campos).Norm();
      VECTOR xaxis = Cross (sky,zaxis).Norm();
      VECTOR yaxis = Cross (zaxis,xaxis);
      TRANSFORM trans;
      trans.Make_General_Transform
        (TRANSMAT(xaxis.x,xaxis.y,xaxis.z,
                  yaxis.x,yaxis.y,yaxis.z,
                  zaxis.x,zaxis.y,zaxis.z,
                  0,0,0));
      trans.Invert();
      TRANSFORM trans2;
      trans2.Make_Translation_Transform(-campos);
      trans2.Add_Transform(trans);
      trans=trans2;
      //  trans.Add_Transform(trans2);

      for (long m=0; m<npart; ++m)
        {
        VECTOR v(p[m].x,p[m].y,p[m].z);
        v=trans.TransPoint(v);
        p[m].x=v.x; p[m].y=v.y; p[m].z=v.z;
        }
#ifdef PROJECTION_OFF
      float64 dist= (campos-lookat).Length();
      float64 xfac=1./(fovfct*dist);
      cout << "Field of fiew: " << 1./xfac*2. << endl;
#endif

      for (long m=0; m<npart; ++m)
        {
#ifdef PROJECTION_OFF
        p[m].x = res*.5 * (p[m].x+fovfct*dist)*xfac;
        p[m].y = res*.5 * (p[m].y+fovfct*dist)*xfac;
#else
        float64 xfac=1./(fovfct*p[m].z);
        p[m].x = res*.5 * (p[m].x+fovfct*p[m].z)*xfac;
        p[m].y = res*.5 * (p[m].y+fovfct*p[m].z)*xfac;
#endif
        p[m].ro = p[m].r;
        p[m].r = p[m].r *res*.5*xfac;
#ifdef MINHSMLPIXEL
        if ((p[m].r <= 0.5) && (p[m].r >= 0.0))
          {
          p[m].r = 0.5;
          p[m].ro = p[m].r/(res*.5*xfac);
          }
#endif
        }


#ifdef USE_MPI
      if (master)                            // ----------- Sorting ------------
        cout << "applying local sort ..." << endl;
#else
      cout << "applying sort (" << npart << ") ..." << endl;
#endif
      times[4] = myTime();

      int sort_type = params.find<int>("sort_type",1);

      switch(sort_type)
        {
        case 0: cout << "skipped sorting ..." << endl;
          break;
        case 1: cout << "sorting by z ..." << endl;
          sort(p.begin(), p.end(), zcmp());
          break;
        case 2: cout << "sorting by value ..." << endl;
          sort(p.begin(), p.end(), vcmp1());
          break;
        case 3: cout << "reverse sorting by value ..." << endl;
          sort(p.begin(), p.end(), vcmp2());
          break;
        case 4: cout << "sorting by size ..." << endl;
          sort(p.begin(), p.end(), hcmp());
          break;
        default:
          planck_fail("unknown sorting choosen ...");
          break;
        }

      if (master)                        // ----------- Coloring ---------------
        cout << "calculating colors (" << npart_all << ") ..." << endl;

     times[5] = myTime();

      int ycut0 = params.find<int>("ycut0",0);
      int ycut1 = params.find<int>("ycut1",res);
      float zmaxval = params.find<float>("zmax",1.e23);
      float zminval = params.find<float>("zmin",0.0);
      float64 brightness = params.find<double>("brightness",1.);
      float64 grayabsorb = params.find<float>("gray_absorption",0.2);

      float64 rfac=1.5;
      vector<particle_splotch> p2;
      for (int m=0; m<npart; ++m)
        {
        if (p[m].z<=0) continue;
        if (p[m].z<=zminval) continue;
        if (p[m].z>=zmaxval) continue;

        float64 r=p[m].r;
        float64 posx=p[m].x, posy=p[m].y;

        float64 rfacr=rfac*r;

        int minx=int(posx-rfacr+1);
        if (minx>=res) continue;
        minx=max(minx,0);
        int maxx=int(posx+rfacr+1);
        if (maxx<=0) continue;
        maxx=min(maxx,res);
        if (minx>=maxx) continue;
        int miny=int(posy-rfacr+1);
        if (miny>=ycut1) continue;
        miny=max(miny,ycut0);
        int maxy=int(posy+rfacr+1);
        if (maxy<=ycut0) continue;
        maxy=min(maxy,ycut1);
        if (miny>=maxy) continue;

        float64 col1=p[m].C1,col2=p[m].C2,col3=p[m].C3;
        clamp (0.0000001,0.9999999,col1);
        if (col_vector)
          {
          clamp (0.0000001,0.9999999,col2);
          clamp (0.0000001,0.9999999,col3);
          }
        float64 intensity=p[m].I;
        clamp (0.0000001,0.9999999,intensity);

        COLOUR e;
        if (col_vector)
          {
          e.r=col1*intensity*brightness;
          e.g=col2*intensity*brightness;
          e.b=col3*intensity*brightness;
          }
        else
          e=amap.Get_Colour(col1)*intensity*brightness;

        COLOUR a=e;

        p2.push_back(particle_splotch(p[m].x, p[m].y,p[m].r,p[m].ro,a,e));
        }
      long nsplotch=p2.size();
      long nsplotch_all=nsplotch;
      mpiMgr.allreduce_sum (nsplotch_all);

      if (master)                          // ----------- Rendering ------------
        cout << "rendering (" << nsplotch_all << "/" << npart << ")..." << endl;

      times[6] = myTime();

      bool a_eq_e = params.find<bool>("a_eq_e",true);

      arr2<COLOUR> pic(res,res);
      pic.fill(COLOUR(0,0,0));

      splotch_renderer renderer;
      renderer.render(p2,pic,a_eq_e,grayabsorb);

      bool colbar = params.find<bool>("colorbar",false);

      times[7] = myTime();

      if (master)
        {
        if (colbar)
          {
          cout << "adding colour bar ..." << endl;
          for (int x=0; x<res; x++)
            {
            float64 temp=x/float64(res);
            COLOUR e=amap.Get_Colour(temp);
            for (int y=0; y<10; y++)
              {
              pic[x][res-1-y].r = e.r;
              pic[x][res-1-y].g = e.g;
              pic[x][res-1-y].b = e.b;
              }
            }
          }
        }

      if (master)                             // ----------- Saving ------------
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

  times[8] = myTime();

  if (master)
    {
    cout << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "Read Data (secs)           : " << times[1]-times[0] << endl;
    cout << "Ranging Data (secs)        : " << times[2]-times[1] << endl;
    cout << "Constructing Palette (secs): " << times[3]-times[2] << endl;
    cout << "Transforming Data (secs)   : " << times[4]-times[3] << endl;
    cout << "Sorting Data (secs)        : " << times[5]-times[4] << endl;
    cout << "Coloring Sub-Data (secs)   : " << times[6]-times[5] << endl;
    cout << "Rendering Sub-Data (secs)  : " << times[7]-times[6] << endl;
    cout << "Write Data (secs)          : " << times[8]-times[7] << endl;
    cout << "Total (secs)               : " << times[8]-times[0] << endl;
    }

  mpiMgr.shutdown();
  }
