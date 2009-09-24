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
 * addings for test
 *
 */
#ifdef USEMPI
#include <mpi.h>
#endif
#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>

#include "arr.h"
#include "cxxutils.h"
#include "paramfile.h"
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"

using namespace std;
using namespace RAYPP;

#ifdef USEMPI
int ThisTask,NTask;
double times[100];
#endif

int main (int argc, char **argv)
{
#ifdef USEMPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  times[0] = MPI_Wtime();

  if(ThisTask == 0)
    {
#endif
      announce ("splotch");
#ifdef USEMPI
      cout << "Application was compiled with MPI support," << endl;
      cout << "runing with " << NTask << " MPI tasks." << endl << endl;
    }
#endif

  planck_assert (argc==2,"usage: splotch <parameter file>");
  paramfile params (argv[1],false);
  vector<particle_sim> p;
  long npart;

#ifdef USEMPI                      // ----------- Reading ---------------
  if(ThisTask==0)
#endif
    cout << "reading data ..." << endl;
  
  int simtype = params.find<int>("simtype");

  switch(simtype)
    {
    case 0: // npart=bin_reader(params,p);
      break;
    case 1: npart=gadget_reader(params,p);
      break;
    case 2: // npart=enzo_reader(params,p);
      break;
    default: 
#ifdef USEMPI
      if(ThisTask==0)
#endif
	cout << "No valid file type given ..." << endl;
      exit(1);
      break;
    }

#ifdef USEMPI                      // ----------- Ranging ---------------
  times[1] = MPI_Wtime();

  long npart_all;
  MPI_Allreduce(&npart, &npart_all, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask==0)
    cout << "ranging values (" << npart_all << ") ..." << endl;
#else
    cout << "ranging values (" << npart << ") ..." << endl;
#endif

  bool log_int = params.find<bool>("log_intensity0",true);
  bool log_col = params.find<bool>("log_color0",true);
  bool asinh_col = params.find<bool>("asinh_color0",false);
  bool col_vector = params.find<bool>("color_is_vector0",false);
  float32 mincol=1e30, maxcol=-1e30,minint=1e30, maxint=-1e30;

  for (int m=0; m<npart; ++m)
    {
      if(log_int) 
	p[m].I = log(p[m].I);
      minint=min(minint,p[m].I);
      maxint=max(maxint,p[m].I);
      if(log_col) 
	p[m].C1 = log(p[m].C1);
      if(asinh_col)
	p[m].C1 = log((double)p[m].C1+
		      sqrt((double)1+(double)p[m].C1*(double)p[m].C1)); 
      mincol=min(mincol,p[m].C1);
      maxcol=max(maxcol,p[m].C1);
      if(col_vector)
	{
	  if(log_col) 
	    {
	      p[m].C2 = log(p[m].C2);
	      p[m].C3 = log(p[m].C3);
	    }
	  if(asinh_col)
	    {
	      p[m].C2 = log((double)p[m].C2+
			    sqrt((double)1+(double)p[m].C2*(double)p[m].C2)); 
	      p[m].C3 = log((double)p[m].C3+
			    sqrt((double)1+(double)p[m].C3*(double)p[m].C3));
	    } 
	  mincol=min(mincol,p[m].C2);
	  mincol=min(mincol,p[m].C3);
	  maxcol=max(maxcol,p[m].C2);
	  maxcol=max(maxcol,p[m].C3);
	}

    }
#ifdef USEMPI
  float32 mincol_all=1e30, maxcol_all=-1e30,minint_all=1e30, maxint_all=-1e30;
  MPI_Allreduce(&mincol, &mincol_all, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&minint, &minint_all, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxcol, &maxcol_all, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxint, &maxint_all, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  mincol=mincol_all;
  minint=minint_all;
  maxcol=maxcol_all;
  maxint=maxint_all;
#endif

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
      if(minval_int != maxval_int) p[m].I = (p[m].I-minval_int)/(maxval_int-minval_int);
      if(minval_col != maxval_col) 
	{
	  p[m].C1 = (p[m].C1-minval_col)/(maxval_col-minval_col);
	  if(col_vector)
	    {
	      p[m].C2 = (p[m].C2-minval_col)/(maxval_col-minval_col);
	      p[m].C3 = (p[m].C3-minval_col)/(maxval_col-minval_col);
	    }
	}
    }


#ifdef USEMPI                      // ----------- Loading Color Maps ---------------
  times[2] = MPI_Wtime();

  if(ThisTask==0)
#endif
    cout << "building color maps ..." << endl;

  COLOURMAP amap,emap;

  FILE * pFile;
  string palname;
  int nColours;
  char foo[100];
  float dlut,startlut,rrr,ggg,bbb,rrr_old,ggg_old,bbb_old;

  palname = params.find<string>("palette0");
  pFile=fopen(palname.c_str(), "r");
  fscanf(pFile,"%s\n", foo ); // skip lut ide     == JASC-PAL
  fscanf(pFile,"%s\n", foo ); // skip lut version == 0100
  fscanf(pFile,"%d\n", &nColours ); 
  startlut = 0.0;
  dlut = 1.0/(float)(nColours-1);
  cout << "Loading " << nColours << " entries of color table"  << endl;             
  fscanf(pFile,"%f %f %f\n", &rrr_old, &ggg_old, &bbb_old );
  for (int i=1; i<nColours; i++)
    {
      fscanf(pFile,"%f %f %f\n", &rrr, &ggg, &bbb );
      amap.Add_Entry(new LINEAR_CMAP_ENTRY(startlut,startlut+dlut,
					   COLOUR(rrr_old/255,ggg_old/255,bbb_old/255),
					   COLOUR(rrr/255,ggg/255,bbb/255)));
      rrr_old=rrr;
      ggg_old=ggg;
      bbb_old=bbb;  
      startlut += dlut;
    }
  fclose(pFile);
  if(startlut < 1)
    amap.Add_Entry(new LINEAR_CMAP_ENTRY(startlut,1,
					  COLOUR(rrr_old/255,ggg_old/255,bbb_old/255),
					  COLOUR(rrr/255,ggg/255,bbb/255)));
  emap=amap;


#ifdef USEMPI                      // ----------- Transforming ------------
  times[3] = MPI_Wtime();

  if(ThisTask==0)
    cout << "applying geometry (" << npart_all << ") ..." << endl;
#else
    cout << "applying geometry (" << npart << ") ..." << endl;
#endif

  int res = params.find<int>("resolution",200);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);

#ifdef GEOMETRY_FILE
  vector<particle_sim> p_orig=p;

  double cam_x,cam_y,cam_z,lat_x,lat_y,lat_z,sky_x,sky_y,sky_z;
  string geometry_file = params.find<string>("geometry_file");
  string line;
  ifstream inp(geometry_file.c_str());
  int linecount=10000;
  int geometry_skip = params.find<int>("geometry_start",0);
  int geometry_incr = params.find<int>("geometry_incr",1);
    
  for(int i=0; i<geometry_skip; i++, linecount++)
    getline(inp, line);
  
  while (getline(inp, line))
    {
      for (long m=0; m<npart; ++m)
	{
	  p[m].x=p_orig[m].x;
	  p[m].y=p_orig[m].y;
	  p[m].z=p_orig[m].z;
	  p[m].r=p_orig[m].r;
	  p[m].ro=p_orig[m].ro;
	  p[m].I=p_orig[m].I;
	  p[m].C1=p_orig[m].C1;
	  p[m].C2=p_orig[m].C2;
	  p[m].C3=p_orig[m].C3;
	  p[m].type=p_orig[m].type;
	}
      sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &cam_x,&cam_y,&cam_z,
	     &lat_x,&lat_y,&lat_z,
	     &sky_x,&sky_y,&sky_z);
      VECTOR campos(cam_x,cam_y,cam_z);
      VECTOR lookat(lat_x,lat_y,lat_z);
      VECTOR sky(sky_x,sky_y,sky_z);
      cout << "Camera: (" << cam_x << "," << cam_y << "," << cam_z << ")" << endl;
      cout << "Lookat: (" << lat_x << "," << lat_y << "," << lat_z << ")" << endl;
      cout << "Sky:    (" << sky_x << "," << sky_y << "," << sky_z << ")" << endl;
      printf("%lf %lf %lf\n",cam_x,cam_y,cam_z);
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
      float64 dist=sqrt((campos.x-lookat.x) * (campos.x-lookat.x)
			+(campos.y-lookat.y) * (campos.y-lookat.y)
			+(campos.z-lookat.z) * (campos.z-lookat.z));
      float64 xfac=1./(fovfct*dist);
      cout << "Field of fiew: " << 1./xfac*2. << endl;
#else
      float64 xfac;
#endif
	
      for (long m=0; m<npart; ++m)
	{
#ifdef PROJECTION_OFF
	  p[m].x = res*.5 * (p[m].x+fovfct*dist)*xfac;
	  p[m].y = res*.5 * (p[m].y+fovfct*dist)*xfac;
#else
	  xfac=1./(fovfct*p[m].z);
	  p[m].x = res*.5 * (p[m].x+fovfct*p[m].z)*xfac;
	  p[m].y = res*.5 * (p[m].y+fovfct*p[m].z)*xfac;
#endif
	  p[m].ro = p[m].r;
	  p[m].r = p[m].r *res*.5*xfac;
#ifdef MINHSMLPIXEL
	  if ((p[m].r <= 0.5) && (p[m].r >= 0.0))
	    {
	      p[m].r = 0.5;
	      p[m].ro = p[m].r/res/0.5/xfac;
	    }
#endif
	}


#ifdef USEMPI                      // ----------- Sorting ------------
      times[4] = MPI_Wtime();

      if(ThisTask==0)
	cout << "applying local sort ..." << endl;
#else
	cout << "applying sort (" << npart << ") ..." << endl;
#endif

      int sort_type = params.find<int>("sort_type",1);

      switch(sort_type)
	{
	case 0: cout << "skiped sorting ..." << endl;
	  break;
	case 1: cout << "sorting by z ..." << endl;
	  sort(&p[0], &p[npart], zcmp());
	  break;
	case 2: cout << "sorting by value ..." << endl;
	  sort(&p[0], &p[npart], vcmp1());
	  break;
	case 3: cout << "reverse sorting by value ..." << endl;
	  sort(&p[0], &p[npart], vcmp2());
	  break;
	case 4: cout << "sorting by size ..." << endl;
	  sort(&p[0], &p[npart], hcmp());
	  break;
	default: cout << "unknown sorting choosen ..." << endl;
	  exit(1);
	  break;
	}


#ifdef USEMPI                      // ----------- Coloring ---------------
     times[5] = MPI_Wtime();

     if(ThisTask==0)
	cout << "calculating colors (" << npart_all << ") ..." << endl;
#else
	cout << "calculating colors (" << npart << ") ..." << endl;
#endif

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
	  col1=max(float64(0.0000001),min(float64(0.9999999),col1));
	  if(col_vector)
	    {
	      col2=max(float64(0.0000001),min(float64(0.9999999),col2));
	      col3=max(float64(0.0000001),min(float64(0.9999999),col3));
	    }
	  float64 intensity=p[m].I;
	  intensity=max(float64(0.0000001),min(float64(0.9999999),intensity));

	  COLOUR e;
	  if(col_vector)
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


#ifdef USEMPI                      // ----------- Rendering ------------
      times[6] = MPI_Wtime();

      long nsplotch_all;
      MPI_Allreduce(&nsplotch, &nsplotch_all, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      if(ThisTask==0)
	cout << "rendering (" << nsplotch_all << "/" << npart << ")..." << endl;
#else
	cout << "rendering (" << nsplotch << "/" << npart << ")..." << endl;
#endif

      bool a_eq_e = params.find<bool>("a_eq_e",true);

      arr2<COLOUR> pic(res,res);
      pic.fill(COLOUR(0,0,0));

      splotch_renderer renderer;
      renderer.render(p2,pic,a_eq_e,grayabsorb);

      bool colbar = params.find<bool>("colorbar",false);
#ifdef USEMPI                   
      times[7] = MPI_Wtime();

      if(ThisTask==0)
	{
#endif
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
#ifdef USEMPI                   
	}
                                // ----------- Saving ------------

      if(ThisTask==0)
#endif
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
#ifdef USEMPI                   
	  if(ThisTask==0)
#endif
	    write_tga(params,pic,res,outfile);
	  break;
	default: 
#ifdef USEMPI
	  if(ThisTask==0)
#endif
	    cout << "No valid image file type given ..." << endl;
	  exit(1);
	  break;
	}


#ifdef GEOMETRY_FILE
      for(int i=1; i<geometry_incr; i++, linecount++)
	getline(inp, line);
    }
#endif

#ifdef USEMPI
  times[8] = MPI_Wtime();

  if(ThisTask==0)
    {
      cout << "--------------------------------------------" << endl; 
      cout << "Summary of timings" << endl;
      cout << "Read Data (secs)           :" << times[1]-times[0] << endl;
      cout << "Ranging Data (secs)        :" << times[2]-times[1] << endl;
      cout << "Constructing Pallet (secs) :" << times[3]-times[2] << endl;
      cout << "Transforming Data (secs)   :" << times[4]-times[3] << endl;
      cout << "Sorting Data (secs)        :" << times[5]-times[4] << endl;
      cout << "Coloring Sub-Data (secs)   :" << times[6]-times[5] << endl;
      cout << "Rendering Sub-Data (secs)  :" << times[7]-times[6] << endl;
      cout << "Write Data (secs)          :" << times[8]-times[7] << endl;
      cout << "Total (secs)               :" << times[8]-times[0] << endl;
   }

  MPI_Finalize();
#endif

  exit(0);
}


