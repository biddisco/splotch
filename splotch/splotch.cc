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

#ifdef VS
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/arr.h"
#include "cxxsupport/cxxutils.h"
#include "cxxsupport/paramfile.h"
#else
#include "mpi_support.h"
#include "arr.h"
#include "cxxutils.h"
#include "paramfile.h"
#endif
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"
//#include "reader/bin_reader.h"

#ifdef VS	//Jin: for compiling under Windows/Visual Studio
#include "reader/gadget_reader.cc"
#include "writer/write_tga.cc"
#endif

#ifdef CUDA
#include "cuda/splotch_cuda.h"
void GoldComparePData
(vector<particle_sim> particle_data, cu_particle_sim* d_particle_data);
cu_color	k_get_color(int ptype, float val, cu_color_map_entry *map, //will move to kernel
			int	mapSize, int *ptype_points, int ptypes);
#endif

#ifdef USE_MPI
#include "mpi.h"
#else
	#ifdef VS
	#include "cuda/VTimer.h"
	#else
	#include <sys/time.h>
	#endif
#endif

using namespace std;
using namespace RAYPP;

double last_time,times[100];


double myTime()
  {
#ifdef USE_MPI
  return MPI_Wtime();
#else
#ifdef VS

  static	VTimer t;
  static	bool bFirstTimeRun =true;

  if (bFirstTimeRun)
  {
	  bFirstTimeRun =false;
	  t.start();
  }

  return t.getTime(); //in seconds

#else
  using namespace std;
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + 1e-6*t.tv_usec;
//  return time(0);
#endif
#endif
  }


int main (int argc, char **argv)
{
  mpiMgr.startup (argc, argv); ///does nothing actually
  last_time = myTime();
  times[0] = last_time;
  bool master = mpiMgr.master(); ///return true

  module_startup ("splotch",argc,2,"splotch <parameter file>",master);
  if (mpiMgr.num_ranks()>1 && master) ///as mpi manager always returns 1
    {
    cout << "Application was compiled with MPI support," << endl;
    cout << "running with " << mpiMgr.num_ranks() << " MPI tasks." << endl << endl;
    }

#ifdef INTERPOLATE
#ifndef GEOMETRY_FILE
  planck_fail("splotch: interpolation without geometry file makes no sense !");
#endif
#endif

// -----------------------------------
// ----------- Needed Data -----------
// -----------------------------------

  //  paramfile params (argv[1],master);
  paramfile params (argv[1],false);
  vector<particle_sim> particle_data; ///row data from file
  vector<particle_splotch> particle_col; ///used for screen coordinates, x, y, r, ro, a, e
  VECTOR campos, lookat, sky; ///A 3D vector class, designed for high efficiency.
  vector<COLOURMAP> amap,emap;
  int ptypes = params.find<int>("ptypes",1); ///each particle type has a color map
#ifdef INTERPOLATE
  vector<particle_sim> particle_data1,particle_data2;
  int snr_start = params.find<int>("snap_start",10);
  int ninterpol = params.find<int>("snap_interpol",8);  
  int snr1=0,snr2=0;
#endif

// ----------------------------------------------
// ----------- Loading Color Maps ---------------
// ----------------------------------------------

  amap.resize(ptypes);

  if (master)
    cout << "building color maps (" << ptypes << ")..." << endl;
  for(int itype=0;itype<ptypes;itype++)
    {
      if(params.find<bool>("color_is_vector"+dataToString(itype),false))
	{
	  if (master)
	    cout << " color of ptype " << itype << " is vector, so no colormap to load ..." << endl;
	}
      else
	{
	  float rrr,ggg,bbb,rrr_old,ggg_old,bbb_old;
	  ifstream infile (params.find<string>("palette"+dataToString(itype)).c_str());
          planck_assert (infile,"could not open palette file  <" + 
                                params.find<string>("palette"+dataToString(itype)) + ">");
	  string dummy;
	  int nColours;
	  infile >> dummy >> dummy >> nColours;
	  if(master) 
	    cout << " loading " << nColours << " entries of color table of ptype " << itype << endl;
	  infile >> rrr_old >> ggg_old >> bbb_old;
	  double step = 1./(nColours-1);
	  for (int i=1; i<nColours; i++)
	    {
	      infile >> rrr >> ggg >> bbb;
	      amap[itype].Add_Entry(new LINEAR_CMAP_ENTRY((i-1)*step,i*step,
							  COLOUR(rrr_old/255,ggg_old/255,bbb_old/255),
							  COLOUR(rrr/255,ggg/255,bbb/255)));
	      rrr_old=rrr; ggg_old=ggg; bbb_old=bbb;
	    }
	}
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

  times[1] += myTime() - last_time;
  last_time = myTime();


#ifdef GEOMETRY_FILE

// -------------------------------------
// -- Looping over a flight path -------
// -------------------------------------
  vector<particle_sim> p_orig;

  ifstream inp(params.find<string>("geometry_file").c_str());
  int linecount=0;
  int geometry_skip = params.find<int>("geometry_start",0);
  int geometry_incr = params.find<int>("geometry_incr",1);

  string line;
  for(int i=0; i<geometry_skip; i++, linecount++)
    getline(inp, line);

  while (getline(inp, line))
    {
      sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &campos.x,&campos.y,&campos.z,
	     &lookat.x,&lookat.y,&lookat.z,
	     &sky.x,&sky.y,&sky.z);
      if(master)
	{
	  cout << endl << "Next entry <" << linecount << "> in geometry file ..." << endl;
	  cout << " Camera: " << campos << endl;
	  cout << " Lookat: " << lookat << endl;
	  cout << " Sky:    " << sky << endl;
	}
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
      if (linecount == geometry_skip)
	{
#endif
	  if (master)
	    cout << endl << "reading data ..." << endl;
	  int simtype = params.find<int>("simtype"); ///2:Gadget2
	  float maxr, minr;
#ifdef INTERPOLATE
	  int ifrac = linecount % ninterpol;
	  int snr1_this = snr_start + (linecount / ninterpol);
	  int snr2_this = snr1_this + 1;
          double frac=(double)ifrac/(double)ninterpol;
#endif
	  switch(simtype)
	    {
	    case 0:
//	            bin_reader_tab(particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
	      break;
	    case 1: 
//	            bin_reader_block(particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
	      break;
	    case 2: 
#ifdef INTERPOLATE
	      cout << "File1: " << snr1_this << " , File2: " << snr2_this << " , interpol fac: " << frac << endl; 
	      cout << " (old files : " << snr1 << " , " << snr2 << ")" << endl; 
              if(snr2 == snr1_this)
		{
  	          cout << " old2 = new1 !" << endl; 
		  particle_data1=particle_data2;
		  snr1 = snr2;
		}
	      if(snr1_this != snr1)
		{
  	          cout << " reading new1 " << snr1_this << endl; 
		  gadget_reader(params,particle_data1,snr1_this);
		  snr1 = snr1_this;
		}
	      if(snr2_this != snr2)
		{
  	          cout << " reading new2 " << snr2_this << endl; 
		  gadget_reader(params,particle_data2,snr2_this);
		  snr2 = snr2_this;
		}
	      if (master)
		cout << "Interpolating between " << particle_data1.size() << " and " << 
		        particle_data2.size() << " particles ..." << endl; 
	      particle_interpolate(particle_data,particle_data1,particle_data2,frac);
#else
	      gadget_reader(params,particle_data,0); ///vector<particle_sim> particle_data;
#ifdef GEOMETRY_FILE
	      p_orig = particle_data;
#endif
#endif
	      break;
	    case 3: //enzo_reader(params,particle_data);
	      break;
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
      long npart=particle_data.size();
      long npart_all=npart;
      mpiMgr.allreduce_sum (npart_all); ///does nothing
      times[2] += myTime() - last_time;
      last_time = myTime();

// -----------------------------------
// ----------- Run Cuda -------------
// -----------------------------------
// After reading, ranging with device
#ifdef CUDA
	VTimer	timer;
	float	time;
	cu_particle_sim	*d_particle_data;
	d_particle_data =new cu_particle_sim[particle_data.size()];

timer.reset();
timer.start();
	cu_init();
timer.stop();
time =timer.getTime();
	cout << endl << "cu_init() cost time:" << time << "s" <<endl;

timer.reset();
timer.start();
	//copy data to local C-like array d_particle_data, in mid of developing only
	for (int i=0; i<particle_data.size(); i++)
		memcpy( &(d_particle_data[i]), &(particle_data[i]), sizeof(cu_particle_sim));
timer.stop();
time =timer.getTime();
	cout << endl << "Copying particles to device cost time:" << time << "s" <<endl;


timer.reset();
timer.start();
	//call cuda range
	cu_range(params, d_particle_data, particle_data.size());
timer.stop();
time =timer.getTime();
	cout << endl << "Ranging with device cost time:" << time << "s" <<endl;

	//then copy particles back, in mid of developing only
	for (int i=0; i<particle_data.size(); i++)
		memcpy( &(particle_data[i]),&(d_particle_data[i]), sizeof(cu_particle_sim));

	//then continue to do transformation...
timer.reset();
timer.start();
	double	c[3]={campos.x, campos.y, campos.z}, 
			l[3]={lookat.x, lookat.y, lookat.z}, 
			s[3]={sky.x, sky.y,	sky.z};
	cu_transform(params, particle_data.size(),c, l, s,d_particle_data);
timer.stop();
time =timer.getTime();

	cout << endl << "Transforming with device cost time:" << time<< "s" <<endl;

	//then copy particles back, in mid of developing only
	for (int i=0; i<particle_data.size(); i++)
		memcpy( &(particle_data[i]),&(d_particle_data[i]), sizeof(cu_particle_sim));
#endif

//#define NO_HOST_RANGING
#ifndef NO_HOST_RANGING
// -----------------------------------
// ----------- Ranging ---------------
// -----------------------------------
      if (master)
		cout << endl << "ranging values (" << npart_all << ") ..." << endl;
	  particle_normalize(params,particle_data,true); ///does log calculations and clamps data
      times[3] += myTime() - last_time;
      last_time = myTime();
#endif

//#define NO_HOST_TRANSFORM
#ifndef NO_HOST_TRANSFORM
// -------------------------------------
// ----------- Transforming ------------
// -------------------------------------
      if (master)
		cout << endl << "applying geometry (" << npart_all << ") ..." << endl;
      particle_project(params, particle_data, campos, lookat, sky);
      times[4] += myTime() - last_time;
      last_time = myTime();
#endif

// --------------------------------
// ----------- Sorting ------------
// --------------------------------
#ifdef USE_MPI
      if (master)
        cout << endl << "applying local sort ..." << endl;
#else
      cout << endl << "applying sort (" << npart << ") ..." << endl;
#endif
      int sort_type = params.find<int>("sort_type",1);
      particle_sort(particle_data,sort_type,true);
      times[5] += myTime() - last_time;
      last_time = myTime();

#ifdef CUDA //2009-10-22: sorting is ignored now to simplified things.
	  if ( sort_type !=0) //means needing sort
		  cout<< endl << "CAUSION: SORTING NEEDED!"<<endl;
#endif


//#define NO_HOST_COLORING
#ifndef NO_HOST_COLORING
// ------------------------------------
// ----------- Coloring ---------------
// ------------------------------------
      if (master)                        
        cout << endl << "calculating colors (" << npart_all << ") ..." << endl;
      particle_colorize(params, particle_data, particle_col, amap, emap);///things goes to array particle_col now
      times[6] += myTime() - last_time;
      last_time = myTime();
#endif

#ifdef CUDA//colorize with device
	//init C style colormap 
	cu_color_map_entry	*amapD;//amap for Device. emap is not used currently
	int	*amapDTypeStartPos; //begin indexes of ptypes in the linear amapD[]
	amapDTypeStartPos =new int[ptypes];
	int	curPtypeStartPos =0;
	int	size =0;
	//first we need to count all the entries to get colormap size
	for(int i=0; i<amap.size(); i++)
	{
	    vector<HANDLE_RAYPP<CMAP_ENTRY> > e;
		e = amap[i].Entry;
		for (int j=0; j<e.size(); j++)
			size++;
	}
	//then fill up the colormap amapD
	amapD =new cu_color_map_entry[size];
	int	index =0;
	for(int i=0; i<amap.size(); i++)
	{
	    vector<HANDLE_RAYPP<CMAP_ENTRY> > e;
		e = amap[i].Entry;
		int j;
		for (j=0; j<e.size(); j++)
		{
			HANDLE_RAYPP<CMAP_ENTRY> h;
			h =e[j];
			COLOUR	clr1, clr2;
			amapD[index].min =h->minval; 
			amapD[index].max =h->maxval;
			clr1 =h->Get_Colour(h->minval);
			clr2 =h->Get_Colour(h->maxval);
			amapD[index].color1.r =clr1.r;
			amapD[index].color1.g =clr1.g;
			amapD[index].color1.b =clr1.b;
			amapD[index].color2.r =clr2.r;
			amapD[index].color2.g =clr2.g;
			amapD[index].color2.b =clr2.b;
			index++;
		}
		amapDTypeStartPos[i] =curPtypeStartPos;
		curPtypeStartPos += j;
	}

	//debug code to test if C-style colormap works
	int	tmp=0;
	for (int i=0; i<size_inout_buffer; i++)
	{
		float	type=inout_buffer[i][0], val=inout_buffer[i][1], 
			r=inout_buffer[i][2], g=inout_buffer[i][3], 
			b=inout_buffer[i][4];
		
		cu_color	clr;
		clr =k_get_color(int(type), val, 
			amapD, size, amapDTypeStartPos, ptypes);
		if (i-tmp>1000)//for locate a bug
		{
			printf("\n%d\n", i);//set conditianal breakpoint here
			tmp =i;
			if ( r!=clr.r || g!=clr.g || b!=clr.b )//if not match, hold the screen
			{
				printf("\n%dth getcolor not match!\n", i);
				getchar();
			}
		}
	}

	//call device

	//delete C style colormap 
	delete	[]amapD;
	delete	[]amapDTypeStartPos;
#endif


// -----------------------------------
// ----------Gold Comparation --------
// -----------------------------------
#ifdef	CUDA 
	//compare to gold result
//	GoldComparePData( particle_data, d_particle_data);
	  //GoldCompareCData();//compare particle_col data
#endif

// -----------------------------------
// ----------- End Cuda --------------
// -----------------------------------
#ifdef CUDA
	cu_end();
	if (d_particle_data)
		delete []d_particle_data;
#endif


// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------
      int res = params.find<int>("resolution",200);
      long nsplotch=particle_col.size();
      long nsplotch_all=nsplotch;
      mpiMgr.allreduce_sum (nsplotch_all);
      if (master)
        cout << endl << "rendering (" << nsplotch_all << "/" << npart_all << ")..." << endl;
      arr2<COLOUR> pic(res,res);
      float64 grayabsorb = params.find<float>("gray_absorption",0.2);
      bool a_eq_e = params.find<bool>("a_eq_e",true);

#ifdef CUDA
//	  cu_render();
#endif

      render(particle_col,pic,a_eq_e,grayabsorb);
      times[7] += myTime() - last_time;
      last_time = myTime();


// ---------------------------------
// ----------- Colorbar ------------
// ---------------------------------
      if (master)
	{
	  bool colbar = params.find<bool>("colorbar",false);
	  if (colbar)
	    {
	      cout << endl << "creating color bar ..." << endl;
	      add_colorbar(params,pic,amap);
	    }
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

      times[8] += myTime() - last_time;
      last_time = myTime();

#ifdef GEOMETRY_FILE
      for(int i=1; i<geometry_incr; i++, linecount++)
        getline(inp, line);
    }
#endif

  

// -------------------------------
// ----------- Timings -----------
// -------------------------------
  times[9] = myTime();
  if (master)
    {
      cout << endl << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "Setup Data (secs)          : " << times[1] << endl;
    cout << "Read Data (secs)           : " << times[2] << endl;
    cout << "Ranging Data (secs)        : " << times[3] << endl;
    cout << "Transforming Data (secs)   : " << times[4] << endl;
    cout << "Sorting Data (secs)        : " << times[5] << endl;
    cout << "Coloring Sub-Data (secs)   : " << times[6] << endl;
    cout << "Rendering Sub-Data (secs)  : " << times[7] << endl;
    cout << "Write Data (secs)          : " << times[8] << endl;
    cout << "Total (secs)               : " << times[9]-times[0] << endl;
    }

  mpiMgr.shutdown();

#ifdef VS
  //Just to hold the screen to read the messages
  cout << endl << "Press any key to end..." ;
  getchar();
#endif
}


////////////////CUDA HELP FUNCTION//////////////////////////////
#ifdef CUDA
void GoldComparePData
(vector<particle_sim> particle_data, cu_particle_sim* d_particle_data)
{
	for (int i=0; i< particle_data.size(); i++)
	{
		if ( particle_data[i].x == d_particle_data[i].x &&
			particle_data[i].y == d_particle_data[i].y &&
			particle_data[i].z == d_particle_data[i].z &&
			particle_data[i].r == d_particle_data[i].r &&
			particle_data[i].ro == d_particle_data[i].ro &&
			particle_data[i].I == d_particle_data[i].I &&
			particle_data[i].C1 == d_particle_data[i].C1 &&
			particle_data[i].C2 == d_particle_data[i].C2 &&
			particle_data[i].C3 == d_particle_data[i].C3
			)
			continue;

		//else
		cout << endl << "Gold compare at particle " << i << endl;
		cout << "gold" << endl;
		cout << particle_data[i].x << ",";
		cout << particle_data[i].y << ",";
		cout << particle_data[i].z << ",";
		cout << particle_data[i].r << ",";
		cout << particle_data[i].ro << ",";
		cout << particle_data[i].I << ",";
		cout << particle_data[i].C1 << ",";
		cout << particle_data[i].C2 << ",";
		cout << particle_data[i].C3<< endl;

		cout << "cuda" << endl;
		cout << d_particle_data[i].x << ",";
		cout << d_particle_data[i].y << ",";
		cout << d_particle_data[i].z << ",";
		cout << d_particle_data[i].r << ",";
		cout << d_particle_data[i].ro << ",";
		cout << d_particle_data[i].I << ",";
		cout << d_particle_data[i].C1 << ",";
		cout << d_particle_data[i].C2 << ",";
		cout << d_particle_data[i].C3 << endl;

		//hold screen for reading
		cout << endl << "Press q to quit, other to continue..." <<endl;
		char c;
		//cin >> c;
		c =getchar();
		if (c == 'q')
			break;
	}
}

//get color from a kernel colormap. should move to kernel later
cu_color	k_get_color(int ptype, float val, cu_color_map_entry *map, //will move to kernel
			int	mapSize, int *ptype_points, int ptypes)
{
	cu_color	clr;
	clr.r =clr.g =clr.b =0.0;
	
	//first find the right entry for this ptype
	if (ptype>=ptypes)
		return clr; //invalid input
	int	start, end;
	start =ptype_points[ptype];
	if ( ptype == ptypes-1)//the last type
		end =mapSize-1;
	else 
		end =ptype_points[ptype+1]-1;
	
	//seach the section of this type to find the val
	for(int i=start; i<=end; i++)
	{
		if ( val>=map[i].min && val<=map[i].max)//if val falls into this entry, set clr
		{
			float	fract = (val-map[i].min)/(map[i].max-map[i].min);
			cu_color	clr1=map[i].color1, clr2=map[i].color2;
			clr.r =clr1.r + fract*(clr2.r-clr1.r);
			clr.g =clr1.g + fract*(clr2.g-clr1.g);
			clr.b =clr1.b + fract*(clr2.b-clr1.b);
			break;
		}
	}

	return clr;
}

#endif
