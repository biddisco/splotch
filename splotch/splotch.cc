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
//switch off host operations. may turn them on sometime for testing.
#define NO_HOST_RANGING
#define NO_HOST_TRANSFORM
#define NO_HOST_COLORING
#define NO_HOST_RENDER
//#define	CUDA_DEVICE_COMBINE
//include head files
#include "cuda/splotch_cuda.h"
//function definitions for cuda/testing use
void GoldComparePData
(vector<particle_sim> particle_data, cu_particle_sim* d_particle_data);
void GoldCompareSData
(vector<particle_splotch> host_data, cu_particle_splotch* device_data);
cu_color	C_get_color(int ptype, float val, cu_color_map_entry *map, //will move to kernel
			int	mapSize, int *ptype_points, int ptypes);
void GoldCompareFBuf(cu_fragment_AeqE *goldBuf, cu_fragment_AeqE *buf, int n);
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
#ifdef INTERPOLATE          // Here only the tow datasets are prepared, interpolation will be done later
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

#ifdef INTERPOLATE
      if (master)
	cout << "Interpolating between " << particle_data1.size() << " and " << 
	  particle_data2.size() << " particles ..." << endl; 
      particle_interpolate(particle_data,particle_data1,particle_data2,frac);
#endif

      long npart=particle_data.size();
      long npart_all=npart;
      mpiMgr.allreduce_sum (npart_all); ///does nothing
      times[2] += myTime() - last_time;
      last_time = myTime();

/////////////Here comes CUDA code//////////////////////////
#ifndef CUDA
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

//2009-10-22: sorting is ignored now to simplified things.
//if sorted, it's needed to copy particle sim back to device
#ifdef CUDA 
	  if ( sort_type !=0) //means needing sort
		  cout<< endl << "CAUSION: SORTING NEEDED!"<<endl;
#endif

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
#ifndef NO_HOST_RENDER
      render(particle_col,pic,a_eq_e,grayabsorb);
      times[7] += myTime() - last_time;
      last_time = myTime();
#endif//NO_HOST_RENDER

#endif //ifndef CUDA


// -----------------------------------
// ----------- Run Cuda -------------
// -----------------------------------
// After reading, ranging with device
#ifdef CUDA
	//CUDA test. for developing only. cut short particle_data
	vector<particle_sim>::iterator it;
	int	testPCount =100000/4, n=0;
	it =particle_data.begin();
	particle_data.erase(it+testPCount, particle_data.end());

	VTimer	timer;
	float	time;
	cu_particle_sim	*d_particle_data;
	d_particle_data =new cu_particle_sim[particle_data.size()];

	//CUDA Init 
timer.reset();
timer.start();
	cu_init();
timer.stop();
time =timer.getTime();
	cout << endl << "cu_init() cost time:" << time << "s" <<endl;

	//Copy particle sim into C-style object d_particle_data
timer.reset();
timer.start();
	//copy data to local C-like array d_particle_data, in mid of developing only
	for (int i=0; i<particle_data.size(); i++)
		memcpy( &(d_particle_data[i]), &(particle_data[i]), sizeof(cu_particle_sim));
timer.stop();
time =timer.getTime();
	cout << endl << "Copying particles to device cost time:" << time << "s" <<endl;

	//CUDA Ranging
timer.reset();
timer.start();
	//call cuda range
	cu_range(params, d_particle_data, particle_data.size());
timer.stop();
time =timer.getTime();
	cout << endl << "Ranging with device cost time:" << time << "s" <<endl;

/*	//then copy particles back, in mid of developing only
	for (int i=0; i<particle_data.size(); i++)
		memcpy( &(particle_data[i]),&(d_particle_data[i]), sizeof(cu_particle_sim));
*/	
	//CUDA RANGING DONE!

	//CUDA Transformation
timer.reset();
timer.start();
	double	c[3]={campos.x, campos.y, campos.z}, 
			l[3]={lookat.x, lookat.y, lookat.z}, 
			s[3]={sky.x, sky.y,	sky.z};
	cu_transform(params, particle_data.size(),c, l, s,d_particle_data);
timer.stop();
time =timer.getTime();
	cout << endl << "Transforming with device cost time:" << time<< "s" <<endl;

	//then copy particles back to host for sorting
	for (int i=0; i<particle_data.size(); i++)
		memcpy( &(particle_data[i]),&(d_particle_data[i]), sizeof(cu_particle_sim));

// --------------------------------
// ----------- Sorting ------------
// --------------------------------
      cout << endl << "applying sort (" << npart << ") ..." << endl;

      int sort_type = params.find<int>("sort_type",1);
      particle_sort(particle_data,sort_type,true);
      times[5] += myTime() - last_time;
      last_time = myTime();

//2009-10-22: sorting is ignored now to simplified things.
//if sorted, it's needed to copy particle sim back to device
//	  if ( sort_type !=0) //means needing sort
//		  cout<< endl << "CAUSION: SORTING NEEDED!"<<endl;
	  //copy sorted data back to device
	  //first to C-style object
	  for(int i=0; i<particle_data.size(); i++)
		memcpy( &(d_particle_data[i]), &(particle_data[i]), sizeof(cu_particle_sim));
	  cu_copy_particle_sim_to_device(d_particle_data, particle_data.size());

	  //CUDA Coloring
timer.reset();
timer.start();
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
	//now let cuda init colormap on device
	cu_colormap_info info;
	info.map =amapD;
	info.mapSize =size;
	info.ptype_points =amapDTypeStartPos;
	info.ptypes =ptypes;
	cu_init_colormap(info);

#ifdef CUDA_TEST_COLORMAP
	//debug code to test if C-style colormap works
	int	tmp=0;
	for (int i=0; i<size_inout_buffer; i++)
	{
		float	type=inout_buffer[i][0], val=inout_buffer[i][1], 
			r=inout_buffer[i][2], g=inout_buffer[i][3], 
			b=inout_buffer[i][4];
		
		if (type==0) continue;

		cu_color	clr;
//		clr =C_get_color(int(type), val, 
//			amapD, size, amapDTypeStartPos, ptypes);
		clr =cu_get_color((int)type, val);
		if ( r!=clr.r || g!=clr.g || b!=clr.b )//if not match, hold the screen
		{
			printf("\n%dth getcolor not match!\n", i);
			printf("\nHOST:   (%d,%f) (%f, %f, %f)\n",(int)type, val, r, g, b);
			printf("\nDEVICE: (%d,%f) (%f, %f, %f)\n",(int)type, val, clr.r, clr.g, clr.b);
			getchar();
		}
		if (i-tmp>1000)//output current position
		{
			printf("\n%d\n", i);//set conditianal breakpoint here
			tmp =i;
		}
	}
#endif //ifdef CUDA_TEST_COLORMAP
	
	//init particle_splotch array memeory
	cu_particle_splotch	*cu_ps;
	size =particle_data.size();
	cu_ps =new cu_particle_splotch[size];
	memset(cu_ps, 0, size);

	//Colorize with device
	cu_colorize(params, cu_ps, size);

	//filter particle_splotch array to a cu_ps_filtered
	cu_particle_splotch	*cu_ps_filtered;
	size =particle_data.size();
	cu_ps_filtered =new cu_particle_splotch[size];
	int	pFiltered=0;
	//do filtering
	unsigned long	posInFragBuf =0;//, countFragments;
	int		minx=1e6,miny=1e6, maxx=-1,maxy=-1;
	for (int i=0; i<size; i++)
	{
		if ( cu_ps[i].isValid )
		{
			//IMPORTANT: set the start position of this particle
			//in fragment buffer
			cu_ps[i].posInFragBuf =posInFragBuf;
			posInFragBuf += (cu_ps[i].maxx -cu_ps[i].minx)*
				(cu_ps[i].maxy -cu_ps[i].miny);
			memcpy(&cu_ps_filtered[pFiltered], 
				&cu_ps[i], sizeof(cu_particle_splotch));
			pFiltered++;

			minx=min(minx,cu_ps[i].minx);
			miny=min(miny,cu_ps[i].miny);
			maxx=max(maxx,cu_ps[i].maxx);
			maxy=max(maxy,cu_ps[i].maxy);
		}
	}

/*	//do gold comparation
//	GoldCompareSData(particle_col, cu_ps_filtered);

	//debug only: then copy data back to host p2 and continue
	//working only with NO_HOST_COLORING
	for (int i=0; i<pFiltered; i++)
	{
		cu_color	a=cu_ps_filtered[i].a, e =cu_ps_filtered[i].e;
		COLOUR	A(a.r,a.g,a.b), E(e.r,e.g,e.b);
		particle_col.push_back( particle_splotch
			(cu_ps_filtered[i].x, cu_ps_filtered[i].y,
			 cu_ps_filtered[i].r,cu_ps_filtered[i].ro,
			 A,E));
	}
*/
timer.stop();
time =timer.getTime();
	cout << endl << "Coloring with device cost time:" << time << "s" <<endl;

	//here is the point that particle sim array on deivce is useless

	//just for rendering time coming next
    last_time = myTime();

/*	//compare to gold result, test only
//	GoldComparePData( particle_data, d_particle_data);
*/	  //GoldCompareCData();//compare particle_col data

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------

	//get parameters for rendering
	int res = params.find<int>("resolution",200);
	long nsplotch=pFiltered;
	long nsplotch_all=nsplotch;
	mpiMgr.allreduce_sum (nsplotch_all);
	if (master)
		cout << endl << "rendering (" << nsplotch_all << "/" << npart_all << ")..." << endl;
	arr2<COLOUR> pic(res,res);
	float64 grayabsorb = params.find<float>("gray_absorption",0.2);
	bool a_eq_e = params.find<bool>("a_eq_e",true);


/*code for test fragment buffer only
timer.reset();
timer.start();
	  //render with render_cu_test()
	  //it now writes to fragBufWrittenByHost!!!
//	  render_cu_test(cu_ps_filtered,pFiltered, pic,a_eq_e,grayabsorb);
timer.stop();
time =timer.getTime();
	cout << endl << "Render to fragment buffer cost time:" << time << "s" <<endl;
*/

	//CUDA Rendering with device
	//must initialize pic first!
    pic.fill(COLOUR(0,0,0));

//here's the point of multi-go loop starts
#ifndef CUDA_DEVICE_COMBINE //combined by host
	//prepare fragment buffer memory space first
	cu_fragment_AeqE	*fragBuf;
	int	nFBufInByte =cu_get_fbuf_size();
	int	nFBufInCell =nFBufInByte/sizeof(cu_fragment_AeqE);
	fragBuf =new cu_fragment_AeqE[nFBufInCell];
	

timer.reset();
timer.start();

	cu_prepare_render(cu_ps_filtered,pFiltered);

	bool bFinished=false;
	int	 startP, endP;
	startP =endP =0;
	int	 nFragments2Render=0;
	//cu_ps_filtered:	the particle array
	//pFiltered:		length of particle array
	//pPosInFragBuf:	length of the whole fragment buffer counted after filterign
	while (!bFinished)
	{
		//find a chunk
		nFragments2Render =0;
		for (startP =endP; endP<pFiltered; endP++)
		{
			int	inc;//increasment
			inc =(cu_ps_filtered[endP].maxx-cu_ps_filtered[endP].minx) *
				(cu_ps_filtered[endP].maxy -cu_ps_filtered[endP].miny);
			if ( nFragments2Render +inc > nFBufInCell)
				break;
			else
				nFragments2Render +=inc;
		}
		if (endP == pFiltered)
			bFinished =true;
		
		//render it
		cu_render1(startP, endP, a_eq_e, grayabsorb);

		//collect result
		cu_get_fbuf(fragBuf, nFragments2Render);

		//combine it
		cu_fragment_AeqE	*bufWrite;//just for test
		bufWrite=fragBuf;
		//bufWrite =fragBufWrittenByHost;
		for (int pPos=startP, fPos=0; pPos<endP && fPos<nFragments2Render; pPos++)
		{
		  for (int x =cu_ps_filtered[pPos].minx; 
			   x <cu_ps_filtered[pPos].maxx; x++)
		  {
			 for( int y =cu_ps_filtered[pPos].miny;
				  y <cu_ps_filtered[pPos].maxy; y++)
			 {
				 pic[x][y].r += bufWrite[fPos].deltaR;
				 pic[x][y].g += bufWrite[fPos].deltaG;
				 pic[x][y].b += bufWrite[fPos].deltaB;
				 fPos ++;
			 }				  
		  }
		}
	}
timer.stop();
time =timer.getTime();
		cout << endl << "Render to fragment buffer cost time:" << time << "s" <<endl;
/////////////////////////////////////////////////////////////////////

	//post-process
timer.reset();
timer.start();
     exptable	xexp(MAX_EXP);
     for(int ix=0; ix<pic.size1(); ix++)
        for(int iy=0; iy<pic.size2(); iy++)
          {
			  pic[ix][iy].r=1-xexp(pic[ix][iy].r);
			  pic[ix][iy].g=1-xexp(pic[ix][iy].g);
			  pic[ix][iy].b=1-xexp(pic[ix][iy].b);
          }
timer.stop();
time =timer.getTime();
	cout << endl << "Post-process pic[] cost time:" << time << "s" <<endl;
#endif //CUDA_DEVICE_COMBINE
	
#ifdef	CUDA_DEVICE_COMBINE
timer.reset();
timer.start();
	cu_prepare_render();
	cu_init_pic(pic.size1(), pic.size2());
	cu_render(cu_ps_filtered, pFiltered,/* pic.size1(),
		pic.size2(),*/a_eq_e, grayabsorb);

	cu_param_combine	cInfo;
	cInfo.xres =pic.size1();
	cInfo.yres =pic.size2();
	cInfo.minx =minx;	cInfo.miny =miny;
	cInfo.maxx =maxx;	cInfo.maxy =maxy;
	cInfo.lenFBuf =posInFragBuf;
	cInfo.pStart =0;
	cInfo.pEnd =pFiltered-1;//as index
	
	cu_combine(cInfo);

	cu_post_process(pic.size1(),pic.size2());
	
	cu_color	*cu_pic=0;
	cu_pic =new cu_color[pic.size1()*pic.size2()];
	memset(cu_pic, 0, sizeof(cu_color)*pic.size1()*pic.size2());
	cu_get_pic(cu_pic, pic.size1(), pic.size2());
	//then dump 1-d cu_pic to 2-d pic(host)
	for (int x=0; x<pic.size1(); x++)
		for (int y=0; y<pic.size2(); y++)
		{
			int idx =x*pic.size2()+y;
			pic[x][y].r =cu_pic[idx].r;
			pic[x][y].g =cu_pic[idx].g;
			pic[x][y].b =cu_pic[idx].b;
		}

/*	//post-process, just for test
     exptable	xexp(MAX_EXP);
     for(int ix=0; ix<pic.size1(); ix++)
        for(int iy=0; iy<pic.size2(); iy++)
          {
			  pic[ix][iy].r=1-xexp(pic[ix][iy].r);
			  pic[ix][iy].g=1-xexp(pic[ix][iy].g);
			  pic[ix][iy].b=1-xexp(pic[ix][iy].b);
          }*/
timer.stop();
time =timer.getTime();
		cout << endl << "Render on device cost time:" << time << "s" <<endl;
#endif	CUDA_DEVICE_COMBINE

//#define CUDA_TEST_EXP
#ifdef	CUDA_TEST_EXP //for testing exponent table
	//init exp table
	cu_init_exptab(MAX_EXP);
	//debug code to test if device exptable works
	printf("\nCheck exp table...\n");
	for (int i=0; i<size_inout_buffer; i++)
	{
		float	f;
		f =cu_get_exp(inout_buffer[i][0]);

		if ( f!=inout_buffer[i][1] )//if not match, hold the screen
		{
			printf("\n%dth get_exp not match!\n", i);
			printf("\nHOST:    (%f, %f)\n",inout_buffer[i][0],inout_buffer[i][1]);
			printf("\nDEVICE:  (%f, %f)\n",inout_buffer[i][0],f);
			getchar();
		}
	}
#endif//if CUDA_TEST_EXP

// -----------------------------------
// ----------- End Cuda --------------
// -----------------------------------
	cu_end();
	//delete things now
	if (d_particle_data)
		delete []d_particle_data;
#ifdef	CUDA_DEVICE_COMBINE
	if (cu_pic)
		delete	[]cu_pic;
#endif //ifdef CUDA_DEVICE_COMBINE
	//delete C style colormap 
	delete	[]amapD;
	delete	[]amapDTypeStartPos;
	//delete cu_particle_splotch objects
	delete	[]cu_ps;
	delete	[]cu_ps_filtered;
	//delete fragment buffer
	delete  []fragBuf;

#endif //ifdef CUDA
//////////////////////CUDA codes end here/////////////////////////////


/////////////////////////////////////////////////////////////////////
//Codes merge here whether using CUDA or not

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
void 
GoldCompareFBuf(cu_fragment_AeqE *goldBuf, cu_fragment_AeqE *buf, int n)
{
	for (int i=0; i<n; i++)
	{
		if ( abs(goldBuf[i].deltaR - buf[i].deltaR) < 1.0e-5 &&
			 abs(goldBuf[i].deltaG - buf[i].deltaG) < 1.0e-5 &&
			 abs(goldBuf[i].deltaB - buf[i].deltaB) < 1.0e-5
			 )
			 continue;

		//else
		cout << endl << "Gold compare fragment " << i << endl;
		cout << "gold" << endl;
		cout << goldBuf[i].deltaR << ",";
		cout << goldBuf[i].deltaG << ",";
		cout << goldBuf[i].deltaB << "," << endl;

		cout << "cuda" << endl;
		cout << buf[i].deltaR << ",";
		cout << buf[i].deltaG << ",";
		cout << buf[i].deltaB << endl;

		//hold screen for reading
		cout << endl << "Press q to quit, other to continue..." <<endl;
		char c;
		//cin >> c;
		c =getchar();
		if (c == 'q')
			break;
	}
}

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

//compare particle_splotche datas
void GoldCompareSData
(vector<particle_splotch> host_data, cu_particle_splotch* device_data)
{
	for (int i=0; i< host_data.size(); i++)
	{
		if ( host_data[i].x == device_data[i].x &&
			host_data[i].y == device_data[i].y &&
			host_data[i].r == device_data[i].r &&
			host_data[i].ro == device_data[i].ro 
			)
			continue;

		//else
		cout << endl << "Gold compare at particle " << i << endl;
		cout << "gold" << endl;
		cout << host_data[i].x << ",";
		cout << host_data[i].y << ",";
		cout << host_data[i].r << ",";
		cout << host_data[i].ro << ",";
		cout << "("<< host_data[i].a.r <<", ";
		cout << host_data[i].a.g <<", " <<host_data[i].a.b;
		cout << ")" << endl;
		cout << "("<< host_data[i].e.r <<", ";
		cout << host_data[i].e.g <<", " <<host_data[i].e.b;
		cout << ")" << endl;

		cout << "cuda" << endl;
		cout << device_data[i].x << ",";
		cout << device_data[i].y << ",";
		cout << device_data[i].r << ",";
		cout << device_data[i].ro << ",";
		cout << "("<< device_data[i].a.r <<", ";
		cout << device_data[i].a.g <<", " <<device_data[i].a.b;
		cout << ")" << endl;
		cout << "("<< device_data[i].e.r <<", ";
		cout << device_data[i].e.g <<", " <<device_data[i].e.b;
		cout << ")" << endl;


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
cu_color	C_get_color(int ptype, float val, cu_color_map_entry *map, //will move to kernel
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
