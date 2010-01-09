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
#ifdef USE_MPI
#include "mpi.h"
#else
#ifdef VS
#include "cuda/VTimer.h"
#else

#include <sys/time.h>
#endif
#endif

#include "splotch/splotchutils.h"
#include "writer/writer.h"
#include "reader/reader.h"
#include "cxxsupport/walltimer.h"

#ifdef VS
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/arr.h"
#include "cxxsupport/cxxutils.h"
#include "cxxsupport/paramfile.h"
#include "windows.h"
#include "reader/gadget_reader.cc"
#include "reader/bin_reader.cc"
#include "writer/write_tga.cc"
#endif

#ifdef CUDA
#include "cuda/splotch_cuda.h"
#include <string.h>
#endif

using namespace std;

#ifdef CUDA
#ifndef VS
#define DWORD long
#define WINAPI
#endif

//function definitions for cuda/testing use
void	GoldComparePData
(vector<particle_sim> particle_data, cu_particle_sim* d_particle_data);
void	GoldCompareSData
(vector<particle_splotch> host_data, cu_particle_splotch* device_data);
cu_color	C_get_color(int ptype, float val, cu_color_map_entry *map, //will move to kernel
			int	mapSize, int *ptype_points, int ptypes);
void	GoldCompareFBuf(cu_fragment_AeqE *goldBuf, cu_fragment_AeqE *buf, int n);

//things for combination with host threads
struct	param_combine_thread//for host combine thread
{
//	bool	bFinished; used in thread combine. not working.
	bool	a_eq_e;
	void	*fbuf;
	int		combineStartP, combineEndP;
	cu_particle_splotch	*ps;
	float	timeUsed;
	arr2<COLOUR>	*pPic;
};
#ifndef NO_WIN_THREAD
DWORD WINAPI combine(void	*param);
DWORD WINAPI TestThreadCombineTime(void	*p);
#else
DWORD WINAPI cu_thread_func(void *pinfo);
DWORD WINAPI cu_draw_chunk(void *pinfo);
#endif

//for record times
enum TimeRecords{
	CUDA_INIT,
	COPY2C_LIKE,
	RANGE,
	TRANSFORMATION,
	COLORIZE,
	FILTER,
	SORT,
	RENDER,
	COMBINE,
	THIS_THREAD,
	TIME_RECORDS,	//to indicate number of times
};

#ifdef CUDA_THREADS
	//struct containing thread task info
	struct thread_info{		
		int	devID;						//index of the device selected
		int	startP, endP;				//start and end particles to handle
		arr2	<COLOUR>	*pPic;		//the output image
		float	times[TIME_RECORDS];	//carry out times of computing
	};

	//some global info shared by all threads
	paramfile	*g_params;
	vector<particle_sim> particle_data; ///row data from file
	VECTOR campos, lookat, sky; ///A 3D vector class, designed for high efficiency.
	vector<COLOURMAP> amap,emap;
	int ptypes = 0;
//	arr2<COLOUR> *g_ppic;	//for testing only

	void	DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread);
#ifndef NO_WIN_THREAD
	DWORD	WINAPI	cu_draw_chunk(void *p);
	DWORD	WINAPI cu_thread_func(void *p);
	DWORD	WINAPI host_thread_func(void *p);
#endif


#endif //if CUDA_THREADS defined 

/*for temp testing
cu_color	pic1[800][800];//for host combine
cu_color	pic2[800][800];//for host rendering
struct region_cmp
{
	int operator()(const cu_particle_splotch &p1, const cu_particle_splotch &p2)
	{
		int	rgn1, rgn2;
		rgn1 =(p1.maxx -p1.minx)*(p1.maxy -p1.miny);
		rgn2 =(p2.maxx -p2.minx)*(p2.maxy -p2.miny);
		return rgn1>rgn2;
	}
};
*/
#endif //ifdef CUDA

int main (int argc, const char **argv)
{
  wallTimer.start("full");
  wallTimer.start("setup");
  bool master = mpiMgr.master(); ///return true

  module_startup ("splotch",argc,argv,2,"<parameter file>",master);
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

  //paramfile params (argv[1],master);
  paramfile params (argv[1],false);
#ifndef CUDA_THREADS
  vector<particle_sim> particle_data; ///row data from file
  VECTOR campos, lookat, sky; ///A 3D vector class, designed for high efficiency.
  vector<COLOURMAP> amap,emap;
  int ptypes = params.find<int>("ptypes",1); ///each particle type has a color map
#else //if CUDA_THREADS defined 
  ptypes = params.find<int>("ptypes",1); ///each particle type has a color map
  g_params =&params;
#endif	//if def CUDA_THREADS they will be a global vars

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
      if(linecount == nextfile)
	{
	  nextfile=linecount+ninterpol;
	  snr1=snr2;
	  snr2=snr2+1;
	}
#endif
    }

  while (getline(inp, line))
    {
      sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i",
	     &campos.x,&campos.y,&campos.z,
	     &lookat.x,&lookat.y,&lookat.z,
	     &sky.x,&sky.y,&sky.z,&ninterpol);
      if(master)
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
      if (linecount == geometry_skip)      // read only once if no interpolation is choosen
	{
#endif
	  if (master)
	    cout << endl << "reading data ..." << endl;
	  int simtype = params.find<int>("simtype"); ///2:Gadget2
	  float maxr, minr;
#ifdef INTERPOLATE
          double frac=(double)(linecount-(nextfile-ninterpol))/(double)ninterpol;
#endif
	  switch(simtype)
	    {
	    case 0:
	      bin_reader_tab(params,particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
	      break;
	    case 1: 
	      bin_reader_block(params,particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
	      break;
	    case 2: 
#ifdef INTERPOLATE          // Here only the tow datasets are prepared, interpolation will be done later
	      cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fac: " << frac << endl;
	      cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl; 
	      cout << " (pos: " << linecount << " , " << nextfile << " , " << ninterpol << ")" << endl; 
              if(snr1 == snr2_now)
		{
  	          cout << " old2 = new1 !" << endl; 
		  particle_data1=particle_data2;
		  snr1_now = snr1;
                  time1 = time2;
		}
	      if(snr1_now != snr1)
		{
  	          cout << " reading new1 " << snr1 << endl; 
		  gadget_reader(params,particle_data1,snr1,&time1);
		  snr1_now = snr1;
		}
	      if(snr2_now != snr2)
		{
  	          cout << " reading new2 " << snr2 << endl; 
		  gadget_reader(params,particle_data2,snr2,&time2);
		  snr2_now = snr2;
		}
#else
	      gadget_reader(params,particle_data,0,&time); ///vector<particle_sim> particle_data;
#ifdef GEOMETRY_FILE
	      p_orig = particle_data;
#endif
#endif
	      break;
	    case 3: //enzo_reader(params,particle_data);
	      break;
            case 4:
              gadget_millenium_reader(params,particle_data,0,&time);
              break;
            case 5:
#if defined(USE_MPI) && defined(USE_MPIIO)
              bin_reader_block_mpi(params,particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
#else
	      planck_fail("mpi reader not available in non MPI compiled version !");
#endif
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
      particle_interpolate(params,particle_data,particle_data1,particle_data2,frac,time1,time2);
#endif

      long npart=particle_data.size();
      long npart_all=npart;
      mpiMgr.allreduce (npart_all,MPI_Manager::Sum); ///does nothing
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
#ifdef USE_MPI
      if (master)
        cout << endl << "applying local sort ..." << endl;
#else
      cout << endl << "applying sort (" << npart << ") ..." << endl;
#endif
      int sort_type = params.find<int>("sort_type",1);
      particle_sort(particle_data,sort_type,true);
      wallTimer.stop("sort");

//2009-10-22: sorting is ignored now to simplified things.
//if sorted, it's needed to copy particle sim back to device
#ifdef CUDA 
	  if ( sort_type !=0) //means needing sort
		  cout<< endl << "CAUTION: SORTING NEEDED!"<<endl;
#endif

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


/////////////Here comes CUDA code//////////////////////////
// -----------------------------------
// ----------- Run Cuda -------------
// -----------------------------------
// After reading, ranging with device
#ifdef CUDA
	//test assign, should remove later
    vector<particle_sim>	tmp;
	int	n1,n2,n3;
	n1 =particle_data.size();
	tmp.assign(particle_data.begin(), particle_data.begin()+1000);
	n2 =particle_data.size();
	n3 =tmp.size();
	bool b;
	b =particle_data.empty();

	//prepare the parameters for the cuda thread
	//the final image
	int res = params.find<int>("resolution",200);
	arr2<COLOUR> pic(res,res);//, pic1(res,res);//pic1 for debug only

	//test host threading
/*	thread_info	ti;
	ti.startP=0;
	ti.endP =25000;
	ti.pPic =&pic;
	HANDLE h=CreateThread( NULL, 0, 
			(LPTHREAD_START_ROUTINE)host_thread_func,&ti, 0, NULL );
	WaitForSingleObject(h, INFINITE);
	goto out;
*/
	//new arrays of thread_info and HANDLE
	int	nDev;
	nDev =params.find<int>("gpu_number",0);
	//see if use host as a working thread
	bool bHostThread;
	bHostThread =params.find<bool>("use_host_as_thread", false);
	int nThread = bHostThread? nDev+1: nDev;
	//init objects for threads control
	thread_info	*tInfo =new thread_info[nThread];
#ifndef NO_WIN_THREAD
	HANDLE		*tHandle =new HANDLE[nThread];
#endif
	//fill in thread_info
	tInfo[0].pPic =&pic;
	for (int i=0; i<nDev; i++)
	{
		tInfo[i].devID =i;
		//local var pic is assigned to the first device
		if ( i!=0 )
			tInfo[i].pPic =new arr2<COLOUR>(res, res);
	}
	//make the last one work for host thread
	if (bHostThread )
	{
		tInfo[nThread-1].devID =-1;
		if ( nThread-1 != 0 )
			tInfo[nThread-1].pPic =new arr2<COLOUR>(res, res);
	}
	//decide how to devide task by another function
	DevideThreadsTasks(tInfo, nThread, bHostThread);

#ifndef NO_WIN_THREAD //to let it compiled in Linux, just for now, 2 Dec 2009.
	//issue the threads
	for (int i=0; i<nDev; i++)
		tHandle[i] =CreateThread( NULL, 0, 
			(LPTHREAD_START_ROUTINE)cu_thread_func,&(tInfo[i]), 0, NULL );
	//issue the host thread too
	if (bHostThread)
		tHandle[nDev] =CreateThread( NULL, 0, 
			(LPTHREAD_START_ROUTINE)host_thread_func,&(tInfo[nDev]), 0, NULL );
	//and wait for them to finish
	WaitForMultipleObjects(nThread, tHandle, true, INFINITE);

#else //do not use thread which is now Windows code
	cu_thread_func (&(tInfo[0])); //just call it as normal function
//	host_thread_func ( &(tInfo[nDev]) );
#endif	//if not NO_WIN_THREAD

	//post-process
//	VTimer	timer;
//	timer.reset();
//	timer.start();
	//combine the results to pic
	if (1)//a_eq_e)
	{
		for (int i=1; i<nThread; i++)
			for (int x=0; x<res; x++) //  error when x =1,
				for (int y=0; y<res; y++)
					pic[x][y] =pic[x][y] + (*tInfo[i].pPic)[x][y];

	}
	else
		;//to be done later...

	if(1)//a_eq_e)
	{
		exptable	xexp(MAX_EXP);
		for(int ix =0; ix <pic.size1(); ix++)
			for(int iy =0; iy <pic.size2(); iy++)
			  {	
				  pic[ix][iy].r=1-xexp(pic[ix][iy].r);
				  pic[ix][iy].g=1-xexp(pic[ix][iy].g);
				  pic[ix][iy].b=1-xexp(pic[ix][iy].b);
			  }
	}
	else
	{
	}
//	timer.stop();
//	cout << endl << "Post-process pic[] cost time:" << timer.getTime() << "s" <<endl;


	//now output the time records
	for (int i=0; i<nThread; i++)
	{
		if ( tInfo[i].devID!= -1)
		{
			cout<< endl <<"Times of GPU" << i << ":" <<endl;
			cout<< "CUDA_INIT:		" << tInfo[i].times[CUDA_INIT] <<endl;
			cout<< "COPY2C_LIKE:		" << tInfo[i].times[COPY2C_LIKE] <<endl;
			cout<< "RANGE:			" << tInfo[i].times[RANGE] <<endl;
			cout<< "TRANSFORMATION:		" << tInfo[i].times[TRANSFORMATION] <<endl;
			cout<< "COLORIZE:		" << tInfo[i].times[COLORIZE] <<endl;
			cout<< "FILTER:			" << tInfo[i].times[FILTER] <<endl;
			cout<< "SORT:			" << tInfo[i].times[SORT] <<endl;
			cout<< "RENDER:			" << tInfo[i].times[RENDER] <<endl;
			cout<< "COMBINE:		" << tInfo[i].times[COMBINE] <<endl;
			cout<< "THIS_THREAD:		" << tInfo[i].times[THIS_THREAD] <<endl;
			cout<<endl;
		}
		else
		{
			cout<< endl <<"Times of CPU as a thread:" <<endl;
			cout<< "RANGE:			" << tInfo[i].times[RANGE] <<endl;
			cout<< "TRANSFORMATION:		" << tInfo[i].times[TRANSFORMATION] <<endl;
			cout<< "COLORIZE:		" << tInfo[i].times[COLORIZE] <<endl;
			cout<< "RENDER:			" << tInfo[i].times[RENDER] <<endl;
			cout<< "THIS_THREAD:		" << tInfo[i].times[THIS_THREAD] <<endl;
			cout<<endl;
		}
	}

	//delete pics that were created
	for (int i=1; i<nThread; i++)
	  delete tInfo[i].pPic;
	//delete thread_info and HANDLE arrays
	delete [] tInfo;
#ifndef NO_WIN_THREAD
	delete [] tHandle;
#endif

#endif	//if def CUDA

//18-11-2009, old code here went to save_cuda_in_main_before_multithread.cpp

//old code
/////////NO CUDA WE TEST THE TIME FOR COMBINE
//////////////////////CUDA codes end here/////////////////////////////


/////////////////////////////////////////////////////////////////////
//Codes merge here whether using CUDA or not

      wallTimer.start("write");
	
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
#ifdef INTERPOLATE
      if(linecount == nextfile)
	{
	  nextfile=linecount+ninterpol;
	  snr1=snr2;
	  snr2=snr2+1;
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
      for(int i=1; i<geometry_incr; i++)
	{
	  getline(inp, line);
	  linecount++;
#ifdef INTERPOLATE
	  sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i",
		 &campos.x,&campos.y,&campos.z,
		 &lookat.x,&lookat.y,&lookat.z,
		 &sky.x,&sky.y,&sky.z,&ninterpol);
	  if(linecount == nextfile)
	    {
	      nextfile=linecount+ninterpol;
	      snr1=snr2;
	      snr2=snr2+1;
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
  //Just to hold the screen to read the messages when debuging
  cout << endl << "Press any key to end..." ;
  getchar();
#endif

}


////////////////CUDA FUNCTIONS AND HELPER FUNCTIONS ////////////////////////////
#ifdef CUDA
DWORD WINAPI host_thread_func(void *p)
{
	printf("\nHost Thread Start!\n");

	thread_info	*tInfo = (thread_info*)p;

	paramfile	params =*g_params;

	vector<particle_sim>	particles;
	vector<particle_sim>::iterator i1,i2;
	i1 =particle_data.begin() +tInfo->startP;
	i2 =particle_data.begin() +tInfo->endP;
	particles.assign(i1, i2);

//	VTimer t, t1;
//	t1.start();
//	t.start();

	particle_normalize(params,particles,false); 

//	t.stop();
//	tInfo->times[RANGE] =t.getTime();
//	t.reset();
//	t.start();

//	vector<particle_splotch> particle_col; not used anymore since Dec 09.
        particle_project(params, particles, campos, lookat, sky);

//	t.stop();
//	tInfo->times[TRANSFORMATION] =t.getTime();

	// ----------- Sorting ------------
	// NO SORING FOR NOW
	//int sort_type = params.find<int>("sort_type",1);
	//particle_sort(particle_data,sort_type,true);

//	t.reset();
//	t.start();

//	particle_colorize(params, particles, particle_col, amap, emap);
	particle_colorize(params, particles, amap, emap); //new calling

//	t.stop();
//	tInfo->times[COLORIZE]=t.getTime();
//	t.reset();
//	t.start();

	int res = params.find<int>("resolution",200);
	float64 grayabsorb = params.find<float>("gray_absorption",0.2);
	bool a_eq_e = params.find<bool>("a_eq_e",true);
//	render_as_thread(particle_col,*(tInfo->pPic),a_eq_e,grayabsorb);
	render_as_thread1(particles,*(tInfo->pPic),a_eq_e,grayabsorb);

//	t.stop();
//	tInfo->times[RENDER] =t.getTime();
//	t1.stop();
//	tInfo->times[THIS_THREAD] =t1.getTime();

	printf("\nHost Thread End!\n");
	return 1;
}


DWORD WINAPI combine(void	*param1)
{
	static	int enter=-1;

	param_combine_thread	*param =(param_combine_thread*)param1;

//	printf("\ncombine in  %d, %d-%d", ++enter, param->combineStartP, param->combineEndP);

	//combine it
	cu_particle_splotch	*ps	=param->ps;
	arr2<COLOUR>	*pPic =param->pPic;

//	VTimer t;
//	t.start();

	if (param->a_eq_e)
	{
		cu_fragment_AeqE	*bufWrite =(cu_fragment_AeqE*)param->fbuf;

		for (int pPos=param->combineStartP, fPos=0; 
			pPos<param-> combineEndP; pPos++)
		{
		  for (int x =ps[pPos].minx; x <ps[pPos].maxx; x++)
		  {
			 for( int y =ps[pPos].miny; y <ps[pPos].maxy; y++)
			 {
				 (*pPic)[x][y].r += bufWrite[fPos].deltaR;
				 (*pPic)[x][y].g += bufWrite[fPos].deltaG;
				 (*pPic)[x][y].b += bufWrite[fPos].deltaB;
				 fPos ++;
			 }				  
		  }
		}
	}
	else
	{
		cu_fragment_AneqE	*bufWrite =(cu_fragment_AneqE*)param->fbuf;

		for (int pPos=param->combineStartP, fPos=0; 
			pPos<param-> combineEndP; pPos++)
		{
		  for (int x =ps[pPos].minx; x <ps[pPos].maxx; x++)
		  {
			 for( int y =ps[pPos].miny; y <ps[pPos].maxy; y++)
			 {
#ifndef	CUDA_THREAD
				 (*pPic)[x][y].r = (*pPic)[x][y].r *bufWrite[fPos].factorR +bufWrite[fPos].deltaR;
				 (*pPic)[x][y].g = (*pPic)[x][y].g *bufWrite[fPos].factorG +bufWrite[fPos].deltaG;
				 (*pPic)[x][y].b = (*pPic)[x][y].b *bufWrite[fPos].factorB +bufWrite[fPos].deltaB;
#else	//if def CUDA_THREAD

#endif	//ifndef CUDA_THREAD
				 fPos ++;
			 }				  
		  }
		}
	}

//	t.stop();
//	param->timeUsed +=t.getTime();

//	printf("\ncombine out %d", enter);
	return 1;
}

DWORD WINAPI cu_thread_func(void *pinfo)
{
//	VTimer	timer, timer1;
	//a new thread info object that will carry each chunk's drawing
	thread_info	ti, *pInfoOutput;
	pInfoOutput =(thread_info*) pinfo;
	ti =*pInfoOutput;

	//do some cleaning for final thread_info
	pInfoOutput->pPic->fill(COLOUR(0.0, 0.0, 0.0));
	memset(pInfoOutput->times, 0, sizeof(float)*TIME_RECORDS);

	//a new pic object residing in ti that will carry the result
	arr2<COLOUR> pic(pInfoOutput->pPic->size1(), pInfoOutput->pPic->size2());
	ti.pPic =&pic;

	//set startP and end P of ti
	int len = cu_get_chunk_particle_count(*g_params);
	if ( len==-1)
	{
		printf("\nGraphics memory setting error\n");
		return -1;
	}

	int	curEnd=0;
	int endP =ti.endP;
	ti.endP =ti.startP;

#ifdef DEBUG
   cout << "cu_thread_func1\n";
#endif
	while( ti.endP<endP )	
	{
		//set range
		ti.endP =ti.startP +len -1;
		if (ti.endP >endP)
			ti.endP =endP;

		//draw chunks one by one
		cu_draw_chunk(&ti);
		//collect image to result
#ifdef DEBUG
   cout << "ti.endP " << ti.endP<< "\n";
   cout << "pic.size1() " << pic.size1()<< "\n";
   cout << "pic.size2() " << pic.size2()<< "\n";
#endif
		for (int x=0; x<pic.size1(); x++)
		  for (int y=0; y<pic.size2(); y++)
		    (*(pInfoOutput->pPic))[x][y] =  (*(pInfoOutput->pPic))[x][y] + pic[x][y];
		//collect times to output
		for (int i=0; i<TIME_RECORDS; i++)
			pInfoOutput->times[i] +=ti.times[i];

		//set range
		ti.startP =ti.endP +1;
	}

	//test 2-goes pased...

#ifdef DEBUG
   cout << "cu_thread_func2\n";
#endif
	return 1;
}

DWORD WINAPI cu_draw_chunk(void *pinfo)
{
//	VTimer	timer, timer1;
	float	time;

//	timer1.reset(); //for the whole thread
//	timer1.start();

	//get the input info
	thread_info	*tInfo = (thread_info*)pinfo;	//	if (!tInfo->devID) return 0;
	int	nParticle =tInfo->endP -tInfo->startP +1;
	//prepare for recording times
//	memset(tInfo->times, 0, sizeof(float)*TIME_RECORDS);

	paramfile	params =*g_params;
	//CUDA test. for developing only. cut short particle_data
	vector<particle_sim>::iterator it;
	int	testPCount =10000, n=0;
	it =particle_data.begin();
//	particle_data.erase(it+testPCount, particle_data.end());

	cu_particle_sim	*d_particle_data;
//	d_particle_data =new cu_particle_sim[particle_data.size()];
	d_particle_data =new cu_particle_sim[nParticle];

	//for each gpu/thread a varible pack is needed
	cu_gpu_vars	gv;
	memset(&gv, 0, sizeof(cu_gpu_vars) );

	//CUDA Init 
//	timer.reset();
//	timer.start();
	cu_init(params, tInfo->devID, &gv);
//	timer.stop();
//	time =timer.getTime();
//	cout << endl << "cu_init() cost time:" << time << "s" <<endl;
//	tInfo->times[CUDA_INIT]=time;

	//Copy particle sim into C-style object d_particle_data
//	timer.reset();
//	timer.start();
	//copy data to local C-like array d_particle_data, in mid of developing only
	for (int i=tInfo->startP,j=0; i<=tInfo->endP; i++, j++)
		memcpy( &(d_particle_data[j]), &(particle_data[i]), sizeof(cu_particle_sim));
//	timer.stop();
//	time =timer.getTime();
//	cout << endl << "Copying particles to device cost time:" << time << "s" <<endl;
	tInfo->times[COPY2C_LIKE]=time;

	//here we analysis how to divide the whole task for large data handling

	//CUDA Ranging
//	timer.reset();
//	timer.start();
	//call cuda range
	cu_range(params, d_particle_data, nParticle, &gv);
//	timer.stop();
//	time =timer.getTime();
//	cout << endl << "Ranging with device cost time:" << time << "s" <<endl;
	tInfo->times[RANGE]=time;

	//old code then copy particles back, in mid of developing only
	//CUDA RANGING DONE!

	//CUDA Transformation
//	timer.reset();
//	timer.start();
	double	c[3]={campos.x, campos.y, campos.z}, 
			l[3]={lookat.x, lookat.y, lookat.z}, 
			s[3]={sky.x, sky.y,	sky.z};
	cu_transform(params, nParticle,c, l, s,d_particle_data, &gv);
//	timer.stop();
//	time =timer.getTime();
//	cout << endl << "Transforming with device cost time:" << time<< "s" <<endl;
//	tInfo->times[TRANSFORMATION]=time;

/*	temporarily ignore sorting 191109.
	it becomes comlicated when using multiple threads with sorting
	//then copy particles back to host for sorting
	for (int i=0; i<nParticle; i++)
		memcpy( &(particle_data[iWRONG]),&(d_particle_data[i]), sizeof(cu_particle_sim));
PROBLEM HERE!

// --------------------------------
// ----------- Sorting ------------
// --------------------------------
//      cout << endl << "applying sort (" << npart << ") ..." << endl;

      int sort_type = params.find<int>("sort_type",1);
      particle_sort(particle_data,sort_type,true);
	  //we can sort by size(r) for balancing with cuda
	  //particle_sort(particle_data,4,true);

	  //copy sorted data back to device, not a must!
	  //first to C-style object
	  for(int i=0; i<particle_data.size(); i++)
		memcpy( &(d_particle_data[i]), &(particle_data[i]), sizeof(cu_particle_sim));
	  cu_copy_particle_sim_to_device(d_particle_data, particle_data.size());
*/

	  //CUDA Coloring
//	timer.reset();
//	timer.start();
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
        cout << "E.SIZE ..... " << e.size() << "\n";
		for (j=0; j<e.size(); j++)
		{
			HANDLE_RAYPP<CMAP_ENTRY> h;
                        cout << "--->><< e[j]      " << e[j]->minval << "\n";
                        cout << "--->><< e[j]      " << e[j]->maxval << "\n";
			h =e[j];
			COLOUR	clr1, clr2;
			amapD[index].min =h->minval; 
			amapD[index].max =h->maxval;
			//clr1 =h->Get_Colour(h->minval);
			//clr2 =h->Get_Colour(h->maxval);
                        cout <<  clr1.r <<"<<<<<<----\n";
                        cout << h->Get_Colour(0.5).r  <<"<<<<<<----\n";
                        cout << h->Get_Colour(1).g << " " << clr1.g <<"<<<<<<----\n";
                        cout << h->Get_Colour(1).b << " " << clr1.b <<"<<<<<<----\n";
                        clr1.r =h->Get_Colour(h->minval).r;
                        clr1.g =h->Get_Colour(h->minval).g;
                        clr1.b =h->Get_Colour(h->minval).b;
                        clr2.r =h->Get_Colour(h->maxval).r;
                        clr2.g =h->Get_Colour(h->maxval).g;
                        clr2.b =h->Get_Colour(h->maxval).b;
			amapD[index].color1.r =clr1.r;
			amapD[index].color1.g =clr1.g;
			amapD[index].color1.b =clr1.b;
			amapD[index].color2.r =clr2.r;
			amapD[index].color2.g =clr2.g;
			amapD[index].color2.b =clr2.b;
                        cout << "pippo\n";
			index++;
		}
		amapDTypeStartPos[i] =curPtypeStartPos;
		curPtypeStartPos += j;
	}
	//now let cuda init colormap on device
	cu_colormap_info clmp_info;
	clmp_info.map =amapD;
	clmp_info.mapSize =size;
	clmp_info.ptype_points =amapDTypeStartPos;
	clmp_info.ptypes =ptypes;
	cu_init_colormap(clmp_info, &gv);

	//old code test color map
	
	//init cu_particle_splotch array memeory
	cu_particle_splotch	*cu_ps;
//	size =particle_data.size();
	size =nParticle;
	cu_ps =new cu_particle_splotch[size];
	memset(cu_ps, 0, size);

	//Colorize with device
	cu_colorize(params, cu_ps, size, &gv);
//      timer.stop();
//      time =timer.getTime();
//	cout << endl << "Coloring with device cost time:" << time << "s" <<endl;
	tInfo->times[COLORIZE]=time;


//////////////////////////////////////////////////////////////////////////////////
//      timer.reset();
//      timer.start();
	//filter particle_splotch array to a cu_ps_filtered
	int	pFiltered=0;
	//for sorting and splitting
	vector <cu_particle_splotch> v_ps;	
	cu_particle_splotch	p;
	//do filtering
	unsigned long	posInFragBuf =0;//, countFragments;
	int		minx=1e6,miny=1e6, maxx=-1,maxy=-1;

//old code observ size
	//selecte valid ones
	for (int i=0; i<size; i++)
	{
		if ( cu_ps[i].isValid )
		{
			//IMPORTANT: set the start position of this particle
			//in fragment buffer
//			cu_ps[i].posInFragBuf =posInFragBuf;
//			posInFragBuf += (cu_ps[i].maxx -cu_ps[i].minx)*
//				(cu_ps[i].maxy -cu_ps[i].miny); SHOULD DO AFTER SORTING!
//			memcpy(&cu_ps_filtered[pFiltered], 
//				&cu_ps[i], sizeof(cu_particle_splotch));
                        
		        memcpy(&p, &cu_ps[i], sizeof(cu_particle_splotch));
			v_ps.push_back(p);
			pFiltered++;

			minx=min(minx,(int)cu_ps[i].minx);
			miny=min(miny,(int)cu_ps[i].miny);
			maxx=max(maxx,(int)cu_ps[i].maxx);
			maxy=max(maxy,(int)cu_ps[i].maxy);
//old code observ size
		}
	}
//old code observ size
//      timer.stop();
//      time =timer.getTime();
//	cout << endl << "Filtering 1 costs time:" << time << "s" <<endl;
	tInfo->times[FILTER] =time;

//      timer.reset();
//      timer.start();

	//split large ones
	int	maxRegion =cu_get_max_region(&gv);
	vector<cu_particle_splotch> v_ps1;//result goes to it
	for (vector<cu_particle_splotch>::iterator i=v_ps.begin();
		i<v_ps.end(); i++)
	{
		int	h, w;
		cu_particle_splotch	p, pNew;
		p =(*i);
		h = p.maxy - p.miny;
		w = p.maxx - p.minx;

		if (h*w <maxRegion)
//		if ( 1)//no splicting test 
		{
			v_ps1.push_back(p);
			continue;
		}

		//now we split
		int	w1;
		w1 = (maxRegion %h==0) ? (maxRegion /h):(maxRegion /h +1);
		//insert new cells
		pNew =p;
		//minx,maxx of pNew need to be set
		for (int minx =p.minx; minx<p.maxx; minx+=w1)
		{
			pNew.minx =minx;
			pNew.maxx =( minx+w1 >=p.maxx) ? p.maxx : minx+w1;
			v_ps1.push_back(pNew);
		}
	}
//      timer.stop();
//      time =timer.getTime();
	//cout << endl << "Filtering 2 costs time:" << time << "s" <<endl;
	tInfo->times[FILTER] +=time;

//      timer.reset();
//      timer.start();
	v_ps.clear();//not useful any more
	//sort the filtered,splitted v_ps
//	sort(v_ps1.begin(), v_ps1.end(), region_cmp());
//      timer.stop();
//      time =timer.getTime();
	//cout << endl << "Sorting costs time:" << time << "s" <<endl;
	tInfo->times[SORT]=time;

//      timer.reset();
//      timer.start();
	//copy to C-style array cu_ps_filtered
	cu_particle_splotch	*cu_ps_filtered;
	size =v_ps1.size();
	cu_ps_filtered =new cu_particle_splotch[size];
	pFiltered =v_ps1.size();
	for(int i=0; i<pFiltered; i++)
	{
		v_ps1[i].posInFragBuf =posInFragBuf;
		int	region =(v_ps1[i].maxx -v_ps1[i].minx)*(v_ps1[i].maxy -v_ps1[i].miny); 
		posInFragBuf +=region;
		cu_ps_filtered[i] =v_ps1[i];
	}
//      timer.stop();
//      time =timer.getTime();
	//cout << endl << "Filtering 3 costs time:" << time << "s" <<endl;
	tInfo->times[FILTER] +=time;


	//old code
	//do gold comparation
	//debug only: then copy data back to host p2 and continue
	//working only with NO_HOST_COLORING
	//here is the point that particle sim array on deivce is useless

	//just for rendering time coming next

	//old codecompare to gold result, test only

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------

	//get parameters for rendering
	int res = params.find<int>("resolution",200);
	long nsplotch=pFiltered;
	long nsplotch_all=nsplotch;
	mpiMgr.allreduce_sum (nsplotch_all);
//	if (master)
//		cout << endl << "rendering (" << nsplotch_all << "/" << npart_all << ")..." << endl;
//	arr2<COLOUR> pic(res,res);
//	arr2<COLOUR> pic =*g_ppic;
	float64 grayabsorb = params.find<float>("gray_absorption",0.2);
	bool a_eq_e = params.find<bool>("a_eq_e",true);


//old code test fragment buffer

	//CUDA Rendering with device

//here's the point of multi-go loop starts
#ifndef CUDA_DEVICE_COMBINE //combined by host
//      timer.reset();//for rendering
//      timer.start();
	//prepare fragment buffer memory space first
	cu_fragment_AeqE	*fragBufAeqE;
	cu_fragment_AneqE	*fragBufAneqE;
	int	nFBufInByte =cu_get_fbuf_size(&gv);
	int	nFBufInCell;
	if (a_eq_e)
	{
		nFBufInCell =nFBufInByte/sizeof(cu_fragment_AeqE);
		fragBufAeqE =new cu_fragment_AeqE[nFBufInCell];
	}
	else
	{
		nFBufInCell =nFBufInByte/sizeof(cu_fragment_AneqE);
		fragBufAneqE =new cu_fragment_AneqE[nFBufInCell];
	}
	


	cu_prepare_render(cu_ps_filtered,pFiltered, &gv);

//old code test AUTOMIC_ADD
	//clear the output pic
	tInfo->pPic ->fill(COLOUR(0.0, 0.0, 0.0));

	//initialize combine vars
	//cu_ps_filtered:	the particle array
	//pFiltered:		length of particle array
	//pPosInFragBuf:	length of the whole fragment buffer counted after filterign
	bool bFinished=false;
	int	 renderStartP, renderEndP, combineStartP, combineEndP;
	renderStartP =renderEndP =0;
	combineStartP = combineEndP=0;
	int	 nFragments2Render=0;
	//some timers
//	VTimer	t1, t2, t3;
//	t1.reset();
//	t2.reset();
//	memset(pic1, 0, sizeof(cu_color)*800*800); this is just for temp


	//now prepare the parameters for combination
	param_combine_thread	param;
	param.pPic =tInfo->pPic;
	param.a_eq_e =a_eq_e;
	if (a_eq_e)
		param.fbuf =(void*)fragBufAeqE;
	else
		param.fbuf =(void*)fragBufAneqE;
	param.ps =cu_ps_filtered;
	param.timeUsed =0.0;

//old code test a_eq_e

//#define HOST_THREAD_RENDER
#ifdef HOST_THREAD_RENDER
	bFinished=true;//let device not working

	param_render_thread	param_render;
	param_render.p	=cu_ps_filtered;
	param_render.start =0;
	param_render.end   =pFiltered;
	param_render.a_eq_e =a_eq_e;
	param_render.grayabsorb =grayabsorb;
	param_render.pic	=(cu_color [][800])pic1;

	render_thread(&param_render);
#endif


	while (!bFinished)
	{
		//find a chunk
		nFragments2Render =0;
		for (renderStartP =renderEndP; renderEndP<pFiltered; renderEndP++)
		{
			int	inc;//increasment
			inc =(cu_ps_filtered[renderEndP].maxx-cu_ps_filtered[renderEndP].minx) *
				(cu_ps_filtered[renderEndP].maxy -cu_ps_filtered[renderEndP].miny);
			if ( nFragments2Render +inc > nFBufInCell)
				break;
			else
				nFragments2Render +=inc;
		}
		if (renderEndP == pFiltered)
			bFinished =true;

		//render it
//		printf("\nThread%d, cu_render1, %d-%d of %d", tInfo->devID, renderStartP, renderEndP, pFiltered);
		cu_render1(renderStartP, renderEndP, 
			a_eq_e, grayabsorb, &gv); 

		//see if it's the first chunk
		if ( renderStartP!=0 )
		{
			//combine it 
//			printf("\nThread%d, combine, %d-%d", tInfo->devID, combineStartP, combineEndP);			
			param.combineStartP =combineStartP;
			param.combineEndP =combineEndP;
			DWORD	id;
			combine(&param);
//			CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)combine, 
//			    &param, 0, &id );

		}		

		//collect result
		cu_get_fbuf(fragBufAeqE, a_eq_e, nFragments2Render, &gv);
		combineStartP=renderStartP;
		combineEndP =renderEndP;

		//see if last chunk
		if (bFinished)
		{
//			printf("\nThread%d, last chunk to combine",tInfo->devID);
//			printf("\nThread%d, combine, %d-%d",tInfo->devID, combineStartP, combineEndP);
			param.combineStartP =combineStartP;
			param.combineEndP =combineEndP;
			combine(&param);
			// param.timeUsed +=t3.getTime();
		}
	}

//      timer.stop();
//      time =timer.getTime();
//	cout << endl << "Render to fragment buffer cost time:" << time << "s" <<endl;
//	cout << endl << "Time for combination:" << param.timeUsed <<endl;
	tInfo->times[RENDER]=time;
	tInfo->times[COMBINE]=param.timeUsed;

//....
//give a test for combination without device involved

/////////////////////////////////////////////////////////////////////

#endif //if not def CUDA_DEVICE_COMBINE
	
//old code device combine

//old code test exp
// -----------------------------------
// ----------- End Cuda --------------
// -----------------------------------
	cu_end(&gv);
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
	if (a_eq_e)
		delete  []fragBufAeqE;
	else
		delete  []fragBufAneqE;

	printf("\nThread %d finished!\n", tInfo->devID);

//	tInfo->times[THIS_THREAD] =timer1.getTime();

	return 1;
}

void	DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread)
{
	bool	bTestLoadBalancing;
	bTestLoadBalancing =g_params->find<bool>("test_load_balancing", false);
	int		hostLoad =g_params->find<int>("host_load",0);

	int	curStart =0, onePercent, averageDevLen, nDev;
	nDev =bHostThread? nThread-1: nThread;
	onePercent =particle_data.size()/100;
	averageDevLen = (nDev!=0)? onePercent *(100-hostLoad)/nDev :
		0;

	for (int i=0; i<nThread; i++)
	{
		tInfo[i].startP =curStart;
		
		if (tInfo[i].devID != -1) //not a host
		{
			if ( bTestLoadBalancing )
			{
				int	gpuLoad;
				gpuLoad =g_params->find<int>("gpu_load"+dataToString(i),0);
				tInfo[i].endP =curStart +gpuLoad* onePercent;
			}
			else
				tInfo[i].endP =curStart +averageDevLen;
		}
		else //if this is a host
		{
			tInfo[i].endP =curStart +hostLoad *onePercent;
		}

		curStart =tInfo[i].endP +1;
	}

	tInfo[nThread-1].endP =particle_data.size()-1;
}

DWORD WINAPI TestThreadCombineTime(void	*p)
{
//	VTimer	t;
//      t.start();
	for(unsigned int i=0; i<545211005; i++)
	{
		int j;
		j=i+1;
//		pic1[0][0].r +=i;
//		pic1[0][0].g +=i;
//		pic1[0][0].b +=i;
	}
//      t.stop();
//	cout << endl << "Time for thread combination:" << t.getTime() <<endl;

	return 1;
}

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
//////////////////////////////////////////////////////////////////////
