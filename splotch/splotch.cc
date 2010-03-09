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
#error  splotch: interpolation without geometry file makes no sense!
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

#if 0
// ----------------------------------------------
// ------- How to build Parameter structure ------
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
  c3=COLOUR(0,0,1);           // blue
  amap.addVal(0,c1);
  amap.addVal(0.5,c2));
  amap.addVal(1.,c3));
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
    if (linecount==geometry_skip)      // read only once if no interpolation is choosen
      {
#endif
      if (master)
        cout << endl << "reading data ..." << endl;
      int simtype = params.find<int>("simtype"); ///2:Gadget2
      float maxr, minr;
#ifdef INTERPOLATE
      double frac=(double)(linecount-(nextfile-ninterpol))/(double)ninterpol;
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
#ifdef INTERPOLATE          // Here only the tow datasets are prepared, interpolation will be done later
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
    if (master)
      (mpiMgr.num_ranks()>1) ?
        cout << endl << "applying local sort ..." << endl :
        cout << endl << "applying sort (" << npart << ") ..." << endl;
    int sort_type = params.find<int>("sort_type",1);
    particle_sort(particle_data,sort_type,true);
    wallTimer.stop("sort");

//2009-10-22: sorting is ignored now to simplified things.
//if sorted, it's needed to copy particle sim back to device
#ifdef CUDA
    if (sort_type!=0) //means needing sort
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
    vector<particle_sim>        tmp;
        int     n1,n2,n3;
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
/*      thread_info     ti;
        ti.startP=0;
        ti.endP =25000;
        ti.pPic =&pic;
        HANDLE h=CreateThread( NULL, 0,
                        (LPTHREAD_START_ROUTINE)host_thread_func,&ti, 0, NULL );
        WaitForSingleObject(h, INFINITE);
        goto out;
*/
        //new arrays of thread_info and HANDLE
        int     nDev;
        nDev =params.find<int>("gpu_number",0);
        //see if use host as a working thread
        bool bHostThread;
        bHostThread =params.find<bool>("use_host_as_thread", false);
        int nThread = bHostThread? nDev+1: nDev;
        //init objects for threads control
        thread_info     *tInfo =new thread_info[nThread];
#ifndef NO_WIN_THREAD
        HANDLE          *tHandle =new HANDLE[nThread];
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
//      host_thread_func ( &(tInfo[nDev]) );
#endif  //if not NO_WIN_THREAD

        //post-process
//      VTimer  timer;
//      timer.reset();
//      timer.start();
        //combine the results to pic
        if (1)//a_eq_e)
        {
                for (int i=1; i<nThread; i++)
                        for (int x=0; x<res; x++) //  error when x =1,
                                for (int y=0; y<res; y++)
                                        pic[x][y] =pic[x][y] + (*tInfo[i].pPic)[x][y];

        }
        else
        {} //to be done later...

        if(1)//a_eq_e)
        {
                exptable        xexp(MAX_EXP);
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
//      timer.stop();
//      cout << endl << "Post-process pic[] cost time:" << timer.getTime() << "s" <<endl;


        //now output the time records
        for (int i=0; i<nThread; i++)
        {
                if ( tInfo[i].devID!= -1)
                {
                        cout<< endl <<"Times of GPU" << i << ":" <<endl;
                        cout<< "CUDA_INIT:              " << tInfo[i].times[CUDA_INIT] <<endl;
                        cout<< "COPY2C_LIKE:            " << tInfo[i].times[COPY2C_LIKE] <<endl;
                        cout<< "RANGE:                  " << tInfo[i].times[RANGE] <<endl;
                        cout<< "TRANSFORMATION:         " << tInfo[i].times[TRANSFORMATION] <<endl;
                        cout<< "COLORIZE:               " << tInfo[i].times[COLORIZE] <<endl;
                        cout<< "FILTER:                 " << tInfo[i].times[FILTER] <<endl;
                        cout<< "SORT:                   " << tInfo[i].times[SORT] <<endl;
                        cout<< "RENDER:                 " << tInfo[i].times[RENDER] <<endl;
                        cout<< "COMBINE:                " << tInfo[i].times[COMBINE] <<endl;
                        cout<< "THIS_THREAD:            " << tInfo[i].times[THIS_THREAD] <<endl;
                        cout<<endl;
                }
                else
                {
                        cout<< endl <<"Times of CPU as a thread:" <<endl;
                        cout<< "RANGE:                  " << tInfo[i].times[RANGE] <<endl;
                        cout<< "TRANSFORMATION:         " << tInfo[i].times[TRANSFORMATION] <<endl;
                        cout<< "COLORIZE:               " << tInfo[i].times[COLORIZE] <<endl;
                        cout<< "RENDER:                 " << tInfo[i].times[RENDER] <<endl;
                        cout<< "THIS_THREAD:            " << tInfo[i].times[THIS_THREAD] <<endl;
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

#endif  //if def CUDA

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
  //Just to hold the screen to read the messages when debuging
  cout << endl << "Press any key to end..." ;
  getchar();
#endif
  }
