/*
 * Copyright (c) 2010-2014
 *              Marzia Rivi (1), Tim Dykes (2)
 *               (1) University of Oxford
 *               (2) University of Portsmouth
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

#include "cuda/cuda_splotch.h"

using namespace std;

#ifdef SPLOTCH_PARAVIEW
//
// paraview version of rendering sets up colour map independently
// and leaves particles on GPU between frames if they have not been modified
//
cu_gpu_vars gv = {0};
long int LEN   = 0;
int nTasksDev  = 1;
int    mydevID = 0;

int cuda_paraview_init(arr2<COLOUR> &pic, vector<particle_sim> &particle, const vec3 &campos, const vec3 &lookat, vec3 &sky, float b_brightness, paramfile &g_params)
{
  long int nP = particle.size();
  pic.fill(COLOUR(0.0, 0.0, 0.0));
  int xres = pic.size1();
  int yres = pic.size2();
 // cout << "resolution = " << xres << " x " << yres << endl;
  int ptypes = g_params.find<int>("ptypes",1);

  // CUDA Init
  // if data was allocated previously, delete it
  if (gv.d_pd!=NULL) {
    cu_end(&gv);
  }
  memset(&gv, 0, sizeof(cu_gpu_vars));
  // Initialize policy class
  CuPolicy *policy = new CuPolicy(xres, yres, g_params); 
  gv.policy = policy;

/*
REMOVED TILED IMPLEMENTATION
#ifndef CUDA_FULL_ATOMICS
  int ntiles = policy->GetNumTiles();
#else
*/
  // For atomic implementation just a placeholder
  int ntiles = 0;
  /*
#endif
  */

  // num particles to manage at once
  float factor = g_params.find<float>("particle_mem_factor", 4);
  LEN = cu_paraview_get_chunk_particle_count(&gv, nTasksDev, sizeof(cu_particle_sim), ntiles, factor, nP);
  if (LEN <= 0)
    {
    cout << "Graphics memory setting error" << endl;
    MPI_Manager::GetInstance()->abort();
    }

  // Set cuda paraview param
  g_params.setParam("cuda_paraview_splotch", true);

  // enable device and allocate arrays
  int error=0;
  tstack_push("Device setup");
  error = cu_init(mydevID, LEN, ntiles, &gv);
  tstack_pop("Device setup");

  return error;
}

void cuda_paraview_rendering(int mydevID, int nTasksDev, arr2<COLOUR> &pic, vector<particle_sim> &particle, const vec3 &campos, const vec3 &lookat, vec3 &sky, float b_brightness, paramfile &g_params, void *gpudata)
{
  tstack_push("CUDA");
  
  int ntiles = gv.policy->GetNumTiles();
  int xres = pic.size1();
  int yres = pic.size2();
  long int nP = particle.size();
  bool doLogs = true;

  tstack_push("Device setup");
  int error = cu_init_params(&gv, g_params, campos, lookat, sky, b_brightness, doLogs);
  tstack_pop("Device setup");
  {
    //a new linear pic object that will carry the result
    float64 grayabsorb = g_params.find<float>("gray_absorption",0.2);
    bool a_eq_e = g_params.find<bool>("a_eq_e",true);
 
    int endP = 0;
    int startP = 0;
    int nPR = 0;

    // We only render a single chunk in the paraview version.
    endP = startP + LEN;   //set range
    if (endP > nP) endP = nP;
    nPR += cu_draw_chunk(mydevID, (cu_particle_sim *) &(particle[startP]), endP-startP, pic, &gv, a_eq_e, grayabsorb, xres, yres, doLogs, gpudata);
    // combine host results of chunks
    cout << "Rank " << MPI_Manager::GetInstance()->rank() << ": Rendered " << nPR << "/" << nP << " particles" << endl << endl;
    startP = endP;
    if (endP<nP) {
      cout << "Not enough memory on GPU to render all particles " << std::endl;
    }

    add_device_image(pic, &gv, xres, yres);
    tstack_pop("CUDA");
//     cu_end(&gv);
  }

 }

#endif
 
void cuda_rendering(int mydevID, int nTasksDev, arr2<COLOUR> &pic, vector<particle_sim> &particle, const vec3 &campos, const vec3 &lookat, vec3 &sky, vector<COLOURMAP> &amap, float b_brightness, paramfile &g_params)
{
  tstack_push("CUDA");
  tstack_push("Device setup");
  long int nP = particle.size();

  // Fill our output image with black
  pic.fill(COLOUR(0.0, 0.0, 0.0));
  int xres = pic.size1();
  int yres = pic.size2();

/*
REMOVED TILED IMPLEMENTATION
#ifndef CUDA_FULL_ATOMICS
  // Create our host image buffer - to be incrementally filled by looped draw call
  arr2<COLOUR> Pic_host(xres,yres);
#else
*/
  // For atomic implementation just a placeholder
  arr2<COLOUR> Pic_host;
  /*
#endif
*/
  // CUDA Init
  // Initialize policy class
  CuPolicy *policy = new CuPolicy(xres, yres, g_params); 

/*
REMOVED TILED IMPLEMENTATION
#ifndef CUDA_FULL_ATOMICS
  int ntiles = policy->GetNumTiles();
#else
*/
  // For atomic implementation just a placeholder
  int ntiles = 0;
  /*
#endif
*/

  // Initialise struct to hold gpu-destined variables
  cu_gpu_vars gv;
  memset(&gv, 0, sizeof(cu_gpu_vars));
  gv.policy = policy;

   // Setup gpu colormap
  int ptypes = g_params.find<int>("ptypes",1);
  setup_colormap(ptypes, amap, &gv);

  // Calculate how many particles to manage at once dependant on GPU memory size
  // Particle mem factor modifies how much memory we request
  float factor = g_params.find<float>("particle_mem_factor", 4);
  long int len = cu_get_chunk_particle_count(&gv, nTasksDev, sizeof(cu_particle_sim), ntiles, factor);
  if (len <= 0)
  {
    cout << "Graphics memory setting error" << endl;
    // Splotch paraview doesnt use mpimgr
    #ifndef SPLOTCH_PARAVIEW
    mpiMgr.abort();
    #endif
  }

  // Initialise device and allocate arrays
  bool doLogs;
  int error = cu_init(mydevID, len, ntiles, &gv);
  tstack_pop("Device setup");
  if (!error)
  {
    // Now we start
    float64 grayabsorb = g_params.find<float>("gray_absorption",0.2);
    bool a_eq_e = g_params.find<bool>("a_eq_e",true);
 
    int endP = 0;
    int startP = 0;
    int nPR = 0;

    // Loop over chunks of particles as big as we can fit in dev mem
    while(endP < nP)
    {
      // Set range and draw first chunk
      endP = startP + len;
      if (endP > nP) endP = nP;
      #ifdef SPLOTCH_PARAVIEW
      // Splotch paraview uses an extra void pointer.
      // This is only to allow compilation with/without paraview, splotch paraview wouldnt actually be in this code block
      void* dummy;
      nPR += cu_draw_chunk(mydevID, (cu_particle_sim *) &(particle[startP]), endP-startP, Pic_host, &gv, a_eq_e, grayabsorb, xres, yres, doLogs, dummy);
      #else
      nPR += cu_draw_chunk(mydevID, (cu_particle_sim *) &(particle[startP]), endP-startP, Pic_host, &gv, a_eq_e, grayabsorb, xres, yres, doLogs);
      #endif
      /*
      REMOVED TILED IMPLEMENTATION
#ifndef CUDA_FULL_ATOMICS
      // Combine host render of large particles to final image
      // No need to do this for atomic implementation
      tstack_push("combine images");
      for (int x=0; x<xres; x++)
        for (int y=0; y<yres; y++)
          pic[x][y] += Pic_host[x][y];
      tstack_pop("combine images");
#endif
*/
    // Splotch paraview doesnt use mpimgr
    #ifndef SPLOTCH_PARAVIEW
      cout << "Rank " << mpiMgr.rank() << ": Rendered " << nPR << "/" << nP << " particles" << endl << endl;
    #endif
      startP = endP;
    }
    // Get device image and combine with final image
    tstack_push("add_device_image()");
    add_device_image(pic, &gv, xres, yres);
    tstack_pop("add_device_image()");
    tstack_pop("CUDA");
    cu_end(&gv);
  }

 }

 
void setup_colormap(int ptypes, vector<COLOURMAP> &amap, cu_gpu_vars* gv)
{
//init C style colormap
  cu_color_map_entry *amapD;//amap for Device
  int *amapDTypeStartPos; //begin indexes of ptypes in the linear amapD[]
  amapDTypeStartPos =new int[ptypes];
  int curPtypeStartPos =0;
  int size =0;
  //first we need to count all the entries to get colormap size
  for (int i=0; i<amap.size(); i++)
    size += amap[i].size();

  //then fill up the colormap amapD
  amapD =new cu_color_map_entry[size];
  int j,index =0;
  for(int i=0; i<amap.size(); i++)
    {
    for (j=0; j<amap[i].size(); j++)
      {
      amapD[index].val = amap[i].getX(j);
      COLOUR c (amap[i].getY(j));
      amapD[index].color.r = c.r;
      amapD[index].color.g = c.g;
      amapD[index].color.b = c.b;
      index++;
      }
    amapDTypeStartPos[i] = curPtypeStartPos;
    curPtypeStartPos += j;
    }
  //now let cuda init colormap on device
  cu_colormap_info clmp_info;
  clmp_info.map =amapD;
  clmp_info.mapSize =size;
  clmp_info.ptype_points =amapDTypeStartPos;
  clmp_info.ptypes =ptypes;
  cu_init_colormap(clmp_info, gv);

  delete []amapD;
  delete []amapDTypeStartPos;
}
