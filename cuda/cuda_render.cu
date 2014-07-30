/*
 * Copyright (c) 2011-2014
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
 
// Implementation only for the case A=E

#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/remove.h>
#include <thrust/copy.h>
#include <thrust/extrema.h>
#include <thrust/scan.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/system_error.h>

#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"

#include "cuda/cuda_render.h"
#include "cuda/cuda_policy.h"
#include "cuda/cuda_kernel.cuh"

using namespace std;

// FOR NON SEPERABLE COMPILATION
// // check for non-active and big particles to remove from the device
// struct particle_notValid
//   {
//     __host__ __device__ 
//     bool operator()(const int flag)
//     {
//       return (flag < 0);
//     }
//   };

// // check for active big particles to copy back to the host
// struct reg_notValid
//   {
//     __host__ __device__
//     bool operator()(const int flag)
//     {
//       return (flag==-2);
//     }
//   };

// struct sum_op
// {
//   __host__ __device__
//   cu_particle_sim operator()(cu_particle_sim& p1, cu_particle_sim& p2) const{

//     cu_particle_sim sum;
//     sum = p1;
//     sum.e.r = p1.e.r + p2.e.r;
//     sum.e.g = p1.e.g + p2.e.g;
//     sum.e.b = p1.e.b + p2.e.b;

//     return sum; 
//    } 
// };


#ifdef SPLOTCH_PARAVIEW
int cu_draw_chunk(int mydevID, cu_particle_sim *d_particle_data, int nParticle, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, int xres, int yres, bool doLogs, void *gpudata)
#else 
int cu_draw_chunk(int mydevID, cu_particle_sim *d_particle_data, int nParticle, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, int xres, int yres, bool doLogs)
#endif
{
  cudaError_t error;

  // Copy data particle to device memory
  // Paraview version might not need to copy
  tstack_push("Data copy");

#ifdef SPLOTCH_PARAVIEW
  if (gpudata) {
    cu_copy_particles_from_gpubuffer(gpudata, nParticle, gv);
  }
  else {
    cu_copy_particles_to_device(d_particle_data, nParticle, gv);
  }
#else 
  cu_copy_particles_to_device(d_particle_data, nParticle, gv);
#endif

  tstack_pop("Data copy");

  // Get parameters for rendering
  int tile_sidex, tile_sidey, width, nxtiles, nytiles;
  gv->policy->GetTileInfo(&tile_sidex, &tile_sidey, &width, &nxtiles, &nytiles);

  tstack_push("do logs");
  if(doLogs)
  {
    cu_range(nParticle, gv);
    cudaThreadSynchronize();
  }
  tstack_pop("do logs");
 
  //--------------------------------------
  //  particle projection and coloring
  //--------------------------------------

  tstack_push("Particle projection & coloring");
  // Project and color particles, set activity 
  // active: either the tile to which it belongs, -2 for big particle, -1 for clipped or ntx*nty for pointlike particles
  cu_process(nParticle, gv, tile_sidex, tile_sidey, width, nxtiles, nytiles);
  cudaThreadSynchronize();
  //cout << cudaGetErrorString(cudaGetLastError()) << endl;
  tstack_pop("Particle projection & coloring");

  int new_ntiles, newParticle, nHostPart=0;
  int nC3 = 0;
  particle_sim *host_part = 0;
  // Thrust dev pointer to particle data
  thrust::device_ptr<cu_particle_sim> dev_ptr_pd((cu_particle_sim *) gv->d_pd);
  try
  { 
   tstack_push("Particle Filtering");
   thrust::device_ptr<int> dev_ptr_flag((int *) gv->d_active);
   // Select big particles to be processed by the host
   // First create buffer to store large particles
   thrust::device_vector<cu_particle_sim> d_host_part(nParticle);
   // Copy particles from device array to big particle array if active flag == -2 
   thrust::device_vector<cu_particle_sim>::iterator end = thrust:: copy_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, d_host_part.begin(), reg_notValid()); 
   // Check how many big particles we had
   nHostPart = end - d_host_part.begin();

   // Copy back big particles
   if (nHostPart > 0)
   {

    cu_particle_sim *d_host_part_ptr = thrust::raw_pointer_cast(&d_host_part[0]);
    error = cudaHostAlloc((void**) &host_part, nHostPart*sizeof(cu_particle_sim), cudaHostAllocDefault);
    if (error != cudaSuccess) cout << "cudaHostAlloc error!" << endl;
    else
    {
      error = cudaMemcpyAsync(host_part, d_host_part_ptr, nHostPart*sizeof(cu_particle_sim), cudaMemcpyDeviceToHost, 0);
      if (error != cudaSuccess) cout << "Big particles Memcpy error!" << endl;
    }
   }

   //Remove non-active and host particles
   thrust::device_ptr<cu_particle_sim> new_end = thrust::remove_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, particle_notValid());
   newParticle = new_end.get() - dev_ptr_pd.get();
   bool pLeft = true;
   if( newParticle != nParticle )
   {
    cout << endl << "Eliminating inactive particles..." << endl;
    long nLeft = newParticle+nHostPart;
    if(!nLeft)
    {
        cout << "No active particles left!" << endl;
        pLeft = false;
    }
    else
      cout << newParticle+nHostPart << "/" << nParticle << " particles left" << endl; 
    // Remove if active < 0 (i.e. -1 for clipped or -2 for big)
    thrust::remove_if(dev_ptr_flag, dev_ptr_flag+nParticle, particle_notValid());
   }
  
   tstack_pop("Particle Filtering");
   
   // We didnt have to render anything, all particles were inactive.
   // FIXME: We should still be rendering large particles on the host...
   if(!pLeft)
   {  
      return 0;
   }  

   tstack_push("Particle Distribution");

   //sort particles according to their tile id
   thrust::sort_by_key(dev_ptr_flag, dev_ptr_flag + newParticle, dev_ptr_pd);

   //compute number of particles for each tile and their starting position
   cudaMemset(gv->d_tiles,0,nxtiles*nytiles+1);
   thrust::device_ptr<int> dev_ptr_nT((int *) gv->d_tiles);
      thrust::device_ptr<int> dev_ptr_tileID((int *) gv->d_tileID);
   thrust::pair< thrust::device_ptr<int>,thrust::device_ptr<int> > end_tiles = thrust::reduce_by_key(dev_ptr_flag, dev_ptr_flag + newParticle, thrust::make_constant_iterator(1), dev_ptr_tileID, dev_ptr_nT);
   new_ntiles = end_tiles.second.get() - dev_ptr_nT.get();

   if (dev_ptr_tileID[new_ntiles-1] == nxtiles*nytiles) nC3 = dev_ptr_nT[new_ntiles-1];
//   cout << nC3 << " of them are point-like particles" << endl;

   thrust::inclusive_scan(dev_ptr_nT, dev_ptr_nT + new_ntiles, dev_ptr_nT);
   tstack_pop("Particle Distribution");

  }
  catch(thrust::system_error &e)
  {
    // output an error message and exit
    std::cerr << "Error accessing vector element: " << e.what() << std::endl;
    exit(-1);
  }
  catch(std::bad_alloc &e)
  {
    std::cerr << "Couldn't allocate vector" << std::endl;
    exit(-1);
  }

  // ----------------------------
  //   particle proper rendering 
  // ----------------------------

  tstack_push("CUDA Rendering");

  // SMALL (point-like) particles rendering on the device
  if (nC3)
  {
    tstack_push("point-like particles rendering");
    cu_indexC3(newParticle, nC3, gv); 
    //cout << cudaGetErrorString(cudaGetLastError()) << endl;
    thrust::device_ptr<int> dev_ptr_Index(gv->d_index);
    thrust::equal_to<int> binary_pred;
    thrust::pair< thrust::device_ptr<int>,thrust::device_ptr<cu_particle_sim> >  new_end_C3;
    thrust::sort_by_key(dev_ptr_Index, dev_ptr_Index + nC3, dev_ptr_pd + newParticle - nC3);
    new_end_C3 = thrust::reduce_by_key(dev_ptr_Index, dev_ptr_Index + nC3, dev_ptr_pd + newParticle - nC3, dev_ptr_Index, dev_ptr_pd + newParticle - nC3,  binary_pred, sum_op());
    int nAffectedPixels = new_end_C3.first.get() - dev_ptr_Index.get();
    cu_renderC3(newParticle, nC3, nAffectedPixels, gv);
    new_ntiles--; 
    tstack_pop("point-like particles rendering");
  }

  // MEDIUM particles rendering on the device
  // 1 block ----> loop on chunk of particles, 1 thread ----> 1 pixel of the particle
  // number of threads in each block = max number of pixels to be rendered for a particle
  int block_size = 4*width*width;
  int dimGrid = new_ntiles;    // number of blocks = number of tiles
  cout << "number of tiles = " << new_ntiles << endl;

//  cudaEvent_t start, stop;
  if (new_ntiles > 0)
  {
//    cudaEventCreate(&start);
//    cudaEventCreate(&stop);
//    cudaEventRecord(start,0);
    cu_renderC2(newParticle, dimGrid, block_size, a_eq_e, (float) grayabsorb, gv, tile_sidex, tile_sidey, width, nytiles);
    //cout << cudaGetErrorString(cudaGetLastError()) << endl;
//    cudaEventRecord(stop,0);
    cout << "Rank " << MPI_Manager::GetInstance()->rank() << " : Device rendering on " << newParticle << " particles" << endl;
  }

  // LARGE particles rendering on the host 
  if (nHostPart > 0)
  {
     cout << "Rank " << MPI_Manager::GetInstance()->rank() << " : Host rendering on " << nHostPart << " particles" << endl;
     host_funct::render_new(host_part, nHostPart, Pic_host, a_eq_e, grayabsorb);
  }

  if (new_ntiles > 0)
  {
//    cudaEventSynchronize(stop);
    //cout << cudaGetErrorString(cudaGetLastError()) << endl;
//    float elapsedTime;
//    cudaEventElapsedTime(&elapsedTime, start, stop);
//    cout << "Device Rendering Time = " << elapsedTime/1000.0 << endl;
//    cudaEventDestroy(start);
//    cudaEventDestroy(stop);

    cudaThreadSynchronize();
  //  cout << "Rank " << MPI_Manager::GetInstance()->rank() << cudaGetErrorString(cudaGetLastError()) << endl;
  }

  tstack_push("point-like particles rendering");

  tstack_pop("point-like particles rendering");
  //cout << cudaGetErrorString(cudaGetLastError()) << endl;

  tstack_pop("CUDA Rendering");

  if (host_part) cudaFreeHost(host_part);
  return nHostPart+newParticle;
}

int add_device_image(arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, int xres, int yres)
{
  int res = xres*yres;

  #ifndef CUDA_USE_ATOMICS
    // add images on the device: pic+pic1+pic2+pic3
    cu_add_images(res, gv);
  #endif
  // cout << "Rank " << MPI_Manager::GetInstance()->rank() << cudaGetErrorString(cudaGetLastError()) << endl;

  COLOUR *Pic = new COLOUR [res];
  // copy back the image
  tstack_push("Data copy");
  cudaError_t error = cudaMemcpy(Pic, gv->d_pic, res * sizeof(cu_color), cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) 
  {
    cout << "Rank " << MPI_Manager::GetInstance()->rank() << " Image copy: Device Memcpy error!" << endl;
    return 0;
  }
  tstack_pop("Data copy");
  tstack_push("combine images");
  for (int x=0; x<xres; x++)
   for (int y=0; y<yres; y++)
      Pic_host[x][y] += Pic[x*yres+y];
  tstack_pop("combine images");
  delete[] Pic;
  return 1;
}