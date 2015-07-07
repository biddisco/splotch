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
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
//#include <thrust/system_error.h>

#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"

#include "cuda/cuda_render.h"
#include "cuda/cuda_policy.h"
#include "cuda/cuda_kernel.cuh"

#define DUMP_BIG_PARTICLES

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
int cu_draw_chunk(int mydevID, cu_particle_sim *d_particle_data, int nParticle, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, int xres, int yres, bool doLogs, int prf, void *gpudata)
#else 
int cu_draw_chunk(int mydevID, cu_particle_sim *d_particle_data, int nParticle, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, int xres, int yres, bool doLogs)
#endif
{

 cudaError_t error;

  // Copy data particle to device memory
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
  
  // Getting logs etc if necessary
  tstack_push("do logs");
  if(doLogs)
  {
    cu_range(nParticle, gv);
    cudaDeviceSynchronize();
  }
  tstack_pop("do logs");
 
  //--------------------------------------
  //  particle projection and coloring
  //--------------------------------------

  tstack_push("Particle projection & coloring");
  // Project and color particles
  cu_process(nParticle, gv);
  cudaDeviceSynchronize();
  //cout << cudaGetErrorString(cudaGetLastError()) << endl;
  tstack_pop("Particle projection & coloring");

  // Clip particles here?

  //--------------------------------------
  //  particle rendering
  //--------------------------------------

  tstack_push("CUDA Rendering");
  
  cu_render(nParticle, gv);
  cudaDeviceSynchronize();
  tstack_pop("CUDA Rendering");

  return nParticle;

}
/*
#endif
*/
int add_device_image(arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, int xres, int yres)
{
  int res = xres*yres;

  // copy back the image
  tstack_push("Data copy");
  cudaError_t error = cudaMemcpy(&Pic_host[0][0], gv->d_pic, res * sizeof(cu_color), cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) 
  {
    cout << "Rank " << MPI_Manager::GetInstance()->rank() << " Image copy: Device Memcpy error!" << endl;
    return 0;
  }
  tstack_pop("Data copy");

   //error = cudaMemset(gv->d_pic, 0, res * sizeof(cu_color));
  if (error != cudaSuccess) 
  {
    cout << "Rank " << MPI_Manager::GetInstance()->rank() << " memset error to clear device image" << endl;
    return 0;
  }
  return 1;
}
