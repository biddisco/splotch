
#ifndef CUDA_KERNEL_H
#define CUDA_KERNEL_H

#include <thrust/random.h>
#include <thrust/tuple.h>
#include "cuda/cuda_utils.h"
#include <cstdio>

//MACROs
#define Pi 3.141592653589793238462643383279502884197
#define MAXSIZE 1000

__device__ __forceinline__ void clamp (float minv, float maxv, float &val);
__device__ __forceinline__ double my_asinh (double val);
__device__ __forceinline__ cu_color get_color(int ptype, float val, int map_size, int map_ptypes);
__global__ void k_process(cu_particle_sim *p, int *p_active, int n, int mapSize, int types, int tile_sidex, int tile_sidey, int width, int nxtiles, int nytiles);
__global__ void k_range(int nP, cu_particle_sim *p);
__device__ int  pixelLocalToGlobal(int lpix, int xo, int yo, int width, int tile_sidey);
__global__ void k_renderC2(int nP, cu_particle_sim *part, int *tileId, int *tilepart, cu_color *pic, cu_color *pic1, cu_color *pic2, cu_color *pic3, int tile_sidex, int tile_sidey, int width, int nytiles);
__global__ void k_indexC3(int n, cu_particle_sim *part, int *index);
__global__ void k_renderC3(int nC3, int *index, cu_particle_sim *part, cu_color *pic);

#ifndef CUDA_USE_ATOMICS
  __global__ void k_add_images(int n, cu_color *pic, cu_color *pic1, cu_color *pic2, cu_color *pic3);
#endif

// check for non-active and big particles to remove from the device
struct particle_notValid
  {
    __host__ __device__    
    bool operator()(const int flag)
    {
      return (flag < 0);
    }
  };
  
__host__ __device__
inline unsigned int hash(unsigned int a)
{
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}
  

// check for active big particles to copy back to the host
struct reg_notValid
  {
    __host__ __device__
    reg_notValid(int prf) 
    {
      _prf = prf/100.0f;
      unsigned int seed = hash(prf);
      // seed a random number generator
      rng = thrust::default_random_engine(seed);

      // create a mapping from random numbers to [0,1)
      u = thrust::uniform_real_distribution<float>(0,1);
    }
    
    __host__ __device__
    bool operator()(int &flag, int index)
    {
      bool copyback = (flag==-2);
      rng.discard(index);
      float x = u(rng);
      printf("%7.3f ",x);
      if (x>_prf) {
        copyback = false;
      }
      return copyback;
    }
    __host__ __device__
    bool operator()(thrust::tuple<int &, int> ref)
    {
      bool copyback = (thrust::get<0>(ref)==-2);
      rng.discard(thrust::get<1>(ref));
      float x = u(rng);
      printf("%7.3f ",x);
      if (x>_prf) {
        copyback = false;
      }
      return copyback;
    }
    
    
    float _prf;
    // variables needed to create random numbers
    thrust::default_random_engine rng;
    thrust::uniform_real_distribution<float> u;
  };

struct sum_op
{
  __host__ __device__
  cu_particle_sim operator()(cu_particle_sim& p1, cu_particle_sim& p2) const{

    cu_particle_sim sum;
    sum = p1;
    sum.e.r = p1.e.r + p2.e.r;
    sum.e.g = p1.e.g + p2.e.g;
    sum.e.b = p1.e.b + p2.e.b;

    return sum; 
   } 
};
#endif