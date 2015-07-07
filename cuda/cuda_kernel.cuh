
#ifndef CUDA_KERNEL_H
#define CUDA_KERNEL_H

#include <thrust/random.h>
#include <thrust/tuple.h>
#include "cuda/cuda_utils.h"
#include <cstdio>

//#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
//#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
//#endif

//MACROs
#define Pi 3.141592653589793238462643383279502884197
#define MAXSIZE 1000

__device__ __forceinline__ void clamp (float minv, float maxv, float &val);
__device__ __forceinline__ double my_asinh (double val);
__device__ __forceinline__ cu_color get_color(int ptype, float val, int map_size, int map_ptypes);

__global__ void k_range(int nP, cu_particle_sim *p);
__global__ void k_process(cu_particle_sim *p, int n, int mapSize, int types);
__global__ void k_render(int nP, cu_particle_sim *part, cu_color *pic);

#endif
