#ifndef SPLOTCH_KERNEL_H
#define SPLOTCH_KERNEL_H
/*
Try accelating splotch with CUDA. July 2009. 
*/

#include "splotch_cuda.h"

__global__ void k_range(d_particle_sim *pd, int n)
{
#ifdef _DEVICEEMU
//    printf("\nk_rnage()");
#endif

    int idx;
    idx =blockIdx.x *blockDim.x + threadIdx.x;
    if (idx >n)
        idx =n;

    pd[idx].x = idx;
    
}

#endif // #ifndef SPLOTCH_KERNEL_H
