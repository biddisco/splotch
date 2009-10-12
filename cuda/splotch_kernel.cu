#ifndef SPLOTCH_KERNEL_H
#define SPLOTCH_KERNEL_H
/*
Try accelating splotch with CUDA. July 2009. 
*/

#include "splotch_cuda.h"

__global__ void k_range()
{
#ifdef _DEVICEEMU
    printf("k_rnage()");
#endif
}

#endif // #ifndef SPLOTCH_KERNEL_H
