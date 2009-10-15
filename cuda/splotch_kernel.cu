#ifndef SPLOTCH_KERNEL_H
#define SPLOTCH_KERNEL_H
/*
Try accelating splotch with CUDA. July 2009. 
*/

#include "splotch_cuda.h"

//Function definitions
__device__ void     d_dump_pr(cu_param_range *pr);
__device__ float    my_asinh(float val);

//Range by kernel step 1
__global__ void k_range1(cu_param_range *pr, cu_particle_sim *p, int n)
{
#ifdef _DEVICEEMU
//    printf("\nk_rnage()");
#endif

    //first get the index m of this thread
    int m;
    m =blockIdx.x *blockDim.x + threadIdx.x;
    if (m >n)
        m =n;

    //now do the calc
    //I, minint, maxint
    if (pr->log_int[p[m].type])
	    p[m].I = log(p[m].I);
//    get_minmax(pr->minint[p[m].type], pr->maxint[p[m].type], p[m].I);

    //C1, mincol, maxcol
    if (pr->log_col[p[m].type])
	    p[m].C1 = log(p[m].C1);
    if(pr->asinh_col[p[m].type])
	    p[m].C1 = my_asinh(p[m].C1);
//    get_minmax(pr->mincol[p[m].type], pr->maxcol[p[m].type], p[m].C1);
      
    //C2, C3, mincol, maxcol
    if (pr->col_vector[p[m].type])
	{
	  if (pr->log_col[p[m].type])
      {
          p[m].C2 = log(p[m].C2);
          p[m].C3 = log(p[m].C3);
      }   
	  if (pr->asinh_col[p[m].type])
      {
	      p[m].C2 = my_asinh(p[m].C2);
	      p[m].C3 = my_asinh(p[m].C3);
	  }
//	  get_minmax(pr->mincol[p[m].type], pr->maxcol[p[m].type], p[m].C2);
//	  get_minmax(pr->mincol[p[m].type], pr->maxcol[p[m].type], p[m].C3);
	}

    //after all threads calculating p[m]
    _syncthreads()
    // now find the min-maxes in this block
    if (m=...) //only one thread should do tihs
    {
        for ...
        {
            
        }
    }
}

//Range by kernel step 2
__global__ void k_range2
(cu_particle_sim *pd, int n, float minval_int, 
float maxval_int,float minval_col,float maxval_col)
{

}

#ifdef _DEVICEEMU
__device__ void d_dump_pr(cu_param_range *pr)
{
    printf("\ndump_pr:\n");
    printf("col_vector, log_int,log_col,asinh_col,");
    printf("mincol, maxcol, minint,maxint\n");

    for (int i=0; i<pr->ptypes; i++)
    {
        printf("%d, %d, %d, %d, %f, %f, %f, %f\n",
            pr->col_vector[i], pr->log_int[i], 
            pr->log_col[i],	pr->asinh_col[i], pr->mincol[i],
    		pr->maxcol[i], pr->minint[i], pr->maxint[i]);
    }                     
}                         
#endif

#ifdef NOT_USED_ANYMORE
__device__ void     get_minmax(float &minv, float &maxv, float val);
__device__ void get_minmax(float &minv, float &maxv, float val)
{
  minv=min(minv,val);
  maxv=max(maxv,val);
}

__device__ float    my_asinh(float val)
{
      return log(val+sqrt(1.+val*val));
}
#endif



#endif // #ifndef SPLOTCH_KERNEL_H