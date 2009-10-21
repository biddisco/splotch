#ifndef SPLOTCH_KERNEL_H
#define SPLOTCH_KERNEL_H
/*
Try accelating splotch with CUDA. July 2009. 
*/

#include "splotch_cuda.h"

//Function definitions
__device__ void     d_dump_pr(cu_param_range *pr);
__device__ float    my_asinh(float val);
__device__ void     my_normalize(float minv, float maxv, float &val);

__device__ void     
dump_transform(cu_param_transform *ptrans);

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

    //C1, mincol, maxcol
    if (pr->log_col[p[m].type])
	    p[m].C1 = log(p[m].C1);
    if(pr->asinh_col[p[m].type])
	    p[m].C1 = my_asinh(p[m].C1);
      
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
	}


/*  //after all threads calculating p[m]
    __syncthreads()
    // now find the min-maxes in this block
    if (m=...) //only one thread should do tihs
    {
        for ...
        {
            
        }
    }
*/
}

//Range by kernel step 2
__global__ void k_range2
(cu_param_range *pr, cu_particle_sim *p, int n, int itype, 
 float minval_int, float maxval_int, 
 float minval_col, float maxval_col)
{
    //first get the index m of this thread
    int m;
    m =blockIdx.x *blockDim.x + threadIdx.x;
    if (m >n)
        m =n;

    //do the calculation
    if(p[m].type == itype)///clamp into (min,max)
    {
       my_normalize(minval_int,maxval_int,p[m].I);
       my_normalize(minval_col,maxval_col,p[m].C1);
       if (pr->col_vector[p[m].type])
       {
            my_normalize(minval_col,maxval_col,p[m].C2);
            my_normalize(minval_col,maxval_col,p[m].C3);
       }
    }        
}

//Transform by kernel
__global__ void k_transform
(cu_particle_sim *p, int n, cu_param_transform *ptrans)
{
    //first get the index m of this thread
    int m;
    m =blockIdx.x *blockDim.x + threadIdx.x;
    if (m >n)
        m =n;
#ifdef _DEVICEEMU //for debug only
    if (m==0) 
    {
        dump_transform(ptrans);
        printf("\np0 = (%f,  %f,  %f)", p[0].x, p[0].y, p[0].z);
        printf("\npn = (%f,  %f,  %f)", p[n].x, p[n].y, p[n].z);
    }
#endif

    //copy parameters to __share__ local memory? later

    //now do x,y,z
    float x,y,z;
    x =p[m].x*ptrans->p[0] + p[m].y*ptrans->p[1] + p[m].z*ptrans->p[2] + ptrans->p[3];
    y =p[m].x*ptrans->p[4] + p[m].y*ptrans->p[5] + p[m].z*ptrans->p[6] + ptrans->p[7];
    z =p[m].x*ptrans->p[8] + p[m].y*ptrans->p[9] + p[m].z*ptrans->p[10]+ ptrans->p[11];
    p[m].x =x;
    p[m].y =y;
    p[m].z =z;
#ifdef _DEVICEEMU //for debug only
    if (m==0) 
    {
        printf("\nresult: = (%f,  %f,  %f)", p[0].x, p[0].y, p[0].z);
    }
#endif

    //do ro and r
    float   xfac =ptrans->xfac;
    if(!ptrans->projection)
    {
      p[m].x = ptrans->res*.5 * (p[m].x+ptrans->fovfct*ptrans->dist)*xfac;
      p[m].y = ptrans->res*.5 * (p[m].y+ptrans->fovfct*ptrans->dist)*xfac;
    }
    else
    {
      xfac=1./(ptrans->fovfct*p[m].z);
      p[m].x = ptrans->res*.5 * (p[m].x+ptrans->fovfct*p[m].z)*xfac;
      p[m].y = ptrans->res*.5 * (p[m].y+ptrans->fovfct*p[m].z)*xfac;
    }

    p[m].ro = p[m].r;
    p[m].r = p[m].r *ptrans->res*.5*xfac;

    if(ptrans->minhsmlpixel)
        if ((p[m].r <= 0.5) && (p[m].r >= 0.0))
        {
             p[m].r = 0.5;
             p[m].ro = p[m].r/(ptrans->res*.5*xfac);
        }
#ifdef _DEVICEEMU //for debug only
    if (m==0) 
    {
        printf("\n\nk_transform(), p[0], xfac: = %f\n", xfac);
    }
#endif

}

/////////help functions///////////////////////////////////
__device__ float    my_asinh(float val)
{
      return log(val+sqrt(1.+val*val));
}

__device__ void my_normalize(float minv, float maxv, float &val)
{
  if (minv!=maxv) val =  (val-minv)/(maxv-minv);
}


#ifdef _DEVICEEMU
//dump the parameters for ranging
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

//dump the parameters for transformation
__device__ void     
dump_transform(cu_param_transform *ptrans)
{
    printf("\ndump transfrom parameter from kernel\n");
    printf("matrix\n");
    printf("%f,  %f,  %f, %f\n", ptrans->p[0],ptrans->p[1],ptrans->p[2],ptrans->p[3]);
    printf("%f,  %f,  %f, %f\n", ptrans->p[4],ptrans->p[5],ptrans->p[6],ptrans->p[7]);
    printf("%f,  %f,  %f, %f\n", ptrans->p[8],ptrans->p[9],ptrans->p[10],ptrans->p[11]);
    printf("\nprojection=%d,  res=%d, fovfct=%f, dist=%f, xfac=%f, minhsmlpixel=%d\n",
        ptrans->projection, ptrans->res, ptrans->fovfct, 
        ptrans->dist, ptrans->xfac, ptrans->minhsmlpixel);
}
                    
#endif

#ifdef NOT_USED_ANYMORE
__device__ void     get_minmax(float &minv, float &maxv, float val);
__device__ void get_minmax(float &minv, float &maxv, float val)
{
  minv=min(minv,val);
  maxv=max(maxv,val);
}
#endif



#endif // #ifndef SPLOTCH_KERNEL_H