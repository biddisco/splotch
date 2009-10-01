/*
Try accelating splotch with CUDA. July 2009. 
*/
//#define s__DEVICE_EMULATION__ not working!

#ifndef _CPP_SPLOTCH_KERNEL_H_
#define _CPP_SPLOTCH_KERNEL_H_

#include "splotch_cu_data.h"

__global__ void k_preCalc( G_VARS *g_vars, PARTICLE *p)
{
	int m;
    m = blockIdx.x *512 + threadIdx.x;
    
/*	if (p[m].z<=0) return;
	if (p[m].z<=g_vars->zminval) return;
	if (p[m].z>=g_vars->zmaxval) return;
*/
	 float  ro=p[m].ro;
	 float  sigma=ro*g_vars->sigma0;
	 float  fac0=-0.5/(sigma*sigma);
	 float  i0=g_vars->i00/ro;
	 float  prefac1 = ro/p[m].r * ro/p[m].r * fac0;
	 float  prefac2 = -0.5*i0;

	p[m].prefac1 =prefac1;
	p[m].prefac2 =prefac2;
	
	 float  posx=p[m].x, posy=p[m].y;
	 int minx, maxx, miny, maxy;
	 minx = max(int(posx-g_vars->rfac*p[m].r),0); 
	 maxx = min(int(posx+g_vars->rfac*p[m].r+1),g_vars->res);
	 miny = max(int(posy-g_vars->rfac*p[m].r),g_vars->ycut0); 
	 maxy = min(int(posy+g_vars->rfac*p[m].r+1),g_vars->ycut1);

	p[m].minx =minx;
	p[m].maxx =maxx;
	p[m].miny =miny;
	p[m].maxy =maxy;
    
}

__global__ void k_initFBuf(FRAGMENT *fBuf)
{
/*This is for one frame buffer. Issued as <<<(800*2),400>>>
    unsigned long index;
    index =blockIdx.x *800 +blockIdx.y*400 +threadIdx.x;

    //initialize, c is set to 0 automatically
    fBuf[index].k_red =1.0;
    fBuf[index].k_green=1.0;
    fBuf[index].k_blue=1.0;
*/

    //for multiple buffers, like 64, issued as <<<(800*800),64>>>
    unsigned int index;
    index =threadIdx.x*800*800 + blockIdx.y *800 +blockIdx.x;
    fBuf[index].k_red =1.0;
    fBuf[index].k_green =1.0;
    fBuf[index].k_blue =1.0;
//    fBuf[index].c_red =(float)index; only for debug
}

__device__ float d_xexp( float* exptab, float arg)
{
    //pre set
    int nexp =10000;
    float   max =-20.0;

   //retrieve an element from exp table to result, all on device
    float expfac;
    expfac =nexp/max;
    
    arg *= expfac;
    int iarg=int(arg);
    if (iarg<0) 
        return 1.0;
    
    if (iarg>=nexp-1) 
        return 0.;
        
    float frac=arg-iarg;
    return (1.0-frac)*exptab[iarg]+frac*exptab[iarg+1];    
}

__global__ void k_shadeA(int nStart, int nEnd, PARTICLE *p, FRAGMENT *fragBuf,  G_VARS *g_vars, float* exptab)
{
    int index = 64/gridDim.x *blockIdx.x + threadIdx.x; //from 0 to 63

//  fragBuf += sizeof(FRAGMENT)*800*800*index; this is wrong!
//  fragBuf = &(fragBuf[index*800*800]);
    unsigned int    bufMov =index*800*800;

    int l = (nEnd -nStart)/64;
    int start, end;
    start =nStart +index*l;
    end =start+l;
    if (index == 63)
        end =nEnd;    
    
    float	k,c;
//	for (int m=index*l; m<(index+1)*l; ++m)
    for (int m=start; m<=end; ++m)
	{
		if (p[m].z<=0) continue;
		if (p[m].z<= g_vars->zminval) continue;
		if (p[m].z>= g_vars->zmaxval) continue;

		for (int y=p[m].miny; y<p[m].maxy; ++y)
		{
		  float  ysq, radsq;
          ysq =(y-p[m].y)*(y-p[m].y);
		  radsq = g_vars->rfac*g_vars->rfac*p[m].r*p[m].r;

		  for (int x=p[m].minx; x<p[m].maxx; ++x)
		  {
			float  dsq = (x-p[m].x)*(x-p[m].x) + ysq;
			if (dsq<radsq)
			{
			  float  fac = p[m].prefac2*d_xexp(exptab,p[m].prefac1*dsq);
              
              //below all fragBuf[x][y] are replaces by fragBuf[y*800+x]

			  k =d_xexp(exptab, fac*p[m].a[0]);
//              c =p[m].q[0] -p[m].q[0]*d_xexp(exptab,fac*p[m].a[0]);
              c =p[m].q[0] * (1.0-k);
			  fragBuf[bufMov+y*800+x].k_red =fragBuf[bufMov+y*800+x].k_red *k;
			  fragBuf[bufMov+y*800+x].c_red =k *fragBuf[bufMov+y*800+x].c_red +c;

			  k =d_xexp(exptab, fac*p[m].a[1]);
//			  c =p[m].q[1] -p[m].q[1]*d_xexp(exptab,fac*p[m].a[1]);
              c =p[m].q[1] * (1.0-k);
              fragBuf[bufMov+y*800+x].k_green =fragBuf[bufMov+y*800+x].k_green *k;
			  fragBuf[bufMov+y*800+x].c_green =k *fragBuf[bufMov+y*800+x].c_green +c;

			  k =d_xexp(exptab, fac*p[m].a[2]);
//			  c =p[m].q[2] -p[m].q[2]*d_xexp(exptab,fac*p[m].a[2]);
              c =p[m].q[2] * (1.0-k);
              fragBuf[bufMov+y*800+x].k_blue =fragBuf[bufMov+y*800+x].k_blue *k;
			  fragBuf[bufMov+y*800+x].c_blue =k *fragBuf[bufMov+y*800+x].c_blue +c;
			}
			//else out of circle, already initialized as (k=1, c=0)
		  }
		}
	}
}

__global__ void k_combineA(FRAGMENT *fBuf)
{
    unsigned int index_first, index_n, distance;
    distance =800*800;
    //it is issued as <<<(m,n),k>>>
//    index_first = 800*blockIdx.y + blockDim.x *blockIdx.x +threadIdx.x; seems not right
    //(800,800),1
    index_first =800 *blockIdx.y +blockIdx.x;
    index_n =index_first;

//    fBuf[index_first].c_red =(float)index_first; debug only
    
    //to every pixel...
    for (int i=1; i<64; i++)
    {
        index_n +=distance;
        fBuf[index_first].k_red =fBuf[index_first].k_red *fBuf[index_n].k_red ;
        fBuf[index_first].c_red =fBuf[index_n].k_red *fBuf[index_first].c_red 
            +fBuf[index_n].c_red ;
        fBuf[index_first].k_green =fBuf[index_first].k_green *fBuf[index_n].k_green ;
        fBuf[index_first].c_green =fBuf[index_n].k_green *fBuf[index_first].c_green 
            +fBuf[index_n].c_green ;
        fBuf[index_first].k_blue =fBuf[index_first].k_blue *fBuf[index_n].k_blue ;
        fBuf[index_first].c_blue =fBuf[index_n].k_blue *fBuf[index_first].c_blue 
            +fBuf[index_n].c_blue ;
    }
}

/*
__global__ void k_xexp( int nexp, float max, float arg, float *result, float *exptab)
{
    // write data to global memory
//    const unsigned int tid = threadIdx.x;
    
    //retrieve an element from exp table to result, all on device
    float expfac;
    expfac =nexp/max;
    
    arg *= expfac;
    int iarg=int(arg);
    if (iarg<0) 
    {
        *result =1.;
        return ;
    }
    
    if (iarg>=nexp-1) 
    {
        *result =0.;
        return ;
    }
    
    float frac=arg-iarg;
    *result = (1.-frac)*exptab[iarg]+frac*exptab[iarg+1];    
    return;
}

__global__ void k_testDouble(float *d)
{
    *d =9.9;
}

__global__ void test_kernel( float *table, int n, float *result)
{
    // write data to global memory
//    const unsigned int tid = threadIdx.x;
    
    //retrieve an element from exp table to result, all on device
    *result =table[n];
}

__global__ void test_kernel2(int *dl)
{
    int k; 
    k = blockIdx.x *512 + threadIdx.x;
    dl[k] =k;
}
*/
#endif // #ifndef _CPP_SPLOTCH_KERNEL_H_
