#ifndef __KERNEL__
#define __KERNEL__

#include "cuda/splotch_cuda.h"

//#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
//#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
//#endif


//MACROs
#define Pi 3.14159265358979323846264338327950288
#define get_xy_from_sn(sn, xmin, ymin, ymax, x, y)\
        {int x1 =sn/(ymax-ymin); int y1 =sn-x1*(ymax-ymin);\
         x  =x1 +xmin; y  =y1 +ymin;}
#define get_sn_from_xy(x,y,maxy,miny, sn)\
    {sn =x*(maxy-miny) +y;}

#define get_minmax(minv, maxv, val) \
         minv=min(minv,val); \
         maxv=max(maxv,val);
#define MAXSIZE 1000

/////////constant memory declaration /////////////////////

__constant__ cu_color_map_entry dmap[MAXSIZE];
__constant__ int ptype_points[10];
__constant__ cu_param dparams;

//help functions

__device__ __forceinline__ void clamp (float minv, float maxv, float &val)
  {
  val = min(maxv, max(minv, val));
  }

//Transform by kernel
__global__ void k_transform
  (cu_particle_sim *p, unsigned long *p_region, bool *p_active, int n, int MaxBlock)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n) return;

  //now do x,y,z
  float x,y,z;
  x =p[m].x*dparams.p[0] + p[m].y*dparams.p[1] + p[m].z*dparams.p[2] + dparams.p[3];
  y =p[m].x*dparams.p[4] + p[m].y*dparams.p[5] + p[m].z*dparams.p[6] + dparams.p[7];
  z =p[m].x*dparams.p[8] + p[m].y*dparams.p[9] + p[m].z*dparams.p[10]+ dparams.p[11];

  //do r
  float xfac = dparams.xfac;
  const float   res2 = 0.5f*dparams.xres;
  const float   ycorr = 0.5f*(dparams.yres-dparams.xres);
  if (!dparams.projection)
    {
    x = res2 * (x+dparams.fovfct*dparams.dist)*xfac;
    y = res2 * (y+dparams.fovfct*dparams.dist)*xfac + ycorr;
    }
  else
    {
    xfac=1.0f/(dparams.fovfct*z);
    x = res2 * (x+dparams.fovfct*z)*xfac;
    y = res2 * (y+dparams.fovfct*z)*xfac + ycorr;
    }

  float r = p[m].r;
  p[m].I /= r;
  r *= res2*xfac;

  const float rfac= sqrt(r*r + 0.25f*dparams.minrad_pix*dparams.minrad_pix)/r;
  r *= rfac;
  p[m].I /= rfac;

  p[m].x = x;
  p[m].y = y;
  p[m].r = r;

  p[m].active = false;
  p_active[m] = false;

  // compute region occupied by the partile
  const float rfacr=dparams.rfac*r;
  int minx=int(x-rfacr+1.f);
  if (minx>=dparams.xres) return;
  minx=max(minx,0);

  int maxx=int(x+rfacr+1.f);
  if (maxx<=0) return;
  maxx=min(maxx,dparams.xres);
  if (minx>=maxx) return;

  int miny=int(y-rfacr+1.f);
  if (miny>=dparams.yres) return;
  miny=max(miny,0);

  int maxy=int(y+rfacr+1.f);
  if (maxy<=0) return;
  maxy=min(maxy,dparams.yres);
  if (miny>=maxy) return;
 
  p[m].active = true;
  p_region[m] = (unsigned long) (maxx-minx)*(maxy-miny);
  if (p_region[m] > (unsigned long) MaxBlock)   p_active[m] = false;
  else  p_active[m] = true;
  }

//fetch a color from color table on device
__device__ __forceinline__ cu_color get_color(int ptype, float val, int mapSize, int ptypes)
  {
  __shared__ int map_size;
  __shared__ int map_ptypes;

  map_size = mapSize;
  map_ptypes = ptypes;
  //first find the right entry for this ptype
  int     start, end;
  start =ptype_points[ptype];
  if ( ptype == map_ptypes-1)//the last type
    end = map_size-1;
  else
    end = ptype_points[ptype+1]-1;

  //search the section of this type to find the val
  int i=start;
  while ((val>dmap[i+1].val) && (i<end)) ++i;

  const float fract = (val-dmap[i].val)/(dmap[i+1].val-dmap[i].val);
  cu_color clr1=dmap[i].color, clr2=dmap[i+1].color;
  cu_color        clr;
  clr.r =clr1.r + fract*(clr2.r-clr1.r);
  clr.g =clr1.g + fract*(clr2.g-clr1.g);
  clr.b =clr1.b + fract*(clr2.b-clr1.b);

  return clr;
  }

//colorize by kernel
__global__ void k_colorize
  (cu_particle_sim *p2, int mapSize, int types, int n)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >= n) return; 
  
  int ptype = p2[m].type;
  float col1=p2[m].e.r,col2=p2[m].e.g,col3=p2[m].e.b;
  clamp (0.0000001,0.9999999,col1);
  if (dparams.col_vector[ptype])
    {
    clamp (0.0000001,0.9999999,col2);
    clamp (0.0000001,0.9999999,col3);
    }
  float intensity = p2[m].I;
  clamp (0.0000001,0.9999999,intensity);
  intensity *= dparams.brightness[ptype];

  cu_color e;
  if (dparams.col_vector[ptype])   // color from file
    {
    e.r=col1*intensity;
    e.g=col2*intensity;
    e.b=col3*intensity;
    }
  else   // get color, associated from physical quantity contained in e.r, from lookup table
    {
  //first find the right entry for this ptype
      if (ptype<types)
      {
        e = get_color(ptype, col1, mapSize, types);
        e.r *= intensity;
        e.g *= intensity;
        e.b *= intensity;
      }
      else
      { e.r = e.g = e.b = 0.0f; }
    }
  p2[m].e.r = e.r;
  p2[m].e.g = e.g;
  p2[m].e.b = e.b;
 }

//device render function k_render1
__global__ void k_render1
  (unsigned long *pos, int nP, unsigned long FragRendered, cu_particle_sim *part,
   void *buf, int *index, bool a_eq_e, float grayabsorb)
{
  //first get the index m of this thread
  int m = blockIdx.x;    // particle id
//  m =blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=nP)//m goes from 0 to nP-1
    return;

  int npix = threadIdx.x;   // pixel number
  
  __shared__ int minx, maxx, miny, maxy, reg;
  __shared__ float r, posx, posy;
  __shared__ cu_color e,q;
  __shared__ unsigned long ppos;
  __shared__ float sigma0;
  if (threadIdx.x == 0)
  {
   e.r = part[m].e.r; e.g = part[m].e.g; e.b = part[m].e.b;
   r = part[m].r; posx = part[m].x; posy = part[m].y;
   float rfacr = dparams.rfac*r;

   minx=int(posx-rfacr+1.f);
   minx=max(minx,0);
   maxx=int(posx+rfacr+1.f);
   maxx=min(maxx,dparams.xres);
   miny=int(posy-rfacr+1.f);
   miny=max(miny,0);
   maxy=int(posy+rfacr+1.f);
   maxy=min(maxy,dparams.yres);
   reg = (maxx-minx)*(maxy-miny);
 
   if (!a_eq_e)
   {
     q.r = e.r/(e.r+grayabsorb);
     q.g = e.g/(e.g+grayabsorb);
     q.b = e.b/(e.b+grayabsorb);
   }

   const float powtmp = __powf(Pi,1.0f/3.0f);
   sigma0 = powtmp/sqrt(2.f*Pi);
   const float intens = -0.5f/(2.f*sqrt(Pi)*powtmp);
   e.r*=intens; e.g*=intens; e.b*=intens;

   ppos = pos[m];     //particle m absolute end position
 //  printf("part = %d, size = %d, ppos =%u \n",m,reg,ppos); 
  }
  __syncthreads();
  if (npix > reg) return;

 //now do the rendering
  float radsq = 2.25f*r*r;
  float stp = -0.5f/(r*r*sigma0*sigma0);

  // relative starting position:
  unsigned long fpos = ppos - FragRendered - (unsigned long) (reg - npix); 
  //if (threadIdx.x == 0) printf("fpos = %u\n",fpos); 
  int x = npix/(maxy-miny) + minx;
  int y = npix%(maxy-miny) + miny;

  //make fbuf the right type
  cu_fragment_AeqE        *fbuf;
  cu_fragment_AneqE       *fbuf1;
  if (a_eq_e)
    fbuf =(cu_fragment_AeqE*) buf;
  else
    fbuf1 =(cu_fragment_AneqE*)buf;

  if (a_eq_e)
  {
 //   for (int x=minx; x<maxx; ++x)
 //   {
     float dxsq=(x-posx)*(x-posx);
 //    for (int y=miny; y<maxy; ++y)
 //     {
        float dsq = (y-posy)*(y-posy) + dxsq;
        if (dsq<radsq)
        {
          float att = __expf(stp*dsq);
          fbuf[fpos].aR = att*e.r;
          fbuf[fpos].aG = att*e.g;
          fbuf[fpos].aB = att*e.b;
        }
        else
        {
          fbuf[fpos].aR = 0.0f;
          fbuf[fpos].aG = 0.0f;
          fbuf[fpos].aB = 0.0f;
        }
        index[fpos] = y+dparams.yres*x;  // pixel index in the image
  //      fpos++;
  //    }//y
  //  }//x
  }
  else
  {
  //  for (int x=minx; x<maxx; ++x)
  //  {
     float dxsq=(x-posx)*(x-posx);
  //   for (int y=miny; y<maxy; ++y)
  //    {
        float dsq = (y-posy)*(y-posy) + dxsq;
        if (dsq<radsq)
        {
          float att = __expf(stp*dsq);
          float expm1;
          expm1 =__expf(att*e.r)-1.0f;
          fbuf1[fpos].aR = expm1;
          fbuf1[fpos].qR = q.r;
          expm1 =__expf(att*e.g)-1.0f;
          fbuf1[fpos].aG = expm1;
          fbuf1[fpos].qG = q.g;
          expm1 =__expf(att*e.b)-1.0f;
          fbuf1[fpos].aB = expm1;
          fbuf1[fpos].qB = q.b;
        }
        else
        {
          fbuf1[fpos].aR = 0.0f;
          fbuf1[fpos].aG = 0.0f;
          fbuf1[fpos].aB = 0.0f;
          fbuf1[fpos].qR = 1.0f;
          fbuf1[fpos].qG = 1.0f;
          fbuf1[fpos].qB = 1.0f;
        }
        index[fpos] = x*dparams.yres+y;
   //     fpos++;
   //   }//y
   // }//x
  }
}


__global__ void k_clear(int n, cu_color *pic)
{
   //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n) m =n;

   pic[m].r = 0.0f;
   pic[m].g = 0.0f;
   pic[m].b = 0.0f;
}

__global__ void k_combine(int n, bool a_eq_e, cu_color *pic, int *index, void *fragBuf)
{
   //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n) return;
  
  if (a_eq_e)
  {
    cu_fragment_AeqE *fragBufAeqE = (cu_fragment_AeqE *)fragBuf;
    pic[index[m]].r += fragBufAeqE[m].aR;
    pic[index[m]].g += fragBufAeqE[m].aG;
    pic[index[m]].b += fragBufAeqE[m].aB;
  }
  else
  {
    cu_fragment_AneqE *fragBufAneqE = (cu_fragment_AneqE *)fragBuf;
    pic[index[m]].r += fragBufAneqE[m].aR *
				   (pic[index[m]].r - fragBufAneqE[m].qR);
    pic[index[m]].g += fragBufAneqE[m].aG *
				   (pic[index[m]].g - fragBufAneqE[m].qG); 
    pic[index[m]].b += fragBufAneqE[m].aB *
                                   (pic[index[m]].b - fragBufAneqE[m].qB);
  }
}


struct sum_op
{
  __host__ __device__
  cu_fragment_AeqE operator()(cu_fragment_AeqE& p1, cu_fragment_AeqE& p2) const{

    cu_fragment_AeqE sum;
    sum.aR = p1.aR + p2.aR;
    sum.aG = p1.aG + p2.aG;
    sum.aB = p1.aB + p2.aB;
    return sum; 
   } 
};

struct reg_notValid
  {
    __host__ __device__
    bool operator()(const bool flag)
    {
      return !flag;
    }
  };

#endif

