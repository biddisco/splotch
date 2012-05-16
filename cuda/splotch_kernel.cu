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
  (cu_particle_sim *p, unsigned long *p_region, char *p_active, int n, int MaxSize)
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
  float xfac2 = dparams.xfac;
  const float   res2 = 0.5f*dparams.xres;
  const float   ycorr = 0.5f*(dparams.yres-dparams.xres);
  if (!dparams.projection)
    {
    x = res2 * (x+dparams.fovfct*dparams.dist)*xfac2;
    y = res2 * (y+dparams.fovfct*dparams.dist)*xfac2 + ycorr;
    }
  else
    {
    xfac2=1.0f/(dparams.fovfct*z);
    x = res2 * (x+dparams.fovfct*z)*xfac2;
    y = res2 * (y+dparams.fovfct*z)*xfac2 + ycorr;
    }

  float r = p[m].r;
  float I = p[m].I;
#ifdef SPLOTCH_CLASSIC
  I *= 0.5f*dparams.bfak/r;
  r*= sqrt(2.f)*dparams.sigma0/dparams.h2sigma;  
#else
  I *= 8.f/(Pi*r*r*r);  //SPH kernel normalization
  I *= dparams.h2sigma*sqrt(Pi)*r;  //integral through the center
#endif

  r *= res2*xfac2;
  const float rcorr= sqrt(r*r + dparams.minrad_pix*dparams.minrad_pix)/r;
  r *= rcorr;
#ifdef SPLOTCH_CLASSIC
  I /= rcorr;
#else
  I /= rcorr*rcorr;
#endif

  p[m].x = x;
  p[m].y = y;
  p[m].r = r;
  p[m].I = I;

  p[m].active = false;
  p_active[m] = -1;	// non active particle

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
  p_active[m] = 1;	// active particle
  p_region[m] = (maxx-minx)*(maxy-miny);
  if (p_region[m] > MaxSize) p_active[m] = 2; // particle to be removed and copied back to the host

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
  cu_color e;  
    e.r=p2[m].e.r;
    e.g=p2[m].e.g;
    e.b=p2[m].e.b;

  float intensity = p2[m].I;
  intensity *= dparams.brightness[ptype];

// get color, associated from physical quantity contained in e.r, from lookup table
  if (!dparams.col_vector[ptype])
     e = get_color(ptype, e.r, mapSize, types);

  p2[m].e.r = e.r*intensity;;
  p2[m].e.g = e.g*intensity;;
  p2[m].e.b = e.b*intensity;;
 }

//device render function k_render1
__global__ void k_render1
  (unsigned long *pos, int nP, unsigned long FragRendered, cu_particle_sim *part,
   void *buf, int *index, bool a_eq_e, float grayabsorb)
{
 // global index of the particle to load on the shared memory
  int m = blockIdx.x *blockDim.x; // + threadIdx.x;
 // if (m >=nP) return;

  __shared__ int minx,maxx,miny,maxy,reg;
  __shared__ float rfacr,radsq,stp, posx,posy;
  __shared__ cu_color e,q;

  int local_chunk_length = blockDim.x;
  if (blockIdx.x == gridDim.x -1) local_chunk_length = nP - blockIdx.x*blockDim.x;
  const int np = threadIdx.x; // local index of the particle to load on the shared memory
			      // and pixel number to process for each particle

  // load chunk of particles on the shared memory: each thread loads a particle (NOT CONVENIENT)
/*  if(m < nP)
  {
     ppos[np] = pos[m];     //particle m absolute end position
     p[np] = part[m];
  }
   __syncthreads();
*/

 //now do the rendering: each thread processes a pixel of particle i
  int x,y;
  for (int i=0; i<local_chunk_length; i++)
  {
    if (threadIdx.x == 0)
    {
      cu_particle_sim p = part[m+i];
      e.r = -p.e.r;   e.g = -p.e.g;  e.b = -p.e.b;
      posx = p.x; posy = p.y;
      rfacr = dparams.rfac*p.r;
      radsq = rfacr*rfacr;
      float sigma = dparams.h2sigma*p.r;
      stp = -1.f/(sigma*sigma);

      minx = int(posx-rfacr+1.f);
      minx=max(minx,0);
      maxx=int(posx+rfacr+1.f);
      maxx=min(maxx,dparams.xres);
      miny=int(posy-rfacr+1.f);
      miny=max(miny,0);
      maxy=int(posy+rfacr+1.f);
      maxy=min(maxy,dparams.yres);
      reg = (maxx-minx)*(maxy-miny);
    }
    __syncthreads();

    // render pixel np of particle i
    if (np < reg)
    {
      // relative starting position:
      unsigned long fpos = pos[m+i] - FragRendered - (unsigned long) (reg - np); 
      x = np/(maxy-miny) + minx;
      y = np%(maxy-miny) + miny;

      if (a_eq_e)
      {
        cu_fragment_AeqE        *fbuf;
        fbuf =(cu_fragment_AeqE*) buf;

        float dxsq = (x-posx)*(x-posx);
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
     }
     else
     {
       cu_fragment_AneqE       *fbuf1;
       fbuf1 =(cu_fragment_AneqE*)buf;
       q.r = e.r/(e.r+grayabsorb);
       q.g = e.g/(e.g+grayabsorb);
       q.b = e.b/(e.b+grayabsorb);
 
       float dxsq=(x-posx)*(x-posx);
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
     }
    }
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

// check for non-active and big particles to remove from the device
struct particle_notValid
  {
    __host__ __device__
    bool operator()(const char flag)
    {
      return ((flag==-1)||(flag==2));
    }
  };

// check for active big particles to copy back to the host
struct reg_notValid
  {
    __host__ __device__
    bool operator()(const char flag)
    {
      return (flag==2);
    }
  };

#endif

