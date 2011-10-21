#ifndef __KERNEL__
#define __KERNEL__

#include "cuda/splotch_cuda.h"

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
    end =map_size-1;
  else
    end =ptype_points[ptype+1]-1;

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
__device__ void k_colorize
  (cu_particle_splotch *p2, int mapSize, int types)
  {
  //first get the index m of this thread
  //int m=blockIdx.x *blockDim.x + threadIdx.x;
  //if (m >n) m =n;

  int ptype = p2->type;
  float col1=p2->e.r,col2=p2->e.g,col3=p2->e.b;
  clamp (0.0000001,0.9999999,col1);
  if (dparams.col_vector[ptype])
    {
    clamp (0.0000001,0.9999999,col2);
    clamp (0.0000001,0.9999999,col3);
    }
  float intensity=p2->I;
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
      { e.r =e.g =e.b =0.0; }
    }
  p2->e=e;

  }

//device render function k_render1
__global__ void k_render1
  (cu_particle_splotch *p, int nP, void *buf, int *index, 
   bool a_eq_e, float grayabsorb, int mapSize, int types, int yres)
{
  //first get the index m of this thread
  int m;
  m =blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=nP)//m goes from 0 to nP-1
    return;

  // coloring
  k_colorize(&(p[m]), mapSize, types);

  //make fbuf the right type
  cu_fragment_AeqE        *fbuf;
  cu_fragment_AneqE       *fbuf1;
  if (a_eq_e)
    fbuf =(cu_fragment_AeqE*) buf;
  else
    fbuf1 =(cu_fragment_AneqE*)buf;

  //now do the rendering
  const float powtmp = pow(Pi,1./3.);
  const float sigma0 = powtmp/sqrt(2*Pi);

  const float r = p[m].r;
  const float radsq = 2.25*r*r;
  const float stp = -0.5/(r*r*sigma0*sigma0);

  cu_color q, e=p[m].e;
  if (!a_eq_e)
   {
     q.r = e.r/(e.r+grayabsorb);
     q.g = e.g/(e.g+grayabsorb);
     q.b = e.b/(e.b+grayabsorb);
   }
  const float intens = -0.5/(2*sqrt(Pi)*powtmp);
  e.r*=intens; e.g*=intens; e.b*=intens;

  const float posx=p[m].x, posy=p[m].y;
  unsigned int fpos =p[m].posInFragBuf;

  if (a_eq_e)
  {
    for (int x=p[m].minx; x<p[m].maxx; ++x)
    {
     float dxsq=(x-posx)*(x-posx);
     for (int y=p[m].miny; y<p[m].maxy; ++y)
      {
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
          fbuf[fpos].aR =0.0;
          fbuf[fpos].aG =0.0;
          fbuf[fpos].aB =0.0;
        }
        index[fpos] = x*yres+y;
        fpos++;
      }//y
    }//x
  }
  else
  {
    for (int x=p[m].minx; x<p[m].maxx; ++x)
    {
     float dxsq=(x-posx)*(x-posx);
     for (int y=p[m].miny; y<p[m].maxy; ++y)
      {
        float dsq = (y-posy)*(y-posy) + dxsq;
        if (dsq<radsq)
        {
          float att = __expf(stp*dsq);
          float   expm1;
          expm1 =__expf(att*e.r)-1.0;
          fbuf1[fpos].aR = expm1;
          fbuf1[fpos].qR = q.r;
          expm1 =__expf(att*e.g)-1.0;
          fbuf1[fpos].aG = expm1;
          fbuf1[fpos].qG = q.g;
          expm1 =__expf(att*e.b)-1.0;
          fbuf1[fpos].aB = expm1;
          fbuf1[fpos].qB = q.b;
        }
        else
        {
          fbuf1[fpos].aR =0.0;
          fbuf1[fpos].aG =0.0;
          fbuf1[fpos].aB =0.0;
          fbuf1[fpos].qR =1.0;
          fbuf1[fpos].qG =1.0;
          fbuf1[fpos].qB =1.0;
        }
        index[fpos] = x*yres+y;
        fpos++;
      }//y
    }//x
  }
}

//Transform by kernel
__global__ void k_transform
  (cu_particle_sim *p, cu_particle_splotch *p2, int n)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n) return;

  //now do x,y,z
  float x,y,z;
  x =p[m].x*dparams.p[0] + p[m].y*dparams.p[1] + p[m].z*dparams.p[2] + dparams.p[3];
  y =p[m].x*dparams.p[4] + p[m].y*dparams.p[5] + p[m].z*dparams.p[6] + dparams.p[7];
  z =p[m].x*dparams.p[8] + p[m].y*dparams.p[9] + p[m].z*dparams.p[10]+ dparams.p[11];

  //do r
  float xfac = dparams.xfac;
  const float   res2 = 0.5*dparams.xres;
  const float   ycorr = .5f*(dparams.yres-dparams.xres);
  if (!dparams.projection)
    {
    x = res2 * (x+dparams.fovfct*dparams.dist)*xfac;
    y = res2 * (y+dparams.fovfct*dparams.dist)*xfac + ycorr;
    }
  else
    {
    xfac=1./(dparams.fovfct*z);
    x = res2 * (x+dparams.fovfct*z)*xfac;
    y = res2 * (y+dparams.fovfct*z)*xfac + ycorr;
    }

  float r = p[m].r;
  p[m].I /= r;
  r *= res2*xfac;

  const float rfac= sqrt(r*r + 0.25*dparams.minrad_pix*dparams.minrad_pix)/r;
  r *= rfac;
  p2[m].I = p[m].I/rfac;

  p2[m].x = x;
  p2[m].y = y;
  p2[m].r = r;

  p2[m].isValid = false;

  // compute region occupied by the partile
  const float rfacr=dparams.rfac*r;
  int minx=int(x-rfacr+1);
  if (minx>=dparams.xres) return;
  minx=max(minx,0);

  int maxx=int(x+rfacr+1);
  if (maxx<=0) return;
  maxx=min(maxx,dparams.xres);
  if (minx>=maxx) return;

  int miny=int(y-rfacr+1);
  if (miny>=dparams.yres) return;
  miny=max(miny,0);

  int maxy=int(y+rfacr+1);
  if (maxy<=0) return;
  maxy=min(maxy,dparams.yres);
  if (miny>=maxy) return; 

  p2[m].minx =minx;  p2[m].miny =miny;
  p2[m].maxx =maxx;  p2[m].maxy =maxy;


  p2[m].e.r = (float) p[m].e.r;
  p2[m].e.g = (float) p[m].e.g;
  p2[m].e.b = (float) p[m].e.b;
  p2[m].type = p[m].type;
  p2[m].isValid = true;
  }

__global__ void k_clear(int n, cu_color *pic)
{
   //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n) m =n;

   pic[m].r = 0.0;
   pic[m].g = 0.0;
   pic[m].b = 0.0;
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

#endif

