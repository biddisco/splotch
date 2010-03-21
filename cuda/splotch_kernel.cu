#ifndef SPLOTCH_KERNEL_H
#define SPLOTCH_KERNEL_H
/*
Try accelerating splotch with CUDA. July 2009.
Copyright things go here.
*/

#include "splotch_cuda.h"

//MACROs
#define Pi 3.14159265358979323846264338327950288
#define get_xy_from_sn(sn, xmin, ymin, ymax, x, y)\
        {int x1 =sn/(ymax-ymin); int y1 =sn-x1*(ymax-ymin);\
         x  =x1 +xmin; y  =y1 +ymin;}
#define get_sn_from_xy(x,y,maxy,miny, sn)\
    {sn =x*(maxy-miny) +y;}

/////////help functions///////////////////////////////////
__device__ float    my_asinh(float val)
  {
  return log(val+sqrt(1.+val*val));
  }

__device__ void my_normalize(float minv, float maxv, float &val)
  {
  if (minv!=maxv) val =  (val-minv)/(maxv-minv);
  }

__device__ void clamp (float minv, float maxv, float &val)
  {
  val = min(maxv, max(minv, val));
  }

//fetch a color from color table on device
__device__ cu_color get_color
  (int ptype, float val, cu_colormap_info info)
  {
  //copy things to local block memory
  __shared__ cu_color_map_entry *map;
  __shared__ int      mapSize;
  __shared__ int *ptype_points;
  __shared__ int ptypes;

  map =info.map;
  mapSize =info.mapSize;
  ptype_points =info.ptype_points;
  ptypes  =info.ptypes;

  cu_color        clr;
  clr.r =clr.g =clr.b =0.0;

  //first find the right entry for this ptype
  if (ptype>=ptypes)
    return clr; //invalid input
  int     start, end;
  start =ptype_points[ptype];
  if ( ptype == ptypes-1)//the last type
    end =mapSize-1;
  else
    end =ptype_points[ptype+1]-1;

  //search the section of this type to find the val
  int i=start;
  while ((val>map[i+1].val) && (i<end))
    ++i;

  float fract = (val-map[i].val)/(map[i+1].val-map[i].val);
  cu_color clr1=map[i].color, clr2=map[i+1].color;
  clr.r =clr1.r + fract*(clr2.r-clr1.r);
  clr.g =clr1.g + fract*(clr2.g-clr1.g);
  clr.b =clr1.b + fract*(clr2.b-clr1.b);

  return clr;
  }

__device__  float get_exp(float arg, cu_exptable_info d_exp_info)
  {
#if 0
  return exp(arg);
#else
  //fetch things to local
  __shared__  float   expfac;
  __shared__  float   *tab1, *tab2;
  __shared__  int     mask1, mask3, nbits;
  expfac  =d_exp_info.expfac;
  tab1    =d_exp_info.tab1;
  tab2    =d_exp_info.tab2;
  mask1   =d_exp_info.mask1;
  mask3   =d_exp_info.mask3;
  nbits   =d_exp_info.nbits;

  int iarg= (int)(arg*expfac);
  //  for final device code
  if (iarg&mask3)
    return (iarg<0) ? 1. : 0.;
  return tab1[iarg>>nbits]*tab2[iarg&mask1];
#endif
  }

__global__ void k_post_process(cu_color *pic, int n, cu_exptable_info exp_info)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n)
    m =n;

  //each pic[m] should do the same calc, so sequence does not matter!
  pic[m].r =1.0 -get_exp( pic[m].r, exp_info);
  pic[m].g =1.0 -get_exp( pic[m].g, exp_info);
  pic[m].b =1.0 -get_exp( pic[m].b, exp_info);
  }

__global__ void k_combine
  (int minx, int miny, int maxx, int maxy, int xres, int yres,
  cu_particle_splotch *p, int pStart, int pEnd, cu_fragment_AeqE *fbuf, cu_color *pic)
  {
  int m =blockIdx.x *blockDim.x + threadIdx.x;
  int n =(maxx-minx)*(maxy-miny);
  if (m >=n)
    m =n;

  //get global coordinate point(x,y) of this thread
  int point_x, point_y;
  get_xy_from_sn(m, minx, miny, maxy, point_x, point_y);

  //go through all particles, for each particle p if point(x,y) is in its region
  //p(minx,miny, maxx,maxy) do the following.
  //find the sequencial number sn1 in p(minx,miny, maxx,maxy), the fragment we are looking
  //for in fragment buffer is fragBuf[ sn1+p.posInFBuf ]
  //grab the fragment f(deltaR,deltaG,deltaB)
  //find the sequencial number sn2 of point(x,y) in the output pic.
  //pic[sn2] += f
  int sn1, sn2, local_x, local_y, fpos;
  for (int i=pStart; i<=pEnd; i++)
    {
    if ( point_x >=p[i].minx && point_x<p[i].maxx &&
         point_y >=p[i].miny && point_y<p[i].maxy)
      {
      local_x =point_x -p[i].minx;
      local_y =point_y -p[i].miny;
      get_sn_from_xy(local_x, local_y, p[i].maxy, p[i].miny,sn1);
      fpos =sn1 +p[i].posInFragBuf;

      get_sn_from_xy(point_x, point_y, yres,0, sn2);
      pic[sn2].r +=fbuf[fpos].deltaR;
      pic[sn2].g +=fbuf[fpos].deltaG;
      pic[sn2].b +=fbuf[fpos].deltaB;
      }
    }
  }

//device render function k_render1
__global__ void k_render1
  (cu_particle_splotch *p,  int startP,  int endP,
  void *buf, bool a_eq_e, float grayabsorb,
  cu_exptable_info d_exp_info)//, int xres, int yres)
  {
  //first get the index m of this thread
  int m, n=endP-startP;
  m =blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n)//m goes from 0 to n-1
    return;
  m +=startP;

  //make fbuf the right type
  cu_fragment_AeqE        *fbuf;
  cu_fragment_AneqE       *fbuf1;
  if (a_eq_e)
    fbuf =(cu_fragment_AeqE*) buf;
  else
    fbuf1 =(cu_fragment_AneqE*)buf;

  //now do the calc
  const float rfac=1.5;
  const float powtmp = pow(Pi,1./3.);
  const float sigma0=powtmp/sqrt(2*Pi);
  const float bfak=1./(2*sqrt(Pi)*powtmp);

  int x0s=0, y0s=0;
  float r=p[m].r;
  float posx=p[m].x, posy=p[m].y;
  posx-=x0s; posy-=y0s;
  float rfacr=rfac*r;

  cu_color a=p[m].a, e, q;
  if (!a_eq_e)
    {
    e=p[m].e;
    q.r=e.r/(a.r+grayabsorb);
    q.g=e.g/(a.g+grayabsorb);
    q.b=e.b/(a.b+grayabsorb);
    }

  float radsq = rfacr*rfacr;
  float prefac1 = -0.5/(r*r*sigma0*sigma0);
  float prefac2 = -0.5*bfak/p[m].ro;
  int minx, miny, maxx, maxy;
  minx =p[m].minx;    miny =p[m].miny;
  maxx =p[m].maxx;    maxy =p[m].maxy;
  unsigned int    fpos;
  fpos =p[m].posInFragBuf -p[startP].posInFragBuf;

  for (int x=minx; x<maxx; ++x)
    {
    float xsq=(x-posx)*(x-posx);
    for (int y=miny; y<maxy; ++y)
      {
      float dsq = (y-posy)*(y-posy) + xsq;
      if (dsq<radsq)
        {
        float fac = prefac2*get_exp(prefac1*dsq, d_exp_info);
        if (a_eq_e)
          {
          fbuf[fpos].deltaR = (fac*a.r);
          fbuf[fpos].deltaG = (fac*a.g);
          fbuf[fpos].deltaB = (fac*a.b);
          }
        else
          {
          float   exp;
          exp =get_exp(fac*a.r, d_exp_info);
          fbuf1[fpos].factorR =exp;
          fbuf1[fpos].deltaR  =q.r*(1.0-exp);
          exp =get_exp(fac*a.g, d_exp_info);
          fbuf1[fpos].factorG =exp;
          fbuf1[fpos].deltaG  =q.g*(1.0-exp);
          exp =get_exp(fac*a.b, d_exp_info);
          fbuf1[fpos].factorB =exp;
          fbuf1[fpos].deltaB  =q.b*(1.0-exp);
          }//if a_eq_e
        }//if dsq<radsq
      else
        {
        if (a_eq_e)
          {
          fbuf[fpos].deltaR =0.0;
          fbuf[fpos].deltaG =0.0;
          fbuf[fpos].deltaB =0.0;
          }
        else
          {
          fbuf1[fpos].deltaR =0.0;
          fbuf1[fpos].deltaG =0.0;
          fbuf1[fpos].deltaB =0.0;
          fbuf1[fpos].factorR =1.0;
          fbuf1[fpos].factorG =1.0;
          fbuf1[fpos].factorB =1.0;
          }
        }
      //for each (x,y)
      fpos++;
      }//y
    }//x
  }

//colorize by kernel
__global__ void k_colorize
  (cu_param_colorize *params, cu_particle_sim *p, int n, cu_particle_splotch *p2,
  cu_colormap_info info)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n)
    m =n;

  //now do the calc, p[m]--->p2[m]
  p2[m].isValid=false;
  if (p[m].z<=0 || p[m].z<=params->zminval || p[m].z>=params->zmaxval)
    return;

  float r=p[m].r;
  float posx=p[m].x, posy=p[m].y;

  float rfacr=params->rfac*r;

  int minx=int(posx-rfacr+1);
  if (minx>=params->res) return;
  minx=max(minx,0);

  int maxx=int(posx+rfacr+1);
  if (maxx<=0) return;
  maxx=min(maxx,params->res);
  if (minx>=maxx) return;

  int miny=int(posy-rfacr+1);
  if (miny>=params->ycut1) return;
  miny=max(miny,params->ycut0);

  int maxy=int(posy+rfacr+1);
  if (maxy<=params->ycut0) return;
  maxy=min(maxy,params->ycut1);
  if (miny>=maxy) return;

  //set region info to output the p2
  p2[m].minx =minx;  p2[m].miny =miny;
  p2[m].maxx =maxx;  p2[m].maxy =maxy;

  float col1=p[m].C1,col2=p[m].C2,col3=p[m].C3;
  clamp (0.0000001,0.9999999,col1);
  if (params->col_vector[p[m].type])
    {
    clamp (0.0000001,0.9999999,col2);
    clamp (0.0000001,0.9999999,col3);
    }
  float intensity=p[m].I;
  clamp (0.0000001,0.9999999,intensity);
  intensity *= params->brightness[p[m].type];

  cu_color e;
  if (params->col_vector[p[m].type])
    {
    e.r=col1*intensity;
    e.g=col2*intensity;
    e.b=col3*intensity;
    }
  else
    {
    e=get_color(p[m].type, col1, info);
    e.r *=intensity;
    e.g *=intensity;
    e.b *=intensity;
    }

  cu_color a=e;

  p2[m].isValid =true;
  p2[m].x =p[m].x;
  p2[m].y =p[m].y;
  p2[m].r =p[m].r;
  p2[m].ro=p[m].ro;
  p2[m].a=a;
  p2[m].e=e;
  }

//Range by kernel step 1
__global__ void k_range1(cu_param_range *pr, cu_particle_sim *p, int n)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n)
    m =n;

  //now do the calc
  //I, minint, maxint
  if (pr->log_int[p[m].type]) //could access invalid address under EMULATION
    p[m].I = log10(p[m].I);

  //C1, mincol, maxcol
  if (pr->log_col[p[m].type])
  p[m].C1 = log10(p[m].C1);
  if (pr->asinh_col[p[m].type])
    p[m].C1 = my_asinh(p[m].C1);

  //C2, C3, mincol, maxcol
  if (pr->col_vector[p[m].type])
    {
    if (pr->log_col[p[m].type])
      {
      p[m].C2 = log10(p[m].C2);
      p[m].C3 = log10(p[m].C3);
      }
    if (pr->asinh_col[p[m].type])
      {
      p[m].C2 = my_asinh(p[m].C2);
      p[m].C3 = my_asinh(p[m].C3);
      }
    }
  }

//Range by kernel step 2
__global__ void k_range2
  (cu_param_range *pr, cu_particle_sim *p, int n, int itype,
  float minval_int, float maxval_int,
  float minval_col, float maxval_col)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
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
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n)
    m =n;

  //copy parameters to __share__ local memory? later

  //now do x,y,z
  float x,y,z;
  x =p[m].x*ptrans->p[0] + p[m].y*ptrans->p[1] + p[m].z*ptrans->p[2] + ptrans->p[3];
  y =p[m].x*ptrans->p[4] + p[m].y*ptrans->p[5] + p[m].z*ptrans->p[6] + ptrans->p[7];
  z =p[m].x*ptrans->p[8] + p[m].y*ptrans->p[9] + p[m].z*ptrans->p[10]+ ptrans->p[11];
  p[m].x =x;
  p[m].y =y;
  p[m].z =z;

  //do ro and r
  float   xfac =ptrans->xfac;
  if (!ptrans->projection)
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

  if (ptrans->minhsmlpixel)
    if ((p[m].r <= 0.5) && (p[m].r >= 0.0))
      {
      p[m].r = 0.5;
      p[m].ro = p[m].r/(ptrans->res*.5*xfac);
      }
  }

#endif // #ifndef SPLOTCH_KERNEL_H
