/*
 * Copyright (c) 2010-2014
 *              Marzia Rivi (1), Tim Dykes (2)
 *               (1) University of Oxford
 *               (2) University of Portsmouth
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */



#ifndef __KERNEL__
#define __KERNEL__

 #include "cuda_kernel.cuh"
//help functions
//#define Pi 3.141592653589793238462643383279502884197
//#define MAXSIZE 1000

__constant__ cu_param dparams;
__constant__ cu_color_map_entry dmap[MAXSIZE];
__constant__ int ptype_points[10];

__device__ __forceinline__ void clamp (float minv, float maxv, float &val)
{
  val = min(maxv, max(minv, val));
}

__device__ __forceinline__   double my_asinh (double val)
{ return log(val+sqrt(1.+val*val)); }

//fetch a color from color table on device
__device__ __forceinline__ cu_color get_color(int ptype, float val, int map_size, int map_ptypes)
{
  //first find the right entry for this ptype
  int     start, end;
  start = ptype_points[ptype];
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

//Transform+coloring by kernel
__global__ void k_process(cu_particle_sim *p, int n, int mapSize, int types)
{
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n) return;

  int ptype = p[m].type;
  float r = p[m].r;
  float I = p[m].I;
  
#ifndef SPLOTCH_PARAVIEW

   float er = p[m].e.r;
   float eg = p[m].e.g;
   float eb = p[m].e.b;   

  // Normalization and clamping 
#ifndef NO_I_NORM
  // Norm and clamp I
    if (dparams.inorm_maxs[ptype]==dparams.inorm_mins[ptype])
      I = 1;
    else
      I = (max(dparams.inorm_mins[ptype],min(dparams.inorm_maxs[ptype],I))-dparams.inorm_mins[ptype])/(dparams.inorm_maxs[ptype]-dparams.inorm_mins[ptype]);
#endif
  // Norm and clamp er
    if (dparams.cnorm_maxs[ptype]==dparams.cnorm_mins[ptype])
      er = 1;
    else
      er = (max(dparams.cnorm_mins[ptype],min(dparams.cnorm_maxs[ptype],er))-dparams.cnorm_mins[ptype])/(dparams.cnorm_maxs[ptype]-dparams.cnorm_mins[ptype]);
  
  // If col_vector[t]
  // norm and clamp eg and eb
    if(dparams.col_vector[ptype])
    {
      if (dparams.cnorm_maxs[ptype]==dparams.cnorm_mins[ptype])
        eg = 1;
      else
        eg = (max(dparams.cnorm_mins[ptype],min(dparams.cnorm_maxs[ptype],er))-dparams.cnorm_mins[ptype])/(dparams.cnorm_maxs[ptype]-dparams.cnorm_mins[ptype]);

      if (dparams.cnorm_maxs[ptype]==dparams.cnorm_mins[ptype])
        eb = 1;
      else
        eb = (max(dparams.cnorm_mins[ptype],min(dparams.cnorm_maxs[ptype],er))-dparams.cnorm_mins[ptype])/(dparams.cnorm_maxs[ptype]-dparams.cnorm_mins[ptype]);
    }
#endif

  
  float x,y,z;
  x =p[m].x*dparams.p[0] + p[m].y*dparams.p[1] + p[m].z*dparams.p[2] + dparams.p[3];
  y =p[m].x*dparams.p[4] + p[m].y*dparams.p[5] + p[m].z*dparams.p[6] + dparams.p[7];
  z =p[m].x*dparams.p[8] + p[m].y*dparams.p[9] + p[m].z*dparams.p[10]+ dparams.p[11];

 // if(-z <= 0.0f){p[m].active = false;return;}; <-----------
 // if(-z >= 1e23){p[m].active = false;return;}; <---------- Z IS REVERSED FOR SOME REASON!?

  //do r
  float xfac2 = dparams.xfac;
  //const float   res2 = 0.5f*dparams.xres;
  //const float   ycorr = 0.5f*(dparams.yres-dparams.xres);
  if (!dparams.projection)
  {
    x = 0.5f*dparams.xres * (x+dparams.fovfct*dparams.dist)*xfac2;
    y = 0.5f*dparams.xres * (y+dparams.fovfct*dparams.dist)*xfac2 + 0.5f*(dparams.yres-dparams.xres);
  }
  else
  {
    xfac2=1.f/(dparams.fovfct*z);
    x = 0.5f*dparams.xres * (x+dparams.fovfct*z)*xfac2;
    y = 0.5f*dparams.xres * (y+dparams.fovfct*z)*xfac2 +  0.5f*(dparams.yres-dparams.xres);
  }

#ifdef SPLOTCH_CLASSIC
  I *= 0.5f*dparams.bfak/r;
  r*= sqrtf(2.f)*dparams.sigma0/dparams.h2sigma;  
#else
  //I *= 8.f/(Pi*r*r*r);  //SPH kernel normalization
  //I *= dparams.h2sigma*sqrtf(Pi)*r;  //integral through the center
  I *= 8.f*dparams.h2sigma/(sqrtf(Pi)*r*r);
#endif

  r *= 0.5f*dparams.xres*xfac2;
  const float rcorr= sqrtf(r*r + dparams.minrad_pix*dparams.minrad_pix)/r;
  r *= rcorr;
#ifdef SPLOTCH_CLASSIC
  I /= rcorr;
#else
  I /= rcorr*rcorr;
#endif
  I *= dparams.brightness[ptype];

  p[m].active = false;

  // compute region occupied by the partile
  //float raux=dparams.rfac;
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
  p[m].x = x;
  p[m].y = y;
  p[m].r = r;
  p[m].I = I;

//coloring
// get color, associated from physical quantity contained in e.r, from lookup table
#ifndef SPLOTCH_PARAVIEW
 cu_color e;
 e.r=er;
 e.g=eg;
 e.b=eb;

 if (!dparams.col_vector[ptype])
    e = get_color(ptype, e.r, mapSize, types);
#endif

 p[m].e.r *= I;
 p[m].e.g *= I;
 p[m].e.b *= I; 
}
 
 
// Calculates logs on device, asinh is commented out because if it is used
// it is done on host
__global__ void k_range(int nP, cu_particle_sim *p)
{

  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=nP) return;

  // Get current particle type
  int ptype = p[m].type;

  // Check if we need to log10 intensity
  if (dparams.log_int[ptype])
  { 
    if(p[m].I > 0)
        p[m].I = log10(p[m].I);
    else
        p[m].I = -38;
  }

  if (dparams.log_col[ptype])
  {
    if(p[m].e.r > 0)
      {
      p[m].e.r = log10(p[m].e.r);
      }
    else
      p[m].e.r =-38;
  }

  if (dparams.col_vector[ptype])
  {
    if (dparams.log_col[ptype])
    {
      p[m].e.g = log10(p[m].e.g);
      p[m].e.b = log10(p[m].e.b);
    }
  }

}

// --------------------------------------
// Render for full atomic implementation
// --------------------------------------
__global__ void k_render(int nP, cu_particle_sim *part, cu_color *pic)
{
  // Get index, double check its not out of bounds 
  // (launch parameters mean it shouldnt be...)
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >=nP) return;

  cu_particle_sim p = part[idx];
  if(p.active)
  {
    // Work out radial factor
    float rfacr = dparams.rfac*p.r;
    float radsq = rfacr*rfacr;
    float stp = -1.f/(dparams.h2sigma*dparams.h2sigma*p.r*p.r);

    // Get min and max pixels affected, clamp to image boundary
    int minx=int(p.x-rfacr+1.f);
    minx=max(minx,0);
    int maxx=int(p.x+rfacr+1.f);
    maxx=min(maxx,dparams.xres); 
    int miny=int(p.y-rfacr+1.f);
    miny=max(miny,0);
    int maxy=int(p.y+rfacr+1.f);
    maxy=min(maxy,dparams.yres);


    // For each pixel on x
    for(int x = minx; x < maxx; ++x)
    {
      
      // Work out x dist from centre and new yminmax 
      float dxsq = (x - p.x)*(x-p.x);
      float dy = sqrt(radsq-dxsq);
      int miny2=max(miny,int(p.y-dy+1)),
          maxy2=min(maxy,int(p.y+dy+1));
      float pre2 = __expf(stp*dxsq);
      // For each pixel on y
      for(int y = miny2; y < maxy2; ++y)
      {
          // Work out y dist from centre  
          float dysq = (y - p.y) * (y - p.y);
          float att = __expf(stp*dysq);
          // Update global image
          atomicAdd(&(pic[x*dparams.yres+y].r),-att*p.e.r*pre2);
          atomicAdd(&(pic[x*dparams.yres+y].g),-att*p.e.g*pre2);
          atomicAdd(&(pic[x*dparams.yres+y].b),-att*p.e.b*pre2);      

      }
    }
  }
}

#endif


