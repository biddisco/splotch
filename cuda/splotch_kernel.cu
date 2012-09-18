#ifndef __KERNEL__
#define __KERNEL__

#include "cuda/splotch_cuda.h"

//#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
//#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
//#endif


//MACROs
#define Pi 3.14159265358979323846264338327950288
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
__global__ void k_process
  (cu_particle_sim *p, int *p_active, int n, int mapSize, int types, int tile_sidex, int tile_sidey, int width, int nxtiles, int nytiles)
  {
  //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >=n) return;

  float r = p[m].r;
  float I = p[m].I;
  int ptype = p[m].type;
  cu_color e;
  e.r=p[m].e.r;
  e.g=p[m].e.g;
  e.b=p[m].e.b;

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

#ifdef SPLOTCH_CLASSIC
  I *= 0.5f*dparams.bfak/r;
  r*= sqrtf(2.f)*dparams.sigma0/dparams.h2sigma;  
#else
  I *= 8.f/(Pi*r*r*r);  //SPH kernel normalization
  I *= dparams.h2sigma*sqrtf(Pi)*r;  //integral through the center
#endif

  r *= res2*xfac2;
  const float rcorr= sqrtf(r*r + dparams.minrad_pix*dparams.minrad_pix)/r;
  r *= rcorr;
#ifdef SPLOTCH_CLASSIC
  I /= rcorr;
#else
  I /= rcorr*rcorr;
#endif
  I *= dparams.brightness[ptype];

  p[m].x = x;
  p[m].y = y;
  p[m].r = r;
  p[m].I = I;

//coloring
// get color, associated from physical quantity contained in e.r, from lookup table
  if (!dparams.col_vector[ptype])
     e = get_color(ptype, e.r, mapSize, types);

  p[m].e.r = e.r*I;
  p[m].e.g = e.g*I;
  p[m].e.b = e.b*I;

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
  // active particle = tile_id to which it belongs to
  p_active[m] = int(y)/tile_sidey + int(x)/tile_sidex*nytiles; 
  //if (p_active[m] < 0 || p_active[m] > nxtiles*nytiles) {printf("x=%f, y=%f, flag=%d\n",x,y,p_active[m]);}
  if ((maxx-minx)*(maxy-miny) <= 1) p_active[m] = nxtiles*nytiles; // point-like particles 
  if (int(rfacr+1.f)>width) 
  {
      p_active[m] = -2; // particle to be removed and copied back to the host 
     // printf("x=%f, y=%f, rfacr=%d\n",x,y,int(rfacr));
  }
}

//colorize by kernel
/*__global__ void k_colorize
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
*/


// change of linear coordinate: from tile to global image
// lpix -> (x,y) -> (X,Y) -> gpix 
__device__ int pixelLocalToGlobal(int lpix, int xo, int yo, int width, int tile_sidey)
{
  // global 2D coordinates
  int x = xo + lpix/(tile_sidey+2*width);
  int y = yo + lpix%(tile_sidey+2*width);

  return x*dparams.yres+y;
}

#define NPSIZE 256

//device render function k_render1
// a_eq_e = false is not supported
__global__ void k_render1
  (int nP, cu_particle_sim *part, int *tileId, int *tilepart, cu_color *pic, cu_color *pic1, cu_color *pic2, cu_color *pic3, int tile_sidex, int tile_sidey, int width, int nytiles)
{
   extern __shared__ cu_color Btile[];
   __shared__ int local_chunk_length, end;
   __shared__ cu_color e[NPSIZE];
   __shared__ float radsq[NPSIZE], stp[NPSIZE];
   __shared__ float posx[NPSIZE], posy[NPSIZE];
   __shared__ int minx[NPSIZE], maxx[NPSIZE], miny[NPSIZE], maxy[NPSIZE];
   
   int tileBsize = (tile_sidex+2*width)*(tile_sidey+2*width);
   int tile = tileId[blockIdx.x];	// tile number 

   if (threadIdx.x == 0)
   {
      end = tilepart[blockIdx.x];
      if (blockIdx.x == 0) local_chunk_length = end;
      else local_chunk_length = end - tilepart[blockIdx.x-1];
      end--;
   }
   __syncthreads();

   int xo = (tile/nytiles)*tile_sidex - width;  // Btile origin x
   int yo = (tile%nytiles)*tile_sidey - width;  // Btile origin y

  //inizialise Btile
  for (int i=threadIdx.x; i<tileBsize; i=i+blockDim.x) 
  {
     Btile[i].r = 0.0f;  Btile[i].g = 0.0f;   Btile[i].b = 0.0f;
  }
   __syncthreads();

  int x,y,k;
  int j = 0;
  int last = blockDim.x;
  //now do the rendering: each thread processes a pixel of particle i

  while (j < local_chunk_length) 
  {
      k = threadIdx.x; 
      if(j+k < local_chunk_length)
      {
        cu_particle_sim p = part[end-k-j];
        e[k] = p.e;
        posx[k] = p.x; posy[k] = p.y;
        float rfacr = dparams.rfac*p.r;
        radsq[k] = rfacr*rfacr;
        stp[k] = -1.f/(dparams.h2sigma*dparams.h2sigma*p.r*p.r);

        minx[k]=int(p.x-rfacr+1.f);
        minx[k]=max(minx[k],0);
        maxx[k]=int(p.x+rfacr+1.f);
        maxx[k]=min(maxx[k],dparams.xres); 
        miny[k]=int(p.y-rfacr+1.f);
        miny[k]=max(miny[k],0);
        maxy[k]=int(p.y+rfacr+1.f);
        maxy[k]=min(maxy[k],dparams.yres);
      }
      __syncthreads(); 

      j += blockDim.x;
      if (j > local_chunk_length) last = local_chunk_length%blockDim.x;
      for (int i=0; i<last; i++)
      {
         int reg = (maxx[i]-minx[i])*(maxy[i]-miny[i]);
         // render pixel threadIdx.x of particle i
         if (threadIdx.x < reg)
         {
           // global pixel coordinates
           x = threadIdx.x/(maxy[i]-miny[i]) + minx[i];
           y = threadIdx.x%(maxy[i]-miny[i]) + miny[i];
           // global pixel index = x*dparams.yres+y
           // localx = x-xo,   localy = y-yo 
           int lp = (x-xo)*(tile_sidey+2*width) + y-yo;  //local pixel index
       //    if (lp >= tileBsize) printf("lp = %d, tile=%d, x=%d, y=%d xr =%f \n",lp,tile,x,y,posx[i]);
           float dsq = (y-posy[i])*(y-posy[i]) + (x-posx[i])*(x-posx[i]);
           if (dsq<radsq[i])
           {
             float att = __expf(stp[i]*dsq);
             Btile[lp].r += -att*e[i].r;
             Btile[lp].g += -att*e[i].g;
             Btile[lp].b += -att*e[i].b;
           }
           else
           {
             Btile[lp].r += 0.0f;
             Btile[lp].g += 0.0f;
             Btile[lp].b += 0.0f;
           }
          }
      }
      __syncthreads();  
  }

  //update inner tile in the global image
  int k0 = width*(tile_sidey+2*width) + width; // starting point
  for (int i=threadIdx.x; i<tile_sidex*tile_sidey; i=i+blockDim.x) 
  {
     j = k0 + i + (i/tile_sidey)*2*width; //add correction due to the boundary
     pic[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
  }
  __syncthreads();

// update boundary in 3 steps: 
// 1. columns

  int ymax = yo + tile_sidey+2*width;
  int xmax = xo + tile_sidex+2*width;
  int step = blockDim.x/2;

  if ((threadIdx.x < step)  && (yo > 0))
  {
    k0 = width*(tile_sidey+2*width);
    for (int i = threadIdx.x; i<tile_sidex*width; i=i+step) 
    {
      j = k0 + i + (i/width)*(tile_sidey+width); //add correction due to the boundary
      pic1[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
    }
  }
  else if ((threadIdx.x >= step)  && (ymax < dparams.yres))
  {
    k0 = width*(tile_sidey+2*width) + width + tile_sidey; 
    for (int i = threadIdx.x - step; i<tile_sidex*width; i=i+step) 
    {
      j = k0 + i + (i/width)*(tile_sidey+width); //add correction due to the boundary
      pic1[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
    }
  }
  __syncthreads();

// 2. rows
  if ((threadIdx.x < step) && (xo > 0))
  {
    k0 = width; 
    for (int i=threadIdx.x; i<tile_sidey*width; i=i+step) 
    {
      j = k0 + i + (i/tile_sidey)*2*width; //add correction due to the boundary
      pic2[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
    }
  }
  else if ((threadIdx.x >= step)  && (xmax < dparams.xres))
  {
    k0 = width + (width+tile_sidex)*(tile_sidey+2*width); // starting point
    for (int i=threadIdx.x - step; i<tile_sidey*width; i=i+step) 
    {
      j = k0 + i + (i/tile_sidey)*2*width; //add correction due to the boundary
      pic2[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
    }
  }
  __syncthreads();

// 3. corners
// dimension corners = 1/4 dimension blocks
  int i;
  if ((threadIdx.x < blockDim.x/4) && (xo > 0) && (yo > 0))
  {
     j = threadIdx.x + (threadIdx.x/width)*(tile_sidey+width);
     pic3[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
  }
  else if ((threadIdx.x >= blockDim.x/4 && threadIdx.x < blockDim.x/2) && (xo > 0) && (ymax < dparams.yres))
  {
     k0 = width + tile_sidey; 
     i = threadIdx.x - blockDim.x/4; 
     j = k0 + i + (i/width)*(tile_sidey+width);
     pic3[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
  }
  else if ((threadIdx.x >= blockDim.x/2 && threadIdx.x < 3*blockDim.x/4) && (xmax < dparams.xres) && (yo > 0))
  {
     k0 = (width + tile_sidex)*(tile_sidey+2*width);
     i = threadIdx.x - blockDim.x/2; 
     j = k0 + i + (i/width)*(tile_sidey+width);
     pic3[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
  }
  else if ((threadIdx.x >= 3*blockDim.x/4) && (xmax < dparams.xres) && (ymax < dparams.yres))
  {
     k0 = (width + tile_sidex)*(tile_sidey+2*width) + width + tile_sidey;
     i = threadIdx.x - 3*blockDim.x/4; 
     j = k0 + i + (i/width)*(tile_sidey+width);
     pic3[pixelLocalToGlobal(j,xo,yo,width,tile_sidey)] = Btile[j];
  }

}


//device render function k_renderC3
// each thread render a particle
__global__ void k_renderC3
  (int n, cu_particle_sim *part, int *index)
{
   //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >= n) return;
 
  cu_particle_sim p = part[m];
  float rfacr = dparams.rfac*p.r;
  int x = int(p.x-rfacr+1.f);
  x = max(x,0);
  int y = int(p.y-rfacr+1.f);
  y = max(y,0);
  index[m] = x*dparams.yres+y; 
  //pixel = -p.e
}

__global__ void k_add_images(int n, cu_color *pic, cu_color *pic1, cu_color *pic2, cu_color *pic3)
{
   //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >n) return;

   pic[m].r += pic1[m].r + pic2[m].r + pic3[m].r;
   pic[m].g += pic1[m].g + pic2[m].g + pic3[m].g;
   pic[m].b += pic1[m].b + pic2[m].b + pic3[m].b;
}

__global__ void k_addC3(int nC3, int *index, cu_particle_sim *part, cu_color *pic)
{
   //first get the index m of this thread
  int m=blockIdx.x *blockDim.x + threadIdx.x;
  if (m >= nC3) return;

  pic[index[m]].r += - part[m].e.r;
  pic[index[m]].g += - part[m].e.g;
  pic[index[m]].b += - part[m].e.b;
  
}

// check for non-active and big particles to remove from the device
struct particle_notValid
  {
    __host__ __device__
    bool operator()(const int flag)
    {
      return (flag < 0);
    }
  };

// check for active big particles to copy back to the host
struct reg_notValid
  {
    __host__ __device__
    bool operator()(const int flag)
    {
      return (flag==-2);
    }
  };

struct sum_op
{
  __host__ __device__
  cu_particle_sim operator()(cu_particle_sim& p1, cu_particle_sim& p2) const{

    cu_particle_sim sum;
    sum = p1;
    sum.e.r = p1.e.r + p2.e.r;
    sum.e.g = p1.e.g + p2.e.g;
    sum.e.b = p1.e.b + p2.e.b;

    return sum; 
   } 
};

#endif

