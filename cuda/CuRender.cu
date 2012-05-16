// Implementation only for the case A=E

#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/remove.h>
#include <thrust/copy.h>
#include <thrust/extrema.h>
#include <thrust/scan.h>

#include "splotch/splotchutils.h"
#include "cuda/CuRender.h"
#include "cuda/CuPolicy.h"

using namespace std;

void cu_draw_chunk(wallTimerSet *times, int mydevID, int startP, int endP, COLOUR *Pic, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, float b_brightness)
{
  cudaError_t error;

  times->start("gpu_thread");
  int nParticle = endP - startP;

  //copy data particle to device memory
  times->start("gcopy");
  cu_particle_sim *d_particle_data = (cu_particle_sim *) &((*particle_data)[startP]);
  cu_copy_particles_to_device(d_particle_data, nParticle, gv);
  times->stop("gcopy");

  //get parameters for rendering
  pair<int, int> res = gv->policy->GetResolution();

  //CUDA Transformation
  times->start("gtransform");
  cu_transform(nParticle, gv);
  times->stop("gtransform");

  times->start("gfilter"); 
  thrust::device_ptr<cu_particle_sim> dev_ptr_pd((cu_particle_sim *) gv->d_pd);
  thrust::device_ptr<char> dev_ptr_flag((char *) gv->d_active);

  //Copy big particles to be processed by the host
  thrust::device_vector<cu_particle_sim> d_host_part(nParticle);
  thrust::device_vector<cu_particle_sim>::iterator end = thrust::copy_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, d_host_part.begin(), reg_notValid()); 
  int nHostPart = end - d_host_part.begin();

// DEBUG ---------------------------------------------------------------------------
/*  unsigned long *region = new unsigned long[nParticle];
  cudaMemcpy(region, gv->d_posInFragBuf, nParticle*sizeof(unsigned long), cudaMemcpyDeviceToHost);
  int cont = 0;
  for (int i=0; i < nParticle; i++) if (region[i] > 1024) cont++;
  cout << cont << " big particles" << endl;
  delete[] region; 
  int cont=0;
  for (int i=0; i < nHostPart; i++) if (host_part[i].active) cont++;
  cout << cont << " big active particles" << endl;
*/
//------------------------------------------------------------------------------------

  //Remove non-active and big particles
  thrust::device_ptr<cu_particle_sim> new_end = thrust::remove_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, particle_notValid());
  int newParticle = new_end.get() - dev_ptr_pd.get();
  if( newParticle != nParticle )
  {
     cout << "Eliminating inactive particles..." << endl;
     cout << newParticle+nHostPart << " particles left" << endl; 
     cout << nHostPart << " of them are processed by the host" << endl; 
  }
  thrust::device_ptr<unsigned long> dev_ptr_reg((unsigned long *) gv->d_posInFragBuf);
  thrust::remove_if(dev_ptr_reg, dev_ptr_reg+nParticle, dev_ptr_flag, particle_notValid());
  
  // max size in pixels of the particles
  thrust::device_ptr<unsigned long> max_size = thrust::max_element(dev_ptr_reg, dev_ptr_reg + newParticle);
  unsigned long *raw_ptr_max = thrust::raw_pointer_cast(max_size);
  unsigned int max_reg_size;
  error = cudaMemcpy(&max_reg_size, raw_ptr_max, sizeof(unsigned long), cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) cout << "max_size Memcpy error!" << endl; 
  else cout << "max particle size = " << max_reg_size << endl;

  //compute the starting position of each particle in the fragment buffer
  thrust::inclusive_scan(dev_ptr_reg, dev_ptr_reg + newParticle, dev_ptr_reg);
  times->stop("gfilter");

// Overlap big particles copy and coloring kernel
  cudaStream_t stream0, stream1;
  cudaStreamCreate(&stream0);
  cudaStreamCreate(&stream1);

  times->start("gcolor");
  cu_particle_sim *host_part = 0;
  if (nHostPart > 0)
  {
    cu_particle_sim *d_host_part_ptr = thrust::raw_pointer_cast(&d_host_part[0]);
    error = cudaHostAlloc((void**) &host_part, nHostPart*sizeof(cu_particle_sim), cudaHostAllocDefault);
    if (error != cudaSuccess) cout << "cudaHostAlloc error!" << endl;
    else
    {
      error = cudaMemcpyAsync(host_part, d_host_part_ptr, nHostPart*sizeof(cu_particle_sim), cudaMemcpyDeviceToHost,stream1);
      if (error != cudaSuccess) cout << "Big particles Memcpy error!" << endl;
    }
  }
  //CUDA Coloring
  dim3 Grid, Block;
  gv->policy->GetDimsBlockGrid(newParticle, &Grid, &Block);
  k_colorize<<<Grid,Block,0,stream0>>>(gv->d_pd, gv->colormap_size, gv->colormap_ptypes, newParticle);
//  cu_colorize(newParticle, gv);
  cudaStreamSynchronize(stream0);
  cout << "Rank " << mpiMgr.rank() << " - GPU " << mydevID << " :  Coloring " << newParticle << "/"<< newParticle << " particles" << endl;
  cudaStreamSynchronize(stream1);
  times->stop("gcolor");

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------

  //clear the device image
  size_t size_Im = res.first * res.second * sizeof(cu_color);
  error = cudaMemset(gv->d_pic,0,size_Im);

  // allocate fragment and index buffers
  int fbsize = gv->policy->GetFBufSize();
  // number of threads in each block = max number of pixels to be rendered for a particle
  int block_size = (int) max_reg_size;
  // each block of threads process n = block_size particles: 
  int nB = (int) fbsize/(sizeof(cu_color)*max_reg_size*block_size);
  int dimGrid = gv->policy->GetMaxGridSize();
  if (dimGrid > nB) dimGrid = nB;

  int chunk_length = dimGrid*block_size; //number of rendered particles at each iteration
  if (newParticle < chunk_length) 
  {
	dimGrid = newParticle/block_size;
        if(newParticle%block_size) dimGrid++;
        chunk_length = newParticle;
  }
  cout << "dimGrid = " << dimGrid << endl;

  // number of fragments = n particles * max_particle_size
  long maxNFrag = chunk_length*max_reg_size;
  cu_allocateFragmentBuffer(maxNFrag, gv);
  cout << "Fragment buffer size = " << maxNFrag << endl;

  // wrap raw pointer with a device_ptr 
  thrust::device_ptr<cu_fragment_AeqE> dev_ptr_FragBuf((cu_fragment_AeqE *) gv->d_fbuf);
  thrust::device_ptr<int> dev_ptr_Index(gv->d_pixel);
  thrust::pair< thrust::device_ptr<int>,thrust::device_ptr<cu_fragment_AeqE> >  new_end_frag;
  thrust::equal_to<int> binary_pred;

  int End_cu_ps = 0; 
  unsigned long nFragments2RenderOld, nFragments2RenderNew;
  nFragments2RenderOld = 0;
  while (End_cu_ps < newParticle)
  {
   // Render particles
   // 1 block ----> 1 particle , 1 thread ----> 1 pixel
   times->start("grender");
   // device rendering
   cu_render1(dimGrid, block_size, chunk_length, End_cu_ps, nFragments2RenderOld, a_eq_e,
             (float) grayabsorb, gv);
   cudaThreadSynchronize();
   //cout << cudaGetErrorString(cudaGetLastError()) << endl;
   End_cu_ps += chunk_length;
   times->stop("grender");
 
  times->start("gcopy");
  error = cudaMemcpy(&nFragments2RenderNew, gv->d_posInFragBuf + End_cu_ps-1, sizeof(unsigned long), cudaMemcpyDeviceToHost);
  times->stop("gcopy");
  if (error != cudaSuccess) cout << "nFragments Memcpy error!" << endl; 

   times->start("gsort");
   thrust::sort_by_key(dev_ptr_Index, dev_ptr_Index + nFragments2RenderNew - nFragments2RenderOld, dev_ptr_FragBuf);
   times->stop("gsort");

   times->start("greduce");
   new_end_frag = thrust::reduce_by_key(dev_ptr_Index, dev_ptr_Index + nFragments2RenderNew - nFragments2RenderOld, dev_ptr_FragBuf, dev_ptr_Index, dev_ptr_FragBuf, binary_pred, sum_op());
   int npixels = new_end_frag.first.get() - dev_ptr_Index.get(); 
   times->stop("greduce");

   times->start("gcombine");
   cu_update_image(npixels, a_eq_e, gv);
   
   nFragments2RenderOld = nFragments2RenderNew;
   if (newParticle - End_cu_ps < chunk_length) 
   {
      chunk_length  = newParticle - End_cu_ps;
      dimGrid = chunk_length/block_size;
      if(chunk_length%block_size) dimGrid++;
   }
   cudaThreadSynchronize();
   //cout << cudaGetErrorString(cudaGetLastError()) << endl;
   times->stop("gcombine");
  }
  cout << "Rank " << mpiMgr.rank() << " - GPU " << mydevID << " : Rendered " << newParticle << "/"<< newParticle << " particles" << endl;

  // copy back the image
  times->start("gcopy");
  error = cudaMemcpy(Pic, gv->d_pic, size_Im, cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) cout << "Device Memcpy error!" << endl; 
  times->start("gcopy");

  times->stop("gpu_thread");

  // host rendering 
  times->start("host_rendering");
  Pic_host.fill(COLOUR(0,0,0));
  if (nHostPart > 0)
  {
     cout << "Rank " << mpiMgr.rank() << ": host Coloring+Rendering " << nHostPart << " particles" << endl;
     host_particle_colorize(*g_params, host_part, nHostPart, amap, b_brightness);
     host_render_new (host_part, nHostPart, Pic_host, a_eq_e, grayabsorb);
  }
  times->stop("host_rendering");

  cu_endChunk(gv);
  if (host_part) cudaFreeHost(host_part);
  cudaStreamDestroy(stream0);
  cudaStreamDestroy(stream1);
}



#define SPLOTCH_CLASSIC

const float32 h2sigma = 0.5*pow(pi,-1./6.);

#ifdef SPLOTCH_CLASSIC
const float32 powtmp = pow(pi,1./3.);
const float32 sigma0 = powtmp/sqrt(2*pi);
const float32 rfac=1.5*h2sigma/(sqrt(2.)*sigma0);
#else
const float32 rfac=1.;
#endif
 
void host_particle_colorize(paramfile &params, cu_particle_sim *p, int npart,
  vector<COLOURMAP> &amap, float b_brightness)
  {
  int nt = params.find<int>("ptypes",1);
  arr<bool> col_vector(nt);
  arr<float32> brightness(nt);

  for(int t=0;t<nt;t++)
    {
    brightness[t] = params.find<float32>("brightness"+dataToString(t),1.f);
    brightness[t] *= b_brightness;
    col_vector[t] = params.find<bool>("color_is_vector"+dataToString(t),false);
    }

#pragma omp parallel
{
  int m;
  COLOUR e;
#pragma omp for schedule(guided,1000)
  for (m=0; m<npart; ++m)
    {
      if (!col_vector[p[m].type])
      {
        e= amap[p[m].type].getVal_const(p[m].e.r);
        p[m].e.r = e.r;
 	p[m].e.g = e.g;
 	p[m].e.b = e.b;
      }
      p[m].e.r *= p[m].I * brightness[p[m].type];
      p[m].e.g *= p[m].I * brightness[p[m].type];
      p[m].e.b *= p[m].I * brightness[p[m].type];
    }
}
  }

const int chunkdim=100;

void host_render_new (cu_particle_sim *p, int npart, arr2<COLOUR> &pic, bool a_eq_e, float32 grayabsorb)
  {
  int xres=pic.size1(), yres=pic.size2();
  int ncx=(xres+chunkdim-1)/chunkdim, ncy=(yres+chunkdim-1)/chunkdim;

  arr2<vector<uint32> > idx(ncx,ncy);
  float32 rcell=sqrt(2.f)*(chunkdim*0.5f-0.5f);
  float32 cmid0=0.5f*(chunkdim-1);

  exptable<float32> xexp(-20.);
#ifdef PLANCK_HAVE_SSE
  const float32 taylorlimit=xexp.taylorLimit();
#endif

//  pic.fill(COLOUR(0,0,0));

//  tstack_push("Chunk preparation");

  for (int i=0; i<npart; ++i)
  {
    cu_particle_sim pp = p[i];
    float32 rfacr = rfac*pp.r;

    int minx=max(0,int(pp.x-rfacr+1)/chunkdim);
    int maxx=min(ncx-1,int(pp.x+rfacr)/chunkdim);
    int miny=max(0,int(pp.y-rfacr+1)/chunkdim);
    int maxy=min(ncy-1,int(pp.y+rfacr)/chunkdim);
    float32 sumsq=(rcell+rfacr)*(rcell+rfacr);
    for (int ix=minx; ix<=maxx; ++ix)
    {
        float32 cx=cmid0+ix*chunkdim;
        for (int iy=miny; iy<=maxy; ++iy)
          {
          float32 cy=cmid0+iy*chunkdim;
          float32 rtot2 = (pp.x-cx)*(pp.x-cx) + (pp.y-cy)*(pp.y-cy);
          if (rtot2<sumsq)
            idx[ix][iy].push_back(i);
          }
     }

   }

//  tstack_replace("Chunk preparation","Rendering proper");

  work_distributor wd (xres,yres,chunkdim,chunkdim);
#pragma omp parallel
{
  arr<float32> pre1(chunkdim);
#ifdef PLANCK_HAVE_SSE
  arr2_align<V4sf,16> lpic(chunkdim,chunkdim);
#else
  arr2<COLOUR> lpic(chunkdim,chunkdim);
#endif
  int chunk;
#pragma omp for schedule(dynamic,1)
  for (chunk=0; chunk<wd.nchunks(); ++chunk)
    {
    int x0, x1, y0, y1;
    wd.chunk_info(chunk,x0,x1,y0,y1);
    int x0s=x0, y0s=y0;
    x1-=x0; x0=0; y1-=y0; y0=0;
    lpic.fast_alloc(x1-x0,y1-y0);
#ifdef PLANCK_HAVE_SSE
    lpic.fill(V4sf(0.));
#else
    lpic.fill(COLOUR(0,0,0));
#endif
    int cx, cy;
    wd.chunk_info_idx(chunk,cx,cy);
    const vector<uint32> &v(idx[cx][cy]);

    for (tsize m=0; m<v.size(); ++m)
      {
      const cu_particle_sim pp=p[v[m]];
      float32 rfacr=pp.r*rfac;
      float32 posx=pp.x, posy=pp.y;
      posx-=x0s; posy-=y0s;
      int minx=int(posx-rfacr+1);
      minx=max(minx,x0);
      int maxx=int(posx+rfacr+1);
      maxx=min(maxx,x1);
      int miny=int(posy-rfacr+1);
      miny=max(miny,y0);
      int maxy=int(posy+rfacr+1);
      maxy=min(maxy,y1);

      float32 radsq = rfacr*rfacr;
      float32 sigma = h2sigma*pp.r;
      float32 stp = -1.f/(sigma*sigma);

      COLOUR a(-pp.e.r,-pp.e.g,-pp.e.b);
#ifdef PLANCK_HAVE_SSE
      V4sf va(a.r,a.g,a.b,0.f);
#endif

      for (int y=miny; y<maxy; ++y)
        pre1[y]=xexp(stp*(y-posy)*(y-posy));

      if (a_eq_e)
        {
        for (int x=minx; x<maxx; ++x)
          {
          float32 dxsq=(x-posx)*(x-posx);
          float32 dy=sqrt(radsq-dxsq);
          int miny2=max(miny,int(posy-dy+1)),
              maxy2=min(maxy,int(posy+dy+1));
          float32 pre2 = xexp(stp*dxsq);
          for (int y=miny2; y<maxy2; ++y)
            {
            float32 att = pre1[y]*pre2;
#ifdef PLANCK_HAVE_SSE
            lpic[x][y]+=va*att;
#else
            lpic[x][y].r += att*a.r;
            lpic[x][y].g += att*a.g;
            lpic[x][y].b += att*a.b;
#endif
            }
          }
        }
      else
        {
        COLOUR q(pp.e.r/(pp.e.r+grayabsorb),
                 pp.e.g/(pp.e.g+grayabsorb),
                 pp.e.b/(pp.e.b+grayabsorb));
#ifdef PLANCK_HAVE_SSE
        float32 maxa=max(abs(a.r),max(abs(a.g),abs(a.b)));
        V4sf vq(q.r,q.g,q.b,0.f);
#endif

        for (int x=minx; x<maxx; ++x)
          {
          float32 dxsq=(x-posx)*(x-posx);
          float32 dy=sqrt(radsq-dxsq);
          int miny2=max(miny,int(posy-dy+1)),
              maxy2=min(maxy,int(posy+dy+1));
          float32 pre2 = xexp(stp*dxsq);
          for (int y=miny2; y<maxy2; ++y)
            {
            float32 att = pre1[y]*pre2;
#ifdef PLANCK_HAVE_SSE
            if ((maxa*att)<taylorlimit)
              lpic[x][y]+=(lpic[x][y]-vq)*va*att;
            else
              {
              V4sf::Tu tmp;
              tmp.v=lpic[x][y].v;
              tmp.d[0] += xexp.expm1(att*a.r)*(tmp.d[0]-q.r);
              tmp.d[1] += xexp.expm1(att*a.g)*(tmp.d[1]-q.g);
              tmp.d[2] += xexp.expm1(att*a.b)*(tmp.d[2]-q.b);
              lpic[x][y]=tmp.v;
              }
#else
            lpic[x][y].r += xexp.expm1(att*a.r)*(lpic[x][y].r-q.r);
            lpic[x][y].g += xexp.expm1(att*a.g)*(lpic[x][y].g-q.g);
            lpic[x][y].b += xexp.expm1(att*a.b)*(lpic[x][y].b-q.b);
#endif
            }
          }
        }
      } // for particle

    for (int ix=0;ix<x1;ix++)
      for (int iy=0;iy<y1;iy++)
#ifdef PLANCK_HAVE_SSE
        {
        COLOUR &c(pic[ix+x0s][iy+y0s]); float32 dum;
        lpic[ix][iy].writeTo(c.r,c.g,c.b,dum);
        }
#else
        pic[ix+x0s][iy+y0s]=lpic[ix][iy];
#endif
    } // for this chunk
} // #pragma omp parallel

//  tstack_pop("Rendering proper");
  }

