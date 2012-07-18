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
#include "splotch/splotch_host.h"
#include "cuda/CuRender.h"
#include "cuda/CuPolicy.h"

#define TILE_SIDEX 16	// x side dimension of the image tile, in terms of pixels
#define TILE_SIDEY 16	// y side dimension of the image tile, in terms of pixels
#define WIDTH_BOUND 8   // width of the boundary around the image tile 

using namespace std;

int cu_draw_chunk(wallTimerSet *times, int mydevID, cu_particle_sim *d_particle_data, int nParticle, COLOUR *Pic, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, float b_brightness, vector<COLOURMAP> &amap, paramfile &g_params)
{
  cudaError_t error;

  //copy data particle to device memory
  times->start("gcopy");
  cu_copy_particles_to_device(d_particle_data, nParticle, gv);
  times->stop("gcopy");

  //get parameters for rendering
  pair<int, int> res = gv->policy->GetResolution();

  //CUDA Transformation
  times->start("gprocess");
  cu_process(nParticle, gv);
  cudaThreadSynchronize();
  times->stop("gprocess");

  times->start("gfilter"); 
  thrust::device_ptr<cu_particle_sim> dev_ptr_pd((cu_particle_sim *) gv->d_pd);
  thrust::device_ptr<int> dev_ptr_flag((int *) gv->d_active);

  //Copy big particles to be processed by the host
  thrust::device_vector<cu_particle_sim> d_host_part(nParticle);
  thrust::device_vector<cu_particle_sim>::iterator end = thrust::copy_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, d_host_part.begin(), reg_notValid()); 
  int nHostPart = end - d_host_part.begin();

  //Remove non-active and big particles
  thrust::device_ptr<cu_particle_sim> new_end = thrust::remove_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, particle_notValid());
  int newParticle = new_end.get() - dev_ptr_pd.get();
  if( newParticle != nParticle )
  {
     cout << endl << "Eliminating inactive particles..." << endl;
     cout << newParticle+nHostPart << " particles left" << endl; 
     cout << nHostPart << " of them are processed by the host" << endl; 
  }
  thrust::device_ptr<unsigned long> dev_ptr_reg((unsigned long *) gv->d_posInFragBuf);
  thrust::remove_if(dev_ptr_reg, dev_ptr_reg+nParticle, dev_ptr_flag, particle_notValid());

  //sort particles according to their tile id
  times->start("gsort");
  //thrust::sort_by_key(dev_ptr_flag, dev_ptr_flag + nParticle, dev_ptr_pd);
  times->stop("gsort");
  
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

// Copy back big particles
  times->start("gcopy");
  particle_sim *host_part = 0;
  if (nHostPart > 0)
  {
    cu_particle_sim *d_host_part_ptr = thrust::raw_pointer_cast(&d_host_part[0]);
    error = cudaHostAlloc((void**) &host_part, nHostPart*sizeof(cu_particle_sim), cudaHostAllocDefault);
    if (error != cudaSuccess) cout << "cudaHostAlloc error!" << endl;
    else
    {
      error = cudaMemcpyAsync(host_part, d_host_part_ptr, nHostPart*sizeof(cu_particle_sim), cudaMemcpyDeviceToHost);
      if (error != cudaSuccess) cout << "Big particles Memcpy error!" << endl;
    }
  }
  times->stop("gcopy");
  //cudaFuncSetCacheConfig(k_colorize, cudaFuncCachePreferL1);

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------

  //clear the device image
  size_t size_Im = res.first * res.second * sizeof(cu_color);
  error = cudaMemset(gv->d_pic,0,size_Im);
  //int ntiles = (res.first/TILE_SIDEX) * (res.second/TILE_SIDEY);

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
             (float) grayabsorb, gv, TILE_SIDEX, TILE_SIDEY, WIDTH_BOUND);
   cudaThreadSynchronize();
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
  cout << "Rank " << mpiMgr.rank() << " : Device rendering on " << newParticle << " particles" << endl;

  // copy back the image
  times->start("gcopy");
  error = cudaMemcpy(Pic, gv->d_pic, size_Im, cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) cout << "Device Memcpy error!" << endl; 
  times->stop("gcopy");

  // host rendering 
  times->start("host_rendering");
  if(nHostPart > 0)
  {
     cout << "Rank " << mpiMgr.rank() << " : Host rendering on " << nHostPart << " particles" << endl;
     host_funct::render_new(host_part, nHostPart, Pic_host, a_eq_e, grayabsorb);
  }
  times->stop("host_rendering");

  cu_endChunk(gv);
  if (host_part) cudaFreeHost(host_part);
  return nHostPart+newParticle;
}
