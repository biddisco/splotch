// Implementation only for the case A=E

#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/remove.h>
#include <thrust/extrema.h>
#include <thrust/scan.h>

#include "cuda/CuRender.h"
#include "cuda/CuPolicy.h"

using namespace std;

void cu_draw_chunk(void *pinfo, COLOUR *Pic, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb)
{
  //get the input info
  thread_info *tInfo = (thread_info*)pinfo;
  tInfo->times.start("gpu_thread");

  int nParticle = tInfo->endP - tInfo->startP + 1;
  cout << "Rank " << mpiMgr.rank() << " - GPU " << tInfo->devID << " : Processing " << nParticle << " particles" << endl;

  //copy data particle to device memory
  tInfo->times.start("gcopy");
  cu_particle_sim *d_particle_data = (cu_particle_sim *) &((*particle_data)[tInfo->startP]);
  cu_copy_particles_to_device(d_particle_data, nParticle, gv);
  tInfo->times.stop("gcopy");

  //CUDA Transformation
  tInfo->times.start("gtransform");
  cu_transform(nParticle, gv);
  tInfo->times.stop("gtransform");

// DEBUG ---------------------------------------------------------------------------
/*  bool *active = new bool[nParticle];
  cudaMemcpy(active, gv->d_active, nParticle*sizeof(bool), cudaMemcpyDeviceToHost);
  long int cont=0;
  for (int i=0; i < nParticle; i++) if (active[i]) cont++;
  cout << " Number active particles = " << cont << endl;
  delete[] active;
  unsigned long *region = new unsigned long[nParticle];
  cudaMemcpy(region, gv->d_posInFragBuf, nParticle*sizeof(unsigned long), cudaMemcpyDeviceToHost);
  unsigned long cont = 0;
  for (int i=0; i < nParticle; i++) cont+=region[i];
  cout << "Tot fragm= " << cont << endl;
*/
//------------------------------------------------------------------------------------


  //Remove non-active particles
  tInfo->times.start("gfilter");
  thrust::device_ptr<bool> dev_ptr_flag((bool *) gv->d_active);
  thrust::device_ptr<cu_particle_sim> dev_ptr_pd((cu_particle_sim *) gv->d_pd);

  thrust::device_ptr<cu_particle_sim> new_end = thrust::remove_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, reg_notValid());
  int newParticle = new_end.get() - dev_ptr_pd.get();
  if( newParticle != nParticle )
  {
     cout << "Removed " << nParticle - newParticle << " particles " << endl;
  }

  thrust::device_ptr<unsigned long> dev_ptr_reg((unsigned long *) gv->d_posInFragBuf);
  thrust::remove_if(dev_ptr_reg, dev_ptr_reg+nParticle, dev_ptr_flag, reg_notValid());
  
  // max size in pixels of the particles
  thrust::device_ptr<unsigned long> max_size = thrust::max_element(dev_ptr_reg, dev_ptr_reg + newParticle);
  unsigned long *raw_ptr_max = thrust::raw_pointer_cast(max_size);
  unsigned int block_size;  
  cudaError_t error = cudaMemcpy(&block_size, raw_ptr_max, sizeof(unsigned long), cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) cout << "max_size Memcpy error!" << endl; 
  else cout << "max particle size = " << block_size << endl;

  unsigned int nB = (unsigned int) gv->policy->GetMaxBlockSize();
  if (block_size > nB) block_size = nB;

  //compute the starting position of each particle in the fragment buffer
  thrust::inclusive_scan(dev_ptr_reg, dev_ptr_reg + newParticle, dev_ptr_reg);
  tInfo->times.stop("gfilter");

//  tInfo->times.start("gcopy");
//  unsigned long *region = new unsigned long[nParticle];
//  cudaMemcpy(region, gv->d_posInFragBuf, nParticle*sizeof(unsigned long), cudaMemcpyDeviceToHost);
//  cout << "Tot fragments = " << region[newParticle-1] << endl;
//  tInfo->times.stop("gcopy");

  //CUDA Coloring
  tInfo->times.start("gcolor");
  cu_colorize(newParticle, gv);
  tInfo->times.stop("gcolor");

/* temporarily ignore sorting 191109.
// --------------------------------
// ----------- Sorting ------------
// --------------------------------
   it becomes complicated when using multiple threads with sorting
 
*/

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------

  //get parameters for rendering
  pair<int, int> res = gv->policy->GetResolution();
//  size_t nFBufInByte = gv->policy->GetFBufSize();

  //clear the device image
  size_t size_Im = res.first * res.second * sizeof(cu_color);
  error = cudaMemset(gv->d_pic,0,size_Im);
  //cout << cudaGetErrorString(cudaGetLastError()) << endl;

  // allocate fragment and index buffers
  int dimGrid = gv->policy->GetMaxGridSize(); 
  if (newParticle < dimGrid) dimGrid = newParticle;

  long maxNFrag = dimGrid*block_size;
  cu_allocateFragmentBuffer(maxNFrag, gv);
  // cout << "Fragment buffer size (num fragments) = " << maxNFrag << endl;

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
   tInfo->times.start("grender");
   cu_render1(dimGrid, block_size, dimGrid, End_cu_ps, nFragments2RenderOld, a_eq_e,
             (float) grayabsorb, gv);
   //cout << cudaGetErrorString(cudaGetLastError()) << endl;
   End_cu_ps += dimGrid;	// number of rendered particles = grid dimension 
   tInfo->times.stop("grender");
 
  tInfo->times.start("gcopy");
  cudaError_t error = cudaMemcpy(&nFragments2RenderNew, gv->d_posInFragBuf + End_cu_ps-1, sizeof(unsigned long), cudaMemcpyDeviceToHost);
  tInfo->times.stop("gcopy");
  if (error != cudaSuccess) cout << "nFragments Memcpy error!" << endl; 
  //else cout << "nFragments = " << nFragments2RenderNew - nFragments2RenderOld << endl;

   tInfo->times.start("gsort");
   thrust::sort_by_key(dev_ptr_Index, dev_ptr_Index + nFragments2RenderNew - nFragments2RenderOld, dev_ptr_FragBuf);
   tInfo->times.stop("gsort");

   tInfo->times.start("greduce");
   new_end_frag = thrust::reduce_by_key(dev_ptr_Index, dev_ptr_Index + nFragments2RenderNew - nFragments2RenderOld, dev_ptr_FragBuf, dev_ptr_Index, dev_ptr_FragBuf, binary_pred, sum_op());
   int npixels = new_end_frag.first.get() - dev_ptr_Index.get(); 
   tInfo->times.stop("greduce");
   //cout << "npixels = " << npixels << endl;

   tInfo->times.start("gcombine");
   cu_update_image(npixels, a_eq_e, gv);
 
  // cout << "Rank " << mpiMgr.rank() << " - GPU " << tInfo->devID << " : Rendered " << End_cu_ps << "/"<< newParticle << " particles" << endl;
   
   nFragments2RenderOld = nFragments2RenderNew;
   if (newParticle - End_cu_ps < dimGrid) dimGrid  = newParticle - End_cu_ps;
   cudaThreadSynchronize();
   //cout << cudaGetErrorString(cudaGetLastError()) << endl;
   tInfo->times.stop("gcombine");
  }
// delete[] region;

  // copy back the image
//  size_t size = size_Im*sizeof(cu_color);
  error = cudaMemcpy(Pic, gv->d_pic, size_Im, cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) cout << "Device Memcpy error!" << endl; 

  cu_endChunk(gv);
  tInfo->times.stop("gpu_thread");
}


