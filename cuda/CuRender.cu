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
#include <thrust/iterator/constant_iterator.h>
#include <thrust/system_error.h>

#include "splotch/splotchutils.h"
#include "splotch/splotch_host.h"
#include "cuda/CuRender.h"
#include "cuda/CuPolicy.h"

using namespace std;

int cu_draw_chunk(int mydevID, cu_particle_sim *d_particle_data, int nParticle, COLOUR *Pic, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, int xres, int yres)
{
  cudaError_t error;

  //copy data particle to device memory
  tstack_push("Data copy");
  cu_copy_particles_to_device(d_particle_data, nParticle, gv); 
  tstack_pop("Data copy");

  //get parameters for rendering
  int tile_sidex, tile_sidey, width, nxtiles, nytiles;
  gv->policy->GetTileInfo(&tile_sidex, &tile_sidey, &width, &nxtiles, &nytiles);
  
  //--------------------------------------
  //  particle projection and coloring
  //--------------------------------------

  tstack_push("Particle projection & coloring");
  cu_process(nParticle, gv, tile_sidex, tile_sidey, width, nxtiles, nytiles);
  cudaThreadSynchronize();
  //cout << cudaGetErrorString(cudaGetLastError()) << endl;
  tstack_pop("Particle projection & coloring");

  int new_ntiles, newParticle, nHostPart;
  particle_sim *host_part = 0;
  try
  { 
   tstack_push("Particle Filtering");

   thrust::device_ptr<cu_particle_sim> dev_ptr_pd((cu_particle_sim *) gv->d_pd);
   thrust::device_ptr<int> dev_ptr_flag((int *) gv->d_active);

   // Select big particles to be processed by the host
   thrust::device_vector<cu_particle_sim> d_host_part(nParticle);
   thrust::device_vector<cu_particle_sim>::iterator end = thrust:: copy_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, d_host_part.begin(), reg_notValid()); 
   nHostPart = end - d_host_part.begin();

   // Copy back big particles
   if (nHostPart > 0)
   {
    cu_particle_sim *d_host_part_ptr = thrust::raw_pointer_cast(&d_host_part[0]);
    error = cudaHostAlloc((void**) &host_part, nHostPart*sizeof(cu_particle_sim), cudaHostAllocDefault);
    if (error != cudaSuccess) cout << "cudaHostAlloc error!" << endl;
    else
    {
      error = cudaMemcpy(host_part, d_host_part_ptr, nHostPart*sizeof(cu_particle_sim), cudaMemcpyDeviceToHost);
      if (error != cudaSuccess) cout << "Big particles Memcpy error!" << endl;
    }
   }

   //Remove non-active and host particles
   thrust::device_ptr<cu_particle_sim> new_end = thrust::remove_if(dev_ptr_pd, dev_ptr_pd+nParticle, dev_ptr_flag, particle_notValid());
   newParticle = new_end.get() - dev_ptr_pd.get();
   if( newParticle != nParticle )
   {
     cout << endl << "Eliminating inactive particles..." << endl;
     cout << newParticle+nHostPart << " particles left" << endl; 
     cout << nHostPart << " of them are processed by the host" << endl; 
     thrust::remove_if(dev_ptr_flag, dev_ptr_flag+nParticle, particle_notValid());
   }
   tstack_pop("Particle Filtering");

   tstack_push("Particle Distribution");

   //sort particles according to their tile id
   thrust::sort_by_key(dev_ptr_flag, dev_ptr_flag + newParticle, dev_ptr_pd);

   //compute number of particles for each tile and their starting position
   cudaMemset(gv->d_tiles,0,nxtiles*nytiles);
   thrust::device_ptr<int> dev_ptr_nT((int *) gv->d_tiles);
   thrust::pair< thrust::device_ptr<int>,thrust::device_ptr<int> > end_tiles = thrust::reduce_by_key(dev_ptr_flag, dev_ptr_flag + newParticle, thrust::make_constant_iterator(1), dev_ptr_flag, dev_ptr_nT);
   new_ntiles = end_tiles.second.get() - dev_ptr_nT.get();
   cout << "number of tiles = " << new_ntiles << endl;
   thrust::inclusive_scan(dev_ptr_nT, dev_ptr_nT + new_ntiles, dev_ptr_nT);

   tstack_pop("Particle Distribution");
  }
  catch(thrust::system_error &e)
  {
    // output an error message and exit
    std::cerr << "Error accessing vector element: " << e.what() << std::endl;
    exit(-1);
  }
  catch(std::bad_alloc &e)
  {
    std::cerr << "Couldn't allocate vector" << std::endl;
    exit(-1);
  }

  // ----------------------------
  //   particle proper rendering 
  // ----------------------------

  //clear the device image
  size_t size_Im = xres * yres * sizeof(cu_color);
  cudaMemset(gv->d_pic,0,size_Im);
  cudaMemset(gv->d_pic1,0,size_Im);
  cudaMemset(gv->d_pic2,0,size_Im);
  cudaMemset(gv->d_pic3,0,size_Im);

  // number of threads in each block = max number of pixels to be rendered for a particle
  int block_size = 4*width*width;
  int dimGrid = new_ntiles;    // number of blocks = number of tiles

  // Device rendering
  // 1 block ----> loop on chunk of particles, 1 thread ----> 1 pixel of the particle
  tstack_push("CUDA Rendering");

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start,0);
  cu_render1(dimGrid, block_size, a_eq_e, (float) grayabsorb, gv, tile_sidex, tile_sidey, width, nytiles);
  cout << "Rank " << mpiMgr.rank() << " : Device rendering on " << newParticle << " particles" << endl;
  cudaEventRecord(stop,0);

  // Host rendering 
  if (nHostPart > 0)
  {
     cout << "Rank " << mpiMgr.rank() << " : Host rendering on " << nHostPart << " particles" << endl;
     host_funct::render_new(host_part, nHostPart, Pic_host, a_eq_e, grayabsorb);
  }

  cudaEventSynchronize(stop);
  float elapsedTime;
  cudaEventElapsedTime(&elapsedTime, start, stop);
  cout << "Device Rendering Time = " << elapsedTime/1000.0 << endl;
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  cudaThreadSynchronize();
  //cout << cudaGetErrorString(cudaGetLastError()) << endl;
  cu_combine(xres * yres, gv);
  cudaThreadSynchronize();
  //cout << cudaGetErrorString(cudaGetLastError()) << endl;
  tstack_pop("CUDA Rendering");

  // copy back the image
  tstack_push("Data copy");
  error = cudaMemcpy(Pic, gv->d_pic, size_Im, cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) cout << "Device Memcpy error!" << endl; 
  tstack_pop("Data copy");

  if (host_part) cudaFreeHost(host_part);
  return nHostPart+newParticle;
}
