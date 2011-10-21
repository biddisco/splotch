/*
Try accelerating splotch with CUDA. July 2009.
Copyright things go here.
*/

#include <stdlib.h>
#include <math.h>

#include <cuda.h>
#include <thrust/version.h>

#include "cxxsupport/lsconstants.h"
#include "cxxsupport/string_utils.h"
#include "splotch/splotchutils.h"
#include "kernel/transform.h"

#include "cuda/splotch_cuda.h"
#include "cuda/CuPolicy.h"

#include "cuda/splotch_kernel.cu"
#include "cuda/CuRender.cu"

using namespace std;

template<typename T> T findParamWithoutChange
  (paramfile *param, std::string &key, T &deflt)
  {
  return param->param_present(key) ? param->find<T>(key) : deflt;
  }

#define CLEAR_MEM(p) if(p) {cudaFree(p); p=0;}


void getCuTransformParams(cu_param &para_trans,
paramfile &params, vec3 &campos, vec3 &lookat, vec3 &sky)
  {
  int xres = params.find<int>("xres",800),
      yres = params.find<int>("yres",xres);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);
  float64 xfac=0.0, dist=0.0;

  sky.Normalize();
  vec3 zaxis = (lookat-campos).Norm();
  vec3 xaxis = crossprod (sky,zaxis).Norm();
  vec3 yaxis = crossprod (zaxis,xaxis);
  TRANSFORM trans;
  trans.Make_General_Transform
        (TRANSMAT(xaxis.x,xaxis.y,xaxis.z,
                  yaxis.x,yaxis.y,yaxis.z,
                  zaxis.x,zaxis.y,zaxis.z,
                  0,0,0));
  trans.Invert();
  TRANSFORM trans2;
  trans2.Make_Translation_Transform(-campos);
  trans2.Add_Transform(trans);
  trans=trans2;
  bool projection = params.find<bool>("projection",true);

  if (!projection)
    {
    float64 dist= (campos-lookat).Length();
    float64 xfac=1./(fovfct*dist);
    cout << " Field of fiew: " << 1./xfac*2. << endl;
    }

  float minrad_pix = params.find<float>("minrad_pix",1.);

  //retrieve the parameters for transformation
  for (int i=0; i<12; i++)
    para_trans.p[i] =trans.Matrix().p[i];
  para_trans.projection=projection;
  para_trans.xres=xres;
  para_trans.yres=yres;
  para_trans.fovfct=fovfct;
  para_trans.dist=dist;
  para_trans.xfac=xfac;
  para_trans.minrad_pix=minrad_pix;
  }


int cu_init(int devID, int nP, cu_gpu_vars* pgv, paramfile &fparams, vec3 &campos, vec3 &lookat, vec3 &sky, float brightness)
  {
  cudaError_t error;
  cudaSetDevice (devID); // initialize cuda runtime
  
  //now prepare memory for d_particle_splotch.
  //one more for dums
  size_t s = nP* sizeof(cu_particle_splotch);
  error = cudaMalloc((void**) &pgv->d_ps_render, s+sizeof(cu_particle_splotch));
  if (error != cudaSuccess) return 1;

  // fragment buffer
  size_t size = pgv->policy->GetFBufSize() <<20;
  error = cudaMalloc((void**) &pgv->d_fbuf, size); 
  if (error != cudaSuccess) return 1;
  // index fragment buffer (index pixels)
  size = pgv->policy->GetIndexSize() <<20;
  error = cudaMalloc((void**) &pgv->d_pixel, size); 
  if (error != cudaSuccess) return 1;
  // image
  size = pgv->policy->GetImageSize() <<20;
  error = cudaMalloc((void**) &pgv->d_pic, size); 
  if (error != cudaSuccess) return 1;

  //retrieve parameters
  cu_param tparams;
  getCuTransformParams(tparams,fparams,campos,lookat,sky);

  tparams.zmaxval   = fparams.find<float>("zmax",1.e23);
  tparams.zminval   = fparams.find<float>("zmin",0.0);
  tparams.ptypes    = fparams.find<int>("ptypes",1);

  for(int itype=0; itype<tparams.ptypes; itype++)
    {
    tparams.brightness[itype] = fparams.find<double>("brightness"+dataToString(itype),1.);
    tparams.brightness[itype] *= b_brightness;
    tparams.col_vector[itype] = fparams.find<bool>("color_is_vector"+dataToString(itype),false);
    }
  tparams.rfac=1.5;

  //dump parameters to device
  error = cudaMemcpyToSymbol(dparams, &tparams, sizeof(cu_param));
  if (error != cudaSuccess) return 1;
  return 0;
  }


int cu_copy_particles_to_device(cu_particle_sim* h_pd, unsigned int n, cu_gpu_vars* pgv)
  {
  cudaError_t error;
  size_t size = (n+1)*sizeof(cu_particle_sim);
  //one more space allocated for the dumb
  error = cudaMalloc((void**) &pgv->d_pd, size);
  if (error != cudaSuccess) 
  {
   cout << "Device Memory: allocation particle data error!" << endl;
   return 1;
  }

  //copy particle data to device
  error = cudaMemcpy(pgv->d_pd, h_pd, size, cudaMemcpyHostToDevice);
  if (error != cudaSuccess)
  {
   cout << "Device Memory: copy particle data error!" << endl;
   return 1;
  }
  return 0;
  }


int cu_transform (unsigned int n, cu_particle_splotch *h_ps, cu_gpu_vars* pgv)
  {
  cudaError_t error;
  //Get block dim and grid dim from pgv->policy object
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  //call device transformation
  k_transform<<<dimGrid,dimBlock>>>(pgv->d_pd, pgv->d_ps_render, n);

  //copy the result out
  size_t size = n* sizeof(cu_particle_splotch);
  error = cudaMemcpy(h_ps, pgv->d_ps_render, size, cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) 
   {
     cout << "Device Memcpy: particle data copy error!" << endl;
     return 1;
   } 
  return 0;
  }

void cu_init_colormap(cu_colormap_info h_info, cu_gpu_vars* pgv)
  {
  //allocate memories for colormap and ptype_points and dump host data into it
  size_t size =sizeof(cu_color_map_entry)*h_info.mapSize;
  cudaMemcpyToSymbol(dmap, h_info.map, size);
  //type
  size =sizeof(int)*h_info.ptypes;
  cudaMemcpyToSymbol(ptype_points, h_info.ptype_points, size);

  //set fields of global variable pgv->d_colormap_info
  pgv->colormap_size   = h_info.mapSize;
  pgv->colormap_ptypes = h_info.ptypes;
  }

/*
void cu_colorize(int n, cu_gpu_vars* pgv)
  {

  //fetch grid dim and block dim and call device
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord( start, 0);

  k_colorize<<<dimGrid,dimBlock>>>(n, pgv->d_ps_render, pgv->colormap_size, pgv->colormap_ptypes);

  cudaEventRecord( stop, 0);
  cudaEventSynchronize(stop);
 // float elapsedTime;
 // cutilSafeCall( cudaEventElapsedTime(&elapsedTime,start,stop));
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  //particle_splotch memory on device will be freed in cu_end
  }
*/

int cu_copy_particles_to_render(cu_particle_splotch *p,
  int n, cu_gpu_vars* pgv)
  {
  cudaError_t error;
  //copy filtered particles into device
  size_t size = n *sizeof(cu_particle_splotch);
  error = cudaMemcpy(pgv->d_ps_render, p,size, cudaMemcpyHostToDevice);
  if (error != cudaSuccess) return 1;
  return 0;
  }

void cu_clear(int n, cu_gpu_vars* pgv)
  {

  //fetch grid dim and block dim and call device
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  k_clear<<<dimGrid,dimBlock>>>(n, pgv->d_pic);
  }

void cu_update_image(int n, bool a_eq_e, cu_gpu_vars* pgv)
{
  //fetch grid dim and block dim and call device
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  k_combine<<<dimGrid,dimBlock>>>(n, a_eq_e, pgv->d_pic, pgv->d_pixel, pgv->d_fbuf);

}


void cu_render1
  (int nP, bool a_eq_e, float grayabsorb, cu_gpu_vars* pgv)
  {
//  cudaEvent_t start,stop;
//  cudaEventCreate(&start);
//  cudaEventCreate(&stop);
//  cudaEventRecord( start, 0);
 
  //get dims from pgv->policy object first
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(nP, &dimGrid, &dimBlock);
  int yres = (pgv->policy->GetResolution()).second;

  //call device
  k_render1<<<dimGrid, dimBlock>>>(pgv->d_ps_render, nP,
    pgv->d_fbuf, pgv->d_pixel, a_eq_e, grayabsorb, pgv->colormap_size, pgv->colormap_ptypes, yres);

//  cudaEventRecord( stop, 0);
//  cudaEventSynchronize(stop);
//  float elapsedTime;
//  cutilSafeCall( cudaEventElapsedTime(&elapsedTime,start,stop));
//  cudaEventDestroy(start);
//  cudaEventDestroy(stop);
  }

/*
void cu_get_fbuf
  (void *h_fbuf, bool a_eq_e, unsigned long n, cu_gpu_vars* pgv)
  {
  size_t size;
  if (a_eq_e)
    size =n* sizeof(cu_fragment_AeqE);
  else
    size =n* sizeof(cu_fragment_AneqE);

  cudaMemcpy(h_fbuf, pgv->d_fbuf,size,cudaMemcpyDeviceToHost);
  }
*/

void cu_end(cu_gpu_vars* pgv)
  {
  CLEAR_MEM((pgv->d_ps_render));
  CLEAR_MEM((pgv->d_fbuf));
  CLEAR_MEM((pgv->d_pixel));
  CLEAR_MEM((pgv->d_pic));

  delete pgv->policy;
  cudaThreadExit();
  }

int cu_get_chunk_particle_count(CuPolicy* policy, size_t psize, float pfactor)
  {
   int gMemSize = policy->GetGMemSize();
   int fBufSize = policy->GetFBufSize();
   int ImSize = policy->GetImageSize();
   int IndexSize = policy->GetIndexSize();
   if (gMemSize <= fBufSize) return 0;

   int spareMem = 20;
   int arrayParticleSize = gMemSize - 3*fBufSize -ImSize - spareMem;
   int len = (int) (arrayParticleSize/(psize*pfactor))*(1<<20); 

   return len;
  }


