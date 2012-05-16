
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

  float dist= (campos-lookat).Length();
  float xfac=1./(fovfct*dist);
  if (!projection)
    cout << " Field of fiew: " << 1./xfac*2. << endl;

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

  float h2sigma = 0.5*pow(Pi,-1./6.);
  para_trans.h2sigma=h2sigma;
#ifdef SPLOTCH_CLASSIC
  float powtmp = pow(Pi,1./3.);
  float sigma0 = powtmp/sqrt(2*Pi);
  para_trans.sigma0=sigma0;
  para_trans.bfak=1./(2.*sqrt(Pi)*powtmp);
  float rfac=1.5*h2sigma/(sqrt(2.)*sigma0);
#else
  float rfac=1.;
#endif
  para_trans.rfac=rfac;
  }


int cu_init(int devID, long int nP, cu_gpu_vars* pgv, paramfile &fparams, vec3 &campos, vec3 &lookat, vec3 &sky, float b_brightness)
  {
  cudaError_t error;
  cudaSetDevice (devID); // initialize cuda runtime
 
  // particle vector  
  size_t size = nP * sizeof(cu_particle_sim);
  error = cudaMalloc((void**) &pgv->d_pd, size);
  if (error != cudaSuccess) 
  {
   cout << "Device Memory: particle data allocation error!" << endl;
   return 1;
  }
  
  // particle size/position vector
  size = nP * sizeof(unsigned long);
  error = cudaMalloc((void**) &pgv->d_posInFragBuf, size);
  if (error != cudaSuccess)
   {
     cout << "Device Malloc " << size << " bytes: size particles vector allocation error!" << endl;
     return 1;
   }
  // particle active flag vector
  size = nP * sizeof(char);
  error = cudaMalloc((void**) &pgv->d_active, size);
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: active flag vector allocation error!" << endl;
     return 1;
   }

  // image
  size = pgv->policy->GetImageSize();
  error = cudaMalloc((void**) &pgv->d_pic, size); 
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: image allocation error!" << endl;
     return 1;
   }

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
//  tparams.rfac=1.5;

  //dump parameters to device
  error = cudaMemcpyToSymbol(dparams, &tparams, sizeof(cu_param));
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: parameters allocation error!" << endl;
     return 1;
   }
  return 0;
  }


int cu_copy_particles_to_device(cu_particle_sim* h_pd, unsigned int n, cu_gpu_vars* pgv)
  {
  cudaError_t error;
  size_t size = n*sizeof(cu_particle_sim);
  error = cudaMemcpy(pgv->d_pd, h_pd, size, cudaMemcpyHostToDevice);
  if (error != cudaSuccess)
  {
   cout << "Device Memory: particle data copy error!" << endl;
   return 1;
  }
  return 0;
  }

int cu_allocateFragmentBuffer(long n, cu_gpu_vars* pgv)
{
  cudaError_t error;
  // fragment buffer
  size_t size = n*sizeof(cu_fragment_AeqE);
  error = cudaMalloc((void**) &pgv->d_fbuf, size); 
  if (error != cudaSuccess) 
   {
     cout << "Device Malloc: fragment buffer allocation error!" << endl;
     return 1;
   }
  // index fragment buffer (index pixels)
  size = n*sizeof(int);
  error = cudaMalloc((void**) &pgv->d_pixel, size); 
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: pixel indexes allocation error!" << endl;
     return 1;
   }
  return 0;
}


int cu_transform (int n, cu_gpu_vars* pgv)
  {
  //Get block dim and grid dim from pgv->policy object
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);
  int maxPartSize = pgv->policy->GetMaxPartSize();

  //call device transformation
  k_transform<<<dimGrid,dimBlock>>>(pgv->d_pd, pgv->d_posInFragBuf, pgv->d_active, n, maxPartSize);
  cudaThreadSynchronize();
 
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
  k_colorize<<<dimGrid,dimBlock>>>(pgv->d_pd, pgv->colormap_size, pgv->colormap_ptypes, n);
  cudaThreadSynchronize();
  }

void cu_clear(int n, cu_gpu_vars* pgv)
  {

  //fetch grid dim and block dim and call device
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  k_clear<<<dimGrid,dimBlock>>>(n, pgv->d_pic);
  }
*/

void cu_update_image(int n, bool a_eq_e, cu_gpu_vars* pgv)
{
  //fetch grid dim and block dim and call device
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  k_combine<<<dimGrid,dimBlock>>>(n, a_eq_e, pgv->d_pic, pgv->d_pixel, pgv->d_fbuf);
}


void cu_render1
  (int grid, int block, int nP, int End_cu_ps, unsigned long FragRendered,
   bool a_eq_e, float grayabsorb, cu_gpu_vars* pgv)
  {
  //get dims from pgv->policy object first
  dim3 dimGrid = dim3(grid); 
  dim3 dimBlock = dim3(block);
 // size_t sizeSharedMem = block*sizeof(cu_particle_sim) ;

  k_render1<<<dimGrid, dimBlock>>>(pgv->d_posInFragBuf+End_cu_ps, nP, FragRendered, pgv->d_pd+End_cu_ps, pgv->d_fbuf, pgv->d_pixel, a_eq_e, grayabsorb);
  }

void cu_endChunk(cu_gpu_vars* pgv)
  {
  CLEAR_MEM((pgv->d_fbuf));
  CLEAR_MEM((pgv->d_pixel));
  }

void cu_endThread(cu_gpu_vars* pgv)
  {
  CLEAR_MEM((pgv->d_pd));
  CLEAR_MEM((pgv->d_posInFragBuf));
  CLEAR_MEM((pgv->d_active));
  CLEAR_MEM((pgv->d_pic));

  delete pgv->policy;
  cudaThreadExit();
  }

long int cu_get_chunk_particle_count(CuPolicy* policy, size_t psize, float pfactor)
  {
   size_t gMemSize = policy->GetGMemSize();
   size_t fBufSize = policy->GetFBufSize();
   size_t ImSize = policy->GetImageSize();
   size_t IndexSize = policy->GetIndexSize();
  // if (gMemSize <= fBufSize) return 0;

   size_t spareMem = 20*(1<<20);
   long int arrayParticleSize = gMemSize - 2*fBufSize - 2*IndexSize - ImSize - spareMem;
   long int len = (long int) (arrayParticleSize/(psize*pfactor)); 

   return len;
  }


