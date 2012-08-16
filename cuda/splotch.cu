
#include <stdlib.h>
#include <math.h>

#include <cuda.h>

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
paramfile &params, const vec3 &campos, const vec3 &lookat, vec3 &sky)
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


int cu_init(int devID, long int nP, int ntiles, cu_gpu_vars* pgv, paramfile &fparams, const vec3 &campos, const vec3 &lookat, vec3 &sky, float b_brightness)
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

  // particle active flag vector
  size = nP * sizeof(int);
  error = cudaMalloc((void**) &pgv->d_active, size);
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: active flag vector allocation error!" << endl;
     return 1;
   }

  // index buffer for C3 particles
  size = nP * sizeof(int);
  error = cudaMalloc((void**) &pgv->d_index, size);
  if (error != cudaSuccess) 
  {
   cout << "Device Memory: index buffer allocation error!" << endl;
   return 1;
  }

  // image + 3 copies
  size = pgv->policy->GetImageSize();
  error = cudaMalloc((void**) &pgv->d_pic, size); 
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: image allocation error!" << endl;
     return 1;
   }
  error = cudaMalloc((void**) &pgv->d_pic1, size); 
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: image 1 allocation error!" << endl;
     return 1;
   }
  error = cudaMalloc((void**) &pgv->d_pic2, size); 
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: image 2 allocation error!" << endl;
     return 1;
   }
  error = cudaMalloc((void**) &pgv->d_pic3, size); 
  if (error != cudaSuccess)
   {
     cout << "Device Malloc: image 3 allocation error!" << endl;
     return 1;
   }

  // tiles
  size = (ntiles+1)*sizeof(int);
  error = cudaMalloc((void**) &pgv->d_tiles, size); 
  if (error != cudaSuccess) 
  {
     cout << "Device Malloc: array tiles allocation error!" << endl;
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

int cu_process (int n, cu_gpu_vars* pgv, int tile_sidex, int tile_sidey, int width, int nxtiles, int nytiles)
  {
  //Get block dim and grid dim from pgv->policy object
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  cudaFuncSetCacheConfig(k_process, cudaFuncCachePreferL1);
  k_process<<<dimGrid,dimBlock>>>(pgv->d_pd, pgv->d_active, n, pgv->colormap_size, pgv->colormap_ptypes, tile_sidex, tile_sidey, width, nxtiles, nytiles);
 
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


void cu_combine(int nP, int nC3, int res, cu_gpu_vars* pgv)
{
  //fetch grid dim and block dim and call device
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(res, &dimGrid, &dimBlock);

  cudaFuncSetCacheConfig(k_add_images, cudaFuncCachePreferL1);
  k_add_images<<<dimGrid,dimBlock>>>(res, pgv->d_pic, pgv->d_pic1, pgv->d_pic2, pgv->d_pic3);

  if (nC3 > 0)
  {
    cudaThreadSynchronize();
    pgv->policy->GetDimsBlockGrid(nC3, &dimGrid, &dimBlock);
    k_addC3<<<dimGrid,dimBlock>>>(nC3, pgv->d_index, pgv->d_pd+nP-nC3, pgv->d_pic);
  }
}


void cu_render1
  (int grid, int block, bool a_eq_e, float grayabsorb, cu_gpu_vars* pgv, int tile_sidex, int tile_sidey, int width, int nytiles)
  {
  //get dims from pgv->policy object first
  dim3 dimGrid = dim3(grid); 
  dim3 dimBlock = dim3(block);
  size_t SharedMem = (tile_sidex+2*width)*(tile_sidey+2*width)*sizeof(cu_color);

  cudaFuncSetCacheConfig(k_render1, cudaFuncCachePreferShared);
  k_render1<<<dimGrid, dimBlock, SharedMem>>>(pgv->d_pd, pgv->d_active, pgv->d_tiles, pgv->d_pic,
  pgv->d_pic1, pgv->d_pic2, pgv->d_pic3, tile_sidex, tile_sidey, width, nytiles);
  }

void cu_indexC3(int nP, int nC3, cu_gpu_vars* pgv)
  {
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(nC3, &dimGrid, &dimBlock);
 
  k_renderC3<<<dimGrid, dimBlock>>>(nC3, pgv->d_pd+nP-nC3, pgv->d_index);
  }


void cu_end(cu_gpu_vars* pgv)
  {
  CLEAR_MEM((pgv->d_pd));
  CLEAR_MEM((pgv->d_active));
  CLEAR_MEM((pgv->d_index));
  CLEAR_MEM((pgv->d_pic));
  CLEAR_MEM((pgv->d_pic1));
  CLEAR_MEM((pgv->d_pic2));
  CLEAR_MEM((pgv->d_pic3));
  CLEAR_MEM((pgv->d_tiles));

  delete pgv->policy;
  cudaThreadExit();
  }

long int cu_get_chunk_particle_count(CuPolicy* policy, size_t psize, int ntiles, float pfactor)
  {
   size_t gMemSize = policy->GetGMemSize();
   size_t ImSize = policy->GetImageSize();
   size_t tiles = ntiles*sizeof(int);

   size_t spareMem = 20*(1<<20);
   long int arrayParticleSize = gMemSize - 4*ImSize - tiles - spareMem;
   long int len = (long int) (arrayParticleSize/((psize+2*sizeof(int))*pfactor)); 
   long int maxlen = policy->GetMaxGridSize() * policy->GetBlockSize();
   if (len > maxlen) len = maxlen;
   return len;
  }


