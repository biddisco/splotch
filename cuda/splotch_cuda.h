/*
Copyright things go here.
*/

#ifndef SPLOTCH_CUDA_H
#define SPLOTCH_CUDA_H

#ifdef VS
#define THREADFUNC DWORD WINAPI
#else
#define THREADFUNC static void*
#endif

#define SPLOTCH_CLASSIC

#include <cstring>
#include "cxxsupport/paramfile.h"
#include "kernel/colour.h"
#include "splotch/splotchutils.h"
class CuPolicy;

//data structs for using on device
//'d_' means device
//typedef particle_sim cu_particle_sim;

struct cu_color
  {
  float r,g,b;
  };

struct cu_particle_sim
  {
    cu_color e;
    float x,y,z,r,I;
    unsigned short type;
    bool active;
  };

#define MAX_P_TYPE 8//('XXXX','TEMP','U','RHO','MACH','DTEG','DISS','VEL')
                                        //in mid of developing only
#define MAX_EXP -20.0

struct cu_param
  {
  float p[12];
  bool  projection;
  int   xres, yres;
  float fovfct, dist, xfac;
  float minrad_pix;
  int ptypes;
  float zmaxval, zminval;
  bool col_vector[MAX_P_TYPE];
  float brightness[MAX_P_TYPE];
  float bfak, h2sigma, sigma0, rfac;
  };


struct cu_color_map_entry
  {
  float val;
  cu_color color;
  };

struct cu_colormap_info
  {
  cu_color_map_entry *map;
  int mapSize;
  int *ptype_points;
  int ptypes;
  };


struct cu_fragment_AeqE
  {
  float aR, aG, aB;
  };

struct cu_fragment_AneqE
  {
  float aR, aG, aB;
  float qR, qG, qB;
  };


struct cu_gpu_vars //variables used by each gpu
  {
  CuPolicy            *policy;
  cu_particle_sim     *d_pd;             //device_particle_data
  unsigned long	      *d_posInFragBuf;
  char		      *d_active; // 0=non-active particle, 1=active particle, 2=active big particle
  void                *d_fbuf;
  int                 *d_pixel;
  cu_color            *d_pic;
  int                 colormap_size;
  int                 colormap_ptypes;
  int 		      maxsize;	//max size of particles to be processed on the device
  };

//functions

int cu_init(int devID, long int nP, cu_gpu_vars* pgv, paramfile &fparams, vec3 &campos, vec3 &lookat, vec3 &sky, float b_brightness);
int cu_copy_particles_to_device(cu_particle_sim* h_pd, unsigned int n, cu_gpu_vars* pgv);
int cu_transform (int n, cu_gpu_vars* pgv);
int cu_allocateFragmentBuffer(long n,cu_gpu_vars* pgv);
void cu_init_colormap(cu_colormap_info info, cu_gpu_vars* pgv);
//void cu_colorize(int n, cu_gpu_vars* pgv);
void cu_render1 (int grid, int block, int nP, int End_cu_ps, 
   unsigned long FragRendered, bool a_eq_e, float grayabsorb, cu_gpu_vars* pgv);
void cu_endChunk (cu_gpu_vars* pgv);
void cu_endThread (cu_gpu_vars* pgv);
long int cu_get_chunk_particle_count(CuPolicy* policy, size_t psize, float pfactor);
void getCuTransformParams(cu_param &para_trans,
      paramfile &params, vec3 &campos, vec3 &lookat, vec3 &sky);
//void cu_clear(int n, cu_gpu_vars* pgv);
void cu_update_image(int n, bool a_eq_e,cu_gpu_vars* pgv);

#endif
