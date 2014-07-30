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

#ifndef SPLOTCH_CUDA_H
#define SPLOTCH_CUDA_H

#define SPLOTCH_CLASSIC

#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>

#include "cxxsupport/paramfile.h"
#include "cxxsupport/lsconstants.h"
#include "cxxsupport/string_utils.h"
#include "kernel/colour.h"
#include "kernel/transform.h"
#include "splotch/splotchutils.h"

#include "cuda/cuda_policy.h"

// Use atomics for image update instead of layered update in k_renderC2
#define CUDA_USE_ATOMICS

//  Data structs for using on device
// 'd_' means device (for device pointers)

// RGB triplet to represent a color
struct cu_color
{
  float               r,g,b;
};

// Struct to represent a splotch particle on the device
struct cu_particle_sim
{
  cu_color            e;
  float               x,y,z,r,I;
  unsigned short      type;
  bool                active;
};

// Max N types as available in GADGET format
// TODO: Allow as many types as requested
#define MAX_P_TYPE 8//('XXXX','TEMP','U','RHO','MACH','DTEG','DISS','VEL')

struct cu_param
{
  float               p[12];
  bool                projection;
  int                 xres, yres;
  float               fovfct, dist, xfac;
  float               minrad_pix;
  int                 ptypes;
  float               zmaxval, zminval;
  bool                col_vector[MAX_P_TYPE];
  float               brightness[MAX_P_TYPE];
  bool                log_int[MAX_P_TYPE];
  bool                log_col[MAX_P_TYPE];
  bool                asinh_col[MAX_P_TYPE];
  float               inorm_mins[MAX_P_TYPE];
  float               inorm_maxs[MAX_P_TYPE];
  float               cnorm_mins[MAX_P_TYPE];
  float               cnorm_maxs[MAX_P_TYPE];
  bool                do_logs;
  float               bfak, h2sigma, sigma0, rfac;
};

// Map element for device color map
struct cu_color_map_entry
{
  float               val;
  cu_color            color;
};

// Info about colormap to pass to device
struct cu_colormap_info
{
  cu_color_map_entry* map;
  int                mapSize;
  int*               ptype_points;
  int                ptypes;
};

// Variables used by each gpu
struct cu_gpu_vars 
{
  CuPolicy            *policy;
  cu_particle_sim     *d_pd;             // Device ptr to device particle data
  int		              *d_active;         // -1=non-active particle, -2=active big particle, n=number of tile, n+1=C3 particles 
  int                 *d_index;
  cu_color            *d_pic;
  cu_color            *d_pic1;
  cu_color            *d_pic2;
  cu_color            *d_pic3;
  int	      	        *d_tiles;	         // Number of particles per tile
  int	      	        *d_tileID;         // Tile ID
  int                 colormap_size;
  int                 colormap_ptypes;
};


int       cu_init(int devID, long int nP, int ntiles, cu_gpu_vars* pgv);
void      cu_get_trans_params(cu_param &para_trans, paramfile &params, const vec3 &campos, const vec3 &lookat, vec3 &sky);
int       cu_copy_particles_to_device(cu_particle_sim* h_pd, unsigned int n, cu_gpu_vars* pgv);
int       cu_process(int n, cu_gpu_vars* pgv, int tile_sidex, int tile_sidey, int width, int nxtiles, int nytiles);
int       cu_range(int nP, cu_gpu_vars* pgv);
void      cu_init_colormap(cu_colormap_info info, cu_gpu_vars* pgv);
void      cu_renderC2(int nP, int grid, int block, bool a_eq_e, float grayabsorb, cu_gpu_vars* pgv, int tile_sidex, int tile_sidey, int width, int nxtiles);
void      cu_indexC3(int nP, int nC3, cu_gpu_vars* pgv);
void      cu_renderC3(int nP, int nC3, int res, cu_gpu_vars* pgv);
void      cu_end(cu_gpu_vars* pgv);
long int  cu_get_chunk_particle_count(cu_gpu_vars* pgv, int nTasksDev, size_t psize, int ntiles, float pfactor);

#ifdef SPLOTCH_PARAVIEW
int cu_init_params(cu_gpu_vars* pgv, paramfile &fparams, const vec3 &campos, const vec3 &lookat, vec3 &sky, float b_brightness, bool& doLogs);
int cu_copy_particles_from_gpubuffer(void *gpubuffer, unsigned int n, cu_gpu_vars* pgv);
long int cu_paraview_get_chunk_particle_count(cu_gpu_vars* pgv, int nTasksDev, size_t psize, int ntiles, float pfactor, int nP);
#endif

#ifndef CUDA_USE_ATOMICS
  void      cu_add_images(int res, cu_gpu_vars* pgv);
#endif



#endif