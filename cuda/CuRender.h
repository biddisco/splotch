#ifndef CURENDER_H
#define CURENDER_H

#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"


using namespace std;

void cu_draw_chunk(wallTimerSet *times, int mydevID, cu_particle_sim *d_particle_data, int nParticle, COLOUR *Pic, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, float b_brightness, vector<COLOURMAP> &amap, paramfile &g_params);

void host_particle_colorize(paramfile &params, cu_particle_sim *p, int npart,
  vector<COLOURMAP> &amap, float b_brightness);
void host_render_new (cu_particle_sim *p, int npart, arr2<COLOUR> &pic, bool a_eq_e, float32 grayabsorb);
#endif
