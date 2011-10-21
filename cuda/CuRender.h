#ifndef CURENDER_H
#define CURENDER_H

#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"


using namespace std;

//void cu_draw_chunk(void *pinfo, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb);
void cu_draw_chunk(void *pinfo, COLOUR *Pic, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb);
int filter_chunk(int StartP, int nParticle, int maxRegion,
                 int nFBufInCell, cu_particle_splotch *cu_ps,
                 cu_particle_splotch *cu_ps_filtered, int *End_cu_ps, 
                 int *nFragments2Render);
void render_chunk(int EndP, int nFBufInCell, cu_particle_splotch *cu_ps_filtered,                           void *fragBuf, cu_gpu_vars *gv, bool a_eq_e, float64 grayabsorb,
                  arr2<COLOUR> &pPic, wallTimerSet &times);

#endif
