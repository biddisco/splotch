#ifndef CURENDER_H
#define CURENDER_H

#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"


using namespace std;

//void cu_draw_chunk(void *pinfo, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb);
void cu_draw_chunk(void *pinfo, COLOUR *Pic, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb);

#endif
