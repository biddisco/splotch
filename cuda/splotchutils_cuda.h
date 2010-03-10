#ifndef SPLOTCHUTILS_CUDA_H
#define SPLOTCHUTILS_CUDA_H

#ifdef CUDA

#include <vector>
#include "cxxsupport/datatypes.h"
#include "kernel/colour.h"
#include "kernel/transform.h"
#include "splotch/splotchutils.h"
#include "cuda/splotch_cuda.h"

#ifdef HOST_THREAD_RENDER
struct  param_render_thread
  {
  cu_particle_splotch     *p;
  int start, end;
  bool a_eq_e;
  double grayabsorb;
  cu_color pic[][800];
  };
#endif //ifdef HOST_THREAD_RENDER

void render_as_thread1 (const vector<particle_sim> &p, arr2<COLOUR> &pic,
                        bool a_eq_e,double grayabsorb);

#ifdef HOST_THREAD_RENDER
DWORD WINAPI render_thread (param_render_thread *param);
#endif //ifdef HOST_THREAD_RENDER

extern "C" void getCuTransformParams(cu_param_transform &para_trans,
                          paramfile &params, double c[3], double l[3], double s[3]);

#endif  //CUDA

#endif
