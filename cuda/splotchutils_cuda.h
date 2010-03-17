#ifndef SPLOTCHUTILS_CUDA_H
#define SPLOTCHUTILS_CUDA_H

#include <vector>
#include "cxxsupport/datatypes.h"
#include "kernel/colour.h"
#include "kernel/transform.h"
#include "splotch/splotchutils.h"
#include "cuda/splotch_cuda.h"

extern "C" void getCuTransformParams(cu_param_transform &para_trans,
                          paramfile &params, double c[3], double l[3], double s[3]);

#endif
