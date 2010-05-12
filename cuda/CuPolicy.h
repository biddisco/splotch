/*
Copyright things go here.
*/
/*
CuPolicy is the class that knows the overall state of cuda application.
All 'magic numbers' are out of this class.
*/
#ifndef CUPOLICY_H
#define CUPOLICY_H

#include "cxxsupport/paramfile.h"
#include "cuda/splotch_cuda.h"
#include <cuda.h>
#include <cutil_inline.h>

class CuPolicy
  {
  private:
    unsigned int m_blockSize, maxregion, fbsize;

  public:
    CuPolicy(paramfile &Param);

    int GetMaxRegion();
    int GetFBufSize();
    int GetSizeDPD(unsigned int n);
    void GetDimsBlockGrid(unsigned int n, dim3 *dimGrid, dim3 *dimBlock);
    unsigned int GetFBufSize(bool a_eq_e);
  };

#endif //CUPOLICY_H
