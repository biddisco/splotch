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

#ifdef __CUDACC__
#include <cuda.h>
#else
struct dim3;
#endif

using namespace std;

class CuPolicy
  {
  private:
    int m_blockSize, maxregion, fbsize;
    pair <int,int> res;
    size_t gmsize;
  public:
    CuPolicy(paramfile &Param);

    pair <int,int> GetResolution();
    int GetMaxRegion();
    int GetFBufSize();
    int GetIndexSize();
    int GetGMemSize();
    int GetImageSize();
    void GetDimsBlockGrid(int n, dim3 *dimGrid, dim3 *dimBlock);
  };

#endif //CUPOLICY_H
