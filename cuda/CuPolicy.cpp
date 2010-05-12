

#include "CuPolicy.h"

CuPolicy::CuPolicy(paramfile &Param)
  {
    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, 0));
    m_blockSize = deviceProp.maxThreadsPerBlock;
    maxregion = Param.find<int>("max_region", 1024);
    fbsize = (Param.find<int>("fragment_buffer_size", 100))<<20;
  }

//Get size of device particle data
int CuPolicy::GetSizeDPD(unsigned int n)
  { return n* sizeof(cu_particle_sim); }

int CuPolicy::GetMaxRegion()
  { return maxregion; }

int CuPolicy::GetFBufSize()
  { return fbsize; }

void CuPolicy::GetDimsBlockGrid(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
  {
    *dimBlock = dim3(m_blockSize);
    unsigned int nBlock = (n+m_blockSize-1)/m_blockSize;
    *dimGrid =dim3(nBlock); 
  }
