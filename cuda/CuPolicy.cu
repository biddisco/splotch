

#include "cuda/CuPolicy.h"

CuPolicy::CuPolicy(paramfile &Param)
  {
    res.first = Param.find<int>("xres",800);
    res.second = Param.find<int>("yres",res.first);

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    p_blockSize = deviceProp.maxThreadsPerBlock;
    m_gridSize = deviceProp.maxGridSize[0];
    gmsize = deviceProp.totalGlobalMem;

    size_t fbsize_def = (size_t) p_blockSize; //min(res.first/20, p_blockSize);
    fbsize_def *= m_gridSize*sizeof(cu_color);
    if ((8*fbsize_def) > gmsize) fbsize_def = gmsize/8;
    fbsize = Param.find<int>("fragment_buffer_size", fbsize_def);
    pix_blockSize = fbsize/(m_gridSize*sizeof(cu_color));
  }

pair<int,int> CuPolicy::GetResolution()
  {
    return res;
  }

size_t CuPolicy::GetFBufSize() // return dimension in terms of bytes
  {
     return fbsize; 
  }

size_t CuPolicy::GetIndexSize() // return dimension in terms of bytes
  {
     int npixels = m_gridSize*pix_blockSize; // (int) fbsize/sizeof(cu_color);
     size_t size = npixels*sizeof(int);
     return size; 
  }

size_t CuPolicy::GetGMemSize() // return dimension in terms of bytes
  { 
   // int MB = 1<<20;
   // int size = gmsize/MB;
    return gmsize; 
  }

int CuPolicy::GetMaxGridSize() 
  { 
    return m_gridSize; 
  }

int CuPolicy::GetMaxBlockSize()
  { 
    return pix_blockSize; 
  }

size_t CuPolicy::GetImageSize()
{
   // int MB = 1<<20;
    size_t size = (res.first)*(res.second)*sizeof(cu_color);
    return size;
}

void CuPolicy::GetDimsBlockGrid(int n, dim3 *dimGrid, dim3 *dimBlock)
  {
    *dimBlock = dim3(p_blockSize);
    int nBlock = (n + p_blockSize - 1)/p_blockSize;
    *dimGrid =dim3(nBlock); 
  }
