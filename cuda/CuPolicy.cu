

#include "cuda/CuPolicy.h"

CuPolicy::CuPolicy(paramfile &Param)
  {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    m_blockSize = deviceProp.maxThreadsPerBlock;
    gmsize = deviceProp.totalGlobalMem;
    maxregion = Param.find<int>("max_region", 1024);
    fbsize = Param.find<int>("fragment_buffer_size", 1024);
//    fbsize = GetGMemSize()/4;
    res.first = Param.find<int>("xres",800),
    res.second = Param.find<int>("yres",res.first);
  }

pair<int,int> CuPolicy::GetResolution()
  {
    return res;
  }

int CuPolicy::GetMaxRegion()
  { return maxregion; }

int CuPolicy::GetFBufSize() // return dimension in terms of Megabytes
  {
     return fbsize; 
  }

int CuPolicy::GetIndexSize() // return dimension in terms of Megabytes
  {
     int npixels = (int) fbsize/sizeof(cu_color);
     int size = npixels*sizeof(int);
     return size; 
  }

int CuPolicy::GetGMemSize() // return dimension in terms of Megabytes
  { 
    int MB = 1<<20;
    int size = gmsize/MB;
    return size; 
  }

int CuPolicy::GetImageSize()
{
    int MB = 1<<20;
    int size = (res.first)*(res.second)*sizeof(cu_color);
    return (int) size/MB + 1;
}

void CuPolicy::GetDimsBlockGrid(int n, dim3 *dimGrid, dim3 *dimBlock)
  {
    *dimBlock = dim3(m_blockSize);
    int nBlock = (n+m_blockSize-1)/m_blockSize;
    *dimGrid =dim3(nBlock); 
  }
