/*
Copyright things go here.
*/

#include "CuPolicy.h"

CuPolicy::CuPolicy(paramfile &Param)
  {
  m_blockSize =Param.find<int>("block_size", 256);
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

void CuPolicy::GetDims1(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
  {
  *dimBlock = dim3(m_blockSize);
  unsigned int nBlock = (n+m_blockSize-1)/m_blockSize;
  *dimGrid =dim3(nBlock);
  }

void CuPolicy::GetDimsPostProcess(int xres, int yres,
  dim3 *dimGrid, dim3 *dimBlock)
  { GetDims1 (xres*yres, dimGrid, dimBlock); }

void CuPolicy::GetDimsCombine(unsigned int minx, unsigned int miny,
  unsigned int maxx, unsigned int maxy,dim3 *dimGrid, dim3 *dimBlock)
  { GetDims1 ((maxx-minx)*(maxy-miny), dimGrid, dimBlock); }

void CuPolicy::GetDimsRange(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
  { GetDims1(n, dimGrid, dimBlock); }

void CuPolicy::GetDimsColorize(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
  { GetDims1(n, dimGrid, dimBlock); }

void CuPolicy::GetDimsRender(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
  { GetDims1(n, dimGrid, dimBlock); }
