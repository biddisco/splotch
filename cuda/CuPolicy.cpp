/*
Copyright things go here.
*/

#include "CuPolicy.h"
#include <cmath>

CuPolicy::CuPolicy()
  : m_pParam(0), m_blockSize (256) {}

CuPolicy::CuPolicy(paramfile *pParam)
  {
  m_pParam =pParam;
  m_blockSize =pParam->find<int>("block_size", 256);
  }

//Get size of device particle data
int CuPolicy::GetSizeDPD(unsigned int n)
  { return n* sizeof(cu_particle_sim); }

int CuPolicy::GetMaxRegion()
  { return m_pParam->find<int>("max_region", 1024); }

int CuPolicy::GetFBufSize()
  {
  //just for test on my own pc now
  int size=m_pParam->find<int>("fragment_buffer_size", 100);
  unsigned int M = 1<<20;
  return int (size*M);
  }

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
