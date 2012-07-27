#include "cuda/CuPolicy.h"

CuPolicy::CuPolicy(int xres, int yres)
  {
    res.first = xres;
    res.second = yres;
    tile_size.first = 12;  // x side dimension of the image tile, in terms of pixels
    tile_size.second = 12; // y side dimension of the image tile, in terms of pixels
    boundary_width = 4;    // width of the boundary around the image tile = max particle radius
    x_num_tiles = xres/tile_size.first;
    y_num_tiles = yres/tile_size.second;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
  //  p_blockSize = deviceProp.maxThreadsPerBlock;
    p_blockSize = 512;
    m_gridSize = deviceProp.maxGridSize[0];
    gmsize = deviceProp.totalGlobalMem;
  }

void CuPolicy::GetTileInfo(int *tile_sidex, int *tile_sidey, int *width, int *nxtiles, int *nytiles)
  {
    *tile_sidex = tile_size.first;
    *tile_sidey = tile_size.second;
    *width = boundary_width;
    *nxtiles = x_num_tiles;
    *nytiles = y_num_tiles;
  }

int CuPolicy::GetNumTiles()
{
   return x_num_tiles*y_num_tiles;
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

int CuPolicy::GetBlockSize() 
  { 
    return p_blockSize; 
  }

size_t CuPolicy::GetImageSize()
{
    size_t size = (res.first)*(res.second)*sizeof(cu_color);
    return size;
}

void CuPolicy::GetDimsBlockGrid(int n, dim3 *dimGrid, dim3 *dimBlock)
  {
    *dimBlock = dim3(p_blockSize);
    int nBlock = (n + p_blockSize - 1)/p_blockSize;
    *dimGrid =dim3(nBlock); 
    if (nBlock > m_gridSize)
      cout << "Error: dim grid = " << nBlock << "too large!" << endl;
  }
