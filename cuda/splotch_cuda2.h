#ifndef SPLOTCH_CUDA2_H
#define SPLOTCH_CUDA2_H

#include "cuda/splotch_cuda.h"
#include "cxxsupport/walltimer.h"

// struct containing pthread task info
/*struct thread_info
  {
  int devID;                  //index of the device selected
  int startP, endP;           //start and end particles to handle
  long npart_all;             //total number of particles
  arr2<COLOUR> *pPic;         //output image computed 
  wallTimerSet times;
  };
*/

extern wallTimerSet cuWallTimers;

int check_device(int rank);
void print_device_info(int rank, int dev);

void cuda_rendering(int mydevID, arr2<COLOUR> &pic, std::vector<particle_sim> &particle, const vec3 &campos, const vec3 &lookat, vec3 &sky, std::vector<COLOURMAP> &amap, float b_brightness, paramfile &g_params);
//void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread);
//void cu_draw_chunk(void *pinfo, cu_gpu_vars* gv);
void setup_colormap(int ptypes, std::vector<COLOURMAP> &amap, cu_gpu_vars* gv);

void GPUReport(wallTimerSet &cuTimers);
void cuda_timeReport();

#endif
