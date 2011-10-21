#ifndef SPLOTCH_CUDA2_H
#define SPLOTCH_CUDA2_H

#include "cuda/splotch_cuda.h"
#include "splotch/splotch_host.h"
#include "cxxsupport/walltimer.h"

THREADFUNC host_thread_func(void *pinfo);
THREADFUNC cu_thread_func(void *pinfo);


// struct containing pthread task info
struct thread_info
  {
  int devID;                  //index of the device selected
  int startP, endP;           //start and end particles to handle
  long npart_all;             //total number of particles
  arr2<COLOUR> *pPic;         //output image computed 
  wallTimerSet times;
  };

//some global info shared by all threads
extern paramfile       *g_params;
extern vec3 campos, lookat, sky;
extern std::vector<COLOURMAP> amap;
extern std::vector<particle_sim> *particle_data;
extern int ptypes;
extern float b_brightness;
extern wallTimerSet cuWallTimers;

int check_device(int rank);
void print_device_info(int rank, int dev);

void cuda_rendering(int mydevID, int nDev, arr2<COLOUR> &pic, std::vector<particle_sim> &particle, float brightness);
void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread);
//void cu_draw_chunk(void *pinfo, cu_gpu_vars* gv);
void setup_colormap(int ptypes, cu_gpu_vars* gv);

void GPUReport(wallTimerSet &cuTimers);
void cuda_timeReport(paramfile &params);

#endif
