#ifndef SPLOTCH_CUDA2_H
#define SPLOTCH_CUDA2_H

#include "cuda/splotch_cuda.h"
#include "splotch/splotch_host.h"
#include "cxxsupport/walltimer.h"

//things for combination with host threads
struct  param_combine_thread //for host combine thread
  {
  // bool    bFinished; used in thread combine. not working.
  bool a_eq_e;
  void *fbuf;
  int combineStartP, combineEndP;
  cu_particle_splotch *ps;
  float timeUsed;
  arr2<COLOUR> *pPic;
  };

//#ifndef NO_WIN_THREAD
THREADFUNC combine(void *param);
//#endif
THREADFUNC cu_thread_func(void *pinfo);
THREADFUNC cu_draw_chunk(void *pinfo);
THREADFUNC host_thread_func(void *pinfo);

// struct containing thread task info
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
extern std::vector<particle_sim> particle_data; //raw data from file
extern vec3 campos, lookat, sky;
extern std::vector<COLOURMAP> amap;
extern int ptypes;
extern wallTimerSet cuWallTimers;

int check_device(int rank);
void device_info(int rank, int dev);
void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread);

void cuda_rendering(int mydevID, int nDev, int res, arr2<COLOUR> &pic);

void GPUReport(wallTimerSet &cuTimers);
void cuda_timeReport(paramfile &params);

#endif
