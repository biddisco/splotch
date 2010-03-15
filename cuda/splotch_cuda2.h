#ifndef SPLOTCH_CUDA2_H
#define SPLOTCH_CUDA2_H

#ifdef CUDA
#include "cuda/splotch_cuda.h"
#include "cuda/splotchutils_cuda.h"
#include <string.h>
#endif

#ifdef CUDA
#ifndef VS
#define DWORD long
#define WINAPI
#endif

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

#ifndef NO_WIN_THREAD
DWORD WINAPI combine(void *param);
#else
DWORD WINAPI cu_thread_func(void *pinfo);
DWORD WINAPI cu_draw_chunk(void *pinfo);
#endif

//for record times
enum TimeRecords {
  CUDA_INIT,
  COPY2C_LIKE,
  RANGE,
  TRANSFORMATION,
  COLORIZE,
  FILTER,
  SORT,
  RENDER,
  COMBINE,
  THIS_THREAD,
  TIME_RECORDS   //to indicate number of times
  };

// struct containing thread task info
struct thread_info
  {
  int devID;                  //index of the device selected
  int startP, endP;           //start and end particles to handle
  arr2<COLOUR> *pPic;         //the output image
  float times[TIME_RECORDS];  //carry out times of computing
  };

//some global info shared by all threads
extern paramfile       *g_params;
extern vector<particle_sim> particle_data; //raw data from file
extern vec3 campos, lookat, sky;
extern vector<COLOURMAP> amap,emap;
extern int ptypes;

void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread);

#ifndef NO_WIN_THREAD
DWORD WINAPI cu_draw_chunk(void *p);
DWORD WINAPI cu_thread_func(void *p);
DWORD WINAPI host_thread_func(void *p);
#endif

void render_cuda(paramfile &params, int &res, arr2<COLOUR> &pic);

#endif //ifdef CUDA

#endif
