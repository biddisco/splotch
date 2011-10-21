#ifndef NO_WIN_THREAD 
#include <pthread.h>
#endif

#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"
#include "cxxsupport/string_utils.h"
#include "cuda/CuPolicy.h"
#include "cuda/CuRender.h"

using namespace std;

paramfile *g_params;
float b_brightness;
int ptypes = 0;
vec3 campos, lookat, sky;
vector<COLOURMAP> amap;
vector<particle_sim> *particle_data;
wallTimerSet cuWallTimers;

void cuda_rendering(int mydevID, int nDev, arr2<COLOUR> &pic, vector<particle_sim> &particle, float brightness)
  {
  //see if host must be a working thread
  b_brightness = brightness;
  particle_data = &particle;
  bool bHostThread = g_params->find<bool>("use_host_as_thread", false);
  int nThread = bHostThread? nDev+1: nDev;
  //init array info for threads control
  thread_info *tInfo = new thread_info[nThread];
  tInfo[0].pPic = &pic;      //local var pic is assigned to the first thread
  tInfo[0].devID = mydevID;
//  tInfo[0].npart_all = npart_all;
  for (int i=1; i<nDev; i++)
    {
    tInfo[i].devID = mydevID+i;
//    tInfo[i].npart_all = npart_all;
    tInfo[i].pPic = new arr2<COLOUR>(pic.size1(), pic.size2());
    }
  //make the last one work for host thread
  if (bHostThread)
    {
    tInfo[nThread-1].devID =-1;
//    tInfo[nThread-1].npart_all = npart_all;
    if (nThread-1 != 0)
      tInfo[nThread-1].pPic = new arr2<COLOUR>(pic.size1(), pic.size2());
    }
  //decide how to divide particles range by another function
  DevideThreadsTasks(tInfo, nThread, bHostThread);

#ifndef NO_WIN_THREAD // create cuda threads on Windows using CreateThread function
  HANDLE *tHandle = new HANDLE[nThread];
  //issue the threads
  for (int i=0; i<nDev; i++)
    tHandle[i] = CreateThread( NULL, 0,
      (LPTHREAD_START_ROUTINE)cu_thread_func,&(tInfo[i]), 0, NULL );
  //issue the host thread too
  if (bHostThread)
    tHandle[nDev] = CreateThread( NULL, 0,
      (LPTHREAD_START_ROUTINE)host_thread_func,&(tInfo[nDev]), 0, NULL );
  WaitForMultipleObjects(nThread, tHandle, true, INFINITE);

#else // create cuda threads on Linux using pthread_create function

//  planck_assert(nDev <= 1, "can't have multiple cuda threads on Linux (yet), so 'gpu_number' must be 1");
  pthread_t *tHandle = new pthread_t[nThread];
  for (int i=0; i<nDev; i++)
     pthread_create(&(tHandle[i]), NULL, cu_thread_func, (void *) &(tInfo[i]) );
  if (bHostThread)
     pthread_create(&(tHandle[nDev]), NULL, host_thread_func, (void *) &(tInfo[nDev]) );
  void *status[nThread];
  for (int i=0; i <nThread; ++i) pthread_join(tHandle[i], &status[i]);
//  cu_thread_func (&(tInfo[0])); //just call it as normal function
#endif  //if not NO_WIN_THREAD

  // combine the results of multiple threads(devices + host) to pic
  for (int i=1; i<nThread; i++)
      for (int x=0; x<pic.size1(); x++)
        for (int y=0; y<pic.size2(); y++)
              pic[x][y] = pic[x][y] + (*tInfo[i].pPic)[x][y];

  if (g_params->getVerbosity())
   for (int i=0; i<nThread; i++)
    {
    if (tInfo[i].devID!=-1)
      {
      cout<< endl <<"Rank " << mpiMgr.rank() << ": Times of GPU" << i << ":" <<endl;
      GPUReport(tInfo[i].times);
      cout<<endl;
      }
    }
  
  if (mpiMgr.master()) cuWallTimers = tInfo[0].times;

  for (int i=1; i<nThread; i++)
    delete tInfo[i].pPic;
  delete [] tInfo;
  delete [] tHandle;
  }


void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread)
  {
  bool bTestLoadBalancing = g_params->find<bool>("test_load_balancing", false);
  unsigned int curStart = 0;
  int hostLoad = bHostThread? g_params->find<int>("host_load",0): 0;
  int nDev = bHostThread? nThread-1: nThread;
  int onePercent = particle_data->size()/100;
  int averageDevLen = (nDev!=0)? onePercent *(100-hostLoad)/nDev : 0;

  for (int i=0; i<nThread; i++)
    {
    tInfo[i].startP = curStart;
    if (tInfo[i].devID != -1) //not a host
      {
      if (bTestLoadBalancing)
        {
        int gpuLoad = g_params->find<int>("gpu_load"+dataToString(i),0);
        tInfo[i].endP = curStart + gpuLoad * onePercent - 1;
        }
      else
        tInfo[i].endP = curStart + averageDevLen - 1;
      }
    else //if this is a host
      {
      tInfo[i].endP = curStart + hostLoad * onePercent - 1;
      }
    curStart = tInfo[i].endP + 1;
    }

  tInfo[nThread-1].endP = particle_data->size()-1;
  }



THREADFUNC cu_thread_func(void *pinfo)
 {
  //a new thread info object that will carry each chunk's drawing
  thread_info *pInfoOutput = (thread_info*) pinfo;
  thread_info ti = *pInfoOutput;

  //do some cleaning for final thread_info
  pInfoOutput->pPic->fill(COLOUR(0.0, 0.0, 0.0));
  int xres = pInfoOutput->pPic->size1();
  int yres = pInfoOutput->pPic->size2();

  //a new linear pic object that will carry the result
  COLOUR Pic_chunk[xres*yres];
  //arr2<COLOUR> pic(pInfoOutput->pPic->size1(), pInfoOutput->pPic->size2());
  //ti.pPic = &pic;

  // Initialize policy class
  CuPolicy *policy = new CuPolicy(*g_params); 

  // num particles to manage at once
  float factor = g_params->find<float>("particle_mem_factor", 3);
  int len = cu_get_chunk_particle_count(policy, sizeof(cu_particle_sim), factor); 
  if (len == 0)
    {
    printf("\nGraphics memory setting error\n");
    mpiMgr.abort();
    }

  //CUDA Init
  cu_gpu_vars gv; //for each gpu a variable pack is needed
  memset(&gv, 0, sizeof(cu_gpu_vars));
  gv.policy = policy;
  // enable device and allocate arrays
  int error;
  error = cu_init(pInfoOutput->devID, len, &gv, *g_params, campos, lookat, sky, b_brightness);
  if (error)
  {
    cout << "Device Memory: allocation error!" << endl;
  }
  else
  {
    setup_colormap(ptypes, &gv);

    float64 grayabsorb = g_params->find<float>("gray_absorption",0.2);
    bool a_eq_e = g_params->find<bool>("a_eq_e",true);

    int endP = ti.endP;
    ti.endP = ti.startP;
    while(ti.endP < endP)
    {
     ti.endP = ti.startP + len - 1;   //set range
     if (ti.endP > endP) ti.endP = endP; 
     cu_draw_chunk(&ti, Pic_chunk, &gv, a_eq_e, grayabsorb);
     // combine results of chunks
     pInfoOutput->times.start("gcombine");
     for (int x=0; x<xres; x++)
      for (int y=0; y<yres; y++)
        (*(pInfoOutput->pPic))[x][y] += Pic_chunk[x*yres+y];
     pInfoOutput->times.stop("gcombine");
     ti.startP = ti.endP + 1;
    }
    pInfoOutput->times = ti.times;
    cu_end(&gv);
  }
 }


void setup_colormap(int ptypes, cu_gpu_vars* gv)
{
//init C style colormap
  cu_color_map_entry *amapD;//amap for Device
  int *amapDTypeStartPos; //begin indexes of ptypes in the linear amapD[]
  amapDTypeStartPos =new int[ptypes];
  int curPtypeStartPos =0;
  int size =0;
  //first we need to count all the entries to get colormap size
  for (int i=0; i<amap.size(); i++)
    size += amap[i].size();

  //then fill up the colormap amapD
  amapD =new cu_color_map_entry[size];
  int j,index =0;
  for(int i=0; i<amap.size(); i++)
    {
    for (j=0; j<amap[i].size(); j++)
      {
      amapD[index].val = amap[i].getX(j);
      COLOUR c (amap[i].getY(j));
      amapD[index].color.r = c.r;
      amapD[index].color.g = c.g;
      amapD[index].color.b = c.b;
      index++;
      }
    amapDTypeStartPos[i] = curPtypeStartPos;
    curPtypeStartPos += j;
    }
  //now let cuda init colormap on device
  cu_colormap_info clmp_info;
  clmp_info.map =amapD;
  clmp_info.mapSize =size;
  clmp_info.ptype_points =amapDTypeStartPos;
  clmp_info.ptypes =ptypes;
  cu_init_colormap(clmp_info, gv);

  delete []amapD;
  delete []amapDTypeStartPos;
}



THREADFUNC host_thread_func(void *p)
  {
  thread_info *tInfo = (thread_info*)p;

  vector<particle_sim>::iterator i1,i2;
  i1 =particle_data->begin() + tInfo->startP;
  i2 =particle_data->begin() + tInfo->endP + 1;
  vector<particle_sim> particles(i1,i2);

  host_rendering(*g_params, particles, *(tInfo->pPic), campos, lookat, sky, amap, b_brightness);
  }


void GPUReport(wallTimerSet &cuTimers)
  {
    cout << "Copy  (secs)               : " << cuTimers.acc("gcopy") << endl;
    cout << "Transforming Data (secs)   : " << cuTimers.acc("gtransform") << endl;
    cout << "Sorting Fragments (secs)   : " << cuTimers.acc("gsort") << endl;
    cout << "Reduce Fragments (secs)    : " << cuTimers.acc("greduce") << endl;
    cout << "Filter Sub-Data (secs)     : " << cuTimers.acc("gfilter") << endl;
    cout << "Rendering Sub-Data (secs)  : " << cuTimers.acc("grender") << endl;
    cout << "Cuda thread (secs)         : " << cuTimers.acc("gpu_thread") << endl;
  }

void cuda_timeReport(paramfile &params)
  {
  if (mpiMgr.master())
    {
    cout << endl << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "--------------------------------------------" << endl;
    cout<< endl <<"Times of GPU:" <<endl;
    GPUReport (cuWallTimers);
    cout <<  "--------------------------------------------" << endl;
    }
  }

