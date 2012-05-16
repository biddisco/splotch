// each mpi_task manage 1 GPU
#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"
#include "cxxsupport/string_utils.h"
#include "cuda/CuPolicy.h"
#include "cuda/CuRender.h"

using namespace std;

paramfile *g_params;
int ptypes = 0;
vec3 campos, lookat, sky;
vector<COLOURMAP> amap;
vector<particle_sim> *particle_data;
wallTimerSet cuWallTimers;

void cuda_rendering(int mydevID, arr2<COLOUR> &pic, vector<particle_sim> &particle, float b_brightness)
  {
  wallTimerSet times;
  particle_data = &particle;
  long int nP = particle_data->size();
  //do some cleaning for final thread_info
  pic.fill(COLOUR(0.0, 0.0, 0.0));
  int xres = pic.size1();
  int yres = pic.size2();
  arr2<COLOUR> Pic_host(xres,yres);

  // Initialize policy class
  CuPolicy *policy = new CuPolicy(*g_params); 

  // num particles to manage at once
  float factor = g_params->find<float>("particle_mem_factor", 4);
  long int len = cu_get_chunk_particle_count(policy, sizeof(cu_particle_sim), factor);

  if (len <= 0)
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
  error = cu_init(mydevID, len, &gv, *g_params, campos, lookat, sky, b_brightness);
  if (!error)
  {
    //a new linear pic object that will carry the result
    COLOUR Pic_chunk[xres*yres];
    setup_colormap(ptypes, &gv);

    float64 grayabsorb = g_params->find<float>("gray_absorption",0.2);
    bool a_eq_e = g_params->find<bool>("a_eq_e",true);

    int endP = 0;
    int startP = 0;
    while(endP < nP)
    {
     endP = startP + len;   //set range
     if (endP > nP) endP = nP; 
     cu_draw_chunk(&times, mydevID, startP, endP, Pic_chunk, Pic_host, &gv, a_eq_e, grayabsorb, b_brightness);
     // combine results of chunks
     times.start("gcombine");
     for (int x=0; x<xres; x++)
      for (int y=0; y<yres; y++)
        pic[x][y] += Pic_chunk[x*yres+y]+Pic_host[x][y];
     times.stop("gcombine");
     startP = endP;
    }
    cu_endThread(&gv);
  }

  if (g_params->getVerbosity())
  { 
      cout<< endl <<"Rank " << mpiMgr.rank() << ": Times of GPU" << mydevID << ":" <<endl;
      GPUReport(times);
      cout<<endl;
   }
  if (mpiMgr.master()) cuWallTimers = times;
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


void GPUReport(wallTimerSet &cuTimers)
  {
    cout << "Copy  (secs)               : " << cuTimers.acc("gcopy") << endl;
    cout << "Transforming Data (secs)   : " << cuTimers.acc("gtransform") << endl;
    cout << "Filter Sub-Data (secs)     : " << cuTimers.acc("gfilter") << endl;
    cout << "Colorize Sub-Data (secs)   : " << cuTimers.acc("gcolor") << endl;
    cout << "Rendering Sub-Data (secs)  : " << cuTimers.acc("grender") << endl;
    cout << "Sorting Fragments (secs)   : " << cuTimers.acc("gsort") << endl;
    cout << "Reduce Fragments (secs)    : " << cuTimers.acc("greduce") << endl;
    cout << "Combine images (secs)      : " << cuTimers.acc("gcombine") << endl;
    cout << "Cuda thread (secs)         : " << cuTimers.acc("gpu_thread") << endl << endl;
    cout << "Host rendering (secs)      : " << cuTimers.acc("host_rendering") << endl;
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

