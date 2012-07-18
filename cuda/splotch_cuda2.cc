// each mpi_task manage 1 GPU
#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"
#include "cxxsupport/string_utils.h"
#include "cuda/CuPolicy.h"
#include "cuda/CuRender.h"

using namespace std;

wallTimerSet cuWallTimers;

void cuda_rendering(int mydevID, arr2<COLOUR> &pic, vector<particle_sim> &particle, const vec3 &campos, const vec3 &lookat, vec3 &sky, vector<COLOURMAP> &amap, float b_brightness, paramfile &g_params)
  {
  long int nP = particle.size();
  pic.fill(COLOUR(0.0, 0.0, 0.0));
  int xres = pic.size1();
  int yres = pic.size2();
  cout << "resolution = " << xres << " x " << yres << endl;
  arr2<COLOUR> Pic_host(xres,yres);
  int ptypes = g_params.find<int>("ptypes",1);

  // Initialize policy class
  CuPolicy *policy = new CuPolicy(g_params); 

  // num particles to manage at once
  float factor = g_params.find<float>("particle_mem_factor", 4);
  long int len = cu_get_chunk_particle_count(policy, sizeof(cu_particle_sim), factor);

  if (len <= 0)
    {
    printf("\nGraphics memory setting error\n");
    mpiMgr.abort();
    }

  //CUDA Init
  cu_gpu_vars gv; 
  memset(&gv, 0, sizeof(cu_gpu_vars));
  gv.policy = policy;
  // enable device and allocate arrays
  int error;
  error = cu_init(mydevID, len, &gv, g_params, campos, lookat, sky, b_brightness);
  if (!error)
  {
    //a new linear pic object that will carry the result
    COLOUR Pic_chunk[xres*yres];
    setup_colormap(ptypes, amap, &gv);

    float64 grayabsorb = g_params.find<float>("gray_absorption",0.2);
    bool a_eq_e = g_params.find<bool>("a_eq_e",true);
 
    int endP = 0;
    int startP = 0;
    int nPR = 0;
    cuWallTimers.start("gpu_draw");
    while(endP < nP)
    {
     endP = startP + len;   //set range
     if (endP > nP) endP = nP; 
     nPR += cu_draw_chunk(&cuWallTimers, mydevID, (cu_particle_sim *) &(particle[startP]), endP-startP, Pic_chunk, Pic_host, &gv, a_eq_e, grayabsorb, b_brightness, amap, g_params);
     // combine results of chunks
     cuWallTimers.start("gcombine");
     for (int x=0; x<xres; x++)
      for (int y=0; y<yres; y++)
        pic[x][y] += Pic_chunk[x*yres+y]+Pic_host[x][y];
     cuWallTimers.stop("gcombine");
     cout << "Rank " << mpiMgr.rank() << ": Rendered " << nPR << "/" << nP << " particles" << endl;
     startP = endP;
    }
    cuWallTimers.stop("gpu_draw");
    cu_endThread(&gv);
  }

  if (g_params.getVerbosity())
  { 
      cout<< endl <<"Rank " << mpiMgr.rank() << ": Times of GPU" << mydevID << ":" <<endl;
      GPUReport(cuWallTimers);
      cout<<endl;
   }
  //if (mpiMgr.master()) cuWallTimers = cu Timers;
 }


void setup_colormap(int ptypes, vector<COLOURMAP> &amap, cu_gpu_vars* gv)
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
    cout << "Project+coloring Data(secs): " << cuTimers.acc("gprocess") << endl;
    cout << "Filter Sub-Data (secs)     : " << cuTimers.acc("gfilter") << endl;
    cout << "Rendering Sub-Data (secs)  : " << cuTimers.acc("grender") << endl;
    cout << "Sorting Fragments (secs)   : " << cuTimers.acc("gsort") << endl;
    cout << "Reduce Fragments (secs)    : " << cuTimers.acc("greduce") << endl;
    cout << "Combine images (secs)      : " << cuTimers.acc("gcombine") << endl;
    cout << "Cuda draw   (secs)         : " << cuTimers.acc("gpu_draw") << endl << endl;
    cout << "Host rendering (secs)      : " << cuTimers.acc("host_rendering") << endl;
  }

void cuda_timeReport()
  {
  if (mpiMgr.master())
    {
    cout << endl << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "--------------------------------------------" << endl;
    cout<< endl <<"Times of GPU:" <<endl;
    GPUReport(cuWallTimers);
    cout <<  "--------------------------------------------" << endl;
    }
  }

