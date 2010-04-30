#ifndef NO_WIN_THREAD 
#include <pthread.h>
#endif

#include "cuda/splotch_cuda2.h"

using namespace std;

paramfile *g_params;
int ptypes = 0;
vector<particle_sim> particle_data;   //raw data from file
vec3 campos, lookat, sky;
vector<COLOURMAP> amap;
wallTimerSet cuWallTimers;

void cuda_rendering(int mydevID, int nDev, int res, arr2<COLOUR> &pic)
  {
  //see if host must be a working thread
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
    tInfo[i].pPic = new arr2<COLOUR>(res, res);
    }
  //make the last one work for host thread
  if (bHostThread)
    {
    tInfo[nThread-1].devID =-1;
//    tInfo[nThread-1].npart_all = npart_all;
    if (nThread-1 != 0)
      tInfo[nThread-1].pPic = new arr2<COLOUR>(res, res);
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

  // combine the results of multiple threads to pic
  for (int i=1; i<nThread; i++)
      for (int x=0; x<res; x++) //  error when x =1,
        for (int y=0; y<res; y++)
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
/*    else
      {
      cout<< endl <<"Times of CPU " << mpiMgr.rank() << " as a thread:" <<endl;
      cout<<endl;
      } */
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
  int curStart = 0;
  int hostLoad = bHostThread? g_params->find<int>("host_load",0): 0;
  int nDev = bHostThread? nThread-1: nThread;
  int onePercent = particle_data.size()/100;
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

 // tInfo[nThread-1].endP = particle_data.size()-1;
  }



THREADFUNC cu_thread_func(void *pinfo)
  {
  //a new thread info object that will carry each chunk's drawing
  thread_info *pInfoOutput = (thread_info*) pinfo;
  thread_info ti = *pInfoOutput;

  //do some cleaning for final thread_info
  pInfoOutput->pPic->fill(COLOUR(0.0, 0.0, 0.0));

  //a new pic object residing in ti that will carry the result
  arr2<COLOUR> pic(pInfoOutput->pPic->size1(), pInfoOutput->pPic->size2());
  ti.pPic = &pic;

  // num particles to manage at once
  int len = cu_get_chunk_particle_count(*g_params); 
  if (len == -1)
    {
    printf("\nGraphics memory setting error\n");
    mpiMgr.abort();
    }

  //CUDA Init
  cu_init(pInfoOutput->devID);

  int curEnd = 0;
  int endP = ti.endP;

  ti.endP = ti.startP;
  while(ti.endP < endP)
    {
    ti.endP =ti.startP + len - 1;   //set range
    if (ti.endP > endP) ti.endP = endP;
    cu_draw_chunk(&ti);
    // combine results of chunks
    for (int x=0; x<pic.size1(); x++)
      for (int y=0; y<pic.size2(); y++)
        (*(pInfoOutput->pPic))[x][y] += pic[x][y];
    ti.startP = ti.endP +1;
    pInfoOutput->times = ti.times;
    }
  }


THREADFUNC cu_draw_chunk(void *pinfo)
  {

  //get the input info
  thread_info *tInfo = (thread_info*)pinfo;
  tInfo->times.start("gpu_thread");

  int nParticle =tInfo->endP -tInfo->startP +1;
  printf("Rank %d - GPU %d : Processing %d particles\n", mpiMgr.rank(), tInfo->devID, nParticle);

  paramfile &params (*g_params);

  //for each gpu/thread a variable pack is needed
  cu_gpu_vars gv;
  memset(&gv, 0, sizeof(cu_gpu_vars));

  //copy data to local C-like array d_particle_data
  tInfo->times.start("gcopy");
  cu_particle_sim *d_particle_data = new cu_particle_sim[nParticle];
  memcpy( &(d_particle_data[0]), &(particle_data[tInfo->startP]),
          nParticle*sizeof(cu_particle_sim));
  tInfo->times.stop("gcopy"); 

  //here we analyse how to divide the whole task for large data handling

  //CUDA Ranging
  tInfo->times.start("grange");
  cu_range(params, d_particle_data, nParticle, &gv);
  tInfo->times.stop("grange");

  //CUDA Transformation
  tInfo->times.start("gtransform");
  double  c[3]={campos.x, campos.y, campos.z},
          l[3]={lookat.x, lookat.y, lookat.z},
          s[3]={sky.x, sky.y,     sky.z};
  cu_transform(params, nParticle,c, l, s,d_particle_data, &gv);
  tInfo->times.stop("gtransform");

/* temporarily ignore sorting 191109.
   it becomes complicated when using multiple threads with sorting
   //then copy particles back to host for sorting
   for (int i=0; i<nParticle; i++)
   memcpy( &(particle_data[iWRONG]),&(d_particle_data[i]), sizeof(cu_particle_sim));
PROBLEM HERE!

// --------------------------------
// ----------- Sorting ------------
// --------------------------------
// cout << endl << "applying sort (" << npart << ") ..." << endl;

   int sort_type = params.find<int>("sort_type",1);
   particle_sort(particle_data,sort_type,true);
     //we can sort by size(r) for balancing with cuda
     //particle_sort(particle_data,4,true);

     //copy sorted data back to device, not a must!
     //first to C-style object
     for(int i=0; i<particle_data.size(); i++)
     memcpy( &(d_particle_data[i]), &(particle_data[i]), sizeof(cu_particle_sim));
    cu_copy_particle_sim_to_device(d_particle_data, particle_data.size());
*/

  //CUDA Coloring
  tInfo->times.start("gcoloring");
  //init C style colormap
  cu_color_map_entry *amapD;//amap for Device. emap is not used currently
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
    amapDTypeStartPos[i] =curPtypeStartPos;
    curPtypeStartPos += j;
    }
  //now let cuda init colormap on device
  cu_colormap_info clmp_info;
  clmp_info.map =amapD;
  clmp_info.mapSize =size;
  clmp_info.ptype_points =amapDTypeStartPos;
  clmp_info.ptypes =ptypes;
  cu_init_colormap(clmp_info, &gv);

  //init cu_particle_splotch array memeory
  cu_particle_splotch *cu_ps;
  size =nParticle;
  cu_ps =new cu_particle_splotch[size];
  memset(cu_ps, 0, size);

  //Colorize with device
  cu_colorize(params, cu_ps, size, &gv);
  tInfo->times.stop("gcoloring");

//////////////////////////////////////////////////////////////////////
  tInfo->times.start("gfilter");
  //filter particle_splotch array to a cu_ps_filtered
  int pFiltered=0;
  //for sorting and splitting
  vector <cu_particle_splotch> v_ps;
  cu_particle_splotch p;
  //do filtering
  unsigned long posInFragBuf =0;
  int minx=1e6,miny=1e6, maxx=-1,maxy=-1;

  //old code observ size
  //select valid ones
  for (int i=0; i<size; i++)
    {
    if (cu_ps[i].isValid)
      {
      //IMPORTANT: set the start position of this particle in fragment buffer
      // cu_ps[i].posInFragBuf =posInFragBuf;
      // posInFragBuf += (cu_ps[i].maxx -cu_ps[i].minx)*
      // (cu_ps[i].maxy -cu_ps[i].miny); SHOULD DO AFTER SORTING!
      // memcpy(&cu_ps_filtered[pFiltered],&cu_ps[i], sizeof(cu_particle_splotch));

      memcpy(&p, &cu_ps[i], sizeof(cu_particle_splotch));
      v_ps.push_back(p);
      pFiltered++;

      minx=min(minx,(int)cu_ps[i].minx);
      miny=min(miny,(int)cu_ps[i].miny);
      maxx=max(maxx,(int)cu_ps[i].maxx);
      maxy=max(maxy,(int)cu_ps[i].maxy);
      }
    }
 
  int maxRegion =cu_get_max_region(&gv);
  vector<cu_particle_splotch> v_ps1;//result goes to it
  for (vector<cu_particle_splotch>::iterator i=v_ps.begin(); i<v_ps.end(); i++)
    {
    int h, w;
    cu_particle_splotch     p, pNew;
    p =(*i);
    h = p.maxy - p.miny;
    w = p.maxx - p.minx;

    if (h*w <maxRegion)
    // if ( 1)//no splitting test
      {
      v_ps1.push_back(p);
      continue;
      }

    //now we split
    int w1 = (maxRegion %h==0) ? (maxRegion /h):(maxRegion /h +1);
    //insert new cells
    pNew =p;
    //minx,maxx of pNew need to be set
    for (int minx =p.minx; minx<p.maxx; minx+=w1)
      {
      pNew.minx =minx;
      pNew.maxx =( minx+w1 >=p.maxx) ? p.maxx : minx+w1;
      v_ps1.push_back(pNew);
      }
    }
  tInfo->times.stop("gfilter");

  tInfo->times.start("gsort");
  v_ps.clear();//not useful any more
  //sort the filtered,splitted v_ps
  // sort(v_ps1.begin(), v_ps1.end(), region_cmp());
  tInfo->times.stop("gsort");

  tInfo->times.start("gfilter");
  //copy to C-style array cu_ps_filtered
  cu_particle_splotch *cu_ps_filtered;
  size =v_ps1.size();
  cu_ps_filtered =new cu_particle_splotch[size];
  pFiltered =v_ps1.size();
  for(int i=0; i<pFiltered; i++)
    {
    v_ps1[i].posInFragBuf =posInFragBuf;
    int     region =(v_ps1[i].maxx -v_ps1[i].minx)*(v_ps1[i].maxy -v_ps1[i].miny);
    posInFragBuf +=region;
    cu_ps_filtered[i] =v_ps1[i];
    }
  tInfo->times.stop("gfilter");

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------

  //get parameters for rendering
  int res = params.find<int>("resolution",200);
  long nsplotch=pFiltered;
  long nsplotch_all=nsplotch;
  mpiMgr.allreduce(nsplotch_all,MPI_Manager::Sum);
  float64 grayabsorb = params.find<float>("gray_absorption",0.2);
  bool a_eq_e = params.find<bool>("a_eq_e",true);

  //CUDA Rendering with device

  //here's the point of multi-go loop starts
  tInfo->times.start("grender");
  //prepare fragment buffer memory space first
  cu_fragment_AeqE  *fragBufAeqE;
  cu_fragment_AneqE *fragBufAneqE;
  int nFBufInByte =cu_get_fbuf_size(&gv);
  int nFBufInCell;
  if (a_eq_e)
    {
    nFBufInCell =nFBufInByte/sizeof(cu_fragment_AeqE);
    fragBufAeqE =new cu_fragment_AeqE[nFBufInCell];
    }
  else
    {
    nFBufInCell =nFBufInByte/sizeof(cu_fragment_AneqE);
    fragBufAneqE =new cu_fragment_AneqE[nFBufInCell];
    }

  cu_prepare_render(cu_ps_filtered,pFiltered, &gv);

  //clear the output pic
  tInfo->pPic ->fill(COLOUR(0.0, 0.0, 0.0));

  //initialize combine vars
  //cu_ps_filtered:       the particle array
  //pFiltered:            length of particle array
  //pPosInFragBuf:        length of the whole fragment buffer counted after filterign

  int renderStartP, renderEndP;
  renderEndP = 0;
  do {
    //find a chunk: compute number of fragments needed to render this chunk
    int nFragments2Render = 0;
    renderStartP =renderEndP;
    while (renderEndP<pFiltered && nFragments2Render < nFBufInCell)
      {
       int inc =(cu_ps_filtered[renderEndP].maxx-cu_ps_filtered[renderEndP].minx) *
           (cu_ps_filtered[renderEndP].maxy -cu_ps_filtered[renderEndP].miny);
       nFragments2Render +=inc;
       renderEndP++;
      }

    //render it
    cu_render1(renderStartP, renderEndP, a_eq_e, grayabsorb, &gv);

    //collect result
    cu_get_fbuf(fragBufAeqE, a_eq_e, nFragments2Render, &gv);

    //combine chunks
    tInfo->times.start("gcombine");
    if (a_eq_e)
    {
      for (int pPos=renderStartP, fPos=0; pPos<renderEndP; pPos++)
      {
        for (int x =cu_ps_filtered[pPos].minx; x <cu_ps_filtered[pPos].maxx; x++)
        {
          for (int y =cu_ps_filtered[pPos].miny; y <cu_ps_filtered[pPos].maxy; y++)
          {
           (*tInfo->pPic)[x][y].r += fragBufAeqE[fPos].deltaR;
           (*tInfo->pPic)[x][y].g += fragBufAeqE[fPos].deltaG;
           (*tInfo->pPic)[x][y].b += fragBufAeqE[fPos].deltaB;
           fPos ++;
          }
        }
      }
    }
    else
    {
      for (int pPos=renderStartP, fPos=0; pPos<renderEndP; pPos++)
      {
        for (int x =cu_ps_filtered[pPos].minx; x <cu_ps_filtered[pPos].maxx; x++)
        {
          for (int y =cu_ps_filtered[pPos].miny; y <cu_ps_filtered[pPos].maxy; y++) 
          {
           (*tInfo->pPic)[x][y].r = (*tInfo->pPic)[x][y].r *fragBufAneqE[fPos].factorR 
                                   +fragBufAneqE[fPos].deltaR;
           (*tInfo->pPic)[x][y].g = (*tInfo->pPic)[x][y].g *fragBufAneqE[fPos].factorG 
                                   +fragBufAneqE[fPos].deltaG;
           (*tInfo->pPic)[x][y].b = (*tInfo->pPic)[x][y].b *fragBufAneqE[fPos].factorB 
                                   +fragBufAneqE[fPos].deltaB;
           fPos ++;
          }
        }
      }
    }
    tInfo->times.stop("gcombine");

  } while (renderEndP < pFiltered);

  tInfo->times.stop("grender");

/////////////////////////////////////////////////////////////////////

  cu_end(&gv);
  delete []d_particle_data;
  delete []amapD;
  delete []amapDTypeStartPos;
  delete []cu_ps;
  delete []cu_ps_filtered;
  if (a_eq_e)
    delete []fragBufAeqE;
  else
    delete []fragBufAneqE;

  tInfo->times.stop("gpu_thread");
  }


THREADFUNC host_thread_func(void *p)
  {
  thread_info *tInfo = (thread_info*)p;

  vector<particle_sim>::iterator i1,i2;
  i1 =particle_data.begin() + tInfo->startP;
  i2 =particle_data.begin() + tInfo->endP + 1;
  vector<particle_sim> particles(i1,i2);

  host_rendering(*g_params, particles, *(tInfo->pPic), campos, lookat, sky, amap);
  }


void GPUReport(wallTimerSet &cuTimers)
  {
    cout << "Copy2C_like (secs)         : " << cuTimers.acc("gcopy") << endl;
    cout << "Ranging Data (secs)        : " << cuTimers.acc("grange") << endl;
    cout << "Transforming Data (secs)   : " << cuTimers.acc("gtransform") << endl;
//    cout << "Sorting Data (secs)        : " << cuTimers.acc("gsort") << endl;
    cout << "Coloring Sub-Data (secs)   : " << cuTimers.acc("gcoloring") << endl;
    cout << "Filter Sub-Data (secs)     : " << cuTimers.acc("gfilter") << endl;
    cout << "Rendering Sub-Data (secs)  : " << cuTimers.acc("grender") << endl;
    cout << "Cuda thread (secs)         : " << cuTimers.acc("gpu_thread") << endl;
  }

void cuda_timeReport(paramfile &params)
  {
  if (mpiMgr.master())
    {
    wallTimers.stop("full");
    cout << endl << "--------------------------------------------" << endl;
    cout << "Summary of timings" << endl;
    cout << "--------------------------------------------" << endl;
    cout<< endl <<"Times of GPU:" <<endl;
    GPUReport (cuWallTimers);
    cout <<  "--------------------------------------------" << endl;

    if (params.find<bool>("use_host_as_thread", false))
      {
      cout<< endl <<"Times of CPU HOST as threads:" <<endl;
      hostTimeReport(wallTimers);
      cout << "--------------------------------------------" << endl;
      }
    cout << "Setup Data (secs)          : " << wallTimers.acc("setup") << endl;
    cout << "Read Data (secs)           : " << wallTimers.acc("read") << endl;
    cout << "Postprocessing (secs)      : " << wallTimers.acc("postproc") << endl;
    cout << "Write Data (secs)          : " << wallTimers.acc("write") << endl;
    cout << "Total (secs)               : " << wallTimers.acc("full") << endl;
    }
  }

