#include "cuda/splotch_cuda2.h"
#include "cxxsupport/walltimer.h"

using namespace std;

paramfile       *g_params;
vector<particle_sim> particle_data; ///raw data from file
vec3 campos, lookat, sky;
vector<COLOURMAP> amap,emap;
int ptypes = 0;

THREADFUNC host_thread_func(void *p)
  {
  printf("\nHost Thread Start!\n");

  thread_info *tInfo = (thread_info*)p;

  paramfile params =*g_params;

  vector<particle_sim> particles;
  vector<particle_sim>::iterator i1,i2;
  i1 =particle_data.begin() +tInfo->startP;
  i2 =particle_data.begin() +tInfo->endP;
  particles.assign(i1, i2);

  wallTimer t, t1;
  t1.start();
  t.start();

  particle_normalize(params,particles,false);

  t.stop();
  tInfo->times[RANGE] =t.acc();
  t.reset();
  t.start();

  particle_project(params, particles, campos, lookat, sky);

  t.stop();
  tInfo->times[TRANSFORMATION] =t.acc();

  // ----------- Sorting ------------
  // NO SORING FOR NOW
  //int sort_type = params.find<int>("sort_type",1);
  //particle_sort(particle_data,sort_type,true);

  t.reset();
  t.start();

  // particle_colorize(params, particles, particle_col, amap, emap);
  particle_colorize(params, particles, amap, emap); //new calling

  t.stop();
  tInfo->times[COLORIZE]=t.acc();
  t.reset();
  t.start();

  int res = params.find<int>("resolution",200);
  float64 grayabsorb = params.find<float>("gray_absorption",0.2);
  bool a_eq_e = params.find<bool>("a_eq_e",true);
  render(particles,*(tInfo->pPic),a_eq_e,grayabsorb,true);

  t.stop();
  tInfo->times[RENDER] =t.acc();
  t1.stop();
  tInfo->times[THIS_THREAD] =t1.acc();

  printf("\nHost Thread End!\n");
  return 1;
  }


THREADFUNC combine (void *param1)
  {
  static int enter=-1;

  param_combine_thread *param = (param_combine_thread*)param1;

  cu_particle_splotch *ps = param->ps;
  arr2<COLOUR> *pPic = param->pPic;

  wallTimer t;
  t.start();

  if (param->a_eq_e)
    {
    cu_fragment_AeqE *bufWrite =(cu_fragment_AeqE*)param->fbuf;

    for (int pPos=param->combineStartP, fPos=0;
         pPos<param->combineEndP; pPos++)
      {
      for (int x =ps[pPos].minx; x <ps[pPos].maxx; x++)
        {
        for (int y =ps[pPos].miny; y <ps[pPos].maxy; y++)
          {
          (*pPic)[x][y].r += bufWrite[fPos].deltaR;
          (*pPic)[x][y].g += bufWrite[fPos].deltaG;
          (*pPic)[x][y].b += bufWrite[fPos].deltaB;
          fPos ++;
          }
        }
      }
    }
  else
    {
    cu_fragment_AneqE *bufWrite =(cu_fragment_AneqE*)param->fbuf;

    for (int pPos=param->combineStartP, fPos=0;
         pPos<param->combineEndP; pPos++)
      {
      for (int x =ps[pPos].minx; x <ps[pPos].maxx; x++)
        {
        for (int y =ps[pPos].miny; y <ps[pPos].maxy; y++)
          {
          (*pPic)[x][y].r = (*pPic)[x][y].r *bufWrite[fPos].factorR +bufWrite[fPos].deltaR;
          (*pPic)[x][y].g = (*pPic)[x][y].g *bufWrite[fPos].factorG +bufWrite[fPos].deltaG;
          (*pPic)[x][y].b = (*pPic)[x][y].b *bufWrite[fPos].factorB +bufWrite[fPos].deltaB;
          fPos ++;
          }
        }
      }
    }

  t.stop();
  param->timeUsed +=t.acc();

  return 1;
  }

THREADFUNC cu_thread_func(void *pinfo)
  {
  //a new thread info object that will carry each chunk's drawing
  thread_info *pInfoOutput=(thread_info*) pinfo,
              ti =*pInfoOutput;

  //do some cleaning for final thread_info
  pInfoOutput->pPic->fill(COLOUR(0.0, 0.0, 0.0));
  memset(pInfoOutput->times, 0, sizeof(float)*TIME_RECORDS);

  //a new pic object residing in ti that will carry the result
  arr2<COLOUR> pic(pInfoOutput->pPic->size1(), pInfoOutput->pPic->size2());
  ti.pPic =&pic;

  //set startP and end P of ti
  int len = cu_get_chunk_particle_count(*g_params);
  if (len==-1)
    {
    printf("\nGraphics memory setting error\n");
    return -1;
    }

  int curEnd=0;
  int endP=ti.endP;
  ti.endP=ti.startP;

  while (ti.endP<endP)
    {
    //set range
    ti.endP =ti.startP +len -1;
    if (ti.endP>endP)
      ti.endP=endP;
    cu_draw_chunk(&ti);
    for (int x=0; x<pic.size1(); x++)
      for (int y=0; y<pic.size2(); y++)
        (*(pInfoOutput->pPic))[x][y] += pic[x][y];
    for (int i=0; i<TIME_RECORDS; i++)
      pInfoOutput->times[i] +=ti.times[i];

    ti.startP =ti.endP +1;
    }

  return 1;
  }

THREADFUNC cu_draw_chunk(void *pinfo)
  {
  wallTimer timer, timer1;
  float time;

  timer1.reset(); //for the whole thread
  timer1.start();

  //get the input info
  thread_info *tInfo = (thread_info*)pinfo;
  int nParticle =tInfo->endP -tInfo->startP +1;
  //prepare for recording times
  memset(tInfo->times, 0, sizeof(float)*TIME_RECORDS);

  paramfile params =*g_params;

  cu_particle_sim *d_particle_data =new cu_particle_sim[nParticle];

  //for each gpu/thread a varible pack is needed
  cu_gpu_vars gv;
  memset(&gv, 0, sizeof(cu_gpu_vars));

  //CUDA Init
  timer.reset();
  timer.start();
  cu_init(params, tInfo->devID, &gv);
  timer.stop();
  time =timer.acc();
  tInfo->times[CUDA_INIT]=time;

  //Copy particle sim into C-style object d_particle_data
  timer.reset();
  timer.start();
  //copy data to local C-like array d_particle_data, in mid of developing only
  for (int i=tInfo->startP,j=0; i<=tInfo->endP; i++, j++)
    memcpy( &(d_particle_data[j]), &(particle_data[i]), sizeof(cu_particle_sim));
  timer.stop();
  time =timer.acc();
  tInfo->times[COPY2C_LIKE]=time;

  //here we analyse how to divide the whole task for large data handling

  //CUDA Ranging
  timer.reset();
  timer.start();

  cu_range(params, d_particle_data, nParticle, &gv);
  timer.stop();
  time =timer.acc();
  tInfo->times[RANGE]=time;

  //CUDA Transformation
  timer.reset();
  timer.start();
  double  c[3]={campos.x, campos.y, campos.z},
          l[3]={lookat.x, lookat.y, lookat.z},
          s[3]={sky.x, sky.y,     sky.z};
  cu_transform(params, nParticle,c, l, s,d_particle_data, &gv);
  timer.stop();
  time =timer.acc();
  tInfo->times[TRANSFORMATION]=time;

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
  timer.reset();
  timer.start();
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
  timer.stop();
  time =timer.acc();
  tInfo->times[COLORIZE]=time;


//////////////////////////////////////////////////////////////////////
  timer.reset();
  timer.start();
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

  timer.stop();
  time =timer.acc();
  tInfo->times[FILTER] =time;

  timer.reset();
  timer.start();

  //split large ones
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
  timer.stop();
  time =timer.acc();
  tInfo->times[FILTER] +=time;

  timer.reset();
  timer.start();
  v_ps.clear();//not useful any more
  //sort the filtered,splitted v_ps
  // sort(v_ps1.begin(), v_ps1.end(), region_cmp());
  timer.stop();
  time =timer.acc();
  tInfo->times[SORT]=time;

  timer.reset();
  timer.start();
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
  timer.stop();
  time =timer.acc();
  tInfo->times[FILTER] +=time;

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
  timer.reset();//for rendering
  timer.start();
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
  bool bFinished=false;
  int renderStartP, renderEndP, combineStartP, combineEndP;
  renderStartP =renderEndP =0;
  combineStartP = combineEndP=0;
  int nFragments2Render=0;
  //some timers
  wallTimer t1, t2, t3;
  t1.reset();
  t2.reset();

  //now prepare the parameters for combination
  param_combine_thread param;
  param.pPic =tInfo->pPic;
  param.a_eq_e =a_eq_e;
  if (a_eq_e)
    param.fbuf =(void*)fragBufAeqE;
  else
    param.fbuf =(void*)fragBufAneqE;
  param.ps =cu_ps_filtered;
  param.timeUsed =0.0;

  while (!bFinished)
    {
    //find a chunk
    nFragments2Render =0;
    for (renderStartP =renderEndP; renderEndP<pFiltered; renderEndP++)
      {
      int inc;//increasment
      inc =(cu_ps_filtered[renderEndP].maxx-cu_ps_filtered[renderEndP].minx) *
           (cu_ps_filtered[renderEndP].maxy -cu_ps_filtered[renderEndP].miny);
      if ( nFragments2Render +inc > nFBufInCell)
        break;
      else
        nFragments2Render +=inc;
      }
    if (renderEndP == pFiltered)
      bFinished =true;

    //render it
    cu_render1(renderStartP, renderEndP, a_eq_e, grayabsorb, &gv);

    //see if it's the first chunk
    if ( renderStartP!=0 )
      {
      param.combineStartP =combineStartP;
      param.combineEndP =combineEndP;
      //DWORD   id;
      combine(&param);
      // CreateThread( NULL, 0, (LPTHREAD_START_ROUTINE)combine, &param, 0, &id );
      }

    //collect result
    cu_get_fbuf(fragBufAeqE, a_eq_e, nFragments2Render, &gv);
    combineStartP=renderStartP;
    combineEndP =renderEndP;

    //see if last chunk
    if (bFinished)
      {
      param.combineStartP =combineStartP;
      param.combineEndP =combineEndP;
      combine(&param);
      // param.timeUsed +=t3.acc();
      }
    }

  timer.stop();
  time =timer.acc();
  tInfo->times[RENDER]=time;
  tInfo->times[COMBINE]=param.timeUsed;

/////////////////////////////////////////////////////////////////////

  cu_end(&gv);
  delete []d_particle_data;
  delete  []amapD;
  delete  []amapDTypeStartPos;
  delete  []cu_ps;
  delete  []cu_ps_filtered;
  if (a_eq_e)
    delete  []fragBufAeqE;
  else
    delete  []fragBufAneqE;

  tInfo->times[THIS_THREAD] =timer1.acc();

  return 1;
  }

void DevideThreadsTasks(thread_info *tInfo, int nThread, bool bHostThread)
  {
  bool bTestLoadBalancing=g_params->find<bool>("test_load_balancing", false);
  int curStart =0;
  int hostLoad =bHostThread? g_params->find<int>("host_load",0): 0;
  int nDev =bHostThread? nThread-1: nThread;
  int onePercent =particle_data.size()/100;
  int averageDevLen = (nDev!=0)? onePercent *(100-hostLoad)/nDev : 0;

  for (int i=0; i<nThread; i++)
    {
    tInfo[i].startP =curStart;
    if (tInfo[i].devID != -1) //not a host
      {
      if (bTestLoadBalancing)
        {
        int gpuLoad=g_params->find<int>("gpu_load"+dataToString(i),0);
        tInfo[i].endP =curStart +gpuLoad* onePercent;
        }
      else
        tInfo[i].endP =curStart +averageDevLen;
      }
    else //if this is a host
      {
      tInfo[i].endP =curStart +hostLoad *onePercent;
      }
    curStart =tInfo[i].endP +1;
    }

  tInfo[nThread-1].endP =particle_data.size()-1;
  }

void render_cuda(paramfile &params, int &res, arr2<COLOUR> &pic)
  {
  //prepare the parameters for the cuda thread
  //the final image
  res = params.find<int>("resolution",200);
  pic.alloc(res,res);

  //new arrays of thread_info and HANDLE
  int nDev;
  nDev =params.find<int>("gpu_number",0);
  //see if use host as a working thread
  bool bHostThread=params.find<bool>("use_host_as_thread", false);
  int nThread = bHostThread? nDev+1: nDev;
  //init objects for threads control
  thread_info *tInfo =new thread_info[nThread];
  //fill in thread_info
  tInfo[0].pPic =&pic;
  for (int i=0; i<nDev; i++)
    {
    tInfo[i].devID =i;
    //local var pic is assigned to the first device
    if (i!=0)
      tInfo[i].pPic =new arr2<COLOUR>(res, res);
    }
  //make the last one work for host thread
  if (bHostThread )
    {
    tInfo[nThread-1].devID =-1;
    if (nThread-1 != 0)
      tInfo[nThread-1].pPic =new arr2<COLOUR>(res, res);
    }
  //decide how to divide task by another function
  DevideThreadsTasks(tInfo, nThread, bHostThread);

#ifndef NO_WIN_THREAD //to let it compiled in Linux, just for now, 2 Dec 2009.
  HANDLE *tHandle =new HANDLE[nThread];
  //issue the threads
  for (int i=0; i<nDev; i++)
    tHandle[i] =CreateThread( NULL, 0,
      (LPTHREAD_START_ROUTINE)cu_thread_func,&(tInfo[i]), 0, NULL );
  //issue the host thread too
  if (bHostThread)
    tHandle[nDev] =CreateThread( NULL, 0,
      (LPTHREAD_START_ROUTINE)host_thread_func,&(tInfo[nDev]), 0, NULL );
  WaitForMultipleObjects(nThread, tHandle, true, INFINITE);

#else //do not use thread which is now Windows code
  planck_assert(nThread==1,
    "can't have multiple threads on Linux (yet), so 'gpu_number' must be 1");
  cu_thread_func (&(tInfo[0])); //just call it as normal function
  //host_thread_func (&(tInfo[0]));
#endif  //if not NO_WIN_THREAD

  // post-process
  wallTimer  timer;
  timer.reset();
  timer.start();
  // combine the results to pic
  if (1)//a_eq_e)
    {
    for (int i=1; i<nThread; i++)
      for (int x=0; x<res; x++) //  error when x =1,
        for (int y=0; y<res; y++)
          pic[x][y] = pic[x][y] + (*tInfo[i].pPic)[x][y];

    }
  else
    {} //to be done later...

  if (1)//a_eq_e)
    {
    exptable xexp(MAX_EXP);
    for (int ix =0; ix <pic.size1(); ix++)
      for (int iy =0; iy <pic.size2(); iy++)
        {
        pic[ix][iy].r=1-xexp(pic[ix][iy].r);
        pic[ix][iy].g=1-xexp(pic[ix][iy].g);
        pic[ix][iy].b=1-xexp(pic[ix][iy].b);
        }
    }
  else
    {
    }
  timer.stop();

  for (int i=0; i<nThread; i++)
    {
    if (tInfo[i].devID!=-1)
      {
      cout<< endl <<"Times of GPU" << i << ":" <<endl;
      cout<< "CUDA_INIT:              " << tInfo[i].times[CUDA_INIT] <<endl;
      cout<< "COPY2C_LIKE:            " << tInfo[i].times[COPY2C_LIKE] <<endl;
      cout<< "RANGE:                  " << tInfo[i].times[RANGE] <<endl;
      cout<< "TRANSFORMATION:         " << tInfo[i].times[TRANSFORMATION] <<endl;
      cout<< "COLORIZE:               " << tInfo[i].times[COLORIZE] <<endl;
      cout<< "FILTER:                 " << tInfo[i].times[FILTER] <<endl;
      cout<< "SORT:                   " << tInfo[i].times[SORT] <<endl;
      cout<< "RENDER:                 " << tInfo[i].times[RENDER] <<endl;
      cout<< "COMBINE:                " << tInfo[i].times[COMBINE] <<endl;
      cout<< "THIS_THREAD:            " << tInfo[i].times[THIS_THREAD] <<endl;
      cout<<endl;
      }
    else
      {
      cout<< endl <<"Times of CPU as a thread:" <<endl;
      cout<< "RANGE:                  " << tInfo[i].times[RANGE] <<endl;
      cout<< "TRANSFORMATION:         " << tInfo[i].times[TRANSFORMATION] <<endl;
      cout<< "COLORIZE:               " << tInfo[i].times[COLORIZE] <<endl;
      cout<< "RENDER:                 " << tInfo[i].times[RENDER] <<endl;
      cout<< "THIS_THREAD:            " << tInfo[i].times[THIS_THREAD] <<endl;
      cout<<endl;
      }
    }

  for (int i=1; i<nThread; i++)
    delete tInfo[i].pPic;
  delete [] tInfo;
#ifndef NO_WIN_THREAD
  delete [] tHandle;
#endif
  }
