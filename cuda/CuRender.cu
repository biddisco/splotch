#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sort.h>
#include <thrust/unique.h>


#include "cuda/CuRender.h"
#include "cuda/CuPolicy.h"

using namespace std;

void cu_draw_chunk(void *pinfo, COLOUR *Pic, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb)
{
  //get the input info
  thread_info *tInfo = (thread_info*)pinfo;
  tInfo->times.start("gpu_thread");

  int nParticle = tInfo->endP - tInfo->startP + 1;
  printf("Rank %d - GPU %d : Processing %d particles\n", mpiMgr.rank(), tInfo->devID, nParticle); fflush(stdout);

  //allocate and copy data particle to device memory
  tInfo->times.start("gcopy");
  cu_particle_sim *d_particle_data = &((*particle_data)[tInfo->startP]);
  cu_copy_particles_to_device(d_particle_data, nParticle, gv);
  tInfo->times.stop("gcopy");

  //init cu_particle_splotch array memory
  cu_particle_splotch *cu_ps;
  cu_ps = new cu_particle_splotch[nParticle];
  memset(cu_ps, 0, nParticle);

  //CUDA Transformation
  tInfo->times.start("gtransform");
  int err;
  err = cu_transform(nParticle, cu_ps, gv);
  if (err) return;
  tInfo->times.stop("gtransform");
  
  cudaError_t error;
  error = cudaFree(gv->d_pd);
  if (error != cudaSuccess) cout << "Device Free error!" << endl; 

/* temporarily ignore sorting 191109.
// --------------------------------
// ----------- Sorting ------------
// --------------------------------
   it becomes complicated when using multiple threads with sorting
 
*/

// ----------------------------------
// ----------- Rendering ------------
// ----------------------------------

  //get parameters for rendering
  int maxRegion = gv->policy->GetMaxRegion();
  pair<int, int> res = gv->policy->GetResolution();

  //prepare fragment buffer memory space first
  void *fragBuf;
  size_t nFBufInByte = gv->policy->GetFBufSize() <<20;

  int nFBufInCell;
  if (a_eq_e)
    {
    nFBufInCell = nFBufInByte/sizeof(cu_fragment_AeqE);
    fragBuf = new cu_fragment_AeqE[nFBufInCell];
    }
  else
    {
    nFBufInCell = nFBufInByte/sizeof(cu_fragment_AneqE);
    fragBuf = new cu_fragment_AneqE[nFBufInCell];
    }

  // new array of particles produced after filter and splitting
  cu_particle_splotch *cu_ps_filtered;
  cu_ps_filtered = new cu_particle_splotch[nParticle];

  //clear the device image
  int size_Im = res.first * res.second;
  cu_clear(size_Im, gv);

  // wrap raw pointer with a device_ptr 
  thrust::device_ptr<cu_fragment_AeqE> dev_ptr_FragBuf((cu_fragment_AeqE *) gv->d_fbuf);
  thrust::device_ptr<int> dev_ptr_Index(gv->d_pixel);

  int End_cu_ps = 0, Start_cu_ps=0;
  while (End_cu_ps < nParticle)
  {
   int nFragments2Render = 0;
   //filter and split particles to a cu_ps_filtered
   tInfo->times.start("gfilter");
   int pFiltered = filter_chunk(Start_cu_ps, nParticle, maxRegion,
                                nFBufInCell, cu_ps, cu_ps_filtered, &End_cu_ps,
                                &nFragments2Render);
   tInfo->times.stop("gfilter");

   tInfo->times.start("gcopy");
   cu_copy_particles_to_render(cu_ps_filtered, pFiltered, gv);
   tInfo->times.stop("gcopy");

   // colorize and render chunks of pFiltered particles
   tInfo->times.start("grender");
   cu_render1(pFiltered, a_eq_e, (float) grayabsorb, gv);
   tInfo->times.stop("grender");

   tInfo->times.start("gsort");
   thrust::sort_by_key(dev_ptr_Index, dev_ptr_Index+nFragments2Render, dev_ptr_FragBuf);
   tInfo->times.stop("gsort");

   tInfo->times.start("greduce");
   thrust::pair< thrust::device_ptr<int>,thrust::device_ptr<cu_fragment_AeqE> >  new_end;
   thrust::equal_to<int> binary_pred;
   new_end = thrust::reduce_by_key(dev_ptr_Index, dev_ptr_Index+nFragments2Render, dev_ptr_FragBuf, dev_ptr_Index, dev_ptr_FragBuf, binary_pred, sum_op());
   int npixels = new_end.first.get() - dev_ptr_Index.get();
   tInfo->times.stop("greduce");

 //  int npixels = nFragments2Render; 
   tInfo->times.start("gcombine");
   cu_update_image(npixels, a_eq_e, gv);
   tInfo->times.stop("gcombine");

   printf("Rank %d - GPU %d : Rendered %d/%d particles\n",  mpiMgr.rank(), tInfo->devID, End_cu_ps, nParticle); 

   Start_cu_ps = End_cu_ps;
  }

  // copy back the image
  size_t size = size_Im*sizeof(cu_color);
  error = cudaMemcpy(Pic, gv->d_pic, size, cudaMemcpyDeviceToHost);
  if (error != cudaSuccess) cout << "Device Memcpy error!" << endl; 
 
  delete []cu_ps;
  delete []cu_ps_filtered;
  if (a_eq_e)
    delete [] ((cu_fragment_AeqE *) fragBuf);
  else
    delete [] ((cu_fragment_AneqE *) fragBuf);

  tInfo->times.stop("gpu_thread");
}


//filter and split particles to a cu_ps_filtered of size nParticles
int filter_chunk(int StartP, int nParticle, int maxRegion, 
                 int nFBufInCell, cu_particle_splotch *cu_ps, 
                 cu_particle_splotch *cu_ps_filtered, int *End_cu_ps, 
                 int *nFragments2Render)
{
  cu_particle_splotch p, pNew;
  int region, nsplit;
  bool finished = false;
  int chunk_dim = nParticle;

  unsigned long posInFragBuf = 0; 
  int pFiltered = 0;
  int i=StartP;  // start chunk position in cu_ps

  // filter particles until cu_ps is finished or cu_ps_filtered array is full
  while(!finished && (i < nParticle))
   {
     //select valid ones
     p = cu_ps[i];
     if (p.isValid)
     {
       int h = p.maxy - p.miny;
       int w = p.maxx - p.minx;
       region = h*w;

       if (region <= maxRegion)
//       if(p.r <= 32.0)
       {
         // set the start position of the particle in fragment buffer
         if (posInFragBuf+region < nFBufInCell) 
         {
           p.posInFragBuf = posInFragBuf;
           cu_ps_filtered[pFiltered] = p;
           pFiltered++;
           posInFragBuf += region;
         }
         else finished = true; 
       }
       else
       { //particle too big -> split along y direction
         pNew = p;
         int w1 = (maxRegion%h == 0) ? (maxRegion/h):(maxRegion/h + 1);
         nsplit = w/w1 + 1;
         if ((pFiltered + nsplit <= chunk_dim) && (posInFragBuf+region < nFBufInCell))
         {
           for (int minx = p.minx; minx < p.maxx; minx += w1)
           {
             pNew.minx = minx;  //minx,maxx of pNew need to be set
             pNew.maxx = (minx+w1 >= p.maxx) ? p.maxx : minx+w1; 
             // set the start position of the particle in fragment buffer
             pNew.posInFragBuf = posInFragBuf; 
             cu_ps_filtered[pFiltered] = pNew;

             pFiltered++;
             int newregion = (pNew.maxx - pNew.minx) * (pNew.maxy - pNew.miny);
             posInFragBuf += newregion;
           }
         }
         else finished = true; 
       }
      }
      i++;
    }

   *End_cu_ps = i;
   *nFragments2Render = posInFragBuf;
   return pFiltered;  // return chunk position reached in cu_ps
}


/*
void combine_chunk(int StartP, int EndP, cu_particle_splotch *cu_ps_filtered, 
                  void *fragBuf, bool a_eq_e, float64 grayabsorb, arr2<COLOUR> &pPic)
{

    if (a_eq_e)
    {
      cu_fragment_AeqE *fragBufAeqE = (cu_fragment_AeqE *)fragBuf;
      for (int pPos=StartP, fPos=0; pPos<EndP; pPos++)
      {
        for (int x =cu_ps_filtered[pPos].minx; x <cu_ps_filtered[pPos].maxx; x++)
        {
          for (int y =cu_ps_filtered[pPos].miny; y <cu_ps_filtered[pPos].maxy; y++)
          {
            pPic[x][y].r += fragBufAeqE[fPos].aR;
            pPic[x][y].g += fragBufAeqE[fPos].aG;
            pPic[x][y].b += fragBufAeqE[fPos].aB;
            fPos ++;
          }
        }
      }
    }
    else
    {
      cu_fragment_AneqE *fragBufAneqE = (cu_fragment_AneqE *)fragBuf;
      for (int pPos=StartP, fPos=0; pPos<EndP; pPos++)
      {
        for (int x =cu_ps_filtered[pPos].minx; x <cu_ps_filtered[pPos].maxx; x++)
        {
          for (int y =cu_ps_filtered[pPos].miny; y <cu_ps_filtered[pPos].maxy; y++) 
          {
            pPic[x][y].r += fragBufAneqE[fPos].aR *
				   (pPic[x][y].r - fragBufAneqE[fPos].qR);
            pPic[x][y].g += fragBufAneqE[fPos].aG *
				   (pPic[x][y].g - fragBufAneqE[fPos].qG); 
            pPic[x][y].b += fragBufAneqE[fPos].aB *
                                   (pPic[x][y].b - fragBufAneqE[fPos].qB);
            fPos ++;
          }
        }
      }
    }
}
*/

