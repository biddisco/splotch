#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef VS
#include "cxxsupport/arr.h"
#include "cxxsupport/cxxutils.h"
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/paramfile.h"
#else
#include "arr.h"
#include "cxxutils.h"
#include "mpi_support.h"
#include "paramfile.h"
#endif
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"

using namespace std;
using namespace RAYPP;

#include "splotch/splotchutils.h"

#define TAG_POSX        11
#define TAG_POSY        12
#define TAG_POSZ        13
#define TAG_SIZE        14
#define TAG_INT         15
#define TAG_COL1        16
#define TAG_COL2        17
#define TAG_COL3        18
#define TAG_TYPE        19
#define TAG_ID          20


void gadget_plain_read_header(bifstream &file, int *npart, double *time)
{
  double massarr[6],redshift;
  int npartall[6];

  file.clear();
  file.rewind();

  file.skip(4);
  file.get(npart,6);
  file.get(massarr,6);
  file >> *time >> redshift;
  file.skip(8);
  file.get(npartall,6);

  file.skip(256 - 6*4 - 6*8 - 3*8 - 6*4);
  file.skip(4);

  }

void gadget_millenium_reader(paramfile &params, vector<particle_sim> &p, int snr, double *time)
{
  int numfiles = params.find<int>("numfiles",1);
  bool doswap = params.find<bool>("swap_endian",true);
  string infilename = params.find<string>("infile");
  int readparallel = params.find<int>("readparallel",1);
  int ptypes = params.find<int>("ptypes",1);

  string filename;
  bifstream infile;

  int ThisTask=mpiMgr.rank(),NTasks=mpiMgr.num_ranks();

  vector<int>	ThisTaskReads,DataFromTask;
  vector<long>	NPartThisTask;
  ThisTaskReads.resize(NTasks);
  DataFromTask.resize(NTasks);
  NPartThisTask.resize(NTasks);

#ifdef USE_MPI
  MPI_Status status;
#endif

  if(mpiMgr.master())
    {
      planck_assert(numfiles >= readparallel,
		    "Number of files must be larger or equal number of parallel reads ...");
      planck_assert(numfiles%readparallel == 0,
		    "Number of files must be a multiple of number of parallel reads ...");
      planck_assert(NTasks >= readparallel,
		    "Number of tasks must be larger or equal number of parallel reads ...");
      planck_assert(NTasks%readparallel == 0,
		    "Number of tasks must be a multiple of number of parallel reads ...");
    }

  // Figure out who will read what 
  int NTaskPerRead=NTasks/readparallel;
  int NFilePerRead=numfiles/readparallel;
  int SendingTask=-1;

  for(int i=0;i<NTasks;i++)
    {
      if(i%NTaskPerRead == 0 && (readparallel > 1 || i == 0))
	{
	  ThisTaskReads[i]=i/NTaskPerRead*NFilePerRead;
	  DataFromTask[i]=-1;
	  SendingTask=i;
	}
      else
	{
	  ThisTaskReads[i]=-1;
	  DataFromTask[i]=SendingTask;
	}
    }

  if(mpiMgr.master())
    {
      int itask=0;
      for(int rt=0;rt<readparallel;rt++)
	{
	  long NPartThisReadTask = 0;
	  for(int f=0;f<NFilePerRead;f++)
	    {
	      int file=rt*NFilePerRead+f;
	      int npartthis[6];
	      if(numfiles>1) filename=infilename+"."+dataToString(file);
	      else           filename=infilename;
	      infile.open(filename.c_str(),doswap);
	      planck_assert (infile,"could not open input file! <" + filename + ">");
	      gadget_plain_read_header(infile,npartthis,time);
	      infile.close();
              cout << " Timestamp from file : " << *time << endl;
	      for(int itype=0;itype<ptypes;itype++)
		{
		  int type = params.find<int>("ptype"+dataToString(itype),0);
                  NPartThisReadTask += npartthis[type];
		}
	    }
	  long SumPartThisReadTask = 0;
	  for(int t=0;t<NTaskPerRead-1;t++)
	    {
	      NPartThisTask[itask] = NPartThisReadTask / NTaskPerRead;
	      SumPartThisReadTask += NPartThisReadTask / NTaskPerRead;
	      itask++;
	    }
	  NPartThisTask[itask] = NPartThisReadTask - SumPartThisReadTask;
	  itask++;
	}
    }
#ifdef USE_MPI
  MPI_Bcast(&NPartThisTask, NTasks, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

  if(mpiMgr.master())
    {
      cout << " Reading " << numfiles << " files by " << readparallel << " tasks ... " << endl;
      
      cout << " Task " << ThisTask << "/" << NTasks << endl;

      cout << " NTaskPerRead/NFilePerRead " << NTaskPerRead << "," << NFilePerRead << endl;

      cout << " ThisTaskReads";
      for(int i=0;i<NTasks;i++)
	cout << ',' << ThisTaskReads[i];
      cout << endl;
      
      cout << " DataFromTask";
      for(int i=0;i<NTasks;i++)
	cout << ',' << DataFromTask[i];
      cout << endl;

      cout << " NPartThis";
      for(int i=0;i<NTasks;i++)
	cout  << ',' << NPartThisTask[i];
      cout << endl;
    }

   //   planck_assert(mpiMgr.num_ranks()==1,
   //		 "sorry, reading Gadget files is not yet MPI parellelized ...");

  long npart=NPartThisTask[ThisTask],nmax=0;

  p.resize(npart);

  for(int i=0;i<NTasks;i++)
    if(NPartThisTask[i] > nmax)
      nmax = NPartThisTask[i];

  float *v1_tmp,*v2_tmp,*v3_tmp;
  int *i1_tmp;
  v1_tmp=new float[nmax];
  v2_tmp=new float[nmax];
  v3_tmp=new float[nmax];
  i1_tmp=new int[nmax];

  if(mpiMgr.master())
    cout << " Reading positions ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
      int ToTask=ThisTask;
      int NPartThis=NPartThisTask[ThisTask];
      long ncount=0;

      for(int f=0;f<NFilePerRead;f++)
	{
	  int npartthis[6];
	  int present=1+2+4+8+16+32;
	  int LastType=0;
	  if(numfiles>1) filename=infilename+"."+dataToString(ThisTaskReads[ThisTask]+f);
	  else           filename=infilename;
          cout << " Task: " << ThisTask << " reading file " << filename << endl;
	  infile.open(filename.c_str(),doswap);
	  planck_assert (infile,"could not open input file! <" + filename + ">");
	  gadget_plain_read_header(infile,npartthis,time);
	  infile.skip(4);
	  for(int itype=0;itype<ptypes;itype++)
	    {
	      int type = params.find<int>("ptype"+dataToString(itype),1);
	      for(int s=LastType+1; s<type; s++)
		if(npartthis[s]>0 && (1<<s & present))
		  {
		    infile.skip(4*3*npartthis[s]);
		  }
	      for(int m=0; m<npartthis[type]; ++m)
		{
		  if(ThisTask == ToTask)
		    {
		      infile >> p[ncount].x >> p[ncount].y >> p[ncount].z;
		      p[ncount].type=itype;
		      ncount++;
		      if(ncount == NPartThisTask[ToTask])
			{
			  ToTask++;
			  ncount=0;
			}
		    }
		  else
		    {
#ifdef USE_MPI
		      infile >> v1_tmp[ncount] >> v2_tmp[ncount] >> v3_tmp[ncount];
		      i1_tmp[ncount] = itype;
		      ncount++;
		      if(ncount == NPartThisTask[ToTask])
			{
			  MPI_Ssend(v1_tmp, NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSX, MPI_COMM_WORLD); 
			  MPI_Ssend(v2_tmp, NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSY, MPI_COMM_WORLD); 
			  MPI_Ssend(v3_tmp, NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSZ, MPI_COMM_WORLD);
			  MPI_Ssend(i1_tmp, NPartThisTask[ToTask], MPI_INT, ToTask, TAG_TYPE, MPI_COMM_WORLD);
			  ToTask++;
			  ncount=0;
			}
#else
		      planck_assert(false,"Should not be executed without MPI support !!!");
#endif
		    }
		}
	      LastType=type;
	    }
	  infile.close();
	}
      planck_assert(ncount == 0,"Some Particles where left when reading Positions ...");
    }
  else
    {
#ifdef USE_MPI
      MPI_Recv(v1_tmp, NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSX, MPI_COMM_WORLD, &status);
      MPI_Recv(v2_tmp, NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSY, MPI_COMM_WORLD, &status);
      MPI_Recv(v3_tmp, NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSZ, MPI_COMM_WORLD, &status);
      MPI_Recv(i1_tmp, NPartThisTask[ThisTask], MPI_INT, DataFromTask[ThisTask], TAG_TYPE, MPI_COMM_WORLD, &status);
      for (int m=0; m<NPartThisTask[ThisTask]; ++m)
	{
	  p[m].x=v1_tmp[m];
	  p[m].y=v2_tmp[m];
	  p[m].z=v3_tmp[m];
	  p[m].type=i1_tmp[m];
	}
#else
      planck_assert(false,"Should not be executed without MPI support !!!");
#endif
    }

  cout << "   -> " << p[0].x << "," << p[0].y << "," << p[0].z << endl;
  cout << "   -> " << p[npart-1].x << "," << p[npart-1].y << "," << p[npart-1].z << endl;

  if(mpiMgr.master())
    cout << " Reading smoothing ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
      int ToTask=ThisTask;
      int NPartThis=NPartThisTask[ThisTask];
      long ncount=0;

      for(int f=0;f<NFilePerRead;f++)
	{
          int npartthis[6];
	  int LastType=0;
	  if(numfiles>1) filename=infilename+"."+dataToString(ThisTaskReads[ThisTask]+f);
	  else           filename=infilename;
	  infile.open(filename.c_str(),doswap);
	  planck_assert (infile,"could not open input file! <" + filename + ">");
	  gadget_plain_read_header(infile,npartthis,time);

	  for(int itype=0;itype<ptypes;itype++)
	    {
	      int type = params.find<int>("ptype"+dataToString(itype),1);
	      float fix_size = params.find<float>("size_fix"+dataToString(itype),1.0);
	      float size_fac = params.find<float>("size_fac"+dataToString(itype),1.0);
	      string label_size = params.find<string>("size_label"+dataToString(itype),"XXXX");
	      if (fix_size == 0.0)
		{
                  planck_assert(false,"Should not be executed for millenium format !!!");
		  // gadget_find_block(infile,label_size);
		  infile.skip(4);
		  int present = params.find<int>("size_present"+dataToString(itype),type);
		  for(int s=LastType+1; s<type; s++)
		    if(npartthis[s]>0 && (1<<s & present))
		      {
			infile.skip(4*npartthis[s]);
		      }
		}
	      for (int m=0; m<npartthis[type]; ++m)
		{
		  if(ThisTask == ToTask)
		    {
		      if (fix_size == 0.0)
			{
			  infile >> p[ncount].r;
			  p[ncount].r *= size_fac;
			}
		      else
			p[ncount].r = fix_size;
		      ncount++;
		      if(ncount == NPartThisTask[ToTask])
			{
			  ToTask++;
			  ncount=0;
			}
		    }
		  else
		    {
#ifdef USE_MPI
		      if (fix_size == 0.0)
			{
			  infile >> v1_tmp[ncount];
			  v1_tmp[ncount] *= size_fac;
			}
		      else
			v1_tmp[ncount] = fix_size;
		      ncount++;
		      if(ncount == NPartThisTask[ToTask])
			{
			  MPI_Ssend(v1_tmp, NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_SIZE, MPI_COMM_WORLD); 
			  ToTask++;
			  ncount=0;
			}
#else
		      planck_assert(false,"Should not be executed without MPI support !!!");
#endif
		    }
		}
	      LastType=type;
	    }
	  infile.close();
	}
      planck_assert(ncount == 0,"Some Particles where left when reading Sizes ...");
    }
  else
    {
#ifdef USE_MPI
      MPI_Recv(v1_tmp, NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_SIZE, MPI_COMM_WORLD, &status);
      for (int m=0; m<NPartThisTask[ThisTask]; ++m)
	p[m].r=v1_tmp[m];
#else
      planck_assert(false,"Should not be executed without MPI support !!!");
#endif
    }

  cout << "   -> " << p[0].r << endl;
  cout << "   -> " << p[npart-1].r << endl;


  if(mpiMgr.master())
    cout << " Reading colors ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
      int ToTask=ThisTask;
      int NPartThis=NPartThisTask[ThisTask];
      long ncount=0;

      for(int f=0;f<NFilePerRead;f++)
	{
          int npartthis[6];
	  int LastType=0;
	  if(numfiles>1) filename=infilename+"."+dataToString(ThisTaskReads[ThisTask]+f);
	  else           filename=infilename;
	  infile.open(filename.c_str(),doswap);
	  planck_assert (infile,"could not open input file! <" + filename + ">");
	  gadget_plain_read_header(infile,npartthis,time);

	  for(int itype=0;itype<ptypes;itype++)
	    {
	      int type = params.find<int>("ptype"+dataToString(itype),1);
	      bool col_vector = params.find<bool>("color_is_vector"+dataToString(itype),true);
	      float col_fac = params.find<float>("color_fac"+dataToString(itype),1.0);
	      int read_col = 1;
              infile.skip(4);
	      infile.skip(npartthis[type]*3*4);
	      infile.skip(4);

	      infile.skip(4);
	      int present = 1+2+4+8+16+32;
	      for(int s=LastType+1; s<type; s++)
		if(npartthis[s]>0 && (1<<s & present))
		  {
		    int nskip=npartthis[s];
		    if(col_vector)
		      nskip *=3;
		    infile.skip(4*nskip);
		  }
	      for (int m=0; m<npartthis[type]; ++m)
		{
		  if(ThisTask == ToTask)
		    {
		      if (read_col > 0)
			{
			  infile >> p[ncount].C1;
			  p[ncount].C1 *= col_fac;
			  if(col_vector)
			    {
			      infile >> p[ncount].C2 >> p[ncount].C3;
			      p[ncount].C2 *= col_fac;
			      p[ncount].C3 *= col_fac;
                              p[ncount].I = 1;
			    }
			}
		      else
			{
			  p[ncount].C1 = 1;
			  p[ncount].C2 = 1;
			  p[ncount].C3 = 1;
                          p[ncount].I = 1;
			}
		      ncount++;
		      if(ncount == NPartThisTask[ToTask])
			{
			  ToTask++;
			  ncount=0;
			}
		    }
		  else
		    {
#ifdef USE_MPI
		      if (read_col > 0)
			{
			  infile >> v1_tmp[ncount];
			  v1_tmp[ncount] *= col_fac;
			  if(col_vector)
			    {
			      infile >> v2_tmp[ncount] >> v3_tmp[ncount];
			      v2_tmp[ncount] *= col_fac;
			      v3_tmp[ncount] *= col_fac;
			    }
			}
		      else
			{
			  v1_tmp[ncount] = 1;
			  v2_tmp[ncount] = 1;
			  v3_tmp[ncount] = 1;
			}
		      ncount++;
		      if(ncount == NPartThisTask[ToTask])
			{
			  MPI_Ssend(v1_tmp, NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_COL1, MPI_COMM_WORLD); 
			  MPI_Ssend(v2_tmp, NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_COL2, MPI_COMM_WORLD); 
			  MPI_Ssend(v3_tmp, NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_COL3, MPI_COMM_WORLD); 
			  ToTask++;
			  ncount=0;
			}
#else
		      planck_assert(false,"Should not be executed without MPI support !!!");
#endif
		    }
		}
	      LastType=type;
	    }
	  infile.close();
	}
      planck_assert(ncount == 0,"Some Particles where left when reading Colors ...");
    }
  else
    {
#ifdef USE_MPI
      MPI_Recv(v1_tmp, NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_COL1, MPI_COMM_WORLD, &status);
      MPI_Recv(v2_tmp, NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_COL2, MPI_COMM_WORLD, &status);
      MPI_Recv(v3_tmp, NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_COL3, MPI_COMM_WORLD, &status);
      for (int m=0; m<NPartThisTask[ThisTask]; ++m)
	{
	  p[m].C1=v1_tmp[m];
	  p[m].C2=v2_tmp[m];
	  p[m].C3=v3_tmp[m];
          p[m].I=1.0;
	}
#else
      planck_assert(false,"Should not be executed without MPI support !!!");
#endif
    }

  cout << "   -> " << p[0].C1 << "," << p[0].C2 << "," << p[0].C3 << "," << p[0].I << endl;
  cout << "   -> " << p[npart-1].C1 << "," << p[npart-1].C2 << "," << p[npart-1].C3 << "," << p[npart-1].I << endl;

}
