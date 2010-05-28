#ifdef USE_MPI
#include "mpi.h"
#endif
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "cxxsupport/arr.h"
#include "cxxsupport/cxxutils.h"
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/paramfile.h"
#include "kernel/bstream.h"
#include "splotch/splotchutils.h"

using namespace std;

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

#ifdef INTERPOLATE

struct idcmp
  {
  bool operator()(const particle_sim &p1, const particle_sim &p2)
    { return p1.id<p2.id; }
  };

#endif

int gadget_find_block (bifstream &file,const string &label)
  {
  int i;
  int32 blocksize=0, blksize;
  char blocklabel[5]={"    "};

  file.clear();
  file.rewind();

  while(!file.eof() && blocksize == 0)
    {
    file >> blksize;
    if(file.eof())
      {
      blksize=-1;
      break;
      }
    planck_assert(blksize==8,
      "wrong structure in GADGET file: " +dataToString(blksize));

    file.get(blocklabel,4);
    blocklabel[4] = '\0';
    i=3;
    while (blocklabel[i]==' ')
      {
      blocklabel[i]='\0';
      i--;
      }
    file >> blocksize;
    file.skip(4);
    if (label!=blocklabel)
      {
      file.skip(blocksize);
      blocksize=0;
      }
    }
  if(file.eof()) file.clear();
  return(blocksize-8);
  }

int gadget_read_header(bifstream &file, int32 *npart, double *time)
  {
  int blocksize = gadget_find_block (file,"HEAD");
  planck_assert (blocksize>0, "Header block not found");
  file.skip(4);
  file.get(npart,6);
  file.skip(6*8);
  file >> *time;
  file.skip(8+8+6*4);

  return blocksize;
  }

#ifdef INTERPOLATE
void gadget_reader(paramfile &params, vector<particle_sim> &p, int snr, double *time)
#else
void gadget_reader(paramfile &params, vector<particle_sim> &p, int /*snr*/, double *time)
#endif
  {
  int numfiles = params.find<int>("numfiles",1);
  bool doswap = params.find<bool>("swap_endian",true);
  string infilename = params.find<string>("infile");
  int readparallel = params.find<int>("readparallel",1);
  int ptypes = params.find<int>("ptypes",1);

  string filename;
  bifstream infile;

  int ThisTask=mpiMgr.rank(),NTasks=mpiMgr.num_ranks();
  arr<int> ThisTaskReads(NTasks), DataFromTask(NTasks);
  arr<long> NPartThisTask(NTasks);

#ifdef USE_MPI
  MPI_Status status;
#endif

#ifdef INTERPOLATE
  infilename += intToString(snr,3);
#endif

  planck_assert(numfiles >= readparallel,
                "Number of files must be larger or equal number of parallel reads ...");
  planck_assert(numfiles%readparallel == 0,
                "Number of files must be a multiple of number of parallel reads ...");
  planck_assert(NTasks >= readparallel,
                "Number of tasks must be larger or equal number of parallel reads ...");
  planck_assert(NTasks%readparallel == 0,
                "Number of tasks must be a multiple of number of parallel reads ...");

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
        filename=infilename;
        if(numfiles>1) filename+="."+dataToString(file);
        infile.open(filename.c_str(),doswap);
        planck_assert (infile,"could not open input file! <" + filename + ">");
        gadget_read_header(infile,npartthis,time);
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

  mpiMgr.bcastRaw(&NPartThisTask[0],NTasks,0);

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

  long npart=NPartThisTask[ThisTask],nmax=0;
  p.resize(npart);

  for(int i=0;i<NTasks;i++)
    if(NPartThisTask[i] > nmax)
      nmax = NPartThisTask[i];

  arr<float> v1_tmp(nmax), v2_tmp(nmax), v3_tmp(nmax);
  arr<int> i1_tmp(nmax);

  if(mpiMgr.master())
    cout << " Reading positions ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
    int ToTask=ThisTask;
    long ncount=0;

    for(int f=0;f<NFilePerRead;f++)
      {
      int npartthis[6];
      int present=1+2+4+8+16+32;
      int LastType=-1;
      filename=infilename;
      if(numfiles>1) filename+="."+dataToString(ThisTaskReads[ThisTask]+f);
      cout << " Task: " << ThisTask << " reading file " << filename << endl;
      infile.open(filename.c_str(),doswap);
      planck_assert (infile,"could not open input file! <" + filename + ">");
      gadget_read_header(infile,npartthis,time);
      gadget_find_block(infile,"POS");
      infile.skip(4);
      for(int itype=0;itype<ptypes;itype++)
        {
        int type = params.find<int>("ptype"+dataToString(itype),0);
        for(int s=LastType+1; s<type; s++)
          if(npartthis[s]>0 && (1<<s & present))
            infile.skip(4*3*npartthis[s]);
        arr<float32> ftmp(3*npartthis[type]);
        infile.get(&ftmp[0],ftmp.size());
        for(int m=0; m<npartthis[type]; ++m)
          {
          if(ThisTask == ToTask)
            {
            p[ncount].x=ftmp[3*m];
            p[ncount].y=ftmp[3*m+1];
            p[ncount].z=ftmp[3*m+2];
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
            v1_tmp[ncount]=ftmp[3*m];
            v2_tmp[ncount]=ftmp[3*m+1];
            v3_tmp[ncount]=ftmp[3*m+2];
            i1_tmp[ncount] = itype;
            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              MPI_Ssend(&v1_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSX, MPI_COMM_WORLD);
              MPI_Ssend(&v2_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSY, MPI_COMM_WORLD);
              MPI_Ssend(&v3_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSZ, MPI_COMM_WORLD);
              MPI_Ssend(&i1_tmp[0], NPartThisTask[ToTask], MPI_INT, ToTask, TAG_TYPE, MPI_COMM_WORLD);
              ToTask++;
              ncount=0;
              }
#else
            planck_fail("Should not be executed without MPI support !!!");
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
    MPI_Recv(&v1_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSX, MPI_COMM_WORLD, &status);
    MPI_Recv(&v2_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSY, MPI_COMM_WORLD, &status);
    MPI_Recv(&v3_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSZ, MPI_COMM_WORLD, &status);
    MPI_Recv(&i1_tmp[0], NPartThisTask[ThisTask], MPI_INT, DataFromTask[ThisTask], TAG_TYPE, MPI_COMM_WORLD, &status);
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

#ifdef INTERPOLATE
#ifdef HIGH_ORDER_INTERPOLATION
  if(mpiMgr.master())
    cout << " Reading velocities ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
    int ToTask=ThisTask;
    long ncount=0;

    for(int f=0;f<NFilePerRead;f++)
      {
      int npartthis[6];
      int present=1+2+4+8+16+32;
      int LastType=-1;
      filename=infilename;
      if(numfiles>1) filename+="."+dataToString(ThisTaskReads[ThisTask]+f);
      infile.open(filename.c_str(),doswap);
      planck_assert (infile,"could not open input file! <" + filename + ">");
      gadget_read_header(infile,npartthis,time);
      gadget_find_block(infile,"VEL");
      infile.skip(4);
      for(int itype=0;itype<ptypes;itype++)
        {
        int type = params.find<int>("ptype"+dataToString(itype),0);
        for(int s=LastType+1; s<type; s++)
          if(npartthis[s]>0 && (1<<s & present))
            infile.skip(4*3*npartthis[s]);
        arr<float32> ftmp(3*npartthis[type]);
        infile.get(&ftmp[0],ftmp.size());
        for(int m=0; m<npartthis[type]; ++m)
          {
          if(ThisTask == ToTask)
            {
            p[ncount].vx=ftmp[3*m];
            p[ncount].vy=ftmp[3*m+1];
            p[ncount].vz=ftmp[3*m+2];
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
            v1_tmp[ncount]=ftmp[3*m];
            v2_tmp[ncount]=ftmp[3*m+1];
            v3_tmp[ncount]=ftmp[3*m+2];
            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              MPI_Ssend(&v1_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSX, MPI_COMM_WORLD);
              MPI_Ssend(&v2_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSY, MPI_COMM_WORLD);
              MPI_Ssend(&v3_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_POSZ, MPI_COMM_WORLD);
              ToTask++;
              ncount=0;
              }
#else
            planck_fail("Should not be executed without MPI support !!!");
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
    MPI_Recv(&v1_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSX, MPI_COMM_WORLD, &status);
    MPI_Recv(&v2_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSY, MPI_COMM_WORLD, &status);
    MPI_Recv(&v3_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_POSZ, MPI_COMM_WORLD, &status);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      {
      p[m].vx=v1_tmp[m];
      p[m].vy=v2_tmp[m];
      p[m].vz=v3_tmp[m];
      }
#else
    planck_fail("Should not be executed without MPI support !!!");
#endif
    }
#endif

  if(mpiMgr.master())
    cout << " Reading ids ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
    int ToTask=ThisTask;
    long ncount=0;

    for(int f=0;f<NFilePerRead;f++)
      {
      int npartthis[6];
      int present=1+2+4+8+16+32;
      int LastType=-1;
      filename=infilename;
      if(numfiles>1) filename+="."+dataToString(ThisTaskReads[ThisTask]+f);
      infile.open(filename.c_str(),doswap);
      planck_assert (infile,"could not open input file! <" + filename + ">");
      gadget_read_header(infile,npartthis,time);
      string label_id = params.find<string>("id_label","ID");
      gadget_find_block(infile,label_id);
      infile.skip(4);
      for(int itype=0;itype<ptypes;itype++)
        {
        int type = params.find<int>("ptype"+dataToString(itype),0);
        for(int s=LastType+1; s<type; s++)
          if(npartthis[s]>0 && (1<<s & present))
            infile.skip(4*npartthis[s]);
        arr<unsigned int> ftmp(npartthis[type]);
        infile.get(&ftmp[0],ftmp.size());
        for(int m=0; m<npartthis[type]; ++m)
          {
          if(ThisTask == ToTask)
            {
            p[ncount].id=ftmp[m];
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
            i1_tmp[ncount]=ftmp[m];
            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              MPI_Ssend(&i1_tmp[0], NPartThisTask[ToTask], MPI_INT, ToTask, TAG_ID, MPI_COMM_WORLD);
              ToTask++;
              ncount=0;
              }
#else
            planck_fail("Should not be executed without MPI support !!!");
#endif
            }
          }
        LastType=type;
        }
        infile.close();
      }
    planck_assert(ncount == 0,"Some Particles where left when reading IDs ...");
    }
  else
    {
#ifdef USE_MPI
    MPI_Recv(&i1_tmp[0], NPartThisTask[ThisTask], MPI_INT, DataFromTask[ThisTask], TAG_ID, MPI_COMM_WORLD, &status);
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
#endif

  if(mpiMgr.master())
    cout << " Reading smoothing ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
    int ToTask=ThisTask;
    // int NPartThis=NPartThisTask[ThisTask];
    long ncount=0;

    for(int f=0;f<NFilePerRead;f++)
      {
      int npartthis[6];
      int LastType=-1;
      filename=infilename;
      if(numfiles>1) filename+="."+dataToString(ThisTaskReads[ThisTask]+f);
      infile.open(filename.c_str(),doswap);
      planck_assert (infile,"could not open input file! <" + filename + ">");
      gadget_read_header(infile,npartthis,time);

      for(int itype=0;itype<ptypes;itype++)
        {
        int type = params.find<int>("ptype"+dataToString(itype),0);
        float fix_size = params.find<float>("size_fix"+dataToString(itype),1.0);
        float size_fac = params.find<float>("size_fac"+dataToString(itype),1.0);
        string label_size = params.find<string>("size_label"+dataToString(itype),"XXXX");
        if (fix_size == 0.0)
          {
          gadget_find_block(infile,label_size);
          infile.skip(4);
          int present = params.find<int>("size_present"+dataToString(itype),type);
          for(int s=LastType+1; s<type; s++)
            if(npartthis[s]>0 && (1<<s & present))
              infile.skip(4*npartthis[s]);
          }
        arr<float32> ftmp(npartthis[type]);
        if (fix_size == 0.0) infile.get(&ftmp[0],ftmp.size());
        for (int m=0; m<npartthis[type]; ++m)
          {
          if(ThisTask == ToTask)
            {
            if (fix_size == 0.0)
              {
              p[ncount].r=ftmp[m];
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
              v1_tmp[ncount]=ftmp[m];
              v1_tmp[ncount] *= size_fac;
              }
            else
              v1_tmp[ncount] = fix_size;
            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              MPI_Ssend(&v1_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_SIZE, MPI_COMM_WORLD);
              ToTask++;
              ncount=0;
              }
#else
            planck_fail("Should not be executed without MPI support !!!");
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
    MPI_Recv(&v1_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_SIZE, MPI_COMM_WORLD, &status);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      p[m].r=v1_tmp[m];
#else
    planck_assert(false,"Should not be executed without MPI support !!!");
#endif
    }


  if(mpiMgr.master())
    cout << " Reading colors ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
    int ToTask=ThisTask;
    // int NPartThis=NPartThisTask[ThisTask];
    long ncount=0;

    for(int f=0;f<NFilePerRead;f++)
      {
      int npartthis[6];
      int LastType=-1;
      filename=infilename;
      if(numfiles>1) filename+="."+dataToString(ThisTaskReads[ThisTask]+f);
      infile.open(filename.c_str(),doswap);
      planck_assert (infile,"could not open input file! <" + filename + ">");
      gadget_read_header(infile,npartthis,time);

      for(int itype=0;itype<ptypes;itype++)
        {
        int type = params.find<int>("ptype"+dataToString(itype),0);
        string label_col = params.find<string>("color_label"+dataToString(itype),"XXXX");
        bool col_vector = params.find<bool>("color_is_vector"+dataToString(itype),false);
        float col_fac = params.find<float>("color_fac"+dataToString(itype),1.0);
        int read_col=gadget_find_block(infile,label_col);
        if (read_col > 0)
          {
          infile.skip(4);
          int present = params.find<int>("color_present"+dataToString(itype),type);
          for(int s=LastType+1; s<type; s++)
            if(npartthis[s]>0 && (1<<s & present))
              {
              int nskip=npartthis[s];
              if(col_vector)
                nskip *=3;
              infile.skip(4*nskip);
              }
          }
        else
          if(mpiMgr.master())
            cout << " Cannot find color field <" << label_col << "> ..." << endl;
        tsize fnread=0;
        if (read_col>0) fnread=npartthis[type];
        if ((read_col>0) && col_vector) fnread=3*npartthis[type];
        arr<float32> ftmp(fnread);
        infile.get(&ftmp[0],ftmp.size());
        for (int m=0; m<npartthis[type]; ++m)
          {
          if(ThisTask == ToTask)
            {
            if (read_col > 0)
              {
              tsize ofs = col_vector ? 3*m : m;
              p[ncount].C1 = ftmp[ofs];
              p[ncount].C1 *= col_fac;
              if(col_vector)
                {
                p[ncount].C2 = ftmp[ofs+1];
                p[ncount].C3 = ftmp[ofs+2];
                p[ncount].C2 *= col_fac;
                p[ncount].C3 *= col_fac;
                }
              }
            else
              {
              p[ncount].C1 = 1;
              p[ncount].C2 = 1;
              p[ncount].C3 = 1;
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
              tsize ofs = col_vector ? 3*m : m;
              v1_tmp[ncount] = ftmp[ofs];
              v1_tmp[ncount] *= col_fac;
              if(col_vector)
                {
                v2_tmp[ncount] = ftmp[ofs+1];
                v3_tmp[ncount] = ftmp[ofs+2];
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
              MPI_Ssend(&v1_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_COL1, MPI_COMM_WORLD);
              MPI_Ssend(&v2_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_COL2, MPI_COMM_WORLD);
              MPI_Ssend(&v3_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_COL3, MPI_COMM_WORLD);
              ToTask++;
              ncount=0;
              }
#else
            planck_fail("Should not be executed without MPI support !!!");
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
    MPI_Recv(&v1_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_COL1, MPI_COMM_WORLD, &status);
    MPI_Recv(&v2_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_COL2, MPI_COMM_WORLD, &status);
    MPI_Recv(&v3_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_COL3, MPI_COMM_WORLD, &status);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      {
      p[m].C1=v1_tmp[m];
      p[m].C2=v2_tmp[m];
      p[m].C3=v3_tmp[m];
      }
#else
    planck_assert(false,"Should not be executed without MPI support !!!");
#endif
    }


  if(mpiMgr.master())
    cout << " Reading intensity ..." << endl;
  if(ThisTaskReads[ThisTask] >= 0)
    {
    int ToTask=ThisTask;
    long ncount=0;

    for(int f=0;f<NFilePerRead;f++)
      {
      int npartthis[6];
      int LastType=-1;
      filename=infilename;
      if(numfiles>1) filename+="."+dataToString(ThisTaskReads[ThisTask]+f);
      infile.open(filename.c_str(),doswap);
      planck_assert (infile,"could not open input file! <" + filename + ">");
      gadget_read_header(infile,npartthis,time);
      for(int itype=0;itype<ptypes;itype++)
        {
        int type = params.find<int>("ptype"+dataToString(itype),0);
        string label_int = params.find<string>("intensity_label"+dataToString(itype),"XXXX");
        int read_int=gadget_find_block(infile,label_int);
        if (read_int > 0)
          {
          infile.skip(4);
          int present = params.find<int>("intensity_present"+dataToString(itype),type);
          for(int s=LastType+1; s<type; s++)
            if(npartthis[s]>0 && (1<<s & present))
              infile.skip(4*npartthis[s]);
          }
        else
          if(mpiMgr.master() && itype==0 && f==0)
            cout << " Cannot find intensity field <" << label_int << "> ..." << endl;
        arr<float32> ftmp(npartthis[type]);
        if (read_int>0) infile.get(&ftmp[0],ftmp.size());
        for (int m=0; m<npartthis[type]; ++m)
          {
          if(ThisTask == ToTask)
            {
            if (read_int > 0)
              p[ncount].I=ftmp[m];
            else
              p[ncount].I = 1;
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
            if (read_int > 0)
              v1_tmp[ncount]=ftmp[m];
            else
              v1_tmp[ncount] = 1;
            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              MPI_Ssend(&v1_tmp[0], NPartThisTask[ToTask], MPI_FLOAT, ToTask, TAG_INT, MPI_COMM_WORLD);
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
    MPI_Recv(&v1_tmp[0], NPartThisTask[ThisTask], MPI_FLOAT, DataFromTask[ThisTask], TAG_INT, MPI_COMM_WORLD, &status);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      p[m].I=v1_tmp[m];
#else
    planck_assert(false,"Should not be executed without MPI support !!!");
#endif
    }

#ifdef INTERPOLATE
  //   parallel_sort(p, npart, sizeof(struct particle_sim), io_compare_P_ID);

  planck_assert(mpiMgr.num_ranks()==1,
     "sorry, interpolating between files is not yet MPI parellelized ...");
  if(mpiMgr.master())
    cout << " sorting particles by ID ..." << endl;

  sort(p.begin(), p.end(), idcmp());
#endif
  }
