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

int gadget_read_header(bifstream &file, int32 *npart, double &time)
  {
  int blocksize = gadget_find_block (file,"HEAD");
  planck_assert (blocksize>0, "Header block not found");
  file.skip(4);
  file.get(npart,6);
  file.skip(6*8);
  file >> time;
  file.skip(8+8+6*4);

  return blocksize;
  }

void gadget_reader(paramfile &params, int interpol_mode,
  vector<particle_sim> &p, vector<uint32> &id, vector<vec3f> &vel, int snr,
  double &time)
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

  if (interpol_mode>0)
    infilename += intToString(snr,3);

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
        cout << " Timestamp from file : " << time << endl;
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

  mpiMgr.bcast(NPartThisTask,0);

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
  if (interpol_mode>0)
    id.resize(npart);
  if (interpol_mode>1)
    vel.resize(npart);

  for(int i=0;i<NTasks;i++)
    if(NPartThisTask[i] > nmax)
      nmax = NPartThisTask[i];

  arr<float> v1_tmp(nmax), v2_tmp(nmax), v3_tmp(nmax);
  arr<uint32> i1_tmp(nmax);

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
            v1_tmp[ncount]=ftmp[3*m];
            v2_tmp[ncount]=ftmp[3*m+1];
            v3_tmp[ncount]=ftmp[3*m+2];
            i1_tmp[ncount] = itype;
            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              mpiMgr.sendRaw(&v1_tmp[0], NPartThisTask[ToTask], ToTask);
              mpiMgr.sendRaw(&v2_tmp[0], NPartThisTask[ToTask], ToTask);
              mpiMgr.sendRaw(&v3_tmp[0], NPartThisTask[ToTask], ToTask);
              mpiMgr.sendRaw(&i1_tmp[0], NPartThisTask[ToTask], ToTask);
              ToTask++;
              ncount=0;
              }
            }
          }
        LastType=type;
        }
      infile.close();
      }
    planck_assert(ncount == 0,"Some particles were left when reading positions ...");
    }
  else
    {
    mpiMgr.recvRaw(&v1_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    mpiMgr.recvRaw(&v2_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    mpiMgr.recvRaw(&v3_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    mpiMgr.recvRaw(&i1_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      {
      p[m].x=v1_tmp[m];
      p[m].y=v2_tmp[m];
      p[m].z=v3_tmp[m];
      p[m].type=i1_tmp[m];
      }
    }

  if (interpol_mode>0)
    {
    if (interpol_mode>1)
      {
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
                vel[ncount].x=ftmp[3*m];
                vel[ncount].y=ftmp[3*m+1];
                vel[ncount].z=ftmp[3*m+2];
                ncount++;
                if(ncount == NPartThisTask[ToTask])
                  {
                  ToTask++;
                  ncount=0;
                  }
                }
              else
                {
                v1_tmp[ncount]=ftmp[3*m];
                v2_tmp[ncount]=ftmp[3*m+1];
                v3_tmp[ncount]=ftmp[3*m+2];
                ncount++;
                if(ncount == NPartThisTask[ToTask])
                  {
                  mpiMgr.sendRaw(&v1_tmp[0], NPartThisTask[ToTask], ToTask);
                  mpiMgr.sendRaw(&v2_tmp[0], NPartThisTask[ToTask], ToTask);
                  mpiMgr.sendRaw(&v3_tmp[0], NPartThisTask[ToTask], ToTask);
                  ToTask++;
                  ncount=0;
                  }
                }
              }
            LastType=type;
            }
          infile.close();
          }
        planck_assert(ncount == 0,"Some particles were left when reading positions ...");
        }
      else
        {
        mpiMgr.recvRaw(&v1_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
        mpiMgr.recvRaw(&v2_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
        mpiMgr.recvRaw(&v3_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
        for (int m=0; m<NPartThisTask[ThisTask]; ++m)
          {
          vel[m].x=v1_tmp[m];
          vel[m].y=v2_tmp[m];
          vel[m].z=v3_tmp[m];
          }
        }
      }

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
              id[ncount]=ftmp[m];
              ncount++;
              if(ncount == NPartThisTask[ToTask])
                {
                ToTask++;
                ncount=0;
                }
              }
            else
              {
              i1_tmp[ncount]=ftmp[m];
              ncount++;
              if(ncount == NPartThisTask[ToTask])
                {
                mpiMgr.sendRaw(&i1_tmp[0], NPartThisTask[ToTask], ToTask);
                ToTask++;
                ncount=0;
                }
              }
            }
          LastType=type;
          }
          infile.close();
        }
      planck_assert(ncount == 0,"Some particles were left when reading IDs ...");
      }
    else
      {
      mpiMgr.recvRaw(&i1_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
      for (int m=0; m<NPartThisTask[ThisTask]; ++m)
        id[m]=i1_tmp[m];
      }
    }

  if(mpiMgr.master())
    cout << " Reading smoothing ..." << endl;
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
            p[ncount++].r = (fix_size==0.0) ? ftmp[m]*size_fac : fix_size;
            if(ncount == NPartThisTask[ToTask])
              {
              ToTask++;
              ncount=0;
              }
            }
          else
            {
            v1_tmp[ncount++] = (fix_size==0.0) ? ftmp[m]*size_fac : fix_size;
            if(ncount == NPartThisTask[ToTask])
              {
              mpiMgr.sendRaw(&v1_tmp[0], NPartThisTask[ToTask], ToTask);
              ToTask++;
              ncount=0;
              }
            }
          }
          LastType=type;
        }
        infile.close();
      }
    planck_assert(ncount == 0,"Some particles were left when reading sizes ...");
    }
  else
    {
    mpiMgr.recvRaw(&v1_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      p[m].r=v1_tmp[m];
    }


  if(mpiMgr.master())
    cout << " Reading colors ..." << endl;
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
              p[ncount].C1 = ftmp[ofs]*col_fac;
              if(col_vector)
                {
                p[ncount].C2 = ftmp[ofs+1]*col_fac;
                p[ncount].C3 = ftmp[ofs+2]*col_fac;
                }
              }
            else
              p[ncount].C1 = p[ncount].C2 = p[ncount].C3 = 1;

            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              ToTask++;
              ncount=0;
              }
            }
          else
            {
            if (read_col > 0)
              {
              tsize ofs = col_vector ? 3*m : m;
              v1_tmp[ncount] = ftmp[ofs]*col_fac;
              if(col_vector)
                {
                v2_tmp[ncount] = ftmp[ofs+1]*col_fac;
                v3_tmp[ncount] = ftmp[ofs+2]*col_fac;
                }
              }
            else
              v1_tmp[ncount] = v2_tmp[ncount] = v3_tmp[ncount] = 1;

            ncount++;
            if(ncount == NPartThisTask[ToTask])
              {
              mpiMgr.sendRaw(&v1_tmp[0], NPartThisTask[ToTask], ToTask);
              mpiMgr.sendRaw(&v2_tmp[0], NPartThisTask[ToTask], ToTask);
              mpiMgr.sendRaw(&v3_tmp[0], NPartThisTask[ToTask], ToTask);
              ToTask++;
              ncount=0;
              }
            }
          }
        LastType=type;
        }
      infile.close();
      }
    planck_assert(ncount == 0,"Some particles were left when reading colors ...");
    }
  else
    {
    mpiMgr.recvRaw(&v1_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    mpiMgr.recvRaw(&v2_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    mpiMgr.recvRaw(&v3_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      {
      p[m].C1=v1_tmp[m]; p[m].C2=v2_tmp[m]; p[m].C3=v3_tmp[m];
      }
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
            p[ncount++].I = (read_int>0) ? ftmp[m] : 1;
            if(ncount == NPartThisTask[ToTask])
              {
              ToTask++;
              ncount=0;
              }
            }
          else
            {
            v1_tmp[ncount++] = (read_int>0) ? ftmp[m] : 1;
            if(ncount == NPartThisTask[ToTask])
              {
              mpiMgr.sendRaw(&v1_tmp[0], NPartThisTask[ToTask], ToTask);
              ToTask++;
              ncount=0;
              }
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
    mpiMgr.recvRaw(&v1_tmp[0], NPartThisTask[ThisTask], DataFromTask[ThisTask]);
    for (int m=0; m<NPartThisTask[ThisTask]; ++m)
      p[m].I=v1_tmp[m];
    }
  }
