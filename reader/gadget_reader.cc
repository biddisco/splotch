#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>

#include "arr.h"
#include "cxxutils.h"
#include "mpi_support.h"
#include "paramfile.h"
#include "kernel/bstream.h"
#include "kernel/colour.h"
#include "config/config.h"
#include "utils/colourmap.h"

using namespace std;
using namespace RAYPP;

#include "splotch/splotchutils.h"

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
    blocklabel[4] = 0;
    i=3;
    while (blocklabel[i]==32)
      {
      blocklabel[i]=0;
      i--;
      }
    // cout << "blocklabel: <" << blocklabel << "> --> <" << label << ">" << endl;
    file >> blocksize;
    // cout << "blocksize: " << blocksize << endl;
    file.skip(4);
    if (label!=blocklabel)
      {
      file.skip(blocksize);
      blocksize=0;
      }
    }
  return(blocksize-8);
  }

int gadget_read_header(bifstream &file, int *npart,double *massarr,
  double &time,double &redshift, int *npartall)
  {
  int blocksize = gadget_find_block (file,"HEAD");
  planck_assert (blocksize>0, "Header block not found");
  file.skip(4);
  file.get(npart,6);
  file.get(massarr,6);
  file >> time >> redshift;
  file.skip(8);
  file.get(npartall,6);

  return blocksize;
  }

void gadget_reader(paramfile &params, vector<particle_sim> &p)
  {
  int numfiles = params.find<int>("numfiles",1);
  bool doswap = params.find<bool>("swap_endian",true);
  string infilename = params.find<string>("infile");

  string label_int = params.find<string>("intensity0","XXXX");
  string label_col = params.find<string>("color0","XXXX");
  string label_size = params.find<string>("size0","XXXX");
  bool col_vector = params.find<bool>("color_is_vector0",false);
  float fix_size = params.find<float>("fix_size0",1.0);
  float hsml_fac = params.find<float>("hsml_fac0",1.0);
  float col_fac = params.find<float>("col_fac0",1.0);

  int npartthis[6],npartall[6];
  double massarr[6],time,redshift;

  planck_assert(mpiMgr.num_ranks()==1,
    "sorry, reading Gadget files is not yet MPI parellelized ...");

  string filename;
  if(numfiles>1) filename=infilename+"."+dataToString(0);
  else           filename=infilename;
  bifstream infile(filename.c_str(),doswap);
  planck_assert (infile,"could not open input file! <" + filename + ">");
  gadget_read_header(infile,npartthis,massarr,time,redshift,npartall);
  infile.close();

  long npart=npartall[0],ncount=0;
  p.resize(npart);

  for(int f=0;f<numfiles;f++)
    {
    if(numfiles>1) filename=infilename+"."+dataToString(f);
    else           filename=infilename;
    bifstream infile(filename.c_str(),doswap);
    planck_assert (infile,"could not open input file! <" + filename + ">");

    gadget_read_header(infile,npartthis,massarr,time,redshift,npartall);

    cout << "Reading positions ..." << endl;
    gadget_find_block(infile,"POS");
    infile.skip(4);
    for (int m=0; m<npartthis[0]; ++m)
      {
      infile >> p[ncount+m].x >> p[ncount+m].y >> p[ncount+m].z;
      p[ncount+m].type=0;
      }

    if (fix_size == 0.0)
      {
      cout << "Reading smoothing (" << label_size << ")..." << endl;
      gadget_find_block(infile,label_size);
      infile.skip(4);
      for (int m=0; m<npartthis[0]; ++m)
        {
        infile >> p[ncount+m].r;
        p[ncount+m].r *= hsml_fac;
        }
      }
    else
      {
      cout << "Setting smoothing to" << fix_size << " ..." << endl;
      for (int m=0; m<npartthis[0]; ++m)
        p[ncount+m].r = fix_size;
      }

    if (gadget_find_block(infile,label_col) > 0)
      {
      cout << "Reading colors (" << label_col << ")..." << endl;
      infile.skip(4);
      for (int m=0; m<npartthis[0]; ++m)
        {
        infile >> p[ncount+m].C1;
	p[ncount+m].C1 *= col_fac;
        if (col_vector)
	  {
	    infile >> p[ncount+m].C2 >> p[ncount+m].C3;
	    p[ncount+m].C2 *= col_fac;
	    p[ncount+m].C3 *= col_fac;
	  }
        }
      }
    else
      {
      cout << "Setting colors to 0 ..." << endl;
      for (int m=0; m<npartthis[0]; ++m)
        p[ncount+m].C1=1.0;
      }
    if (gadget_find_block(infile,label_int) > 0)
      {
      cout << "Reading intensity (gas) ..." << endl;
      infile.skip(4);
      for (int m=0; m<npartthis[0]; ++m)
      infile >> p[ncount+m].I;
      }
    else
      {
      cout << "Setting intensity to 1 ..." << endl;
      for (int m=0; m<npartthis[0]; ++m)
        p[ncount+m].I=1.0;
      }

    ncount+=npartthis[0];
    }

  planck_assert(npart==ncount,
    "Something went wrong when reading the file !");
  }
