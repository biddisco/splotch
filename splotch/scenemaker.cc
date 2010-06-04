/*
 * Copyright (c) 2004-2010
 *              Martin Reinecke (1), Klaus Dolag (1)
 *               (1) Max-Planck-Institute for Astrophysics
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include "splotch/scenemaker.h"
#include "cxxsupport/walltimer.h"
#include "reader/reader.h"

using namespace std;

// Higher order interpolation would be:
// Time between snapshots (cosmology!)
//    dt=(z2t(h1.redshift)-z2t(h0.redshift))*0.7
// Velocity factors:
//    v_unit1=v_unit/l_unit/sqrt(h1.time)*dt
//    v_unit0=v_unit/l_unit/sqrt(h0.time)*dt
// Delta_v (cosmology)
//    vda=2*(x1-x0)-(v0*v_unit0+v1*v_unit1)
// Delta_t (0..1)
//     t=FLOAT(k)/FLOAT(nint) == frac (!)
// Interpolated positions:
//    x=x0+v0*v_unit0*t+0.5*(v1*v_unit1-v0*v_unit0+vda)*t^2
// Interpolated velocities:
//    v=v0+t*(v1-v0)

void sceneMaker::particle_interpolate(vector<particle_sim> &p, double frac)
  {
  cout << " Time1/2 = " << time1 << "," << time2 << endl;

  double v_unit1, v_unit2;
  if (interpol_mode>1)
    {
    double h = params.find<double>("hubble",0.7);
    double O = params.find<double>("omega",0.3);
    double L = params.find<double>("lambda",0.7);
    double mparsck = 3.0856780e+24;
    double l_unit = params.find<double>("l_unit",3.0856780e+21);
    double v_unit = params.find<double>("v_unit",100000.00);
    double t1 = log(sqrt(L/O*time1*time1*time1)+sqrt((L/O*time1*time1*time1)+1))/1.5/sqrt(L)/h/1e7*mparsck;
    double t2 = log(sqrt(L/O*time2*time2*time2)+sqrt((L/O*time2*time2*time2)+1))/1.5/sqrt(L)/h/1e7*mparsck;
    double dt = (t2 - t1) * h;
    v_unit1=v_unit/l_unit/sqrt(time1)*dt;
    v_unit2=v_unit/l_unit/sqrt(time2)*dt;
    }

  vector<pair<uint32,uint32> > v;
  v.reserve(min(p1.size(),p2.size()));
  {
  tsize i1=0,i2=0;
  while(i1<p1.size() && i2<p2.size())
    {
    if (id1[idx1[i1]]==id2[idx2[i2]])
      {
      v.push_back(pair<uint32,uint32>(idx1[i1],idx2[i2]));
      i1++; i2++;
      }
    else if (id1[idx1[i1]]<id2[idx2[i2]])
      i1++;
    else if (id1[idx1[i1]]>id2[idx2[i2]])
      i2++;
    }
  }

  tsize npart=v.size();
  p.resize(npart);

#pragma omp parallel
{
  tsize i;
#pragma omp for schedule(guided,1000)
  for (i=0; i<npart; ++i)
    {
    tsize i1=v[i].first, i2=v[i].second;
    planck_assert (p1[i1].type==p2[i2].type,
      "interpolate: can not interpolate between different types !");

    vec3f pos;
    if (interpol_mode>1)
      {
      double vda_x = 2 * (p2[i2].x-p1[i1].x) - (vel1[i1].x*v_unit1 + vel2[i2].x*v_unit2);
      double vda_y = 2 * (p2[i2].y-p1[i1].y) - (vel1[i1].y*v_unit1 + vel2[i2].y*v_unit2);
      double vda_z = 2 * (p2[i2].z-p1[i1].z) - (vel1[i1].z*v_unit1 + vel2[i2].z*v_unit2);
      pos.x = p1[i1].x + vel1[i1].x * v_unit1 * frac
           + 0.5 * (vel2[i2].x * v_unit2 - vel1[i1].x * v_unit1 + vda_x) * frac * frac;
      pos.y = p1[i1].y + vel1[i1].y * v_unit1 * frac
           + 0.5 * (vel2[i2].y * v_unit2 - vel1[i1].y * v_unit1 + vda_y) * frac * frac;
      pos.z = p1[i1].z + vel1[i1].z * v_unit1 * frac
           + 0.5 * (vel2[i2].z * v_unit2 - vel1[i1].z * v_unit1 + vda_z) * frac * frac;
      }
    else
      {
      pos.x = (1-frac) * p1[i1].x  + frac*p2[i2].x;
      pos.y = (1-frac) * p1[i1].y  + frac*p2[i2].y;
      pos.z = (1-frac) * p1[i1].z  + frac*p2[i2].z;
      }

    p[i]=particle_sim(
         COLOUR((1-frac) * p1[i1].e.r + frac*p2[i2].e.r,
                (1-frac) * p1[i1].e.g + frac*p2[i2].e.g,
                (1-frac) * p1[i1].e.b + frac*p2[i2].e.b),
         pos.x,pos.y,pos.z,
         (1-frac) * p1[i1].r  + frac*p2[i2].r,
         (1-frac) * p1[i1].I  + frac*p2[i2].I,
         p1[i1].type,p1[i1].active);
    }
}

  if(mpiMgr.master())
    cout << " found " << p.size() << " common particles ..." << endl;
  }

sceneMaker::sceneMaker (paramfile &par)
  : params(par)
  {
  string geometry_file = params.find<string>("geometry_file","");
  geomfile = (geometry_file!="");
  done=false;

  interpol_mode = params.find<int>("interpolation_mode",0);
  if (interpol_mode>0)
    {
    planck_assert(mpiMgr.num_ranks()==1,
       "Sorry, interpolating between files is not yet MPI parallelized ...");

    int snr_start = params.find<int>("snap_start",10);
    snr1=snr_start;
    snr2=snr_start+1;
    snr1_now=snr2_now=-1;
    }

  if (geomfile)
    {
    inp.open(geometry_file.c_str());
    planck_assert (inp, "could not open geometry file '" + geometry_file +"'");
    linecount=ninterpol=nextfile=0;
    geometry_skip = params.find<int>("geometry_start",0);
    geometry_incr = params.find<int>("geometry_incr",1);

    string line;
    for(int i=0; i<geometry_skip; i++, linecount++)
      {
      getline(inp, line);
      if (interpol_mode>0)
        {
        if (linecount==nextfile)
          {
          nextfile=linecount+ninterpol;
          snr1=snr2;
          snr2++;
          }
        }
      }
    }
  }

void sceneMaker::fetchFiles(vector<particle_sim> &particle_data, double &frac)
  {
  if (mpiMgr.master())
    cout << endl << "reading data ..." << endl;
  int simtype = params.find<int>("simtype");
  frac=0;
  if (interpol_mode>0) frac = (linecount-(nextfile-ninterpol))/double(ninterpol);
  switch (simtype)
    {
    case 0:
      bin_reader_tab(params,particle_data);
      break;
    case 1:
      bin_reader_block(params,particle_data);
      break;
    case 2:
      if (interpol_mode>0) // Here only the two datasets are prepared, interpolation will be done later
        {
        cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fac: " << frac << endl;
        cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl;
        cout << " (pos: " << linecount << " , " << nextfile << " , " << ninterpol << ")" << endl;
        if (snr1==snr2_now)
          {
          cout << " old2 = new1!" << endl;
          p1.swap(p2);
          id1.swap(id2);
          idx1.swap(idx2);
          vel1.swap(vel2);

          snr1_now = snr1;
          time1 = time2;
          }
        if (snr1_now!=snr1)
          {
          cout << " reading new1 " << snr1 << endl;
          gadget_reader(params,interpol_mode,p1,id1,vel1,snr1,time1);
          buildIndex(id1.begin(),id1.end(),idx1);
          snr1_now = snr1;
          }
        if (snr2_now!=snr2)
          {
          cout << " reading new2 " << snr2 << endl;
          gadget_reader(params,interpol_mode,p2,id2,vel2,snr2,time2);
          buildIndex(id2.begin(),id2.end(),idx2);
          snr2_now = snr2;
          }
        }
      else
        {
        double dummy;
        gadget_reader(params,interpol_mode,particle_data,id1,vel1,0,dummy);
        if (geomfile)
          p_orig = particle_data;
        }
      break;
#if 0
    case 3:
      enzo_reader(params,particle_data);
      break;
#endif
    case 4:
      {
      double dummy;
      gadget_millenium_reader(params,particle_data,0,&dummy);
      break;
      }
    case 5:
#if defined(USE_MPIIO)
      {
      float maxr, minr;
      bin_reader_block_mpi(params,particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
      }
#else
      planck_fail("mpi reader not available in non MPI compiled version !");
#endif
      break;
    case 6:
      mesh_reader(params,particle_data);
      break;
#ifdef HDF5
    case 7:
      hdf5_reader(params,particle_data);
      break;
#endif
    default:
      planck_fail("No valid file type given ...");
      break;
    }
  }

bool sceneMaker::getNextScene (vector<particle_sim> &particle_data,
  vec3 &campos, vec3 &lookat, vec3 &sky, string &outfile)
  {
  wallTimers.start("read");
  if ((!geomfile)&&done) return false;

  bool master = mpiMgr.master();

  if (geomfile)
    {
    string line;
    if (!getline(inp, line)) return false;
    sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %i",
      &campos.x,&campos.y,&campos.z,
      &lookat.x,&lookat.y,&lookat.z,
      &sky.x,&sky.y,&sky.z,&ninterpol);
    if (master)
      {
      cout << endl << "Next entry <" << linecount << "> in geometry file ..." << endl;
      cout << " Camera:    " << campos << endl;
      cout << " Lookat:    " << lookat << endl;
      cout << " Sky:       " << sky << endl;
      if (interpol_mode>0)
        cout << " ninterpol: " << ninterpol << endl;
      }
    if ((interpol_mode>0) && (linecount==0) && (nextfile==0))
      nextfile=linecount+ninterpol;
    }
  else
    {
    campos.x=params.find<double>("camera_x");
    campos.y=params.find<double>("camera_y");
    campos.z=params.find<double>("camera_z");
    lookat.x=params.find<double>("lookat_x");
    lookat.y=params.find<double>("lookat_y");
    lookat.z=params.find<double>("lookat_z");
    sky.x=params.find<double>("sky_x",0);
    sky.y=params.find<double>("sky_y",0);
    sky.z=params.find<double>("sky_z",0);
    }

// -----------------------------------
// ----------- Reading ---------------
// -----------------------------------

  double frac;
  if (interpol_mode==0)
    {
    if ((!geomfile) || (geomfile&&(linecount==geometry_skip))) // read only once if no interpolation is chosen
      fetchFiles(particle_data,frac);
    else
      if (geomfile) particle_data = p_orig;
    }
  else
    fetchFiles(particle_data,frac);

  if (interpol_mode>0)
    {
    if (master)
      cout << "Interpolating between " << p1.size() << " and " <<
        p2.size() << " particles ..." << endl;
    particle_interpolate(particle_data,frac);
    }
  outfile = params.find<string>("outfile");
  if (geomfile)
    {
    outfile += intToString(linecount,4) + ".tga";
    linecount++;
    if ((interpol_mode>0) && (linecount==nextfile))
      {
      nextfile=linecount+ninterpol;
      snr1=snr2;
      snr2++;
      }

    string line;
    for (int i=1; i<geometry_incr; i++)
      {
      getline(inp, line);
      linecount++;
      if ((interpol_mode>0) && (linecount==nextfile))
        {
        nextfile=linecount+ninterpol;
        snr1=snr2;
        snr2++;
        }
      }
    }

  done=true;

  wallTimers.stop("read");
  return true;
  }
