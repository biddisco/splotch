#include "splotch/scenemaker.h"
#include "cxxsupport/walltimer.h"
#include "reader/reader.h"

using namespace std;

#ifdef INTERPOLATE

namespace {

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

void particle_interpolate(paramfile &params, vector<particle_sim> &p,
  const vector<particle_sim> &p1, const vector<particle_sim> &p2, double frac,
  double time1, double time2)
  {
  cout << " Time1/2 = " << time1 << "," << time2 << endl;

#ifdef HIGH_ORDER_INTERPOLATION
  double h = params.find<double>("hubble",0.7);
  double O = params.find<double>("omega",0.3);
  double L = params.find<double>("lambda",0.7);
  double mparsck = 3.0856780e+24;
  double l_unit = params.find<double>("l_unit",3.0856780e+21);
  double v_unit = params.find<double>("v_unit",100000.00);
  double t1 = log(sqrt(L/O*time1*time1*time1)+sqrt((L/O*time1*time1*time1)+1))/1.5/sqrt(L)/h/1e7*mparsck;
  double t2 = log(sqrt(L/O*time2*time2*time2)+sqrt((L/O*time2*time2*time2)+1))/1.5/sqrt(L)/h/1e7*mparsck;
  double dt = (t2 - t1) * h;
  double v_unit1=v_unit/l_unit/sqrt(time1)*dt;
  double v_unit2=v_unit/l_unit/sqrt(time2)*dt;
#endif

  vector<pair<uint32,uint32> >v;
  v.reserve(min(p1.size(),p2.size()));
  tsize i1=0,i2=0;
  while(i1<p1.size() && i2<p2.size())
    {
    if (p1[i1].id==p2[i2].id)
      {
      v.push_back(pair<uint32,uint32>(i1,i2));
      i1++; i2++;
      }
    else if (p1[i1].id<p2[i2].id)
      i1++;
    else if (p1[i1].id>p2[i2].id)
      i2++;
    }

  tsize npart=v.size();
  p.resize(npart);

#pragma omp parallel
{
  int i;
#pragma omp for schedule(guided,10000)
  for (i=0; i<npart; ++i)
    {
    tsize i1=v[i].first, i2=v[i].second;
    planck_assert (p1[i1].type==p2[i2].type,
      "interpolate: can not interpolate between different types !");
#ifdef HIGH_ORDER_INTERPOLATION
    double vda_x = 2 * (p2[i2].x-p1[i1].x) - (p1[i1].vx*v_unit1 + p2[i2].vx*v_unit2);
    double vda_y = 2 * (p2[i2].y-p1[i1].y) - (p1[i1].vy*v_unit1 + p2[i2].vy*v_unit2);
    double vda_z = 2 * (p2[i2].z-p1[i1].z) - (p1[i1].vz*v_unit1 + p2[i2].vz*v_unit2);
#endif
    p[i]=particle_sim(
#ifdef HIGH_ORDER_INTERPOLATION
         p1[i1].x + p1[i1].vx * v_unit1 * frac
           + 0.5 * (p2[i2].vx * v_unit2 - p1[i1].vx * v_unit1 + vda_x) * frac * frac,
         p1[i1].y + p1[i1].vy * v_unit1 * frac
           + 0.5 * (p2[i2].vy * v_unit2 - p1[i1].vy * v_unit1 + vda_y) * frac * frac,
         p1[i1].z + p1[i1].vz * v_unit1 * frac
           + 0.5 * (p2[i2].vz * v_unit2 - p1[i1].vz * v_unit1 + vda_z) * frac * frac,
#else
         (1-frac) * p1[i1].x  + frac*p2[i2].x,
         (1-frac) * p1[i1].y  + frac*p2[i2].y,
         (1-frac) * p1[i1].z  + frac*p2[i2].z,
#endif
         (1-frac) * p1[i1].r  + frac*p2[i2].r,
         (1-frac) * p1[i1].ro + frac*p2[i2].ro,
         (1-frac) * p1[i1].I  + frac*p2[i2].I,
         (1-frac) * p1[i1].C1 + frac*p2[i2].C1,
         (1-frac) * p1[i1].C2 + frac*p2[i2].C2,
         (1-frac) * p1[i1].C3 + frac*p2[i2].C3,
         p1[i1].type,p1[i1].active,p1[i1].e,p1[i1].id
#ifdef HIGH_ORDER_INTERPOLATION
        ,(1-frac) * p1[i1].vx  + frac*p2[i2].vx,
         (1-frac) * p1[i1].vy  + frac*p2[i2].vy,
         (1-frac) * p1[i1].vz  + frac*p2[i2].vz
#endif
         );
    }
}

  if(mpiMgr.master())
    cout << " found " << p.size() << " common particles ..." << endl;
  }

} // unnamed namespace

#endif

sceneMaker::sceneMaker (paramfile &par)
  : params(par)
  {
  string geometry_file = params.find<string>("geometry_file","");
  geomfile = (geometry_file!="");
  done=false;
#ifdef INTERPOLATE
  snr_start = params.find<int>("snap_start",10);
  snr1=snr_start;
  snr2=snr_start+1;
  snr1_now=snr2_now=-1;
#endif

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
#ifdef INTERPOLATE
      if (linecount==nextfile)
        {
        nextfile=linecount+ninterpol;
        snr1=snr2;
        snr2++;
        }
#endif
      }
    }
  }

bool sceneMaker::getNextScene (vector<particle_sim> &particle_data,
  vec3 &campos, vec3 &lookat, vec3 &sky, string &outfile)
  {
  wallTimers.start("read");
  if ((!geomfile)&&done) return false;

  bool master = mpiMgr.master();
  double dummy;

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
#ifdef INTERPOLATE
      cout << " ninterpol: " << ninterpol << endl;
#endif
      }
#ifdef INTERPOLATE
    if(linecount == 0 && nextfile == 0)
      nextfile=linecount+ninterpol;
#endif
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

#ifndef INTERPOLATE
  if ((!geomfile) || (geomfile&&(linecount==geometry_skip))) // read only once if no interpolation is chosen
    {
#endif
    if (master)
      cout << endl << "reading data ..." << endl;
    int simtype = params.find<int>("simtype"); // 2:Gadget2
    float maxr, minr;
#ifdef INTERPOLATE
    double frac=(linecount-(nextfile-ninterpol))/double(ninterpol);
#endif
    switch (simtype)
      {
      case 0:
        bin_reader_tab(params,particle_data);
        break;
      case 1:
        bin_reader_block(params,particle_data);
        break;
      case 2:
#ifdef INTERPOLATE // Here only the two datasets are prepared, interpolation will be done later
        cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fac: " << frac << endl;
        cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl;
        cout << " (pos: " << linecount << " , " << nextfile << " , " << ninterpol << ")" << endl;
        if (snr1==snr2_now)
          {
          cout << " old2 = new1!" << endl;
          particle_data1=particle_data2;
          snr1_now = snr1;
          time1 = time2;
          }
        if (snr1_now!=snr1)
          {
          cout << " reading new1 " << snr1 << endl;
          gadget_reader(params,particle_data1,snr1,&time1);
          snr1_now = snr1;
          }
        if (snr2_now!=snr2)
          {
          cout << " reading new2 " << snr2 << endl;
          gadget_reader(params,particle_data2,snr2,&time2);
          snr2_now = snr2;
          }
#else
        gadget_reader(params,particle_data,0,&dummy);
        if (geomfile)
          p_orig = particle_data;
#endif
        break;
#if 0
      case 3:
        enzo_reader(params,particle_data);
        break;
#endif
      case 4:
        gadget_millenium_reader(params,particle_data,0,&dummy);
        break;
      case 5:
#if defined(USE_MPIIO)
        bin_reader_block_mpi(params,particle_data, &maxr, &minr, mpiMgr.rank(), mpiMgr.num_ranks());
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
#ifndef INTERPOLATE
    }
  else
    if (geomfile) particle_data = p_orig;
#endif

#ifdef INTERPOLATE
  if (master)
    cout << "Interpolating between " << particle_data1.size() << " and " <<
      particle_data2.size() << " particles ..." << endl;
  particle_interpolate(params,particle_data,particle_data1,particle_data2,frac,time1,time2);
#endif
  outfile = params.find<string>("outfile");
  if (geomfile)
    {
    outfile += intToString(linecount,4) + ".tga";
    linecount++;
#ifdef INTERPOLATE
    if (linecount==nextfile)
      {
      nextfile=linecount+ninterpol;
      snr1=snr2;
      snr2++;
      }
#endif

    string line;
    for (int i=1; i<geometry_incr; i++)
      {
      getline(inp, line);
      linecount++;
#ifdef INTERPOLATE
      if (linecount==nextfile)
        {
        nextfile=linecount+ninterpol;
        snr1=snr2;
        snr2++;
        }
#endif
      }
    }

  done=true;

  wallTimers.stop("read");
  return true;
  }
