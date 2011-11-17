/*
 * Copyright (c) 2004-2011
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
#include <iostream>
#include <fstream>

#include "splotch/scenemaker.h"
#include "cxxsupport/lsconstants.h"
#include "cxxsupport/walltimer.h"
#include "cxxsupport/cxxutils.h"
#include "reader/reader.h"
//boost
#include "booster/mesh_vis.h"

using namespace std;

void sceneMaker::particle_normalize(std::vector<particle_sim> &p, bool verbose)
{
  int nt = params.find<int>("ptypes",1);
  arr<bool> col_vector(nt),log_int(nt),log_col(nt),asinh_col(nt);
  arr<Normalizer<float32> > intnorm(nt), colnorm(nt);

  for(int t=0; t<nt; t++)
  {
    log_int[t] = params.find<bool>("intensity_log"+dataToString(t),false);
    log_col[t] = params.find<bool>("color_log"+dataToString(t),false);
    asinh_col[t] = params.find<bool>("color_asinh"+dataToString(t),false);
    col_vector[t] = params.find<bool>("color_is_vector"+dataToString(t),false);
  }

  int npart=p.size();
  tstack_push("minmax");
  #pragma omp parallel
  {
    // FIXME: the "+20" serves as protection against false sharing,
    // should be done more elegantly
    arr<Normalizer<float32> > inorm(nt+20), cnorm(nt+20);
    int m;
    #pragma omp for schedule(guided,1000)
    for (m=0; m<npart; ++m) // do log calculations if requested
    {
      int t=p[m].type;

      if (log_int[t])
	{
	  if(p[m].I > 0)
	    {
	      p[m].I = log10(p[m].I);
	      inorm[t].collect(p[m].I);
	    }
	  else
	    {
	      p[m].I = -38;
	    }
	}
      else
	inorm[t].collect(p[m].I);

      if (log_col[t])
	{
	  if(p[m].e.r > 0)
	    {
	      p[m].e.r = log10(p[m].e.r);
	      cnorm[t].collect(p[m].e.r);
	    }
	  else
	    {
	      p[m].e.r =-38;
	    }
	}
      else 
	{
	  if (asinh_col[t])
	    p[m].e.r = my_asinh(p[m].e.r);
	  cnorm[t].collect(p[m].e.r);
	}

      if (col_vector[t])
      {
        if (log_col[t])
        {
          p[m].e.g = log10(p[m].e.g);
          p[m].e.b = log10(p[m].e.b);
        }
        if (asinh_col[t])
        {
          p[m].e.g = my_asinh(p[m].e.g);
          p[m].e.b = my_asinh(p[m].e.b);
        }
        cnorm[t].collect(p[m].e.g);
        cnorm[t].collect(p[m].e.b);
      }
    }
    #pragma omp critical
    for(int t=0; t<nt; t++)
    {
      intnorm[t].collect(inorm[t]);
      colnorm[t].collect(cnorm[t]);
    }
  }

  for(int t=0; t<nt; t++)
  {
    mpiMgr.allreduce(intnorm[t].minv,MPI_Manager::Min);
    mpiMgr.allreduce(colnorm[t].minv,MPI_Manager::Min);
    mpiMgr.allreduce(intnorm[t].maxv,MPI_Manager::Max);
    mpiMgr.allreduce(colnorm[t].maxv,MPI_Manager::Max);

    if (verbose && mpiMgr.master())
    {
      cout << " For particles of type " << t << ":" << endl;
      cout << " From data: " << endl;
      cout << " Color Range:     " << colnorm[t].minv << " (min) , " <<
           colnorm[t].maxv << " (max) " << endl;
      cout << " Intensity Range: " << intnorm[t].minv << " (min) , " <<
           intnorm[t].maxv << " (max) " << endl;
    }

    if(params.param_present("intensity_min"+dataToString(t)))
      intnorm[t].minv = params.find<float>("intensity_min"+dataToString(t));

    if(params.param_present("intensity_max"+dataToString(t)))
      intnorm[t].maxv = params.find<float>("intensity_max"+dataToString(t));

    if(params.param_present("color_min"+dataToString(t)))
      colnorm[t].minv = params.find<float>("color_min"+dataToString(t));

    if(params.param_present("color_max"+dataToString(t)))
      colnorm[t].maxv = params.find<float>("color_max"+dataToString(t));

    if (verbose && mpiMgr.master())
    {
      cout << " Restricted to: " << endl;
      cout << " Color Range:     " << colnorm[t].minv << " (min) , " <<
           colnorm[t].maxv << " (max) " << endl;
      cout << " Intensity Range: " << intnorm[t].minv << " (min) , " <<
           intnorm[t].maxv << " (max) " << endl;
    }
  }
  tstack_pop("minmax");

  // set variables for the output to the image log file
  if(params.param_present("color_min0"))
    colmin0=params.find<float>("color_min0");
  else
    colmin0=colnorm[0].minv;

  if(params.param_present("color_max0"))
    colmax0=params.find<float>("color_max0");
  else
    colmax0=colnorm[0].maxv;

  if(params.param_present("brightness0"))
    bright0=params.find<float>("brightness0");
  else
    bright0=1.0;

  tstack_push("clamp");

  #pragma omp parallel
  {
    int m;
    #pragma omp for schedule(guided,1000)
    for(m=0; m<npart; ++m)
    {
      int t=p[m].type;
      intnorm[t].normAndClamp(p[m].I);
      colnorm[t].normAndClamp(p[m].e.r);
      if (col_vector[t])
      {
        colnorm[t].normAndClamp(p[m].e.g);
        colnorm[t].normAndClamp(p[m].e.b);
      }
    }
  }
  tstack_pop("clamp");
} // END particle_normalize

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

// booster main variables

Mesh_vis * Mesh = NULL;
Mesh_dim MeshD;

void sceneMaker::particle_interpolate(vector<particle_sim> &p, double frac)
{
  cout << " Time1/2 = " << time1 << "," << time2 << endl;

  releaseMemory(p);

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

  vector<pair<MyIDType,MyIDType> > v;
  v.reserve(min(p1.size(),p2.size()));
  {
    tsize i1=0,i2=0;
    while(i1<p1.size() && i2<p2.size())
    {
      if (id1[idx1[i1]]==id2[idx2[i2]])
      {
        //  if(p1[idx1[i1]].type==p2[idx2[i2]].type)
        v.push_back(pair<MyIDType,MyIDType>(idx1[i1],idx2[i2]));
        i1++;
        i2++;
      }
      else if (id1[idx1[i1]]<id2[idx2[i2]])
        i1++;
      else if (id1[idx1[i1]]>id2[idx2[i2]])
        i2++;
    }
  }

  tsize npart=v.size();
  p.resize(npart);

  bool periodic = params.find<bool>("periodic",true);
  double boxhalf = boxsize / 2;

  #pragma omp parallel
  {
    int i;
    #pragma omp for schedule(guided,1000)
    for (i=0; i<npart; ++i)
    {
      tsize i1=v[i].first, i2=v[i].second;
      /*
      planck_assert (p1[i1].type==p2[i2].type,
        "interpolate: cannot interpolate between different particle types!");
      */
      vec3f pos;
      double x1,x2,y1,y2,z1,z2;

      x1 = p1[i1].x;
      x2 = p2[i2].x;
      y1 = p1[i1].y;
      y2 = p2[i2].y;
      z1 = p1[i1].z;
      z2 = p2[i2].z;

      if (periodic)
      {
        if(abs(x2 - x1) > boxhalf)
          (x2 > x1) ? x2 -= boxsize : x2 += boxsize;

        if(abs(y2 - y1) > boxhalf)
          (y2 > y1) ? y2 -= boxsize : y2 += boxsize;

        if(abs(z2 - z1) > boxhalf)
          (z2 > z1) ? z2 -= boxsize : z2 += boxsize;
      }
      if (interpol_mode>1)
      {
        double vda_x = 2 * (x2-x1) - (vel1[i1].x*v_unit1 + vel2[i2].x*v_unit2);
        double vda_y = 2 * (y2-y1) - (vel1[i1].y*v_unit1 + vel2[i2].y*v_unit2);
        double vda_z = 2 * (z2-z1) - (vel1[i1].z*v_unit1 + vel2[i2].z*v_unit2);
        pos.x = x1 + vel1[i1].x * v_unit1 * frac
                + 0.5 * (vel2[i2].x * v_unit2 - vel1[i1].x * v_unit1 + vda_x) * frac * frac;
        pos.y = y1 + vel1[i1].y * v_unit1 * frac
                + 0.5 * (vel2[i2].y * v_unit2 - vel1[i1].y * v_unit1 + vda_y) * frac * frac;
        pos.z = z1 + vel1[i1].z * v_unit1 * frac
                + 0.5 * (vel2[i2].z * v_unit2 - vel1[i1].z * v_unit1 + vda_z) * frac * frac;
      }
      else
      {
        pos.x = (1-frac) * x1  + frac*x2;
        pos.y = (1-frac) * y1  + frac*y2;
        pos.z = (1-frac) * z1  + frac*z2;
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

  // Set the time variables for the output to the image log file
  // which is written in getNextScene()
  intTime     = (1.-frac)*time1     + frac*time2;
  intRedshift = (1.-frac)*redshift1 + frac*redshift2;

  if(mpiMgr.master())
    cout << " found " << p.size() << " common particles ..." << endl;

}


sceneMaker::sceneMaker (paramfile &par)
  : cur_scene(-1), params(par), snr1_now(-1), snr2_now(-1)
{
  double fidx = params.find<double>("fidx",0);

  std::map<std::string,std::string> sceneParameters;
  sceneParameters.clear();

  string outfile = params.find<string>("outfile","demo");

  // do nothing if we are only analyzing ...
  if (params.find<bool>("AnalyzeSimulationOnly",false))
  {
    string outfilen = outfile+intToString(0,4);
    scenes.push_back(scene(sceneParameters,outfilen,false,false));
    return;
  }

  string geometry_file = params.find<string>("scene_file","");
  interpol_mode = params.find<int>("interpolation_mode",0);
  if (geometry_file=="")
  {
    string outfilen = outfile+intToString(0,4);
    scenes.push_back(scene(sceneParameters,outfilen,false,false));
  }
  else
  {
    if (interpol_mode>0)
      planck_assert(mpiMgr.num_ranks()==1,
                    "Sorry, interpolating between files is not yet MPI parallelized ...");
    ifstream inp(geometry_file.c_str());
    planck_assert (inp, "could not open scene file '" + geometry_file +"'");
    int current_scene = params.find<int>("scene_start",0);
    int scene_incr = params.find<int>("scene_incr",1);
    //
    string line;
    getline(inp, line);
    double tmpDbl;
    if (sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
          &tmpDbl,&tmpDbl,&tmpDbl,&tmpDbl,&tmpDbl,&tmpDbl,&tmpDbl,&tmpDbl,&tmpDbl,&tmpDbl)==10)
    {
      //      cerr << "DEBUG: old geometry file format detected." << endl;
      line.assign("camera_x camera_y camera_z lookat_x lookat_y lookat_z sky_x sky_y sky_z fidx");
      inp.seekg(0, ios_base::beg);
    }
    else
    {
      //      cerr << "DEBUG: new geometry file format detected." << endl;
    }
    //
    std::vector<std::string> sceneParameterKeys, sceneParameterValues;
    split(line, sceneParameterKeys);

    cout << endl << "Found following parameters to be modulated in the scene file: " << endl;
    for(unsigned int i=0;i<sceneParameterKeys.size();i++)
      cout << sceneParameterKeys[i] << endl;

    for (int i=0; i<current_scene; ++i)
      getline(inp, line);

    while (getline(inp, line))
    {
      sceneParameters.clear();
      string outfilen = outfile+intToString(current_scene,4);
      split(line, sceneParameterValues);

      planck_assert(sceneParameterKeys.size()==sceneParameterValues.size(),
        "ERROR in scene file detected, please check!  Quitting.");
      for (unsigned int i=0; i<sceneParameterKeys.size(); i++)
        sceneParameters[sceneParameterKeys[i]] = sceneParameterValues[i];
      bool reuse=false;

      fidx = params.find<double>("fidx",0);

      if (scenes.size()>0)
        if (approx(fidx,stringToData<double>(scenes[scenes.size()-1].sceneParameters.find("fidx")->second)))
          scenes[scenes.size()-1].keep_particles=reuse=true;

      scenes.push_back(scene(sceneParameters,outfilen,false,reuse));
      current_scene += scene_incr;
      for (int i=0; i<scene_incr-1; ++i)
        getline(inp, line);
    }
  }
  double eye_separation = degr2rad * params.find<double>("EyeSeparation",0);
  if (eye_separation>0)
  {
    vector<scene> sc_orig;
    sc_orig.swap(scenes);
    for (tsize i=0; i<sc_orig.size(); ++i)
    {
      scenes.push_back(sc_orig[i]);
      scenes.push_back(sc_orig[i]);
      scene &sa = scenes[scenes.size()-2], &sb = scenes[scenes.size()-1];
      sa.keep_particles=true;
      sb.keep_particles=false;
      sb.reuse_particles=true;

      vec3 lookat(params.find<double>("lookat_x"),params.find<double>("lookat_y"),params.find<double>("lookat_z"));
      vec3 campos(params.find<double>("campos_x"),params.find<double>("campos_y"),params.find<double>("campos_z"));
      vec3 sky(params.find<double>("sky_x",0),params.find<double>("sky_y",0),params.find<double>("sky_z",1));

      vec3 view = lookat - campos;

      // Real sky vector 'sky_real' is the given sky vector 'sky' projected into the plane
      // which lies orthogonal to the looking vector 'view', which connects the
      // camera 'campos' with the lookat point 'look'

      double cosa = dotprod (view,sky) / (view.Length() * sky.Length());

      vec3 sky_real = sky - view * cosa * sky.Length() / view.Length();
      vec3 right = crossprod (sky,view);

      double distance = eye_separation * view.Length();

      vec3 campos_r = campos - right / right.Length() * distance*0.5;
      sa.sceneParameters["campos_x"] = dataToString(campos_r.x);
      sa.sceneParameters["campos_y"] = dataToString(campos_r.y);
      sa.sceneParameters["campos_z"] = dataToString(campos_r.z);

      vec3 campos_l = campos + right / right.Length() * distance*0.5;
      sb.sceneParameters["campos_x"] = dataToString(campos_l.x);
      sb.sceneParameters["campos_y"] = dataToString(campos_l.y);
      sb.sceneParameters["campos_z"] = dataToString(campos_l.z);

      sa.outname = "left_"+sa.outname;
      sb.outname = "right_"+sb.outname;
    }
  }
}

void sceneMaker::fetchFiles(vector<particle_sim> &particle_data, double fidx)
{
  if (scenes[cur_scene].reuse_particles)
  {
    particle_data=p_orig;
    return;
  }

  tstack_push("Input");
  if (mpiMgr.master())
    cout << endl << "reading data ..." << endl;
  int simtype = params.find<int>("simtype");
  int spacing = params.find<double>("snapshot_spacing",1);
  int snr1 = int(fidx/spacing)*spacing, snr2=snr1+spacing;
  double frac=(fidx-snr1)/spacing;

  switch (simtype)
  {
  case 0:
    bin_reader_tab(params,particle_data);
    break;
  case 1:
    bin_reader_block(params,particle_data);
    break;
  case 2:
    if (interpol_mode>0) // Here only the two data sets are prepared, interpolation will be done later
    {
      cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fraction: " << frac << endl;
      cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl;
      if (snr1==snr2_now)
      {
        cout << " old2 = new1!" << endl;
        p1.swap(p2);
        id1.swap(id2);
        idx1.swap(idx2);
        vel1.swap(vel2);

        snr1_now = snr1;
        time1 = time2;
        redshift1 = redshift2;
      }
      if (snr1_now!=snr1)
      {
        cout << " reading new1 " << snr1 << endl;
        gadget_reader(params,interpol_mode,p1,id1,vel1,snr1,time1,boxsize);
	mpiMgr.barrier();
        tstack_replace("Input","Particle index generation");
        buildIndex(id1.begin(),id1.end(),idx1);
        tstack_replace("Particle index generation","Input");
        snr1_now = snr1;
      }
      if (snr2_now!=snr2)
      {
        cout << " reading new2 " << snr2 << endl;
        gadget_reader(params,interpol_mode,p2,id2,vel2,snr2,time2,boxsize);
	mpiMgr.barrier();
        tstack_replace("Input","Particle index generation");
        buildIndex(id2.begin(),id2.end(),idx2);
        tstack_replace("Particle index generation","Input");
        snr2_now = snr2;
      }
    }
    else
    {
      intRedshift=-1.0;  // (khr) Is redshift available by default from Gadget data?
      gadget_reader(params,interpol_mode,particle_data,id1,vel1,snr1,intTime,boxsize);
    }
    break;
  case 3:
#if 0
    enzo_reader(params,particle_data);
#else
    planck_fail("Enzo reader not available in this version!");
#endif
    break;
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
  planck_fail("mpi reader not available in non MPI compiled version!");
#endif
  break;
  case 6:
    mesh_reader(params,particle_data);
    break;
#ifdef HDF5
  case 7:
    hdf5_reader(params,particle_data);
    break;
  case 8:
    // GADGET HDF5 READER
    //
    if (interpol_mode>0) // Here only the two data sets are prepared, interpolation will be done later
    {
      cout << "Loaded file1: " << snr1_now << " , file2: " << snr2_now << " , interpol fraction: " << frac << endl;
      cout << " (needed files : " << snr1 << " , " << snr2 << ")" << endl;
      if (snr1==snr2_now)
      {
        cout << " old2 = new1!" << endl;
        p1.swap(p2);
        id1.swap(id2);
        idx1.swap(idx2);
        vel1.swap(vel2);

        snr1_now = snr1;
        time1 = time2;
        redshift1 = redshift2;
      }
      if (snr1_now!=snr1)
      {
        cout << " reading new1 " << snr1 << endl;
        gadget_hdf5_reader(params,interpol_mode,p1,id1,vel1,snr1,time1,redshift1,boxsize);
	mpiMgr.barrier();
        tstack_replace("Input","Particle index generation");
        buildIndex(id1.begin(),id1.end(),idx1);
        tstack_replace("Particle index generation","Input");
        snr1_now = snr1;
      }
      if (snr2_now!=snr2)
      {
        cout << " reading new2 " << snr2 << endl;
        gadget_hdf5_reader(params,interpol_mode,p2,id2,vel2,snr2,time2,redshift2,boxsize);
	mpiMgr.barrier();
        tstack_replace("Input","Particle index generation");
        buildIndex(id2.begin(),id2.end(),idx2);
        tstack_replace("Particle index generation","Input");
        snr2_now = snr2;
      }
    }
    else
    {
      gadget_hdf5_reader(params,interpol_mode,particle_data,id1,vel1,0,intTime,intRedshift,boxsize);
    }
    break;
#endif
#ifdef SPLVISIVO
  case 10:
    if(!visivo_reader(params,particle_data,opt))
      planck_fail("Invalid read data ...");
    break;
#endif
    case 11:
      tipsy_reader(params,particle_data);
	break;
    default:
      planck_fail("No valid file type given ...");
      break;
    }

  mpiMgr.barrier();
  tstack_pop("Input");

  if (interpol_mode>0)
  {
    if (mpiMgr.master())
      cout << "Interpolating between " << p1.size() << " and " <<
           p2.size() << " particles ..." << endl;
    tstack_push("Time interpolation");
    particle_interpolate(particle_data,frac);
    tstack_pop("Time interpolation");
  }
  tstack_push("Particle ranging");
  tsize npart_all = particle_data.size();
  mpiMgr.allreduce (npart_all,MPI_Manager::Sum);
  if (mpiMgr.master())
    cout << endl << "host: ranging values (" << npart_all << ") ..." << endl;
  particle_normalize(particle_data, true);
  tstack_pop("Particle ranging");

  if (scenes[cur_scene].keep_particles) p_orig = particle_data;

// boost initialization

  bool boost = params.find<bool>("boost",false);
  if(boost)
  {
    cout << "Boost setup..." << endl;
    mesh_creator(particle_data, &Mesh, &MeshD);
    randomizer(particle_data, Mesh, MeshD);
  }

}

bool sceneMaker::getNextScene (vector<particle_sim> &particle_data, vector<particle_sim> &r_points,
                               vec3 &campos, vec3 &lookat, vec3 &sky, string &outfile)
{
  if (tsize(++cur_scene) >= scenes.size()) return false;

  const scene &scn=scenes[cur_scene];

  // patch the params object with parameter values relevant to the current scene
  std::map<std::string,std::string> sceneParameters=scn.sceneParameters;
  for (std::map<std::string,std::string>::iterator it=sceneParameters.begin(); it!=sceneParameters.end(); ++it)
  {
    params.setParam(it->first, it->second);
  }

  // Fetch the values from the param object which may have been altered by the scene file or copied from the opt object.
  campos=vec3(params.find<double>("camera_x"),params.find<double>("camera_y"),params.find<double>("camera_z"));
  lookat=vec3(params.find<double>("lookat_x"),params.find<double>("lookat_y"),params.find<double>("lookat_z"));
  sky   =vec3(params.find<double>("sky_x",0), params.find<double>("sky_y",0), params.find<double>("sky_z",1));

  outfile=scn.outname;
  double fidx=params.find<double>("fidx",0);

  fetchFiles(particle_data,fidx);

  if (params.find<bool>("periodic",true))
  {
    int npart = particle_data.size();
    double boxhalf = boxsize / 2;

    if(mpiMgr.master())
      cout << " doing parallel box wrap " << boxsize << endl;
    #pragma omp parallel
    {
      int m;
      #pragma omp for schedule(guided,1000)
      for (m=0; m<npart; ++m)
      {
        if(particle_data[m].x - lookat.x > boxhalf)
          particle_data[m].x -= boxsize;
        if(lookat.x - particle_data[m].x > boxhalf)
          particle_data[m].x += boxsize;
        if(particle_data[m].y - lookat.y > boxhalf)
          particle_data[m].y -= boxsize;
        if(lookat.y - particle_data[m].y > boxhalf)
          particle_data[m].y += boxsize;
        if(particle_data[m].z - lookat.z > boxhalf)
          particle_data[m].z -= boxsize;
        if(lookat.z - particle_data[m].z > boxhalf)
          particle_data[m].z += boxsize;
      }
    }
  }

  // Let's try to boost!!!
  bool boost = params.find<bool>("boost",false);
  if(boost)
  {
    cout << "Boost!!!" << endl;
    m_rotation(params, &Mesh, MeshD, campos, lookat, sky);
    p_selector(particle_data, Mesh, MeshD, r_points);
  }

  // Dump information on the next rendered image into a similarly named file
  // which makes image postprocessing with additional scripts *much* easier.
  // For the moment, this is possible with Gadget data only.
  int simtype = params.find<int>("simtype");
  if ( mpiMgr.master() && ((simtype==2)||(simtype==8)))
  {
    string logFileName;
    logFileName.assign(outfile);
    logFileName.append(".log");

    ofstream logFile(logFileName.c_str());
    //
    logFile << "campos "   << campos; // endl is included in the overloaded "<<" operator
    logFile << "lookat "   << lookat;
    logFile << "sky "      << sky;
    //
    if (interpol_mode > 0)
    {
      logFile << "fidx "   << fidx    << endl;
    }
    logFile << "redshift " << intRedshift << endl;
    logFile << "time "     << intTime     << endl;
    //
    logFile << "colmin0 "  << colmin0 << endl;
    logFile << "colmax0 "  << colmax0 << endl;
    logFile << "bright0 "  << bright0 << endl;
    //
    logFile.close();
  }

  return true;
}
