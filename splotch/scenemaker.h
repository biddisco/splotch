#ifndef SPLOTCH_SCENE_MAKER_H
#define SPLOTCH_SCENE_MAKER_H

#include <vector>
#include <string>
#include <fstream>

#include "cxxsupport/paramfile.h"
#include "splotch/splotchutils.h"

class sceneMaker
  {
  private:
    paramfile &params;

    bool geomfile;
    int interpol_mode;

// only used if interpol_mode>0
    std::vector<particle_sim> p1,p2;
    std::vector<uint32> id1,id2,idx1,idx2;
    int snr1,snr2,snr1_now,snr2_now;
    double time1,time2;
// only used if interpol_mode>1
    std::vector<vec3f> vel1,vel2;

// only used if interpol_mode>0
    void particle_interpolate(std::vector<particle_sim> &p, double frac);

    // only used if geomfile==true
    std::vector<particle_sim> p_orig;
    std::ifstream inp;
    int linecount,ninterpol,nextfile;
    int geometry_incr, geometry_skip;
    // only used if geomfile==false
    bool done;

    void fetchFiles(std::vector<particle_sim> &particle_data);

  public:
    sceneMaker (paramfile &par);

    bool getNextScene (std::vector<particle_sim> &particle_data, vec3 &campos,
      vec3 &lookat, vec3 &sky, std::string &outfile);
  };

#endif
