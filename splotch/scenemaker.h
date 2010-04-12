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
#ifdef INTERPOLATE
    std::vector<particle_sim> particle_data1,particle_data2;
    int snr_start;
    int snr1,snr2,snr1_now,snr2_now;
    double time1,time2;
#endif

    // only used if geomfile==true
    std::vector<particle_sim> p_orig;
    std::ifstream inp;
    int linecount,ninterpol,nextfile;
    int geometry_incr, geometry_skip;
    // only used if geomfile==false
    bool done;

  public:
    sceneMaker (paramfile &par);

    bool getNextScene (std::vector<particle_sim> &particle_data, vec3 &campos,
      vec3 &lookat, vec3 &sky, std::string &outfile);
  };

#endif
