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
#ifdef INTERPOLATE
#ifndef GEOMETRY_FILE
#error Splotch: interpolation without geometry file makes no sense!
#endif
#endif
    paramfile &params;
#ifdef INTERPOLATE
    std::vector<particle_sim> particle_data1,particle_data2;
    int snr_start;
    int snr1,snr2,snr1_now,snr2_now;
    double time1,time2;
#endif
#ifdef GEOMETRY_FILE
    std::vector<particle_sim> p_orig;
    std::ifstream inp;
    int linecount,ninterpol,nextfile;
    int geometry_incr, geometry_skip;
#else
    bool done;
#endif

  public:
    sceneMaker (paramfile &par);

    bool getNextScene (std::vector<particle_sim> &particle_data, vec3 &campos,
      vec3 &lookat, vec3 &sky, std::string &outfile);
  };

#endif
