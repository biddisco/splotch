#ifndef SPLOTCH_HOST_H
#define SPLOTCH_HOST_H

#include "splotch/splotchutils.h"
#include "cxxsupport/walltimer.h"

void host_rendering(bool master, paramfile &params, long npart_all, 
                    arr2<COLOUR> &pic, std::vector<particle_sim> &particle_data,
                    vec3 &campos, vec3 &lookat, vec3 &sky, std::vector<COLOURMAP> &amap);

void render_new (const std::vector<particle_sim> &p, arr2<COLOUR> &pic,
  bool a_eq_e,double grayabsorb, bool nopostproc);
void render_classic (const std::vector<particle_sim> &p, arr2<COLOUR> &pic,
  bool a_eq_e,double grayabsorb, bool nopostproc);

void particle_normalize(paramfile &params, std::vector<particle_sim> &p, bool verbose);
void particle_project(paramfile &params, std::vector<particle_sim> &p,
  const vec3 &campos, const vec3 &lookat, vec3 sky);
void particle_colorize(paramfile &params, std::vector<particle_sim> &p,
  std::vector<COLOURMAP> &amap);
void particle_sort(std::vector<particle_sim> &p, int sort_type, bool verbose);

#endif
