#ifndef SPLOTCH_HOST_H
#define SPLOTCH_HOST_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "splotch/splotchutils.h"

#ifdef SPLVISIVO
void host_rendering (paramfile &params, std::vector<particle_sim> &particle_data,
  arr2<COLOUR> &pic, const vec3 &campos, const vec3 &lookat, const vec3 &sky,
  std::vector<COLOURMAP> &amap, float b_brightness, VisIVOServerOptions &opt);
#else
void host_rendering (paramfile &params, std::vector<particle_sim> &particle_data,
  arr2<COLOUR> &pic, const vec3 &campos, const vec3 &lookat, const vec3 &sky,
  std::vector<COLOURMAP> &amap, float b_brightness);
#endif

void particle_colorize(paramfile &params, std::vector<particle_sim> &particle_data,
  std::vector<COLOURMAP> &amap);

void particle_project(paramfile &params, std::vector<particle_sim> &particle_data,
  const vec3 &campos, const vec3 &lookat, vec3 sky);

void particle_sort(std::vector<particle_sim> &particle_data, int sort_type, bool verbose);

void particle_eliminate(paramfile &params, std::vector<particle_sim> &particle_data, tsize &npart_all);

void render_new (std::vector<particle_sim> &particle_data, arr2<COLOUR> &pic,
  bool a_eq_e, float32 grayabsorb);

#endif
