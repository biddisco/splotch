#ifndef SPLOTCH_HOST_H
#define SPLOTCH_HOST_H

#include "splotch/splotchutils.h"

void host_processing(bool master, paramfile &params, long npart_all, 
                     std::vector<particle_sim> &particle_data, vec3 &campos, vec3 &lookat, vec3 &sky, std::vector<COLOURMAP> &amap);

#endif
