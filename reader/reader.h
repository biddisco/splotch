#ifndef READER_H
#define READER_H

#include "splotch/splotchutils.h"

#ifdef INTERPOLATE
void gadget_reader(paramfile &params, std::vector<particle_sim> &p,
  std::vector<uint32> &id, int snr, double &time);
#ifdef HIGH_ORDER_INTERPOLATION
void gadget_reader(paramfile &params, std::vector<particle_sim> &p,
  std::vector<uint32> &id, std::vector<vec3f> &vel, int snr, double &time);
#endif
#else
void gadget_reader(paramfile &params, std::vector<particle_sim> &p, double &time);
#endif
void gadget_millenium_reader(paramfile &params, std::vector<particle_sim> &p, int snr, double *time);
void bin_reader_tab (paramfile &params, std::vector<particle_sim> &points);
void bin_reader_block (paramfile &params, std::vector<particle_sim> &points);
long bin_reader_block_mpi (paramfile &params, std::vector<particle_sim> &points, float *maxr, float *minr, int mype, int npes);
void mesh_reader(paramfile &params, std::vector<particle_sim> &points);
void hdf5_reader(paramfile &params, std::vector<particle_sim> &points);
#endif
