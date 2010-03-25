#ifndef READER_H
#define READER_H

#include "splotch/splotchutils.h"

void gadget_reader(paramfile &params, std::vector<particle_sim> &p, int snr, double *time);
void gadget_millenium_reader(paramfile &params, std::vector<particle_sim> &p, int snr, double *time);
void bin_reader_tab (paramfile &params, std::vector<particle_sim> &points, float &maxr, float &minr);
void bin_reader_block (paramfile &params, std::vector<particle_sim> &points, float &maxr, float &minr);
long bin_reader_block_mpi (paramfile &params, std::vector<particle_sim> &points, float *maxr, float *minr, int mype, int npes);
void mesh_reader(paramfile &params, std::vector<particle_sim> &points, float &maxr, float &minr);
void hdf5_reader(paramfile &params, std::vector<particle_sim> &points, float &maxr, float &minr);
#endif
