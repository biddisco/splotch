/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2009 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_MPI_SUPPORT_H
#define PLANCK_MPI_SUPPORT_H

#ifdef USE_MPI
#include "mpi.h"
#include "arr.h"

template <typename T> inline MPI_Datatype getMpiType();
template<> inline MPI_Datatype getMpiType<char>()
  { return MPI_CHAR; }
template<> inline MPI_Datatype getMpiType<int>()
  { return MPI_INT; }
template<> inline MPI_Datatype getMpiType<unsigned int>()
  { return MPI_UNSIGNED; }
template<> inline MPI_Datatype getMpiType<long>()
  { return MPI_LONG; }
template<> inline MPI_Datatype getMpiType<unsigned long>()
  { return MPI_UNSIGNED_LONG; }
template<> inline MPI_Datatype getMpiType<long long>()
  { return MPI_LONG_LONG; }
template<> inline MPI_Datatype getMpiType<unsigned long long>()
  { return MPI_UNSIGNED_LONG_LONG; }
template<> inline MPI_Datatype getMpiType<float>()
  { return MPI_FLOAT; }
template<> inline MPI_Datatype getMpiType<double>()
  { return MPI_DOUBLE; }
template<> inline MPI_Datatype getMpiType<long double>()
  { return MPI_LONG_DOUBLE; }

class MPI_Manager
  {
  private:
    void gatherv_helper1_m (int nval_loc, arr<int> &nval, arr<int> &offset,
      int &nval_tot)
      {
      gather_m (nval_loc, nval);
      nval_tot=0;
      for (tsize i=0; i<nval.size(); ++i)
        nval_tot+=nval[i];
      offset.alloc(num_ranks());
      offset[0]=0;
      for (tsize i=1; i<offset.size(); ++i)
        offset[i]=offset[i-1]+nval[i-1];
      }

  public:
    ~MPI_Manager() {}

    void startup (int &argc, char **&argv)
      {
      MPI_Init(&argc,&argv);
      MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
      }
    void shutdown ()
      { MPI_Finalize(); }
    int num_ranks() const
      { int res; MPI_Comm_size(MPI_COMM_WORLD, &res); return res; }
    int rank() const
      { int res; MPI_Comm_rank(MPI_COMM_WORLD, &res); return res; }
    bool master() const
      { return (rank() == 0); }
    void barrier() const
      { MPI_Barrier (MPI_COMM_WORLD); }

#if 0
    void gather_m_rawVoid (const void *in, tsize insize, void *out, PDT type)
      {
      MPI_Gather (const_cast<void *>(in),
        insize,pdt2mpi(type),out,insize,pdt2mpi(type),0,MPI_COMM_WORLD);
      }
    void gather_s_rawVoid (const void *in, tsize insize, PDT type)
      {
      MPI_Gather (const_cast<void *>(in),
        insize,pdt2mpi(type),0,insize,pdt2mpi(type),0,MPI_COMM_WORLD);
      }
#endif
    template<typename T> void gather_m (const T &in, arr<T> &out)
      {
      MPI_Datatype dtype = getMpiType<T>();
      out.alloc(num_ranks());
      MPI_Gather(const_cast<T *>(&in),1,dtype,&out[0],1,dtype,0,MPI_COMM_WORLD);
      }
    template<typename T> void gather_s (const T &in)
      {
      MPI_Datatype dtype = getMpiType<T>();
      MPI_Gather(const_cast<T *>(&in),1,dtype,0,1,dtype,0,MPI_COMM_WORLD);
      }

    template<typename T> void gatherv_m (const arr<T> &in, arr<T> &out)
      {
      int nval_loc = in.size(), nval_tot;
      MPI_Datatype dtype = getMpiType<T>();
      arr<int> nval, offset;
      gatherv_helper1_m (nval_loc,nval,offset,nval_tot);
      out.alloc(nval_tot);
      MPI_Gatherv(const_cast<T *>(&in[0]),nval_loc,dtype,&out[0],&nval[0],
        &offset[0],dtype,0,MPI_COMM_WORLD);
      }
    template<typename T> void gatherv_s (const arr<T> &in)
      {
      int nval_loc = in.size();
      MPI_Datatype dtype = getMpiType<T>();
      gather_s (nval_loc);
      MPI_Gatherv(const_cast<T *>(&in[0]),nval_loc,dtype,0,0,0,dtype,0,
        MPI_COMM_WORLD);
      }
    template<typename T> void gatherv (const arr<T> &in, arr<T> &out)
      { master() ? gatherv_m(in,out) : gatherv_s(in); }

    template<typename T> void gatherv_m (const arr2<T> &in, arr2<T> &out)
      {
      int nval_loc = in.size(), nval_tot;
      MPI_Datatype dtype = getMpiType<T>();
      arr<int> nval, offset;
      gatherv_helper1_m (nval_loc, nval, offset, nval_tot);
      out.alloc(nval_tot/in.size2(),in.size2());
      MPI_Gatherv(const_cast<T *>(&in[0][0]),nval_loc,dtype,&out[0][0],&nval[0],
        &offset[0],dtype,0,MPI_COMM_WORLD);
      }
    template<typename T> void gatherv_s (const arr2<T> &in)
      {
      int nval_loc = in.size();
      MPI_Datatype dtype = getMpiType<T>();
      gather_s (nval_loc);
      MPI_Gatherv(const_cast<T *>(&in[0][0]),nval_loc,dtype,0,0,0,dtype,0,
        MPI_COMM_WORLD);
      }

    template<typename T> void broadcast_raw (T *ptr, tsize nvals)
      { MPI_Bcast (ptr,nvals,getMpiType<T>(),0,MPI_COMM_WORLD); }
    template<typename T> void broadcast (T &val)
      { MPI_Bcast (&val,1,getMpiType<T>(),0,MPI_COMM_WORLD); }
    template<typename T> void allreduce_sum (T &val)
      {
      T val2;
      MPI_Allreduce (&val,&val2,1,getMpiType<T>(),MPI_SUM,MPI_COMM_WORLD);
      val=val2;
      }
    template<typename T> void allreduce_sum_raw (T *ptr, tsize nvals)
      {
      arr<T> v2(nvals);
      MPI_Allreduce (ptr,&v2[0],nvals,getMpiType<T>(),MPI_SUM,MPI_COMM_WORLD);
      for (tsize i=0; i<nvals; ++i) ptr[i]=v2[i];
      }
    template<typename T> void allreduce_min_raw (T *ptr, tsize nvals)
      {
      arr<T> v2(nvals);
      MPI_Allreduce (ptr,&v2[0],nvals,getMpiType<T>(),MPI_MIN,MPI_COMM_WORLD);
      for (tsize i=0; i<nvals; ++i) ptr[i]=v2[i];
      }
    template<typename T> void allreduce_min (T &val)
      {
      T val2;
      MPI_Allreduce (&val,&val2,1,getMpiType<T>(),MPI_MIN,MPI_COMM_WORLD);
      val=val2;
      }
    template<typename T> void allreduce_max_raw (T *ptr, tsize nvals)
      {
      arr<T> v2(nvals);
      MPI_Allreduce (ptr,&v2[0],nvals,getMpiType<T>(),MPI_MAX,MPI_COMM_WORLD);
      for (tsize i=0; i<nvals; ++i) ptr[i]=v2[i];
      }
    template<typename T> void allreduce_max (T &val)
      {
      T val2;
      MPI_Allreduce (&val,&val2,1,getMpiType<T>(),MPI_MAX,MPI_COMM_WORLD);
      val=val2;
      }
  };

#else

class MPI_Manager
  {
  public:
    void startup (int &, char **&) {}
    void shutdown () {}
    int num_ranks() const { return 1; }
    int rank() const { return 0; }
    bool master() const { return true; }
    void barrier() const {}

    template<typename T> void gatherv (const arr<T> &in, arr<T> &out)
      { out=in; }
    template<typename T> void gatherv_m (const arr2<T> &in, arr2<T> &out)
      { out=in; }
    template<typename T> void gatherv_s (const arr2<T> &) {}
    template<typename T> void broadcast_raw (T *, tsize) {}
    template<typename T> void broadcast (T &) {}
    template<typename T> void allreduce_sum_raw (T *, tsize) {}
    template<typename T> void allreduce_sum (T &) {}
    template<typename T> void allreduce_min_raw (T *, tsize) {}
    template<typename T> void allreduce_min (T &) {}
    template<typename T> void allreduce_max_raw (T *, tsize) {}
    template<typename T> void allreduce_max (T &) {}
  };

#endif

extern MPI_Manager mpiMgr;

#endif
