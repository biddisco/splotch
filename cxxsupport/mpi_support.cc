#include "mpi_support.h"

MPI_Manager mpiMgr;

#ifdef USE_MPI
#include "mpi.h"
#else
#include <cstring>
#endif

using namespace std;

#ifdef USE_MPI

namespace {

MPI_Datatype ndt2mpi (NDT type)
  {
  switch (type)
    {
    case NAT_CHAR: return MPI_CHAR;
    case NAT_INT: return MPI_INT;
    case NAT_UINT: return MPI_UNSIGNED;
    case NAT_LONG: return MPI_LONG;
    case NAT_ULONG: return MPI_UNSIGNED_LONG;
    case NAT_LONGLONG: return MPI_LONG_LONG;
    case NAT_ULONGLONG: return MPI_UNSIGNED_LONG_LONG;
    case NAT_FLOAT: return MPI_FLOAT;
    case NAT_DOUBLE: return MPI_DOUBLE;
    case NAT_LONGDOUBLE: return MPI_LONG_DOUBLE;
    default: planck_fail ("Unsupported type");
    }
  }
MPI_Op op2mop (MPI_Manager::redOp op)
  {
  switch (op)
    {
    case MPI_Manager::Min: return MPI_MIN;
    case MPI_Manager::Max: return MPI_MAX;
    case MPI_Manager::Sum: return MPI_SUM;
    default: planck_fail ("unsupported reduction operation");
    }
  }

} // unnamed namespace

#endif

#ifdef USE_MPI

void MPI_Manager::gatherv_helper1_m (int nval_loc, arr<int> &nval,
  arr<int> &offset, int &nval_tot) const
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

#else

void MPI_Manager::gatherv_helper1_m (int nval_loc, arr<int> &nval,
  arr<int> &offset, int &nval_tot) const
  {
  nval.alloc(1);
  nval[0]=nval_tot=nval_loc;
  offset.alloc(1);
  offset[0]=0;
  }

#endif

#ifdef USE_MPI

MPI_Manager::MPI_Manager ()
  {
  MPI_Init(0,0);
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  }
MPI_Manager::~MPI_Manager ()
  { MPI_Finalize(); }

int MPI_Manager::num_ranks() const
  { int res; MPI_Comm_size(MPI_COMM_WORLD, &res); return res; }
int MPI_Manager::rank() const
  { int res; MPI_Comm_rank(MPI_COMM_WORLD, &res); return res; }
bool MPI_Manager::master() const
  { return (rank() == 0); }

#else

MPI_Manager::MPI_Manager () {}
MPI_Manager::~MPI_Manager () {}

int MPI_Manager::num_ranks() const { return 1; }
int MPI_Manager::rank() const { return 0; }
bool MPI_Manager::master() const { return true; }

#endif

#ifdef USE_MPI

void MPI_Manager::gatherRawVoid (const void *in, tsize num, void *out, NDT type)
  const
  {
  MPI_Datatype dtype = ndt2mpi(type);
  MPI_Gather(const_cast<void *>(in),1,dtype,out,num,dtype,0,MPI_COMM_WORLD);
  }
void MPI_Manager::gathervRawVoid (const void *in, tsize num, void *out,
  const int *nval, const int *offset, NDT type) const
  {
  MPI_Datatype dtype = ndt2mpi(type);
  MPI_Gatherv(const_cast<void *>(in),num,dtype,out,const_cast<int *>(nval),
    const_cast<int *>(offset),dtype,0,MPI_COMM_WORLD);
  }

void MPI_Manager::allreduceRawVoid (const void *in, void *out, NDT type,
  tsize num, redOp op) const
  {
  MPI_Allreduce (const_cast<void *>(in),out,num,ndt2mpi(type),op2mop(op),
    MPI_COMM_WORLD);
  }
void MPI_Manager::allreduceRawVoid_inplace (void *data, NDT type, tsize num,
  redOp op) const
  {
  arr<char> tmp(num);
  MPI_Allreduce (data,&tmp[0],num,ndt2mpi(type),op2mop(op),MPI_COMM_WORLD);
  memcpy (data,&tmp[0],num*ndt2size(type));
  }
void MPI_Manager::reduceRawVoid (const void *in, void *out, NDT type, tsize num,
  redOp op, int root) const
  {
  MPI_Reduce (const_cast<void *>(in),out,num,ndt2mpi(type),op2mop(op), root,
    MPI_COMM_WORLD);
  }

void MPI_Manager::bcastRawVoid (void *data, NDT type, tsize num, int root) const
  { MPI_Bcast (data,num,ndt2mpi(type),root,MPI_COMM_WORLD); }

#else

void MPI_Manager::gatherRawVoid (const void *in, tsize num, void *out, NDT type)
  const
  { memcpy (out, in, num*ndt2size(type)); }
void MPI_Manager::gathervRawVoid (const void *in, tsize num, void *out,
  const int *, const int *, NDT type) const
  { memcpy (out, in, num*ndt2size(type)); }

void MPI_Manager::allreduceRawVoid (const void *in, void *out, NDT type,
  tsize num, redOp) const
  { memcpy (out, in, num*ndt2size(type)); }
void MPI_Manager::allreduceRawVoid_inplace (void *, NDT, tsize, redOp) const
  {}
void MPI_Manager::reduceRawVoid (const void *in, void *out, NDT type, tsize num,
  redOp, int) const
  { memcpy (out, in, num*ndt2size(type)); }

void MPI_Manager::bcastRawVoid (void *, NDT, tsize, int) const
  {}

#endif
