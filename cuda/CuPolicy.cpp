#include ".\cupolicy.h"
//#using <mscorlib.dll>

CuPolicy::CuPolicy(void)
{
	pParam =0;
}

CuPolicy::CuPolicy(paramfile *pParam)
{
	this->pParam =pParam;
}

CuPolicy::~CuPolicy(void)
{
}

int CuPolicy::GetSizeDPD(unsigned int n)
{
	int m =sizeof(cu_particle_sim);
	return n* sizeof(cu_particle_sim);
}

void	CuPolicy::GetDimsRange(dim3 *dimGrid, dim3 *dimBlock)
{
	*dimGrid =dim3(1449);
	*dimBlock =dim3(256);
}