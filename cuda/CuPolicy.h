/*
CuPlicy is the class that knows the overall state of cuda application.
All 'magic nunbers' are out of this class.
*/

#pragma once

#include "cxxsupport/paramfile.h"
#include "splotch_cuda.h"
#include "cuda.h"
#include <cutil_inline.h>

class CuPolicy
{
private:
	paramfile	*pParam;

public:
	CuPolicy(void);
	CuPolicy(paramfile *pParam);
	~CuPolicy(void);

	int		GetSizeDPD(unsigned int n);
	void	GetDimsRange(dim3 *dimGrid, dim3 *dimBlock);
};
