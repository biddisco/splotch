/*
CuPlicy is the class that knows the overall state of cuda application.
All 'magic nunbers' are out of this class.
*/

#pragma once

#include "cxxsupport/paramfile.h"

class CuPolicy
{
private:
	paramfile	*pParam;

public:
	CuPolicy(void);
	CuPolicy(paramfile *pParam);
	~CuPolicy(void);
};
