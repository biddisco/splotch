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
