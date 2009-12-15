/*
Copyright things go here.
*/

#include ".\cupolicy.h"
#include <math.h>
//#using <mscorlib.dll>

CuPolicy::CuPolicy(void)
{
	m_pParam =0;
	m_blockSize =256;
//	m_fragmentCount =30000000; //30M, for the first 25000 particles
}

CuPolicy::CuPolicy(paramfile *pParam)
{
	this->m_pParam =pParam;
	m_blockSize =pParam->find<int>("block_size", 256);
}

CuPolicy::~CuPolicy(void)
{
}

//Get size of device particle datas
int CuPolicy::GetSizeDPD(unsigned int n)
{
	int m =sizeof(cu_particle_sim);//debug
	return n* sizeof(cu_particle_sim);
}

int		CuPolicy::GetMaxRegion()
{
	return m_pParam->find<int>("max_region", 1024);
}

int		CuPolicy::GetFBufSize()
{
	//just for test on my own pc now
	int	size;
	size =m_pParam->find<int>("fragment_buffer_size", 100);
	unsigned int	M =int( pow(2.0, 20) );
	return int ( size*M );
}

/*only used in test version
unsigned int CuPolicy::GetFBufSize(bool a_eq_e)
{
	if (a_eq_e)
		return m_fragmentCount *sizeof(cu_fragment_AeqE);
	else
		return m_fragmentCount *sizeof(cu_fragment_AneqE);
}
*/
void	CuPolicy::GetDims1(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
{
	*dimBlock =dim3(m_blockSize);
	unsigned int nBlock;
	if ( n %  m_blockSize == 0) //like 512
		nBlock =n/m_blockSize;
	else
		nBlock =n/m_blockSize +1;
	*dimGrid =dim3(nBlock);
	
}

void	CuPolicy::GetDimsPostProcess(int xres, int yres,
									 dim3 *dimGrid, dim3 *dimBlock)
{
	int	n =xres*yres;
	GetDims1 (n, dimGrid, dimBlock);
}

void	CuPolicy::GetDimsCombine(unsigned int minx, unsigned int miny, 
			unsigned int maxx, unsigned int maxy,dim3 *dimGrid, dim3 *dimBlock)
{
	int n =(maxx-minx)*(maxy-miny);
	GetDims1 (n, dimGrid, dimBlock);
}


void	CuPolicy::GetDimsRange(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
{
	GetDims1(n, dimGrid, dimBlock);
}

void	CuPolicy::GetDimsColorize(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
{
	GetDims1(n, dimGrid, dimBlock);}

void	CuPolicy::GetDimsRender(unsigned int n, dim3 *dimGrid, dim3 *dimBlock)
{
	GetDims1(n, dimGrid, dimBlock);
}


